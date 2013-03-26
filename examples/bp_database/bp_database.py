import tempfile
import subprocess
import gzip
import urllib
import os.path
import os
import glob
import numpy as np
import scipy.cluster.vq as vq
import matplotlib.pyplot as plt

# Download PDBs #
pdb_list = '../../database/pdb_list/DNA_2.0_noprot.txt'
server_name = 'http://www.rcsb.org/pdb/files/'
for pdb_id in open(pdb_list) :
    filename = '%s.pdb.gz' % pdb_id[:-1].lower()
    filename_pdb = '%s.pdb' % pdb_id[:-1].lower()
    addr = '%s%s' %( server_name, filename )
    if  ( (not os.path.exists(filename)) and
        (not os.path.exists(filename_pdb)) ):
        print 'Downloading', filename
        urllib.urlretrieve( addr, filename ) # Download the files

        #Uncompression
    if os.path.exists(filename):
        f = gzip.open(filename)
        out = open(filename_pdb, 'w')
        out.write( f.read() )
        out.close()
        f.close()
        os.remove(filename)

# Run 3-DNA on the pdb files (need to get 3DNA installed) #
def run_3dna( name ):
    temp = tempfile.NamedTemporaryFile()
    subprocess.check_call('find_pair %s %s' % (name, temp.name), shell=True)
    try :
        subprocess.check_call('analyze %s' % temp.name, shell=True)
    except :
        pass

for filename in glob.glob('*.pdb'):
    if not os.path.exists(filename.replace('.pdb', '.out')):
        run_3dna( filename )

# Load the bp-step data from 3DNA output #
seq_list = ['AA','AT','AG','AC','GA','GC','GT','GG','TT','TA','TC','TG','CC','CG','CA','CT']
params = [[] for i in xrange(16)]
for out_3dna in glob.glob('*.out'):
    is_readline1 = False
    is_readline2 = False
    wc_bp = []
    for line in open(out_3dna):
        if is_readline1:
            if len(line) < 4 or line[:4] == '****':
                is_readline1 = False
            else :
                bp_num = int( line[:4] )
                bp_tag = line[34:39]
                if bp_tag in ['-----', 'x----', '----x']:
                    wc_bp.append( bp_num )
        if 'Strand I                    Strand II' in line:
            is_readline1 = True

        if is_readline2:
            if ('~~~~' in line) or ('****' in line):
                break
            if ('----' not in line):
                bp_num = int( line[:4] )
                if bp_num in wc_bp and (bp_num + 1) in wc_bp:
                    try:
                        bp_type = line[5:7].upper()
                        idx = seq_list.index(bp_type)
                        params[idx].append( np.fromstring( line[10:], sep=' ' ) )
                    except:
                        pass
        if 'Shift     Slide' in line:
            is_readline2 = True

# data filtering #
params_filtered = []
for params_i in params:
    params_i = np.array( params_i )
    params_i = params_i[ params_i[:,5] > 5 ]
    params_i = params_i[ params_i[:,2] < 5.5 ]
    data_avg = np.average(params_i, axis=0)
    data_std = np.std(params_i, axis=0)
    data_min = data_avg - 4 * data_std
    data_max = data_avg + 4 * data_std
    params_i = params_i[ np.all( params_i >= data_min[np.newaxis,:], axis=1 ) ]
    params_i = params_i[ np.all( params_i <= data_max[np.newaxis,:], axis=1 ) ]
    params_filtered.append(params_i)

params_name = ['shift', 'slide', 'rise', 'tilt', 'roll', 'twist']
all_params = np.vstack( params_filtered ) #list without sequence dependence

# Plot the histogram before clustering #
fig = plt.figure(figsize=(12,8))
for i in xrange(6):
    subplt = fig.add_subplot(2,3,i+1)
    plt.hist( all_params[:,i], 30, histtype='step' )
    plt.title(params_name[i])

# Clustering #
n_cluster = 3
whitened = vq.whiten( all_params )
avg1 = np.average( whitened[:,1] )
whitened[:,1] = (whitened[:,1] - avg1) * 1.2 + avg1
centr,f = vq.kmeans(whitened, n_cluster)
reject_cluster = centr[ centr[:,1] < -1 ] # cluster with center slide < -1 -> A-DNA
acpt_cluster = centr[ centr[:,1] >= -1 ] # the other two B-DNA clusters
centr = np.vstack( (reject_cluster, acpt_cluster) )

# Plot after clustering #
idx,_ = vq.vq(whitened, centr)
fig = plt.figure(figsize=(12,8))
for i in xrange(6):
    subplt = fig.add_subplot(2,3,i+1)
    for j in xrange(n_cluster):
        plt.hist( all_params[idx==j][:,i], 30, histtype='step' )
    plt.title(params_name[i])

# Remove the A-DNA cluster from the sequence-dependent params #
start = 0
final_params = []
for params in params_filtered:
    end = start + params.shape[0]
    new_idx = idx[start:end]
    params = params[new_idx != 0] # remove data belonging to the reject cluster
    params[:,3:] = np.radians( params[:,3:] ) # Convert angles to radians
    final_params.append( params )
    start = end

# Output the database file #
expr = "np.savez( 'DNA_2.0_noprot',"
for i in xrange(16):
    expr += " %s = final_params[%d], " % (seq_list[i], i)
expr += ')'
eval(expr)



plt.show()
