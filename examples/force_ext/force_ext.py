# -*- coding: UTF-8 -*-

import subprocess
import numpy as np
import matplotlib.pyplot as plt
from helixmc.pose import HelixPose

# Run HelixMC #
cmdline  = 'helixmc-run '
cmdline += '-params_file ../../database/DNA_default.npz '
cmdline += '-n_bp 100 '
cmdline += '-n_step 10 '
cmdline += '-seq GCCG '
cmdline += '-force 5 '
cmdline += '-compute_fuller_link True '
cmdline += ' -out_frame test_run'

print 'HelixMC command line:', cmdline
subprocess.check_call( cmdline.split() )

# Data Analysis #
data = np.load('MC_data.npz')
print 'Data content:', data.files
print 'Z-coord:', data['coord_terminal'][:,2] #2 for the z-elements
print 'Avg. Z:', np.average(data['coord_terminal'][:,2]) #avg. z-extension in Å
print 'Avg. Link:', np.average(data['twist'] + data['writhe']) #avg. link in radian

# Helix Plotting #
pose = HelixPose('test_run.npz')
pose.plot_centerline(show=False) #plot the centerline
plt.title('Center line')


pose.plot_helix(show=False) #plot the entire helix
plt.title('Helix Plot')
plt.show()
