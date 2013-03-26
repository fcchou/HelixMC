# -*- coding: UTF-8 -*-

import subprocess
import numpy as np
import matplotlib.pyplot as plt

# Run HelixMC #
cmdline  = 'helixmc-run '
cmdline += '-params_file ../../database/DNA_default.npz '
cmdline += '-n_bp 1000 '
cmdline += '-n_step 50 '
cmdline += '-force 0 '
cmdline += '-check_fuller fuller_check.out '

print 'Fuller-Check command line:', cmdline
subprocess.check_call( cmdline.split() )

# Data Plotting #
data = np.loadtxt('fuller_check.out')
data /= 2 * np.pi #Convert from radians to # of turns
plt.figure()
plt.plot(data[:,1], data[:,0], 'rx')
plt.ylabel("Fuller's Writhe (turn)" )
plt.xlabel("Exact Writhe (turn)" )
plt.show()
