# -*- coding: UTF-8 -*-

import subprocess
import numpy as np

# Quick prerun to get the avg. link at relaxed state #
cmdline = 'helixmc-run '
cmdline += '-params DNA_default.npz '
cmdline += '-n_bp 500 '
cmdline += '-n_step 20 '
cmdline += '-force 5 '
cmdline += '-compute_fuller_link '
cmdline += '-out prerun'

print 'Prerun command line:', cmdline
subprocess.check_call(cmdline.split())

data = np.load('prerun.npz')  # Load the prerun output
avg_link = np.average(data['twist'] + data['writhe'])  # avg. link in radian
# avg. z-extension in Å
print 'Avg. Z relaxed:', np.average(data['coord_terminal'][:, 2])
print 'Avg. link relaxed:', avg_link

# Link constrained run #
target_link = avg_link + 8 * (2 * np.pi)  # set the center at +8 turns
print 'Target link:', target_link

cmdline = 'helixmc-run '
cmdline += '-params DNA_default.npz '
cmdline += '-n_bp 500 '
cmdline += '-n_step 20 '
cmdline += '-force 5 '
cmdline += '-compute_fuller_link '
cmdline += '-target_link %f ' % target_link
cmdline += '-torsional_stiffness 2000 '
cmdline += '-out link_cst'

print 'Link contrained command line:', cmdline
subprocess.check_call(cmdline.split())

# Data Analysis #
data = np.load('link_cst.npz')
# avg. z-extension in Å
print 'Avg. Z (8 turns):', np.average(data['coord_terminal'][:, 2])
# avg. link in radian
print 'Avg. Link (8 turns):', np.average(data['twist'] + data['writhe'])
