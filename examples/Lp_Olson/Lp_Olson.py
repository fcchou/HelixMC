import numpy as np
from helixmc.random_step import RandomStepAgg
import Nucleic_Acids_MC as namc

# Set up the RandomStep #
database = '../../database/DNA_default.npz'
random_step = RandomStepAgg( database )

# Compute the average R and o #
n_step = 1000000
avg_R = np.zeros( (3,3) )
avg_o = np.zeros(3)
for i in xrange(n_step):
    params, o, R = random_step()
    avg_o += o
    avg_R += R
avg_o /= n_step
avg_R /= n_step

# Get the bending persistence length #
A = np.eye(4)
A[:3,:3] = avg_R
A[3,:3] = avg_o
Lp = np.linalg.matrix_power(A, 10000)[3,2]
print 'Bending persistence length:', (Lp / 10), 'nm'

