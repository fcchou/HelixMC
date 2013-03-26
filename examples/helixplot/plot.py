from helixmc.pose import HelixPose
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3

# Here we have the final frames of 3 independent HelixMC runs,
# named data0.npz, data1.npz and data3.npz. We will plot these
# data in this example.

# Load data #
pose_list = []
pose_list.append( HelixPose( 'data0.npz' ) )
pose_list.append( HelixPose( 'data1.npz' ) )
pose_list.append( HelixPose( 'data2.npz' ) )

coords = [ pose.coord / 10 for pose in pose_list ]
colors = 'kbrgycm'

# Top View #
plt.figure()
for i, coord in enumerate(coords):
    plt.plot( coord[:,0], coord[:,1], color = colors[i] )
plt.ylabel('Y (nm)')
plt.xlabel('X (nm)')
plt.title( 'Top view' )

# Side View #
plt.figure()
for i, coord in enumerate(coords):
    plt.plot( coord[:,0], coord[:,2], color = colors[i] )
plt.ylabel('Z (nm)')
plt.xlabel('X (nm)')
plt.title( 'side view' )

# 3D Plot #
fig = plt.figure()
ax = p3.Axes3D(fig)
for i, pose in enumerate(pose_list):
    pose.plot_centerline( color = colors[i], show = False, fig_ax = ax )

plt.show()

