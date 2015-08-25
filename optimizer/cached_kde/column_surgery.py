import iopro
import numpy as np

a_0p0 = iopro.genfromtxt('cached_kde.ada.0p0.csv')
a_0p25 = iopro.genfromtxt('cached_kde.ada.0p25.csv')
a_0p50 = iopro.genfromtxt('cached_kde.ada.0p5.csv')
N, D = a_0p50.shape

#a = np.hstack((a0[:,0:2], a5[:,2:]))
#np.savetxt("cached_kde.ada.hack1.csv", a, delimiter=" ")

#a = np.hstack((a5[:,0:2], a0[:,2:]))
#np.savetxt("cached_kde.ada.hack2.csv", a, delimiter=" ")

#a = np.hstack((a5[:,0].reshape((N, 1)), a0[:,1].reshape((N,1)), a5[:,2:]))
#np.savetxt("cached_kde.ada.hack3.csv", a, delimiter=" ")

# !!
# mildly adapting c0?
a = np.hstack((a_0p0[:,0].reshape((N, 1)),
               a_0p50[:,1].reshape((N, 1)),
               a_0p50[:,2].reshape((N, 1)),
               a_0p50[:,3].reshape((N, 1)),
               a_0p50[:,4].reshape((N, 1))))
np.savetxt("cached_kde.ada.hack4.csv", a, delimiter=" ")

# over-adaptation = bad
#a = np.hstack((a5[:,0].reshape((N, 1)),
#               a5[:,1].reshape((N, 1)),
#               a10[:,2].reshape((N, 1)),
#               a5[:,3].reshape((N, 1)),
#               a5[:,4].reshape((N, 1))))
#np.savetxt("cached_kde.ada.hack5.csv", a, delimiter=" ")

#a = np.hstack((a10[:,0].reshape((N, 1)),
#               a10[:,1].reshape((N, 1)),
#               a5[:,2].reshape((N, 1)),
#               a5[:,3].reshape((N, 1)),
#               a5[:,4].reshape((N, 1))))
#np.savetxt("cached_kde.ada.hack6.csv", a, delimiter=" ")

# !!
#a = np.hstack((a25[:,0].reshape((N, 1)),
#               a50[:,1].reshape((N, 1)),
#               a50[:,2].reshape((N, 1)),
#               a50[:,3].reshape((N, 1)),
#               a50[:,4].reshape((N, 1))))
#np.savetxt("cached_kde.ada.hack7.csv", a, delimiter=" ")

# !!
# alpha smaller => c0 larger?
