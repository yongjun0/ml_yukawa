import sys
import numpy as np
import matplotlib.pyplot as plt

kappa = str(sys.argv[1])
gamma = str(sys.argv[2])

N = 10000
Npd = 8000
dt = 0.1
q_max = 30

T = Npd * dt
L = (4*np.pi*N/3)**(1./3)
dq = 2*np.pi/L
Nq = int(q_max/dq)
 
#frequency range
dw = 2*np.pi/T

w = np.arange(-dw*Npd, dw*Npd, dw)

path0 = "../Yukawa/kappa_"+str(kappa)+"/gamma_"+str(gamma)+"/data/"
fn = "sqw_avg.npy"
data = np.load(path0+fn)

sqw = data
print(sqw.shape)
print(Nq, len(w))
print(w)
xmax = 5
xmin = -xmax
if(1):
    plt.figure()
    plt.plot(w, sqw[10, :])
    plt.xlim(xmin, xmax)
    plt.xlabel("$\omega/ \omega_i$", fontsize = 16)
    #plt.ylabel("$S(q, \omega/ \omega_i)$", fontsize = 16)
    plt.show()
