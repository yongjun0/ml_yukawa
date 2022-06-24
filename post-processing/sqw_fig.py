import sys
import numpy as np
import matplotlib.pyplot as plt

kappa = int(sys.argv[1])
gamma = int(sys.argv[2])

N = 10000
Npd = 8001
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
    plt.plot(w, sqw[0, :])
    plt.plot(w, sqw[5, :])
    plt.plot(w, sqw[10, :])
    plt.plot(w, sqw[15, :])
    plt.plot(w, sqw[20, :])
    plt.plot(w, sqw[25, :])
    plt.plot(w, sqw[30, :])
    plt.xlim(xmin, xmax)
    plt.show()
