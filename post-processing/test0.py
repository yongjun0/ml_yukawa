import numpy as np
import numba as nb
import sys
import re
import glob
import matplotlib.pyplot as plt

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

kappa = 1
path0="../Yukawa/kappa_1/"
data1 = np.load(path0+"gamma_11/data/sqw_avg.npy")
data2 = np.load(path0+"gamma_14/data/sqw_avg.npy")
data3 = np.load(path0+"gamma_72/data/sqw_avg.npy")
data4 = np.load(path0+"gamma_144/data/sqw_avg.npy")
data5 = np.load(path0+"gamma_217/data/sqw_avg.npy")
data6 = np.load(path0+"gamma_234/data/sqw_avg.npy")

plt.figure()
plt.subplot(121)
plt.plot(w, data1[0], 'r')
plt.plot(w, data2[0], 'g')
plt.plot(w, data3[0], 'b')
plt.plot(w, data4[0], 'r-.')
plt.plot(w, data5[0], 'g-.')
plt.xlim(0, 0.5)
plt.subplot(122)
plt.plot(w, data1[0], 'b-.')
plt.xlim(0, 0.5)
plt.show()
