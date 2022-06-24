'''
#y1: PRE: 
Ionic Transport in High Energy-Density Matter by LG Stanton, and MS Murillo

#y2: High Energy Density Physics
 Viscosity estimates of liquid metals and warm dense matter using the Yukawa reference system by MS Muriilo
# Yongjun Choi
'''
import matplotlib.pyplot as plt
import numpy as np

x = np.arange(0,3.1,0.1)

y1 = 5.97+1.93*x+1.16*x*x+1.44*x*x*x
y2 = 171.8+82.8*(np.exp(0.565*x**1.38) - 1)
#for i in range(len(x)):
#  print np.arange(y1[i], y2[i], (y2[i]-y1[i])/4)
for i in range(len(x)):
    print(f"{x[i]:.2f}, {y1[i]:.2f}, {y2[i]:.2f}")


kappa = [0, 1, 2, 3]
gamma_eff = [1, 10, 100, 150]
gamma = np.array([[10, 50, 150], [14, 72, 217], \
        [31, 158, 476], [100, 503, 1510]])

#kappa = [0, 0.5, 1, 1.5, 2, 2.5, 3]
#gamma = np.array([[10, 50, 150], [11, 61, 183], [14, 72, 271], [22, 116, 328], \
#        [31, 158, 476], [60, 328, 890], [100, 503, 1510]])

#for i in range(len(kappa)):
#  for j in range(len(gamma_eff)):
#    gamma = gamma_eff[j]*np.exp(kappa[i])

fs = 15
plt.figure(1)
plt.plot(x, y1, 'b')
plt.plot(x, y2, 'b')



#    gamma = gamma_eff[j]*np.exp(kappa[i])
plt.subplot(211)
for i in range(len(kappa)):
  for j in range(3):    
    plt.semilogy(x, y1)
    plt.semilogy(x, y2)
    plt.semilogy(kappa[i], gamma[i, j], marker='o', color='red')
plt.subplot(212)
for i in range(len(kappa)):
  for j in range(3):    
    plt.plot(kappa[i], gamma[i, j], marker='o', color='red')
    plt.plot(x, y1)
    plt.plot(x, y2)
plt.xlabel(r"$\kappa$", fontsize=fs)
plt.ylabel(r"$\Gamma$", fontsize=fs)
plt.show()
