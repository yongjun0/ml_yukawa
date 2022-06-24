import numpy as np
import numba as nb
import sys
import re
import glob
import matplotlib.pyplot as plt

def add_values_in_dict(sample_dict, key, list_of_values):
    ''' Append multiple values to a key in 
        the given dictionary '''
    if key not in sample_dict:
        sample_dict[key] = list()
    sample_dict[key].extend(list_of_values)
    return sample_dict
###########
t_kappa = 0.5
t_gamma = 11
###########

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
w = np.arange(0, dw*Npd, dw)
wmax = 10.0
idx_wmax = int(20/dw)

path0 = "data_sample/sqw/"
file_list = sorted(glob.glob(path0+"*"))
data = np.load(file_list[0])

t_sqw = np.empty_like(data)
#t_sqw = np.zeros((Nq, Npd))

kg = {}
#kg = add_values_in_dict(kg, '1', [22])
#kg = add_values_in_dict(kg, '1', [21])
for ifile in file_list:
    kappa  = re.search('(([0-9]*[.])?[0-9]+)', ifile)
    gamma  = re.search('((?<=-)([0-9]*[.])?[0-9]+)', ifile)
    kg_dict = add_values_in_dict(kg, kappa[0], [gamma[0]])


den_weight = 0.0
nom_weight = 0.0
weight = 0.0

print(t_kappa, t_gamma)
print("======")
for key in kg_dict:
    ikappa = key
    for igamma in kg_dict[key]:
        if(ikappa != str(t_kappa) or igamma != str(t_gamma)):
            #print(ikappa, igamma)
            fn = path0+key+"-"+igamma+".npy" 
            igamma = float(igamma)
            ikappa = float(ikappa)
            
            dis = np.sqrt((t_gamma-igamma)**2 + (t_kappa-ikappa)**2)
            dis_sq = 10*dis**2
            den_weight += np.exp(-dis_sq)
            nom_weight = np.exp(-dis_sq)
            data = np.load(fn)
            t_sqw += nom_weight*data

t_sqw /=den_weight
#data = np.load(path0+'0.5-11.npy')
data = np.load("../Yukawa/kappa_0.5/gamma_11/data/sqw_avg.npy")
plt.figure()
plt.plot(w, t_sqw[0, Npd:], 'r')
plt.plot(w, 2*data[0, Npd:], 'b')
plt.xlim(0, 2)
plt.show()

