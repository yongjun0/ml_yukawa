import numpy as np
import numba as nb
import sys
import re
import glob
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def add_values_in_dict(sample_dict, key, list_of_values):
    ''' Append multiple values to a key in 
        the given dictionary '''
    if key not in sample_dict:
        sample_dict[key] = list()
    sample_dict[key].extend(list_of_values)
    return sample_dict
###########
t_kappa = float(sys.argv[1])
t_gamma = int(sys.argv[2])
#t_kappa = 2.2
#t_gamma = 80 
###########

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
#data = np.load("../Yukawa/kappa_2.1/gamma_80/data/sqw_avg.npy")
data = np.load("data_sample/sqw/"+str(t_kappa)+"-"+str(t_gamma)+".npy")


ax = np.zeros(4, object)
fig, ((ax[0], ax[1]), (ax[2], ax[3])) = plt.subplots(2, 2, sharex = 'col', figsize = (12, 6))

tx = 0.65
ty = 0.45
fl = 16

for ij in range(4):
    nq = (ij)*10
    ax[ij].plot(w, t_sqw[nq, Npd:], 'r', label='Guess')
    ax[ij].plot(w, data[nq, Npd:], 'b', label='True')

    label = f"q = {((nq+1)*dq):3.2f}"
    ax[ij].text(tx, ty, r'$'+label+'$', fontsize = fl, transform=ax[ij].transAxes)
    ax[ij].set_xlim([0, 2])
    ax[ij].tick_params(axis='both', which='major', labelsize = fl)
    ax[ij].legend(loc='upper right', fontsize = fl)


ax[2].set_xlabel("$\omega$", fontsize = fl)
ax[3].set_xlabel("$\omega$", fontsize = fl)
plt.subplots_adjust(hspace=0)
plt.suptitle(f"$S(q, \omega), \kappa = {t_kappa}\; \Gamma ={t_gamma}$")


#plt.show()
#pdf=PdfPages(f'figures/{t_kappa}-{t_gamma}.pdf')
#pdf.savefig(fig)
#pdf.close()
plt.savefig(f'figures/{t_kappa}-{t_gamma}.png')
