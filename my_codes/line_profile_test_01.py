'''
Sythetic line profile generation. Calculates optical depth and 
realtive flux or intensity for a particular line, such as CIV,
OVIII and Fe XXV. Plots a two panel figure showing tau and flux. 
'''

# importing necessary libraries and scripts

import numpy as np
from opacity_table_v3 import *           # reads XSTAR lookup table
from matplotlib import pyplot as plt  
import matplotlib      
from my_athena_reader import *           # reads Athena++ output file

####################################################################

# constants we'll need
mp  = 1.6726231e-24                      # 1.67262178e-24 g
pi  = np.pi
h   = 6.6262e-27                         # plank's constant in erg / s
C   = 3.0e10                             # 2.99792458e+10 cm / s
kb  = 1.380658e-16                       # 1.380658e-16 erg / K
mu  = 0.6
gamma = 1.6666667
nx = 1e8                                 # /cm^3

####################################################################

# reading data from XSTAR lookup table and Athena++ results

Model = 'A'
Line = 'FeXXV'

if(Model=='A'):
    hydro_file   = "Athenatab/agn1.hep18.merged.out1.01200.tab"
    hydro_file1 =  "Athenatab/agn1.hep18.merged.out2.01200.tab"
else:
    hydro_file   = "Athenatab/agn1.hep17.merged.out1.01200.tab"
    hydro_file1 =  "Athenatab/agn1.hep17.merged.out2.01200.tab"
    
if(Line=='FeXXV'):
    opacity_file = "XSTARtab/fe_xxv_6.7keV.dat"
    lines = [1.8505]
    matom = mp * 52       # mass of target atom in g
if(Line=='CIV'):
    opacity_file = "XSTARtab/c_iv_1550.dat"
    lines = [1550.77, 1548.2]
    matom = mp * 12      # mass of target atom in g
if(Line=='OVIII'):
    opacity_file = "XSTARtab/o_viii_19_1.dat"
    lines = [18.9725,18.9671]
    matom = mp * 16      # mass of target atom in g

####################################################################

# setting the basic parameters, velocity distribution space

nu_lab0 = []
for i in range(len(lines)):
    nu_lab0.append(C*1e8/lines[i])      # rest frequency of line

optable   = opacity_table(opacity_file)
hydroData = read_tab(hydro_file)
hydroData1 = read_tab(hydro_file1)

ymin = -7e7/C
ymax = 5e6/C
Nv=1000
nu_dist = np.linspace(ymin,ymax,Nv)
vel_dist = []
v = []
for j in range(Nv):
    vel_dist.append(nu_dist[j]*C)       # velocity from Doppler shift from assumed distribution
    v.append(vel_dist[j]/1e5)           # velocity in units of km/s

####################################################################
    
# extracting variables from Athena data
    
r = hydroData['x1v']
r = r.tolist()
dr = hydroData['x1v'][1:] - hydroData['x1v'][:-1] 
dr = dr.tolist()
nur = C*1e8/lines[0]
rho = hydroData['rho']
rho = rho.tolist()
vel = hydroData['vel1']
vel = vel.tolist()
vel1 = vel
vel = [-j for j in vel]                # sign of velocity needs to be inverted
press = hydroData['press']
press = press.tolist()
xis = hydroData1['user_out_var0']
xis = xis.tolist()

####################################################################

# calculating temperature and number density

temp = []
n_den = []
for k in range(len(press)):
    temp.append(press[k]*mu*mp/kb/rho[k])
    n_den.append(rho[k]/mu/mp)

####################################################################
    
# figure setup
    
fs = 20
fig, axs = plt.subplots(nrows=2,ncols=1,figsize=(7,5))  
colr = ['blue','red']
plt.subplots_adjust(left=0.12, bottom=0.08, right=0.88, top=0.95, wspace=0.1, hspace=0.)
plt.rcParams['font.family'] = 'serif' #'STIXGeneral' #'serif'
matplotlib.rcParams['font.size'] = '16'
matplotlib.rcParams['ps.fonttype'] = 42 #note: fontype 42 compatible with MNRAS style file when saving figs
matplotlib.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['axes.labelsize'] = 16
plt.rcParams['axes.linewidth'] = 1.
plt.rcParams['xtick.major.size'] = 6
plt.rcParams['xtick.minor.size'] = 3
plt.rcParams['ytick.major.size'] = 6
plt.rcParams['ytick.minor.size'] = 3
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 14
plt.rcParams['legend.numpoints'] = 1  # uses 1 symbol instead of 2
plt.rcParams['legend.frameon'] = False 
plt.rcParams['legend.handletextpad'] = 0.3

####################################################################

# starting line profile calculations

for i,snap in enumerate(lines):
    nu_lab = nu_lab0[i]
    print ("starting line profile calculations for ", snap)
    
    # val = sigma*nXSTAR, f0_array = nu_D (Doppler shifted frequency)
    # vth is thermal velocity, vi is thermal velocity normalized to c
    # k0 is the correctly normalized alpha obtained from XSTAR
    
    alpha_array = []
    f0_array = []
    k0 = []
    vth = []
    vi = []
    
    for k in range(len(xis)):
        val = optable.get_opacity(np.log10(xis[k]),np.log10(temp[k]),i)
        alpha_array.append(val*n_den[k]/nx)
        k0.append(alpha_array[k]/np.sqrt(pi))
        f0_array.append( (1. + vel[k]/C) )
        vth.append(np.sqrt(2*kb*temp[k]/matom))
        vi.append(vth[k]/C)
    
    flux = []
    tau = []
    
    for k in range(len(nu_dist)):
        s = 0
        tmp = []
        for j in range(len(dr)):
            y = nu_dist[k]
            y -= vel[j]/C
            y *= 1/vi[j]
            phi = np.exp(-y*y)
            k1 = k0[j]*phi
            tmp.append(k1)
        tmp1 = []
        for j in range(len(tmp)-1):
            tmp1.append(0.5*(tmp[j]+tmp[j+1])*dr[j])
        tau.append(np.sum(tmp1))
        flux.append(np.exp(-tau[k]))
    
    axs[0].plot(v,tau,color=colr[i],lw=2,label='%.4f'%snap)
    axs[1].plot(v,flux,color=colr[i],lw=2)

axs[0].set_ylabel(r"$\tau_\nu$", fontsize=fs,rotation=0,labelpad=14)
axs[0].set_ylim(0,max(tau)+0.2*max(tau))
axs[0].set_xlim(-700,50)
axs[0].legend(frameon=False)
axs[1].set_ylabel(r"$I_\nu$", fontsize=fs,rotation=0,labelpad=14)
axs[1].set_xlabel('v [km/s]',fontsize=fs)
axs[1].set_ylim(min(flux)-0.2*min(flux),1)
axs[1].set_xlim(-700,50)
axs[0].set_title('Model '+Model+', '+Line)

for ax in axs:
    ax.minorticks_on()
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    ax.tick_params(axis='both', direction='in',width=1.) 
    ax.tick_params(which='minor',axis='both',direction='in',width=1.) 
    if(ax==axs[0]):
        ax.tick_params(which='major',axis='both',direction='in',width=1.) 
    
plt.setp(axs[0].get_xticklabels(), visible=False)
plt.setp(axs[0].get_yticklabels()[0],visible=False)
plt.setp(axs[1].get_yticklabels()[-1],visible=False)

fn = 'Imgs/'+'Model'+Model+Line+'.png'
plt.savefig(fn,dpi=300,bbox_inches='tight')
plt.show()