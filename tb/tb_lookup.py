from src.pygmid import Lookup as lk
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# setup mpl
mpl.rcParams['axes.spines.right'] = False
mpl.rcParams['axes.spines.top'] = False
mpl.rcParams.update({"axes.grid" : True})
mpl.use("TkAgg")

#%%
# Specify appropriate path...
NCH = lk('xt018_ne5.pkl')  # load process data into pygmid lookup object

VDSs = NCH['VDS']       # lookup object has pseudo-array access to data
VGSs = np.arange(0.4, 0.6, 0.05)


# %%
# Plot ID versus VDS
ID = NCH.look_up('ID', vds=VDSs, vgs=VGSs)
plt.figure()
plt.plot(VDSs, 1e6*ID.T)
plt.ylabel(r"$I_D$ [$\mu$A]")
plt.xlabel(r"$V_{DS}$ [V]")
plt.title(r'$I_D$ vs. $V_{DS}$ for varying $V_{GS}$')
plt.legend(VGSs)
plt.show()

# %% Gm over Id methodology applied to folded-cascode OTA

common_mode_output_voltage = 0.6
output_swing = 0.8
G = 2
L0 = 50

Vds_min = (common_mode_output_voltage - output_swing/4)/2
# gm_id = 2/Vds_min
gm_id = 15

beta_max = 1/(1+G)
beta = 0.75*beta_max # first-order optimum
kappa = 0.7 # conservative optimum

# Channel length sweep
L23 = L45 = L =  np.linspace(0.06,1,100)
gm_gds2 = NCH.look_up('GM_GDS', gm_id=gm_id, vds=0.2, L=L23)
gm_gds3 = NCH.look_up('GM_GDS', gm_id=gm_id, vds=0.4, L=L23)
gm_gds4 = NCH.look_up('GM_GDS', gm_id=gm_id, vds=0.4, L=L45)
gm_gds5 = NCH.look_up('GM_GDS', gm_id=gm_id, vds=0.2, L=L45)

# Chosen length
L0_swept = beta*kappa/(1/(1+gm_gds5)/gm_gds4+1/(1+gm_gds2/3)/gm_gds3)






ID = NCH.look_up('ID', vds=VDSs, vgs=VGSs)
plt.figure()
plt.plot(VDSs, 1e6*ID.T)
plt.ylabel(r"$I_D$ [$\mu$A]")
plt.xlabel(r"$V_{DS}$ [V]")
plt.title(r'$I_D$ vs. $V_{DS}$ for varying $V_{GS}$')
plt.legend(VGSs)
plt.show()