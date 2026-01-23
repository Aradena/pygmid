from src.pygmid import Lookup as lk
from src.pygmid import interp1
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from prefixed import Float
import scipy as sp

# setup mpl
mpl.rcParams['axes.spines.right'] = False
mpl.rcParams['axes.spines.top'] = False
mpl.rcParams.update({"axes.grid" : True})
mpl.use("TkAgg")

#%% Load process data into pygmid lookup object
# NCH = lk('xt018_ne5.pkl')
# PCH = lk('xt018_pe5.pkl')
NCH = lk('65nch.mat')
PCH = lk('65pch.mat')

# %% Gm over Id methodology applied to folded-cascode OTA

common_mode_output_voltage = 0.6
output_swing = 0.8
class s: pass  # structure to collect specs
s.G = 2

Vds_min = (common_mode_output_voltage - output_swing/4)/2
# gm_id = 2/Vds_min
gm_id = 15

beta_max = 1/(1+s.G)
beta = 0.75*beta_max # first-order optimum
kappa = 0.7 # conservative optimum

# Channel length sweep
# L23 = L45 = L =  np.linspace(0.5,5,100)
L23 = L45 = L =  np.linspace(0.06,1,100)
gm_gds2 = NCH.look_up('GM_GDS', gm_id=gm_id, vds=0.2, L=L23)
gm_gds3 = NCH.look_up('GM_GDS', gm_id=gm_id, vds=0.4, L=L23)
gm_gds4 = PCH.look_up('GM_GDS', gm_id=gm_id, vds=0.4, L=L45).reshape(-1,1)
gm_gds5 = PCH.look_up('GM_GDS', gm_id=gm_id, vds=0.2, L=L45)[:, np.newaxis]
# alternatives for transforming a (horizontal) 1D vector to a vertical vector
# requires to add a dimension

# Chosen length
L0_swept = beta*kappa/(1/(1+gm_gds5)/gm_gds4+1/(1+gm_gds2/3)/gm_gds3)

# This custom formatter removes trailing zeros, e.g. "1.0" becomes "1", and
def fmt(x):
    s = f"{x:.1f}"
    if s.endswith("0"):
        s = f"{x:.0f}"
    # return rf"{s} \%" if plt.rcParams["text.usetex"] else f"{s} %"
    return rf"{s}"

# fig, ax = plt.subplots()
# levels= np.linspace(10,100,10)
# CS = ax.contour(L23, L45, L0_swept) #, levels)
# ax.set_ylabel(r"$L_{4,5}$ [$\mu$m]")
# ax.set_xlabel(r"$L_{2,3}$ [$\mu$m]")
# ax.set_title(r'Low frequency loop gain $L_0$ as a function of channel lengths')
# ax.clabel(CS, CS.levels, fmt=fmt, fontsize=10)
# plt.show()
# plt.close()

# L23 = L45 = 2  # 300 < L0 < 400  (50 with book techno in 65 nm)
L23 = L45 = 0.4  # L0 ~ 60


gmb_gm3 = NCH.look_up('GMB_GM', gm_id=gm_id, vds=0.4, vsb=0.2, L=L23)
gm_css3 = NCH.look_up('GM_CSS', gm_id=gm_id, vds=0.4, vsb=0.2, L=L23)
cdd_css3 = NCH.look_up('CDD_CSS', gm_id=gm_id, vds=0.4, vsb=0.2, L=L23)
cdd_w3 = NCH.look_up('CDD_CSS', gm_id=gm_id, vds=0.4, vsb=0.2, L=L23)
cdd_w2 = NCH.look_up('CDD_CSS', gm_id=gm_id, vds=0.2, L=L23) # vsb=0.2,
fp2 = 1/2/np.pi*gm_css3*(1+gmb_gm3)/(1+2*cdd_css3*2*(cdd_w2/cdd_w3))
print(f'Non dominant pole fp2 = {Float(fp2):!.2h}Hz') #%.2E Hz' % fp2)

# %% Design for a given settling time & noise budget, minimise power dissipation
# s.ts = 100e-9
s.ts = 5e-9
s.ed = 0.1e-2
s.fu1 = 1/2/np.pi*np.log(1/s.ed)/s.ts
print("fp2/s.fu1 = %.2f gives a phase margin of about "%(fp2/s.fu1))

s.vod = 400e-6 # total integrated noise
s.T = 300 # Kelvin
s.fan_out = 0.5

def folded_cascode(s, d):
    m1 = {}

    # Extra hypothesis
    gamma1 = gamma5 = gamma2 = 0.7
    gm_id5 = gm_id2 = gm_id
    gm1_gm3 = 1

    # 2D sweep for finding beta with a physical reality
    for beta in d.beta:
        # Exclude last point for removing issues with interpolation to close to upper boundary

        gm_id1 = d.gm_id1
        # Excess noise factor
        alpha = 2*gamma1*(1+gamma5/gamma1*gm_id5/d.gm_id1+2*gamma2/gamma1*gm_id2/d.gm_id1)

        # Output capacitances
        CLtot = alpha/beta/s.vod**2*sp.constants.Boltzmann*s.T
        # s.fan_out = FO = CL/CS   G = CS/CF   FO.G = CL/CF
        CF = CLtot/(1+d.rself)/(s.fan_out*s.G+1-beta)
        CS = s.G * CF
        # CL = s.fan_out * CS

        # Current division factor kappa
        gm_gds1 = NCH.look_up('GM_GDS', gm_id=d.gm_id1, vsb=0, L=d.L1)
        kappa = 1/(1+gm1_gm3/gm_gds1+2/gm_gds2)

        # Transconductance
        gm1 = 2*np.pi*s.fu1*CLtot/beta/kappa

        # Current
        id1 = gm1/d.gm_id1

        # Transit frequency
        wTi = NCH.look_up('GM_CGG', gm_id=d.gm_id1, vsb=0, L=d.L1)

        # Input capacitances
        Cgg1 = gm1/wTi #= Cgs1
        # Cgd1 = NCH.look_up('CGD_W', gm_id=d.gm_id1, vsb=0, L=d.L1)
        jd1 = NCH.look_up('ID_W', gm_id=d.gm_id1, vsb=0, L=d.L1)
        W = id1/jd1
        Cgd1 = W*NCH.look_up('CGD_W', gm_id=d.gm_id1, vsb=0, L=d.L1)
        Cin = Cgg1 + Cgd1*gm1_gm3

        # Actual beta & physical design point
        beta_actual = CF/(CF+CS+Cin)
        gm_id1_physical = interp1(beta_actual, d.gm_id1, kind='nearest', fill_value="extrapolate")(beta)
        physical_index = np.where(d.gm_id1==gm_id1_physical)


        # Save
        for var in ['gm_id1', 'id1', 'wTi', 'Cgg1', 'Cgd1', 'kappa', 'CLtot', 'CF', 'CS', 'W', 'alpha']:
            if var not in m1.keys():
                m1[var] = np.array([])
            m1[var] = np.append(m1[var], locals()[var][physical_index])

        pass
        # m1.gm_id1.append(gm_id1_physical)
        # m1.id1.append(id1[physical_index])
        # m1.wTi.append(wTi[physical_index])
        # m1.Cgg1.append(Cgg1[physical_index])
        # m1.Cgd1.append(Cgd1[physical_index])
        # m1.kappa.append(kappa[physical_index])
        # m1.CLtot.append(CLtot[physical_index])
        # m1.CF.append(CF[physical_index])
        # m1.CS.append(CS[physical_index])
        # m1.W.append(W[physical_index])

        # Check hypothesis
        # wrong formula ? gamma = NCH.look_up('STH_GM', gm_id=gm_id1_physical, L=d.L1)/4/sp.constants.Boltzmann/s.T  #  vgs=vgs1
        # m1.gamma.append(gamma)
        # plt.plot(np.array(m1.gm_id1), np.array(m1.gm_id1).reshape(-1, 1) * np.array(m1.id1))

    plt.figure()
    plt.plot(m1['gm_id1'], m1['id1'], label='L = %.2f um'%d.L1)
    plt.ylabel(r"$I_D$ [mA]")
    plt.xlabel(r"$g_m/I_D$ [S/A]")
    plt.title(r'$I_D$ vs. $g_m/I_D$ for varying $L$')
    plt.legend() #np.around(d.L1, decimals=2))
    plt.show()

    return m1, p1


L1 = [0.5, 0.6, 0.7, 0.8]
class d: pass # structure to collect design parameters
d.rself = 0
d.gm_id1 = np.linspace(3,25, 100)
d.beta = beta_max*np.linspace(0.2,1,100)

m1 = np.zeros(len(L1))
p1 = np.zeros(len(L1))
for i,L in enumerate(L1):
    d.L1 = L
    m1[i], p1[i] = folded_cascode(s, d)

    # Current minimun, new rself
    # id1 = min(m1.id1)
    # gm_id1

# jd1 = NCH.look_up('ID_W', gm_id=d.gm_id1, vsb=0, L=d.L1)

pass




