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

def load_parameters_65nm():
    #%% Load process data into pygmid lookup object
    class d: pass # structure to collect design parameters
    class s: pass  # structure to collect specs

    NCH = lk('65nch.mat')
    PCH = lk('65pch.mat')

    # Cascode gm/id defined by output swing
    common_mode_output_voltage = 0.6
    output_swing = 0.8
    Vds_min = (common_mode_output_voltage - output_swing/4)/2
    gm_id_cas_max = 2/Vds_min
    d.gm_id_cas = 15
    Lsweep = np.linspace(0.06,1,100)

    # Specifications
    s.G = 2
    s.ts = 5e-9
    s.ed = 0.1e-2
    s.vod = 400e-6  # total integrated noise
    s.T = 300  # Kelvin
    s.fan_out = 0.5

    # Design variable
    d.L1 = [0.1, 0.2, 0.3, 0.4]
    d.Lcas = 0.4  # L0 ~ 60
    # d.rself = 0
    d.cself = 0
    d.gm_id1 = np.linspace(3,27, 100)
    d.gamma = 0.7

    # NCH = lk('xt018_ne5.pkl')
    # PCH = lk('xt018_pe5.pkl')
    # d.Lcas = 2  # 300 < L0 < 400  (50 with book techno in 65 nm)
    # s.ts = 100e-9
    # d.L1 = [0.5, 0.6, 0.7, 0.8]
    # d.gm_id1 = np.linspace(3,25, 100)
    # Lsweep = np.linspace(0.5,5,100)

    return NCH, PCH, s, d, Lsweep

# This custom formatter removes trailing zeros, e.g. "1.0" becomes "1", and
def fmt(x):
    s = f"{x:.1f}"
    if s.endswith("0"):
        s = f"{x:.0f}"
    # return rf"{s} \%" if plt.rcParams["text.usetex"] else f"{s} %"
    return rf"{s}"

def optimise_folded_cascode():
    # %% Gm over Id methodology applied to folded-cascode OTA

    NCH, PCH, s, d, Lsweep = load_parameters_65nm()

    beta_max = 1/(1+s.G)
    beta = 0.75*beta_max # first-order optimum
    kappa = 0.7 # conservative optimum

    # Channel length sweep
    gm_gds2 = NCH.look_up('GM_GDS', gm_id=d.gm_id_cas, vds=0.2, L=Lsweep)
    gm_gds3 = NCH.look_up('GM_GDS', gm_id=d.gm_id_cas, vds=0.4, L=Lsweep)
    gm_gds4 = PCH.look_up('GM_GDS', gm_id=d.gm_id_cas, vds=0.4, L=Lsweep).reshape(-1,1)
    gm_gds5 = PCH.look_up('GM_GDS', gm_id=d.gm_id_cas, vds=0.2, L=Lsweep)[:, np.newaxis]
    # alternatives for transforming a (horizontal) 1D vector to a vertical vector
    # requires to add a dimension

    # Chosen length
    L0_swept = beta*kappa/(1/(1+gm_gds5)/gm_gds4+1/(1+gm_gds2/3)/gm_gds3)


    # fig, ax = plt.subplots()
    # levels= np.linspace(10,100,10)
    # CS = ax.contour(L23, L45, L0_swept) #, levels)
    # ax.set_ylabel(r"$L_{4,5}$ [$\mu$m]")
    # ax.set_xlabel(r"$L_{2,3}$ [$\mu$m]")
    # ax.set_title(r'Low frequency loop gain $L_0$ as a function of channel lengths')
    # ax.clabel(CS, CS.levels, fmt=fmt, fontsize=10)
    # plt.show()
    # plt.close()


    gmb_gm3 = NCH.look_up('GMB_GM', gm_id=d.gm_id_cas, vds=0.4, vsb=0.2, L=d.Lcas)
    gm_css3 = NCH.look_up('GM_CSS', gm_id=d.gm_id_cas, vds=0.4, vsb=0.2, L=d.Lcas)
    cdd_css3 = NCH.look_up('CDD_CSS', gm_id=d.gm_id_cas, vds=0.4, vsb=0.2, L=d.Lcas)
    cdd_w3 = NCH.look_up('CDD_CSS', gm_id=d.gm_id_cas, vds=0.4, vsb=0.2, L=d.Lcas)
    cdd_w2 = NCH.look_up('CDD_CSS', gm_id=d.gm_id_cas, vds=0.2, L=d.Lcas) # vsb=0.2,
    fp2 = 1/2/np.pi*gm_css3*(1+gmb_gm3)/(1+2*cdd_css3*2*(cdd_w2/cdd_w3))
    print(f'Non dominant pole fp2 = {Float(fp2):!.2h}Hz') #%.2E Hz' % fp2)
    fp2 = 2e9

    # %% Design for a given settling time & noise budget, minimise power dissipation
    s.fu1 = 1/2/np.pi*np.log(1/s.ed)/s.ts
    print("fp2/s.fu1 = %.2f gives a phase margin of about "%(fp2/s.fu1))

    d.beta = beta_max*np.linspace(0.2,1,100)
    plt.figure()
    for i,L in enumerate(d.L1):
        d.L1 = L
        r = folded_cascode(NCH, PCH, s, d)

        # plt.plot(r.gm_id1, r.id1*1e3, label='L = %.2f um'%d.L1)
        plt.plot(r.gm_id1, r.rself, label='L = %.2f um'%d.L1)

    plt.ylabel(r"$I_D$ [mA]")
    plt.xlabel(r"$g_m/I_D$ [S/A]")
    plt.title(r'$I_D$ vs. $g_m/I_D$ for varying $L$')
    # plt.xlim([5, 27])
    # plt.ylim([0, 1.2])
    plt.legend() #np.around(d.L1, decimals=2))
    plt.show()

    # jd1 = NCH.look_up('ID_W', gm_id=d.gm_id1, vsb=0, L=d.L1)

    pass

def folded_cascode(NCH, PCH, s, d):

    # Extra hypothesis
    gamma1 = gamma5 = gamma2 = d.gamma
    gm_id5 = gm_id2 = d.gm_id_cas

    # Precomputation
    gm_id1 = d.gm_id1  # simplify saving singular optimum variable among input array
    gm1_gm3 = d.gm_id1 / d.gm_id_cas
    gm_gds1 = NCH.look_up('GM_GDS', gm_id=d.gm_id1, L=d.L1)
    gm_gds2 = NCH.look_up('GM_GDS', gm_id=d.gm_id_cas, L=d.Lcas)
    wT1 = NCH.look_up('GM_CGG', gm_id=d.gm_id1, L=d.L1)  # Transit frequency
    cgd_cgg1 = NCH.look_up('CGD_CGG', gm_id=d.gm_id1, L=d.L1)
    cdd_gm3 = NCH.look_up('CDD_GM', gm_id=d.gm_id_cas, L=d.Lcas)
    cdd_gm4 = PCH.look_up('CDD_GM', gm_id=d.gm_id_cas, L=d.Lcas)

    # cgd_w1 = NCH.look_up('CGD_W', gm_id=d.gm_id1, vsb=0, L=d.L1)
    # jd1 = NCH.look_up('ID_W', gm_id=d.gm_id1, vsb=0, L=d.L1)
    # cgb_cgg1 = NCH.look_up('CGB_CGG', gm_id=d.gm_id1, L=d.L1)

    # store results
    class r: pass

    # 2D sweep for finding beta with a physical reality
    for beta in d.beta:
        # Excess noise factor
        alpha = 2*gamma1*(1+gamma5/gamma1*gm_id5/d.gm_id1+2*gamma2/gamma1*gm_id2/d.gm_id1)

        # Output capacitances
        CLtot = alpha/beta/s.vod**2*sp.constants.Boltzmann*s.T
        # s.fan_out = FO = CL/CS   G = CS/CF   FO.G = CL/CF
        # CF = CLtot/(1+d.rself)/(s.fan_out*s.G+1-beta)
        CF = (CLtot-d.cself)/(s.fan_out*s.G+1-beta)
        CS = s.G * CF
        # CL = s.fan_out * CS

        # Current division factor kappa
        kappa = 1/(1+gm1_gm3/gm_gds1+2/gm_gds2)

        # Transconductance
        gm1 = 2*np.pi*s.fu1*CLtot/beta/kappa

        # Current
        id1 = gm1/d.gm_id1

        # Input capacitances
        Cgs1 = gm1/wT1
        Cgd1 = Cgs1 * cgd_cgg1
        # Alternative to compute Cgd
        # W = id1 / jd1
        # Cgd1 = W * cgd_w1
        Cin = Cgs1 + Cgd1*gm1_gm3 # cgb_cgg1 is neglectable

        # Actual beta & physical design point
        beta_actual = CF/(CF+CS+Cin)
        gm_id1_physical = interp1(beta_actual, d.gm_id1, kind='nearest', fill_value="extrapolate")(beta)
        physical_index = np.where(d.gm_id1==gm_id1_physical)

        # Save
        for var in ['id1', 'gm_id1', 'CLtot']: #['gm_id1', 'id1', 'wTi', 'Cgg1', 'Cgd1', 'kappa', 'CLtot', 'CF', 'CS', 'W', 'alpha']:
            if var not in r.__dict__.keys():
                setattr(r, var, np.array([]))
            setattr(r, var, np.append(r.__dict__[var], locals()[var][physical_index]))
        pass

    # Check hypothesis
    # wrong formula ? gamma = NCH.look_up('STH_GM', gm_id=gm_id1_physical, L=d.L1)/4/sp.constants.Boltzmann/s.T  #  vgs=vgs1
    # m1.gamma.append(gamma)
    # plt.plot(np.array(m1.gm_id1), np.array(m1.gm_id1).reshape(-1, 1) * np.array(m1.id1))

    gm34 = r.id1/d.gm_id_cas
    r.cself = gm34*(cdd_gm3+cdd_gm4)
    r.rself = r.cself/r.CLtot

    return r



if __name__=='__main__':
    optimise_folded_cascode()
