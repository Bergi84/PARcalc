import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline

def calcUmolPerLm(spdFileName):
    lumPerWatt = 683

    # read Luminous efficacy spectrum
    # ISO/CIE 11664:2015
    lumEffX, lumEffY = np.split(np.transpose(np.genfromtxt("luminousEfficiency.csv")), 2)
    lumEffX = np.ravel(lumEffX)
    lumEffY = np.ravel(lumEffY)

    lumEff_Sp = make_interp_spline(lumEffX, lumEffY, k=3)

    # read relative spectral power density
    spd_X, spd_Y = np.split(np.transpose(np.genfromtxt(spdFileName)), 2)
    spd_X = np.ravel(spd_X)
    spd_Y = np.ravel(spd_Y)

    spd_Sp = make_interp_spline(spd_X, spd_Y, k=3)

    # interpolate both cures to fixed X values
    X = np.linspace(np.min(spd_X), np.max(spd_X), 1000)
    lumEff_Y = lumEff_Sp(X)
    spd_Y = spd_Sp(X)

    # normalize relative spectral power density
    intSPD = np.trapz(spd_Y, X)
    spd_Y_norm = spd_Y/intSPD

    # multiply spd with luminous efficacy and lumen per Watt
    # to get the power (photonic power, not electrical) to lumen efficiency for this lamp
    spd_lm_Y = np.multiply(spd_Y_norm, lumEff_Y * lumPerWatt)   # lm/W
    lmPerW = np.trapz(spd_lm_Y, X)

    c = 2.998*10**8      # speed of light m/s
    h = 6.626*10**-34    # planck's constant J*s
    N_A = 6.022*10**23   # avogadro's number 1/mol

    spd_E_Y = h*c/(X*10**-9)*N_A        # J/mol, X is nm so we need to multiply with 10^-9
    spd_molPerW = np.multiply(spd_E_Y, spd_lm_Y)    # lm*s/mol

    spd_molPerW_Sp = make_interp_spline(X, spd_molPerW, k=3)

    X_par = np.linspace(400, 700, 1000)
    y_par = spd_molPerW_Sp(X_par)

    molPerLm = 1/np.trapz(y_par, X_par)

    return molPerLm


k_CREE_J2835_4000K_90CRI = calcUmolPerLm("CREE_J2835_4000K 90CRI.csv")

print("%.4f" % (k_CREE_J2835_4000K_90CRI*10**6), "µmol/(Lm*s)\n")
