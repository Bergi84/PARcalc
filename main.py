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


k_CREE_J2835_2700K_95CRI = calcUmolPerLm("CREE_J2835_2700K_95CRI.csv")
k_CREE_J2835_3000K_95CRI = calcUmolPerLm("CREE_J2835_3000K_95CRI.csv")
k_CREE_J2835_3500K_95CRI = calcUmolPerLm("CREE_J2835_3500K_95CRI.csv")
k_CREE_J2835_4000K_95CRI = calcUmolPerLm("CREE_J2835_4000K_95CRI.csv")
k_CREE_J2835_5000K_95CRI = calcUmolPerLm("CREE_J2835_5000K_95CRI.csv")
k_CREE_J2835_6500K_95CRI = calcUmolPerLm("CREE_J2835_6500K_95CRI.csv")

k_CREE_J2835_2200K_90CRI = calcUmolPerLm("CREE_J2835_2200K_90CRI.csv")
k_CREE_J2835_2700K_90CRI = calcUmolPerLm("CREE_J2835_2700K_90CRI.csv")
k_CREE_J2835_3000K_90CRI = calcUmolPerLm("CREE_J2835_3000K_90CRI.csv")
k_CREE_J2835_3500K_90CRI = calcUmolPerLm("CREE_J2835_3500K_90CRI.csv")
k_CREE_J2835_4000K_90CRI = calcUmolPerLm("CREE_J2835_4000K_90CRI.csv")
k_CREE_J2835_5000K_90CRI = calcUmolPerLm("CREE_J2835_5000K_90CRI.csv")
k_CREE_J2835_5700K_90CRI = calcUmolPerLm("CREE_J2835_5700K_90CRI.csv")
k_CREE_J2835_6500K_90CRI = calcUmolPerLm("CREE_J2835_6500K_90CRI.csv")

k_CREE_J2835_2200K_80CRI = calcUmolPerLm("CREE_J2835_2200K_80CRI.csv")
k_CREE_J2835_2700K_80CRI = calcUmolPerLm("CREE_J2835_2700K_80CRI.csv")
k_CREE_J2835_3000K_80CRI = calcUmolPerLm("CREE_J2835_3000K_80CRI.csv")
k_CREE_J2835_3500K_80CRI = calcUmolPerLm("CREE_J2835_3500K_80CRI.csv")
k_CREE_J2835_4000K_80CRI = calcUmolPerLm("CREE_J2835_4000K_80CRI.csv")
k_CREE_J2835_5000K_80CRI = calcUmolPerLm("CREE_J2835_5000K_80CRI.csv")
k_CREE_J2835_5700K_80CRI = calcUmolPerLm("CREE_J2835_5700K_80CRI.csv")
k_CREE_J2835_6500K_80CRI = calcUmolPerLm("CREE_J2835_6500K_80CRI.csv")

k_CREE_J2835_2700K_70CRI = calcUmolPerLm("CREE_J2835_2700K_70CRI.csv")
k_CREE_J2835_3000K_70CRI = calcUmolPerLm("CREE_J2835_3000K_70CRI.csv")
k_CREE_J2835_4000K_70CRI = calcUmolPerLm("CREE_J2835_4000K_70CRI.csv")
k_CREE_J2835_5000K_70CRI = calcUmolPerLm("CREE_J2835_5000K_70CRI.csv")
k_CREE_J2835_5700K_70CRI = calcUmolPerLm("CREE_J2835_5700K_70CRI.csv")
k_CREE_J2835_6500K_70CRI = calcUmolPerLm("CREE_J2835_6500K_70CRI.csv")

k_LUMILEDS_LUXEON_2835HE_1800K_90CRI = calcUmolPerLm("Lumileds_LUXEON_90CRI_2835HE_1800K.csv")
k_LUMILEDS_LUXEON_2835HE_2200K_90CRI = calcUmolPerLm("Lumileds_LUXEON_90CRI_2835HE_2200K.csv")
k_LUMILEDS_LUXEON_2835HE_2700K_90CRI = calcUmolPerLm("Lumileds_LUXEON_90CRI_2835HE_2700K.csv")
k_LUMILEDS_LUXEON_2835HE_3000K_90CRI = calcUmolPerLm("Lumileds_LUXEON_90CRI_2835HE_3000K.csv")
k_LUMILEDS_LUXEON_2835HE_3500K_90CRI = calcUmolPerLm("Lumileds_LUXEON_90CRI_2835HE_3500K.csv")
k_LUMILEDS_LUXEON_2835HE_4000K_90CRI = calcUmolPerLm("Lumileds_LUXEON_90CRI_2835HE_4000K.csv")
k_LUMILEDS_LUXEON_2835HE_5000K_90CRI = calcUmolPerLm("Lumileds_LUXEON_90CRI_2835HE_5000K.csv")
k_LUMILEDS_LUXEON_2835HE_5700K_90CRI = calcUmolPerLm("Lumileds_LUXEON_90CRI_2835HE_5700K.csv")
k_LUMILEDS_LUXEON_2835HE_6500K_90CRI = calcUmolPerLm("Lumileds_LUXEON_90CRI_2835HE_6500K.csv")


print("CREE 2700K, 95CRI: %.4f" % (k_CREE_J2835_2700K_95CRI*10**6), "µmol/(Lm*s)")
print("CREE 3000K, 95CRI: %.4f" % (k_CREE_J2835_3000K_95CRI*10**6), "µmol/(Lm*s)")
print("CREE 3500K, 95CRI: %.4f" % (k_CREE_J2835_3500K_95CRI*10**6), "µmol/(Lm*s)")
print("CREE 4000K, 95CRI: %.4f" % (k_CREE_J2835_4000K_95CRI*10**6), "µmol/(Lm*s)")
print("CREE 5000K, 95CRI: %.4f" % (k_CREE_J2835_5000K_95CRI*10**6), "µmol/(Lm*s)")
print("CREE 6500K, 95CRI: %.4f" % (k_CREE_J2835_6500K_95CRI*10**6), "µmol/(Lm*s)\n")

print("CREE 2200K, 90CRI: %.4f" % (k_CREE_J2835_2200K_90CRI*10**6), "µmol/(Lm*s)")
print("CREE 2700K, 90CRI: %.4f" % (k_CREE_J2835_2700K_90CRI*10**6), "µmol/(Lm*s)")
print("CREE 3000K, 90CRI: %.4f" % (k_CREE_J2835_3000K_90CRI*10**6), "µmol/(Lm*s)")
print("CREE 3500K, 90CRI: %.4f" % (k_CREE_J2835_3500K_90CRI*10**6), "µmol/(Lm*s)")
print("CREE 4000K, 90CRI: %.4f" % (k_CREE_J2835_4000K_90CRI*10**6), "µmol/(Lm*s)")
print("CREE 5000K, 90CRI: %.4f" % (k_CREE_J2835_5000K_90CRI*10**6), "µmol/(Lm*s)")
print("CREE 5700K, 90CRI: %.4f" % (k_CREE_J2835_5700K_90CRI*10**6), "µmol/(Lm*s)")
print("CREE 6500K, 90CRI: %.4f" % (k_CREE_J2835_6500K_90CRI*10**6), "µmol/(Lm*s)\n")

print("CREE 2200K, 80CRI: %.4f" % (k_CREE_J2835_2200K_80CRI*10**6), "µmol/(Lm*s)")
print("CREE 2700K, 80CRI: %.4f" % (k_CREE_J2835_2700K_80CRI*10**6), "µmol/(Lm*s)")
print("CREE 3000K, 80CRI: %.4f" % (k_CREE_J2835_3000K_80CRI*10**6), "µmol/(Lm*s)")
print("CREE 3500K, 80CRI: %.4f" % (k_CREE_J2835_3500K_80CRI*10**6), "µmol/(Lm*s)")
print("CREE 4000K, 80CRI: %.4f" % (k_CREE_J2835_4000K_80CRI*10**6), "µmol/(Lm*s)")
print("CREE 5000K, 80CRI: %.4f" % (k_CREE_J2835_5000K_80CRI*10**6), "µmol/(Lm*s)")
print("CREE 5700K, 80CRI: %.4f" % (k_CREE_J2835_5700K_80CRI*10**6), "µmol/(Lm*s)")
print("CREE 6500K, 80CRI: %.4f" % (k_CREE_J2835_6500K_80CRI*10**6), "µmol/(Lm*s)\n")

print("CREE 2700K, 70CRI: %.4f" % (k_CREE_J2835_2700K_70CRI*10**6), "µmol/(Lm*s)")
print("CREE 3000K, 70CRI: %.4f" % (k_CREE_J2835_3000K_70CRI*10**6), "µmol/(Lm*s)")
print("CREE 4000K, 70CRI: %.4f" % (k_CREE_J2835_4000K_70CRI*10**6), "µmol/(Lm*s)")
print("CREE 5000K, 70CRI: %.4f" % (k_CREE_J2835_5000K_70CRI*10**6), "µmol/(Lm*s)")
print("CREE 5700K, 70CRI: %.4f" % (k_CREE_J2835_5700K_70CRI*10**6), "µmol/(Lm*s)")
print("CREE 6500K, 70CRI: %.4f" % (k_CREE_J2835_6500K_70CRI*10**6), "µmol/(Lm*s)\n")

print("LUMILEDS 1800K, 90CRI: %.4f" % (k_LUMILEDS_LUXEON_2835HE_1800K_90CRI*10**6), "µmol/(Lm*s)")
print("LUMILEDS 2200K, 90CRI: %.4f" % (k_LUMILEDS_LUXEON_2835HE_2200K_90CRI*10**6), "µmol/(Lm*s)")
print("LUMILEDS 2700K, 90CRI: %.4f" % (k_LUMILEDS_LUXEON_2835HE_2700K_90CRI*10**6), "µmol/(Lm*s)")
print("LUMILEDS 3000K, 90CRI: %.4f" % (k_LUMILEDS_LUXEON_2835HE_3000K_90CRI*10**6), "µmol/(Lm*s)")
print("LUMILEDS 3500K, 90CRI: %.4f" % (k_LUMILEDS_LUXEON_2835HE_3500K_90CRI*10**6), "µmol/(Lm*s)")
print("LUMILEDS 4000K, 90CRI: %.4f" % (k_LUMILEDS_LUXEON_2835HE_4000K_90CRI*10**6), "µmol/(Lm*s)")
print("LUMILEDS 5000K, 90CRI: %.4f" % (k_LUMILEDS_LUXEON_2835HE_5000K_90CRI*10**6), "µmol/(Lm*s)")
print("LUMILEDS 5700K, 90CRI: %.4f" % (k_LUMILEDS_LUXEON_2835HE_5700K_90CRI*10**6), "µmol/(Lm*s)")
print("LUMILEDS 6500K, 90CRI: %.4f" % (k_LUMILEDS_LUXEON_2835HE_6500K_90CRI*10**6), "µmol/(Lm*s)\n")