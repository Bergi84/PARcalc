import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline

lumPerWatt = 683

# read Luminous efficacy spectrum
# ISO/CIE 11664:2015
lumEffX, lumEffY = np.split(np.transpose(np.genfromtxt("luminousEfficiency.csv")), 2)
lumEffX = np.ravel(lumEffX)
lumEffY = np.ravel(lumEffY)

lumEff_Sp = make_interp_spline(lumEffX, lumEffY, k=3)

# read relative spectral power density
spd_X, spd_Y = np.split(np.transpose(np.genfromtxt("CREE_J2835_4000K 90CRI.csv")), 2)
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
spd_lm_Y = np.multiply(spd_Y_norm, lumEff_Y * lumPerWatt)
lmPerW = np.trapz(spd_lm_Y, X)

print("%.2f" % lmPerW, " lm/W")

plt.subplot(3,1,1)
plt.plot(X, lumEff_Y)

plt.subplot(3,1,2)
plt.plot(X, spd_Y)

plt.subplot(3,1,3)
plt.plot(X, spd_lm_Y)

plt.show()
