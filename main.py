import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline

lumPerWatt = 683

# read Luminous efficacy spectrum
# ISO/CIE 11664:2015
lumEffX, lumEffY = np.split(np.transpose(np.genfromtxt("luminousEfficiency.csv")), 2)
lumEffX = np.ravel(lumEffX)
lumEffY = np.ravel(lumEffY)



# read relative spectral power density
spd_X, lm_spd_Y = np.split(np.transpose(np.genfromtxt("CREE_J2835_4000K 90CRI.csv")), 2)
spd_X = np.ravel(spd_X)
lm_spd_Y = np.ravel(lm_spd_Y)

lm_spd_Sp = make_interp_spline(spd_X, lm_spd_Y, k=3)
lm_spd_Y = lm_spd_Sp(lumEffX)

# normalize relative spectral power density
intSPD = np.trapz(lm_spd_Y, lumEffX)
lm_spd_Y_norm = lm_spd_Y/intSPD

# multiply spd with luminous efficacy and lumen per Watt
# to get the power (photonic power, not electrical) to lumen efficiency for this lamp
p_spd_Y = np.multiply(lm_spd_Y_norm, lumEffY * lumPerWatt)
lmPerW = np.trapz(p_spd_Y, lumEffX)

print("%.2f" % lmPerW, " lm/W")


plt.plot(lumEffX, p_spd_Y)
plt.show()
