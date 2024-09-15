"""
    Example of the usage of the module for a 2 mm parallel-plate 
    ionization chamber exposed to a 0.1415 Gy pulse of 2.5 us duration.
"""

import PPICSimulation as PPICS
import matplotlib.pylab as plt
import time
import numpy as np

# dpp           : Dose per pulse in Gy.
# alpha         : Volume recombination parameter in m^{3} s^{-1}. Usually we 
#                 are using a value between 0.9E-12 m^{3} s^{-1} up to 1.5E-12 m^{3} s^{-1}.
# pulseDuration : Pulse duration of the beam in s.
# temperature   : Temperature of the air in the ionization chamber in degree Celsius.
# pressure      : Pressure of the air inside the ionization chamber in hPa.
# rHumidity     : Relative humidity of the air in the ionization chamber in %.
# voltage       : Bias applied voltage (always positive) in V.
# gap           : Distance between electrodes of the ionization chamber in m.
# Ndw           : Calibration coefficient of the ionization chamber in Gy C^{-1}.
#                 The calibration coefficient must have applied all the factor
#                 related to the charge released in the medium but not the
#                 temperature and pressure correction which the simulation will
#                 take into account.
# area          : Area of the sensitive volume of the ionization chamber in m^{2}.
# n             : Number of discretization steps in position. A reasonable number
#                 is around 1000. You may increase to have lower numerical error.

dpp            = 0.1415  # Gy
alpha          = 1.1E-12 # m^{3}s^{-1}
pulseDuration  = 2.5E-6  # s
temperature    = 18.77   # ºC
pressure       = 1013.0  # hPa
rHumidity      = 50      # %
voltage        = 100     # V
gap            = 2E-3    # m
Ndw            = 8.18E7  # Gy C^{-1}
kQ             = 0.8954
area           = 2.01E-4 # m^{2}
n              = 2000

inputParameters = [dpp, pulseDuration, alpha, voltage, temperature, pressure,
                   rHumidity, gap, n, Ndw * kQ, area]

t0 = time.time()
CCE, FEF0, FEF1, Q_coll, intensity = PPICS.pulsedSimulation(*inputParameters)

dt = time.time() - t0

# CCE    : Charge collection efficiency
# FEF0   : Free electron fraction in relation to the released charge in the 
#          medium.
# FEF1   : Free electron fraction in relation to the collected charge.
# Q_coll : Collected charge per pulse referenced to 20 ºC and 1013.25 hPa.

print(f"CCE    = {CCE:.4f}")
print(f"FEF0   = {FEF0:.4f}")
print(f"FEF1   = {FEF1:.4f}")
print(f"Q_coll = {Q_coll * 1E9:.4f} nC\n")

print(f"Elapsed time = {dt:.4f} s")


## Plot of the instantaneous current from the simulation.
intensity = np.array(intensity)

fig, ax = plt.subplots(figsize = (8, 6))

ax.plot(intensity[:, 0] * 1E6, intensity[:, 1] * 1E3, "-r", label = "Electrons")
ax.plot(intensity[:, 0] * 1E6, intensity[:, 2] * 1E3, "-b", label = "Positive ions")
ax.plot(intensity[:, 0] * 1E6, intensity[:, 3] * 1E3, "-m", label = "Negative ions")

ax.set_xlabel(r"Time (us)")
ax.set_ylabel(r"Intensity (mA)")

ax.set_xscale("log")
ax.set_yscale("log")

ax.set_ylim([1E-4, 1E-1])
ax.set_xlim([1E-6, 1E3])

fig.tight_layout()

plt.show()

# This script with the above parameters print:
# CCE    = 0.4020
# FEF0   = 0.1122
# FEF1   = 0.2791
# Q_coll = 0.7766 nC
