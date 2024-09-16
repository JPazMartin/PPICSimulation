"""
    Python extension to simulate the charge collection efficiency for parallel
    plate ionization chambers. 
"""

import numpy as np
import os

from utils cimport mobPos, mobNeg, linearInterpolation
from libc.math cimport exp

def pulsedSimulation(double dosepp, double pulseDuration, double alpha,
                     double voltage, double temperature, double pressure,
                     double humidity, double gap, unsigned int n, double Ndw,
                     double area, double[::1] timeStruct = np.array([]), 
                     double[::1] doseStruct = np.array([])):
    
    """
    This function calculates the charge collection efficiency of a parallel
    plate ionization chamber exposed to a given dose per pulse in pulsed beams.
    
    *-- Inputs:
        
        :dosepp:         Dose per pulse in Gy.
        :pulseDuration:  Pulse duration in s.
        :alpha_ref:      Volume recombination constant in m^{3} s^{-1}.
        :voltage:        Voltage in V.
        :temperature:    Temperature in degree Celsius.
        :pressure:       Pressure in hPa.
        :humidity:       Relative humidity in %.
        :gap:            Ionization chamber gap in m.
        :Ndw:            Calibration coefficient in Gy C^{-1}. The calibration
                         coefficient must include all the correction factors
                         that affect the released charge excluding k_TP.
        :area:           Area of the sensitive volume in m^{2}.
        :timeStruct:     Optional, array with the time in s for a custom pulse
                         structure.
        :DoseStruct:     Optional, array with the dose rate in arbitrary units
                         for a custom pulse structure.
    
    *-- Returns:
        
        :CCE:       Charge collection efficiency.
        :FEF0:      Free electron fraction in relation to the released charge.
        :FEF1:      Free electron fraction in relation to the collected charge.
        :Q_coll:    Collected charge at 20 degree Celsius and 1013.25 hPa (STP).
        :Intensity: Array with the time in s and intensity in A for each of the
                    three species considered.

                    . Intensity[:, 0] -> Array with the times in s.
                    . Intensity[:, 1] -> Array with the instantaneous current 
                                         from electrons in A.
                    . Intensity[:, 2] -> Array with the instantaneous current
                                         from positive ions in A.
                    . Intensity[:, 3] -> Array with the instantaneous current 
                                         from negative ions in A.
    """
    
    cdef:
        
        ## Basic physical constants.
        double e  = 1.6021766208E-19 # C
        double e0 = 8.854187817E-12  # F m^{-1}
        double er = 1.000589
        
        ## Reference temperature and pressure
        double refTemperature = 293.15  # K
        double refPressure    = 1013.25 # hPa

        ## Negative and positive ion mobilities
        double k_neg = mobNeg(humidity, pressure, temperature)
        double k_pos = mobPos(humidity, pressure, temperature)

        double dx = gap / n

        ## -- Load electron properties:
        ## Path where the module is placed
        str c_dir = os.path.dirname(__file__)
        ## Paths of the electron velocity and ionization rate
        str p_data = f"{c_dir}/Data/dataElectrons.txt"
        ## -- Load Electron properties --
        double[:, :] dataElectrons = np.loadtxt(p_data)
        double[::1] eFieldTable    = np.ascontiguousarray(dataElectrons[:, 0])
        double[::1] eVelocityTable = np.ascontiguousarray(dataElectrons[:, 1])
        double[::1] eIRateTable    = np.ascontiguousarray(dataElectrons[:, 2])
        
        unsigned int nData = len(eFieldTable)

        ## Temperature and pressure correction
        double k_TP = (273.15 + temperature) / refTemperature * refPressure / pressure
        
        ## Variables needed for the simulation
        double time, tStep, chargePos, chargeNeg, chargeE, tMax, n0, recomb
        double attachment, multiplication, eFieldSum, nSum, eMax, eSum, vMax 
        double rel, Ie, Ineg, Ipos

        ## List to hold the intensities of each carrier specie.
        list Intensity = []
        
        ## Arrays for electron parameters
        double[::1] eVelocity   = np.empty(n)
        double[::1] eAttachment = np.empty(n)
        double[::1] eIRate      = np.empty(n)

        ## Electric field
        double[::1] eField = np.ones(n) * (abs(voltage) / gap)

        ## Charge densities.
        double[::1] nE    = np.zeros(n)
        double[::1] nNeg  = np.zeros(n)
        double[::1] nPos  = np.zeros(n)

        ## To avoid using a linear interpolation in each step the parameters
        ## are interpolated into small steps and the parameters are considered
        ## to be constant along this steps. This improves the speed fo the code.
        int nBin = 100000
        double[::1] eFieldBin     = np.zeros(nBin)
        double[::1] eVelocityBin  = np.zeros(nBin)
        double[::1] eIRateBin     = np.zeros(nBin)

        ## Other needed parameters
        double[::1] totalCharge = np.empty(n)

        double stepE = (eFieldTable[nData - 1] - eFieldTable[0]) / nBin

        ## If 0 a constant charge release is used during the pulse duration. If 
        ## 1 use the custom dosepp rate.
        int structure = 0
    
        unsigned int i, idxMaxE, index, idxBin, nStruct = 0

        double[::1] n0Struct = np.empty(0)


    n0 = dosepp / (e * Ndw * area * gap * pulseDuration * k_TP)

    ## Compute the charge released per unit of volume and time:
    if len(timeStruct) != 0:

        structure = 1

        nStruct  = len(timeStruct)
        n0Struct = np.empty(nStruct)

        ## Normalize the pulse structure just in case:
        doseStruct *= (dosepp / np.trapz(doseStruct, timeStruct))

        for i in range(nStruct):
            n0Struct[i] = doseStruct[i] / (e * Ndw * area * gap * k_TP)

        pulseDuration  = timeStruct[nStruct - 1]

    ## Take into account the temperature and pressure for the electron
    ## transport properties.
    for i in range(nData):
        
        eFieldTable[i] = eFieldTable[i] / k_TP
        eIRateTable[i] = eIRateTable[i] / k_TP

    for i in range(1, nBin + 1):

        eFieldBin[i] = stepE * i

    for i in range(nBin):
        eVelocityBin[i] = linearInterpolation(eFieldTable, eVelocityTable,
                                              eFieldBin[i], nData)
        eIRateBin[i]    = linearInterpolation(eFieldTable, eIRateTable,
                                              eFieldBin[i], nData)

    ## -- Initialize variables --
    nSum = 0; index = 0; time = 0; chargePos = 0; chargeNeg = 0; chargeE = 0

    ## Dafault initial time step.
    tStep = 1E-13

    ## Maximum time step. Just to avoid problems when the pulse duration becomes
    ## really short in comparison with the computed time step.
    tMax  = pulseDuration / 2000
    if structure: tMax = np.min(np.diff(timeStruct)) / 2

    ## -- Simulation starts --
    while (nSum > 1) or (time < pulseDuration):
        
        ## Update variables
        rel     = tStep / dx
        time   += tStep
        index  += 1

        ## -- Zeroth step: Calculate or update transport properties --        
        for i in range(n):
            
            idxBin = int((eField[i] - eFieldBin[0]) / stepE + 0.5)

            if idxBin >= nBin:
                print("Error: Electric field is higher than the provided from tables.")
                return 0, 0, 0, 0, np.array([0, 0, 0, 0])
            
            ## Electron attachment from Boissonnat et al. (arXiv:1609.03740v1)
            eAttachment[i] = 95.24E-9 * (1 - exp(-eField[i] / 258.5E3))
            eVelocity[i]   = eVelocityBin[idxBin]
            eIRate[i]      = eIRateBin[idxBin]
        
        ## -- First step: Irradiation --
        if time <= pulseDuration:

            if structure:
                n0 = linearInterpolation(timeStruct, n0Struct, time, nStruct)

            for i in range(n):
         
                nE[i]   += n0 * tStep
                nPos[i] += n0 * tStep
         
        ### -- Second step: Ion-ion Recombination --
        for i in range(n):
    
            recomb = alpha * nNeg[i] * nPos[i] * tStep
        
            nNeg[i] -= recomb
            nPos[i] -= recomb
        
          
        # -- Third step: Set to zero --
        for i in range(n):
    
            if nNeg[i] < 1E-30: nNeg[i] = 0.0
            if nPos[i] < 1E-30: nPos[i] = 0.0
            if nE[i] < 1E-30:   nE[i]   = 0.0
            
        ## -- Forth step: Electric field perturbation --
        eFieldSum = 0
        totalCharge[0] = nPos[0] - nE[0] - nNeg[0]

        for i in range(1, n):
        
            totalCharge[i] = totalCharge[i - 1] + nPos[i] - nE[i] - nNeg[i]
        
        for i in range(n):
            
            eField[i] = (e / (2 * er * e0)) * (gap / n) * (2 * totalCharge[i]\
                        - totalCharge[n - 1])
            eFieldSum += eField[i]
            
        for i in range(n):
            
            eField[i] += abs(voltage) / gap - eFieldSum / n
        
        ## -- Fifth step: Propagation --
        
        # *-- Positive carriers:
        chargePos += nPos[n - 1] * k_pos * eField[n - 1] * rel * e * dx

        for i in range(1, n):
    
            nPos[n - i] += (nPos[n - i - 1] * eField[n - i - 1] - nPos[n - i] \
                             * eField[n - i]) * rel * k_pos
    
        nPos[0] -= nPos[0] * k_pos * eField[0] * rel      
        
        # *-- Negative carriers:
        chargeNeg += nNeg[0] * k_neg * eField[0] * rel * e * dx

        for i in range(n - 1):
                
            nNeg[i] += (nNeg[i + 1] * eField[i + 1] - nNeg[i] * eField[i]) * rel * k_neg
            
        nNeg[n - 1] -= nNeg[n - 1] * k_neg * eField[n - 1] * rel
        
        # *-- Electrons:
        chargeE += nE[0] * eVelocity[0] * rel * e * dx
    
        for i in range(n - 1):
                
            nE[i] += (nE[i + 1] * eVelocity[i + 1] - nE[i] * eVelocity[i]) * rel
            
        nE[n - 1] += - nE[n - 1] * eVelocity[n - 1] * rel
        
        ## -- Sixth step: Attachment --
        for i in range(n):
    
            attachment = nE[i] * tStep / eAttachment[i]
            nNeg[i] += attachment
            nE[i]   -= attachment

        # -- Seven step: Multiplication --
        for i in range(n):
            
            multiplication = nE[i] * eIRate[i] * eVelocity[i] * tStep
            
            nE[i]   += multiplication
            nPos[i] += multiplication
        
        
        ## -- Eighth step: Step adjustment --
        eMax = eField[0]; eSum = 0; idxMaxE = 0
        
        for i in range(n):
            
            eSum += nE[i]
            
            if (eField[i] > eMax): 
                eMax     = eField[i] 
                idxMaxE = i
        
        if eVelocity[idxMaxE] > (k_neg * eField[idxMaxE]) and eSum > 1:
            
            vMax = eVelocity[idxMaxE]
            
        else:
            
            vMax = k_neg * eField[idxMaxE]
            
        tStep = 0.4 * (gap / n) / vMax
        
        if (time < pulseDuration) and (tStep > tMax):
            
            tStep = tMax
        
        ## -- Ninth step: Carriers in the volume --
        nSum = 0; Ie = 0; Ipos = 0; Ineg = 0
        for i in range(n):
            
            nSum += nPos[i] + nNeg[i] + nE[i]
        
            Ie   += eVelocity[i] * nE[i] * e * area / n
            Ipos += k_pos * eField[i] * nPos[i] * e * area / n
            Ineg += k_neg * eField[i] * nNeg[i] * e * area / n

        Intensity.append([time, Ie, Ipos, Ineg])

    CCE   = chargePos / (n0 * pulseDuration * e * gap) # CCE
    FEF0  = chargeE / (n0 * pulseDuration * e * gap)   # FEF (Boag definition 0)
    FEF1  = chargeE / chargePos                        # FEF (Boag definition 1)
    Qcoll = chargePos * k_TP * area                    # Charge referenced to STP
    
    return CCE, FEF0, FEF1, Qcoll, Intensity
