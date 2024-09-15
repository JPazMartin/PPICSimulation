"""
    Functions for calculating the positive and negative ions mobilities in air 
    from:
    Zhang et al. Prediction of Average Mobility of Ions from Corona Discharge in
    Air with Respect to Pressure, Humidity and Temperature" IEEE Transactions on
    Dielectrics and Electrical Insulation Vol. 26, No. 5; October 2019
    (10.1109/TDEI.2019.008001)
"""

from libc.math cimport exp, sqrt # C-functions

cpdef double relativeToAbsoluteHumidity(double temperature, double rHumidity):
    
    cdef double temperature2 = temperature + 273.15
    cdef double value
    
    if temperature >= 0:
        
        value = 10**((10.286 * temperature2 - 2148.4909) / (temperature2 - 35.85)) / 1000
        
    else :
        
        value = 10**(12.5633 - 2670.59 / 273.1) / 1000
        
    return rHumidity * value / 100

cpdef double alphaPos(double temperature, double aHumidity):
    
    cdef:

        double temperature2 = temperature + 273.15
        double a1 = 1.438E-4 * temperature2**(0.451)
        double a2 = 1.964E-4 * temperature2**(0.513)
        double a3 = -0.95053 - exp((292.3191 - temperature2) / 12.49447)
        double h0 = 0.05285 + exp((temperature2 - 301.5777) / 18.99722)
    
    return a1 + (a2 - a1) / (1 + 10**(a3 * (h0 - aHumidity)))


cpdef double betaPos(double temperature, double aHumidity):
    
    cdef:
        double temperature2 = temperature + 273.15
        double b1 = -0.55579 - exp((temperature2 - 341.1570) / 21.32626)
        double b2 = 7.84584 - 0.06263 * temperature2 + 1.16655E-4 * temperature**2
        double b3 = 0.0377 + 0.144 / (1 + 10**((278 - temperature2) / 16.8))
    
    
    return b1 + b2 * exp(-aHumidity / b3)


cpdef double alphaNeg(double temperature, double aHumidity):
    
    cdef:

        double temperature2 = temperature + 273.15
        double a1 = 2.476E-4 * temperature2**(0.377)
        double a2 = 2.359E-4 * temperature2**(0.499)
        double a3 = -0.79778 - exp((283.4930 - temperature2) / 9.22739)
        double h0 = -0.04443 + exp((temperature2 - 301.1779) / 27.45497)
    
    return a1 + (a2 - a1) / (1 + 10**(a3 * (h0 - aHumidity)))


cpdef double betaNeg(double temperature, double aHumidity):
    
    
    cdef:
        double temperature2 = temperature + 273.15
        double b1 = -0.57882 - exp((temperature2 - 370.5096) / 46.72047)
        double b2 = 2.31325 - 0.02408 * temperature2 + 4.93039E-5 * temperature2**2
        double b3 = 0.00152 + 0.112 / (1 + 10**((269 - temperature2) / 21.3))
    
    return b1 + b2 * exp(-aHumidity / b3)

cpdef double mobPos(double rHumidity, double pressure, double temperature):
    
    cdef:
    
        double aHumidity = relativeToAbsoluteHumidity(temperature, rHumidity)
        
    pressure = pressure / 10
    
    return 1 / sqrt(temperature + 273.15) * alphaPos(temperature, aHumidity) * (3 * pressure / (temperature + 273.15))**betaPos(temperature, aHumidity)
    
cpdef double mobNeg(double rHumidity, double pressure, double temperature):
    
    cdef:
        
        double aHumidity = relativeToAbsoluteHumidity(temperature, rHumidity)
        
    pressure = pressure / 10
    
    return 1 / sqrt(temperature + 273.15) * alphaNeg(temperature, aHumidity) * (3 * pressure / (temperature + 273.15))**betaNeg(temperature, aHumidity)


cdef double linearInterpolation(double[::1] x_theo, double[::1] y_theo, double x_ex, unsigned int n):
    
    
    cdef:
        
        unsigned int i
        double   y = y_theo[0]
        
    if x_ex > x_theo[n - 1]:
        
        return x_ex * (y_theo[n - 1] - y_theo[n - 2]) / (x_theo[n - 1] - x_theo[n - 2]) + (x_theo[n - 1] * y_theo[n - 2] - x_theo[n - 2] * y_theo[n - 1]) / (x_theo[n - 1] - x_theo[n - 2])
    
    for i in range(1, n):
        
        if x_ex < x_theo[i]:

            if x_ex > x_theo[i - 1]:

                y = y_theo[i - 1] + (x_ex - x_theo[i - 1]) / (x_theo[i] - x_theo[i - 1]) * (y_theo[i] - y_theo[i - 1])
                
                return y
            
        if x_ex == x_theo[i]:
            
            return y_theo[i]
            
    if x_ex < x_theo[0]:
        
        return y_theo[0]
        