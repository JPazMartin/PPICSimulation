cpdef double mobPos(double H, double P, double T)
cpdef double mobNeg(double H, double P, double T)

cdef double linearInterpolation(double[::1] x_theo, double[::1] y_theo, double x_ex, unsigned int n)