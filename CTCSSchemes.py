# Numerical schemes for simulating linear advection for outer code
# linearAdvect.py 

# If you are using Python 2.7 rather than Python 3, import various
# functions from Python 3 such as to use real number division
# rather than integer division. ie 3/2  = 1.5  rather than 3/2 = 1
#from __future__ import absolute_import, division, print_function

# The numpy package for numerical functions and pi
import numpy as np

# FTCS for the first time step
from FTCSSchemes import *

def CTCS(phiOld, c, nt):
    "Linear advection of profile in phiOld using CTCS, Courant number c"
    "for nt time-steps"
    
    nx = len(phiOld)

    phi = FTCS(phiOld, c, 1)

    # new time-step array for phi
    phiNew = phiOld.copy()


    # CTCS for each time-step
    for it in range(nt-1):
        # Loop through all space using remainder after division (%)
        # to cope with periodic boundary conditions
        for j in range(nx):
            phiNew[j] = phiOld[j] - c*\
                     (phi[(j+1)%nx] - phi[(j-1)%nx])
        
        # update arrays for next time-step
        phiOld = phi.copy()
        phi = phiNew.copy()

    return phiNew