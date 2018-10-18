# Numerical schemes for simulating linear advection for outer code
# linearAdvect.py 

# If you are using Python 2.7 rather than Python 3, import various
# functions from Python 3 such as to use real number division
# rather than integer division. ie 3/2  = 1.5  rather than 3/2 = 1
#from __future__ import absolute_import, division, print_function

# The numpy package for numerical functions and pi
import numpy as np

def CIP(phiOld, phiOldDash, c, dx, nt):
    "Linear advection of profile in phiOld using Semi-Lagrangian, Courant number c"
    "for nt time-steps"



    nx = len(phiOld)

    # new time-step array for phi
    phi = phiOld.copy()
    phidash = phiOldDash.copy()

    # CIP for each time-step
    for it in range(nt):
        # Loop through all space using remainder after division (%)
        # to cope with periodic boundary conditions
        for j in range(nx):
            phi[j] = 2*(phiOld[j%nx] - phiOld[(j-1)%nx])*c*c*c - (phiOldDash[j] + phiOldDash[(j-1)%nx])*c*c*c*dx \
             - 3*(phiOld[j] - phiOld[(j-1)%nx])*c*c + (2*phiOldDash[j] + phiOldDash[(j-1)%nx])*c*c*dx \
             - phiOldDash[j]*c*dx + phiOld[j]

            phidash[j] = - 6*c*c*(phiOld[j] - phiOld[(j-1)%nx])/dx + 3*c*c*(phiOldDash[j] + phiOldDash[(j-1)%nx]) \
             + 6*c*(phiOld[j] - phiOld[(j-1)%nx])/dx - 2*c*(2*phiOldDash[j] + phiOldDash[(j-1)%nx]) \
             + phiOldDash[j]

        # update arrays for next time-step
        phiOld = phi.copy()
        phiOldDash = phidash.copy()

    return phi

