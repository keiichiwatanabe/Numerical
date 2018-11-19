# Numerical schemes for simulating linear advection for outer code
# linearAdvect.py 

# If you are using Python 2.7 rather than Python 3, import various
# functions from Python 3 such as to use real number division
# rather than integer division. ie 3/2  = 1.5  rather than 3/2 = 1
#from __future__ import absolute_import, division, print_function

# The numpy package for numerical functions and pi
import numpy as np

def Errorlists(phi, phiDash, phiAnalytic, dx, nt):
    "Linear advection of profile in phiOld using FTBS, Courant number c"
    "for nt time-steps"

    FTBSplots = [] # initialized the list of FTBS L2 error norm
    CTCSplots = [] # initialized the list of CTCS L2 error norm
    CIPplots  = [] # initialized the list of CIP L2 error norm

    # Advect the profile using finite difference for all the time steps
    for i in range(14):
        phiFTBS = FTBS(phi, 0.5*(-7 + i), nt)
        FTBSplots.extend(l2ErrorNorm(phiFTBS, phiAnalytic))

        phiCTCS = CTCS(phi, 0.5*(-7 + i), nt)
        CTCSplots.extend(l2ErrorNorm(phiCTCS, phiAnalytic))

        phiCIP  = CIP(phi, phiDash, 0.5*(-4 + i), dx, nt)
        CIPplots.extend(l2ErrorNorm(phiCIP, phiAnalytic))
    return 