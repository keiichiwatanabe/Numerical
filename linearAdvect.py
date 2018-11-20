#!/usr/bin/python3

# Outer code for setting up the linear advection problem on a uniform
# grid and calling the function to perform the linear advection and plot.


### The matplotlib package contains plotting functions              ###
import matplotlib.pyplot as plt

# read in all the linear advection schemes, initial conditions and other
# code associated with this application
from initialConditions import *
from FTBSSchemes import *
from CTCSSchemes import *
from SMSchemes import *
from diagnostics import *
from Errorlists import *


### The main code is inside a function to avoid global variables    ###
def main():
    "Advect the initial conditions using various advection schemes and"
    "compare results"

    # Parameters
    xmin = 0
    xmax = 4
    nx = 201    # Number of points from x=0 to x=4
    nt = 100    # Number of time-steps
    c = 0.8     # Courant number

    # Derived parameters
    dx = (xmax - xmin)/nx
    
    # spatial points for plotting and for defining initial conditions
    x = np.arange(xmin, xmax, dx)

    # Initial conditions for cosine bell
    BellOld = cosBell(x, 0, 0.75)

    # Initial conditions for square wave
    SquareOld = squareWave(x, 0, 0.6)

    # Exact solution is the initial condition shifted around the domain
    BellAnalytic = cosBell((x - c*nt*dx)%(xmax - xmin), 0, 0.75)
    SquareAnalytic = squareWave((x - c*nt*dx)%(xmax - xmin), 0, 0.6)

    # Advect the profile using finite difference for all the time steps
    BellFTBS = FTBS(BellOld.copy(), c, nt)
    BellCTCS = CTCS(BellOld.copy(), c, nt)
    BellSM = SM(BellOld.copy(), c, nt)

    SquareFTBS = FTBS(SquareOld.copy(), c, nt)
    SquareCTCS = CTCS(SquareOld.copy(), c, nt)
    SquareSM = SM(SquareOld.copy(), c, nt)
    
    # Calculate and print out error norms for cosine bell
    print("Bell FTBS l2 error norm = ", l2ErrorNorm(BellFTBS, BellAnalytic))
    print("Bell FTBS linf error norm = ", lInfErrorNorm(BellFTBS, BellAnalytic))
    print("Bell CTCS l2 error norm = ", l2ErrorNorm(BellCTCS, BellAnalytic))
    print("Bell CTCS linf error norm = ", lInfErrorNorm(BellCTCS, BellAnalytic))
    print("Bell CSM l2 error norm = ", l2ErrorNorm(BellSM, BellAnalytic))
    print("Bell CSM linf error norm = ", lInfErrorNorm(BellSM, BellAnalytic))

    # Calculate and print out error norms for square wave
    print("Square FTBS l2 error norm = ", l2ErrorNorm(SquareFTBS, SquareAnalytic))
    print("Square FTBS linf error norm = ", lInfErrorNorm(SquareFTBS, SquareAnalytic))
    print("Square CTCS l2 error norm = ", l2ErrorNorm(SquareCTCS, SquareAnalytic))
    print("Square CTCS linf error norm = ", lInfErrorNorm(SquareCTCS, SquareAnalytic))
    print("Square CSM l2 error norm = ", l2ErrorNorm(SquareSM, SquareAnalytic))
    print("Square CSM linf error norm = ", lInfErrorNorm(SquareSM, SquareAnalytic))

    # Plot the solutions for cosine bell
    font = {'size'   : 12}
    plt.rc('font', **font)
    plt.figure(1)
    plt.clf()
    plt.ion()
    plt.plot(x, BellOld, label='Initial', color='black')
    plt.plot(x, BellAnalytic, label='Analytic', color='black', 
             linestyle='--', linewidth=2)
    plt.plot(x, BellFTBS, label='FTBS', color='blue')
    plt.plot(x, BellCTCS, label='CTCS', color='green')
    plt.plot(x, BellSM, label='CSM', color='red')
    plt.axhline(0, linestyle=':', color='black')
    plt.ylim([-0.2,1.2])
    plt.legend(bbox_to_anchor=(1, 1))
    plt.xlabel('$x$')
    input('press return to save file for plots of cosine bell and continue')
    plt.savefig('plots/Plots_Bell.pdf')

    # Plot the solutions for square wave
    font = {'size'   : 12}
    plt.rc('font', **font)
    plt.figure(1)
    plt.clf()
    plt.ion()
    plt.plot(x, SquareOld, label='Initial', color='black')
    plt.plot(x, SquareAnalytic, label='Analytic', color='black', 
             linestyle='--', linewidth=2)
    plt.plot(x, SquareFTBS, label='FTBS', color='blue')
    plt.plot(x, SquareCTCS, label='CTCS', color='green')
    plt.plot(x, SquareSM, label='CSM', color='red')
    plt.axhline(0, linestyle=':', color='black')
    plt.ylim([-0.2,1.2])
    plt.legend(bbox_to_anchor=(1 , 1))
    plt.xlabel('$x$')
    input('press return to save file for plots of square wave and continue')
    plt.savefig('plots/Plots_Square.pdf')

    # Plot
    C = [-1.5, -1.25, -1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1, 1.25,\
       1.5, 1.75, 2, 2.25, 2.5, 2.75, 3]
    cSize = len(C)

    # Calculate the error norms for square wave
    FTBSBellplots = np.zeros(cSize)
    for i in range(cSize):
        Bell_FTBS = FTBS(BellOld.copy(), C[0] + 0.25*i, nt)
        FTBSBellplots[i] = np.log(l2ErrorNorm(Bell_FTBS, BellAnalytic))
    

    CTCSBellplots = np.zeros(cSize)
    for i in range(cSize):
        Bell_CTCS = CTCS(SquareOld.copy(), C[0] + 0.25*i, nt)
        CTCSBellplots[i] = np.log(l2ErrorNorm(Bell_CTCS, BellAnalytic))

    SMBellplots = np.zeros(cSize)
    for i in range(cSize):
        Bell_SM = SM(BellOld.copy(), C[0] + 0.25*i, nt)
        SMBellplots[i] = np.log(l2ErrorNorm(Bell_SM, BellAnalytic))

    # Plot the L2 errors for square wave    
    ax = plt.gca()
    ax.set_yscale('log') 
    font = {'size'   : 12}
    plt.rc('font', **font)
    plt.figure(1)
    plt.clf()
    plt.ion()
    plt.plot(C, FTBSBellplots, label='FTBS', color='blue')
    plt.plot(C, CTCSBellplots, label='CTCS', color='green')
    plt.plot(C, SMBellplots, label='CSM', color='red')
    plt.axhline(0, linestyle=':', color='black')
    plt.ylabel('L2 error with log scale')
    plt.legend(loc="upper left")
    plt.xlabel('Courant number')
    input('press return to save file for plots of L2 error norms of'
    ' Bell function and continue')
    plt.savefig('plots/Plots_Error_Bell.pdf')


    # Calculate the error norms for square wave
    FTBSSquareplots = np.zeros(cSize)
    for i in range(cSize):
        Square_FTBS = FTBS(SquareOld.copy(), C[0] + 0.25*i, nt)
        FTBSSquareplots[i] = np.log(l2ErrorNorm(Square_FTBS, SquareAnalytic))
    

    CTCSSquareplots = np.zeros(cSize)
    for i in range(cSize):
        Square_CTCS = CTCS(SquareOld.copy(), C[0] + 0.25*i, nt)
        CTCSSquareplots[i] = np.log(l2ErrorNorm(Square_CTCS, SquareAnalytic))

    SMSquareplots = np.zeros(cSize)
    for i in range(cSize):
        Square_SM = SM(SquareOld.copy(), C[0] + 0.25*i, nt)
        SMSquareplots[i] = np.log(l2ErrorNorm(Square_SM, SquareAnalytic))

    # Plot the L2 errors for square wave    
    ax = plt.gca()
    ax.set_yscale('log') 
    font = {'size'   : 12}
    plt.rc('font', **font)
    plt.figure(1)
    plt.clf()
    plt.ion()
    plt.plot(C, FTBSSquareplots, label='FTBS', color='blue')
    plt.plot(C, CTCSSquareplots, label='CTCS', color='green')
    plt.plot(C, SMSquareplots, label='CSM', color='red')
    plt.axhline(0, linestyle=':', color='black')
    plt.ylabel('L2 error with log scale')
    plt.legend(loc="upper left")
    plt.xlabel('Courant number')
    input('press return to save file for plots of L2 error norms of'
    ' Square wave and continue')
    plt.savefig('plots/Plots_Error_Square.pdf')
    


### Run the function main defined in this file                      ###
main()

