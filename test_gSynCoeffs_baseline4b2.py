import simseedburst_func_25
from pylab import *
import sys
gNoiseCoeff = 1.07
gNoiseCoeffI = 1.07

coeff = float(sys.argv[5])
if True:
  gEE = 1.07*coeff
  gEI = float(sys.argv[4])
  gIE = 1.07*coeff
  gII = float(sys.argv[4])
  rdseed = int(sys.argv[6])
  NTTC = 8#80#50#150
  NIN = 2#20#3#50
  Nsavedcells = NTTC+NIN # or 1.. this takes mem.
  #def simseedburst_func(NTTC=1, NIN=1, tstop=10200,mutID=0,rdSeed=1,Econ=0.0004,Icon=0.001,nseg=5,rateCoeff=1.0,gNoiseCoeff=1.0,gNoiseCoeffI=1.0,gSynEE=1.0,gSynEI=1.0,gSynIE=1.0,gSynII=1.0,Ncells2save=1,sparsedt=1.0,Nsyns2save=1,conn\M=[],gToBlock=[],blockEfficiency=0.0):                                                                                                                                                                                                   
  # data=simseedburst_func.simseedburst_func(NTTC, NIN, 11000,0,rdseed,0.00039,0.0006,5,1.0,gNoiseCoeff,gNoiseCoeffI,gEE,gEI,gIE,gII,1,1.0,1)
  data=simseedburst_func_25.simseedburst_func(NTTC, NIN, 11000,0,rdseed,0.00039,0.0006,5,1.0,gNoiseCoeff,gNoiseCoeffI,gEE,gEI,gIE,gII,Nsavedcells,1.0,1)

