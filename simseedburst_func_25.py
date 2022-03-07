print("simseedburst started")
from mpi4py import MPI
print("MPI loaded")
from neuron import h
print("NEURON loaded")
import matplotlib
matplotlib.use('Agg')
print("matplotlib loaded")
import numpy
from pylab import *
import time
import scipy.io
import pickle
import sys
import mutation_stuff
import approxhaynetstuff
import mytools
import resource

#9.12.2016: Copied from ../approxhaynet. Make use the updated (groupsyn-version) TTC.hoc

def simseedburst_func(NTTC=1, NIN=1, tstop=10200,mutID=0,rdSeed=1,Econ=0.0004,Icon=0.001,nseg=5,rateCoeff=1.0,gNoiseCoeff=1.0,gNoiseCoeffI=1.0,gSynEE=1.0,gSynEI=1.0,gSynIE=1.0,gSynII=1.0,Ncells2save=1,sparsedt=1.0,Nsyns2save=1,connM=[],gToBlock=[],blockEfficiency=0.0, weight_I_NMDA=1.0):
  MEM_used = []
  MEM_used.append(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss); print("MEM now "+str(MEM_used[-1]))
  myrandsavemat = int(100000*gSynEE+1000*NTTC+rdSeed)

  rank = MPI.COMM_WORLD.Get_rank()
  dt = 0.025

  coeffCoeffs = [[0.25,0],[0.125,0],[0.5,0],[0.5,1.0/3],[0.5,2.0/3],[0.5,1.0],[-0.25,0],[-0.125,0],[-0.5,0]]
  MT = mutation_stuff.getMT(False)
  defVals = mutation_stuff.getdefvals()
  keyList = list(defVals.keys())
  for idefval in range(0,len(keyList)):
    if type(defVals[keyList[idefval]]) is not list:
      defVals[keyList[idefval]] = [defVals[keyList[idefval]], defVals[keyList[idefval]]] #make the dictionary values [somatic, apical]
  updatedVars = ['somatic','apical','basal'] # the possible classes of segments that defVals may apply to
  whichDefVal = [0,1,0]                      # use the defVal[0] for somatic and basal segments and defVal[1] for apical segments
  unpicklefile = open('scalings_cs.sav', 'rb')
  unpickledlist = pickle.load(unpicklefile,encoding='bytes')
  unpicklefile.close()
  theseCoeffsAllAll = unpickledlist[0]

  filename = 'pars_withmids_combfs_final'
  unpicklefile = open(filename+".sav", 'rb')
  unpickledlist = pickle.load(unpicklefile,encoding='bytes')
  unpicklefile.close()
  par_names = unpickledlist[0]
  par_values = unpickledlist[1]
  paramdict = {}
  for i in range(0,len(par_names)):
    paramdict[par_names[i].decode("utf-8")] = par_values[i]

  h("""
{load_file("stdlib.hoc")}
{load_file("stdrun.hoc")}

initialization_tstart = startsw()

strdef fileName
objref fileObj

fileObj = new File()

rdSeed = """+str(rdSeed)+"""
NTTC = """+str(NTTC)+""" //Number of thick-tufted cells
NIN = """+str(NIN)+"""   //Number of interneurons
connectivity = 1
noisyst = 0

tstop = """+str(tstop)+"""
rcpWeightFactor = 1.5 // the factor by which reciprocal weights are stronger than unidirectional weights
pT2Tr = 0.06 //probability of reciprocating an existing connection to another L5bPC
pT2T = 0.13 //probability of a L5bPC being connected to another L5bPC
pT2I = 0.04 //probability of a L5PC being connected to inhibitory neuron, estimated from Figure 7B in Markram et al. 2015
pI2T = 0.09 //probability of an inhibitory neuron being connected to L5PC, estimated from Figure 7B in Markram et al. 2015
pI2I = 0.06 //probability of an inhibitory neuron being connected to another inhibitory neuron, estimated from Figure 7B in Markram et al. 2015

Econ = """+str(Econ)+""" //excitatory synaptic conductance
Icon = """+str(Icon)+""" //inhibitory synaptic conductance
NcontE = 5 // number of excitatory synaptic contacts per connection (to L5PC)
NcontI = 20 // number of inhibitory synaptic contacts per connection (to L5PC), estimated from Figure 7A in Markram et al. 2015
NcontE_I = 7 // number of excitatory synaptic contacts per connection (to interneuron), estimated from Figure 7A in Markram et al. 2015
NcontI_I = 17 // number of excitatory synaptic contacts per connection (to interneuron), estimated from Figure 7A in Markram et al. 2015

NsynE = 10000 // number of excitatory synapses
NsynI = 2500 // number of inhibitory synapses
NsynE_IN = 2400 // number of excitatory synapses
NsynI_IN = 600 // number of inhibitory synapses
gNoiseCoeff = """+str(gNoiseCoeff)+""" // scaling of background synaptic conductances
gNoiseCoeffI = """+str(gNoiseCoeffI)+""" // scaling of background synaptic conductances
rateE = """+str(0.72*rateCoeff)+""" // average rate of presynaptic excitatory cells
rateI = """+str(7.0*rateCoeff)+""" // average rate of presynaptic inhibitory cells
mainBifurcation = 650

{Ncells2save = """+str(Ncells2save)+"""}
sparsedt = """+str(sparsedt)+""" // recordings of [Ca], I_SK and vApical are done with low temporal resolution to save memory
{gSynEE = """+str(gSynEE)+"""}
{gSynEI = """+str(gSynEI)+"""}
{gSynIE = """+str(gSynIE)+"""}
{gSynII = """+str(gSynII)+"""}

objref tempvec
tempvec = new Vector()
{tempvec.append(NTTC)}
{Ncells2save = """+str(Ncells2save)+"""}
""")
  print("Params OK!")

  h("""
{load_file("models/TTC.hoc")}
{load_file("models/sIN.hoc")}

objref MC_TTC, MC_sIN
objref sl //synaptic locations list

objref rds1,rds2,rds3
{rds1 = new Random(1000*rdSeed)}
{rds1.uniform(0,1)} //random for microcircuit connectivity and noisyst

objref conMatEE, conMatEI, conMatIE, conMatII
conMatEE = new Matrix(NTTC,NTTC) //E -> E
conMatEI = new Matrix(NTTC,NIN)  //I -> E
conMatIE = new Matrix(NIN,NTTC)  //E -> I
conMatII = new Matrix(NIN,NIN)   //I -> I

for(i=0;i<NTTC;i+=1){
        conMatEE.x[i][i]=0
}
for(i=0;i<NIN;i+=1){
        conMatII.x[i][i]=0
}
""")
  MEM_used.append(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss); print("MEM now "+str(MEM_used[-1]))
  if len(connM) == 0:
    h("""
for(i=0;i<(NTTC-1);i+=1){
        for(j=(i+1);j<NTTC;j+=1){
                if (connectivity){
                        pcon = rds1.repick()
                        if (pcon<pT2Tr){
                                conMatEE.x[i][j]=rcpWeightFactor*gSynEE
                                conMatEE.x[j][i]=rcpWeightFactor*gSynEE
                        } else {
                                if (pcon<(pT2Tr + 0.5*pT2T)){
                                        conMatEE.x[i][j]=gSynEE
                                        conMatEE.x[j][i]=0
                                } else {
                                        if (pcon<(pT2Tr + pT2T)){
                                                conMatEE.x[i][j]=0
                                                conMatEE.x[j][i]=gSynEE
                                        } else {
                                                conMatEE.x[i][j]=0
                                                conMatEE.x[j][i]=0
                                        }
                                }
                        }
                } else {
                        conMatEE.x[i][j]=0
                        conMatEE.x[j][i]=0
                }
        }
}
""")

    h("""
for(i=0;i<NTTC;i+=1){
        for(j=0;j<NIN;j+=1){
                if (connectivity){
                        pcon = rds1.repick()
                        if (pcon<pT2I){
                                conMatEI.x[i][j]=gSynEI
                        } else {
                                conMatEI.x[i][j]=0
                        }
                } else {
                        conMatEI.x[i][j]=0
                }
        }
}
""")

    h("""
for(i=0;i<NIN;i+=1){
        for(j=0;j<NTTC;j+=1){
                if (connectivity){
                        pcon = rds1.repick()
                        if (pcon<pI2T){
                                conMatIE.x[i][j]=gSynIE
                        } else {
                                conMatIE.x[i][j]=0
                        }
                } else {
                        conMatIE.x[i][j]=0
                }
        }
}
""")

    h("""
for(i=0;i<NIN;i+=1){
        for(j=0;j<NIN;j+=1){
                if (connectivity && i!=j){
                        pcon = rds1.repick()
                        if (pcon<pI2I){
                                conMatII.x[i][j]=gSynII
                        } else {
                                conMatII.x[i][j]=0
                        }
                } else {
                        conMatII.x[i][j]=0
                }
        }
}
""")
  else: #TODO: also for conMatEI, IE, II
    for i in range(0,NTTC):
      for j in range(0,NTTC):
        if connM[i][j]:
          h("conMatEE.x["+str(i)+"]["+str(j)+"]="+str(gSynEE*connM[i][j]))
      for j in range(0,NIN):
        if connM[i][NTTC+j]:
          h("conMatEI.x["+str(i)+"]["+str(j)+"]="+str(gSynEI*connM[i][NTTC+j]))
    for i in range(0,NIN):
      for j in range(0,NTTC):
        if connM[NTTC+i][j]:
          h("conMatIE.x["+str(i)+"]["+str(j)+"]="+str(gSynIE*connM[NTTC+i][j]))
      for j in range(0,NIN):
        if connM[NTTC+i][NTTC+j]:
          h("conMatII.x["+str(i)+"]["+str(j)+"]="+str(gSynII*connM[NTTC+i][NTTC+j]))
  print("Connectivity OK!")
  MEM_used.append(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss); print("MEM now "+str(MEM_used[-1]))

  h("""
{load_file("netparmpi.hoc")}
objref epnm

epnm = new ParallelNetManager(NTTC+NIN)
{epnm.round_robin()}

strdef treename
objref NsynsE, NsynsI, preTrainList
objref NsynsE_IN, NsynsI_IN // JFK
""")

  for i in range(0,NTTC):
    h("""
  i = """+str(i)+"""
  if (epnm.gid_exists(i)) {
        print \"rank = """+str(rank)+""", gid \", i, \" exists\\n\"
        MC_TTC = new TTC()
        epnm.register_cell(i,MC_TTC)
        epnm.pc.gid2cell(i).initRand(1000*rdSeed+i)
        epnm.pc.gid2cell(i).setnetworkparameters(rcpWeightFactor,Econ,Icon,NsynE,NsynI,NcontE,NcontI,1.0,1.0,1.0,gNoiseCoeff)
  }""")
  for i in range(0,NIN):
    h("""
  i = """+str(i)+"""
  if (epnm.gid_exists(NTTC+i)) {
        print \"rank = """+str(rank)+""", gid \", NTTC+i, \" exists\\n\"
        MC_sIN = new sIN()
        epnm.register_cell(NTTC+i,MC_sIN)
        epnm.pc.gid2cell(NTTC+i).initRand(1000*rdSeed+NTTC+i)
        epnm.pc.gid2cell(NTTC+i).setnetworkparameters(1.0,Econ,Icon,NsynE_IN,NsynI_IN,1.0,"""+str(weight_I_NMDA)+""",1.0,gNoiseCoeffI,NcontE_I,NcontI_I)
  }""")
    print("ParallelNetManager " + str(i)+ " OK!")

  approxhaynetstuff.setparams(paramdict,NTTC) # TMM should i have one for NIN? and does this contain all now that i have a different INpop?

  for i in range(0,NTTC):
    h("""
  i = """+str(i)+"""
  if (epnm.gid_exists(i)) {
    lengthA = epnm.pc.gid2cell(i).apic[0].L + epnm.pc.gid2cell(i).apic[1].L
    lengthB = epnm.pc.gid2cell(i).dend.L
    //print "NTTC ", i, epnm.pc.gid2cell(i)
    pA = lengthA/(lengthA + lengthB)
    {NsynsE = new List()}
    {NsynsI = new List()}
    for i1 = 0, 2 {
      //print "NTTC i1, nseg ",i1, """+str(nseg)+"""
      {NsynsE.append(new Vector("""+str(nseg)+"""))}
      {NsynsI.append(new Vector("""+str(nseg)+"""))}
    }
    // print "NTTC NsynsE ", NsynsE

    for(i1=0;i1<(NsynE+NsynI);i1+=1){
      if (epnm.pc.gid2cell(i).rd1.repick()<pA){
        treename = "apic"
        compInd = 1
      } else {
        treename = "dend"
        compInd = 0
      }

      sl = epnm.pc.gid2cell(i).locateSites(treename,epnm.pc.gid2cell(i).rd1.repick()*epnm.pc.gid2cell(i).getLongestBranch(treename))
      sitenum = int((sl.count()-1)*epnm.pc.gid2cell(i).rd1.repick())
      compInd = compInd + sl.o[sitenum].x[0] // if we are at apical, and sl.o[sitenum].x[0]=1, then compInd = 2, otherwise 1 at apical, and 0 at basal
      segInd = int(sl.o[sitenum].x[1]*"""+str(nseg)+""")
      if (i1<NsynE) {
        NsynsE.o[compInd].x[segInd] = NsynsE.o[compInd].x[segInd] + 1
      } else {
        NsynsI.o[compInd].x[segInd] = NsynsI.o[compInd].x[segInd] + 1
      }
    }

    {epnm.pc.gid2cell(i).distributeSyn(NsynsE,NsynsI)}
    print "distributeSyn(NsynsE,NsynsI) OK"

    {preTrainList = new List()}
    {rds2 = new Random(1000*rdSeed+i)}//random for presynaptic trains
    {rds3 = new Random(1000*rdSeed+i)}//random for presynaptic trains
    {rds2.negexp(1/rateE)}
    {rds3.negexp(1/rateI)}
    print "rds OK"

    for(compInd=0;compInd<3;compInd+=1){
      for(segInd=0;segInd<"""+str(nseg)+""";segInd+=1){
        {preTrainList.append(new Vector())}
        pst=0 //presynaptic spike time
        if(NsynsE.o[compInd].x[segInd]==0) {
          print "Warning: NsynsE.o[",compInd,"].x[",segInd,"] = 0!!!!"
          pst = 1e6
        }
        while(pst < tstop){
          pst+= 1000*rds2.repick()/NsynsE.o[compInd].x[segInd]
          {preTrainList.o[preTrainList.count()-1].append(pst)}
        }
      }
    }
    for(compInd=0;compInd<3;compInd+=1){
      for(segInd=0;segInd<"""+str(nseg)+""";segInd+=1){
        {preTrainList.append(new Vector())}
        pst=0 //presynaptic spike time
        if(NsynsI.o[compInd].x[segInd]==0) {
          print "Warning: NsynsI.o[",compInd,"].x[",segInd,"] = 0!!!!"
          pst = 1e6
        }
        while(pst < tstop){
          pst+= 1000*rds3.repick()/NsynsI.o[compInd].x[segInd]
          {preTrainList.o[preTrainList.count()-1].append(pst)}
        }
      }
    }
    {epnm.pc.gid2cell(i).setpretrains(preTrainList)}
    {epnm.pc.gid2cell(i).queuePreTrains()}
    print \"distributeSyn(), setpretrains(), queuePreTrains() OK! i=\", i, \".\"
  }
""")


  for i in range(0,NIN):
    h("""
  i = """+str(i)+"""
  if (epnm.gid_exists(NTTC+i)) {
    //print "NIN ", NTTC+i, epnm.pc.gid2cell(NTTC+i)
    pA = 0.5

    {NsynsE_IN = new List()} // TMM I had to make this list, as well as the one in the next line
    {NsynsI_IN = new List()}

    for i1 = 0, 1 { // TMM why 0,2?
      //print "NIN i1,nseg ",i1, """+str(nseg)+"""
      {NsynsE_IN.append(new Vector("""+str(nseg)+"""))}
      {NsynsI_IN.append(new Vector("""+str(nseg)+"""))}
    }

    for(i1=0;i1<(NsynE_IN+NsynI_IN);i1+=1){
      if (epnm.pc.gid2cell(NTTC+i).rd1.repick()<pA){
        treename = "dend[1]"
        compInd = 1
      } else {
        treename = "dend[0]"
        compInd = 0
      }
      sl = epnm.pc.gid2cell(NTTC+i).locateSites(treename,epnm.pc.gid2cell(NTTC+i).rd1.repick()*epnm.pc.gid2cell(NTTC+i).getLongestBranch(treename))
      sitenum = int((sl.count()-1)*epnm.pc.gid2cell(NTTC+i).rd1.repick())
      //compInd = compInd + sl.o[sitenum].x[0] // if we are at apical, and sl.o[sitenum].x[0]=1, then compInd = 2, otherwise 1 at apical, and 0 at basal
      segInd = int(sl.o[sitenum].x[1]*"""+str(nseg)+""")
      if (i1<NsynE_IN) {
        NsynsE_IN.o[compInd].x[segInd] = NsynsE_IN.o[compInd].x[segInd] + 1
      } else {
        NsynsI_IN.o[compInd].x[segInd] = NsynsI_IN.o[compInd].x[segInd] + 1
      }
    }    
    {epnm.pc.gid2cell(NTTC+i).distributeSyn(NsynsE_IN,NsynsI_IN)}
    print "distributeSyn(NsynsE_IN,NsynsI_IN) OK"

    {preTrainList = new List()}
    {rds2 = new Random(1000*rdSeed+i)}//random for presynaptic trains
    {rds3 = new Random(1000*rdSeed+i)}//random for presynaptic trains
    {rds2.negexp(1/rateE)}
    {rds3.negexp(1/rateI)}
    print "rds OK"


    /*
    {preTrainList.append(new Vector())}
    pst=0 //presynaptic spike time
    while(pst < tstop){
      pst+= 1000*rds2.repick()/NsynE_IN
      {preTrainList.o[preTrainList.count()-1].append(pst)}
    }

    {preTrainList.append(new Vector())}
    pst=0 //presynaptic spike time
    while(pst < tstop){
      pst+= 1000*rds3.repick()/NsynI_IN
      {preTrainList.o[preTrainList.count()-1].append(pst)}
    }
    */

    /// JFK
    for(compInd=0;compInd<2;compInd+=1){ //TMM why to 3, also added this whole block
      for(segInd=0;segInd<"""+str(nseg)+""";segInd+=1){
        {preTrainList.append(new Vector())}
        pst=0 //presynaptic spike time
        if(NsynsE_IN.o[compInd].x[segInd]==0) {
          print "Warning: NsynsE_IN.o[",compInd,"].x[",segInd,"] = 0!!!!"
          pst = 1e6
        }
        while(pst < tstop){
          pst+= 1000*rds2.repick()/NsynsE_IN.o[compInd].x[segInd]
          {preTrainList.o[preTrainList.count()-1].append(pst)}
        }
      }
    }
    
    for(compInd=0;compInd<2;compInd+=1){
      for(segInd=0;segInd<"""+str(nseg)+""";segInd+=1){
        {preTrainList.append(new Vector())}
        pst=0 //presynaptic spike time
        if(NsynsI_IN.o[compInd].x[segInd]==0) {
          print "Warning: NsynsI_IN.o[",compInd,"].x[",segInd,"] = 0!!!!"
          pst = 1e6
        }
        while(pst < tstop){
          pst+= 1000*rds3.repick()/NsynsI_IN.o[compInd].x[segInd]
          {preTrainList.o[preTrainList.count()-1].append(pst)}
        }
      }
    }
    /// KFJ


    {epnm.pc.gid2cell(NTTC+i).setpretrains(preTrainList)}
    {epnm.pc.gid2cell(NTTC+i).queuePreTrains()}

    print \"distributeSyn(), setpretrains(), queuePreTrains() OK! i=\", i, \".\"
  }
""")

  print("Spike trains OK!")

  # sys.exit()

  MEM_used.append(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss); print("MEM now "+str(MEM_used[-1])); scipy.io.savemat('MEM_used_'+str(myrandsavemat)+'.mat',{'MEM_used': MEM_used})

  #LOAD MUTATIONS
  counter = -1
  found = 0
  for igene in range(0, len(MT)):
    for imut in range(0, len(MT[igene])):
      theseCoeffs = theseCoeffsAllAll[0][igene][imut]

      nVals = len(MT[igene][imut])*[0]
      thesemutvars = []
      for imutvar in range(0,len(MT[igene][imut])):
        thesemutvars.append(MT[igene][imut][imutvar][0])
        if type(MT[igene][imut][imutvar][1]) is int or type(MT[igene][imut][imutvar][1]) is float:
          MT[igene][imut][imutvar][1] = [MT[igene][imut][imutvar][1]]
        nVals[imutvar] = len(MT[igene][imut][imutvar][1])
      cumprodnVals = cumprod(nVals)
      allmutvars = cumprodnVals[len(MT[igene][imut])-1]*[thesemutvars]
      allmutvals = []
      for iallmutval in range(0,cumprodnVals[len(MT[igene][imut])-1]):
        allmutvals.append([0]*len(thesemutvars))
      for iallmutval in range(0,cumprodnVals[len(MT[igene][imut])-1]):
        for imutvar in range(0,len(MT[igene][imut])):
          if imutvar==0:
            allmutvals[iallmutval][imutvar] = MT[igene][imut][imutvar][1][iallmutval%nVals[imutvar]]
          else:
            allmutvals[iallmutval][imutvar] = MT[igene][imut][imutvar][1][int(iallmutval/cumprodnVals[imutvar-1])%nVals[imutvar]]

      for iallmutval in range(0, len(theseCoeffs)):
        if igene == 0 and imut == 0 and iallmutval == 0:
          iters = [-1,0,2,6,8]
        else:
          iters = [0,2,6,8]
        for iiter in range(0,len(iters)):
          iter = iters[iiter]
          counter = counter + 1
          if counter == mutID:
            found = 1
            break
        if found:
          break
      if found:
        break
    if found:
      break
  if not found:
    print("Not found corresponding mutID, mutID = " + str(mutID))
    sys.exit()

  if iter >= 0:
    thisCoeff = coeffCoeffs[iter][0]*theseCoeffs[iallmutval] + coeffCoeffs[iter][1]*(1.0 - 0.5*theseCoeffs[iallmutval])
  else:
    thisCoeff = 0
  print("iter="+str(iter)+", thisCoeff="+str(thisCoeff))

  mutText = ""
  for imutvar in range(0,len(MT[igene][imut])):
    if imutvar > 0 and imutvar%2==0:
      mutText = mutText+"\n"
    mutvars = allmutvars[iallmutval][imutvar]
    mutvals = allmutvals[iallmutval][imutvar]
    if type(mutvars) is str:
      mutvars = [mutvars]
    mutText = mutText + str(mutvars) + ": "
    for kmutvar in range(0,len(mutvars)):
      mutvar = mutvars[kmutvar]
      if mutvar.find('offm') > -1 or mutvar.find('offh') > -1 or mutvar.find('ehcn') > -1:
        newVal =  [x+mutvals*thisCoeff for x in defVals[mutvar]]
        if mutvals >= 0 and kmutvar==0:
          mutText = mutText + "+" + str(mutvals) +" mV"
        elif kmutvar==0:
          mutText = mutText  + str(mutvals) +" mV"
      else:
        newVal = [x*(mutvals**thisCoeff) for x in defVals[mutvar]]
        if kmutvar==0:
          mutText = mutText + "*" + str(mutvals)
      if kmutvar < len(mutvars)-1:
        mutText = mutText + ", "
      if mutvar.find('_Ih') > -1:
        updateThese = [1,1,1]
      elif mutvar.find('_Ca_HVA') > -1 or mutvar.find('_Ca_LVAst') > -1 or mutvar.find('_SKv3.1') > -1 or mutvar.find('_Ca_HVA') > -1 or mutvar.find('_SK_E2') > -1 or mutvar.find('_NaTa_t') > -1 or mutvar.find('_CaDynamics_E2') > -1:
        updateThese = [1,1,0]
      elif mutvar.find('_K_Pst') > -1 or mutvar.find('_K_Tst') > -1 or mutvar.find('_Nap_Et2') > -1:
        updateThese = [1,0,0]
      elif mutvar.find('_Im') > -1:
        updateThese = [0,1,0]
      else:
        print("Error: str=" + str(mutvar))
        updatedThese = [0,0,0]
      for iupdated in range(0,3):
        if updateThese[iupdated]:
          print("""forsec epnm.pc.gid2cell(i)."""+str(updatedVars[iupdated])+""" {
"""+mutvar+""" = """+str(newVal[whichDefVal[iupdated]])+"""
}""") 
          for i in range(0,NTTC):
            h("""
i = """+str(i)+"""
if (epnm.gid_exists(i)) {
  forsec epnm.pc.gid2cell(i)."""+str(updatedVars[iupdated])+""" {
    """+mutvar+""" = """+str(newVal[whichDefVal[iupdated]])+"""
  }
}""")

  print(mutText)
  h("""
thisCa = 0.0001
for(i=0;i<NTTC;i+=1){
  if (epnm.gid_exists(i)) {
    thisCa = epnm.pc.gid2cell(i).soma.minCai_CaDynamics_E2
  }
}
""")
  thisCa = h.thisCa

  myMechs = ['Ca_HVA','Ca_LVAst','Ih','Im','K_Pst','K_Tst','NaTa_t','Nap_Et2','SK_E2','SKv3_1','']
  myMechToAdd = ""
  for iblock in range(0,len(gToBlock)):
    for iMech in range(0,len(myMechs)):
      if gToBlock[iblock] in myMechs[iMech]:
        break
    if iMech <= 9:
      if type(blockEfficiency) is list:
        for i in range(0,NTTC):
          h("""
i = """+str(i)+"""
if (epnm.gid_exists(i)) {
  epnm.pc.gid2cell(i).soma g"""+str(myMechs[iMech])+"""bar_"""+str(myMechs[iMech])+""" = g"""+str(myMechs[iMech])+"""bar_"""+str(myMechs[iMech])+""" * """+str(blockEfficiency[0])+"""
  forsec epnm.pc.gid2cell(i).apical { g"""+str(myMechs[iMech])+"""bar_"""+str(myMechs[iMech])+""" = g"""+str(myMechs[iMech])+"""bar_"""+str(myMechs[iMech])+""" * """+str(blockEfficiency[1])+""" }
}""")
        print(("""
i = """+str(i)+"""
if (epnm.gid_exists(i)) {
  epnm.pc.gid2cell(i).soma g"""+str(myMechs[iMech])+"""bar_"""+str(myMechs[iMech])+""" = g"""+str(myMechs[iMech])+"""bar_"""+str(myMechs[iMech])+""" * """+str(blockEfficiency[0])+"""
  forsec epnm.pc.gid2cell(i).apical { g"""+str(myMechs[iMech])+"""bar_"""+str(myMechs[iMech])+""" = g"""+str(myMechs[iMech])+"""bar_"""+str(myMechs[iMech])+""" * """+str(blockEfficiency[1])+""" }
}"""))
        myMechToAdd = myMechToAdd+myMechs[iMech]+'x'+str(blockEfficiency[0])+'-'+str(blockEfficiency[1])+'_'        
      else:
        print(("forall if(ismembrane(\""+str(myMechs[iMech])+"\")) { g"+str(myMechs[iMech])+"bar_"+str(myMechs[iMech])+" = g"+str(myMechs[iMech])+"bar_"+str(myMechs[iMech])+" * "+str(blockEfficiency)+" }"))
        h("forall if(ismembrane(\""+str(myMechs[iMech])+"\")) { g"+str(myMechs[iMech])+"bar_"+str(myMechs[iMech])+" = g"+str(myMechs[iMech])+"bar_"+str(myMechs[iMech])+" * "+str(blockEfficiency)+" }")
        myMechToAdd = myMechToAdd+myMechs[iMech]+'x'+str(blockEfficiency)+'_'
    else:
      print("Error: No mechanism recognized")
  if weight_I_NMDA != 1.0:
    myMechToAdd = myMechToAdd+'wI_NMDA'+str(weight_I_NMDA)+'_'
  print("mutstuff OK!")

  h("""
v_init = -80
cai0_ca_ion = thisCa
dt = """+str(dt)+"""
objref syninds, conMatRows
for(i=0;i<NTTC;i+=1){
  if (epnm.gid_exists(i)) {
   epnm.pc.gid2cell(i).insertMCcons(conMatEE.getcol(i), conMatIE.getcol(i))
  }
  print "insertMCcons exc OK!"
}
for(i=0;i<NIN;i+=1){
  print "insertMCcons inh OK!"
  if (epnm.gid_exists(NTTC+i)) {
    epnm.pc.gid2cell(NTTC+i).insertMCcons(conMatEI.getcol(i), conMatII.getcol(i))
  }
}

{syninds = new Vector()}
{conMatRows = new List()}
for(i=0;i<NTTC;i+=1){
        syninds.append(2*3*"""+str(nseg)+""")
}
for(i=0;i<NIN;i+=1){
        syninds.append(2)
}

// appending the microcircuit connections
for(j=0;j<NTTC;j+=1){
        conMatRows.append(new Vector())
        for(i=0;i<NTTC;i+=1){
                conMatRows.o[j].insrt(i,conMatEE.x[i][j])
                if (conMatEE.x[i][j] != 0){
                        //print "E->E Synapse at i=", i, ", j=", j, ", syninds.x[j] = ", syninds.x[j], "... ", syninds.x[j]+NcontE-1
                        for(jj=0;jj<NcontE;jj+=1){
                                epnm.nc_append(i,j,syninds.x[j],1,0.5)
                                syninds.x[j] +=1
                        }
                }
        }
        for(i=0;i<NIN;i+=1){
                conMatRows.o[j].insrt(NTTC+i,conMatIE.x[i][j])
                if (conMatIE.x[i][j] != 0){
                        //print "I->E Synapse at i=", NTTC+i, ", j=", j, ", syninds.x[j] = ", syninds.x[j], "... ", syninds.x[j]+NcontI-1
                        for(jj=0;jj<NcontI;jj+=1){
                                epnm.nc_append(NTTC+i,j,syninds.x[j],1,0.5)
                                syninds.x[j] +=1
                        }
                }
        }
}
for(j=0;j<NIN;j+=1){
        conMatRows.append(new Vector())
        for(i=0;i<NTTC;i+=1){
                conMatRows.o[NTTC+j].insrt(i,conMatEI.x[i][j])
                if (conMatEI.x[i][j] != 0){
                        //print "E->I Synapse at i=", i, ", j=", NTTC+j, ", syninds.x[NTTC+j] = ", syninds.x[NTTC+j]
                        epnm.nc_append(i,NTTC+j,syninds.x[NTTC+j],NcontE_I,0.5)
                        syninds.x[NTTC+j] +=1
                }
        }
        for(i=0;i<NIN;i+=1){
                conMatRows.o[NTTC+j].insrt(NTTC+i,conMatII.x[i][j])
                if (conMatII.x[i][j] != 0){
                        //print "I->I Synapse at i=", NTTC+i, ", j=", NTTC+j, ", syninds.x[NTTC+j] = ", syninds.x[NTTC+j]
                        epnm.nc_append(NTTC+i,NTTC+j,syninds.x[NTTC+j],NcontI_I,0.5)
                        syninds.x[j] +=1
                }
        }
}
""")
  h("forall nseg="+str(nseg))
  for i in range(0,NIN):
    h("if (epnm.gid_exists(NTTC+i)) { epnm.pc.gid2cell(NTTC+i).soma nseg = 1 }")
  print("Syninds OK!")
  MEM_used.append(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss); print("MEM now "+str(MEM_used[-1])); scipy.io.savemat('MEM_used_'+str(myrandsavemat)+'.mat',{'MEM_used': MEM_used})

  h("""
objref vSomaList, tvecList, caSomaList, skSomaList, cahvaSomaList, calvaSomaList
objref natSomaList, napSomaList, ihSomaList, kv31SomaList, ktSomaList, kpSomaList, IList
objref apcvecList, apcList, netcon, nil, spikes, spikedCells
objref spikeList //JFK

{spikes = new Vector()}
{spikedCells = new Vector()}
{spikeList = new Vector()}//JFK

{apcvecList = new List()}
{apcList = new List()}
{vSomaList = new List()}
{caSomaList = new List()}
{skSomaList = new List()}
{cahvaSomaList = new List()}
{calvaSomaList = new List()}
{natSomaList = new List()}
{napSomaList = new List()}
{ihSomaList = new List()}
{kv31SomaList = new List()}
{ktSomaList = new List()}
{kpSomaList = new List()}
{IList = new List()}
{tvecList = new List()}
""")
  MEM_used.append(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss); print("MEM now "+str(MEM_used[-1])); scipy.io.savemat('MEM_used_'+str(myrandsavemat)+'.mat',{'MEM_used': MEM_used})

  Nsyns = numpy.array(h.syninds)-2*3*nseg
  cumpNsyns = numpy.cumsum(Nsyns)/numpy.sum(Nsyns)
  randVec = [rand() for x in range(0,Nsyns2save)]
  if sum(Nsyns) > 0:
    cellsSynRecorded = [next(i for i,x in enumerate(cumpNsyns) if x > randVec[j]) for j in range(0,Nsyns2save)]
    synIndsRecorded = [int(2*3*nseg+rand()*Nsyns[i]) for i in cellsSynRecorded]
  else:
    cellsSynRecorded = []
    synIndsRecorded = []

  for i in range(0,int(h.Ncells2save)):
  # for i in [NTTC-1,NTTC,NTTC+1]:
  # for i in range(NTTC,NTTC+int(h.Ncells2save)):

    h("""
      i = """+str(i)+"""
      if (epnm.gid_exists(i)) {
        {vSomaList.append(new Vector())}
        {caSomaList.append(new Vector())}
        {skSomaList.append(new Vector())}
        {cahvaSomaList.append(new Vector())}
        {calvaSomaList.append(new Vector())}
        {natSomaList.append(new Vector())}
        {napSomaList.append(new Vector())}
        {ihSomaList.append(new Vector())}
        {kv31SomaList.append(new Vector())}
        {ktSomaList.append(new Vector())}
        {kpSomaList.append(new Vector())}
        {tvecList.append(new Vector())}

        access epnm.pc.gid2cell(i).soma
        {vSomaList.o[vSomaList.count()-1].record(&v(0.5),dt)}
        if (i < NTTC) {
          {caSomaList.o[caSomaList.count()-1].record(&cai(0.5),sparsedt)}
          {skSomaList.o[skSomaList.count()-1].record(&ik_SK_E2(0.5),sparsedt)}
          {cahvaSomaList.o[skSomaList.count()-1].record(&ica_Ca_HVA(0.5),sparsedt)}
          {calvaSomaList.o[skSomaList.count()-1].record(&ica_Ca_LVAst(0.5),sparsedt)}
          {natSomaList.o[skSomaList.count()-1].record(&ina_NaTa_t(0.5),sparsedt)}
          {napSomaList.o[skSomaList.count()-1].record(&ina_Nap_Et2(0.5),sparsedt)}
          {ihSomaList.o[skSomaList.count()-1].record(&ihcn_Ih(0.5),sparsedt)}
          {kv31SomaList.o[skSomaList.count()-1].record(&ik_SKv3_1(0.5),sparsedt)}
          {ktSomaList.o[skSomaList.count()-1].record(&ik_K_Tst(0.5),sparsedt)}
          {kpSomaList.o[skSomaList.count()-1].record(&ik_K_Pst(0.5),sparsedt)}
        }
      }
    """)
    indSynIndsRecorded = [ix for ix,x in enumerate(cellsSynRecorded) if x==i]
    for isyn in range(0,len(indSynIndsRecorded)):
      h("""
  if (epnm.gid_exists(i)) {
    {IList.append(new Vector())}
    {IList.o[IList.count()-1].record(&epnm.pc.gid2cell(i).synlist.o["""+str(synIndsRecorded[indSynIndsRecorded[isyn]])+"""].i, sparsedt)}
  }
""")
    print("epnm.gid " + str(i)+ " OK!")

  #for i in range(0,
  MEM_used.append(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss); print("MEM now "+str(MEM_used[-1])); scipy.io.savemat('MEM_used_'+str(myrandsavemat)+'.mat',{'MEM_used': MEM_used})

  h("objref nclist, excvec")
  h("nclist = new List()")
  h("excvec = new Vector()")


  # h("objref tobj,timervecc,idvecc,recncses")
  # h("timervecc = new Vector()")
  # h("idvecc = new Vector()")
  # h("recncses = new List()")

  for i in range(0,NTTC+NIN):
    h("""
  i = """+str(i)+"""
  if (epnm.gid_exists(i)) {
          access epnm.pc.gid2cell(i).soma
          {apcList.append(new APCount(0.5))}
          {apcvecList.append(new Vector())}
          apcList.o[apcList.count()-1].thresh= -40
          {apcList.o[apcList.count()-1].record(apcvecList.o[apcList.count()-1])}
          {netcon = new NetCon(&v(0.5), nil)}
          netcon.threshold = -20
          {netcon.record(spikes, spikedCells)  }
          {nclist.append(netcon)}

         // epnm.pc.gid2cell(i).soma tobj = new NetCon(&v(0.5),nil)
         // {tobj.record(timervecc,idvecc,i+1)}
         // {recncses.append(tobj)}




          if(i < NTTC) {
            excvec.append(1)
          } else {
            excvec.append(0)
          }
}""")
    print("epnm.gid " + str(i)+ " OK!")

  print("Connection matrix:")
  print(str(numpy.array(h.conMatRows)))
  MEM_used.append(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss); print("MEM now "+str(MEM_used[-1])); scipy.io.savemat('MEM_used_'+str(myrandsavemat)+'.mat',{'MEM_used': MEM_used})







  # h("objref clomp")
  # h("""
  #   epnm.pc.gid2cell(3).soma clomp = new IClamp(0.5)
  #   clomp.del = 1000
  #   clomp.dur = 10000
  #   clomp.amp = 1000
  #   """)





  # sys.exit()


  h("""
{epnm.set_maxstep(100)}

stdinit()

if (epnm.gid_exists(0)) {
        print \"\\n\"
        sim_tstart = startsw()
        initializationtime = (sim_tstart-initialization_tstart)/3600
        print \"Initialization completed. Initialization took \", initializationtime, \" hours\\n\"
        print \"Starting simulation\\n\"
        print \"\\n\"
}
""")
  MEM_used.append(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss); print("MEM now "+str(MEM_used[-1])); scipy.io.savemat('MEM_used_'+str(myrandsavemat)+'.mat',{'MEM_used': MEM_used})
  h("""
{epnm.psolve(tstop)}

if (epnm.gid_exists(0)) {
        simruntime = (startsw() - sim_tstart)/3600
        print \"Simulation took \", simruntime, \" hours\\n\"
}

""")
  print("Simulation OK!")
  print(h.tstop)
  print("Printing nclist")
  for i in range(0,int(h.nclist.count())):

    h("""

      //print nclist.count()
      print nclist.o(i)
      //print nclist.browser()
    """)


  #times = numpy.array(h.tvecList)
  vSoma = numpy.array(h.vSomaList)
  spikes = numpy.array(h.spikes)
  spikedCells = numpy.array(h.spikedCells)

  picklelist = [spikes,spikedCells,vSoma,[x.hname() for x in h.nclist],array(h.excvec)]
  file = open('spikes_parallel_'+myMechToAdd+str(nseg)+'_'+str(tstop)+'_'+str(mutID)+'_'+str(Econ)+'_'+str(Icon)+'_'+str(rateCoeff)+'_'+str(gNoiseCoeff)+'_'+str(gNoiseCoeffI)+'_'+str(gSynEE)+'_'+str(gSynEI)+'_'+str(gSynIE)+'_'+str(gSynII)+'_'+str(rdSeed)+'_'+str(rank)+'_of_'+str(NTTC+NIN)+'.sav', 'wb')
  pickle.dump(picklelist,file)
  file.close()
  print('Saved to spikes_parallel_'+myMechToAdd+str(nseg)+'_'+str(tstop)+'_'+str(mutID)+'_'+str(Econ)+'_'+str(Icon)+'_'+str(rateCoeff)+'_'+str(gNoiseCoeff)+'_'+str(gNoiseCoeffI)+'_'+str(gSynEE)+'_'+str(gSynEI)+'_'+str(gSynIE)+'_'+str(gSynII)+'_'+str(rdSeed)+'_'+str(rank)+'_of_'+str(NTTC+NIN)+'.sav')
  caSoma = numpy.array(h.caSomaList)
  skSoma = numpy.array(h.skSomaList)
  cahvaSoma = numpy.array(h.cahvaSomaList)
  calvaSoma = numpy.array(h.calvaSomaList)
  natSoma = numpy.array(h.natSomaList)
  napSoma = numpy.array(h.napSomaList)
  ihSoma = numpy.array(h.ihSomaList)
  kv31Soma = numpy.array(h.kv31SomaList)
  ktSoma = numpy.array(h.ktSomaList)
  kpSoma = numpy.array(h.kpSomaList)
  Is = numpy.array(h.IList)
  if len(caSoma) > 0:
    sparseTimes = [sparsedt*x for x in range(0,len(caSoma[0]))]
  else:
    sparseTimes = []


  if len(vSoma) > 0:
    times = [dt*x for x in range(0,len(vSoma[0]))]
    if len(vSoma)==1:
      spikes = mytools.spike_times(times,vSoma.T)
    else:
      spikeTimes = []
      spikeNeurs = []
      for i in range(0,len(vSoma)):
        theseSpikes = mytools.spike_times(times,vSoma[i])
        spikeTimes = spikeTimes + theseSpikes
        spikeNeurs = spikeNeurs + [i for x in theseSpikes]
      spikes = [spikeTimes, spikeNeurs]
  else:
    times = []


    h("""
      objref g
      proc plotraster() {
        g = new Graph()
        spikedCells.mark(g, spikes, "|")
      }

      proc myrun() {
        run()
        plotraster()
      }
      """)

  h("""
{epnm.pc.runworker()}
{epnm.pc.done()}
""")
  print("runworkers OK!")

  return [times,vSoma,spikes]


