/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*****************************************************************************
 * $Source: /home/cvs/namd/cvsroot/namd2/src/Controller.C,v $
 * $Author: jim $
 * $Date: 2015/12/17 21:22:03 $
 * $Revision: 1.1312 $
 *****************************************************************************/

#include "InfoStream.h"
#include "memusage.h"
#include "Node.h"
#include "Molecule.h"
#include "SimParameters.h"
#include "Controller.h"
#include "ReductionMgr.h"
#include "CollectionMaster.h"
#include "Output.h"
#include "strlib.h"
#include "BroadcastObject.h"
#include "NamdState.h"
#include "ScriptTcl.h"
#include "Broadcasts.h"
#include "LdbCoordinator.h"
#include "Thread.h"
#include <math.h>
#include <signal.h>
#include "NamdOneTools.h"
#include "PatchMap.h"
#include "PatchMap.inl"
#include "Random.h"
#include "imd.h"
#include "IMDOutput.h"
#include "BackEnd.h"
#include <fstream>
#include <iomanip>
#include <errno.h>
#include "qd.h"

#include "ComputeNonbondedMICKernel.h"

#if(CMK_CCS_AVAILABLE && CMK_WEB_MODE)
extern "C" void CApplicationDepositNode0Data(char *);
#endif

#ifndef cbrt
  // cbrt() not in math.h on goneril
  #define cbrt(x)  pow(x,(double)(1.0/3.0))
#endif

//#define DEBUG_PRESSURE
#define MIN_DEBUG_LEVEL 3
//#define DEBUGM
#include "Debug.h"

#define XXXBIGREAL 1.0e32

class PressureProfileReduction {
private:
  RequireReduction *reduction;
  const int nslabs;
  const int freq;
  int nelements;
  int count;
  char *name;

  BigReal *average;

public:
  PressureProfileReduction(int rtag, int numslabs, int numpartitions,
      const char *myname, int outputfreq)
  : nslabs(numslabs), freq(outputfreq) {
    name = strdup(myname);
    nelements = 3*nslabs * numpartitions;
    reduction = ReductionMgr::Object()->willRequire(rtag,nelements);

    average = new BigReal[nelements];
    count = 0;
  }
  ~PressureProfileReduction() {
    delete [] average;
    delete reduction;
    free(name);
  }
  // 
  void getData(int firsttimestep, int step, const Lattice &lattice, 
      BigReal *total) {
    reduction->require();

    int i;
    double inv_volume = 1.0 / lattice.volume();
    // accumulate the pressure profile computed for this step into the average.
    int arraysize = 3*nslabs;
    for (i=0; i<nelements; i++) {
      BigReal val = reduction->item(i) * inv_volume;
      average[i] += val;
      total[i % arraysize] += val;
    }
    count++;
    if (!(step % freq)) {
      // convert NAMD internal virial to pressure in units of bar 
      BigReal scalefac = PRESSUREFACTOR * nslabs;
 
      iout << "PPROFILE" << name << ": " << step << " ";
      if (step == firsttimestep) {
        // output pressure profile for this step
        for (i=0; i<nelements; i++) 
          iout << reduction->item(i)*scalefac*inv_volume << " ";
      } else {
        // output pressure profile averaged over the last Freq steps.
        scalefac /= count;
        for (i=0; i<nelements; i++) 
          iout << average[i]*scalefac << " ";
      }
      iout << "\n" << endi; 
      // Clear the average for the next block
      memset(average, 0, nelements*sizeof(BigReal));
      count = 0;
    }
  }
};


Controller::Controller(NamdState *s) :
	computeChecksum(0), marginViolations(0), pairlistWarnings(0),
	simParams(Node::Object()->simParameters),
	state(s),
	collection(CollectionMaster::Object()),
        startCTime(0),
        firstCTime(CmiTimer()),
        startWTime(0),
        firstWTime(CmiWallTimer()),
        startBenchTime(0),
        stepInFullRun(0),
	ldbSteps(0),
	fflush_count(3)
{
    broadcast = new ControllerBroadcasts;
    min_reduction = ReductionMgr::Object()->willRequire(REDUCTIONS_MINIMIZER,1);
    // for accelMD
    if (simParams->accelMDOn) {
       // REDUCTIONS_BASIC wil contain all potential energies and arrive first
       amd_reduction = ReductionMgr::Object()->willRequire(REDUCTIONS_BASIC);
       // contents of amd_reduction be submitted to REDUCTIONS_AMD
       submit_reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_AMD);
       // REDUCTIONS_AMD will contain Sequencer contributions and arrive second
       reduction = ReductionMgr::Object()->willRequire(REDUCTIONS_AMD);
    } else {
       amd_reduction = NULL;
       submit_reduction = NULL;
       reduction = ReductionMgr::Object()->willRequire(REDUCTIONS_BASIC);
    }
    // pressure profile reductions
    pressureProfileSlabs = 0;
    pressureProfileCount = 0;
    ppbonded = ppnonbonded = ppint = NULL;
    if (simParams->pressureProfileOn || simParams->pressureProfileEwaldOn) {
      int ntypes = simParams->pressureProfileAtomTypes;
      int nslabs = pressureProfileSlabs = simParams->pressureProfileSlabs;
      // number of partitions for pairwise interactions
      int npairs = (ntypes * (ntypes+1))/2;
      pressureProfileAverage = new BigReal[3*nslabs];
      memset(pressureProfileAverage, 0, 3*nslabs*sizeof(BigReal));
      int freq = simParams->pressureProfileFreq;
      if (simParams->pressureProfileOn) {
        ppbonded = new PressureProfileReduction(REDUCTIONS_PPROF_BONDED, 
            nslabs, npairs, "BONDED", freq);
        ppnonbonded = new PressureProfileReduction(REDUCTIONS_PPROF_NONBONDED, 
            nslabs, npairs, "NONBONDED", freq);
        // internal partitions by atom type, but not in a pairwise fashion
        ppint = new PressureProfileReduction(REDUCTIONS_PPROF_INTERNAL, 
            nslabs, ntypes, "INTERNAL", freq);
      } else {
        // just doing Ewald, so only need nonbonded
        ppnonbonded = new PressureProfileReduction(REDUCTIONS_PPROF_NONBONDED, 
            nslabs, npairs, "NONBONDED", freq);
      }
    }
    random = new Random(simParams->randomSeed);
    random->split(0,PatchMap::Object()->numPatches()+1);

    rescaleVelocities_sumTemps = 0;  rescaleVelocities_numTemps = 0;
    rescaleVelocities_count    = 0;
    rescaleVelocities_sbeta    = 0;
    rescaleVelocities_sum1     = 0;
    rescaleVelocities_sumBeta  = 0;
    rescaleVelocities_sumBeta2 = 0;
    rescaleVelocities_sumDbde  = 0;
    berendsenPressure_avg = 0; berendsenPressure_count = 0;
    // strainRate tensor is symmetric to avoid rotation
    langevinPiston_strainRate =
	Tensor::symmetric(simParams->strainRate,simParams->strainRate2);
    if ( ! simParams->useFlexibleCell ) {
      BigReal avg = trace(langevinPiston_strainRate) / 3.;
      langevinPiston_strainRate = Tensor::identity(avg);
    } else if ( simParams->useConstantRatio ) {
#define AVGXY(T) T.xy = T.yx = 0; T.xx = T.yy = 0.5 * ( T.xx + T.yy );\
		 T.xz = T.zx = T.yz = T.zy = 0.5 * ( T.xz + T.yz );
      AVGXY(langevinPiston_strainRate);
#undef AVGXY
    }
    langevinPiston_origStrainRate = langevinPiston_strainRate;
    origLattice = state->lattice;
    smooth2_avg = XXXBIGREAL;
    temp_avg = 0;
    pressure_avg = 0;
    groupPressure_avg = 0;
    avg_count = 0;
    pressure_tavg = 0;
    groupPressure_tavg = 0;
    tavg_count = 0;
    checkpoint_stored = 0;
    drudeBondTemp = 0;
    drudeBondTempAvg = 0;
}

Controller::~Controller(void)
{
    delete broadcast;
    delete reduction;
    delete min_reduction;
    delete amd_reduction;
    delete submit_reduction;
    delete ppbonded;
    delete ppnonbonded;
    delete ppint;
    delete [] pressureProfileAverage;
    delete random;
}

void Controller::threadRun(Controller* arg)
{
    arg->algorithm();
}

void Controller::run(void)
{
    // create a Thread and invoke it
    DebugM(4, "Starting thread in controller on this=" << this << "\n");
    thread = CthCreate((CthVoidFn)&(threadRun),(void*)(this),CTRL_STK_SZ);
    CthSetStrategyDefault(thread);
#if CMK_BLUEGENE_CHARM
    BgAttach(thread);
#endif
    awaken();
}


void Controller::algorithm(void)
{
  int scriptTask;
  int scriptSeq = 0;
  BackEnd::awaken();
  while ( (scriptTask = broadcast->scriptBarrier.get(scriptSeq++)) != SCRIPT_END) {
    switch ( scriptTask ) {
      case SCRIPT_OUTPUT:
        enqueueCollections(FILE_OUTPUT);
        outputExtendedSystem(FILE_OUTPUT);
        break;
      case SCRIPT_FORCEOUTPUT:
        enqueueCollections(FORCE_OUTPUT);
        break;
      case SCRIPT_MEASURE:
        enqueueCollections(EVAL_MEASURE);
        break;
      case SCRIPT_REINITVELS:
        iout << "REINITIALIZING VELOCITIES AT STEP " << simParams->firstTimestep
          << " TO " << simParams->initialTemp << " KELVIN.\n" << endi;
        break;
      case SCRIPT_RESCALEVELS:
        iout << "RESCALING VELOCITIES AT STEP " << simParams->firstTimestep
          << " BY " << simParams->scriptArg1 << "\n" << endi;
        break;
      case SCRIPT_CHECKPOINT:
        iout << "CHECKPOINTING AT STEP " << simParams->firstTimestep
          << "\n" << endi;
        checkpoint_stored = 1;
        checkpoint_lattice = state->lattice;
        checkpoint_state = *(ControllerState*)this;
        break;
      case SCRIPT_REVERT:
        iout << "REVERTING AT STEP " << simParams->firstTimestep
          << "\n" << endi;
        if ( ! checkpoint_stored )
          NAMD_die("Unable to revert, checkpoint was never called!");
        state->lattice = checkpoint_lattice;
        *(ControllerState*)this = checkpoint_state;
        break;
      case SCRIPT_CHECKPOINT_STORE:
        iout << "STORING CHECKPOINT AT STEP " << simParams->firstTimestep
          << " TO KEY " << simParams->scriptStringArg1;
        if ( CmiNumPartitions() > 1 ) iout << " ON REPLICA " << simParams->scriptIntArg1;
        iout << "\n" << endi;
        if ( simParams->scriptIntArg1 == CmiMyPartition() ) {
          if ( ! checkpoints.count(simParams->scriptStringArg1) ) {
            checkpoints[simParams->scriptStringArg1] = new checkpoint;
          }
          checkpoint &cp = *checkpoints[simParams->scriptStringArg1];
          cp.lattice = state->lattice;
          cp.state = *(ControllerState*)this;
        } else {
          Node::Object()->sendCheckpointReq(simParams->scriptIntArg1,simParams->scriptStringArg1,
                                            scriptTask, state->lattice, *(ControllerState*)this);
        }
        break;
      case SCRIPT_CHECKPOINT_LOAD:
        iout << "LOADING CHECKPOINT AT STEP " << simParams->firstTimestep
          << " FROM KEY " << simParams->scriptStringArg1;
        if ( CmiNumPartitions() > 1 ) iout << " ON REPLICA " << simParams->scriptIntArg1;
        iout << "\n" << endi;
        if ( simParams->scriptIntArg1 == CmiMyPartition() ) {
          if ( ! checkpoints.count(simParams->scriptStringArg1) ) {
            NAMD_die("Unable to load checkpoint, requested key was never stored.");
          }
          checkpoint &cp = *checkpoints[simParams->scriptStringArg1];
          state->lattice = cp.lattice;
          *(ControllerState*)this = cp.state;
        } else {
          Node::Object()->sendCheckpointReq(simParams->scriptIntArg1,simParams->scriptStringArg1,
                                            scriptTask, state->lattice, *(ControllerState*)this);
        }
        break;
      case SCRIPT_CHECKPOINT_SWAP:
        iout << "SWAPPING CHECKPOINT AT STEP " << simParams->firstTimestep
          << " FROM KEY " << simParams->scriptStringArg1;
        if ( CmiNumPartitions() > 1 ) iout << " ON REPLICA " << simParams->scriptIntArg1;
        iout << "\n" << endi;
        if ( simParams->scriptIntArg1 == CmiMyPartition() ) {
          if ( ! checkpoints.count(simParams->scriptStringArg1) ) {
            NAMD_die("Unable to swap checkpoint, requested key was never stored.");
          }
          checkpoint &cp = *checkpoints[simParams->scriptStringArg1];
          std::swap(state->lattice,cp.lattice);
          std::swap(*(ControllerState*)this,cp.state);
        } else {
          Node::Object()->sendCheckpointReq(simParams->scriptIntArg1,simParams->scriptStringArg1,
                                            scriptTask, state->lattice, *(ControllerState*)this);
        }
        break;
      case SCRIPT_CHECKPOINT_FREE:
        iout << "FREEING CHECKPOINT AT STEP " << simParams->firstTimestep
          << " FROM KEY " << simParams->scriptStringArg1;
        if ( CmiNumPartitions() > 1 ) iout << " ON REPLICA " << simParams->scriptIntArg1;
        iout << "\n" << endi;
        if ( simParams->scriptIntArg1 == CmiMyPartition() ) {
          if ( ! checkpoints.count(simParams->scriptStringArg1) ) {
            NAMD_die("Unable to free checkpoint, requested key was never stored.");
          }
          delete checkpoints[simParams->scriptStringArg1];
          checkpoints.erase(simParams->scriptStringArg1);
        } else {
          Node::Object()->sendCheckpointReq(simParams->scriptIntArg1,simParams->scriptStringArg1,
                                            scriptTask, state->lattice, *(ControllerState*)this);
        }
        break;
      case SCRIPT_ATOMSENDRECV:
      case SCRIPT_ATOMSEND:
      case SCRIPT_ATOMRECV:
        break;
      case SCRIPT_MINIMIZE:
        minimize();
        break;
      case SCRIPT_RUN:
      case SCRIPT_CONTINUE:
        integrate(scriptTask);
        break;
      default:
        NAMD_bug("Unknown task in Controller::algorithm");
    }
    BackEnd::awaken();
  }
  enqueueCollections(END_OF_RUN);
  outputExtendedSystem(END_OF_RUN);
  terminate();
}


//
// XXX static and global variables are unsafe for shared memory builds.
// The use of global and static vars should be eliminated.
//
extern int eventEndOfTimeStep;

// Handle SIGINT so that restart files get written completely.
static int gotsigint = 0;
static void my_sigint_handler(int sig) {
  if (sig == SIGINT) gotsigint = 1;
}
extern "C" {
  typedef void (*namd_sighandler_t)(int);
}

void Controller::integrate(int scriptTask) {
    char traceNote[24];
  
    int step = simParams->firstTimestep;

    const int numberOfSteps = simParams->N;
    const int stepsPerCycle = simParams->stepsPerCycle;

    const int zeroMomentum = simParams->zeroMomentum;

    nbondFreq = simParams->nonbondedFrequency;
    const int dofull = ( simParams->fullElectFrequency ? 1 : 0 );
    if (dofull)
      slowFreq = simParams->fullElectFrequency;
    else
      slowFreq = simParams->nonbondedFrequency;
    if ( step >= numberOfSteps ) slowFreq = nbondFreq = 1;

    langRescaleFactorPrev = 1;
    tNHCInit();
    rescaleVelocitiesInit();

  if ( scriptTask == SCRIPT_RUN ) {

    reassignVelocities(step);  // only for full-step velecities
    rescaleaccelMD(step);

    receivePressure(step);
    if ( zeroMomentum && dofull && ! (step % slowFreq) )
						correctMomentum(step);
    // they shouldn't be called
    //langRescaleVelocities(step, TRUE);
    //tNHCRescaleVelocities(step, TRUE);
    adaptTempUpdate(step); // Init adaptive tempering;

    printFepMessage(step);
    printTiMessage(step);
    printDynamicsEnergies(step);
    outputFepEnergy(step);
    outputTiEnergy(step);
    if(traceIsOn()){
        traceUserEvent(eventEndOfTimeStep);
        sprintf(traceNote, "s:%d", step);
        traceUserSuppliedNote(traceNote);
    }
    outputExtendedSystem(step);
    rebalanceLoad(step);

  }
    keHistInit();
    if ( simParams->rescaleInitTotal != 0 ) {
      rescaleForTotalEnergy();
    }
    iout << "Potential energy " << totalEnergy - kineticEnergy
         << ", kinetic energy " << kineticEnergy
         << ", total energy " << totalEnergy << "\n" << endi;
    std::ofstream fsEnergyLog;
    if ( simParams->energyLogFreq > 0 ) {
      fsEnergyLog.open(simParams->energyLogFile, std::ios_base::app);
      fsEnergyLog << step << " " << (totalEnergy - kineticEnergy);
      if ( simParams->energyLogTotal )
        fsEnergyLog << " " << totalEnergy;
      if ( simParams->adaptTempOn )
        fsEnergyLog << " " << adaptTempT;
      fsEnergyLog << std::endl;
    }

    // Handling SIGINT doesn't seem to be working on Lemieux, and it
    // sometimes causes the net-xxx versions of NAMD to segfault on exit, 
    // so disable it for now.
    // namd_sighandler_t oldhandler = signal(SIGINT, 
    //  (namd_sighandler_t)my_sigint_handler);
    for ( ++step ; step <= numberOfSteps; ++step )
    {
        rescaleVelocities(step);
	tcoupleVelocities(step);
	langRescaleVelocities(step, FALSE);
	tNHCRescaleVelocities(step, FALSE);
	berendsenPressure(step);
	langevinPiston1(step);
        rescaleaccelMD(step);
	enqueueCollections(step);  // after lattice scaling!
	receivePressure(step); // calling reduction->require()
        if ( zeroMomentum && dofull && ! (step % slowFreq) )
						correctMomentum(step);
	langevinPiston2(step);
        reassignVelocities(step);
	
        langRescaleVelocities(step, TRUE);
        tNHCRescaleVelocities(step, TRUE);
        Bool scaled = adaptTempUpdate(step);
        keHistUpdate(step);
        printDynamicsEnergies(step);
        if ( fsEnergyLog.is_open() && step % simParams->energyLogFreq == 0 ) {
          fsEnergyLog << step << " " << (totalEnergy - kineticEnergy);
          if ( simParams->energyLogTotal )
            fsEnergyLog << " " << totalEnergy;
          if ( simParams->adaptTempOn )
            fsEnergyLog << " " << adaptTempT;
          fsEnergyLog << std::endl;
        }
        outputFepEnergy(step);
        outputTiEnergy(step);
        if(traceIsOn()){
            traceUserEvent(eventEndOfTimeStep);
            sprintf(traceNote, "s:%d", step);
            traceUserSuppliedNote(traceNote);
        }
  // if (gotsigint) {
  //   iout << iINFO << "Received SIGINT; shutting down.\n" << endi;
  //   NAMD_quit();
  // }
        outputExtendedSystem(step);
#if CYCLE_BARRIER
        cycleBarrier(!((step+1) % stepsPerCycle),step);
#elif  PME_BARRIER
        cycleBarrier(dofull && !(step%slowFreq),step);
#elif  STEP_BARRIER
        cycleBarrier(1, step);
#endif
		
	if(Node::Object()->specialTracing || simParams->statsOn){		
		 int bstep = simParams->traceStartStep;
		 int estep = bstep + simParams->numTraceSteps;		
		 if(step == bstep){
			 traceBarrier(1, step);
		 }
		 if(step == estep){			
			 traceBarrier(0, step);
		 }
	 }

#ifdef MEASURE_NAMD_WITH_PAPI
	if(simParams->papiMeasure) {
		int bstep = simParams->papiMeasureStartStep;
		int estep = bstep + simParams->numPapiMeasureSteps;
		if(step == bstep) {
			papiMeasureBarrier(1, step);
		}
		if(step == estep) {
			papiMeasureBarrier(0, step);
		}
	}
#endif
        if ( scaled && (ldbSteps == 1 || LdbCoordinator::Object()->getNumStepsToRun() == 1) ) {
          // collect Hi's if we're about to rebalance load
          collection->enqueueHi(step);
        }
	//if ( ldbSteps == 1 ) CkPrintf("step %d, before rebalanceLoad(), Controller PE %d/%d, thread %p\n", step, CkMyPe(), CkNumPes(), CthSelf());
        rebalanceLoad(step);

#if  PME_BARRIER
        cycleBarrier(dofull && !((step+1)%slowFreq),step);   // step before PME
#endif
    }
    // signal(SIGINT, oldhandler);
    
    adaptTempDone(step);
    rescaleVelocitiesSave(step);
    if ( fsEnergyLog.is_open() ) {
      fsEnergyLog << std::endl;
      fsEnergyLog.close();
    }
    tNHCDone(step);
    keHistDone(step);
}


#define CALCULATE \
  printMinimizeEnergies(step); \
  outputExtendedSystem(step); \
  rebalanceLoad(step); \
  if ( step == numberOfSteps ) return; \
  else ++step;

#define MOVETO(X) \
  if ( step == numberOfSteps ) { \
    if ( minVerbose ) { iout << "LINE MINIMIZER: RETURNING TO " << mid.x << " FROM " << last.x << "\n" << endi; } \
    if ( newDir || (mid.x-last.x) ) { \
      broadcast->minimizeCoefficient.publish(minSeq++,mid.x-last.x); \
    } else { \
      broadcast->minimizeCoefficient.publish(minSeq++,0.); \
      broadcast->minimizeCoefficient.publish(minSeq++,0.); \
      min_reduction->require(); \
      broadcast->minimizeCoefficient.publish(minSeq++,0.); \
    } \
    enqueueCollections(step); \
    CALCULATE \
  } else if ( (X)-last.x ) { \
    broadcast->minimizeCoefficient.publish(minSeq++,(X)-last.x); \
    newDir = 0; \
    last.x = (X); \
    enqueueCollections(step); \
    CALCULATE \
    last.u = min_energy; \
    last.dudx = -1. * min_f_dot_v; \
    last.noGradient = min_huge_count; \
    if ( minVerbose ) { \
      iout << "LINE MINIMIZER: POSITION " << last.x << " ENERGY " << last.u << " GRADIENT " << last.dudx; \
      if ( last.noGradient ) iout << " HUGECOUNT " << last.noGradient; \
      iout << "\n" << endi; \
    } \
  }

struct minpoint {
  BigReal x, u, dudx, noGradient;
  minpoint() : x(0), u(0), dudx(0), noGradient(1) { ; }
};

void Controller::minimize() {
  // iout << "Controller::minimize() called.\n" << endi;

  const int minVerbose = simParams->minVerbose;
  const int numberOfSteps = simParams->N;
  int step = simParams->firstTimestep;
  slowFreq = nbondFreq = 1;
  BigReal linegoal = simParams->minLineGoal;  // 1.0e-3
  const BigReal goldenRatio = 0.5 * ( sqrt(5.0) - 1.0 );

  CALCULATE

  int minSeq = 0;

  // just move downhill until initial bad contacts go away or 100 steps
  int old_min_huge_count = 2000000000;
  int steps_at_same_min_huge_count = 0;
  for ( int i_slow = 0; min_huge_count && i_slow < 100; ++i_slow ) {
    if ( min_huge_count >= old_min_huge_count ) {
      if ( ++steps_at_same_min_huge_count > 10 ) break;
    } else {
      old_min_huge_count = min_huge_count;
      steps_at_same_min_huge_count = 0;
    }
    iout << "MINIMIZER SLOWLY MOVING " << min_huge_count << " ATOMS WITH BAD CONTACTS DOWNHILL\n" << endi;
    broadcast->minimizeCoefficient.publish(minSeq++,1.);
    enqueueCollections(step);
    CALCULATE
  }
  if ( min_huge_count ) {
    iout << "MINIMIZER GIVING UP ON " << min_huge_count << " ATOMS WITH BAD CONTACTS\n" << endi;
  }
  iout << "MINIMIZER STARTING CONJUGATE GRADIENT ALGORITHM\n" << endi;

  int atStart = 2;
  int errorFactor = 10;
  BigReal old_f_dot_f = min_f_dot_f;
  broadcast->minimizeCoefficient.publish(minSeq++,0.);
  broadcast->minimizeCoefficient.publish(minSeq++,0.); // v = f
  int newDir = 1;
  min_f_dot_v = min_f_dot_f;
  min_v_dot_v = min_f_dot_f;
  while ( 1 ) {
    // line minimization
    // bracket minimum on line
    minpoint lo,hi,mid,last;
    BigReal x = 0;
    lo.x = x;
    lo.u = min_energy;
    lo.dudx = -1. * min_f_dot_v;
    lo.noGradient = min_huge_count;
    mid = lo;
    last = mid;
    if ( minVerbose ) {
      iout << "LINE MINIMIZER: POSITION " << last.x << " ENERGY " << last.u << " GRADIENT " << last.dudx;
      if ( last.noGradient ) iout << " HUGECOUNT " << last.noGradient;
      iout << "\n" << endi;
    }
    BigReal tol = fabs( linegoal * min_f_dot_v );
    iout << "LINE MINIMIZER REDUCING GRADIENT FROM " <<
            fabs(min_f_dot_v) << " TO " << tol << "\n" << endi;
    int start_with_huge = last.noGradient;
    min_reduction->require();
    BigReal maxstep = 0.1 / sqrt(min_reduction->item(0));
    x = maxstep; MOVETO(x);
    // bracket minimum on line
    while ( last.u < mid.u ) {
      if ( last.dudx < mid.dudx * (5.5 - x/maxstep)/5. ) {
        x = 6 * maxstep; break;
      }
      lo = mid; mid = last;
      x += maxstep;
      if ( x > 5.5 * maxstep ) break;
      MOVETO(x)
    }
    if ( x > 5.5 * maxstep ) {
      iout << "MINIMIZER RESTARTING CONJUGATE GRADIENT ALGORITHM DUE TO POOR PROGRESS\n" << endi;
      broadcast->minimizeCoefficient.publish(minSeq++,0.);
      broadcast->minimizeCoefficient.publish(minSeq++,0.); // v = f
      newDir = 1;
      old_f_dot_f = min_f_dot_f;
      min_f_dot_v = min_f_dot_f;
      min_v_dot_v = min_f_dot_f;
      continue;
    }
    hi = last;
#define PRINT_BRACKET \
    iout << "LINE MINIMIZER BRACKET: DX " \
         << (mid.x-lo.x) << " " << (hi.x-mid.x) << \
        " DU " << (mid.u-lo.u) << " " << (hi.u-mid.u) << " DUDX " << \
        lo.dudx << " " << mid.dudx << " " << hi.dudx << " \n" << endi;
    PRINT_BRACKET
    // converge on minimum on line
    int itcnt;
    for ( itcnt = 10; itcnt > 0; --itcnt ) {
      int progress = 1;
      // select new position
      if ( mid.noGradient ) {
       if ( ( mid.x - lo.x ) > ( hi.x - mid.x ) ) {  // subdivide left side
	x = (1.0 - goldenRatio) * lo.x + goldenRatio * mid.x;
	MOVETO(x)
	if ( last.u <= mid.u ) {
	  hi = mid; mid = last;
	} else {
	  lo = last;
	}
       } else {  // subdivide right side
	x = (1.0 - goldenRatio) * hi.x + goldenRatio * mid.x;
	MOVETO(x)
	if ( last.u <= mid.u ) {
	  lo = mid; mid = last;
	} else {
	  hi = last;
	}
       }
      } else {
       if ( mid.dudx > 0. ) {  // subdivide left side
        BigReal altxhi = 0.1 * lo.x + 0.9 * mid.x;
        BigReal altxlo = 0.9 * lo.x + 0.1 * mid.x;
        x = mid.dudx*(mid.x*mid.x-lo.x*lo.x) + 2*mid.x*(lo.u-mid.u);
        BigReal xdiv = 2*(mid.dudx*(mid.x-lo.x)+(lo.u-mid.u));
        if ( xdiv ) x /= xdiv; else x = last.x;
        if ( x > altxhi ) x = altxhi;
        if ( x < altxlo ) x = altxlo;
        if ( x-last.x == 0 ) break;
        MOVETO(x)
        if ( last.u <= mid.u ) {
	  hi = mid; mid = last;
	} else {
          if ( lo.dudx < 0. && last.dudx > 0. ) progress = 0;
	  lo = last;
	}
       } else {  // subdivide right side
        BigReal altxlo = 0.1 * hi.x + 0.9 * mid.x;
        BigReal altxhi = 0.9 * hi.x + 0.1 * mid.x;
        x = mid.dudx*(mid.x*mid.x-hi.x*hi.x) + 2*mid.x*(hi.u-mid.u);
        BigReal xdiv = 2*(mid.dudx*(mid.x-hi.x)+(hi.u-mid.u));
        if ( xdiv ) x /= xdiv; else x = last.x;
        if ( x < altxlo ) x = altxlo;
        if ( x > altxhi ) x = altxhi;
        if ( x-last.x == 0 ) break;
        MOVETO(x)
        if ( last.u <= mid.u ) {
	  lo = mid; mid = last;
	} else {
          if ( hi.dudx > 0. && last.dudx < 0. ) progress = 0;
	  hi = last;
	}
       }
      }
      PRINT_BRACKET
      if ( ! progress ) {
        MOVETO(mid.x);
        break;
      }
      if ( fabs(last.dudx) < tol ) break;
      if ( lo.x != mid.x && (lo.u-mid.u) < tol * (mid.x-lo.x) ) break;
      if ( mid.x != hi.x && (hi.u-mid.u) < tol * (hi.x-mid.x) ) break;
    }
    // new direction
    broadcast->minimizeCoefficient.publish(minSeq++,0.);
    BigReal c = min_f_dot_f / old_f_dot_f;
    c = ( c > 1.5 ? 1.5 : c );
    if ( atStart ) { c = 0; --atStart; }
    if ( c*c*min_v_dot_v > errorFactor*min_f_dot_f ) {
      c = 0;
      if ( errorFactor < 100 ) errorFactor += 10;
    }
    if ( c == 0 ) {
      iout << "MINIMIZER RESTARTING CONJUGATE GRADIENT ALGORITHM\n" << endi;
    }
    broadcast->minimizeCoefficient.publish(minSeq++,c); // v = c*v+f
    newDir = 1;
    old_f_dot_f = min_f_dot_f;
    min_f_dot_v = c * min_f_dot_v + min_f_dot_f;
    min_v_dot_v = c*c*min_v_dot_v + 2*c*min_f_dot_v + min_f_dot_f;
  }

}

#undef MOVETO
#undef CALCULATE

void Controller::berendsenPressure(int step)
{
  if ( simParams->berendsenPressureOn ) {
   berendsenPressure_count += 1;
   berendsenPressure_avg += controlPressure;
   const int freq = simParams->berendsenPressureFreq;
   if ( ! (berendsenPressure_count % freq) ) {
    Tensor factor = berendsenPressure_avg / berendsenPressure_count;
    berendsenPressure_avg = 0;
    berendsenPressure_count = 0;
    // We only use on-diagonal terms (for now)
    factor = Tensor::diagonal(diagonal(factor));
    factor -= Tensor::identity(simParams->berendsenPressureTarget);
    factor *= ( ( simParams->berendsenPressureCompressibility / 3.0 ) *
       simParams->dt * freq / simParams->berendsenPressureRelaxationTime );
    factor += Tensor::identity(1.0);
#define LIMIT_SCALING(VAR,MIN,MAX,FLAG) {\
         if ( VAR < (MIN) ) { VAR = (MIN); FLAG = 1; } \
         if ( VAR > (MAX) ) { VAR = (MAX); FLAG = 1; } }
    int limited = 0;
    LIMIT_SCALING(factor.xx,1./1.03,1.03,limited)
    LIMIT_SCALING(factor.yy,1./1.03,1.03,limited)
    LIMIT_SCALING(factor.zz,1./1.03,1.03,limited)
#undef LIMIT_SCALING
    if ( limited ) {
      iout << iERROR << "Step " << step <<
	" cell rescaling factor limited.\n" << endi;
    }
    broadcast->positionRescaleFactor.publish(step,factor);
    state->lattice.rescale(factor);
   }
  } else {
    berendsenPressure_avg = 0;
    berendsenPressure_count = 0;
  }
}

void Controller::langevinPiston1(int step)
{
  if ( simParams->langevinPistonOn )
  {
    Tensor &strainRate = langevinPiston_strainRate;
    int cellDims = simParams->useFlexibleCell ? 1 : 3;
    BigReal dt = simParams->dt;
    BigReal dt_long = slowFreq * dt;
    BigReal kT = BOLTZMANN * simParams->langevinPistonTemp;
    BigReal tau = simParams->langevinPistonPeriod;
    BigReal mass = controlNumDegFreedom * kT * tau * tau * cellDims;

#ifdef DEBUG_PRESSURE
    iout << iINFO << "entering langevinPiston1, strain rate: " << strainRate << "\n";
#endif

    if ( ! ( (step-1) % slowFreq ) )
    {
      BigReal gamma = 1 / simParams->langevinPistonDecay;
      BigReal f1 = exp( -0.5 * dt_long * gamma );
      BigReal f2 = sqrt( ( 1. - f1*f1 ) * kT / mass );
      strainRate *= f1;
      if ( simParams->useFlexibleCell ) {
        // We only use on-diagonal terms (for now)
        if ( simParams->useConstantRatio ) {
	  BigReal r = f2 * random->gaussian();
	  strainRate.xx += r;
	  strainRate.yy += r;
	  strainRate.zz += f2 * random->gaussian();
        } else {
	  strainRate += f2 * Tensor::diagonal(random->gaussian_vector());
        }
      } else {
	strainRate += f2 * Tensor::identity(random->gaussian());
      }

#ifdef DEBUG_PRESSURE
      iout << iINFO << "applying langevin, strain rate: " << strainRate << "\n";
#endif
    }

    // Apply surface tension.  If surfaceTensionTarget is zero, we get
    // the default (isotropic pressure) case.
    
    Tensor ptarget;
    ptarget.xx = ptarget.yy = ptarget.zz = simParams->langevinPistonTarget;
    if ( simParams->surfaceTensionTarget != 0. ) {
      ptarget.xx = ptarget.yy = simParams->langevinPistonTarget - 
        simParams->surfaceTensionTarget / state->lattice.c().z;
    }

    strainRate += ( 0.5 * dt * cellDims * state->lattice.volume() / mass ) *
      ( controlPressure - ptarget );

    if (simParams->fixCellDims) {
      if (simParams->fixCellDimX) strainRate.xx = 0;
      if (simParams->fixCellDimY) strainRate.yy = 0;
      if (simParams->fixCellDimZ) strainRate.zz = 0;
    }


#ifdef DEBUG_PRESSURE
    iout << iINFO << "integrating half step, strain rate: " << strainRate << "\n";
#endif

    if ( simParams->langevinPistonBarrier ) {
    if ( ! ( (step-1-slowFreq/2) % slowFreq ) )
    {
      // We only use on-diagonal terms (for now)
      Tensor factor;
      if ( !simParams->useConstantArea ) {
        factor.xx = exp( dt_long * strainRate.xx );
        factor.yy = exp( dt_long * strainRate.yy );
      } else {
        factor.xx = factor.yy = 1;
      }
      factor.zz = exp( dt_long * strainRate.zz );
      broadcast->positionRescaleFactor.publish(step,factor);
      state->lattice.rescale(factor);
#ifdef DEBUG_PRESSURE
      iout << iINFO << "rescaling by: " << factor << "\n";
#endif
    }
    } else { // langevinPistonBarrier
    if ( ! ( (step-1-slowFreq/2) % slowFreq ) )
    {
      if ( ! positionRescaleFactor.xx ) {  // first pass without MTS
      // We only use on-diagonal terms (for now)
      Tensor factor;
      if ( !simParams->useConstantArea ) {
        factor.xx = exp( dt_long * strainRate.xx );
        factor.yy = exp( dt_long * strainRate.yy );
      } else {
        factor.xx = factor.yy = 1;
      }
      factor.zz = exp( dt_long * strainRate.zz );
      broadcast->positionRescaleFactor.publish(step,factor);
      positionRescaleFactor = factor;
      strainRate_old = strainRate;
      }
      state->lattice.rescale(positionRescaleFactor);
#ifdef DEBUG_PRESSURE
      iout << iINFO << "rescaling by: " << factor << "\n";
#endif
    }
    if ( ! ( (step-slowFreq/2) % slowFreq ) )
    {
      // We only use on-diagonal terms (for now)
      Tensor factor;
      if ( !simParams->useConstantArea ) {
        factor.xx = exp( (dt+dt_long) * strainRate.xx - dt * strainRate_old.xx );
        factor.yy = exp( (dt+dt_long) * strainRate.yy - dt * strainRate_old.yy );
      } else {
        factor.xx = factor.yy = 1;
      }
      factor.zz = exp( (dt+dt_long) * strainRate.zz - dt * strainRate_old.zz );
      broadcast->positionRescaleFactor.publish(step+1,factor);
      positionRescaleFactor = factor;
      strainRate_old = strainRate;
    }
    }

    // corrections to integrator
    if ( ! ( step % nbondFreq ) )
    {
#ifdef DEBUG_PRESSURE
      iout << iINFO << "correcting strain rate for nbond, ";
#endif
      strainRate -= ( 0.5 * dt * cellDims * state->lattice.volume() / mass ) *
		( (nbondFreq - 1) * controlPressure_nbond );
#ifdef DEBUG_PRESSURE
      iout << "strain rate: " << strainRate << "\n";
#endif
    }
    if ( ! ( step % slowFreq ) )
    {
#ifdef DEBUG_PRESSURE
      iout << iINFO << "correcting strain rate for slow, ";
#endif
      strainRate -= ( 0.5 * dt * cellDims * state->lattice.volume() / mass ) *
		( (slowFreq - 1) * controlPressure_slow );
#ifdef DEBUG_PRESSURE
      iout << "strain rate: " << strainRate << "\n";
#endif
    }
    if (simParams->fixCellDims) {
      if (simParams->fixCellDimX) strainRate.xx = 0;
      if (simParams->fixCellDimY) strainRate.yy = 0;
      if (simParams->fixCellDimZ) strainRate.zz = 0;
    }

  }

}

void Controller::langevinPiston2(int step)
{
  if ( simParams->langevinPistonOn )
  {
    Tensor &strainRate = langevinPiston_strainRate;
    int cellDims = simParams->useFlexibleCell ? 1 : 3;
    BigReal dt = simParams->dt;
    BigReal dt_long = slowFreq * dt;
    BigReal kT = BOLTZMANN * simParams->langevinPistonTemp;
    BigReal tau = simParams->langevinPistonPeriod;
    BigReal mass = controlNumDegFreedom * kT * tau * tau * cellDims;

    // corrections to integrator
    if ( ! ( step % nbondFreq ) )
    {
#ifdef DEBUG_PRESSURE
      iout << iINFO << "correcting strain rate for nbond, ";
#endif
      strainRate += ( 0.5 * dt * cellDims * state->lattice.volume() / mass ) *
		( (nbondFreq - 1) * controlPressure_nbond );
#ifdef DEBUG_PRESSURE
      iout << "strain rate: " << strainRate << "\n";
#endif
    }
    if ( ! ( step % slowFreq ) )
    {
#ifdef DEBUG_PRESSURE
      iout << iINFO << "correcting strain rate for slow, ";
#endif
      strainRate += ( 0.5 * dt * cellDims * state->lattice.volume() / mass ) *
		( (slowFreq - 1) * controlPressure_slow );
#ifdef DEBUG_PRESSURE
      iout << "strain rate: " << strainRate << "\n";
#endif
    }

    // Apply surface tension.  If surfaceTensionTarget is zero, we get
    // the default (isotropic pressure) case.
   
    Tensor ptarget;
    ptarget.xx = ptarget.yy = ptarget.zz = simParams->langevinPistonTarget;
    if ( simParams->surfaceTensionTarget != 0. ) {
      ptarget.xx = ptarget.yy = simParams->langevinPistonTarget - 
        simParams->surfaceTensionTarget / state->lattice.c().z;
    }

    strainRate += ( 0.5 * dt * cellDims * state->lattice.volume() / mass ) *
      ( controlPressure - ptarget );

 
#ifdef DEBUG_PRESSURE
    iout << iINFO << "integrating half step, strain rate: " << strainRate << "\n";
#endif

    if ( ! ( step % slowFreq ) )
    {
      BigReal gamma = 1 / simParams->langevinPistonDecay;
      BigReal f1 = exp( -0.5 * dt_long * gamma );
      BigReal f2 = sqrt( ( 1. - f1*f1 ) * kT / mass );
      strainRate *= f1;
      if ( simParams->useFlexibleCell ) {
        // We only use on-diagonal terms (for now)
        if ( simParams->useConstantRatio ) {
	  BigReal r = f2 * random->gaussian();
	  strainRate.xx += r;
	  strainRate.yy += r;
	  strainRate.zz += f2 * random->gaussian();
        } else {
	  strainRate += f2 * Tensor::diagonal(random->gaussian_vector());
        }
      } else {
	strainRate += f2 * Tensor::identity(random->gaussian());
      }
#ifdef DEBUG_PRESSURE
      iout << iINFO << "applying langevin, strain rate: " << strainRate << "\n";
#endif
    }

#ifdef DEBUG_PRESSURE
    iout << iINFO << "exiting langevinPiston2, strain rate: " << strainRate << "\n";
#endif
    if (simParams->fixCellDims) {
      if (simParams->fixCellDimX) strainRate.xx = 0;
      if (simParams->fixCellDimY) strainRate.yy = 0;
      if (simParams->fixCellDimZ) strainRate.zz = 0;
    }
  }


}

void Controller::rescaleVelocities(int step)
{
  const int rescaleFreq = simParams->rescaleFreq;
  if ( rescaleFreq > 0 ) {
    rescaleVelocities_sumTemps += temperature;  ++rescaleVelocities_numTemps;
    BigReal ek = 0.5 * BOLTZMANN * temperature * numDegFreedom;
    BigReal beta = (0.5 * numDegFreedom - 1) / ek;
    // the following quantities are accumulated to compute the block average of beta
    rescaleVelocities_count += 1;
    rescaleVelocities_sbeta += beta;
    // the following quantities are accumulated to compute beta'(E)
    rescaleVelocities_sum1     += 1;
    rescaleVelocities_sumBeta  += beta;
    rescaleVelocities_sumBeta2 += beta * beta;
    rescaleVelocities_sumDbde  += -beta / ek;
    if ( rescaleVelocities_numTemps == rescaleFreq ) {
      BigReal avgTemp = rescaleVelocities_sumTemps / rescaleVelocities_numTemps;
      BigReal rescaleTemp = simParams->rescaleTemp;
      if ( simParams->adaptTempOn && simParams->adaptTempRescale && (step > simParams->adaptTempFirstStep ) && 
                (!(simParams->adaptTempLastStep > 0) || step < simParams->adaptTempLastStep )) {
        rescaleTemp = adaptTempT;
      }
      BigReal factor = sqrt(rescaleTemp/avgTemp);
      if ( simParams->rescaleAdaptiveOn ) {
        // recompute the velocity-rescaling factor
        BigReal bref = 1.0 / (BOLTZMANN * rescaleTemp);
        BigReal dbdk = rescaleVelocities_sumDbde / rescaleVelocities_sum1;
        BigReal bet = rescaleVelocities_sumBeta / rescaleVelocities_sum1;
        BigReal bvar = rescaleVelocities_sumBeta2 / rescaleVelocities_sum1 - bet * bet;
        BigReal dbde;
        if ( simParams->rescaleAdaptiveDKdE > 0 ) { // heuristic method
          dbde = -simParams->rescaleAdaptiveDKdE * bref * bref
               / (0.5 * numDegFreedom + simParams->rescaleAdaptiveDKdE - 2);
        } else { // exact method
          if ( rescaleVelocities_sum1 < 10 ) {
            dbde = dbdk;
          } else {
            dbde = dbdk + bvar;
            //CkPrintf("dbdk %g, dbde %g, bvar %g, threshold %g, count %g\n",
            //    dbdk, dbde, bvar, simParams->rescaleAdaptiveDKdEMin * dbdk, rescaleVelocities_count);
            if ( dbde > simParams->rescaleAdaptiveDKdEMin * dbdk )
              dbde = simParams->rescaleAdaptiveDKdEMin * dbdk;
          }
        }
        BigReal bave = rescaleVelocities_sbeta / rescaleVelocities_count;
        BigReal dbeta = bref - bave;
        BigReal de = dbeta / dbde;
        BigReal s;
        if ( simParams->rescaleAdaptiveMag > 0 ) {
          s = simParams->rescaleAdaptiveMag * (de / ek);
        } else {
          s = simParams->rescaleAdaptiveZoom * (de / ek)
            * simParams->rescaleFreq / rescaleVelocities_sum1;
        }
        if ( s > 0.5 ) s = 0.5;
        else if ( s < -0.5 ) s = -0.5;
        factor = sqrt(1 + s);
        if ( fmod(rescaleVelocities_sum1, simParams->rescaleAdaptiveFileFreq) < 0.5 ) {
          CkPrintf("step %d, factor %g, s %g, dbeta %g, bet %g/%g, dbde %g/%g = %5.3f, bvar %g, delE/K %g, tp %g, dof %ld\n",
              step, factor, s, dbeta, bet, bref, dbde, dbdk, dbde/dbdk, bvar, de/ek, temperature, (long) numDegFreedom); // getchar();
          rescaleVelocitiesSave(step);
        }
      }
      broadcast->velocityRescaleFactor.publish(step,factor);
      //iout << "RESCALING VELOCITIES AT STEP " << step
      //     << " FROM AVERAGE TEMPERATURE OF " << avgTemp
      //     << " TO " << rescaleTemp << " KELVIN.\n" << endi;
      // reset accumulators for block averages
      rescaleVelocities_sumTemps = 0;  rescaleVelocities_numTemps = 0;
      rescaleVelocities_count = 0;
      rescaleVelocities_sbeta = 0;
    }
  }
}

void Controller::rescaleVelocitiesInit(void)
{
  if ( numDegFreedom <= 0 )
    numDegFreedom = Node::Object()->molecule->num_deg_freedom();
  if ( simParams->rescaleFreq <= 0 || !simParams->rescaleAdaptiveOn )
    return;
  rescaleVelocitiesLoad();
}

void Controller::rescaleVelocitiesLoad(void)
{
  if ( simParams->rescaleFreq <= 0 || !simParams->rescaleAdaptiveOn )
    return;
  std::ifstream fs(simParams->rescaleAdaptiveFile);
  if ( !fs ) return;
  double cnt, beta, beta2, dbde;
  fs >> cnt >> beta >> beta2 >> dbde;
  fs.close();
  iout << "LOADED ADAPTIVE VELOCITY-RESCALING DATA FROM "
       << simParams->rescaleAdaptiveFile << ": STEP " << cnt
       << ", BETA " << beta << ", VAR(BETA) " << beta2
       << ", BETA'(E) " << dbde << "\n" << endi;
  beta2 += beta * beta;
  if ( fmod(cnt, simParams->rescaleFreq) > 0.5 ) {
    cnt = cnt / simParams->rescaleFreq * simParams->rescaleFreq;
  }
  rescaleVelocities_sum1 = cnt;
  rescaleVelocities_sumBeta = cnt * beta;
  rescaleVelocities_sumBeta2 = cnt * beta2;
  rescaleVelocities_sumDbde = cnt * dbde;
}

void Controller::rescaleVelocitiesSave(int step)
{
  if ( simParams->rescaleFreq <= 0 || !simParams->rescaleAdaptiveOn )
    return;
  if ( rescaleVelocities_sum1 <= 0 )
    return;
  NAMD_backup_file(simParams->rescaleAdaptiveFile);
  ofstream_namd fs(simParams->rescaleAdaptiveFile);
  if ( !fs ) return;
  BigReal beta = rescaleVelocities_sumBeta / rescaleVelocities_sum1;
  BigReal beta2 = rescaleVelocities_sumBeta2 / rescaleVelocities_sum1
                - beta * beta;
  BigReal dbde = rescaleVelocities_sumDbde / rescaleVelocities_sum1;
  char buf[200];
  sprintf(buf, "%.0f %.10f %.10f %.10f\n", rescaleVelocities_sum1, beta, beta2, dbde);
  fs << buf;
}

// scale velocities to achieve a pre-specified total energy
void Controller::rescaleForTotalEnergy()
{
  BigReal pe = totalEnergy - kineticEnergy;
  BigReal target = simParams->rescaleInitTotal
                 + random->gaussian() * simParams->rescaleInitDev;
  BigReal ke = target - pe;
  BigReal factor = (ke > 0) ? sqrt(ke/kineticEnergy) : 1.0;
  CkPrintf("Info: changing the initial KE from %g to %g, target total energy %g\n",
      kineticEnergy, ke, target);
  broadcast->velocityRescaleFactor.publish(0, factor);
  kineticEnergy = ke;
  totalEnergy = ke + pe;
}

void Controller::correctMomentum(int step) {

    Vector momentum;
    momentum.x = reduction->item(REDUCTION_HALFSTEP_MOMENTUM_X);
    momentum.y = reduction->item(REDUCTION_HALFSTEP_MOMENTUM_Y);
    momentum.z = reduction->item(REDUCTION_HALFSTEP_MOMENTUM_Z);
    const BigReal mass = reduction->item(REDUCTION_MOMENTUM_MASS);

    if ( momentum.length2() == 0. )
      iout << iERROR << "Momentum is exactly zero; probably a bug.\n" << endi;
    if ( mass == 0. )
      NAMD_die("Total mass is zero in Controller::correctMomentum");

    momentum *= (-1./mass);

    broadcast->momentumCorrection.publish(step+slowFreq,momentum);
}

//Modifications for alchemical fep
void Controller::printFepMessage(int step)
{
  if (simParams->alchFepOn) {
    const BigReal alchLambda = simParams->alchLambda;
    const BigReal alchLambda2 = simParams->alchLambda2;
    const BigReal alchTemp = simParams->alchTemp;
    const int alchEquilSteps = simParams->alchEquilSteps;
    iout << "FEP: RESETTING FOR NEW FEP WINDOW "
         << "LAMBDA SET TO " << alchLambda << " LAMBDA2 " << alchLambda2
         << "\nFEP: WINDOW TO HAVE " << alchEquilSteps
         << " STEPS OF EQUILIBRATION PRIOR TO FEP DATA COLLECTION.\n"
         << "FEP: USING CONSTANT TEMPERATURE OF " << alchTemp 
         << " K FOR FEP CALCULATION\n" << endi;
  }
} 
void Controller::printTiMessage(int step)
{
  if (simParams->alchThermIntOn) {
    const BigReal alchLambda = simParams->alchLambda;
    const BigReal alchTemp = simParams->alchTemp;
    const int alchEquilSteps = simParams->alchEquilSteps;
    iout << "TI: RESETTING FOR NEW WINDOW "
         << "LAMBDA SET TO " << alchLambda 
         << "\nTI: WINDOW TO HAVE " << alchEquilSteps 
         << " STEPS OF EQUILIBRATION PRIOR TO TI DATA COLLECTION.\n" << endi;
  }
} 
//fepe

void Controller::reassignVelocities(int step)
{
  const int reassignFreq = simParams->reassignFreq;
  if ( ( reassignFreq > 0 ) && ! ( step % reassignFreq ) ) {
    BigReal newTemp = simParams->reassignTemp;
    newTemp += ( step / reassignFreq ) * simParams->reassignIncr;
    if ( simParams->reassignIncr > 0.0 ) {
      if ( newTemp > simParams->reassignHold && simParams->reassignHold > 0.0 )
        newTemp = simParams->reassignHold;
    } else {
      if ( newTemp < simParams->reassignHold )
        newTemp = simParams->reassignHold;
    }
    iout << "REASSIGNING VELOCITIES AT STEP " << step
         << " TO " << newTemp << " KELVIN.\n" << endi;
  }
}

void Controller::tcoupleVelocities(int step)
{
  if ( simParams->tCoupleOn )
  {
    const BigReal tCoupleTemp = simParams->tCoupleTemp;
    BigReal coefficient = 1.;
    if ( temperature > 0. ) coefficient = tCoupleTemp/temperature - 1.;
    broadcast->tcoupleCoefficient.publish(step,coefficient);
  }
}

// Ref.: Canonical sampling through velocity rescaling
// Bussi, Donadio, and Parrinello, JCP 126, 014101 (2007)
void Controller::langRescaleVelocities(int step, Bool isPrev)
{
  if ( simParams->langRescaleOn ) {
    BigReal tp = simParams->langRescaleTemp;
    // use the temperature from adaptive tempering, if any
    if ( simParams->adaptTempOn && simParams->adaptTempRescale
     && (step > simParams->adaptTempFirstStep )
     && (!(simParams->adaptTempLastStep > 0) || step < simParams->adaptTempLastStep )) {
      tp = adaptTempT;
    }
    tp *= BOLTZMANN;

    BigReal dt = simParams->dt / simParams->langRescaleDt;
    dt *= 0.5; // doing a half time step, called twice in an MD step
    BigReal c = exp(-dt);

    // integrate two half steps
    int dof = numDegFreedom;
    if ( dof <= 0 ) dof = Node::Object()->molecule->num_deg_freedom();
    BigReal r = random->gaussian();
    BigReal r2 = random->chisqr(dof - 1);
    BigReal ek1 = BOLTZMANN * temperature * dof * 0.5;
    BigReal ek2 = ek1 + (1 - c) * ((r2 + r * r) * 0.5 * tp - ek1)
                + 2 * r * sqrt(c * (1 - c) * ek1 * 0.5 * tp);
    if ( ek2 < 0 ) ek2 = 0;
    BigReal fac2 = ek2 / ek1;
    BigReal factor = sqrt( fac2 );
    if ( !isPrev ) {
      // combine with the rescaling factor carried from
      // the second half-step from the previous MD step
      broadcast->langRescaleFactor.publish(step, factor * langRescaleFactorPrev);
    } else {
      // save it for the rescaling in the next step
      langRescaleFactorPrev = factor;
    }
    temperature *= fac2;
    totalEnergy += kineticEnergy * (fac2 - 1);
    kineticEnergy *= fac2;
    kineticEnergyCentered *= fac2;
  }
}

void Controller::tNHCInit(void)
{
  if ( !simParams->tNHCOn ) return;

  int nnhc = simParams->tNHCLen, i;
  tNHCzeta = new BigReal [nnhc];
  tNHCmass = new BigReal [nnhc];

  // Note: numDegFreedom has not been set yet
  int dof = Node::Object()->molecule->num_deg_freedom();
  BigReal per = simParams->tNHCPeriod / (2 * M_PI);
  BigReal kT = BOLTZMANN * simParams->tNHCTemp;
  // reference mass choices, see Appendix B of JCP 97 (4) 2635
  BigReal mass2 = per * per * kT;
  BigReal mass1 = mass2 * dof;

  for ( i = 0; i < nnhc; i++ ) {
    tNHCzeta[i] = 0;
    tNHCmass[i] = ( i == 0 ) ? mass1 : mass2;
  }

  // try to load the chain variables, ok if it fails
  tNHCLoad();
  CkPrintf("NHC %d, mass1 %g, mass2 %g\n", dof, mass1, mass2);
  for ( i = 0; i < nnhc; i++ )
    CkPrintf("NHC %d: zeta %g, mass %g\n", i+1, tNHCzeta[i], tNHCmass[i]);
  tNHCRescaleFactorPrev = 1.0;
}

void Controller::tNHCDone(int step)
{
  if ( simParams->tNHCOn ) {
    tNHCSave(step);
    delete[] tNHCzeta;
    delete[] tNHCmass;
  }
}

// Nose-Hoover chain thermostat
// Ref.: Nose-Hoover chains: The canonical ensemble via continuous dynamics 
// Glenn J. Martyna, Michael L. Klein, and Mark Tuckerman, JCP 97 (4) 2635
void Controller::tNHCRescaleVelocities(int step, Bool isPrev)
{
  if ( simParams->tNHCOn ) {
    BigReal tp = simParams->tNHCTemp;
    // use the temperature from adaptive tempering, if any
    if ( simParams->adaptTempOn && simParams->adaptTempRescale
     && (step > simParams->adaptTempFirstStep )
     && (!(simParams->adaptTempLastStep > 0) || step < simParams->adaptTempLastStep )) {
      tp = adaptTempT;
    }
    tp *= BOLTZMANN;

    Real dt = simParams->dt * 0.5; // only for half step
    int dof = numDegFreedom;
    if ( dof <= 0 ) dof = Node::Object()->molecule->num_deg_freedom();
    int nnhc = simParams->tNHCLen, i, j, k;
    BigReal s, GQ, mvv, factor;

    mvv = BOLTZMANN * temperature * dof;
    for ( j = nnhc - 1; j >= 0; j-- ) {
      s = ( j == nnhc - 1 ) ? 1 : exp(-tNHCzeta[j+1]*dt*0.25);
      if ( j == 0 ) {
        GQ = mvv - dof * tp;
      } else {
        GQ = tNHCmass[j-1] * tNHCzeta[j-1] * tNHCzeta[j-1] - tp;
      }
      tNHCzeta[j] = (tNHCzeta[j] * s + GQ /tNHCmass[j] * dt*0.5) * s;
    }

    // velocity rescaling factor
    factor = exp( -tNHCzeta[0] * dt );
    if ( !isPrev ) {
      // combine with the rescaling factor carried from
      // the second half-step from the previous MD step
      broadcast->tNHCRescaleFactor.publish(step, factor * tNHCRescaleFactorPrev);
    } else {
      tNHCRescaleFactorPrev = factor;
    }

    BigReal fac2 = factor * factor;
    temperature *= fac2;
    totalEnergy += kineticEnergy * (fac2 - 1);
    kineticEnergy *= fac2;
    kineticEnergyCentered *= fac2;
    mvv *= fac2;

    for ( j = 0; j < nnhc; j++ ) {
      s = ( j == nnhc - 1 ) ? 1 : exp(-tNHCzeta[j+1]*dt*0.25);
      if ( j == 0 ) {
        GQ = mvv - dof * tp;
      } else {
        GQ = tNHCmass[j-1] * tNHCzeta[j-1] * tNHCzeta[j-1] - tp;
      }
      tNHCzeta[j] = (tNHCzeta[j] * s + GQ /tNHCmass[j] * dt*0.5) * s;
    }

    if ( simParams->tNHCFileFreq > 0 && step % simParams->tNHCFileFreq == 0 )
      tNHCSave(step);
  }
}

void Controller::tNHCSave(int step)
{
  if ( !simParams->tNHCOn ) return;

  int i, nnhc = simParams->tNHCLen;
  NAMD_backup_file(simParams->tNHCFile);
  ofstream_namd fs(simParams->tNHCFile);
  if ( !fs ) {
    iout << "Error: cannot write " << simParams->tNHCFile << "\n" << endi;
    return;
  }
  fs << nnhc << " " << step << "\n";
  for ( i = 0; i < nnhc; i++ )
    fs << tNHCzeta[i] << " ";
  fs << "\n";
  for ( i = 0; i < nnhc; i++ )
    fs << tNHCmass[i] << " ";
  fs << "\n";
}

void Controller::tNHCLoad(void)
{
  if ( !simParams->tNHCOn ) return;

  int i, nnhc, step;
  std::ifstream fs(simParams->tNHCFile);
  if ( !fs ) {
    iout << "Cannot read " << simParams->tNHCFile << "\n" << endi;
    return;
  }
  fs >> nnhc >> step;
  if ( nnhc != simParams->tNHCLen ) {
    iout << "Error: NH-chain length mismatch " << nnhc
         << " vs. " << simParams->tNHCLen << "\n" << endi;
    return;
  }
  for ( i = 0; i < nnhc; i++ )
    fs >> tNHCzeta[i];

  if ( simParams->tNHCFileReadMass ) {
    for ( i = 0; i < nnhc; i++ )
      fs >> tNHCmass[i];
  }
  fs.close();
}

void Controller::keHistInit(void)
{
  if ( !simParams->keHistOn ) return;
  BigReal keHistTemp = simParams->thermostatTemp();
  BigReal ke = BOLTZMANN * keHistTemp * numDegFreedom / 2;
  keHistBinMax = (int) (5.0 * ke / simParams->keHistBin);
  CkPrintf("keHistInit: temperature %g, dof %d, keHistBinMax %d\n", keHistTemp, numDegFreedom, keHistBinMax);
  keHist = new BigReal [keHistBinMax];
  int i;
  for ( i = 0; i < keHistBinMax; i++ ) keHist[i] = 0;
  keHistLoad(); // try to load the previous histogram
}

void Controller::keHistUpdate(int step)
{
  if ( !simParams->keHistOn ) return;
  BigReal ke = BOLTZMANN * temperature * numDegFreedom / 2;
  if ( simParams->adaptTempOn ) {
    ke *= simParams->thermostatTemp() / adaptTempT;
  }
  int i = (int) ( ke / simParams->keHistBin );
  if ( i < keHistBinMax ) keHist[i] += 1;
  if ( step > 0 && step % simParams->keHistFileFreq == 0 ) {
    keHistSave(step);
  }
}

void Controller::keHistSave(int step)
{
  if ( !simParams->keHistOn ) return;
  NAMD_backup_file(simParams->keHistFile);
  ofstream_namd fs(simParams->keHistFile);
  if ( !fs ) return;
  fs << "# " << numDegFreedom << " " << step << "\n";
  int i;
  BigReal tot = 0;
  for ( i = 0; i < keHistBinMax; i++ )
    tot += keHist[i];

  // normalization of the reference curve
  BigReal norm = (numDegFreedom % 2) ? 0.5 * log(M_PI) : 0;
  for ( i = 2 - numDegFreedom % 2; i < numDegFreedom; i += 2 )
    norm += log(i*0.5);
  BigReal tp = simParams->thermostatTemp() * BOLTZMANN;
  BigReal dk = simParams->keHistBin;
  for ( i = 0; i < keHistBinMax; i++ ) {
    if ( keHist[i] <= 0 ) continue;
    double hist = keHist[i] / ( dk * tot );
    double ke = (i + 0.5) * dk;
    double histref = exp(log(ke/tp) * (numDegFreedom*0.5-1) -ke/tp - norm) / tp;
    fs << (i + 0.5) * dk << "\t" << hist << "\t" << histref << "\t" << keHist[i] << "\n";
  }
}

void Controller::keHistLoad(void)
{
  if ( !simParams->keHistOn ) return;
  std::ifstream fs(simParams->keHistFile);
  if ( !fs ) return;
  std::string buf;

  std::getline(fs, buf);
  while ( fs.good() ) {
    double ke, hist1, hist2, hist;
    std::getline(fs, buf);
    sscanf(buf.c_str(), "%lf%lf%lf%lf", &ke, &hist1, &hist2, &hist);
    int i = (int) (ke / simParams->keHistBin);
    keHist[i] = hist;
  }
  fs.close();
  iout << "Loaded previous histogram from "
       << simParams->keHistFile << ".\n" << endi;
}

void Controller::keHistDone(int step)
{
  if ( !simParams->keHistOn ) return;
  keHistSave(step);
  delete[] keHist;
}

static char *FORMAT(BigReal X)
{
  static char tmp_string[25];
  const double maxnum = 9999999999.9999;
  if ( X > maxnum ) X = maxnum;
  if ( X < -maxnum ) X = -maxnum;
  sprintf(tmp_string," %14.4f",X); 
  return tmp_string;
}

static char *FORMAT(const char *X)
{
  static char tmp_string[25];
  sprintf(tmp_string," %14s",X); 
  return tmp_string;
}

static char *ETITLE(int X)
{
  static char tmp_string[21];
  sprintf(tmp_string,"ENERGY: %7d",X); 
  return  tmp_string;
}

void Controller::receivePressure(int step, int minimize)
{
    Node *node = Node::Object();
    Molecule *molecule = node->molecule;
    SimParameters *simParameters = node->simParameters;
    Lattice &lattice = state->lattice;

    reduction->require();

    Tensor virial;
    Tensor virial_normal;
    Tensor virial_nbond;
    Tensor virial_slow;
    Tensor pressure_amd;
#ifdef ALTVIRIAL
    Tensor altVirial_normal;
    Tensor altVirial_nbond;
    Tensor altVirial_slow;
#endif
    Tensor intVirial;
    Tensor intVirial_normal;
    Tensor intVirial_nbond;
    Tensor intVirial_slow;
    Vector extForce_normal;
    Vector extForce_nbond;
    Vector extForce_slow;
    BigReal volume;

#if 1
    numDegFreedom = molecule->num_deg_freedom();
    int numGroupDegFreedom = molecule->num_group_deg_freedom();
    int numFixedGroups = molecule->num_fixed_groups();
    int numFixedAtoms = molecule->num_fixed_atoms();
#endif
#if 0
    int numAtoms = molecule->numAtoms;
    numDegFreedom = 3 * numAtoms;
    int numGroupDegFreedom = 3 * molecule->numHydrogenGroups;
    int numFixedAtoms =
	( simParameters->fixedAtomsOn ? molecule->numFixedAtoms : 0 );
    int numLonepairs = molecule->numLonepairs;
    int numFixedGroups = ( numFixedAtoms ? molecule->numFixedGroups : 0 );
    if ( numFixedAtoms ) numDegFreedom -= 3 * numFixedAtoms;
    if (numLonepairs) numDegFreedom -= 3 * numLonepairs;
    if ( numFixedGroups ) numGroupDegFreedom -= 3 * numFixedGroups;
    if ( ! ( numFixedAtoms || molecule->numConstraints
	|| simParameters->comMove || simParameters->langevinOn ) ) {
      numDegFreedom -= 3;
      numGroupDegFreedom -= 3;
    }
    if (simParameters->pairInteractionOn) {
      // this doesn't attempt to deal with fixed atoms or constraints
      numDegFreedom = 3 * molecule->numFepInitial;
    }
    int numRigidBonds = molecule->numRigidBonds;
    int numFixedRigidBonds =
	( simParameters->fixedAtomsOn ? molecule->numFixedRigidBonds : 0 );

    // numLonepairs is subtracted here because all lonepairs have a rigid bond
    // to oxygen, but all of the LP degrees of freedom are dealt with above
    numDegFreedom -= ( numRigidBonds - numFixedRigidBonds - numLonepairs);
#endif

    kineticEnergyHalfstep = reduction->item(REDUCTION_HALFSTEP_KINETIC_ENERGY);
    kineticEnergyCentered = reduction->item(REDUCTION_CENTERED_KINETIC_ENERGY);

    BigReal groupKineticEnergyHalfstep = kineticEnergyHalfstep -
	reduction->item(REDUCTION_INT_HALFSTEP_KINETIC_ENERGY);
    BigReal groupKineticEnergyCentered = kineticEnergyCentered -
	reduction->item(REDUCTION_INT_CENTERED_KINETIC_ENERGY);

    BigReal atomTempHalfstep = 2.0 * kineticEnergyHalfstep
					/ ( numDegFreedom * BOLTZMANN );
    BigReal atomTempCentered = 2.0 * kineticEnergyCentered
					/ ( numDegFreedom * BOLTZMANN );
    BigReal groupTempHalfstep = 2.0 * groupKineticEnergyHalfstep
					/ ( numGroupDegFreedom * BOLTZMANN );
    BigReal groupTempCentered = 2.0 * groupKineticEnergyCentered
					/ ( numGroupDegFreedom * BOLTZMANN );

    /*  test code for comparing different temperatures  
    iout << "TEMPTEST: " << step << " " << 
	atomTempHalfstep << " " <<
	atomTempCentered << " " <<
	groupTempHalfstep << " " <<
	groupTempCentered << "\n" << endi;
  iout << "Number of degrees of freedom: " << numDegFreedom << " " <<
    numGroupDegFreedom << "\n" << endi;
     */

    GET_TENSOR(virial_normal,reduction,REDUCTION_VIRIAL_NORMAL);
    GET_TENSOR(virial_nbond,reduction,REDUCTION_VIRIAL_NBOND);
    GET_TENSOR(virial_slow,reduction,REDUCTION_VIRIAL_SLOW);

#ifdef ALTVIRIAL
    GET_TENSOR(altVirial_normal,reduction,REDUCTION_ALT_VIRIAL_NORMAL);
    GET_TENSOR(altVirial_nbond,reduction,REDUCTION_ALT_VIRIAL_NBOND);
    GET_TENSOR(altVirial_slow,reduction,REDUCTION_ALT_VIRIAL_SLOW);
#endif

    GET_TENSOR(intVirial_normal,reduction,REDUCTION_INT_VIRIAL_NORMAL);
    GET_TENSOR(intVirial_nbond,reduction,REDUCTION_INT_VIRIAL_NBOND);
    GET_TENSOR(intVirial_slow,reduction,REDUCTION_INT_VIRIAL_SLOW);

    GET_VECTOR(extForce_normal,reduction,REDUCTION_EXT_FORCE_NORMAL);
    GET_VECTOR(extForce_nbond,reduction,REDUCTION_EXT_FORCE_NBOND);
    GET_VECTOR(extForce_slow,reduction,REDUCTION_EXT_FORCE_SLOW);
    Vector extPosition = lattice.origin();
    virial_normal -= outer(extForce_normal,extPosition);
    virial_nbond -= outer(extForce_nbond,extPosition);
    virial_slow -= outer(extForce_slow,extPosition);


    kineticEnergy = kineticEnergyCentered;
    temperature = 2.0 * kineticEnergyCentered / ( numDegFreedom * BOLTZMANN );

    if (simParameters->drudeOn) {
      BigReal drudeComKE
        = reduction->item(REDUCTION_DRUDECOM_CENTERED_KINETIC_ENERGY);
      BigReal drudeBondKE
        = reduction->item(REDUCTION_DRUDEBOND_CENTERED_KINETIC_ENERGY);
      int g_bond = 3 * molecule->numDrudeAtoms;
      int g_com = numDegFreedom - g_bond;
      temperature = 2.0 * drudeComKE / (g_com * BOLTZMANN);
      drudeBondTemp = (g_bond!=0 ? (2.*drudeBondKE/(g_bond*BOLTZMANN)) : 0.);
    }

    if ( (volume=lattice.volume()) != 0. )
    {

      if (simParameters->LJcorrection && volume) {
        // Apply tail correction to pressure
        //printf("Volume is %f\n", volume);
        //printf("Applying tail correction of %f to virial\n", molecule->tail_corr_virial / volume);
        virial_normal += Tensor::identity(molecule->tail_corr_virial / volume);
      }

      // kinetic energy component included in virials
      pressure_normal = virial_normal / volume;
      groupPressure_normal = ( virial_normal - intVirial_normal ) / volume;

      if (simParameters->accelMDOn) {
        pressure_amd = virial_amd / volume;
        pressure_normal += pressure_amd;
        groupPressure_normal +=  pressure_amd;
      }

      if ( minimize || ! ( step % nbondFreq ) )
      {
        pressure_nbond = virial_nbond / volume;
        groupPressure_nbond = ( virial_nbond - intVirial_nbond ) / volume;
      }

      if ( minimize || ! ( step % slowFreq ) )
      {
        pressure_slow = virial_slow / volume;
        groupPressure_slow = ( virial_slow - intVirial_slow ) / volume;
      }

/*
      iout << "VIRIALS: " << virial_normal << " " << virial_nbond << " " <<
	virial_slow << " " << ( virial_normal - intVirial_normal ) << " " <<
	( virial_nbond - intVirial_nbond ) << " " <<
	( virial_slow - intVirial_slow ) << "\n";
*/

      pressure = pressure_normal + pressure_nbond + pressure_slow; 
      groupPressure = groupPressure_normal + groupPressure_nbond +
						groupPressure_slow;
    }
    else
    {
      pressure = Tensor();
      groupPressure = Tensor();
    }

    if ( simParameters->useGroupPressure )
    {
      controlPressure_normal = groupPressure_normal;
      controlPressure_nbond = groupPressure_nbond;
      controlPressure_slow = groupPressure_slow;
      controlPressure = groupPressure;
      controlNumDegFreedom = molecule->numHydrogenGroups - numFixedGroups;
      if ( ! ( numFixedAtoms || molecule->numConstraints
	|| simParameters->comMove || simParameters->langevinOn ) ) {
        controlNumDegFreedom -= 1;
      }
    }
    else
    {
      controlPressure_normal = pressure_normal;
      controlPressure_nbond = pressure_nbond;
      controlPressure_slow = pressure_slow;
      controlPressure = pressure;
      controlNumDegFreedom = numDegFreedom / 3;
    }

    if (simParameters->fixCellDims) {
      if (simParameters->fixCellDimX) controlNumDegFreedom -= 1;
      if (simParameters->fixCellDimY) controlNumDegFreedom -= 1;
      if (simParameters->fixCellDimZ) controlNumDegFreedom -= 1;
    }

    if ( simParameters->useFlexibleCell ) {
      // use symmetric pressure to control rotation
      // controlPressure_normal = symmetric(controlPressure_normal);
      // controlPressure_nbond = symmetric(controlPressure_nbond);
      // controlPressure_slow = symmetric(controlPressure_slow);
      // controlPressure = symmetric(controlPressure);
      // only use on-diagonal components for now
      controlPressure_normal = Tensor::diagonal(diagonal(controlPressure_normal));
      controlPressure_nbond = Tensor::diagonal(diagonal(controlPressure_nbond));
      controlPressure_slow = Tensor::diagonal(diagonal(controlPressure_slow));
      controlPressure = Tensor::diagonal(diagonal(controlPressure));
      if ( simParameters->useConstantRatio ) {
#define AVGXY(T) T.xy = T.yx = 0; T.xx = T.yy = 0.5 * ( T.xx + T.yy );\
		 T.xz = T.zx = T.yz = T.zy = 0.5 * ( T.xz + T.yz );
        AVGXY(controlPressure_normal);
        AVGXY(controlPressure_nbond);
        AVGXY(controlPressure_slow);
        AVGXY(controlPressure);
#undef AVGXY
      }
    } else {
      controlPressure_normal =
		Tensor::identity(trace(controlPressure_normal)/3.);
      controlPressure_nbond =
		Tensor::identity(trace(controlPressure_nbond)/3.);
      controlPressure_slow =
		Tensor::identity(trace(controlPressure_slow)/3.);
      controlPressure =
		Tensor::identity(trace(controlPressure)/3.);
    }

#ifdef DEBUG_PRESSURE
    iout << iINFO << "Control pressure = " << controlPressure <<
      " = " << controlPressure_normal << " + " <<
      controlPressure_nbond << " + " << controlPressure_slow << "\n" << endi;
#endif
    if ( simParams->accelMDOn && simParams->accelMDDebugOn && ! (step % simParameters->accelMDOutFreq) ) {
        iout << " PRESS " << trace(pressure)*PRESSUREFACTOR/3. << "\n"
             << " PNORMAL " << trace(pressure_normal)*PRESSUREFACTOR/3. << "\n"
             << " PNBOND " << trace(pressure_nbond)*PRESSUREFACTOR/3. << "\n"
             << " PSLOW " << trace(pressure_slow)*PRESSUREFACTOR/3. << "\n"
             << " PAMD " << trace(pressure_amd)*PRESSUREFACTOR/3. << "\n"
             << " GPRESS " << trace(groupPressure)*PRESSUREFACTOR/3. << "\n"
             << " GPNORMAL " << trace(groupPressure_normal)*PRESSUREFACTOR/3. << "\n"
             << " GPNBOND " << trace(groupPressure_nbond)*PRESSUREFACTOR/3. << "\n"
             << " GPSLOW " << trace(groupPressure_slow)*PRESSUREFACTOR/3. << "\n"
             << endi;
   }
}

void Controller::rescaleaccelMD(int step, int minimize)
{
    if ( !simParams->accelMDOn ) return;

    amd_reduction->require();

    // copy all to submit_reduction
    for ( int i=0; i<REDUCTION_MAX_RESERVED; ++i ) {
      submit_reduction->item(i) += amd_reduction->item(i);
    }
    submit_reduction->submit();

    if (step == simParams->firstTimestep) accelMDdVAverage = 0;
//    if ( minimize || ((step < simParams->accelMDFirstStep ) || (step > simParams->accelMDLastStep ))) return;
    if ( minimize || (step < simParams->accelMDFirstStep ) || ( simParams->accelMDLastStep > 0 && step > simParams->accelMDLastStep )) return;

    Node *node = Node::Object();
    Molecule *molecule = node->molecule;
    Lattice &lattice = state->lattice;

    const BigReal accelMDE = simParams->accelMDE;
    const BigReal accelMDalpha = simParams->accelMDalpha;
    const BigReal accelMDTE = simParams->accelMDTE;
    const BigReal accelMDTalpha = simParams->accelMDTalpha;
    const int accelMDOutFreq = simParams->accelMDOutFreq;

    BigReal bondEnergy;
    BigReal angleEnergy;
    BigReal dihedralEnergy;
    BigReal improperEnergy;
    BigReal crosstermEnergy;
    BigReal boundaryEnergy;
    BigReal miscEnergy;
    BigReal amd_electEnergy;
    BigReal amd_ljEnergy;
    BigReal amd_electEnergySlow;
    BigReal amd_groLJEnergy;
    BigReal amd_groGaussEnergy;
    BigReal amd_goTotalEnergy;
    BigReal potentialEnergy;
    BigReal factor_dihe = 1;
    BigReal factor_tot = 1;
    BigReal testV;
    BigReal dV = 0;
    Vector  accelMDfactor;
    Tensor vir;
    Tensor vir_dihe;
    Tensor vir_normal;
    Tensor vir_nbond;
    Tensor vir_slow;

    bondEnergy = amd_reduction->item(REDUCTION_BOND_ENERGY);
    angleEnergy = amd_reduction->item(REDUCTION_ANGLE_ENERGY);
    dihedralEnergy = amd_reduction->item(REDUCTION_DIHEDRAL_ENERGY);
    improperEnergy = amd_reduction->item(REDUCTION_IMPROPER_ENERGY);
    crosstermEnergy = amd_reduction->item(REDUCTION_CROSSTERM_ENERGY);
    boundaryEnergy = amd_reduction->item(REDUCTION_BC_ENERGY);
    miscEnergy = amd_reduction->item(REDUCTION_MISC_ENERGY);

    GET_TENSOR(vir_dihe,amd_reduction,REDUCTION_VIRIAL_AMD_DIHE);
    GET_TENSOR(vir_normal,amd_reduction,REDUCTION_VIRIAL_NORMAL);
    GET_TENSOR(vir_nbond,amd_reduction,REDUCTION_VIRIAL_NBOND);
    GET_TENSOR(vir_slow,amd_reduction,REDUCTION_VIRIAL_SLOW);

    if ( !( step % nbondFreq ) ) {
      amd_electEnergy = amd_reduction->item(REDUCTION_ELECT_ENERGY);
      amd_ljEnergy = amd_reduction->item(REDUCTION_LJ_ENERGY);
      amd_groLJEnergy = amd_reduction->item(REDUCTION_GRO_LJ_ENERGY);
      amd_groGaussEnergy = amd_reduction->item(REDUCTION_GRO_GAUSS_ENERGY);
      BigReal goNativeEnergy = amd_reduction->item(REDUCTION_GO_NATIVE_ENERGY);
      BigReal goNonnativeEnergy = amd_reduction->item(REDUCTION_GO_NONNATIVE_ENERGY);
      amd_goTotalEnergy = goNativeEnergy + goNonnativeEnergy;
    } else {
      amd_electEnergy = electEnergy;
      amd_ljEnergy = ljEnergy;
      amd_groLJEnergy = groLJEnergy;
      amd_groGaussEnergy = groGaussEnergy;
      amd_goTotalEnergy = goTotalEnergy;
    }

    if ( !( step % slowFreq ) ) {
      amd_electEnergySlow = amd_reduction->item(REDUCTION_ELECT_ENERGY_SLOW);
    } else {
      amd_electEnergySlow = electEnergySlow;
    }

    potentialEnergy = bondEnergy + angleEnergy + dihedralEnergy +
        improperEnergy + amd_electEnergy + amd_electEnergySlow + amd_ljEnergy +
        crosstermEnergy + boundaryEnergy + miscEnergy +
        amd_goTotalEnergy + amd_groLJEnergy + amd_groGaussEnergy;

    if (simParams->accelMDdihe) {

        testV = dihedralEnergy + crosstermEnergy;
        if ( testV < accelMDE ) {
           factor_dihe = accelMDalpha/(accelMDalpha + accelMDE - testV);
           factor_dihe *= factor_dihe;
           vir = vir_dihe * (factor_dihe - 1.0);
           dV = (accelMDE - testV)*(accelMDE - testV)/(accelMDalpha + accelMDE - testV);
           accelMDdVAverage += dV;
        }  

    } else if (simParams->accelMDdual) {

        testV = dihedralEnergy + crosstermEnergy;
        if ( testV < accelMDE ) {
           factor_dihe = accelMDalpha/(accelMDalpha + accelMDE - testV);
           factor_dihe *= factor_dihe;
           vir = vir_dihe * (factor_dihe - 1.0) ;
           dV = (accelMDE - testV)*(accelMDE - testV)/(accelMDalpha + accelMDE - testV);
        }
 
        testV = potentialEnergy - dihedralEnergy - crosstermEnergy;
        if ( testV < accelMDTE ) {
           factor_tot = accelMDTalpha/(accelMDTalpha + accelMDTE - testV);
           factor_tot *= factor_tot;
           vir += (vir_normal - vir_dihe + vir_nbond + vir_slow) * (factor_tot - 1.0);
           dV += (accelMDTE - testV)*(accelMDTE - testV)/(accelMDTalpha + accelMDTE - testV);
        }
       accelMDdVAverage += dV;

    } else {

        testV = potentialEnergy;
        if ( testV < accelMDE ) {
           factor_tot = accelMDalpha/(accelMDalpha + accelMDE - testV);
           factor_tot *= factor_tot;
           vir = (vir_normal + vir_nbond + vir_slow) * (factor_tot - 1.0);
           dV = (accelMDE - testV)*(accelMDE - testV)/(accelMDalpha + accelMDE - testV);
           accelMDdVAverage += dV;
        }
    } 
 
    accelMDfactor[0]=factor_dihe;
    accelMDfactor[1]=factor_tot;
    accelMDfactor[2]=1;
    broadcast->accelMDRescaleFactor.publish(step,accelMDfactor);
    virial_amd = vir; 

    if ( factor_tot < 0.001 ) {
       iout << iWARN << "accelMD is using a very high boost potential, simulation may become unstable!"
            << "\n" << endi;
    }

    if ( ! (step % accelMDOutFreq) ) {
        if ( !(step == simParams->firstTimestep) ) {
             accelMDdVAverage = accelMDdVAverage/accelMDOutFreq; 
        }
        iout << "ACCELERATED MD: STEP " << step
             << " dV "   << dV
             << " dVAVG " << accelMDdVAverage 
             << " BOND " << bondEnergy
             << " ANGLE " << angleEnergy
             << " DIHED " << dihedralEnergy+crosstermEnergy
             << " IMPRP " << improperEnergy
             << " ELECT " << amd_electEnergy+amd_electEnergySlow
             << " VDW "  << amd_ljEnergy
             << " POTENTIAL "  << potentialEnergy << "\n"
             << endi;

        accelMDdVAverage = 0;

        if (simParams->accelMDDebugOn) {
           BigReal volume;
           Tensor p_normal;
           Tensor p_nbond;
           Tensor p_slow;
           Tensor p;
           if ( (volume=lattice.volume()) != 0. )  {
                 p_normal = vir_normal/volume;
                 p_nbond  = vir_nbond/volume;
                 p_slow = vir_slow/volume;
                 p = vir/volume;
           }    
           iout << " accelMD Scaling Factor: " << accelMDfactor << "\n" 
                << " accelMD PNORMAL " << trace(p_normal)*PRESSUREFACTOR/3. << "\n"
                << " accelMD PNBOND " << trace(p_nbond)*PRESSUREFACTOR/3. << "\n"
                << " accelMD PSLOW " << trace(p_slow)*PRESSUREFACTOR/3. << "\n"
                << " accelMD PAMD " << trace(p)*PRESSUREFACTOR/3. << "\n" 
                << endi;
        }
   }
}

void Controller::adaptTempInit(int step) {
    if (!simParams->adaptTempOn) return;
    iout << iINFO << "INITIALISING ADAPTIVE TEMPERING\n" << endi;
    adaptTempMCSize = simParams->adaptTempMCSize;
    adaptTempMCTot = adaptTempMCAcc = 0;
    adaptTempMCDAcc = adaptTempMCFail = 0;
    adaptTempDt = simParams->adaptTempDt;
    adaptTempLangTot = adaptTempLangAcc = 0;
    adaptTempLangDAcc = adaptTempLangFail = 0;
    if (simParams->adaptTempInFile[0] != '\0') {
      iout << iINFO << "READING ADAPTIVE TEMPERING RESTART FILE\n" << endi;
      std::ifstream adaptTempRead(simParams->adaptTempInFile);
      if (adaptTempRead) {
        adaptTempRead.exceptions( std::ios::failbit | std::ios::badbit );
        try {
          std::string buf(""); // now read by lines to make extension easier
          std::stringstream ss;
          ss.exceptions( std::ios::failbit | std::ios::badbit );
          int readInt;
          BigReal readReal;
          std::getline(adaptTempRead, buf);
          ss.str(buf);
          ss >> readInt; // step
          ss >> adaptTempT;     // KELVIN
          ss >> adaptTempBetaMin;  // KELVIN
          ss >> adaptTempBetaMax;  // KELVIN
          adaptTempBetaMin = 1./adaptTempBetaMin; // KELVIN^-1
          adaptTempBetaMax = 1./adaptTempBetaMax; // KELVIN^-1
          // In case file is manually edited
          if (adaptTempBetaMin > adaptTempBetaMax)
            std::swap(adaptTempBetaMin, adaptTempBetaMax);
          ss >> adaptTempBins;     
          ss >> adaptTempCg;
          ss >> adaptTempDt;
          ss.clear(); // clear eof
          adaptTempPotEnergyAveNum  = new double[adaptTempBins];
          adaptTempPotEnergyAveDen  = new double[adaptTempBins];
          adaptTempPotEnergySamples = new long[adaptTempBins];
          adaptTempPotEnergyVarNum  = new double[adaptTempBins];
          adaptTempPotEnergyVar     = new double[adaptTempBins];
          adaptTempPotEnergyAve     = new double[adaptTempBins];
          adaptTempBetaN            = new BigReal[adaptTempBins + 1];
          adaptTempDBeta = (adaptTempBetaMax - adaptTempBetaMin)/(adaptTempBins);
          for(int j = 0; j < adaptTempBins; ++j) {
            std::getline(adaptTempRead, buf);
            ss.str(buf);
            ss >> readReal; // 1 / T
            ss >> adaptTempPotEnergyAve[j];
            ss >> adaptTempPotEnergyVar[j];
            ss >> adaptTempPotEnergySamples[j];
            ss >> adaptTempPotEnergyAveNum[j];
            ss >> adaptTempPotEnergyVarNum[j];
            ss >> adaptTempPotEnergyAveDen[j];
            // ss >> readReal; // InvW
            ss.clear(); // clear eof
            if ( simParams->adaptTempEmptyData ) {
              adaptTempPotEnergyVar[j] = 0;
              adaptTempPotEnergySamples[j] = 0;
              adaptTempPotEnergyAveNum[j] = 0;
              adaptTempPotEnergyVarNum[j] = 0;
              adaptTempPotEnergyAveDen[j] = 0;
            }
          }
          for ( int j = 0; j <= adaptTempBins; j++ )
            adaptTempBetaN[j] = adaptTempBetaMin + j * adaptTempDBeta;
          adaptTempBinMinus = new int[adaptTempBins];
          adaptTempBinPlus  = new int[adaptTempBins];
          // read in data for separate accumulators
          if ( simParams->adaptTempSepOn ) {
            char info[256];
            try {
              std::getline(adaptTempRead, buf);
              if ( strncasecmp(buf.c_str(), "SEP BEGIN", 9) != 0 )
                throw std::ios::failure("");
            } catch ( const std::ios::failure& e ) {
              sprintf(info, "No beginning for separator accumulators, file %s\n",
                  simParams->adaptTempInFile);
              NAMD_die(info);
            }
            adaptTempSepAcc = new AdaptTempSepAcc[adaptTempBins];
            // loop over each accumulator
            int i, j, i1, j1;
            double total;
            for ( i = 0; i < adaptTempBins; i++ ) {
              AdaptTempSepAcc *acc = adaptTempSepAcc + i;
              std::getline(adaptTempRead, buf);
              if ( 4 != sscanf(buf.c_str(), "%d%d%d%lf", &i1, &adaptTempBinMinus[i], &adaptTempBinPlus[i], &total) ) {
                sprintf(info, "Broken info for estimator %d, file %s\n%s",
                    i, simParams->adaptTempInFile, buf.c_str());
                NAMD_die(info);
              }
              acc->init(adaptTempBinMinus[i], adaptTempBinPlus[i]);
              acc->total = total;
              for ( j = 0; j < acc->winSize; j++ ) {
                std::getline(adaptTempRead, buf);
                if ( 5 != sscanf(buf.c_str(), "%d%lf%lf%lf%lf\n", &j1,
                      &acc->sumw[j], &acc->ave[j], &acc->var[j], &acc->cnt[j]) ) {
                  sprintf(info, "Broken member %d/%d for estimator %d, file %s\n%s",
                      j, acc->winSize, i, simParams->adaptTempInFile, buf.c_str());
                  NAMD_die(info);
                }
                acc->sumE[j] = acc->sumw[j] * acc->ave[j];
                acc->sumE2[j] = acc->sumw[j] * (acc->var[j] + acc->ave[j] * acc->ave[j]);
              }
              std::getline(adaptTempRead, buf); // blank line
              if ( simParams->adaptTempEmptyData ) acc->empty();
            }
            std::getline(adaptTempRead, buf);
            if ( strncasecmp(buf.c_str(), "SEP END", 7) != 0 ) {
              sprintf(info, "No ending for separator accumulators, file %s\n",
                  simParams->adaptTempInFile);
              NAMD_die(info);
            }
          }
          try { // to read MC or Langevin data
            std::getline(adaptTempRead, buf);
            const char *p = buf.c_str();
            if ( strncasecmp(p, "MC ", 3) == 0 ) {
              sscanf(p + 3, "%lf%lf%lf%lf%lf", &adaptTempMCSize,
                  &adaptTempMCTot, &adaptTempMCAcc, &adaptTempMCDAcc, &adaptTempMCFail);
            } else if ( strncasecmp(p, "LANGEVIN ", 9) == 0 ) {
              sscanf(p + 9, "%lf%lf%lf%lf", &adaptTempLangTot, &adaptTempLangAcc,
                  &adaptTempLangDAcc, &adaptTempLangFail);
            }
            //CkPrintf("%s\n%g %g %g\n", p, adaptTempMCSize, adaptTempMCTot, adaptTempMCAcc); getchar();
          } catch ( const std::ios::failure& e ) {
            iout << "Failed to read MC/Langevin data from " <<  simParams->adaptTempInFile << "\n" << endi;
          }
        } catch ( const std::ios::failure& e ) {
          NAMD_die("Failed to read the ADAPTIVE TEMPERING restart file.\n");
        }
        adaptTempRead.close();
      }
      else NAMD_die("Could not open ADAPTIVE TEMPERING restart file.\n");
    } 
    else {
      adaptTempBins = simParams->adaptTempBins;
      adaptTempPotEnergyAveNum  = new double[adaptTempBins];
      adaptTempPotEnergyAveDen  = new double[adaptTempBins];
      adaptTempPotEnergySamples = new long[adaptTempBins];
      adaptTempPotEnergyVarNum  = new double[adaptTempBins];
      adaptTempPotEnergyVar     = new double[adaptTempBins];
      adaptTempPotEnergyAve     = new double[adaptTempBins];
      adaptTempBetaN            = new BigReal[adaptTempBins + 1];
      adaptTempBetaMax = 1./simParams->adaptTempTmin;
      adaptTempBetaMin = 1./simParams->adaptTempTmax;
      adaptTempCg = simParams->adaptTempCgamma;   
      adaptTempDBeta = (adaptTempBetaMax - adaptTempBetaMin)/(adaptTempBins);
      adaptTempT = simParams->thermostatTemp();
      for(int j = 0; j < adaptTempBins; ++j){
        adaptTempPotEnergyAveNum[j] = 0.;
        adaptTempPotEnergyAveDen[j] = 0.;
        adaptTempPotEnergySamples[j] = 0;
        adaptTempPotEnergyVarNum[j] = 0.;
        adaptTempPotEnergyVar[j] = 0.;
        adaptTempPotEnergyAve[j] = 0.;
      }
      for ( int j = 0; j <= adaptTempBins; j++ )
        adaptTempBetaN[j] = adaptTempBetaMin + j * adaptTempDBeta;
      // compute the window boundaries
      adaptTempBinMinus = new int[adaptTempBins];
      adaptTempBinPlus  = new int[adaptTempBins];
      for ( int j = 0; j < adaptTempBins; j++ ) {
        BigReal betaMid = adaptTempBetaMin + (j + 0.5) * adaptTempDBeta;
        BigReal deltaBeta = betaMid * simParams->adaptTempWindowSize;
        int deltaBins = (int) (deltaBeta / adaptTempDBeta + 0.5);
        // adjust the window size such that the window does not exceed the boundaries
        if ( j - deltaBins < 0 )
          deltaBins = j;
        if ( j + deltaBins + 1 > adaptTempBins )
          deltaBins = adaptTempBins - j - 1;
        adaptTempBinMinus[j] = j - deltaBins;
        adaptTempBinPlus[j] = j + deltaBins + 1;
      }
      // initialize the separate accumulators
      if ( simParams->adaptTempSepOn ) {
        adaptTempSepAcc = new AdaptTempSepAcc[adaptTempBins];
        for ( int i = 0; i < adaptTempBins; ++i )
          adaptTempSepAcc[i].init(adaptTempBinMinus[i], adaptTempBinPlus[i]);
      }
    }
    // wrap adaptTempT within the range
    BigReal beta = 1.0 / adaptTempT;
    if ( beta <  adaptTempBetaMin )
      beta = adaptTempBetaMin * 0.99 + adaptTempBetaMax * 0.01;
    if ( beta >= adaptTempBetaMax )
      beta = adaptTempBetaMin * 0.01 + adaptTempBetaMax * 0.99;
    adaptTempT = 1.0 / beta;
    iout << iINFO << "ADAPTIVE TEMPERING: TEMPERATURE RANGE: [" << 1./adaptTempBetaMax << "," << 1./adaptTempBetaMin << "] KELVIN\n";
    iout << iINFO << "ADAPTIVE TEMPERING: NUMBER OF BINS TO STORE POT. ENERGY: " << adaptTempBins << "\n";
    iout << iINFO << "ADAPTIVE TEMPERING: ADAPTIVE BIN AVERAGING PARAMETER: " << adaptTempCg << "\n";
    iout << iINFO << "ADAPTIVE TEMPERING: SCALING CONSTANT FOR LANGEVIN TEMPERATURE UPDATES: " << adaptTempDt << "\n";
    iout << iINFO << "ADAPTIVE TEMPERING: INITIAL TEMPERATURE: " << adaptTempT<< "\n";
    if (simParams->adaptTempRestartFreq > 0) {
      NAMD_backup_file(simParams->adaptTempRestartFile);
      adaptTempRestartFile.open(simParams->adaptTempRestartFile);
       if(!adaptTempRestartFile)
        NAMD_die("Error opening restart file");
    }
}

void Controller::adaptTempDone(int step) {
    if (!simParams->adaptTempOn) return;
    adaptTempWriteRestart(step);
    delete [] adaptTempPotEnergyAveNum;
    delete [] adaptTempPotEnergyAveDen;
    delete [] adaptTempPotEnergySamples;
    delete [] adaptTempPotEnergyVarNum;
    delete [] adaptTempPotEnergyVar;
    delete [] adaptTempPotEnergyAve;
    delete [] adaptTempBetaN;
    delete [] adaptTempBinMinus;
    delete [] adaptTempBinPlus;
    if ( simParams->adaptTempSepOn )
      delete [] adaptTempSepAcc;
}

void Controller::adaptTempWriteRestart(int step) {
    if ( simParams->adaptTempOn ) {
        char s[1024];
        if ( !simParams->adaptTempFixedAve ) {
          // compute the average values
          for ( int i = 0; i < adaptTempBins; i++ )
            adaptTempPotEnergyAve[i] = adaptTempGetPEAve(i);
        }
        if ( !simParams->adaptTempRestartAppend ) {
          adaptTempRestartFile.seekbegin();
        }
        iout << "ADAPTEMP: WRITING RESTART FILE AT STEP " << step << "\n" << endi;
        adaptTempRestartFile << step << " ";
        adaptTempRestartFile << adaptTempT << " ";     // KELVIN
        adaptTempRestartFile << 1./adaptTempBetaMin << " ";  // KELVIN
        adaptTempRestartFile << 1./adaptTempBetaMax << " ";  // KELVIN
        adaptTempRestartFile << adaptTempBins << " ";     
        adaptTempRestartFile << adaptTempCg << " ";
        adaptTempRestartFile << adaptTempDt;
        adaptTempRestartFile << "\n";
        for(int j = 0; j < adaptTempBins; ++j) {
          BigReal bet = adaptTempBetaMin + (j + 0.5) * adaptTempDBeta;
          // use printf for better precision control
          sprintf(s, "%.12f\t%14.6f\t%14.6f\t%ld\t%22.14e\t%22.14e\t%22.14e\t%10.6f\n",
              bet, adaptTempPotEnergyAve[j], adaptTempPotEnergyVar[j],
              adaptTempPotEnergySamples[j], adaptTempPotEnergyAveNum[j],
              adaptTempPotEnergyVarNum[j], adaptTempPotEnergyAveDen[j],
              adaptTempGetInvW(1.0 / bet));
          adaptTempRestartFile << s;
        }
        // data for separate accumulators
        if ( simParams->adaptTempSepOn ) {
          adaptTempRestartFile << "SEP BEGIN\n";
          // loop over each accumulator
          for ( int i = 0; i < adaptTempBins; i++ ) {
            AdaptTempSepAcc *acc = adaptTempSepAcc + i;
            acc->trim();
            sprintf(s, "%d %d %d %.0f\n", i,
                adaptTempBinMinus[i], adaptTempBinPlus[i], acc->total);
            adaptTempRestartFile << s;
            for ( int j = 0; j < acc->winSize; j++ ) {
              sprintf(s, "%d\t%22.14e\t%22.14e\t%22.14e\t%.0f\n", j,
                  acc->sumw[j], acc->ave[j], acc->var[j], acc->cnt[j]);
              adaptTempRestartFile << s;
            }
            adaptTempRestartFile << "\n";
          }
          adaptTempRestartFile << "SEP END\n";
        }
        if ( simParams->adaptTempMCMove ) {
          sprintf(s, "MC %g %.0f %.0f %22.14e %.0f\n", adaptTempMCSize,
              adaptTempMCTot, adaptTempMCAcc, adaptTempMCDAcc, adaptTempMCFail);
        } else {
          sprintf(s, "LANGEVIN %.0f %.0f %22.14e %.0f\n", adaptTempLangTot,
              adaptTempLangAcc, adaptTempLangDAcc, adaptTempLangFail);
        }
        adaptTempRestartFile << s;
        adaptTempRestartFile.flush(); 
    }
}    

BigReal Controller::adaptTempGetInvW(BigReal tp)
{
    return pow(BOLTZMANN * tp, -simParams->adaptTempWeightExp);
}

BigReal Controller::adaptTempGetPEAve(int i, BigReal def)
{
    const BigReal varCntMin = simParams->adaptTempFreq * 10;
    BigReal potEnergyAverage;

    if ( simParams->adaptTempSepOn ) {
      potEnergyAverage = adaptTempSepAcc[i].iiave(varCntMin, def);
    } else {
      int j, w0, w1;
      // Get Averaging Limits:
      int nMinus = adaptTempBinMinus[i];
      int nPlus  = adaptTempBinPlus[i];
      // Variables for <E(beta)> estimate:
      BigReal den0 = 0.0, den1 = 0.0, ene0 = 0.0, ene1 = 0.0;
      // Integral terms
      BigReal A0 = 0; // Sum_{from beta_minus to beta_{i+1} }
                      //   (beta - beta_minus)/(beta_{i+1} - beta_minus) var(E)
      BigReal A1 = 0; // Sum_{from beta_{i+1} to beta_plus }
                      //   (beta - beta_plus) /(beta_plus  - beta_{i+1}) var(E)
      BigReal A2 = 0; // 0.5 * DBeta * var(E) at bin i
      //A0 phi_s integral for beta_minus < beta < beta_{i+1}

      // compute the default variance for lack of data
      BigReal denMax = 0, defVar = 0;
      for ( j = nMinus; j < nPlus; j++ ) {
        if ( adaptTempPotEnergyAveDen[j] > denMax ) {
          denMax = adaptTempPotEnergyAveDen[j];
          defVar = adaptTempPotEnergyVar[j];
        }
      }
      if ( defVar <= 0 ) defVar = 1.0;

      BigReal var;
      for ( j = nMinus; j <= i; j++ ) {
        den0 += adaptTempPotEnergyAveDen[j];
        ene0 += adaptTempPotEnergyAveNum[j];
        var = (adaptTempPotEnergySamples[j] > varCntMin) ? adaptTempPotEnergyVar[j] : defVar;
        A0 += var * (j - nMinus + 0.5);
      }
      w0 = i + 1 - nMinus;
      den0 /= w0;
      ene0 /= w0;
      A0 /= w0;

      //A2 phi_t integral for beta_i
      var = (adaptTempPotEnergySamples[i] > varCntMin) ? adaptTempPotEnergyVar[i] : defVar;
      A2 = 0.5 * var;

      //A1 phi_s integral for beta_{i+1} < beta < beta_plus
      for ( j = i + 1; j < nPlus; j++ ) {
        den1 += adaptTempPotEnergyAveDen[j];
        ene1 += adaptTempPotEnergyAveNum[j];
        var = (adaptTempPotEnergySamples[j] > varCntMin) ? adaptTempPotEnergyVar[j] : defVar;
        A1 += var * (j - nPlus + 0.5);
      }
      w1 = nPlus - i - 1;
      if ( w1 > 0 ) {
        den1 /= w1;
        ene1 /= w1;
        A1 /= w1;
      }

      if ( den0 + den1 <= 0 ) {
        potEnergyAverage = def;
      } else {
        // Now calculate a+ and a-
        BigReal aplus = ( den0 + den1 > 0 ) ? (A0 - A2) / (A0 - A1) : 0;
        if ( aplus < 0 ) aplus = 0;
        if ( aplus > 1 ) aplus = 1;
        BigReal aminus = 1 - aplus;
        potEnergyAverage = ( aminus * ene0 + aplus * ene1 )
                         / ( aminus * den0 + aplus * den1 );
        if (simParams->adaptTempDebug) {
          iout << "ADAPTEMP DEBUG:"  << "\n"
               << "     adaptTempBin:    " << i << "\n"
               << "     Samples:   " << adaptTempPotEnergySamples[i] << "\n"
               << "     potentialEnergeAverage:  " << potEnergyAverage << "\n"
               << "     adaptTemp:   " << adaptTempT<< "\n"
               << "     betaMin:   " << adaptTempBetaMin << "\n"
               << "     betaMax:   " << adaptTempBetaMax << "\n"
               << "     deltaBeta: " << (i - nMinus) * adaptTempDBeta << "\n"
               << "     betaMinus: " << adaptTempBetaN[nMinus] << "\n"
               << "     betaPlus:  " << adaptTempBetaN[nPlus] << "\n"
               << "     nMinus:    " << nMinus << "\n"
               << "     nPlus:     " << nPlus << "\n"
               << "     A0:        " << A0 << "\n"
               << "     A1:        " << A1 << "\n"
               << "     A2:        " << A2 << "\n"
               << "     a+:        " << aplus << "\n"
               << "     a-:        " << aminus << "\n"
               << "     aveEner:   " << potEnergyAverage << "\n"
               << "     aveEne0:   " << ene0 / (den0 + 1e-16) << "\n"
               << "     aveEne1:   " << ene1 / (den1 + 1e-16) << "\n"
               << "     den0:      " << den0 << "\n"
               << "     den1:      " << den1 << "\n"
               << endi;
        }
      }
    }
    return potEnergyAverage;
}

// return Integral { beta to nbeta } E(beta) d beta = Z(beta) - Z(nbeta)
BigReal Controller::adaptTempGetIntE(BigReal beta, int i, BigReal nbeta, int ni)
{
    double delta = 0, epave = 0, beta_n;
    int j, sgn = ( i <= ni ) ? 1 : -1;
    for ( j = i; ; j += sgn ) {
      if ( !simParams->adaptTempFixedAve ) // recompute the average energy
        adaptTempPotEnergyAve[j] = epave = adaptTempGetPEAve(j, epave);
      beta_n = ( j == ni ) ? nbeta : adaptTempBetaN[j + (sgn > 0)];
      delta += adaptTempPotEnergyAve[j] * (beta_n - beta);
      beta = beta_n;
      if ( j == ni ) break;
    }
    return delta / BOLTZMANN;
}

BigReal Controller::adaptTempMCMove(BigReal tp, BigReal ep)
{
    double kB = BOLTZMANN, x = simParams->adaptTempWeightExp;
    double beta = 1./tp, nbeta, r, delta, nepave;
    int i, ni, acc = 0;
    r = random->gaussian();
    nbeta = beta * exp(adaptTempMCSize * r); // evenly change ln(beta)
    i  = (int) ( (beta  - adaptTempBetaMin) / adaptTempDBeta );
    ni = (int) ( (nbeta - adaptTempBetaMin) / adaptTempDBeta );
    if ( nbeta >= adaptTempBetaMin && ni < adaptTempBins ) {
      delta = adaptTempGetIntE(beta, i, nbeta, ni);
      nepave = adaptTempPotEnergyAve[ni];
      delta = delta - ep * (nbeta - beta)/kB + (x - 1) * log(beta/nbeta);
      //CkPrintf("delta %g cf %g, beta %g, %g, ep %g, %g, bin %d, %d\n", delta,
      //    ((adaptTempPotEnergyAve[i]+adaptTempPotEnergyAve[ni])/2 - ep) * (nbeta - beta) / kB,
      //    beta, nbeta, ep, nepave, i, ni); // getchar();
      acc = ( delta > 0 || random->uniform() < exp(delta) );
      if ( acc ) { // for d(acc. ratio)/d(ln beta)
        if ( delta < 0 )
          adaptTempMCDAcc += -r * ((ep - nepave)*(nbeta/kB) + x - 1);
        double mbeta = beta * exp((adaptTempMCSize + simParams->adaptTempMCSizeInc) * r);
        if ( mbeta < adaptTempBetaMin || mbeta >= adaptTempBetaMax )
          adaptTempMCFail += 1;
      }
    }
    adaptTempMCTot += 1;
    adaptTempMCAcc += acc;
    // adjust the MC move size automatically
    if ( simParams->adaptTempMCAutoAR > 0 && adaptTempMCTot > 100 ) {
      double dacc = adaptTempMCDAcc;
      if ( dacc > 0 ) dacc = 0;
      if ( adaptTempMCFail > 5 ) { // out-of-boundary rejections
        dacc -= adaptTempMCFail / simParams->adaptTempMCSizeInc;
      } else { // approximation
        dacc -= 2 * adaptTempMCAcc;
      }
      double del = ( simParams->adaptTempMCAutoAR - acc ) / dacc;
      if ( del >  0.5 * adaptTempMCSize ) del =  0.5 * adaptTempMCSize;
      if ( del < -0.5 * adaptTempMCSize ) del = -0.5 * adaptTempMCSize;
      adaptTempMCSize += del;
    }
    return acc ? 1.0/nbeta : tp;
}

BigReal Controller::adaptTempLangevin(BigReal tp, BigReal ep)
{
    double kB = BOLTZMANN, x = simParams->adaptTempWeightExp;
    double beta = 1./tp, Beta = beta/kB;
    int i = (int) ( (beta - adaptTempBetaMin) / adaptTempDBeta );
    double epave = adaptTempGetPEAve(i, 0);
    double r = random->gaussian();
    double dt = adaptTempDt, a = sqrt(2. * dt);
    double de = ep - epave + x / Beta;
    double dtp = dt * de / kB + r * a * tp;
    double ntp = tp + dtp;
    double nbeta = 1./ntp, nBeta = nbeta/kB;
    int ni = (int) ( (nbeta - adaptTempBetaMin) / adaptTempDBeta );
    // use Monte Carlo to accept the new temperature
    int acc = 0;
    if ( nbeta >= adaptTempBetaMin && ni < adaptTempBins ) {
      double delta = adaptTempGetIntE(beta, i, nbeta, ni);
      double nepave = adaptTempPotEnergyAve[ni];
      double nde = ep - nepave + x / nBeta;
      double nr = (-dtp/ntp - dt*nde*nBeta) / a;
      delta = (Beta - nBeta) * ep + delta
            + (x - 3) * log(Beta/nBeta) + (r*r - nr*nr)/2;
      acc = ( delta > 0 || random->uniform() < exp(delta) );
      //CkPrintf("acc %d, delta %g, beta %g, %g, ep %g, i %d, %d\n", acc, delta, beta, nbeta, ep, i, ni);
      if ( acc ) { // for d(acc. ratio)/da
        if ( delta < 0 ) {
          BigReal dnB = -nBeta*nBeta/Beta * (r + Beta*de*a); // d(nBeta)/da
          BigReal nvar = adaptTempPotEnergyVar[ adaptTempPotEnergySamples[ni] > 0 ? ni : i ]; // -E'(Beta)
          BigReal dnr = -(nBeta*nde + nr/a) + (1./Beta/a - (ep - nepave + nBeta*nvar)*a/2) * dnB; // d(nr)/da
          adaptTempLangDAcc += -((nde - 3./nBeta) * dnB + nr * dnr); // d(delta)/da
        }
        double ma = a + simParams->adaptTempMCSizeInc;
        double mbeta = 1./(tp + (ma*ma/2) * de/kB + r * ma * tp);
        if ( mbeta < adaptTempBetaMin || mbeta >= adaptTempBetaMax )
          adaptTempLangFail += 1;
      }
    }
    adaptTempLangTot += 1;
    adaptTempLangAcc += acc;
    // adjust the Langevin integration time automatically
    if ( simParams->adaptTempDtAutoAR > 0 && adaptTempLangTot > 100 ) {
      double dacc = adaptTempLangDAcc;
      if ( dacc > 0 ) dacc = 0;
      if ( adaptTempLangFail > 5 ) { // out-of-boundary rejections
        dacc -= adaptTempLangFail / simParams->adaptTempMCSizeInc;
      } else { // approximation
        dacc -= 2 * adaptTempLangAcc;
      }
      double del = ( simParams->adaptTempDtAutoAR - acc ) / dacc;
      if ( del >  0.5 * a ) del =  0.5 * a;
      if ( del < -0.5 * a ) del = -0.5 * a;
      adaptTempDt = (a + del) * (a + del) / 2;
    }
    return acc ? ntp : tp;
}

Bool Controller::adaptTempUpdate(int step, int minimize)
{
    Bool scaled = FALSE;
    //Beta = 1./T
    if ( !simParams->adaptTempOn ) return scaled;
    int j = 0;
    if (step == simParams->firstTimestep) {
        adaptTempInit(step);
        return scaled;
    }
    if ( minimize || (step < simParams->adaptTempFirstStep ) || 
        ( simParams->adaptTempLastStep > 0 && step > simParams->adaptTempLastStep )) return scaled;
    //Calculate Current inverse temperature and bin 
    BigReal adaptTempBeta = 1./adaptTempT;
    adaptTempBin = (int) ( (adaptTempBeta - adaptTempBetaMin) / adaptTempDBeta );

    if ( adaptTempBeta < adaptTempBetaMin || adaptTempBin >= adaptTempBins )
        iout << iWARN << " adaptTempBin out of range: adaptTempBin: " << adaptTempBin  
                      << " adaptTempBeta: " << adaptTempBeta 
                      << " adaptTempDBeta: " << adaptTempDBeta 
                      << " betaMin:" << adaptTempBetaMin 
                      << " betaMax: " << adaptTempBetaMax << "\n";
    if ( adaptTempBin < 0 ) {
      adaptTempBin = 0;
    } else if ( adaptTempBin >= adaptTempBins ) {
      adaptTempBin = adaptTempBins - 1;
    }
    adaptTempPotEnergySamples[adaptTempBin] += 1;
    BigReal gammaAve = 1.-adaptTempCg/adaptTempPotEnergySamples[adaptTempBin];
    if ( gammaAve < 0 ) gammaAve = 0;

    BigReal potentialEnergy;
    BigReal potEnergyAverage;
    BigReal potEnergyVariance;

    // recompute the potential energy
    BigReal bondEnergy = reduction->item(REDUCTION_BOND_ENERGY);
    BigReal angleEnergy = reduction->item(REDUCTION_ANGLE_ENERGY);
    BigReal dihedralEnergy = reduction->item(REDUCTION_DIHEDRAL_ENERGY);
    BigReal improperEnergy = reduction->item(REDUCTION_IMPROPER_ENERGY);
    BigReal crosstermEnergy = reduction->item(REDUCTION_CROSSTERM_ENERGY);
    BigReal boundaryEnergy = reduction->item(REDUCTION_BC_ENERGY);
    BigReal miscEnergy = reduction->item(REDUCTION_MISC_ENERGY);
    if ( step % nbondFreq == 0 ) {
      electEnergy = reduction->item(REDUCTION_ELECT_ENERGY);
      ljEnergy = reduction->item(REDUCTION_LJ_ENERGY);
      BigReal volume = state->lattice.volume();
      if (simParams->LJcorrection && volume) {
        // Apply tail correction to energy
        ljEnergy += Node::Object()->molecule->tail_corr_ener / volume;
      }
      groLJEnergy = reduction->item(REDUCTION_GRO_LJ_ENERGY);
      groGaussEnergy = reduction->item(REDUCTION_GRO_GAUSS_ENERGY);
      goNativeEnergy = reduction->item(REDUCTION_GO_NATIVE_ENERGY);
      goNonnativeEnergy = reduction->item(REDUCTION_GO_NONNATIVE_ENERGY);
      goTotalEnergy = goNativeEnergy + goNonnativeEnergy;
    }
    if ( step % slowFreq == 0 ) {
      electEnergySlow = reduction->item(REDUCTION_ELECT_ENERGY_SLOW);
    }
    potentialEnergy = bondEnergy + angleEnergy + dihedralEnergy
	+ improperEnergy + electEnergy + electEnergySlow + ljEnergy
        + crosstermEnergy + boundaryEnergy + miscEnergy
        + goTotalEnergy + groLJEnergy + groGaussEnergy;
    totalEnergy = potentialEnergy + kineticEnergy;

    BigReal invw = adaptTempGetInvW(adaptTempT);
    //calculate new bin average and variance using adaptive averaging
    adaptTempPotEnergyAveNum[adaptTempBin] = adaptTempPotEnergyAveNum[adaptTempBin]*gammaAve + potentialEnergy * invw;
    adaptTempPotEnergyAveDen[adaptTempBin] = adaptTempPotEnergyAveDen[adaptTempBin]*gammaAve + invw;
    adaptTempPotEnergyVarNum[adaptTempBin] = adaptTempPotEnergyVarNum[adaptTempBin]*gammaAve + potentialEnergy*potentialEnergy * invw;
    
    potEnergyAverage = adaptTempPotEnergyAveNum[adaptTempBin]/adaptTempPotEnergyAveDen[adaptTempBin];
    potEnergyVariance = adaptTempPotEnergyVarNum[adaptTempBin]/adaptTempPotEnergyAveDen[adaptTempBin];
    potEnergyVariance -= potEnergyAverage*potEnergyAverage;

    if ( !simParams->adaptTempFixedAve )
      adaptTempPotEnergyAve[adaptTempBin] = potEnergyAverage;
    adaptTempPotEnergyVar[adaptTempBin] = potEnergyVariance;

    if ( simParams->adaptTempSepOn ) { // update separate accumulators
      for ( j = adaptTempBin; j >= 0; j-- ) { // search downward
        if ( adaptTempBinPlus[j] <= adaptTempBin ) break;
        adaptTempSepAcc[j].add(adaptTempBin, potentialEnergy, invw, adaptTempCg);
      }
      for ( j = adaptTempBin + 1; j < adaptTempBins; j++ ) { // search upward
        if ( adaptTempBinMinus[j] > adaptTempBin ) break;
        adaptTempSepAcc[j].add(adaptTempBin, potentialEnergy, invw, adaptTempCg);
      }
    }

    // Weighted integral of <Delta E^2>_beta dbeta <= Eq 4 of JCP 132 244101
    // Integrals of Eqs 5 and 6 is done as piecewise assuming <Delta E^2>_beta
    // is constant for each bin. This is to estimate <E(beta)> where beta \in
    // (beta_i,beta_{i+1}) using Eq 2 of JCP 132 244101
    if ( ! ( step % simParams->adaptTempFreq ) ) {
      BigReal dT; // dT is the new temperature
      if ( adaptTempPotEnergySamples[adaptTempBin] <= simParams->adaptTempSamplesMin ) {
        // avoid making temperature transitions without enough samples
        dT = adaptTempT;
      } else if ( simParams->adaptTempMCMove ) { // Monte Carlo temperature update
        dT = adaptTempMCMove(adaptTempT, potentialEnergy);
      } else { // Langevin equation temperature update
        dT = adaptTempLangevin(adaptTempT, potentialEnergy);
      }

      BigReal tScale = dT / adaptTempT;
      BigReal vScale = sqrt(tScale);
      // for velocity-rescaling-based thermostats,
      // carry the velocity-rescaling factor to the next step
      if ( simParams->langRescaleOn ) {
        langRescaleFactorPrev *= vScale;
      } else if ( simParams->tNHCOn ) {
        tNHCRescaleFactorPrev *= vScale;
      }
      adaptTempT = dT; 
      broadcast->adaptTemperature.publish(step,adaptTempT);
      scaled = TRUE;
      // temperature is to be used for the Langevin velocity-rescaling
      // and NH-chain thermostats, so it needs to be updated.
      temperature *= tScale;
      totalEnergy += kineticEnergy * (tScale - 1);
      kineticEnergy *= tScale;
      kineticEnergyCentered *= tScale;
    }
    if ( step % simParams->adaptTempRestartFreq == 0 )
      adaptTempWriteRestart(step);
    if ( ! (step % simParams->adaptTempOutFreq) ) {
      char info[256];
      sprintf(info, "ADAPTEMP: STEP %d TEMP %g BIN %d ENERGY %.5f AVG %g VAR %g",
          step, adaptTempT, adaptTempBin, potentialEnergy,
          adaptTempPotEnergyAve[adaptTempBin], adaptTempPotEnergyVar[adaptTempBin]);
      iout << info;
      if ( simParams->adaptTempMCMove ) { // Monte Carlo
        if ( adaptTempMCTot > 0 ) {
          BigReal acc = adaptTempMCAcc / adaptTempMCTot;
          BigReal dacc = adaptTempMCDAcc / adaptTempMCTot;
          if ( dacc > 0 ) dacc = 0;
          dacc -= adaptTempMCFail / simParams->adaptTempMCSizeInc / adaptTempMCTot;
          BigReal ar = simParams->adaptTempMCAutoAR;
          if ( ar <= 0 ) ar = 0.5;
          BigReal newsize = adaptTempMCSize + (ar - acc) / dacc;
          if ( newsize < 0 ) newsize = 0;
          sprintf(info, " MC %.0f(%.0f) ACC. RATIO %.3f%% DAR %g SIZE %g -> %g",
              adaptTempMCTot, adaptTempMCFail, 100*acc, dacc, adaptTempMCSize, newsize);
          iout << info;
        }
      } else { // Langevin equation
        if ( adaptTempLangTot > 0 ) {
          BigReal acc = adaptTempLangAcc / adaptTempLangTot;
          BigReal dacc = adaptTempLangDAcc / adaptTempLangTot;
          if ( dacc > 0 ) dacc = 0;
          dacc -= adaptTempLangFail / simParams->adaptTempMCSizeInc / adaptTempLangTot;
          BigReal oldsize = sqrt(2*adaptTempDt);
          BigReal ar = simParams->adaptTempDtAutoAR;
          if ( ar <= 0 ) ar = 0.5;
          BigReal newsize = oldsize + (ar - acc) / dacc;
          if ( newsize < 0 ) newsize = 0;
          BigReal newdt = newsize * newsize / 2;
          sprintf(info, " LANGEVIN %.0f(%.0f) ACC. RATIO %.3f%% DAR %g SIZE %g -> %g, DT %g -> %g",
              adaptTempLangTot, adaptTempLangFail, 100*acc, dacc, oldsize, newsize, adaptTempDt, newdt);
          iout << info;
        }
      }
      iout << "\n" << endi;
   }
   return scaled;
}


void Controller::compareChecksums(int step, int forgiving) {
    Node *node = Node::Object();
    Molecule *molecule = node->molecule;

    // Some basic correctness checking
    BigReal checksum, checksum_b;

    checksum = reduction->item(REDUCTION_ATOM_CHECKSUM);
    if ( ((int)checksum) != molecule->numAtoms ) {
      char errmsg[256];
      sprintf(errmsg, "Bad global atom count! (%d vs %d)\n",
              (int)checksum, molecule->numAtoms);
      NAMD_bug(errmsg);
    }

    checksum = reduction->item(REDUCTION_COMPUTE_CHECKSUM);
    if ( ((int)checksum) && ((int)checksum) != computeChecksum ) {
      if ( ! computeChecksum ) {
        computesPartitioned = 0;
      } else if ( (int)checksum > computeChecksum && ! computesPartitioned ) {
        computesPartitioned = 1;
      } else {
        NAMD_bug("Bad global compute count!\n");
      }
      computeChecksum = ((int)checksum);
    }

    checksum_b = checksum = reduction->item(REDUCTION_BOND_CHECKSUM);
    if ( checksum_b && (((int)checksum) != molecule->numCalcBonds) ) {
      if ( forgiving && (((int)checksum) < molecule->numCalcBonds) )
        iout << iWARN << "Bad global bond count!\n" << endi;
      else NAMD_bug("Bad global bond count!\n");
    }

    checksum = reduction->item(REDUCTION_ANGLE_CHECKSUM);
    if ( checksum_b && (((int)checksum) != molecule->numCalcAngles) ) {
      if ( forgiving && (((int)checksum) < molecule->numCalcAngles) )
        iout << iWARN << "Bad global angle count!\n" << endi;
      else NAMD_bug("Bad global angle count!\n");
    }

    checksum = reduction->item(REDUCTION_DIHEDRAL_CHECKSUM);
    if ( checksum_b && (((int)checksum) != molecule->numCalcDihedrals) ) {
      if ( forgiving && (((int)checksum) < molecule->numCalcDihedrals) )
        iout << iWARN << "Bad global dihedral count!\n" << endi;
      else NAMD_bug("Bad global dihedral count!\n");
    }

    checksum = reduction->item(REDUCTION_IMPROPER_CHECKSUM);
    if ( checksum_b && (((int)checksum) != molecule->numCalcImpropers) ) {
      if ( forgiving && (((int)checksum) < molecule->numCalcImpropers) )
        iout << iWARN << "Bad global improper count!\n" << endi;
      else NAMD_bug("Bad global improper count!\n");
    }

    checksum = reduction->item(REDUCTION_THOLE_CHECKSUM);
    if ( checksum_b && (((int)checksum) != molecule->numCalcTholes) ) {
      if ( forgiving && (((int)checksum) < molecule->numCalcTholes) )
        iout << iWARN << "Bad global Thole count!\n" << endi;
      else NAMD_bug("Bad global Thole count!\n");
    }

    checksum = reduction->item(REDUCTION_ANISO_CHECKSUM);
    if ( checksum_b && (((int)checksum) != molecule->numCalcAnisos) ) {
      if ( forgiving && (((int)checksum) < molecule->numCalcAnisos) )
        iout << iWARN << "Bad global Aniso count!\n" << endi;
      else NAMD_bug("Bad global Aniso count!\n");
    }

    checksum = reduction->item(REDUCTION_CROSSTERM_CHECKSUM);
    if ( checksum_b && (((int)checksum) != molecule->numCalcCrossterms) ) {
      if ( forgiving && (((int)checksum) < molecule->numCalcCrossterms) )
        iout << iWARN << "Bad global crossterm count!\n" << endi;
      else NAMD_bug("Bad global crossterm count!\n");
    }

    checksum = reduction->item(REDUCTION_EXCLUSION_CHECKSUM);
    if ( ((int)checksum) > molecule->numCalcExclusions &&
         ( ! simParams->mollyOn || step % slowFreq ) ) {
      if ( forgiving )
        iout << iWARN << "High global exclusion count ("
                      << ((int)checksum) << " vs "
                      << molecule->numCalcExclusions << "), possible error!\n"
          << iWARN << "This warning is not unusual during minimization.\n"
          << iWARN << "Decreasing pairlistdist or cutoff that is too close to periodic cell size may avoid this.\n" << endi;
      else {
        char errmsg[256];
        sprintf(errmsg, "High global exclusion count!  (%d vs %d)  System unstable or pairlistdist or cutoff too close to periodic cell size.\n",
                (int)checksum, molecule->numCalcExclusions);
        NAMD_bug(errmsg);
      }
    }

    if ( ((int)checksum) && ((int)checksum) < molecule->numCalcExclusions && ! simParams->goGroPair ) {
      if ( forgiving || ! simParams->fullElectFrequency ) {
        iout << iWARN << "Low global exclusion count!  ("
          << ((int)checksum) << " vs " << molecule->numCalcExclusions << ")";
        if ( forgiving ) {
          iout << "\n"
            << iWARN << "This warning is not unusual during minimization.\n"
            << iWARN << "Increasing pairlistdist or cutoff may avoid this.\n";
        } else {
          iout << "  System unstable or pairlistdist or cutoff too small.\n";
        }
        iout << endi;
      } else {
        char errmsg[256];
        sprintf(errmsg, "Low global exclusion count!  (%d vs %d)  System unstable or pairlistdist or cutoff too small.\n",
                (int)checksum, molecule->numCalcExclusions);
        NAMD_bug(errmsg);
      }
    }

#ifdef NAMD_CUDA
    checksum = reduction->item(REDUCTION_EXCLUSION_CHECKSUM_CUDA);
    if ( ((int)checksum) > molecule->numCalcFullExclusions &&
         ( ! simParams->mollyOn || step % slowFreq ) ) {
      if ( forgiving )
        iout << iWARN << "High global CUDA exclusion count ("
                      << ((int)checksum) << " vs "
                      << molecule->numCalcFullExclusions << "), possible error!\n"
          << iWARN << "This warning is not unusual during minimization.\n"
          << iWARN << "Decreasing pairlistdist or cutoff that is too close to periodic cell size may avoid this.\n" << endi;
      else {
        char errmsg[256];
        sprintf(errmsg, "High global CUDA exclusion count!  (%d vs %d)  System unstable or pairlistdist or cutoff too close to periodic cell size.\n",
                (int)checksum, molecule->numCalcFullExclusions);
        NAMD_bug(errmsg);
      }
    }

    if ( ((int)checksum) && ((int)checksum) < molecule->numCalcFullExclusions && ! simParams->goGroPair ) {
      if ( forgiving || ! simParams->fullElectFrequency ) {
        iout << iWARN << "Low global CUDA exclusion count!  ("
          << ((int)checksum) << " vs " << molecule->numCalcFullExclusions << ")";
        if ( forgiving ) {
          iout << "\n"
            << iWARN << "This warning is not unusual during minimization.\n"
            << iWARN << "Increasing pairlistdist or cutoff may avoid this.\n";
        } else {
          iout << "  System unstable or pairlistdist or cutoff too small.\n";
        }
        iout << endi;
      } else {
        char errmsg[256];
        sprintf(errmsg, "Low global CUDA exclusion count!  (%d vs %d)  System unstable or pairlistdist or cutoff too small.\n",
                (int)checksum, molecule->numCalcFullExclusions);
        NAMD_bug(errmsg);
      }
    }
#endif

    checksum = reduction->item(REDUCTION_MARGIN_VIOLATIONS);
    if ( ((int)checksum) && ! marginViolations ) {
      iout << iERROR << "Margin is too small for " << ((int)checksum) <<
        " atoms during timestep " << step << ".\n" << iERROR <<
      "Incorrect nonbonded forces and energies may be calculated!\n" << endi;
    }
    marginViolations += (int)checksum;

    checksum = reduction->item(REDUCTION_PAIRLIST_WARNINGS);
    if ( simParams->outputPairlists && ((int)checksum) && ! pairlistWarnings ) {
      iout << iINFO <<
        "Pairlistdist is too small for " << ((int)checksum) <<
        " computes during timestep " << step << ".\n" << endi;
    }
    if ( simParams->outputPairlists )  pairlistWarnings += (int)checksum;

    checksum = reduction->item(REDUCTION_STRAY_CHARGE_ERRORS);
    if ( checksum ) {
      if ( forgiving )
        iout << iWARN << "Stray PME grid charges ignored!\n" << endi;
      else NAMD_bug("Stray PME grid charges detected!\n");
    }
}

void Controller::printTiming(int step) {

    if ( simParams->outputTiming && ! ( step % simParams->outputTiming ) )
    {
      const double endWTime = CmiWallTimer() - firstWTime;
      const double endCTime = CmiTimer() - firstCTime;

      // fflush about once per minute
      if ( (((int)endWTime)>>6) != (((int)startWTime)>>6 ) ) fflush_count = 2;

      const double elapsedW = 
	(endWTime - startWTime) / simParams->outputTiming;
      const double elapsedC = 
	(endCTime - startCTime) / simParams->outputTiming;

      const double remainingW = elapsedW * (simParams->N - step);
      const double remainingW_hours = remainingW / 3600;

      startWTime = endWTime;
      startCTime = endCTime;

#ifdef NAMD_CUDA
      if ( simParams->outputEnergies < 60 &&
           step < (simParams->firstTimestep + 10 * simParams->outputTiming) ) {
        iout << iWARN << "Energy evaluation is expensive, increase outputEnergies to improve performance.\n" << endi;
      }
#endif
      if ( step >= (simParams->firstTimestep + simParams->outputTiming) ) {
	CmiPrintf("TIMING: %d  CPU: %g, %g/step  Wall: %g, %g/step"
		  ", %g hours remaining, %f MB of memory in use.\n",
		  step, endCTime, elapsedC, endWTime, elapsedW,
		  remainingW_hours, memusage_MB());
        if ( fflush_count ) { --fflush_count; fflush(stdout); }
      }
    }
}

void Controller::printMinimizeEnergies(int step) {

    rescaleaccelMD(step,1);
    receivePressure(step,1);

    Node *node = Node::Object();
    Molecule *molecule = node->molecule;
    compareChecksums(step,1);

    printEnergies(step,1);

    min_energy = totalEnergy;
    min_f_dot_f = reduction->item(REDUCTION_MIN_F_DOT_F);
    min_f_dot_v = reduction->item(REDUCTION_MIN_F_DOT_V);
    min_v_dot_v = reduction->item(REDUCTION_MIN_V_DOT_V);
    min_huge_count = (int) (reduction->item(REDUCTION_MIN_HUGE_COUNT));

}

void Controller::printDynamicsEnergies(int step) {

    Node *node = Node::Object();
    Molecule *molecule = node->molecule;
    SimParameters *simParameters = node->simParameters;
    Lattice &lattice = state->lattice;

    compareChecksums(step);

    printEnergies(step,0);
}

void Controller::printEnergies(int step, int minimize)
{
    Node *node = Node::Object();
    Molecule *molecule = node->molecule;
    SimParameters *simParameters = node->simParameters;
    Lattice &lattice = state->lattice;

    // Drude model ANISO energy is added into BOND energy
    // and THOLE energy is added into ELECT energy

    BigReal bondEnergy;
    BigReal angleEnergy;
    BigReal dihedralEnergy;
    BigReal improperEnergy;
    BigReal crosstermEnergy;
    BigReal boundaryEnergy;
    BigReal miscEnergy;
    BigReal potentialEnergy;
    BigReal flatEnergy;
    BigReal smoothEnergy;

    Vector momentum;
    Vector angularMomentum;
    BigReal volume = lattice.volume();

    bondEnergy = reduction->item(REDUCTION_BOND_ENERGY);
    angleEnergy = reduction->item(REDUCTION_ANGLE_ENERGY);
    dihedralEnergy = reduction->item(REDUCTION_DIHEDRAL_ENERGY);
    improperEnergy = reduction->item(REDUCTION_IMPROPER_ENERGY);
    crosstermEnergy = reduction->item(REDUCTION_CROSSTERM_ENERGY);
    boundaryEnergy = reduction->item(REDUCTION_BC_ENERGY);
    miscEnergy = reduction->item(REDUCTION_MISC_ENERGY);

    if ( minimize || ! ( step % nbondFreq ) )
    {
      electEnergy = reduction->item(REDUCTION_ELECT_ENERGY);
      ljEnergy = reduction->item(REDUCTION_LJ_ENERGY);

      if (simParameters->LJcorrection && volume) {
        // Apply tail correction to energy
        //printf("Volume is %f\n", volume);
        //printf("Applying tail correction of %f to energy\n", molecule->tail_corr_ener / volume);
        ljEnergy += molecule->tail_corr_ener / volume;
        ljEnergy_f += molecule->tail_corr_ener / volume;
        ljEnergy_f_left += molecule->tail_corr_ener / volume;
      }

      // JLai
      groLJEnergy = reduction->item(REDUCTION_GRO_LJ_ENERGY);
      groGaussEnergy = reduction->item(REDUCTION_GRO_GAUSS_ENERGY);

      goNativeEnergy = reduction->item(REDUCTION_GO_NATIVE_ENERGY);
      goNonnativeEnergy = reduction->item(REDUCTION_GO_NONNATIVE_ENERGY);
      goTotalEnergy = goNativeEnergy + goNonnativeEnergy;

//fepb
      electEnergy_f = reduction->item(REDUCTION_ELECT_ENERGY_F);
      ljEnergy_f = reduction->item(REDUCTION_LJ_ENERGY_F);
      ljEnergy_f_left = reduction->item(REDUCTION_LJ_ENERGY_F_LEFT);

      electEnergy_ti_1 = reduction->item(REDUCTION_ELECT_ENERGY_TI_1);
      ljEnergy_ti_1 = reduction->item(REDUCTION_LJ_ENERGY_TI_1);
      electEnergy_ti_2 = reduction->item(REDUCTION_ELECT_ENERGY_TI_2);
      ljEnergy_ti_2 = reduction->item(REDUCTION_LJ_ENERGY_TI_2);
//fepe
    }

    if ( minimize || ! ( step % slowFreq ) )
    {
      electEnergySlow = reduction->item(REDUCTION_ELECT_ENERGY_SLOW);
//fepb
      electEnergySlow_f = reduction->item(REDUCTION_ELECT_ENERGY_SLOW_F);

      electEnergySlow_ti_1 = reduction->item(REDUCTION_ELECT_ENERGY_SLOW_TI_1);
      electEnergySlow_ti_2 = reduction->item(REDUCTION_ELECT_ENERGY_SLOW_TI_2);
      electEnergyPME_ti_1 = reduction->item(REDUCTION_ELECT_ENERGY_PME_TI_1);
      electEnergyPME_ti_2 = reduction->item(REDUCTION_ELECT_ENERGY_PME_TI_2);
//fepe
    }

    momentum.x = reduction->item(REDUCTION_MOMENTUM_X);
    momentum.y = reduction->item(REDUCTION_MOMENTUM_Y);
    momentum.z = reduction->item(REDUCTION_MOMENTUM_Z);
    angularMomentum.x = reduction->item(REDUCTION_ANGULAR_MOMENTUM_X);
    angularMomentum.y = reduction->item(REDUCTION_ANGULAR_MOMENTUM_Y);
    angularMomentum.z = reduction->item(REDUCTION_ANGULAR_MOMENTUM_Z);

    // Ported by JLai
    potentialEnergy = bondEnergy + angleEnergy + dihedralEnergy +
	improperEnergy + electEnergy + electEnergySlow + ljEnergy +
	crosstermEnergy + boundaryEnergy + miscEnergy + goTotalEnergy + groLJEnergy + groGaussEnergy;
    // End of port
    totalEnergy = potentialEnergy + kineticEnergy;
    flatEnergy = totalEnergy +
        (1.0/3.0)*( kineticEnergyHalfstep - kineticEnergyCentered);
    if ( !(step%slowFreq) ) {
      // only adjust based on most accurate energies
      BigReal s = (4.0/3.0)*( kineticEnergyHalfstep - kineticEnergyCentered);
      if ( smooth2_avg == XXXBIGREAL ) smooth2_avg = s;
      if ( step != simParams->firstTimestep ) {
        smooth2_avg *= 0.9375;
        smooth2_avg += 0.0625 * s;
      }
    }
    smoothEnergy = flatEnergy + smooth2_avg -
        (4.0/3.0)*( kineticEnergyHalfstep - kineticEnergyCentered);

    if ( simParameters->outputMomenta && ! minimize &&
         ! ( step % simParameters->outputMomenta ) )
    {
      iout << "MOMENTUM: " << step 
           << " P: " << momentum
           << " L: " << angularMomentum
           << "\n" << endi;
    }

    if ( simParameters->outputPressure ) {
      pressure_tavg += pressure;
      groupPressure_tavg += groupPressure;
      tavg_count += 1;
      if ( minimize || ! ( step % simParameters->outputPressure ) ) {
        iout << "PRESSURE: " << step << " "
           << PRESSUREFACTOR * pressure << "\n"
           << "GPRESSURE: " << step << " "
           << PRESSUREFACTOR * groupPressure << "\n";
        if ( tavg_count > 1 ) iout << "PRESSAVG: " << step << " "
           << (PRESSUREFACTOR/tavg_count) * pressure_tavg << "\n"
           << "GPRESSAVG: " << step << " "
           << (PRESSUREFACTOR/tavg_count) * groupPressure_tavg << "\n";
        iout << endi;
        pressure_tavg = 0;
        groupPressure_tavg = 0;
        tavg_count = 0;
      }
    }

    // pressure profile reductions
    if (pressureProfileSlabs) {
      const int freq = simParams->pressureProfileFreq;
      const int arraysize = 3*pressureProfileSlabs;
      
      BigReal *total = new BigReal[arraysize];
      memset(total, 0, arraysize*sizeof(BigReal));
      const int first = simParams->firstTimestep;

      if (ppbonded)    ppbonded->getData(first, step, lattice, total);
      if (ppnonbonded) ppnonbonded->getData(first, step, lattice, total);
      if (ppint)       ppint->getData(first, step, lattice, total);
      for (int i=0; i<arraysize; i++) pressureProfileAverage[i] += total[i];
      pressureProfileCount++;

      if (!(step % freq)) {
        // convert NAMD internal virial to pressure in units of bar 
        BigReal scalefac = PRESSUREFACTOR * pressureProfileSlabs;
   
        iout << "PRESSUREPROFILE: " << step << " ";
        if (step == first) {
          // output pressure profile for this step
          for (int i=0; i<arraysize; i++) {
            iout << total[i] * scalefac << " ";
          }
        } else {
          // output pressure profile averaged over the last count steps.
          scalefac /= pressureProfileCount;
          for (int i=0; i<arraysize; i++) 
            iout << pressureProfileAverage[i]*scalefac << " ";
        }
        iout << "\n" << endi; 
       
        // Clear the average for the next block
        memset(pressureProfileAverage, 0, arraysize*sizeof(BigReal));
        pressureProfileCount = 0;
      }
      delete [] total;
    }
  
   if ( step != simParams->firstTimestep || stepInFullRun == 0 ) {  // skip repeated first step
    if ( stepInFullRun % simParams->firstLdbStep == 0 ) {
     int benchPhase = stepInFullRun / simParams->firstLdbStep;
     if ( benchPhase > 0 && benchPhase < 10 ) {
#ifdef NAMD_CUDA
      if ( simParams->outputEnergies < 60 ) {
        iout << iWARN << "Energy evaluation is expensive, increase outputEnergies to improve performance.\n";
      }
#endif
      iout << iINFO;
      if ( benchPhase < 4 ) iout << "Initial time: ";
      else iout << "Benchmark time: ";
      iout << CkNumPes() << " CPUs ";
      {
        BigReal wallPerStep =
		(CmiWallTimer() - startBenchTime) / simParams->firstLdbStep;
	BigReal ns = simParams->dt / 1000000.0;
	BigReal days = 1.0 / (24.0 * 60.0 * 60.0);
	BigReal daysPerNano = wallPerStep * days / ns;
	iout << wallPerStep << " s/step " << daysPerNano << " days/ns ";
        iout << memusage_MB() << " MB memory\n" << endi;
      }
     }
     startBenchTime = CmiWallTimer();
    }
    ++stepInFullRun;
   }

    printTiming(step);

    Vector pairVDWForce, pairElectForce;
    if ( simParameters->pairInteractionOn ) {
      GET_VECTOR(pairVDWForce,reduction,REDUCTION_PAIR_VDW_FORCE);
      GET_VECTOR(pairElectForce,reduction,REDUCTION_PAIR_ELECT_FORCE);
    }

    // callback to Tcl with whatever we can
#ifdef NAMD_TCL
#define CALLBACKDATA(LABEL,VALUE) \
		labels << (LABEL) << " "; values << (VALUE) << " ";
#define CALLBACKLIST(LABEL,VALUE) \
		labels << (LABEL) << " "; values << "{" << (VALUE) << "} ";
    if (step == simParams->N && node->getScript() && node->getScript()->doCallback()) {
      std::ostringstream labels, values;
#if CMK_BLUEGENEL
      // the normal version below gives a compiler error
      values.precision(16);
#else
      values << std::setprecision(16);
#endif
      CALLBACKDATA("TS",step);
      CALLBACKDATA("BOND",bondEnergy);
      CALLBACKDATA("ANGLE",angleEnergy);
      CALLBACKDATA("DIHED",dihedralEnergy);
      CALLBACKDATA("CROSS",crosstermEnergy);
      CALLBACKDATA("IMPRP",improperEnergy);
      CALLBACKDATA("ELECT",electEnergy+electEnergySlow);
      CALLBACKDATA("VDW",ljEnergy);
      CALLBACKDATA("BOUNDARY",boundaryEnergy);
      CALLBACKDATA("MISC",miscEnergy);
      CALLBACKDATA("POTENTIAL",potentialEnergy);
      CALLBACKDATA("KINETIC",kineticEnergy);
      CALLBACKDATA("TOTAL",totalEnergy);
      CALLBACKDATA("TEMP",temperature);
      CALLBACKLIST("PRESSURE",pressure*PRESSUREFACTOR);
      CALLBACKLIST("GPRESSURE",groupPressure*PRESSUREFACTOR);
      CALLBACKDATA("VOLUME",lattice.volume());
      CALLBACKLIST("CELL_A",lattice.a());
      CALLBACKLIST("CELL_B",lattice.b());
      CALLBACKLIST("CELL_C",lattice.c());
      CALLBACKLIST("CELL_O",lattice.origin());
      labels << "PERIODIC "; values << "{" << lattice.a_p() << " "
		<< lattice.b_p() << " " << lattice.c_p() << "} ";
      if (simParameters->drudeOn) {
        CALLBACKDATA("DRUDEBOND",drudeBondTemp);
      }
      if ( simParameters->pairInteractionOn ) {
        CALLBACKLIST("VDW_FORCE",pairVDWForce);
        CALLBACKLIST("ELECT_FORCE",pairElectForce);
      }

      labels << '\0';  values << '\0';  // insane but makes Linux work
      state->callback_labelstring = labels.str();
      state->callback_valuestring = values.str();
      // node->getScript()->doCallback(labelstring.c_str(),valuestring.c_str());
    }
#undef CALLBACKDATA
#endif

    drudeBondTempAvg += drudeBondTemp;

    temp_avg += temperature;
    pressure_avg += trace(pressure)/3.;
    groupPressure_avg += trace(groupPressure)/3.;
    avg_count += 1;

    if ( simParams->outputPairlists && pairlistWarnings &&
				! (step % simParams->outputPairlists) ) {
      iout << iINFO << pairlistWarnings <<
        " pairlist warnings in past " << simParams->outputPairlists <<
	" steps.\n" << endi;
      pairlistWarnings = 0;
    }
    
    // NO CALCULATIONS OR REDUCTIONS BEYOND THIS POINT!!!
    if ( ! minimize &&  step % simParameters->outputEnergies ) return;
    // ONLY OUTPUT SHOULD OCCUR BELOW THIS LINE!!!

    if (simParameters->IMDon && !(step % simParameters->IMDfreq)) {
      IMDEnergies energies;
      energies.tstep = step;
      energies.T = temp_avg/avg_count;
      energies.Etot = totalEnergy;
      energies.Epot = potentialEnergy;
      energies.Evdw = ljEnergy;
      energies.Eelec = electEnergy + electEnergySlow;
      energies.Ebond = bondEnergy;
      energies.Eangle = angleEnergy;
      energies.Edihe = dihedralEnergy + crosstermEnergy;
      energies.Eimpr = improperEnergy;
      Node::Object()->imd->gather_energies(&energies);
    }

    if ( marginViolations ) {
      iout << iERROR << marginViolations <<
        " margin violations detected since previous energy output.\n" << endi;
    }
    marginViolations = 0;

    if ( (step % (10 * (minimize?1:simParameters->outputEnergies) ) ) == 0 )
    {
	iout << "ETITLE:      TS";
	iout << FORMAT("BOND");
	iout << FORMAT("ANGLE");
	iout << FORMAT("DIHED");
	if ( ! simParameters->mergeCrossterms ) iout << FORMAT("CROSS");
	iout << FORMAT("IMPRP");
        iout << "     ";
	iout << FORMAT("ELECT");
	iout << FORMAT("VDW");
	iout << FORMAT("BOUNDARY");
	iout << FORMAT("MISC");
	iout << FORMAT("KINETIC");
        iout << "     ";
	iout << FORMAT("TOTAL");
	iout << FORMAT("TEMP");
	iout << FORMAT("POTENTIAL");
	// iout << FORMAT("TOTAL2");
	iout << FORMAT("TOTAL3");
	iout << FORMAT("TEMPAVG");
	if ( volume != 0. ) {
          iout << "     ";
	  iout << FORMAT("PRESSURE");
	  iout << FORMAT("GPRESSURE");
	  iout << FORMAT("VOLUME");
	  iout << FORMAT("PRESSAVG");
	  iout << FORMAT("GPRESSAVG");
	}
        if (simParameters->drudeOn) {
          iout << "     ";
	  iout << FORMAT("DRUDEBOND");
	  iout << FORMAT("DRBONDAVG");
        }
	// Ported by JLai
	if (simParameters->goGroPair) {
	  iout << "     ";
	  iout << FORMAT("GRO_PAIR_LJ");
	  iout << FORMAT("GRO_PAIR_GAUSS");
	}

	if (simParameters->goForcesOn) {
          iout << "     ";
	  iout << FORMAT("NATIVE");
	  iout << FORMAT("NONNATIVE");
	  //iout << FORMAT("REL_NATIVE");
	  //iout << FORMAT("REL_NONNATIVE");
	  iout << FORMAT("GOTOTAL");
	  //iout << FORMAT("GOAVG");
	}
	// End of port -- JLai
	iout << "\n\n" << endi;
    }

    // N.B.  HP's aCC compiler merges FORMAT calls in the same expression.
    //       Need separate statements because data returned in static array.

    iout << ETITLE(step);
    iout << FORMAT(bondEnergy);
    iout << FORMAT(angleEnergy);
    if ( simParameters->mergeCrossterms ) {
      iout << FORMAT(dihedralEnergy+crosstermEnergy);
    } else {
      iout << FORMAT(dihedralEnergy);
      iout << FORMAT(crosstermEnergy);
    }
    iout << FORMAT(improperEnergy);
    iout << "     ";
    iout << FORMAT(electEnergy+electEnergySlow);
    iout << FORMAT(ljEnergy);
    iout << FORMAT(boundaryEnergy);
    iout << FORMAT(miscEnergy);
    iout << FORMAT(kineticEnergy);
    iout << "     ";
    iout << FORMAT(totalEnergy);
    iout << FORMAT(temperature);
    iout << FORMAT(potentialEnergy);
    // iout << FORMAT(flatEnergy);
    iout << FORMAT(smoothEnergy);
    iout << FORMAT(temp_avg/avg_count);
    if ( volume != 0. )
    {
        iout << "     ";
	iout << FORMAT(trace(pressure)*PRESSUREFACTOR/3.);
	iout << FORMAT(trace(groupPressure)*PRESSUREFACTOR/3.);
	iout << FORMAT(volume);
	iout << FORMAT(pressure_avg*PRESSUREFACTOR/avg_count);
	iout << FORMAT(groupPressure_avg*PRESSUREFACTOR/avg_count);
    }
    if (simParameters->drudeOn) {
        iout << "     ";
	iout << FORMAT(drudeBondTemp);
	iout << FORMAT(drudeBondTempAvg/avg_count);
    }
    // Ported by JLai
    if (simParameters->goGroPair) {
      iout << "     ";
      iout << FORMAT(groLJEnergy);
      iout << FORMAT(groGaussEnergy);
    }

    if (simParameters->goForcesOn) {
      iout << "     ";
      iout << FORMAT(goNativeEnergy);
      iout << FORMAT(goNonnativeEnergy);
      //iout << FORMAT(relgoNativeEnergy);
      //iout << FORMAT(relgoNonnativeEnergy);
      iout << FORMAT(goTotalEnergy);
      //iout << FORMAT("not implemented");
    } // End of port -- JLai
    iout << "\n\n" << endi;

#if(CMK_CCS_AVAILABLE && CMK_WEB_MODE)
     char webout[80];
     sprintf(webout,"%d %d %d %d",(int)totalEnergy,
	     (int)(potentialEnergy),
	     (int)kineticEnergy,(int)temperature);
     CApplicationDepositNode0Data(webout);
#endif

    if (simParameters->pairInteractionOn) {
      iout << "PAIR INTERACTION:";
      iout << " STEP: " << step;
      iout << " VDW_FORCE: ";
      iout << FORMAT(pairVDWForce.x);
      iout << FORMAT(pairVDWForce.y);
      iout << FORMAT(pairVDWForce.z);
      iout << " ELECT_FORCE: ";
      iout << FORMAT(pairElectForce.x);
      iout << FORMAT(pairElectForce.y);
      iout << FORMAT(pairElectForce.z);
      iout << "\n" << endi;
    }
    drudeBondTempAvg = 0;
    temp_avg = 0;
    pressure_avg = 0;
    groupPressure_avg = 0;
    avg_count = 0;

    if ( fflush_count ) {
      --fflush_count;
      fflush(stdout);
    }
}

void Controller::writeExtendedSystemLabels(ofstream_namd &file) {
  Lattice &lattice = state->lattice;
  file << "#$LABELS step";
  if ( lattice.a_p() ) file << " a_x a_y a_z";
  if ( lattice.b_p() ) file << " b_x b_y b_z";
  if ( lattice.c_p() ) file << " c_x c_y c_z";
  file << " o_x o_y o_z";
  if ( simParams->langevinPistonOn ) {
    file << " s_x s_y s_z s_u s_v s_w";
  }
  file << std::endl;
}

void Controller::writeExtendedSystemData(int step, ofstream_namd &file) {
  Lattice &lattice = state->lattice;
  file.precision(12);
  file << step;
    if ( lattice.a_p() ) file << " " << lattice.a().x << " " << lattice.a().y << " " << lattice.a().z;
    if ( lattice.b_p() ) file << " " << lattice.b().x << " " << lattice.b().y << " " << lattice.b().z;
    if ( lattice.c_p() ) file << " " << lattice.c().x << " " << lattice.c().y << " " << lattice.c().z;
    file << " " << lattice.origin().x << " " << lattice.origin().y << " " << lattice.origin().z;
  if ( simParams->langevinPistonOn ) {
    Vector strainRate = diagonal(langevinPiston_strainRate);
    Vector strainRate2 = off_diagonal(langevinPiston_strainRate);
    file << " " << strainRate.x;
    file << " " << strainRate.y;
    file << " " << strainRate.z;
    file << " " << strainRate2.x;
    file << " " << strainRate2.y;
    file << " " << strainRate2.z;
  }
  file << std::endl;
}

void Controller::enqueueCollections(int timestep)
{
  if ( Output::coordinateNeeded(timestep) ){
    collection->enqueuePositions(timestep,state->lattice);
  }
  if ( Output::velocityNeeded(timestep) )
    collection->enqueueVelocities(timestep);
  if ( Output::forceNeeded(timestep) )
    collection->enqueueForces(timestep);
}

//Modifications for alchemical fep
static char *FEPTITLE(int X)
{ 
  static char tmp_string[21];
  sprintf(tmp_string, "FepEnergy: %6d ",X);
  return tmp_string;
}

static char *TITITLE(int X)
{ 
  static char tmp_string[21];
  sprintf(tmp_string, "TI   %6d ",X);
  return tmp_string;
}

void Controller::outputFepEnergy(int step) {
 if (simParams->alchFepOn) {
  const int stepInRun = step - simParams->firstTimestep;
  const int alchEquilSteps = simParams->alchEquilSteps;
  const BigReal alchLambda = simParams->alchLambda;
  const BigReal alchLambda2 = simParams->alchLambda2;
  const bool alchEnsembleAvg = simParams->alchEnsembleAvg; 
  const bool FepWhamOn = simParams->alchFepWhamOn;

  if (alchEnsembleAvg && (stepInRun == 0 || stepInRun == alchEquilSteps)) {
    FepNo = 0;
    exp_dE_ByRT = 0.0;
    net_dE = 0.0;
  }
  BigReal dE = electEnergy_f + electEnergySlow_f + ljEnergy_f
		- (electEnergy + electEnergySlow + ljEnergy);
  BigReal RT = BOLTZMANN * simParams->alchTemp;

  if (alchEnsembleAvg){
  FepNo++;
  exp_dE_ByRT += exp(-dE/RT);
  net_dE += dE;
  }
 
  if (stepInRun == 0) {
    if (!fepFile.is_open()) {
      NAMD_backup_file(simParams->alchOutFile);
      fepFile.open(simParams->alchOutFile);
      iout << "OPENING FEP ENERGY OUTPUT FILE\n" << endi;
      if(alchEnsembleAvg){
      fepSum = 0.0;
      fepFile << "#            STEP                 Elec                            "
              << "vdW                    dE           dE_avg         Temp             dG\n"
              << "#                           l             l+dl      "
              << "       l            l+dl         E(l+dl)-E(l)" << std::endl;
      }
      else{
					if(!FepWhamOn){ 
            fepFile << "#            STEP                 Elec                            "
                    << "vdW                    dE         Temp\n"
                    << "#                           l             l+dl      "
                    << "       l            l+dl         E(l+dl)-E(l)" << std::endl;
          } 
      }
    }
    if(!step){
  		if(!FepWhamOn){ 
        fepFile << "#NEW FEP WINDOW: "
                << "LAMBDA SET TO " << alchLambda << " LAMBDA2 " 
                << alchLambda2 << std::endl;
      }
    }
  }
  if ((alchEquilSteps) && (stepInRun == alchEquilSteps)) {
    fepFile << "#" << alchEquilSteps << " STEPS OF EQUILIBRATION AT "
            << "LAMBDA " << simParams->alchLambda << " COMPLETED\n"
            << "#STARTING COLLECTION OF ENSEMBLE AVERAGE" << std::endl;
  }
  if ((simParams->N) && (simParams->alchOutFreq && ((step%simParams->alchOutFreq)==0))) {
    writeFepEnergyData(step, fepFile);
    fepFile.flush();
  }
  if (alchEnsembleAvg && (step == simParams->N)) {
    fepSum = fepSum + dG;
    fepFile << "#Free energy change for lambda window [ " << alchLambda
	    << " " << alchLambda2 << " ] is " << dG << " ; net change until now is " << fepSum << std::endl;
    fepFile.flush();
  }
 }
}

void Controller::outputTiEnergy(int step) {
 if (simParams->alchThermIntOn) {
  const int stepInRun = step - simParams->firstTimestep;
  const int alchEquilSteps = simParams->alchEquilSteps;
  const BigReal alchLambda = simParams->alchLambda;
  
  if (stepInRun == 0 || stepInRun == alchEquilSteps) {
    TiNo = 0;
    net_dEdl_elec_1 = 0;
    net_dEdl_elec_2 = 0;
    net_dEdl_lj_1 = 0;
    net_dEdl_lj_2 = 0;
  }
  if (stepInRun == 0 || (! ((step - 1) % simParams->alchOutFreq))) {
    // output of instantaneous dU/dl now replaced with running average
    // over last alchOutFreq steps (except for step 0)
    recent_TiNo = 0;
    recent_dEdl_elec_1 = 0;
    recent_dEdl_elec_2 = 0;
    recent_dEdl_lj_1 = 0;
    recent_dEdl_lj_2 = 0;
  }
  TiNo++;
  recent_TiNo++;
  // FB - PME is no longer scaled by global lambda, but by the respective
  // lambda as dictated by elecLambdaStart. All electrostatics now go together.
  net_dEdl_elec_1 += electEnergy_ti_1 + electEnergySlow_ti_1 + electEnergyPME_ti_1;
  net_dEdl_elec_2 += electEnergy_ti_2 + electEnergySlow_ti_2 + electEnergyPME_ti_2;
  net_dEdl_lj_1 += ljEnergy_ti_1;
  net_dEdl_lj_2 += ljEnergy_ti_2;
  recent_dEdl_elec_1 += electEnergy_ti_1 + electEnergySlow_ti_1 + electEnergyPME_ti_1; 
  recent_dEdl_elec_2 += electEnergy_ti_2 + electEnergySlow_ti_2 + electEnergyPME_ti_2; 
  recent_dEdl_lj_1 += ljEnergy_ti_1;
  recent_dEdl_lj_2 += ljEnergy_ti_2;

  if (stepInRun == 0) {
    BigReal alchElecLambdaStart = simParams->alchElecLambdaStart;
    BigReal alchVdwLambdaEnd = simParams->alchVdwLambdaEnd;
    BigReal elec_lambda_1 = (alchLambda <= alchElecLambdaStart)? 0. : \
            (alchLambda - alchElecLambdaStart) / (1. - alchElecLambdaStart);
    BigReal elec_lambda_2 = ((1-alchLambda) <= alchElecLambdaStart)? 0. : \
            ((1-alchLambda) - alchElecLambdaStart) / (1. - alchElecLambdaStart);
    BigReal vdw_lambda_1 =  (alchLambda >= alchVdwLambdaEnd)? 1. : \
            alchLambda / alchVdwLambdaEnd; 
    BigReal vdw_lambda_2 =  ((1-alchLambda) >= alchVdwLambdaEnd)? 1. : \
            (1-alchLambda) / alchVdwLambdaEnd; 
    if (!tiFile.is_open()) {
      //tiSum = 0.0;
      NAMD_backup_file(simParams->alchOutFile);
      tiFile.open(simParams->alchOutFile);
      iout << "OPENING TI ENERGY OUTPUT FILE\n" << endi;
      tiFile << "#       STEP      Elec_dU/dl      Elec_avg        vdW_dU/dl      vdw_avg       Elec_dU/dl      Elec_avg      vdW_dU/dl       vdw_avg       PME_dU/dl      PME_avg\n"
              << "#               <---------------------PARTITION 1------------------------>    <---------------------PARTITION 2--------------------->" 
              << std::endl;
    }
    tiFile << "#NEW TI WINDOW: "
            << "LAMBDA " << alchLambda 
            << "\n#PARTITION 1 VDW LAMBDA " << vdw_lambda_1 
            << "\n#PARTITION 1 ELEC LAMBDA " << elec_lambda_1 
            << "\n#PARTITION 2 VDW LAMBDA " << vdw_lambda_2 
            << "\n#PARTITION 2 ELEC LAMBDA " << elec_lambda_2 
            << "\n" << std::endl;
  }
  if (stepInRun == alchEquilSteps) {
    tiFile << "#" << alchEquilSteps << " STEPS OF EQUILIBRATION AT "
            << "LAMBDA " << simParams->alchLambda << " COMPLETED\n"
            << "#STARTING COLLECTION OF ENSEMBLE AVERAGE" << std::endl;
  }
  if (simParams->alchOutFreq && ((step%simParams->alchOutFreq)==0)) {
    writeTiEnergyData(step, tiFile);
    tiFile.flush();
  }
 }
}

void Controller::writeFepEnergyData(int step, ofstream_namd &file) {
  BigReal eeng = electEnergy+electEnergySlow;
  BigReal eeng_f = electEnergy_f + electEnergySlow_f;
  BigReal dE = eeng_f + ljEnergy_f - eeng - ljEnergy;
  BigReal dE_Left = eeng_f + ljEnergy_f_left - eeng - ljEnergy;
  BigReal RT = BOLTZMANN * simParams->alchTemp;
  const bool alchEnsembleAvg = simParams->alchEnsembleAvg;
  const int stepInRun = step - simParams->firstTimestep;
	const bool FepWhamOn = simParams->alchFepWhamOn;
	const bool WCARepuOn = simParams->alchFepWCARepuOn;
	const BigReal WCArcut1 = simParams->alchFepWCArcut1;
	const BigReal WCArcut2 = simParams->alchFepWCArcut2;
	const BigReal WCArcut3 = simParams->alchFepWCArcut3;
	const BigReal alchLambda = simParams->alchLambda;
	
	const BigReal alchRepLambda = simParams->alchRepLambda;
	const BigReal alchDispLambda = simParams->alchDispLambda;
	const BigReal alchElecLambda = simParams->alchElecLambda;
  if(stepInRun){
		if(!FepWhamOn){
			fepFile << FEPTITLE(step);
			fepFile << FORMAT(eeng);
			fepFile << FORMAT(eeng_f);
			fepFile << FORMAT(ljEnergy);
			fepFile << FORMAT(ljEnergy_f);
		}
		else{	// FepWhamOn = ON
			if(WCARepuOn){
				if(WCArcut1<WCArcut2) {	// [s1,s2]
 					fepFile << "FEP_WCA_REP  ";
					fepFile << FORMAT(WCArcut1);
					fepFile << FORMAT(WCArcut2);
					fepFile << FORMAT(1.0);
		      fepFile << FORMAT(dE_Left);
				}
				if(WCArcut2<WCArcut3) {	// [s2,s3]
					if(WCArcut1<WCArcut2) fepFile << " BREAK ";
 					fepFile << "FEP_WCA_REP  ";
					fepFile << FORMAT(WCArcut2);
					fepFile << FORMAT(WCArcut3);
					fepFile << FORMAT(0.0);
		      fepFile << FORMAT(dE);
				}
        fepFile << std::endl;
			}
			else if(simParams->alchFepWCADispOn)	{
				fepFile << "FEP_WCA_DISP ";
				fepFile << FORMAT(alchDispLambda);
			}
			else if(simParams->alchFepElecOn)	{
				fepFile << "FEP_ELEC     ";
				fepFile << FORMAT(alchElecLambda);
			}
		}
		if( ! WCARepuOn ) {
      fepFile << FORMAT(dE);
		}
    if(alchEnsembleAvg){
      BigReal dE_avg = net_dE/FepNo;
      fepFile << FORMAT(dE_avg);
    }
		if(!FepWhamOn){ 
			fepFile << FORMAT(temperature);
		}
    if(alchEnsembleAvg){
      dG = -(RT * log(exp_dE_ByRT/FepNo));
      fepFile << FORMAT(dG);
    } 
		if( ! WCARepuOn ) {
      fepFile << std::endl;
		}
  }
}
//fepe
void Controller::writeTiEnergyData(int step, ofstream_namd &file) {
  tiFile << TITITLE(step);
  tiFile << FORMAT(recent_dEdl_elec_1 / recent_TiNo);
  tiFile << FORMAT(net_dEdl_elec_1/TiNo);
  tiFile << FORMAT(recent_dEdl_lj_1 / recent_TiNo);
  tiFile << FORMAT(net_dEdl_lj_1/TiNo);
  tiFile << FORMAT(recent_dEdl_elec_2 / recent_TiNo);
  tiFile << FORMAT(net_dEdl_elec_2/TiNo);
  tiFile << FORMAT(recent_dEdl_lj_2 / recent_TiNo);
  tiFile << FORMAT(net_dEdl_lj_2/TiNo);
  tiFile << std::endl;
}

void Controller::outputExtendedSystem(int step)
{

  if ( step >= 0 ) {

    // Write out eXtended System Trajectory (XST) file
    if ( simParams->xstFrequency &&
         ((step % simParams->xstFrequency) == 0) )
    {
      if ( ! xstFile.is_open() )
      {
        iout << "OPENING EXTENDED SYSTEM TRAJECTORY FILE\n" << endi;
        NAMD_backup_file(simParams->xstFilename);
        xstFile.open(simParams->xstFilename);
        while (!xstFile) {
          if ( errno == EINTR ) {
            CkPrintf("Warning: Interrupted system call opening XST trajectory file, retrying.\n");
            xstFile.clear();
            xstFile.open(simParams->xstFilename);
            continue;
          }
          char err_msg[257];
          sprintf(err_msg, "Error opening XST trajectory file %s",simParams->xstFilename);
          NAMD_err(err_msg);
        }
        xstFile << "# NAMD extended system trajectory file" << std::endl;
        writeExtendedSystemLabels(xstFile);
      }
      writeExtendedSystemData(step,xstFile);
      xstFile.flush();
    }

    // Write out eXtended System Configuration (XSC) files
    //  Output a restart file
    if ( simParams->restartFrequency &&
         ((step % simParams->restartFrequency) == 0) &&
         (step != simParams->firstTimestep) )
    {
      iout << "WRITING EXTENDED SYSTEM TO RESTART FILE AT STEP "
		<< step << "\n" << endi;
      char fname[140];
      const char *bsuffix = ".old";
      strcpy(fname, simParams->restartFilename);
      if ( simParams->restartSave ) {
        char timestepstr[20];
        sprintf(timestepstr,".%d",step);
        strcat(fname, timestepstr);
        bsuffix = ".BAK";
      }
      strcat(fname, ".xsc");
      NAMD_backup_file(fname,bsuffix);
      ofstream_namd xscFile(fname);
      while (!xscFile) {
        if ( errno == EINTR ) {
          CkPrintf("Warning: Interrupted system call opening XSC restart file, retrying.\n");
          xscFile.clear();
          xscFile.open(fname);
          continue;
        }
        char err_msg[257];
        sprintf(err_msg, "Error opening XSC restart file %s",fname);
        NAMD_err(err_msg);
      } 
      xscFile << "# NAMD extended system configuration restart file" << std::endl;
      writeExtendedSystemLabels(xscFile);
      writeExtendedSystemData(step,xscFile);
      if (!xscFile) {
        char err_msg[257];
        sprintf(err_msg, "Error writing XSC restart file %s",fname);
        NAMD_err(err_msg);
      } 
    }

  }

  //  Output final coordinates
  if (step == FILE_OUTPUT || step == END_OF_RUN)
  {
    int realstep = ( step == FILE_OUTPUT ?
        simParams->firstTimestep : simParams->N );
    iout << "WRITING EXTENDED SYSTEM TO OUTPUT FILE AT STEP "
		<< realstep << "\n" << endi;
    static char fname[140];
    strcpy(fname, simParams->outputFilename);
    strcat(fname, ".xsc");
    NAMD_backup_file(fname);
    ofstream_namd xscFile(fname);
    while (!xscFile) {
      if ( errno == EINTR ) {
        CkPrintf("Warning: Interrupted system call opening XSC output file, retrying.\n");
        xscFile.clear();
        xscFile.open(fname);
        continue;
      }
      char err_msg[257];
      sprintf(err_msg, "Error opening XSC output file %s",fname);
      NAMD_err(err_msg);
    } 
    xscFile << "# NAMD extended system configuration output file" << std::endl;
    writeExtendedSystemLabels(xscFile);
    writeExtendedSystemData(realstep,xscFile);
    if (!xscFile) {
      char err_msg[257];
      sprintf(err_msg, "Error writing XSC output file %s",fname);
      NAMD_err(err_msg);
    } 
  }

  //  Close trajectory file
  if (step == END_OF_RUN) {
    if ( xstFile.is_open() ) {
      xstFile.close();
      iout << "CLOSING EXTENDED SYSTEM TRAJECTORY FILE\n" << endi;
    }
  }

}


void Controller::recvCheckpointReq(const char *key, int task, checkpoint &cp) {  // responding replica
  switch ( task ) {
    case SCRIPT_CHECKPOINT_STORE:
      if ( ! checkpoints.count(key) ) {
        checkpoints[key] = new checkpoint;
      }
      *checkpoints[key] = cp;
      break;
    case SCRIPT_CHECKPOINT_LOAD:
      if ( ! checkpoints.count(key) ) {
        NAMD_die("Unable to load checkpoint, requested key was never stored.");
      }
      cp = *checkpoints[key];
      break;
    case SCRIPT_CHECKPOINT_SWAP:
      if ( ! checkpoints.count(key) ) {
        NAMD_die("Unable to swap checkpoint, requested key was never stored.");
      }
      std::swap(cp,*checkpoints[key]);
      break;
    case SCRIPT_CHECKPOINT_FREE:
      if ( ! checkpoints.count(key) ) {
        NAMD_die("Unable to free checkpoint, requested key was never stored.");
      }
      delete checkpoints[key];
      checkpoints.erase(key);
      break;
  }
}

void Controller::recvCheckpointAck(checkpoint &cp) {  // initiating replica
  if ( checkpoint_task == SCRIPT_CHECKPOINT_LOAD || checkpoint_task == SCRIPT_CHECKPOINT_SWAP ) {
    state->lattice = cp.lattice;
    *(ControllerState*)this = cp.state;
  }
  CkpvAccess(_qd)->process();
}


void Controller::rebalanceLoad(int step)
{
  if ( ! ldbSteps ) { 
    ldbSteps = LdbCoordinator::Object()->getNumStepsToRun();
  }
  if ( ! --ldbSteps ) {
    //CkPrintf("Controller: rebalancing %d, thread %p\n", step, CthSelf());
    startBenchTime -= CmiWallTimer();
	Node::Object()->outputPatchComputeMaps("before_ldb", step);
    LdbCoordinator::Object()->rebalance(this);	
	startBenchTime += CmiWallTimer();
    fflush_count = 3;
  }
}

void Controller::cycleBarrier(int doBarrier, int step) {
#if USE_BARRIER
	if (doBarrier) {
	  broadcast->cycleBarrier.publish(step,1);
	  CkPrintf("Cycle time at sync Wall: %f CPU %f\n",
		  CmiWallTimer()-firstWTime,CmiTimer()-firstCTime);
	}
#endif
}

void Controller::traceBarrier(int turnOnTrace, int step) {
	CkPrintf("Cycle time at trace sync (begin) Wall at step %d: %f CPU %f\n", step, CmiWallTimer()-firstWTime,CmiTimer()-firstCTime);	
	CProxy_Node nd(CkpvAccess(BOCclass_group).node);
	nd.traceBarrier(turnOnTrace, step);
	CthSuspend();
}

void Controller::resumeAfterTraceBarrier(int step){
	broadcast->traceBarrier.publish(step,1);
	CkPrintf("Cycle time at trace sync (end) Wall at step %d: %f CPU %f\n", step, CmiWallTimer()-firstWTime,CmiTimer()-firstCTime);	
	awaken();
}

#ifdef MEASURE_NAMD_WITH_PAPI
void Controller::papiMeasureBarrier(int turnOnMeasure, int step){
	CkPrintf("Cycle time at PAPI measurement sync (begin) Wall at step %d: %f CPU %f\n", step, CmiWallTimer()-firstWTime,CmiTimer()-firstCTime);	
	CProxy_Node nd(CkpvAccess(BOCclass_group).node);
	nd.papiMeasureBarrier(turnOnMeasure, step);
	CthSuspend();
}

void Controller::resumeAfterPapiMeasureBarrier(int step){
	broadcast->papiMeasureBarrier.publish(step,1);
	CkPrintf("Cycle time at PAPI measurement sync (end) Wall at step %d: %f CPU %f\n", step, CmiWallTimer()-firstWTime,CmiTimer()-firstCTime);	
	awaken();
}
#endif

void Controller::terminate(void) {
  BackEnd::awaken();
  CthFree(thread);
  CthSuspend();
}

