/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef CONTROLLER_H
#define CONTROLLER_H

#include "converse.h"
#include "Node.h"
#include "common.h"
#include "fstream_namd.h"
#include <string>
#include <map>

class ControllerBroadcasts;
class NamdState;
class SimParameters;
class RequireReduction;
class SubmitReduction;

#ifdef MEM_OPT_VERSION
class CollectionMasterHandler;
#else
class CollectionMaster;
#endif

class Random;
class PressureProfileReduction;

struct ControllerState {
    Tensor langevinPiston_strainRate;
    Tensor berendsenPressure_avg;
    int berendsenPressure_count;
    BigReal smooth2_avg;
};

class Controller : protected ControllerState
{
public:
    Controller(NamdState *s);
    virtual ~Controller(void);
    void run(void);             // spawn thread, etc.
    void awaken(void) { CthAwaken(thread); };
    void resumeAfterTraceBarrier(int);
#ifdef MEASURE_NAMD_WITH_PAPI
	void resumeAfterPapiMeasureBarrier(int step);
#endif

protected:
    friend class ScriptTcl;
    friend class Node;
    friend class CheckpointMsg;
    virtual void algorithm(void);	// subclasses redefine this method

    void integrate(int); // Verlet integrator
    void minimize(); // CG minimizer
      RequireReduction *min_reduction;

    void receivePressure(int step, int minimize = 0);
      Tensor pressure_normal;
      Tensor pressure_nbond;
      Tensor pressure_slow;
      Tensor virial_amd;
      Tensor groupPressure_normal;
      Tensor groupPressure_nbond;
      Tensor groupPressure_slow;
      Tensor controlPressure_normal;
      Tensor controlPressure_nbond;
      Tensor controlPressure_slow;
      int nbondFreq;
      int slowFreq;
      BigReal temp_avg;
      BigReal pressure_avg;
      BigReal groupPressure_avg;
      int avg_count;
      Tensor pressure_tavg;
      Tensor groupPressure_tavg;
      int tavg_count;
    void compareChecksums(int,int=0);
      int computeChecksum;
      int marginViolations;
      int pairlistWarnings;
    void printTiming(int);
    void printMinimizeEnergies(int);
      BigReal min_energy;
      BigReal min_f_dot_f;
      BigReal min_f_dot_v;
      BigReal min_v_dot_v;
      int min_huge_count;
    void printDynamicsEnergies(int);
    void printEnergies(int step, int minimize);
      int numDegFreedom;
      int stepInFullRun;
      BigReal totalEnergy;
      BigReal electEnergy;
      BigReal electEnergySlow;
      BigReal ljEnergy;
      BigReal groLJEnergy;
      BigReal groGaussEnergy;
      BigReal goNativeEnergy;
      BigReal goNonnativeEnergy;
      BigReal goTotalEnergy;
//fepb
      BigReal electEnergy_f;
      BigReal electEnergySlow_f;
      BigReal ljEnergy_f;
      BigReal ljEnergy_f_left;	// used by WCA repulsive, [s1,s2]
      BigReal exp_dE_ByRT;
      BigReal net_dE;
      BigReal dG;
      int FepNo;
      void printFepMessage(int);
      BigReal fepSum;
//fepe

      BigReal electEnergy_ti_1;
      BigReal electEnergySlow_ti_1;
      BigReal ljEnergy_ti_1;
      BigReal electEnergy_ti_2;
      BigReal electEnergySlow_ti_2;
      BigReal ljEnergy_ti_2;
      BigReal net_dEdl_elec_1;
      BigReal net_dEdl_elec_2;
      BigReal net_dEdl_lj_1;
      BigReal net_dEdl_lj_2;
      BigReal electEnergyPME_ti_1;
      BigReal electEnergyPME_ti_2;
      int TiNo;
      BigReal recent_dEdl_elec_1;
      BigReal recent_dEdl_elec_2;
      BigReal recent_dEdl_lj_1;
      BigReal recent_dEdl_lj_2;
      int recent_TiNo;
      void printTiMessage(int);

      BigReal drudeBondTemp; // temperature of Drude bonds
      BigReal drudeBondTempAvg;

      BigReal potentialEnergy;
      BigReal kineticEnergy;
      BigReal kineticEnergyHalfstep;
      BigReal kineticEnergyCentered;
      BigReal temperature;
      // BigReal smooth2_avg;
      BigReal smooth2_avg2;  // avoid internal compiler error
      Tensor pressure;
      Tensor groupPressure;
      int controlNumDegFreedom;
      Tensor controlPressure;
    void enqueueCollections(int);
    void correctMomentum(int step);
    void rescaleForTotalEnergy(void);
    void rescaleVelocities(int);
    void rescaleVelocitiesInit(void);
    void rescaleVelocitiesLoad(void);
    void rescaleVelocitiesSave(int);
      BigReal rescaleVelocities_sumTemps;
      int rescaleVelocities_numTemps;
      // the following quantities are used to compute the block average of beta
      BigReal rescaleVelocities_count;
      BigReal rescaleVelocities_sbeta;
      // the following quantities are used to compute beta'(E)
      BigReal rescaleVelocities_sum1;
      BigReal rescaleVelocities_sumBeta;
      BigReal rescaleVelocities_sumBeta2;
      BigReal rescaleVelocities_sumDbde;
    void reassignVelocities(int);
    void tcoupleVelocities(int);
    void langRescaleVelocities(int, Bool);
    BigReal langRescaleFactorPrev;
    void tNHCInit(void);
    void tNHCDone(int);
    void tNHCRescaleVelocities(int, Bool);
    void tNHCSave(int);
    void tNHCLoad(void);
    BigReal tNHCRescaleFactorPrev;
    BigReal *tNHCzeta;
    BigReal *tNHCmass;
    void keHistInit(void);
    void keHistDone(int);
    void keHistUpdate(int);
    // save the kinetic energy to file
    // the first column is the kinetic energy
    // the second and third columns are the
    // normalized histogram and the reference value
    void keHistSave(int);
    void keHistLoad(void);
    BigReal *keHist;
    int keHistBinMax;
    void berendsenPressure(int);
      // Tensor berendsenPressure_avg;
      // int berendsenPressure_count;
    void langevinPiston1(int);
    void langevinPiston2(int);
      Tensor langevinPiston_origStrainRate;
      Tensor strainRate_old;  // for langevinPistonBarrier no
      Tensor positionRescaleFactor;  // for langevinPistonBarrier no

    int ldbSteps;
    void rebalanceLoad(int);
      int fflush_count;
    void cycleBarrier(int,int);	
	
	void traceBarrier(int, int);

#ifdef MEASURE_NAMD_WITH_PAPI
	void papiMeasureBarrier(int, int);
#endif

    // void suspend(void) { CthSuspend(); };
    void terminate(void);

    Random *random;
    SimParameters *const simParams;	// for convenience
    NamdState *const state;		// access data in state
    RequireReduction *reduction;
    RequireReduction *amd_reduction;
    SubmitReduction *submit_reduction;

    // data for pressure profile reductions and output
    PressureProfileReduction *ppbonded;
    PressureProfileReduction *ppnonbonded;
    PressureProfileReduction *ppint;
    int pressureProfileSlabs;
    int pressureProfileCount;
    BigReal *pressureProfileAverage;

    CollectionMaster *const collection;
    
    ControllerBroadcasts * broadcast;
    ofstream_namd xstFile;
    void outputExtendedSystem(int step);
    void writeExtendedSystemLabels(ofstream_namd &file);
    void writeExtendedSystemData(int step, ofstream_namd &file);

//fepb
    ofstream_namd fepFile;
    void outputFepEnergy(int step);
    void writeFepEnergyData(int step, ofstream_namd &file);
//fepe
    ofstream_namd tiFile;
    void outputTiEnergy(int step);
    void writeTiEnergyData(int step, ofstream_namd &file);

    // for checkpoint/revert
    int checkpoint_stored;
    Lattice checkpoint_lattice;
    ControllerState checkpoint_state;

    struct checkpoint {
      Lattice lattice;
      ControllerState state;
    };
    std::map<std::string,checkpoint*> checkpoints;
    int checkpoint_task;
    void recvCheckpointReq(const char *key, int task, checkpoint &cp);
    void recvCheckpointAck(checkpoint &cp);

    Lattice origLattice;

//for accelMD
   void rescaleaccelMD (int step, int minimize = 0);
   BigReal accelMDdVAverage;

//JS for adaptive temperature sampling
   void adaptTempInit(int step);
   void adaptTempDone(int step);
   BigReal adaptTempGetInvW(BigReal tp);
   BigReal adaptTempRegression(void);
   BigReal adaptTempGetPEAve(int i, BigReal def = 0, BigReal beta = 0);
   BigReal adaptTempGetIntE(BigReal beta, int i, BigReal nbeta, int ni);
   BigReal adaptTempMCMove(BigReal tp, BigReal ep);
   BigReal adaptTempLangevin(BigReal tp, BigReal ep);
   Bool adaptTempUpdate(int step, int minimize = 0);
   void adaptTempWriteRestart(int step);
   int *adaptTempBinMinus;
   int *adaptTempBinPlus;
   // separator accumulator
   struct AdaptTempSepAcc {
     int bin0; // first bin
     int winSize;
     double *sumw;
     double *sumE;
     double *sumE2;
     double *ave;
     double *var;
     double *cnt;
     double total;
     double invGamma;

     AdaptTempSepAcc(void) {
       winSize = 0;
     }

     ~AdaptTempSepAcc(void) {
       if ( winSize > 0 ) {
         delete [] sumw;
         delete [] sumE;
         delete [] sumE2;
         delete [] ave;
         delete [] var;
         delete [] cnt;
       }
     }

     void empty(void) {
       for ( int j = 0; j < winSize; j++ ) {
         sumw[j]  = 0;
         sumE[j]  = 0;
         sumE2[j] = 0;
         ave[j]   = 0;
         var[j]   = 0;
         cnt[j]   = 0;
       }
     }

     // initialize the window
     void init(int minus, int plus) {
       bin0 = minus;
       winSize = plus - minus;
       if ( winSize % 2 == 0 ) NAMD_die("Window size should be odd");
       total = 0;
       invGamma = 1;
       sumw  = new double[winSize];
       sumE  = new double[winSize];
       sumE2 = new double[winSize];
       ave   = new double[winSize];
       var   = new double[winSize];
       cnt   = new double[winSize];
       empty();
     }

     // compute the average and variance of each bin
     void trim(void) {
       if ( invGamma > 1.0 ) { // renormalize the running weight
         double gam = 1.0/invGamma;
         for ( int j = 0; j < winSize; j++ ) {
           sumw[j]  *= gam;
           sumE[j]  *= gam;
           sumE2[j] *= gam;
         }
         invGamma = 1;
       }
       
       for ( int j = 0; j < winSize; j++ ) {
         if ( sumw[j] > 0 ) {
           ave[j] = sumE[j] / sumw[j];
           var[j] = sumE2[j] / sumw[j] - ave[j] * ave[j];
         } else {
           ave[j] = 0;
           var[j] = 0;
         }
         //CkPrintf("j %d, %g %g %g %g %g\n", j, sumw[j], sumE[j], sumE2[j], ave[j], var[j]);
       }
     }

     // add a data point from bin i to this accumulator
     void add(int i, BigReal potEne, BigReal invw, BigReal cg) {
       i -= bin0; // convert to the local index
       if ( i < 0 || i >= winSize ) {
         CkPrintf("Bad local index for bin0 %d: i %d, winSize %d\n", bin0, i + bin0, winSize);
         NAMD_die("Adaptive tempering: bad local index.");
       }
       total += 1;
       double gamma = 1 - cg / total;
       if ( gamma < 1e-8 ) gamma = 1e-8;
       invGamma /= gamma;
       invw *= invGamma;
       sumw[i]  += invw;
       sumE[i]  += invw * potEne;
       sumE2[i] += invw * potEne * potEne;
       cnt[i]   += 1;
       //CkPrintf("adding to bin %d with potEne %g\n", i, potEne);
     }

     // compute the average energy from the integral identity
     BigReal iiave(BigReal varCntMin, BigReal def = 0) {
       int j, mid = winSize / 2;
       double den0 = 0.0, den1 = 0.0, ene0 = 0.0, ene1 = 0.0;
       double A0 = 0, A1 = 0, A2 = 0, cntMax = 0, defVar = 0, varj;
       if ( total <= 0 ) return def;
       trim();
       // compute the default variance from the most populated bin
       for ( j = 0; j < winSize; j++ ) {
         if ( sumw[j] > cntMax ) {
           cntMax = sumw[j];
           defVar = var[j];
         }
       }
       if ( defVar <= 0 ) defVar = 1.0;
       // left side
       for ( j = 0; j <= mid; j++ ) {
         den0 += sumw[j];
         ene0 += sumE[j];
         varj = (cnt[j] > varCntMin) ? var[j] : defVar;
         A0 += varj * (j + 0.5);
         //CkPrintf("+ %d %g %g %g\n", j, sumw[j], sumE[j], invw[j]);
       }
       den0 /= mid + 1;
       ene0 /= mid + 1;
       A0 /= mid + 1;
       // middle bin correction
       varj = (cnt[mid] > varCntMin) ? var[mid] : defVar;
       A2 = 0.5 * varj;
       // right side
       for ( j = mid + 1; j < winSize; j++ ) {
         den1 += sumw[j];
         ene1 += sumE[j];
         varj = (cnt[j] > varCntMin) ? var[j] : defVar;
         A1 += varj * (j - winSize + 0.5);
         //CkPrintf("+ %d %g %g %g\n", j, sumw[j], sumE[j], invw[j]);
       }
       if ( mid > 0 ) {
         den1 /= mid;
         ene1 /= mid;
         A1 /= mid;
       }
       if ( den0 + den1 <= 0 ) return def;
       // compute a+ and a-
       double aplus = ( den0 + den1 > 0 ) ? (A0 - A2) / (A0 - A1) : 0;
       if ( aplus < 0 ) aplus = 0;
       if ( aplus > 1 ) aplus = 1;
       double aminus = 1 - aplus;
       //CkPrintf("A0 %g, A1 %g, A2 %g, a- %g, a+ %g, ave %g, %g, %g\n", A0, A1, A2, aminus, aplus, sum0/(den0 + 1e-16), sum1/(den1 + 1e-16), (aminus * sum0 + aplus * sum1)/(aminus * den0 + aplus * den1)); getchar();
       return ( aminus * ene0 + aplus * ene1 )
            / ( aminus * den0 + aplus * den1 );
     }
   };
   AdaptTempSepAcc *adaptTempSepAcc;
   double  adaptTempMCSize;
   double  adaptTempMCTot, adaptTempMCAcc;
   double  adaptTempMCDAcc, adaptTempMCFail; // accumulators for adjusting the MC size
   double  adaptTempLangTot, adaptTempLangAcc;
   double  adaptTempLangDAcc, adaptTempLangFail; // accumulators for adjusting Dt
   double  *adaptTempPotEnergyAveNum;
   double  *adaptTempPotEnergyAveDen;
   double  *adaptTempPotEnergyVarNum;
   double  *adaptTempPotEnergyAve;
   double  *adaptTempPotEnergyVar;
   long    *adaptTempPotEnergySamples;
   BigReal *adaptTempBetaN;
   BigReal adaptTempT;
   BigReal adaptTempBetaMin;
   BigReal adaptTempBetaMax;
   int     adaptTempBin;
   int     adaptTempBins;
   BigReal adaptTempDBeta;
   BigReal adaptTempCg;
   BigReal adaptTempDt;
   Bool    adaptTempAnaDirty;
   BigReal adaptTempAnaSlope;
   BigReal adaptTempAnaIntercept;
   ofstream_namd adaptTempRestartFile;
  
    // special atoms
    void specInit(int scriptTask, int step);

private:
    CthThread thread;
    static void threadRun(Controller*);

    double startCTime;
    double startWTime;
    double firstCTime;
    double firstWTime;
    double startBenchTime;

    int computesPartitioned;
};

#endif // CONTROLLER_H

