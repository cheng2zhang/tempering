# NAMD files

This file can be viewed online

https://github.com/cheng2zhang/tempering/blob/master/prog/NAMD/README.thstat.md

### thstat

This patch contains several modifications to NAMD 2.11:

 + Modifications related to the adaptive tempering module
  * Velocities rescaling after temperature transitions in adaptive tempering (due to Justin).
  * Fixing (hopefully) the one-step mismatch problem in Sequencer.C and Controller.C.
  * Computing the potential energy at the beginning `adaptTempUpdate()` (due to Justin).
  * Holding back temperature transitions until a certain number of samples per bin (due to Justin).
  * Number-of-visits-weighted integral identity.
  * Implementing the Monte Carlo scheme for temperature transitions.
  * Implementing an MC-corrected Langevin equation integration scheme.
  * Implementing the separate accumulator scheme for adaptive averaging.
  * Fine tuning of the overall temperature distribution.
  * Features and options for the restart file.
  * Disabling the code for `adaptTempRandom` and `adaptTempAutoDt`.
  * Miscellaneous modifications.

 + New thermostats
  * Integrating the Langevin-style velocity-rescaling thermostat and Nose-Hoover thermostat.
  * Adaptively rescaling the velocity to approach an asymptotic microcanonical ensemble.

 + Logging
  * Monitoring the distribution of the (reduced) kinetic energy.
  * Logging the potential energy.


#### Velocity rescaling after temperature transitions

Since the temperature is a dynamic variable in adaptive tempering,
the velocities must be rescaled correspondingly after a temperature transition.
Simulated tempering based on the potential energy does not work without this step.
Please see details in doc/vsTmove.pdf, and ../test/Argon_NAMD_ST/cmp.png for an example.

#### One-step mismatch in Controller.C and Sequencer.C

The original implementation places adaptTempUpdate() at the beginning of an MD step in Controller.C,
but the function of the same name at the end of an MD step in Sequencer.C.
This causes the temperature to differ by an MD step in Controller.C and Sequencer.C.
To fix the problem, we move the adaptTempUpdate() in Controller.C to the end of an MD step.

However, this causes a fatal error in the run-time:
as the function broadcasts the adaptive temperature to Sequencers,
if the Sequencers receive the new temperature after the rebalanceLoad() call in Controller.C,
the program will crash.
To fix the problem, we further ask all Sequencers to send a "Hi" message back to the Controller
upon receiving the new adaptive temperature,
and the Controller must collect the Hi messages from all Sequencers before calling rebalanceLoad().
This quick fix appears to prevent the above mentioned problem,
although the underlying reason is still unclear.
The code of receiving and sending Hi messages are implemented
in the CollectionMaster and CollectionMgr modules, respectively.
Please see, e.g., CollectionMgr::submitHi() and CollectionMaster::enqueueHi().


#### Recomputing the potential energy in adaptTempUpdate()

Previously, the potential energy was deduced from the conservation of the total energy
by subtracting the kinetic energy from the total energy.
But that is inaccurate because the MD integrator only approximately
conserves the total energy.
The code now recomputes the potential energy at the beginning of adaptTempUpdate().

#### Holding back temperature transitions until a certain number of samples per bin

In an initial stage, the temperature often drifts to one end of the temperature spectrum.
This is due to lack of equilibration and bad initial estimate of the average energy.
To alleviate this problem, we can hold back temperature transitions until the number of
samples within the current bin exceeds the number given by `adaptTempSamplesMin`
```
adaptTempSamplesMin     2000
```

#### Number-of-visits-weighted integral identity

The integral identity used for computing the average energy can be modified
to use the number of visits as a weighting factor.
This modification hopefully improve the stability during the equilibration stage.

#### Implementing a Monte Carlo scheme for temperature transitions

In addition to the Langevin equation approach,
adaptive tempering can now be done through direct Monte Carlo (MC).
To use this feature, set
```
adaptTempMCMove    on
adaptTempMCSize    0.01
```
The move size, specified by `adaptTempMCSize` is given as a fraction of the current temperature.
The value should be adjusted such that the acceptance ratio, as `ACC. RATIO` printed out on the screen, is roughly 50%.
To automate the adjusting process during equilibration, one can set the target acceptance ratio
```
adaptTempMCAutoAR   0.5
```
The adjustment occurs in every step, but with decreasing magnitude (~ 1/t) to approach a fixed asymptotic limit.
However, for safety, the automatic adjusting feature should only be used in equilibration.

In order to compute how the acceptance ratio changes with the MC move size,
the program needs to virtually increase the MC move size by a little amount specified by
```
adaptTempMCSizeInc  0.0005
```
This value is shared in the MC-corrected Langevin equation

#### Implementing an MC-corrected Langevin equation integration scheme

For the Langevin equation, the value recommended in the paper (JCP, 2010, 132, 244101)
0.0001 appears to too large in that it broadens the potential energy distribution (due to Justin).
The Langevin equation is now corrected by a Monte Carlo scheme
so it is correct no matter the size of time step.
However, the acceptance ratio of the correction step should be greater than 0.5.
To automatically adjust the move size during equilibration,
one can set the target acceptance ratio as
```
adaptTempDtAutoAR   0.5
```
The adjustment occurs in every step, but with decreasing magnitude (~ 1/t).
However, for safety, the option should only be used in equilibration.

The approximate equivalence relationship between `adaptTempDt` and `adaptTempMCSize` is
```
adaptTempDt ~ adaptTempMCSize^2 / 2
```
With the new correction scheme, the Langevin equation may offer a slightly larger
move size than the Monte Carlo scheme.

#### Implementing the separate accumulator scheme

To make the adaptive averaging scheme work most efficiently,
each bin needs a separator accumulator associated with the window.
In this way the parameter Cgamma applies to all bins within the window.
To use the feature, set
```
adaptTempSep    on
```
We usually use an `adaptTempCgamma` around 0.1 with such a scheme, during the beginning of simulation.
Later on, one may choose to set `adaptTempCgamma` to 0.0, or fixing the weight with `adaptTempFixedAve`.
When this feature is turned on, the data for the separate accumulators
will be appended to the restart file.

#### Fine-tuning of the overall temperature distribution

The overall temperature distribution can now be tuned by the new parameter `adaptTempWeightExp`.
In terms of the distribution of the inverse-temperature, beta, this parameter corresponds to x as in

  w(beta) ~ 1/beta^x.

Equivalently, the distribution of temperature T, is proportional to T^(x - 2).
Thus, to achieve a flat-T histogram, we need to set x = 2.
To achieve a flat-beta histogram, we need to set x = 0.
The default value is 1.0, which corresponds to a flat-lnT histogram. 
```
adaptTempWeightExp    1
```

#### Feature and options for the restart file

When writting the restart file, it is now overwritten instead of appended.
If the appending behavior is needed, the user can set the option `adaptTempRestartAppend`.

Several options are added to the output restart file.
  * Adding the inverse temperature as the first column of the restart file (due to Justin).
  * Using the average energy computed from the integral identity as the average energy in the restart file (the second column).
  * Adding the inverse weight to the last column of the restart file (column 8).
  * Increasing the precision of the restart file.
The product of column 4 (histogram) and column 8 (inverse weight) should be roughly a constant after a long run.

Two options are added to use the input restart file. 
The option `adaptTempFixedAve` can fix the average energies read from the input restart file (due to Justin).
During the simulation, the number of visits and other data will still be accumulated,
but the second column of the output restart file regarding the average energy will be fixed.

The option `adaptTempEmptyData` can empty data after loading the input restart file.
This feature can be used in conjugation with `adaptTempFixedAve` to see if
the average energy from the input restart file is able to produce the desired temperature distribution.

The program now throws out an exception when there is a failure reading the input restart file.

#### Disabling the code for `adaptTempRandom` and `adaptTempAutoDt`

The option `adaptTempRandom` is deprecated.
If the Langevin equation drives the temperature out of range,
the scheme triggered by the option `adaptTempRandom` will randomly pick a new temperature in the range.
This strategy is not exact.
In the new version, an invalid temperature transition is abandoned and
the old adaptive temperature is kept,
just as one would do in a failed Monte Carlo move.

The option `adaptTempAutoDt` is deprecated.
A similar functionality is provided by `adaptTempDtAutoAR` as discussed above.

A warning message will be displayed if either option is used.

#### Miscellaneous modifications

  * Adding the option `adaptTempWindowSize` to adjust the window size of integral identity (due to Justin).
  * Allowing `adaptTempInFile` and `adaptTempBins` to be set simultaneously in the configuration file, the former overrides the latter.
  * Reducing the default window size for the integral identity from 0.04 to 0.02.
  * Trying to change how often the `LJcorrection` is computed.
  * Langevin time step adaptTempDt is dimensionless (previously it assumes the unit of femto-second, SimParameters.h).
  * Issuing a warning for using adaptive tempering with the original velocity rescaling scheme. The original velocity rescaling does not sample a canonical distribution, and should not be for use with adaptive tempering for production runs.  As a replacement, use the Langevin-style velocity rescaling thermostat instead.

#### Integrating the Langevin-style velocity rescaling thermostat and Nose-Hoover thermostat

Currently, the only rigorous thermostat for sampling a Boltzmann distribution is the Langevin dynamics.
However, thermostats based on uniformly scaling the velocity are usually more efficient.
Please see the demo, thermostat_ljdemo.html, for a comparison of the autocorrelation functions
of the kinetic energy and velocity.
We implement two global thermostats based on velocity rescaling.

* The velocity-rescaling thermostat [Bussi, Donadio, and Parrinello, JCP 126, 014101 (2007)].
* The Nose-Hoover chain thermostat [Martyna, Klein, and Tuckerman, JCP 97, 2635 (1992)].

The former can be understood as the stochastic version of the latter,
and it has the advantage of no need of saving the state of the chain variables.

To use the velocity-rescaling thermostat, one need to add these lines to the configuration file:
```
langRescale          on
langRescaleTemp      300
langRescaleDt        100.0
```
The parameter `langRescaleDt` is the inverse viscosity in femtoseconds:
in each MD step, the Langevin equation for the total kinetic energy
will be integrated by an effective time step of `timestep / langRescaleDt`,
where `timestep` is the MD step size in femtoseconds (usually 1.0 or 2.0).

To use the Nose-Hoover chain thermostat, one need to add these lines to the configuration file:
```
tNHC                  on
tNHCTemp              300
tNHCLen               5
tNHCPeriod            100.0
tNHCFile              mysimul.nhc
tNHCFileFreq          1000
tNHCFileReadMass      off
```
The parameter `tNHCLen` is the number of chain variables in the Nose-Hoover chains,
and when `tNHCLen = 1` it recovers the original Nose-Hoover thermostat.
The parameter `tNHCPeriod` determines the masses of the chain variables:
the mass is proportional to kB T times `tNHCPeriod` squared divided by 4 Pi^2.
The mass of the first chain variable is further multiplied by the number of degree of freedom.
The masses of the rest of the chain variables are assumed to be the same.
The parameter `tNHCFile` is the file recording the state of the chain variables.
This is a text file of three lines.
The first line gives the number of chain variables, and the MD time step.
The second line records the velocities of the chain variables.
The last line records the masses of the chain variables
(which are usually ignored on loading unless `tNHCFileReadMass` is turned on).
This file, if exists, is automatically loaded at the beginning of a simulation.
If a fresh restart is needed, please delete the file before running NAMD.
The parameter `tNHCFileFreq` specifies how often we should save the state of chain variable.
The parameter `tNHCFileReadMass` is used to instruct NAMD to read masses from `tNHCFile`.
This allows the user to specify the masses of the chain variables explicitly.

Caution.  Thermostats (including the Langevin dynamics and the old velocity rescaling) are exclusive,
in each simulation, only one thermostat can be turned on.

#### Adaptively rescaling the velocity to approach an asymptotic microcanonical ensemble

To better control temperature in the microcanonical ensemble,
we can use the native NAMD mechanism of velocity rescaling,
but gradually reduce the scaling magnitude.
To enable this feature, turn on the adaptive velocity scaling feature (`rescaleAdaptive`)
```
rescaleTemp               300
rescaleFreq               10
rescaleAdaptive           on
rescaleAdaptiveFile       adaptvrescale.dat
rescaleAdaptiveFileFreq   100000
```
In this way, the magnitude of velocity scaling is modified by a factor of 1/t,
where t is the number of times of such scaling so far.
The file specified by `rescaleAdaptiveFile` is automatically reloaded.
In initial runs please delete this file.

By default, we use the exact method to compute dbeta/dE.
But this can be too demanding on the precision,
a workaround is to set `rescaleAdaptiveDedk` to a number greater than 1.0,
which is roughly the change of the total energy divided by the
change of the kinetic energy in response to a temperature change
```
rescaleAdaptiveDedk 1.5
```

#### Monitoring the distribution of the (reduced) kinetic energy

To monitor the distribution of the kinetic energy,
we supply the following options:
```
keHist               on
keHistBin            1.0
keHistFile           ke.dat
keHistFileFreq       1000
```
The parameter `keHistBin` is the bin size for the kinetic energy (the unit is kcal/mol).
The parameter `keHistFile` gives the output histogram file, which is a text file.
The first column is the kinetic energy.
The second column is the normalized observed histogram of the kinetic energy.
The third column is the reference value.
So in gnuplot, we can plot the file using
```
plot "ke.dat" u 1:2 w l t "Observed", "ke.dat" u 1:3 w l t "Reference"
```
This file is saved every `keHistFileFreq` MD steps,
and it is automatically reloaded upon resuming a simulation.
Thus, to disable reloading the previous file, please delete this file before running NAMD.

For adaptive tempering, the temperature is a variable, so the kinetic energy
is multiplied by a factor of (the thermostat reference temperature / the instantaneous temperature).
The instantaneous temperature is defined as twice the kinetic energy
divided by the number of degrees of freedom divided by the Boltzmann constant.
Please see ../test/Argon_NAMD_ST/ke.png for an example.

#### Logging the potential energy

To log the potential energy, set `energyLogFile`.
```
energyLogFile      ene.log
energyLogFreq      1
```
By default the frequency of logging the potential energy is 1.

For a regular simulation, the logging outputs the MD step and the potential energy.
For a simulation using adaptive tempering, the logging also outputs the temperature as the last column.

Once we have the log file, the histogram can be constructed using the Python script `mkhist.py`:
```
./mkhist.py --dE=5 --dT=0.5 -T300 ene1.log
```
The command will produce an output `ene1.his` accordingly,
with the energy bin size being 5 kcal/mol,
and in the adaptive tempering case, collecting frame whose temperature is between 299.5K to 300.5K.
If the input log file is omitted, all `ene*.log` files under the current directory will be processed.


## Applying patches

The patch file is produced by using `diff` on the `thstat_old` and `thstat_new` directory.
See the tutorial in http://www.thegeekstuff.com/2014/12/patch-command-examples/
A convenient `make` command is provided:
```
make thstat.patch
```

To use the patch
```
patch -b -p3 < thstat.patch
```

