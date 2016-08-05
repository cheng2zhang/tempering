# NAMD files


### thstat

This patch contains several modifications to NAMD 2.11:

  * The correction of the velocity-scaling problem (due to Justin) after temperature transitions in adaptive tempering.
  * Fixing (hopefully) the one-step mismatch problem in Sequencer.C and Controller.C.
  * Fixing the bin index overflow problem (due to Justin) in adaptTempUpdate() in Controller.C.
  * Properly overwriting (instead of appending) the restart file (due to Justin).
  * Computing the potential energy at the beginning adaptTempUpdate() (due to Justin).
  * Implementing the separate accumulator scheme.
  * Implementing the Monte Carlo scheme for temperature change.
  * Disabling the code for the ad hoc adaptTempRandom scheme, which lacks theoretically foundation.
  * Issuing a warning for using adaptive tempering with the original velocity rescaling, which does not rigorously sample the Boltzmann distribution.
  * Miscellaneous modifications.
  * Integrating the Langevin-style velocity-rescaling thermostat and Nose-Hoover thermostat.
  * Adaptively rescaling the velocity to approach an asymptotic microcanonical ensemble.
  * Monitoring the distribution of the (reduced) kinetic energy.
  * Logging the potential energy.


#### Velocity-scaling after temperature transition

Since the temperature is a dynamic variable in adaptive tempering,
the velocities must be scaled corresponding after a temperature transition.
The original implementation lacks this step, and thus is incorrect.
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
and the Controller must wait for the Hi messages from all Sequencers before calling rebalanceLoad().
This quick fix appears to prevent the above mentioned problem,
although the underlying reason is still unclear.
The code of receiving and sending Hi messages are implemented in the CollectionMaster and CollectionMgr modules, respectively.
Please see, e.g., CollectionMgr::submitHi() and CollectionMaster::enqueueHi().


#### Possible index overflow and underflow

In adaptTempUpdate(), the temperature index adaptTempBin (line 1997) may be out of range.
We fix this by correcting the index to 0 or adaptTempBins - 1.
This should work if the overflow or underflow is due to rounding error.
But the user should still make sure that the thermostat temperature lying within the range at the beginning.

#### Properly overwriting the the restart file

Currently the restart file is appended instead of overwritten because of a programming mistake.
This is fixed in the patch.  If the appending behavior is desired,
the user can set the option `adaptTempRestartAppend`.

#### Recomputing the potential energy in adaptTempUpdate()

The old code deduce the potential energy from the conservation of the total energy.
But that is inaccurate because the MD integrator only approximately conserves the total energy.
The code now recompute the potential energy at the beginning of adaptTempUpdate().
We also attempt to fix a possible bug of how often the LJcorrection is computed.

#### Implementing the separate accumulator scheme

To make the adaptive averaging scheme work most efficiently,
each bin need a separator accumulator associated with the window.
This feature is now implemented.  To use it, set
```
adaptTempSep    on
```

#### Implementing a Monte Carlo scheme for temperature transition

In addition to the Langevin equation, adaptive tempering can now be done through a Monte Carlo scheme.
To use this feature, set
```
adaptTempMCMove    on
```

#### Disabling the code for adaptTempRandom

If the Langevin equation drives the temperature out of range,
the scheme by adaptTempRandom will randomly pick a new temperature in the range.
This is an ad hoc and theoretically-unfounded strategy.
Our fix is to simply abandon the invalid temperature transition and keep the old adaptive temperature,
just as one would do in a failed Monte Carlo move.

#### Miscellaneous modifications

  * Allowing adaptTempInFile and adaptTempBins to be set simultaneously, the former overrides the latter.
  * Adding the inverse temperature as the first column of the restart file.
  * Using the average energy computed from the integral identity as the average energy in the restart file (the second column).
  * Throwing out an exception when reading from the restart file fails.
  * Adding the option `adaptTempFixedAve` to the fix the average energies from the input restart file (due to Justin).

#### Issuing a warning for using adaptive tempering with the original velocity rescaling

The orginal velocity rescaling does not sample a canonical distribution.
Therefore, it should not be used with adaptive tempering for production runs.

#### Integrating the Langevin-style velocity rescaling thermostat and Nose-Hoover thermostat

Currently, the only rigorous thermostat for sampling a Boltzmann distribution is the Langevin dynamics.
However, thermostats based on uniformly scaling the velocity are usually more efficient.
Please see the demo, thermostat_ljdemo.html, for a comparison of the autocorrelation functions
of the kinetic energy and velocity.
We implement two examples of this type of global scaling thermostats.

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

#### Adaptively rescale the velocity to approach an asymptotic microcanonical ensemble

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
But this can be too demanding on the precision.
A workaround is to set `rescaleAdaptiveDedk` to a number greater than 1.0,
which is roughly the change of the total energy divided by the
change of the kinetic energy in response to a temperature change
```
rescaleAdaptiveDedk 1.5
```

#### Monitoring the distribution of the (reduced) kinetic energy

To monitor the distribution of the kinetic energy,
we further supply the following options
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

For a regular simulation, the logging outputs step and the potential energy.
For a simulation using adaptive tempering, the logging also outputs the temperature as the last column.


## Apply patches

http://www.thegeekstuff.com/2014/12/patch-command-examples/

```
make thstat.patch
```

To use the patch
```
patch -b -p3 < thstat.patch
