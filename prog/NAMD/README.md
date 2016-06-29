# NAMD files


Variances     | Description
--------------|---------------------
minmod        | minimal corrections to the original version
spec          | also report end-to-end distance of C-alpha atoms in every step
thstat        | Langevin-style velocity rescaling and Nose-Hoover thermostats

## Variances

### thstat

This patch contains several modifications to NAMD 2.11

 1. The correction of the velocity-scaling problem (due to Justin) after temperature transitions in adaptive tempering.
 2. Fixing (hopefully) the one-step mismatch in Sequencer.C and Controller.C.
 3. Fixing the bin index overflow problem (due to Justin) in adaptTempUpdate() in Controller.C.
 4. Disabling the code for adaptTempRandom, which lacks theoretically foundation.
 5. Issuing a warning for using adaptive tempering with the original velocity rescaling, which does not sample the Boltzmann distribution.
 6. Integrating the Langevin-style velocity rescaling thermostat and Nose-Hoover thermostat.
 7. Monitoring the distribution of the (reduced) kinetic energy.


#### Velocity-scaling after temperature transition

Since the temperature is a dynamic variable in adaptive tempering, 
the velocities must be scaled corresponding after a temperature transition.
The original implementation lacks this step, and thus is incorrect.

#### One-step mismatch in Controller.C and Sequencer.C

The original implementation places adaptTempUpdate() at the beginning of an MD step in Controller.C,
but the function of the same name at the end of an MD step in Sequencer.C.
This causes the temperature differs by an MD step in Controller.C and Sequencer.C.
To fix the problem, we moved the adaptTempUpdate() in Controller.C to the end of an MD step.

However, this fix sometimes causes a fatal error:
as the function broadcasts the adaptive temperature to Sequencers,
if the sequencers receive the new temperature after the rebalanceLoad() call in Controller.C,
the program will crash.
To fix the problem, we further ask all Sequencers to send a Hi message back to the Controller
upon receiving the new adaptive temperature,
and the Controller must collect the Hi messages from all Sequencers before calling rebalanceLoad().
This quick fix appears to prevent the above mentioned problem,
although the underlying reason is still unclear.
The code of receiving and sending Hi messages are implemented in CollectionMaster and CollectionMgr modules,
e.g., CollectionMgr::submitHi() and CollectionMaster::enqueueHi().


#### Possible index overflow and underflow

In adaptTempUpdate(), the temperature index adaptTempBin (line 1997) may be out of range.
We fix this by correct the index to the 0 or adaptTempBins - 1.
This should work if the overflow or underflow is due to rounding error.

#### Disabling the code for adaptTempRandom

If the Langevin equation drives the temperature out of range,
the scheme by adaptTempRandom will randomly pick a new temperature in the range.
This is dangerous and theoretically-unfounded strategy.
Our fix is to simply abandon the move and keep the old adaptive temperature,
just as one would do in a failed Monte Carlo move.

#### Issuing a warning for using adaptive tempering with the original velocity rescaling

The orginal velocity rescaling does not sample a canonical distribution.
Therefore, it should not be used with adaptive tempering for production runs.

#### Integrating the Langevin-style velocity rescaling thermostat and Nose-Hoover thermostat

Currently, the only rigorous thermostat for sampling a Boltzmann distribution is the Langevin dynamics.
However, thermostats based on uniformly scaling the velocity are usually more efficient.
We implement two examples of this type of global scaling thermostats.

* The velocity rescaling thermostat [Bussi, Donadio, and Parrinello, JCP 126, 014101 (2007)].
* The Nose-Hoover chain thermostat [Martyna, Klein, and Tuckerman, JCP 97, 2635 (1992)].

The former can be understood as the stochastic version of the latter,
and it has the advantage of no need of saving the state of the chain variables.

To use the velocity rescaling thermostat, one need to add these lines to the configuration file:
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
the mass is proportional to kB * T times `tNHCPeriod` squared divided by 4 Pi^2.
The mass of the first chain variable is further multiplied by the number of degree of freedom.
The masses of the rest of the chain variables are assumed to be the same.
The parameter `tNHCFile` is the file recording the state of the chain variables.
This is a text file of three lines.
The first line gives the number of chain variables.
The second line records the velocities of the chain variables.
The last line records the masses of the chain variables
(which are usually ignored on loading unless `tNHCFileReadMass` is turned on).
The parameter `tNHCFileFreq` specifies how often we should save the state of chain variable.
The parameter `tNHCFileReadMass` is used to instruct NAMD read mass from the `tNHCFile`.
This allows the user to specify masses of the chain variables explicitly.

Caution.  Thermostats (including the Langevin dynamics and the old velocity rescaling) are exclusive,
in each simulation, only one thermostat can be turned on.


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
is multipled by a factor of the reference temperature / the instantaneous temperature.
The instantaneous temperature is defined as twice the kinetic energy
divided by the number of degrees of freedom divided by the Boltzmann constant.


## Apply patches

http://www.thegeekstuff.com/2014/12/patch-command-examples/

```
make minmod.patch
```

To use the patch
```
patch -b -p3 < minmod.patch
