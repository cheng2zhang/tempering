Justin's 100 argon test case.

## Comparison of using velocity scaling after the temperature transition

Automatically:
```
make cmp.png
```
Optionally, use -B to regenerate `novscale.txt` and `vscale.txt`.
```
make -B cmp.png
```

Manually

```
path/to/original/namd2 +p3 argon_st.conf
tail -n60 argon_st1.rst > novscale.txt
make clean
path/to/modified/namd2 +p3 argon_st.conf
tail -n60 argon_st1.rst > vscale.txt
make cmp.png
```

If Langevin dynamics is not used, then it appears to be fine.
Probably because velocity scaling is a very strong (but perhaps inexact?) thermostat.


## Tips

Use `make clean` to clean up running files.

## About the units

### From reduced to real

Energies should be multiplied by a factor of eps = 0.238.
Temperatures by a factor of kB/eps = 0.001987191/0.238 = 0.00835.
Length by a factor of sigma = 3.405.

