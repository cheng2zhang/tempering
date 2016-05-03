Justin's 100 argon test case.

Use `make clean` to clean up running files.

The energy should be multiplied by a factor.


## Comparison

```
/path/to/original/namd2 +p2 argon_namd_st.conf
tail -n40 argon_st1.rst > novscale.txt
make clean
/path/to/modified/namd2 +p2 argon_namd_st.conf
tail -n40 argon_st1.rst > vscale.txt
make cmp.png
```

