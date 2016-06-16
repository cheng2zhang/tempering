# NAMD files


Subdirectory  | Description
--------------|---------------------
original      | Original files from the NAMD source on 2016/05/02
minmod        | minimal corrections to the original version


## Patches

http://www.thegeekstuff.com/2014/12/patch-command-examples/

```
make minmod.patch
```

To use the patch
```
patch -b -p3 < minmod.patch
