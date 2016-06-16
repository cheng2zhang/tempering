
dir_orig   = ~/lwork/NAMD_2.11_Source/src
dir_minmod = ~/lwork/NAMD_mods/NAMD_2.11_tempering/src
dir_spec   = ~/lwork/NAMD_mods/NAMD_2.11_spec/src

minmod.patch::
	ln -f `diff -wdBqr --exclude="*~" --exclude=".*" $(dir_orig) $(dir_minmod) | cut -f2 -d" "` minmod_old/
	ln -f `diff -wdBqr --exclude="*~" --exclude=".*" $(dir_orig) $(dir_minmod) | cut -f4 -d" "` minmod_new/
	git add minmod_*/*
	diff -u -wdB minmod_old minmod_new > $@; [ $$? -eq 1 ]

