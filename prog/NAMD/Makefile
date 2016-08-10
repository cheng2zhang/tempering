
dir_orig   = ~/lwork/NAMD_2.11_Source/src
dir_minmod = ~/lwork/NAMD_mods/NAMD_2.11_tempering/src
dir_spec   = ~/lwork/NAMD_mods/NAMD_2.11_spec/src
dir_thstat = ~/lwork/NAMD_mods/NAMD_2.11_thstat/src

all: minmod.patch spec.patch thstat.patch

minmod.patch::
	rm -f $@ minmod_*/*
	ln -f `diff -wdBqr --exclude="*~" --exclude=".*" $(dir_orig) $(dir_minmod) | cut -f2 -d" "` minmod_old/
	ln -f `diff -wdBqr --exclude="*~" --exclude=".*" $(dir_orig) $(dir_minmod) | cut -f4 -d" "` minmod_new/
	git add minmod_*/*
	diff -u -wdB minmod_old minmod_new > $@; [ $$? -eq 1 ]

spec.patch::
	rm -f $@ spec_*/* spec_*/.*~
	ln -f `diff -wdBqr --exclude="*~" --exclude=".*" $(dir_orig) $(dir_spec) | cut -f2 -d" "` spec_old/
	ln -f `diff -wdBqr --exclude="*~" --exclude=".*" $(dir_orig) $(dir_spec) | cut -f4 -d" "` spec_new/
	git add spec_*/*
	diff -u -wdB spec_old spec_new > $@; [ $$? -eq 1 ]

thstat.patch::
	rm -f $@ thstat_*/* thstat_*/.*~
	ln -f `diff -wdBqr --exclude="*~" --exclude=".*" $(dir_orig) $(dir_thstat) | cut -f2 -d" "` thstat_old/
	ln -f `diff -wdBqr --exclude="*~" --exclude=".*" $(dir_orig) $(dir_thstat) | cut -f4 -d" "` thstat_new/
	git add thstat_*/*
	diff -u -wdB thstat_old thstat_new > $@; [ $$? -eq 1 ]

NAMD_thstat.zip::
	$(RM) -f $@
	$(MAKE) -B thstat.patch
	zip -r $@ thstat* *.html \
	  README.thstat.md \
	  doc/vsTmove.tex doc/*.pdf \
	  doc/fig/*.gp doc/fig/*.eps doc/fig/*.pdf \
	  test/Argon_NAMD_ST/*.prm \
	  test/Argon_NAMD_ST/*.conf \
	  test/Argon_NAMD_ST/argon_300K.coor \
	  test/Argon_NAMD_ST/*.psf \
	  test/Argon_NAMD_ST/cmp.png \
	  test/Argon_NAMD_ST/ke.png \
	  test/waterbox/*.conf \
	  test/waterbox/*.inp \
	  test/waterbox/waterbox.* \
	  ./mkhist.py \
	  test/alanin/*.conf

NAMD_spec.zip::
	$(RM) -f $@
	$(MAKE) -B spec.patch
	zip -r $@ spec* *.html \
	  README.thstat.md README.spec.md \
	  doc/vsTmove.tex doc/*.pdf \
	  doc/fig/*.gp doc/fig/*.eps doc/fig/*.pdf \
	  test/Argon_NAMD_ST/*.prm \
	  test/Argon_NAMD_ST/*.conf \
	  test/Argon_NAMD_ST/argon_300K.coor \
	  test/Argon_NAMD_ST/*.psf \
	  test/Argon_NAMD_ST/cmp.png \
	  test/Argon_NAMD_ST/ke.png \
	  test/waterbox/*.conf \
	  test/waterbox/*.inp \
	  test/waterbox/waterbox.* \
	  ./mkhist.py \
	  test/alanin/*.conf

Bossman:
	make NAMD_thstat.zip
	make NAMD_spec.zip
	rsync -avz NAMD_*.zip /Bossman/cz1/tempering/ 
