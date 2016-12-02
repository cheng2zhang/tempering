
dir_orig   = ~/lwork/NAMD_2.11_Source/src
dir_minmod = ~/lwork/NAMD_mods/NAMD_2.11_tempering/src
dir_spec   = ~/lwork/NAMD_mods/NAMD_2.11_spec/src
dir_thstat = ~/lwork/NAMD_mods/NAMD_2.11_thstat/src

all: minmod.patch spec.patch thstat.patch

minmod.patch::
	rm -f $@ minmod_*/* minmod_*/.*~
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

testfiles = test/alanin/*.conf \
	    test/alanin/*.param \
	    test/alanin/alanin.pdb \
	    test/alanin/alanin.psf \
	    test/argon/*.prm \
	    test/argon/*.conf \
	    test/argon/argon_300K.coor \
	    test/argon/*.psf \
	    test/argon/cmp.png \
	    test/argon/ke.png \
	    test/argon/Makefile \
	    test/argon/*.gp \
	    test/waterbox/*.conf \
	    test/waterbox/*.inp \
	    test/waterbox/waterbox.* \
	    test/waterbox/Makefile \
	    test/waterbox/*.gp \
	    test/butane/*.conf \
	    test/butane/*ethers* \
	    test/butane/butane.* \
	    test/butane/Makefile \
	    test/butane/*.gp \
	    test/gly5/*.conf \
	    test/gly5/par*.inp \
	    test/gly5/gly5.psf \
	    test/gly5/gly5.pdb \
	    test/gly5/Makefile \
	    test/gly5/*.gp \
	    test/gly35/*.conf \
	    test/gly35/gly35.psf \
	    test/gly35/gly35.pdb \
	    test/gly35/Makefile

NAMD_thstat.zip::
	$(RM) -f $@
	$(MAKE) -B thstat.patch
	zip -r $@ thstat* *.html \
	  README.thstat.md \
	  doc/vsTmove.tex doc/*.pdf \
	  doc/fig/*.gp doc/fig/*.eps doc/fig/*.pdf \
	  $(testfiles) \
	  ./mkhist.py ./calcspec.py

NAMD_spec.zip::
	$(RM) -f $@
	$(MAKE) -B spec.patch
	zip -r $@ spec* *.html \
	  README.thstat.md README.spec.md \
	  doc/vsTmove.tex doc/*.pdf \
	  doc/fig/*.gp doc/fig/*.eps doc/fig/*.pdf \
	  $(testfiles) \
	  ./mkhist.py ./calcspec.py

Bossman:
	make NAMD_thstat.zip
	make NAMD_spec.zip
	rsync -avz NAMD_*.zip /Bossman/cz1/tempering/

clean:
	rm -f NAMD_*.zip *~ .*~ */*~ */.*~
