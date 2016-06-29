cmp.png: cmp.gp vscale.txt novscale.txt
	./$<
	gnome-open $@

# The following command is used to get the number of bins used in adaptive tempering
# grep adaptTempBins argon_st.conf | awk '{print $2}'
# The grep part gets a line like, "adaptTempBins   60",
# then the awk part takes the second column of 60.
# In makefile $2 should appear as $$2 for the second column.

# use the original NAMD (without velocity scaling after tempering)
novscale.txt:
	make clean
	../../../NAMD_2.11_Source/Linux-x86_64-g++/namd2 argon_st.conf
	tail -n`grep adaptTempBins argon_st.conf | awk '{print $$2}'` argon_st1.rst > $@

# use the modified NAMD (with velocity scaling after tempering)
vscale.txt:
	make clean
	../../../NAMD_mods/NAMD_2.11_tempering/Linux-x86_64-g++/namd2 argon_st.conf
	tail -n`grep adaptTempBins argon_st.conf | awk '{print $$2}'` argon_st1.rst > $@

clean:
	rm -f argon_st[1-9]* FFTW_NAMD_* ke.dat *.nhc *~
