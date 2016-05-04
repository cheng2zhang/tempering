cmp.png: cmp.gp vscale.txt novscale.txt
	./$<
	gnome-open $@

# The following command is used to get the number of bins used in adaptive tempering
# grep adaptTempBins argon_st.conf | awk '{print $$2}'

# use the original NAMD (without velocity scaling after tempering)
novscale.txt:
	make clean
	../../../NAMD_CVS-2016-05-02_Source/Linux-x86_64-g++/namd2 +p2 argon_st.conf
	tail -n`grep adaptTempBins argon_st.conf | awk '{print $$2}'` argon_st1.rst > $@

# use the modified NAMD (with velocity scaling after tempering)
vscale.txt:
	make clean
	../../../NAMD_modified/Linux-x86_64-g++/namd2 +p2 argon_st.conf
	tail -n`grep adaptTempBins argon_st.conf | awk '{print $$2}'` argon_st1.rst > $@

clean:
	rm -f argon_st[1-9]* FFTW_NAMD_* *~
