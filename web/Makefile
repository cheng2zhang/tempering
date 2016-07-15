clean:
	rm -rf a.out *~ .*.un~ */*~ */.*.un~ *.dat *.pos *.log
	rstrip.py -Rlv

Bossman:
	python htmlpack.py
	rsync -avz stlj_pack.html /Bossman/cz1/tempering/web/

Dropbox: clean
	$(MAKE) -C .. $@
