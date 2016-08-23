clean:
	rm -f ene.log dih.log ke.dat *~ .*~ *.old *.dcd *restart* *.xst *.BAK

fnrst=narrow1.rst

ehist:
	../mkhist.py --dE=5 --dT=0.5 --rst=$(fnrst) ene*.log

dhist:
	../mkhist.py --dx=0.01745329252 --dT=0.5 --rst=$(fnrst) --colT=3 --colE=4 dih.log
	../mkhist.py --dx=0.01745329252 --dT=0.5 --rst=$(fnrst)                   dih0.log

hist: ehist dhist

