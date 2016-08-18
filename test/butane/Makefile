clean:
	rm -f ene*.log dih*.log ke.dat *~ .*~ *.old *.dcd *restart* *.xst *.BAK

ehist:
	../mkhist.py --dE=5 --dT=0.5 ene*.log

dhist:
	../mkhist.py --dE=0.0698 --dT=0.5 --colT=3 dih.log
	../mkhist.py --dE=0.0698 --dT=0.5          dih0.log

hist: ehist dhist

