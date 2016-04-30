CC = icc -Wall -Wremarks -O2 -g
LM =
#CC = gcc -O3 -Wall -Wextra -g
#LM = -lm

srcs = $(wildcard *.c)
bins = $(patsubst %.c,%, $(srcs))
deps = $(wildcard *.h) Makefile

all: $(bins)

$(bins): % : %.c $(deps)
	$(CC) $< -o $@ $(LM)

clean:
	rm -rf $(bins) $(bin2ds) $(bins_d) $(bins2ds_d) \
	  a.out *~ .*.un~ */*~ */.*.un~ *.dat *.pos *.log
	rstrip.py -Rlv

zip:
	zip -r temper.zip *.[ch] *.gp README.md Makefile

mddemo: tempmd
	./tempmd 1 > good.dat
	./tempmd 0 > bad.dat
	gnuplot vscale.gp
	gnome-open vscale.png
