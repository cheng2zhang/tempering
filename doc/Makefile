prj=vsTmove
ver=3.0

figgps=$(wildcard fig/*.gp)
pdffigs=$(patsubst %.gp,%.pdf, $(figgps))
epsfigs=$(patsubst %.gp,%.eps, $(figgps))

alltexpdf = $(patsubst %.tex,%.pdf, $(wildcard *.tex))

all: $(alltexpdf)

$(alltexpdf) : %.pdf : %.tex simul.bib $(pdffigs) Makefile
	make -C fig
	pdflatex $<
	bibtex $*
	pdflatex $<
	pdflatex $<

zip: vsTmove.zip

vsTmove.zip: vsTmove.tex vsTmove.pdf \
             simul.bib Makefile README.md \
	     $(epsfigs) $(pdffigs) $(figgps)
	zip -r $@ $^
