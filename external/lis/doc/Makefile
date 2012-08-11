# Makefile for LaTeX Systems

FILE=lis-manual

all: en ja

# English version

en: 
	latex $(FILE)-en.tex
	latex $(FILE)-en.tex
	dvipdfmx $(FILE)-en.dvi

# Japanese version

ja: 
	platex $(FILE)-ja.tex
	platex $(FILE)-ja.tex
	dvipdfmx $(FILE)-ja.dvi

clean: 
	-$(RM) *.dvi *.aux *.toc *.log

distclean: clean
	-$(RM) *.pdf
