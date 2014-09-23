all: ms.pdf

%.pdf : %.tex
	pdflatex $*
	pdflatex $*
	pdflatex $*
clean:
	rm -rf $(foreach suff, pdf aux bbl blg dvi log, *.$(suff))
