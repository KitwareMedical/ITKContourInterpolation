

NAME = MorphologicalContourInterpolation

BIBTEXDATABSENAME = InsightJournal


all: pdf

pdf: ${NAME}.pdf

bbl: ${NAME}.bbl


${NAME}.bbl: ${BIBTEXDATABSENAME}.bib ${NAME}.aux
	bibtex ${NAME}


${NAME}.pdf: ${NAME}.tex  ${BIBTEXDATABSENAME}.bib
	pdflatex ${NAME}.tex
	bibtex ${NAME}
	pdflatex ${NAME}.tex
	pdflatex ${NAME}.tex


${NAME}.aux: ${NAME}.tex  ${BIBTEXDATABSENAME}.bib
	pdflatex ${NAME}.tex


clean:
	rm ${NAME}.pdf ${NAME}.aux ${NAME}.log ${NAME}.bbl ${NAME}.blg ${NAME}.toc ${NAME}.out ${NAME}.brf

