option = # -shell-escape

compiler = pdflatex

ics_and_preds.pdf: ics_and_preds.tex
	${compiler} ${option}  ics_and_preds.tex
	bibtex ics_and_preds
	${compiler} ${option}  ics_and_preds.tex
	${compiler} ${option}  ics_and_preds.tex
	gnome-open ics_and_preds.pdf

clean:
	rm -f *~



