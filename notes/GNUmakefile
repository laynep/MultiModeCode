option = # -shell-escape

compiler = pdflatex

multi-modecodev3.pdf: multi-modecodev3.tex
	${compiler} ${option}  multi-modecodev3.tex
	#bibtex multi-modecodev3
	${compiler} ${option}  multi-modecodev3.tex
	${compiler} ${option}  multi-modecodev3.tex
	gnome-open multi-modecodev3.pdf

clean:
	rm -f *~



