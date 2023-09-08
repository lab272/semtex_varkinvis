# Make a PDF of the semtex user guide and optionally Doxygen html.
###############################################################################

# Make a PDF user guide document (userguide.pdf).
# It is assumed that you have a working latex distribution.

default:
	pdflatex userguide
	bibtex userguide
	pdflatex userguide
	pdflatex userguide

# After running the command below, open ./html/index.html in a web browser.
# It is assumed that you have doxygen and graphviz packages installed.

doxygen:
	doxygen Doxyfile
