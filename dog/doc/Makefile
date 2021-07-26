# Make a PDF of the dog user guide and optionally Doxygen html.
#
# $Id$
###############################################################################

# Make a PDF user guide document (userguide.pdf).
# It is assumed that you have a working latex distribution.

default:
	pdflatex guidedog
	bibtex   guidedog
	pdflatex guidedog
	pdflatex guidedog

# After running the command below, open ./html/index.html in a web browser.
# It is assumed that you have doxygen and graphviz packages installed.

doxygen:
	doxygen Doxyfile
