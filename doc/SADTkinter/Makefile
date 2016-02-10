
LATEX=platex
MAKEINDEX=makeindex
DVIPS=dvips
DVIPDF=dvipdfm

DOCPREFIX=SADTkinter

ps:	${DOCPREFIX}.ps

pdf:	${DOCPREFIX}.pdf

clean:
	rm -f ${DOCPREFIX}.aux ${DOCPREFIX}.log ${DOCPREFIX}.toc \
	${DOCPREFIX}.idx ${DOCPREFIX}.ilg ${DOCPREFIX}.ind \
	${DOCPREFIX}.dvi ${DOCPREFIX}.ps ${DOCPREFIX}.pdf

release:	ps pdf

${DOCPREFIX}.dvi:	${DOCPREFIX}.tex
	rm -f ${DOCPREFIX}.ind
	${LATEX} ${DOCPREFIX}.tex
	${LATEX} ${DOCPREFIX}.tex
	${MAKEINDEX} ${DOCPREFIX}.idx
	${LATEX} ${DOCPREFIX}.tex
	${LATEX} ${DOCPREFIX}.tex

${DOCPREFIX}.ps:	${DOCPREFIX}.dvi
	${DVIPS}  -o $@ ${DOCPREFIX}.dvi

${DOCPREFIX}.pdf:	${DOCPREFIX}.dvi
	${DVIPDF} -o $@ ${DOCPREFIX}.dvi

# End of File
