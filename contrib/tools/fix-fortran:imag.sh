#!/bin/sh

fortran=`find src -name '*.f'`

for f in ${fortran}; do
    sed -E \
	-e 's/([^AaDd])[Ii][Mm][Aa][Gg][ 	]*\(/\1dimag(/g' \
	${f} > ${f}.new
    if cmp ${f} ${f}.new >/dev/null; then
	rm -f ${f}.new
    else
	mv -f ${f}.new ${f}
    fi
done

# End of File
