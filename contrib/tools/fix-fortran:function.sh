#!/bin/sh

fortran=`find src -name '*.f'`

for f in ${fortran}; do
    sed -E \
	-e 's/^(.*[^ ])[ 	]+[Ff][Uu][Nn][Cc][Tt][Ii][Oo][Nn][ 	]+([^ ]+)\*([0-9]+)[ 	]*\(/\1*\3 function \2(/' \
	${f} > ${f}.new
    if cmp ${f} ${f}.new >/dev/null; then
	rm -f ${f}.new
    else
	mv -f ${f}.new ${f}
    fi
done

# End of File
