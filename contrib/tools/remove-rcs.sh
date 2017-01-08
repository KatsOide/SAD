#!/bin/sh

usage () {
    echo "usage: ${cmd_name} [-w] [-l] -t RCS_tag_name files..."
    exit 2
}

cmd_name="`basename $0`"

tag_name=
remove_mode=WORD
while [ $# -gt 0 ]; do
    case $1 in
	-t)
	    if [ $# -gt 1 ]; then
		shift
		tag_name="$1"
	    else
		usage
	    fi
	    ;;
	-w)
	    remove_mode=WORD
	    ;;
	-l)
	    remove_mode=LINE
	    ;;
	-h)
	    usage
	    ;;
	--)
	    break
	    ;;
	-*)
	    usage
	    ;;
	*)
	    break
	    ;;
    esac
    shift
done

if [ "x${tag_name}" = "x" ]; then
    usage
fi

while [ $# -gt 0 ]; do
    if [ -r "$1" ]; then
	case ${remove_mode} in
	    LINE)
		grep -v -e "\\\$${tag_name}:[^\$][^\$]*\\\$" "$1" >"$1.new"
		if cmp "$1" "$1.new" >/dev/null; then
		    rm -f "$1.new"
		else
		    mv -f "$1.new" "$1"
		fi
		;;
	    WORD)
		sed -e "s/\\\$${tag_name}:[^\$][^\$]*\\\$//g" "$1" >"$1.new"
		if cmp "$1" "$1.new" >/dev/null; then
		    rm -f "$1.new"
		else
		    mv -f "$1.new" "$1"
		fi
		;;
	    *)
		;;
	esac
    fi
    shift
done

exit 0

# End of File
