#!/bin/bash
#
# cwd.sh
#
# Call a command in a spcified working directory
#

if [ $# -lt 2 ]; then
	echo "$0: wrong numbe of arguments"
	return 255
fi

DIR=$1
shift
( cd ${DIR} && $* )

