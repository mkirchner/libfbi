#!/bin/bash
#
# sync a new build of the docs with the gh-pages branch
#
# Copyright (c) 2010 Marc Kirchner
#

if [ $# -lt 1 ]; then
    echo "Usage: $0 <libfbi-source-dir> <libfbi-build-dir>"
    exit 1
fi

# the sync target depends on the doc target, hence the doc path should 
# be availale
if [ ! -f $2/doc/html/index.html ]; then
    echo "Please build the doc target before attempting to sync."
    exit 1
fi

# attempt to switch branches
( cd $1 ; git checkout gh-pages )

# make sure we are on gh-pages
if [ `( cd $1; git branch ) | grep '*' | cut -f2 -d' ' | grep '^gh-pages$' | wc -l` -lt 1 ] ; then
    echo "Switching branches failed. Aborting."
    exit 1
fi

# get a copy of the docs
( cd $1 ; rm -rf doc/html/* )
( cd $2 ; tar cf - doc/html | ( cd $1 ; tar xvf - ) )

echo "Synced working copy with current doc."
echo "Please check that everything is ok, then call:"
echo "    git commit -a -m \"doc update\""
echo "    git push gh-pages origin/gh-pages"

# commit the change
#( cd $1 ; git commit -a -m "doc update" )

# switch back to master branch
#( cd $1 ; git checkout master )



