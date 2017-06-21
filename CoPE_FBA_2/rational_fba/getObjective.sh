#!/bin/sh
tail -1 $1 > optTail

cat optTail | gawk '{ for(i=1; i<=NF; i++) if( index($i,"1") != 0) print i-1; }'
