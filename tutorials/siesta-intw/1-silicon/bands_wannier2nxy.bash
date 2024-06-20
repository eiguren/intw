#!/usr/bin/env bash

PREFIX=si


cat ${PREFIX}_band.dat | sed -e's/[[:space:]]*$//' | awk 'BEGIN{k=0} {if ($0=="") {k=0 ; print $0} else {k=k+1 ; print k, $2}}'
