#!/bin/bash

# simple script to plot the output wave functions of denchar and intw with gnuplot

if [ -z "$1" ]; then
  echo "No argument supplied"
  exit
fi

band=$1

SPIN_COLIN=$(cat h2o.save.intw/crystal.dat | grep -A1 "SPIN_COLIN" | tail -n 1 | xargs)
SPIN_NONCOLIN=$(cat h2o.save.intw/crystal.dat | grep -A1 "SPIN_NONCOLIN" | tail -n 1 | xargs)
SPIN_SO=$(cat h2o.save.intw/crystal.dat | grep -A1 "SPIN_SO" | tail -n 1 | xargs)
if [ "${SPIN_COLIN}" == "T" ]; then
  echo "COLINEAR"
  spin=2
elif [ "${SPIN_NONCOLIN}" == "T" ]; then
  echo "NON COLINEAR"
  spin=4
elif [ "${SPIN_SO}" == "T" ]; then
  echo "NON COLINEAR"
  spin=4
else
  echo "NON POLARIZED"
  spin=1
fi


if test h2o.selected.WFSX -nt denchar.out; then
  # h2o.selected.WFSX is newer than denchar.out
  echo "Running denchar..."
  rm -f h2o.CON.*
  cp h2o.selected.WFSX h2o.WFSX
  ../../../../Denchar/Src/denchar < denchar.fdf > denchar.out
fi

if test h2o.selected.WFSX -nt plot_wfc.out; then
  # h2o.selected.WFSX is newer than plot_wfc.out
  echo "Running plot_wfc.x..."
  rm -f fort.*
  /home/haritz/Codes/intw/intw3.3/build/src/utilities/plot_wfc.x < plot_wfc.in > plot_wfc.out
fi


#
# NON-POLARIZED
#
if [ "${spin}" -eq "1" ]; then
  intw_band=$(echo "10000+${band}" | bc)
  # cat h2o.CON.K1.WF.${band}.MOD | sed ':a;N;$!ba;s/\n/newline/g' | sed 's/newline\s*\S*\s*\S*\s*\S*newline newline/newline/g' | sed 's/newline/\n/g' > tmp_denchar
  cat h2o.CON.K1.WF.${band}.MOD | sed -r '/^\s*$/d' > tmp_denchar
  gnuplot -p <<-EOF
    #
    set term x11 1
    set title "h2o.CON.K1.WF.${band}.MOD vs fort.${intw_band}"
    # plot "h2o.CON.K1.WF.${band}.MOD" every 121::::14398:14398 u 1:3 , "fort.${intw_band}" every ::::119:119 u 1:3
    plot "tmp_denchar" every 120::::14279:14279 u 1:3 w l notitle, "fort.${intw_band}" every ::::119:119 u 1:3 notitle
    # set term x11 2
    # plot "tmp_denchar" every 1::::119:119 u 2:3 w l, "fort.${intw_band}" every 120::::14400:14400 u 2:3
EOF
  rm -f tmp_denchar
fi

#
# COLINEAR
#
if [ "${spin}" -eq "2" ]; then
  intw_band_up=$(echo "10000+2*${band}-1" | bc)
  intw_band_down=$(echo "20000+2*${band}" | bc)
  cat h2o.CON.K1.WF.${band}.UP.MOD | sed -r '/^\s*$/d' > tmp_denchar_up
  cat h2o.CON.K1.WF.${band}.DOWN.MOD | sed -r '/^\s*$/d' > tmp_denchar_down
  gnuplot -p <<-EOF
    #
    set term x11 1
    set title "h2o.CON.K1.WF.${band}.UP.MOD vs fort.${intw_band_up}"
    plot "tmp_denchar_up" every 120::::14279:14279 u 1:3 w l notitle, "fort.${intw_band_up}" every ::::119:119 u 1:3 notitle
    set term x11 2
    set title "h2o.CON.K1.WF.${band}.DOWN.MOD vs fort.${intw_band_down}"
    plot "tmp_denchar_down" every 120::::14279:14279 u 1:3 w l notitle, "fort.${intw_band_down}" every ::::119:119 u 1:3 notitle
EOF
  rm -f tmp_denchar_up tmp_denchar_down
fi


#
# NON-COLINEAR
#
if [ "${spin}" -eq "4" ]; then
  intw_band_up=$(echo "10000+${band}" | bc)
  intw_band_down=$(echo "20000+${band}" | bc)
  cat h2o.CON.K1.WF.${band}.UP.MOD | sed -r '/^\s*$/d' > tmp_denchar_up
  cat h2o.CON.K1.WF.${band}.DOWN.MOD | sed -r '/^\s*$/d' > tmp_denchar_down
  gnuplot -p <<-EOF
    #
    set term x11 1
    set title "h2o.CON.K1.WF.${band}.UP.MOD vs fort.${intw_band_up}"
    plot "tmp_denchar_up" every 120::::14279:14279 u 1:3 w l notitle, "fort.${intw_band_up}" every ::::119:119 u 1:3 notitle
    set term x11 2
    set title "h2o.CON.K1.WF.${band}.DOWN.MOD vs fort.${intw_band_down}"
    plot "tmp_denchar_down" every 120::::14279:14279 u 1:3 w l notitle, "fort.${intw_band_down}" every ::::119:119 u 1:3 notitle
EOF
fi
