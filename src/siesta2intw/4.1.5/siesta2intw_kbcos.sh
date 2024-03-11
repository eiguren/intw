#!/bin/sh


if [ "$1" = "" ]; then
  echo "  ERROR: No input file"
  exit
fi

PSEUDO=$1

extension="${PSEUDO#*.}"
if [ ! $extension = "ion.xml" ]; then
  echo "The extension of the file should be ion.xml!"
  exit 1
fi


if [ -f $PSEUDO ]; then
  LABEL=$(grep "<label>" $PSEUDO | sed -e 's/<label>//' -e 's/<\/label>//' -e 's/ //g')
fi


if [ "$2" = "" ]; then
  path=$(dirname "${PSEUDO}")
  files=$(find $path -type f -print | xargs grep "*  WELCOME TO SIESTA  *" | wc -l)
  if [ $files -gt 1 ]; then
    echo "More than one output file!"
    echo "Run the script again with a second argument indicating the output file."
    exit 1
  else
    OUTPUT=$(find $path -type f -print | xargs grep "\*  WELCOME TO SIESTA  \*" | tail -n 1 | sed -e 's/\*  WELCOME TO SIESTA  \*//' -e 's/ //g' -e 's/\://g')
  fi
else
  OUTPUT=$2
fi


sed -n "/atom: Called for $LABEL/,/atom: -------------------------------------------------------------------------/p" $OUTPUT | grep "kbcos"
sed -n "/atom: Called for $LABEL/,/atom: -------------------------------------------------------------------------/p" $OUTPUT | grep "kbcos" | sed -e 's/l=//' -e 's/rc=//' -e 's/el=//' -e 's/Ekb=//' -e 's/kbcos=//' | awk '{print $5}' > $LABEL.kbcos

exit 0
