#!/bin/sh
#script de translation des fichier fortran DEC -> linux
printf "."
echo $1 >> diarytrans.txt
/bin/cp $1 $1.20020719
/bin/cat $1 | /bin/sed  -e 's/delete(/zuicloseone(/' > $1.tmp
mv  $1.tmp $1
exit
