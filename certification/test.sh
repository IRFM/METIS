#!/bin/sh

#source arch.inc 2> /dev/null

#set -x 

trap '
exit ' 2

if [ ! -d certification/output ]; then
   mkdir certification/output
fi

if [ "$1" = "modules" ]; then
  if [ $# = 1 ]; then  
    for i in `ls certification/modules/Cronos_test*.mat*`
    do
      test=`echo $i | awk -F "." '{ print $1 }'`
      lance=` echo matlab -nodesktop -nodisplay -r \"zineb_test\(\'${test}\'\)\;exit\" `
      sh -c "${lance}"
    done
    for i in `ls -d certification/modules/*/ `
    do
      if [ "$i" != "certification/modules/CVS/" ]; then
         time -p ${i}/check.sh 
      fi
    done
  else
    if [ -d certification/${1}/${2} ]; then
      time -p certification/${1}/${2}/check.sh
    else
      test=`echo $2 | awk -F "." '{ print $1 }'` 
      lance=` echo matlab -nodesktop -nodisplay -r \"zineb_test\(\'${test}\'\)\;exit\" `
      sh -c "${lance}"
    fi
  fi
fi
if [ "$1" = "cronosfunctions" ]; then
   if [ $# = 1 ]; then
      for i in `ls certification/cronosfunctions/Cronos_test*.mat*`
      do
        test=`echo $i | awk -F "." '{ print $1 }'`
        lance=` echo matlab -nodesktop -nodisplay -r \"zineb_test\(\'${test}\'\)\;exit\" `
        sh -c "${lance}"
      done
   else
      test=`echo $2 | awk -F "." '{ print $1 }'`
      lance=` echo matlab -nodesktop -nodisplay -r \"zineb_test\(\'${test}\'\)\;exit\" `
      sh -c "${lance}"
   fi
fi
if [ "$1" = "cronostestrun" ]; then
   if [ $# = 1 ]; then
      for i in `ls certification/cronostestrun/*.mat*`
      do
        test=`echo $i | awk -F "." '{ print $1 }'`
        lance=` echo matlab -nodesktop -nodisplay -r \"zineb_test\(\'${test}\'\)\;exit\" `
        sh -c "${lance}"
      done
   else
      test=`echo $2 | awk -F "." '{ print $1 }'`
      lance=` echo matlab -nodesktop -nodisplay -r \"zineb_test\(\'${test}\'\)\;exit\" `
      sh -c "${lance}"
   fi
fi
if [ "$1" = "metis" ]; then
   if [ $# = 1 ]; then
      for i in `ls certification/metis/*.mat*`
      do
        test=`echo $i`
        lance=` echo matlab -nodesktop -nodisplay -r \"metis_test\(\'${test}\'\)\;exit\" `
        sh -c "${lance}"
      done
   else
      test=`echo $2`
      lance=` echo matlab -nodesktop -nodisplay -r \"metis_test\(\'${test}\'\)\;exit\" `
      sh -c "${lance}"
   fi
fi
if [ "$1" = "imasinmetis" ]; then
   if [ $# = 1 ]; then
        lance=` echo matlab -nodesktop -nodisplay -r \"metis_test_imas\;exit\" `
        echo $lance
        sh -c "${lance}"
    fi
fi

