#!/bin/sh 
# usage /path/zcapture.sh "name" "path_to_bin_zwpick" "filename.formatname"
idw=`xwininfo -name "$1" | grep  "Window id:" | awk -F\: '{print $3}'|  awk '{print $1}'` 
echo $idw
exec `$2/xwpick -window $idw  $3
exit
