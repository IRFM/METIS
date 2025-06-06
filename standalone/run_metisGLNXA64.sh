#!/bin/sh
# script for execution of deployed applications
#
# Sets up the MCR environment for the current $ARCH and executes 
# the specified command.
#
exe_name=$0
exe_dir=`dirname "$0"`
echo "------------------------------------------"
if [ "x$1" = "x" ]; then
  echo Usage:
  echo    $0 \<deployedMCRroot\> args
else
  echo Setting up environment variables
  MCRROOT="$1"
  echo ---
  LD_LIBRARY_PATH_MEM=${LD_LIBRARY_PATH};
  LD_LIBRARY_PATH=.:${MCRROOT}/runtime/glnxa64 ;
  LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/bin/glnxa64 ;
  LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/sys/os/glnxa64;
	MCRJRE=${MCRROOT}/sys/java/jre/glnxa64/jre/lib/amd64 ;
	LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}/native_threads ; 
	LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}/server ;
	LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}/client ;
	LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE} ;  
	LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${LD_LIBRARY_PATH_MEM};  
  XAPPLRESDIR=${MCRROOT}/X11/app-defaults ;
  export LD_LIBRARY_PATH;
  export XAPPLRESDIR;
  echo LD_LIBRARY_PATH is ${LD_LIBRARY_PATH};
  shift 1
  args=
  while [ $# -gt 0 ]; do
      token=$1
      args="${args} ${token}" 
      shift
  done
  "${exe_dir}"/metisGLNXA64 $args
fi
exit

# On the target computer, append the following to your LD_LIBRARY_PATH environment variable:  /pfs/home/jfa/public/metisruntime/MCR/v80/runtime/glnxa64:/pfs/home/jfa/public/metisruntime/MCR/v80/bin/glnxa64:/pfs/home/jfa/public/metisruntime/MCR/v80/sys/os/glnxa64:/pfs/home/jfa/public/metisruntime/MCR/v80/sys/java/jre/glnxa64/jre/lib/amd64/native_threads:/pfs/home/jfa/public/metisruntime/MCR/v80/sys/java/jre/glnxa64/jre/lib/amd64/server:/pfs/home/jfa/public/metisruntime/MCR/v80/sys/java/jre/glnxa64/jre/lib/amd64    Next, set the XAPPLRESDIR environment variable to the following value:  /pfs/home/jfa/public/metisruntime/MCR/v80/X11/app-defaults  
