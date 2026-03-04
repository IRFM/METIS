#/bin/bash 
#
# script to test metis4muscle3 using local installation of muslce3:
#    module load tools_dc/14_m2020b
#    pip install --user muscle3
#    export PATH=$PATH:"/Home/JA32999/.local/bin"
#    pip3 install --user ymmsl
# 
# alternative installation at :
#
# 	python3 -m venv muscle3_venv
# 	muscle3_venv/bin/activate
# 	(muscle3_venv)~$ pip install -U pip setuptools wheel
# 	(muscle3_venv)~$ pip install muscle3
# 	(muscle3_venv)~$  pip3 install ymmsl
# 
#

module load mamba/python39 mamba-forge
module load tools_dc/15_m2020b
export PATH=$PATH:"$HOME/.local/bin"
mkdir  -p /dev/shm/$USER/MUSCLE3
mkdir  -p /dev/shm/$USER/MUSCLE3/metis
mkdir  -p /dev/shm/$USER/MUSCLE3/scheduler
muscle_manager --start-all ../YMMSL/metis4muscle3_test.ymmsl 