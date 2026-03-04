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
# specific env
module load mamba/python39 mamba-forge
module load tools_dc/16_m2020b

# launch Muslce3 workflow
DIR_ORIGIN="$( dirname -- "${BASH_SOURCE[0]}"; )";   # Get the directory name
DIR_ORIGIN="$( realpath -ePL -- "$DIR_ORIGIN"; )";    # Resolve its full path if need be
echo "Simulation source files are in $DIR_ORIGIN"
# path to local Muscle3 installation if needed
export PATH=$PATH:"$HOME/.local/bin"
# temporary directory for Matlab IDS serialize and deserialize  (one per Matlab actor)
# useful only with AL before 5.4.2
mkdir  -p /dev/shm/$USER/$$/MUSCLE3
mkdir  -p /dev/shm/$USER/$$/MUSCLE3/metis
mkdir  -p /dev/shm/$USER/$$/MUSCLE3/scheduler
export IMAS_AL_SERIALIZER_TMP_DIR_METIS="/dev/shm/$USER/$$/MUSCLE3/metis"
export IMAS_AL_SERIALIZER_TMP_DIR_SCHEDULER="/dev/shm/$USER/$$/MUSCLE3/scheduler"
# path to metis4muscle3.m file
export DIR_METIS4MUSCLE3="$DIR_ORIGIN/mfile"
echo "path to metis4muscle3 is $DIR_METIS4MUSCLE3"
# launch workflow
muscle_manager --start-all $DIR_ORIGIN/YMMSL/metis4muscle3_dataprovider.ymmsl 
# clean temporary directories
rm -rf /dev/shm/$USER/$$
