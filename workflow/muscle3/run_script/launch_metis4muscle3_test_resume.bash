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

# create resume directory
mkdir -p test_resume

# launch Muslce3 workflow
DIR_ORIGIN="$( dirname -- "${BASH_SOURCE[0]}"; )";   # Get the directory name
DIR_ORIGIN="$( realpath -ePL -- "$DIR_ORIGIN"; )";    # Resolve its full path if need be
echo "Simulation source files are in $DIR_ORIGIN"
# path to local Muscle3 installation if needed
export PATH=$PATH:"$HOME/.local/bin"
# path to metis4muscle3.m file
export DIR_METIS4MUSCLE3="$DIR_ORIGIN/mfile"
echo "path to metis4muscle3 is $DIR_METIS4MUSCLE3"
# launch workflow
timeout --kill-after=3600 --signal=SIGTERM 3600 muscle_manager --log-level DEBUG --run-dir run_METIS4MUSCLE3_FULL_TEST_20250408_102954 --start-all $DIR_ORIGIN/YMMSL/metis4muscle3_full_test.ymmsl run_METIS4MUSCLE3_FULL_TEST_20250408_102954/instances/metis/snapshots/metis_2.pack
# clean temporary directories
rm -rf /dev/shm/$USER/$$
