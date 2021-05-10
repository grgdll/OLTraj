#!/usr/bin/bash
# this script should be run from inside the Data directory

# Source the python env
source /users/rsg/fne/PyVirtualEnv/TAPAS_Python3/bin/activate
echo "Start `date`"
for iyear in {1998..2018}
do
   python download_script.py --year $iyear 
done
echo "End `date`"
