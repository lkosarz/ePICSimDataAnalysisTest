#!/bin/bash

EICSHELL=eic-shell

## Pass commands to eic-shell
ENV_VARS=$(cat <<-END

## Set environment
source /opt/detector/epic-main/bin/thisepic.sh
#source /opt/detector/setup.sh

#export LOCAL_PREFIX='pwd'

#source ${LOCAL_PREFIX}/install/bin/thisepic.sh
#source ${LOCAL_PREFIX}/epic/install/setup.sh

### Export detector libraries
#export LD_LIBRARY_PATH=${LOCAL_PREFIX}/epic/install/lib:$LD_LIBRARY_PATH
\n

END

)

CONDOR_DIR=condorSim
OUT_DIR=output

mkdir ${CONDOR_DIR}
mkdir "${CONDOR_DIR}/${OUT_DIR}"

LISTNAME="files_${5}_${4}.list"
NHEAD=$(($((${4}+1))*${3}))

echo "head -${NHEAD} ${1} | tail -${3} | tee ${LISTNAME}"
head -${NHEAD} ${1} | tail -${3} | tee ${LISTNAME}

echo ${LISTNAME}
cat ${LISTNAME}
	
echo -e "${ENV_VARS} root -l -b -q 'readFrameRoot.C+(\"${LISTNAME}\", \"${2}\")'"
echo -e "${ENV_VARS} root -l -b -q 'readFrameRoot.C+(\"${LISTNAME}\", \"${2}\")'"	| ${EICSHELL} 	

exit
