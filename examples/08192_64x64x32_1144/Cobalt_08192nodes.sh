#!/bin/bash -x

# must be 1144
# --attrs "location=<blockname>"

export QUEUE=default
export NODES=8192
export MAPPING=ACBDET
export GEOMETRY=64x64x32

export BG_SHAREDMEMSIZE=32
export L1P_POLICY=std
export PROG=hacc_tpm
export RANKS_PER_NODE=16
export OMP_NUM_THREADS=4 
export BG_THREADLAYOUT=1   # 1 - default next core first; 2 - my core first

export NPROCS=$((NODES*RANKS_PER_NODE)) 
export OUTPUT=N${NODES}_R${RANKS_PER_NODE}_${NPROCS}_${OMP_NUM_THREADS}

export VARS="VPROF_PROFILE=yes:RUNJOB_MAPPING=${MAPPING}:FAST_WAKEUP=TRUE:PAMID_COLLECTIVES=1:BG_SHAREDMEMSIZE=${BG_SHAREDMEMSIZE}:OMP_NUM_THREADS=${OMP_NUM_THREADS}:L1P_POLICY=${L1P_POLICY}:BG_THREADLAYOUT=${BG_THREADLAYOUT}:XLSMPOPTS=stack=4000000"

rm -f core.* ${OUTPUT}.cobaltlog ${OUTPUT}.error ${OUTPUT}.output

qsub -q ${QUEUE} -A Performance -n ${NODES} --mode c${RANKS_PER_NODE} -t 4:00:00 -O $OUTPUT \
 --env ${VARS} \
 $PROG indat cmbM000.tf m000 INIT ALL_TO_ALL -w -R -N 512 -t ${GEOMETRY}

#touch ${OUTPUT}.output
#tail -f ---disable-inotify ${OUTPUT}.output
