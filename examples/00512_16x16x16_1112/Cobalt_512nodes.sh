#!/bin/bash -x
export BG_SHAREDMEMSIZE=32
export L1P_POLICY=std

export PROG=hacc_tpm
export NODES=512
export RANKS_PER_NODE=8
export OMP_NUM_THREADS=8 
export BG_THREADLAYOUT=1   # 1 - default next core first; 2 - my core first

export NPROCS=$((NODES*RANKS_PER_NODE)) 
export OUTPUT=HACC_N${NODES}_R${RANKS_PER_NODE}_${NPROCS}_${OMP_NUM_THREADS}

export VARS="BG_MAPPING=ABCDET:FAST_WAKEUP=TRUE:VPROF_PROFILE=yes:PAMID_COLLECTIVES=1:BG_SHAREDMEMSIZE=${BG_SHAREDMEMSIZE}:OMP_NUM_THREADS=${OMP_NUM_THREADS}:L1P_POLICY=${L1P_POLICY}:BG_THREADLAYOUT=${BG_THREADLAYOUT}:XLSMPOPTS=stack=4000000"

rm -f core.* ${OUTPUT}.cobaltlog ${OUTPUT}.error ${OUTPUT}.output

qsub -q default -A Performance -n ${NODES} --mode c${RANKS_PER_NODE} -t 1:59:00 -O $OUTPUT \
 --env ${VARS} \
 $PROG indat cmbM000.tf m000 INIT ALL_TO_ALL -w -R -N 512 -a final -f refresh -t 16x16x16

