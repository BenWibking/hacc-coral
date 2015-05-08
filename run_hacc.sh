ln -s ./src/cpu/linux-generic/hacc_tpm .
mpirun -np 8 ./hacc_tpm indat cmbM000.tf m000 INIT ALL_TO_ALL -w -R -N 512 -t 2x2x2
