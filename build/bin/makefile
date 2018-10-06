run :
	(cd /home/export/online1/cpc152/finals/fgn/CPC2018-GROMACS/src/gromacs/mdlib/nbnxn_kernels/sw_subcore && make device) || exit 1
	(cd .. && make -j 20) || exit 1
	GROMACS_STEP=1 bsub -I -b -q q_sw_cpc -cgsp 64 -n 64 -np 4  -share_size 6500 -host_stack 500 -J woca ./mdrun_mpi_d -s /home/export/online1/cpc/ion_channel-st.tpr -v
