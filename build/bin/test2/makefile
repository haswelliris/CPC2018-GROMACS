run :
	bsub -I -b -q q_sw_cpc -cgsp 64 -n 256 -np 4  -share_size 6500 -host_stack 500 -J woca2 ../mdrun_mpi_d -s /home/export/online1/cpc/lignocellulose-rf.BGQ-st.tpr -v
	~/online1/gmx_mpi_d check -e ~/online1/s2.edr -e2 ener.edr
