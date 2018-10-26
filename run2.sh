cd /home/export/online1/cpc152/finals/ljx/code/result2
LOAD_BALANCE_STEP=37 bsub -I -b -q q_sw_cpc -cgsp 64 -n 64 -np 4  -share_size 6500 -host_stack 500 -J woca2 ../build/bin/mdrun_mpi -s ../sample/lignocellulose-rf.BGQ-st.tpr -v
