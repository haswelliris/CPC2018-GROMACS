cd /home/export/online1/cpc152/finals/ljx/code/result3
LOAD_BALANCE_STEP=8 bsub -I -b -q q_sw_cpc -cgsp 64 -n 512 -np 4  -share_size 6500 -host_stack 500 -J woca2 ../build/bin/mdrun_mpi -s ../sample/lignocellulose-rf.BGQ-st.tpr -v
