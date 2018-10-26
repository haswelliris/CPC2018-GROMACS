cd /home/export/online1/cpc152/finals/ljx/code/result1
LOAD_BALANCE_STEP=10 bsub -I -b -q q_sw_expr -cgsp 64 -n 16 -np 4  -share_size 6500 -host_stack 500 -J woca ../build/bin/mdrun_mpi -s ../sample/ion_channel-st.tpr -v
