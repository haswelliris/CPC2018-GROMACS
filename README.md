# CPC2018-GROMACS
cpc2018 final application : gromacs
## 正确性校验
算例1 ：温度，总能量3%的相对误差，压强不做要求。  
算例2.3 温度，总能量3%的相对误差，压强10%的相对误差。  
注意check命令只会输出误差大于1%的数值项
### 算例1
```
~/online1/gmx_mpi_d check -e ~/online1/s1.edr -e2 你的输出.edr
```
### 算例2
```
~/online1/gmx_mpi_d check -e ~/online1/s2.edr -e2 你的输出.edr
```
### 算例3
```
~/online1/gmx_mpi_d check -e ~/online1/s3.edr -e2 你的输出.edr
```
# 神威
## 编译
1. 预准备
```
mkdir build
cd build
```
2. cmake选项配置
注意，把`GMX_CYCLE_SUBCOUNTERS`打开以启用模块计时功能  
默认开启`GMX_DOUBLE`双精度浮点  
默认只编译`mdrun`主程序  
```
LD=mpiCC CC=mpicc CXX=mpiCC cmake .. \
-DGMX_FFT_LIBRARY=fftpack \
-DGMX_MPI=on \
-DGMX_DOUBLE=ON -DGMX_BUILD_MDRUN_ONLY=ON \
-DBUILD_SHARED_LIBS=off \
-DGMX_CYCLE_SUBCOUNTERS=on \
-LH
```
3. 编译
```
make -j
```
## 运行
### 手动设置运行步数
```
export GROMACS_STEP=步数
```
  
请把输入文件`ion_channel-st.tpr`和`lignocellulose-rf.BGQ-st.tpr`放在`./bin/`路径下
### 算例1：
```
bsub -I -b -q q_sw_cpc -cgsp 64 -n 16 -np 4  -share_size 6500 -host_stack 500 -J woca ./bin/mdrun_mpi_d -s ./bin/ion_channel-st.tpr -v
```
### 算例2：
```
bsub -I -b -q q_sw_cpc -cgsp 64 -n 64 -np 4  -share_size 6500 -host_stack 500 -J woca2 ./mdrun_mpi_d -s ./lignocellulose-rf.BGQ-st.tpr -v
```
### 算例3：
```
bsub -I -b -q q_sw_cpc -cgsp 64 -n 512 -np 4  -share_size 6500 -host_stack 500 -J woca3 bin/mdrun_mpi_d -s ./bin/lignocellulose-rf.BGQ-st.tpr -v
```
# intel
# 请切换到intel分支
## 编译环境配置
基本编译环境
```
sudo apt-get install build-essential asciidoc binutils bzip2 gawk gettext git libncurses5-dev libz-dev patch unzip zlib1g-dev gcc-multilib p7zip p7zip-full msmtp libssl-dev ccache make cmake ccmake gcc g++
```
自行安装intel编译器parallel studio xe 2018 update3  
以下是icc需要用到的额外包：  
```
sudo apt install  libpmi0 libibumad-dev libibverbs-dev libudev-dev
wget http://mirrors.kernel.org/ubuntu/pool/main/u/udev/libudev0_175-0ubuntu9_amd64.deb
opkg install libudev0_175-0ubuntu9_amd64.deb
```
## 编译
用原版或者intel分支的代码。

```
LD=mpiicc CC=mpiicc CXX=mpiicpc cmake .. -DGMX_FFT_LIBRARY=fftpack -DGMX_MPI=on -DGMX_DOUBLE=ON -DGMX_BUILD_MDRUN_ONLY=ON -DBUILD_SHARED_LIBS=off -LH
```
剩下的都一样，不写了