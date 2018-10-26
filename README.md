# 第二届国产CPU并行应用挑战赛
# 决赛题：GROMACS在神威·太湖之光上的优化加速
## GROMACS简介
[GROMACS](http://www.gromacs.org/)是一款分子动力学应用程序，旨在模拟包含数百到数百万个粒子的系统的牛顿运动方程。  
GROMACS设计用于模拟具有大量复杂键合相互作用的生物化学分子，包括蛋白质，脂，核酸甚至是一些非生物来源的分子，多聚物等等。

## 神威·太湖之光简介
[神威·太湖之光](https://zh.wikipedia.org/wiki/%E7%A5%9E%E5%A8%81%C2%B7%E5%A4%AA%E6%B9%96%E4%B9%8B%E5%85%89) （Sunway TaihuLight）是由中国国家并行计算机工程技术研究中心研制的超级计算机，2016年6月20日在LINPACK性能测试中以 93 PFLOPS 的测试结果超越同为中国组建的天河二号，成为世界上最快的超级计算机，直到2018年6月8日被美国的超级计算机高峰超越。截止2018年10月，神威·太湖之光为世界上第二快的超级计算机，中国第一快的超级计算机。  

## 本次比赛中用到的优化策略
我们从计算、访存、模型三个方面进行了优化：  
1. 修改任务划分策略，使用athread库实现从核加速PP计算核
2. 使用软件实现的缓存预取优化访存开销
3. 使用预估计算量的负载均衡使任务划分更均匀
4. 对PP计算核进行了循环展开和部分向量化
5. PME与PP主从核异步计算彼此掩盖了部分计算时间

更多技术细节详见[优化报告](./CPC2018-SYSU.pdf)

## 总结
|算例|原始运行时间(秒)|优化后运行时间(秒)|加速比|
|----|---|---|---|
|1|6556.575|785|835.2%|
|2|3781.820|472.551|800.3%|
|3|578.761|119.415|484.7%|
最终，我们在[CPC2018决赛](http://cpc.nsccwx.cn/)中取得第4名的成绩。  

## 编译运行说明
通过bash脚本快速编译运行，注意修改路径，并把[测试样例](https://mp.weixin.qq.com/s/-G66jmIXe6BPLjCA1HwZZA)放在`./sample/`文件夹中

```
init_build.sh ：执行完全重新编译
run1.sh       ：运行算例1
run2.sh       ：运行算例2
run3.sh       ：运行算例3
```
   
# 原版GROMACS编译运行和正确性校验说明
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
## 神威
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
-DGMX_DOUBLE=ON \
-DGMX_BUILD_MDRUN_ONLY=ON \
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
## intel平台
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
