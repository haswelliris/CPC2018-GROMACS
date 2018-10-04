#HOW TO USE
主核代码include SwHost.h
从核代码include SwDevice.h
调整配置可以写在 SwConfig.h

#ABOUT
SwHost.h和SwDevice.h都包含了SwConfig.h，作为共同配置文件。
从核代码的运行入口是device_main()
需要用户自定义的函数是device_run()，在框架中以extern的形式定义了，所以在写从核代码的时候用户入口是device_run()

比如在swlbm中我就让device_run()分别调用1*1*50和2*2*20的计算
