#include "sw/SwDevice.h"

void device_run()
{
    static int called = 0;
    if(called == 0)
    {
        TLOG("Device noticed.\n");
        called = 1;
    }
}