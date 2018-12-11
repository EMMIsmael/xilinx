/*
 * hw_util.c
 *
 *  Created on: 16 May 2017
 *      Author: JULIAN MORTIMER
 */
#include "include.h"

u32  hwRdReg( u32 addr, u32 offs );
void hwWrReg( u32 addr, u32 offs, u32 data );
u64  hwRdReg64( u32 addr, u32 offs );
void hwWrReg64( u32 addr, u32 offs, u64 data );


void hwWriteCtrlReg( INTPTR baseAddress, u32 regOffset, u32 data );

void hwWriteData( INTPTR baseAddress,
		         u32 ctrlRegOffs, u32 ctrlRegWr,
		         u32 addrRegOffs, u32 addr,
		         u32 dataRegOffs, u32 data );

void hwWriteData64( INTPTR baseAddress,
                 u32 ctrlRegOffs, u32 ctrlRegWr,
                 u32 addrRegOffs, u32 addr,
                 u32 dataRegOffs, u64 data );

u32 hwReadData( INTPTR baseAddress,
		       u32 ctrlRegOffs, u32 ctrlRegRd, u32 ctrlRegAck,
		       u32 statRegOffs, u32 statRegRdy,
		       u32 addrRegOffs, u32 addr,
		       u32 dataRegOffs );

u64 hwReadData64( INTPTR baseAddress,
               u32 ctrlRegOffs, u32 ctrlRegRd, u32 ctrlRegAck,
               u32 statRegOffs, u32 statRegRdy,
               u32 addrRegOffs, u32 addr,
               u32 dataRegOffs );

void hwWriteCtrlReg( INTPTR baseAddress, u32 regOffset, u32 data )
{
//    hwWrReg( baseAddress, regOffset, 0 );
    hwWrReg( baseAddress, regOffset, data );
    hwWrReg( baseAddress, regOffset, 0 );
}


void hwWriteData( INTPTR baseAddress,
		         u32 ctrlRegOffs, u32 ctrlRegWr,
		         u32 addrRegOffs, u32 addr,
		         u32 dataRegOffs, u32 data )
{
    hwWrReg( baseAddress, addrRegOffs, addr );
    hwWrReg( baseAddress, dataRegOffs, data );
    hwWriteCtrlReg( baseAddress, ctrlRegOffs, ctrlRegWr );
}

void hwWriteData64( INTPTR baseAddress,
                 u32 ctrlRegOffs, u32 ctrlRegWr,
                 u32 addrRegOffs, u32 addr,
                 u32 dataRegOffs, u64 data )
{
    hwWrReg( baseAddress, addrRegOffs, addr );
    hwWrReg64( baseAddress, dataRegOffs, data );
    hwWriteCtrlReg( baseAddress, ctrlRegOffs, ctrlRegWr );
}

u32 hwReadData( INTPTR baseAddress,
		       u32 ctrlRegOffs, u32 ctrlRegRd, u32 ctrlRegAck,
		       u32 statRegOffs, u32 statRegRdy,
		       u32 addrRegOffs, u32 addr,
		       u32 dataRegOffs )
{
    u32 u;
	hwWrReg( baseAddress, addrRegOffs, addr );
    hwWriteCtrlReg( baseAddress, ctrlRegOffs, ctrlRegRd );
    double f = SysTimer->ReadTimerLongSec( SysTimer, SysTimerId );
    while ( ( hwRdReg( baseAddress, statRegOffs ) & statRegRdy ) == 0 )
    {
       if ( ( SysTimer->ReadTimerLongSec( SysTimer, SysTimerId ) - f ) > 1.1 * HW_READ_BRAM_TIMEOUT ) return 0;
    }
    u = hwRdReg( baseAddress, dataRegOffs );
    hwWriteCtrlReg( baseAddress, ctrlRegOffs, ctrlRegAck );
    return u;
}

u64 hwReadData64( INTPTR baseAddress,
               u32 ctrlRegOffs, u32 ctrlRegRd, u32 ctrlRegAck,
               u32 statRegOffs, u32 statRegRdy,
               u32 addrRegOffs, u32 addr,
               u32 dataRegOffs )
{
    hwWrReg( baseAddress, addrRegOffs, addr );
    hwWriteCtrlReg( baseAddress, ctrlRegOffs, ctrlRegRd );
    double f = SysTimer->ReadTimerLongSec( SysTimer, SysTimerId );
    while ( ( hwRdReg( baseAddress, statRegOffs ) & statRegRdy ) == 0 )
    {
       if ( ( SysTimer->ReadTimerLongSec( SysTimer, SysTimerId ) - f ) > 1.1 * HW_READ_BRAM_TIMEOUT ) return 0;
    }
    u64 u = hwRdReg64( baseAddress, dataRegOffs );
    hwWriteCtrlReg( baseAddress, ctrlRegOffs, ctrlRegAck );
    return ( u );
}

// beware - non atomic
u64 hwRdReg64( u32 addr, u32 offs )
{
    u64 lo = hwRdReg( addr, offs + 0 ) & 0xFFFFFFFF;
    u64 hi = hwRdReg( addr, offs + 4 ) & 0xFFFFFFFF;
    return ( ( hi << 32 ) | lo );
}

// beware - non atomic
void hwWrReg64( u32 addr, u32 offs, u64 data )
{
    u64 dh    = ( data >> 32 ) & 0xFFFFFFFF;
    u32 datah = ( u32 ) dh;
    u64 dl    = data & 0xFFFFFFFF;
    u32 datal = ( u32 ) dl;
    hwWrReg( addr, offs + 0, datal );
    hwWrReg( addr, offs + 4, datah );
}

u32 hwRdReg( u32 addr, u32 offs )
{
    u32 rp = 0;
    if ( addr < SYS_MAX_DEVICES )
    {
        rp = SysDvceSimHwReg( addr, offs >> 2, 0, false );
    }
    else
    {
        rp = hwReadReg( addr, offs );
    }
    return ( rp );
}

void hwWrReg( u32 addr, u32 offs, u32 data )
{
    if ( addr < SYS_MAX_DEVICES )
    {
        ( void ) SysDvceSimHwReg( addr, offs >> 2, data, true );
    }
    else
    {
        hwWriteReg( addr, offs, data );
    }
}

