/*
 * GM_Simple_Timer.c
 *
 *  Created on: 17 May 2017
 *      Author: JULIAN MORTIMER
 */
#include "include.h"
#include "Device_Manager.h"

#define BIT_TIMER_RESET      ( 0 )
#define BIT_TIMER_INT_ACK    ( 1 )
#define BIT_TIMER_START      ( 2 )
#define BIT_TIMER_STOP       ( 3 )
#define BIT_TIMER_CLEAR      ( 4 )
#define BIT_TIMER_READ       ( 5 )
#define BIT_TIMER_READ_CLEAR ( 6 )

#define TIMER_CTRL_REG       ( 0 )
#define TIMER_INFO_REG       ( 1 )

#define TIMER_CTR_LO_REG     ( 3 )
#define TIMER_CTR_HI_REG     ( 4 )


typedef union timerInfo TimerInfo;
union timerInfo {
	u32 slvReg;
	struct {
      u8  nSlvReg;
      u8  RWbits;
      u16 devId;
	};
};

typedef struct sTimer Timer;
struct sTimer {
	char   name[SYS_MAX_OBJ_NAME_LEN + 1];
	u32    baseAddr;
	u32    id;
	u32    (* IsRunning )( Timer * );
    void   (* Reset )( Timer * );
	void   (* IntAck )( Timer * );
	void   (* Start )( Timer * );
	void   (* Stop )( Timer * );
	void   (* Clear )( Timer * );
	u32    (* ReadLoWord )( Timer * );
	u32    (* ReadHiWord )( Timer * );
	u32    (* ReadClear )( Timer * );
	double (* ReadLong )( Timer * );
	void   (* ReadSlvRegs )( Timer *p, u32 *pData );
	u32    (* ReadSlvReg )( Timer *p, u32 reg );
	void   (* WriteSlvReg )( Timer *p, u32 reg, u32 data );
	u32    (* GetDeviceId )( Timer * p );
	u32    (* GetNSlvReg )( Timer * p );
	u32    (* BaseAddr  )( Timer * );
	u32    (* OfsCtrlReg  )( Timer * );
	u32    (* OfsInfoReg  )( Timer * );
	u32    (* OfsCtrLoReg  )( Timer * );
	u32    (* OfsCtrHiReg  )( Timer * );
	u32    (* BitCtrlReset  )( Timer * );
	u32    (* BitCtrlIntAck  )( Timer * );
	u32    (* BitCtrlStart  )( Timer * );
	u32    (* BitCtrlStop  )( Timer * );
	u32    (* BitCtrlRead  )( Timer * );
	u32    (* BitCtrlClear  )( Timer * );
	u32    (* BitCtrlReadClear )( Timer * );
};


// API-accessed management
// API (via parent) and DeviceManager-accessed
static DeviceErrCode InitTimerAPI( TimerAPI *pa );
static Device*       AllocateTimer( DeviceManager *pm );
static u32           InitTimer( TimerAPI *pa, u32 baseAddr, u32 id, const char name[], u32 deviceId );

static DeviceErrCode DevManErr( TimerAPI *pa );
static u32        DevManInitialised( TimerAPI *pa );
static const char *DevManErrMsg( TimerAPI *pa );
static const char *GetTimerName( TimerAPI *pa, u32 id );
static void       ResetTimer( TimerAPI *pa, u32 id );
static void       StartTimer( TimerAPI *pa, u32 id );
static void       StopTimer( TimerAPI *pa, u32 id );
static u32        ReadTimerHi( TimerAPI *pa, u32 id );
static u32        ReadTimerLo( TimerAPI *pa, u32 id );
static void       ClearTimer( TimerAPI *pa, u32 id );
static double     ReadTimerLong( TimerAPI *pa, u32 id );
static double     ReadTimerLongSec( TimerAPI *pa, u32 id );
static u32        ReadClearTimer( TimerAPI *pa, u32 id );
static u32        TimerIsRunning( TimerAPI *pa, u32 id );
static void       ReadTimerSlvRegs( TimerAPI *pa, u32 id, u32 *pData );
static void       WriteTimerSlvReg(  TimerAPI *pa, u32 id, u32 reg, u32 data );
static u32        ReadTimerSlvReg(  TimerAPI *pa, u32 id, u32 reg );
static u32        GetTimerNSlvReg(  TimerAPI *pa, u32 id  );
static u32        GetTimerDeviceId(  TimerAPI *pa, u32 id );
static void       DumpTimerSlvRegs( TimerAPI *pa, u32 id, u32 fmtsel );

// Timer accessed by pointer from this module only
static u32        IsRunning( Timer *p );
static void       Reset( Timer *p );
static void       IntAck( Timer *p );
static void       Start( Timer *p );
static void       Stop( Timer *p );
static void       Clear( Timer *p );
static u32        ReadLoWord( Timer *p );
static u32        ReadHiWord( Timer *p );
static u32        ReadClear( Timer *p );
static double     ReadLong( Timer *p );
static void       ReadSlvRegs( Timer *p, u32 *pData );
static void       WriteSlvReg( Timer *p, u32 reg, u32 data );
static u32        ReadSlvReg( Timer *p, u32 reg );
static u32        GetNSlvReg( Timer * p );
static u32        GetDeviceId( Timer * p );

static inline u32 BitCtrlReset( Timer *p )     { return( BIT_TIMER_RESET);        }
static inline u32 BitCtrlIntAck( Timer *p )    { return( BIT_TIMER_INT_ACK );     }
static inline u32 BitCtrlStart( Timer *p )     { return( BIT_TIMER_START );       }
static inline u32 BitCtrlStop( Timer *p )      { return( BIT_TIMER_STOP );        }
static inline u32 BitCtrlClear( Timer *p )     { return( BIT_TIMER_CLEAR );       }
static inline u32 BitCtrlRead( Timer *p )      { return( BIT_TIMER_READ );        }
static inline u32 BitCtrlReadClear( Timer *p ) { return( BIT_TIMER_READ_CLEAR );  }

static inline u32 BaseAddr( Timer *p )         { return( p->baseAddr );           }
static inline u32 OfsCtrlReg( Timer *p )       { return( TIMER_CTRL_REG   << 2 ); }
static inline u32 OfsInfoReg( Timer *p )       { return( TIMER_INFO_REG   << 2 ); }
static inline u32 OfsCtrLoReg( Timer *p )      { return( TIMER_CTR_LO_REG << 2 ); }
static inline u32 OfsCtrHiReg( Timer *p )      { return( TIMER_CTR_HI_REG << 2 ); }

static Timer LocalTimer;
static Timer *TimerArray[SYS_MAX_DEVICES];
static DeviceManager TimerManager;

static const Timer TIMER_TEMPLATE = {
  "\0",
  0,
  0,
  IsRunning,
  Reset,
  IntAck,
  Start,
  Stop,
  Clear,
  ReadLoWord,
  ReadHiWord,
  ReadClear,
  ReadLong,
  ReadSlvRegs,
  ReadSlvReg,
  WriteSlvReg,
  GetDeviceId,
  GetNSlvReg,
  BaseAddr,
  OfsCtrlReg,
  OfsInfoReg,
  OfsCtrLoReg,
  OfsCtrHiReg,
  BitCtrlReset,
  BitCtrlIntAck,
  BitCtrlStart,
  BitCtrlStop,
  BitCtrlRead,
  BitCtrlClear,
  BitCtrlReadClear
};


static TimerAPI _SysTimerAPI_ = {
	&TimerManager,
	DevManErr,
	DevManErrMsg,
	DevManInitialised,
	InitTimerAPI,
	InitTimer,
	GetTimerName,
	ResetTimer,
	StartTimer,
	StopTimer,
    ClearTimer,
    ReadTimerLo,
    ReadTimerHi,
    ReadClearTimer,
    ReadTimerLong,
    ReadTimerLongSec,
    TimerIsRunning,
    ReadTimerSlvRegs,
    WriteTimerSlvReg,
    ReadTimerSlvReg,
    GetTimerNSlvReg,
    GetTimerDeviceId,
    DumpTimerSlvRegs,
};
TimerAPI *SysTimerAPI = &_SysTimerAPI_;

// =======================
// Debugging Aids
// =======================
static void DumpTimerSlvRegs( TimerAPI *pa, u32 id, u32 fmtsel )
{
	char title[SYS_MAX_OBJ_NAME_LEN + 50];
	u32 buf[SYS_TIMER_N_SLV_REGS];
	pa->ReadTimerSlvRegs( pa, id, buf );
	u32 fmt = (fmtsel == 0 ? DUMP_FMT_HEX_DEC : DUMP_FMT_BIN32 );
	sprintf(title, "Name: %s (Device Id: %08X, Slave Registers: %d)", pa->Name( pa, id ), (unsigned int)pa->DeviceId(pa, id ), (unsigned int)pa->NSlvReg(pa, id ));
	utilScreenDumpData(0, fmt, DUMP_TYP_SLV_REG, ( u32 )buf, pa->NSlvReg(pa, id ), 0, 0, title);
	dispSysTimerRunning();
}

// =======================
// Timer Manager Overloads
// =======================
// API user can name the device, this returns an id.
// User should provide an initialised LUT.

static Device* AllocateTimer( DeviceManager *pm )
{
	Timer *pd  = ( Timer* ) malloc(sizeof( Timer ));
	if( !pd ) return ( Device * )NULL;
    *pd = TIMER_TEMPLATE;
	return ( Device *) pd;
}

static DeviceErrCode InitTimerAPI( TimerAPI *pa )
{
   // Insert user code here
   LocalTimer         = TIMER_TEMPLATE;
   DeviceAPI *pd      = (DeviceAPI *)pa;
   DeviceManager *pm  = (DeviceManager *)pa->pDeviceManager;
   *pm                = DEVICE_MANAGER_TEMPLATE;
   pm->AllocateDevice = AllocateTimer;
   pm->DeviceArray    = (Device **) TimerArray;
   pm->pAPI           = pd;
   return pm->InitDeviceManager( pm );
}


static u32 InitTimer( TimerAPI *pa, u32 baseAddr, u32 id, const char name[], u32 deviceId)
{
	Timer *pd;
	DeviceManager *pm = (DeviceManager *)pa->pDeviceManager;
	if ( !pm->IsInitialised ) pa->InitTimerAPI( pa );
	if (strlen(name) > SYS_MAX_OBJ_NAME_LEN)
    {
        pm->SetDeviceErr( pm, DEVICE_NAME_TOO_LONG );
    	return ( ID_DEVICE_ERROR );
    }


	if( id == ID_NEW_DEVICE ) // Fresh device
	{
	    // Check for existence of correct HW device
		pd           = &LocalTimer;
	    pd->baseAddr = baseAddr;

	    if ( pd->GetDeviceId( pd ) != deviceId )
	    {
	    	pm->SetDeviceErr(pm, DEVICE_BAD_DEVICE_ID );
	    	return ID_DEVICE_ERROR;
	    }

	    pd = (Timer *)pm->AddDevice( pm );
		if ( !pd ) return( ID_DEVICE_ERROR );
	}
	else
	{
        // here to reinitialise
		if((pd = (Timer *)pm->GetDevicePtr( pm, id )) == NULL)
	    {
		   pm->SetDeviceErr( pm, DEVICE_NOT_FOUND );
		   return ( ID_DEVICE_ERROR );
	    }
	}
    pd->baseAddr = baseAddr;

    strcpy( pd->name, name );
    pd->Reset( pd );
    return( pd->id );
}


static DeviceErrCode DevManErr( TimerAPI *pa )
{
	DeviceManager *pm = (DeviceManager *)pa->pDeviceManager;
    return pm->GetDeviceErr( pm );
}

static const char* DevManErrMsg( TimerAPI *pa )
{
	DeviceManager *pm = (DeviceManager *)pa->pDeviceManager;
    return pm->GetDeviceErrMsg( pm );
}

static u32 DevManInitialised( TimerAPI *pa )
{
	DeviceManager *pm = (DeviceManager *)pa->pDeviceManager;
    return pm->IsInitialised( pm );
}

//================================
// API-accessible device functions
//================================

const char *GetTimerName(  TimerAPI *pa, u32 id )
{
	DeviceManager *pm = (DeviceManager *)pa->pDeviceManager;
	Timer *pt =         (Timer *)pm->GetDevicePtr( pm, id );
	if ( pt )
	{
		return( pt->name );
	}
	else
	{
		pm->SetDeviceErr( pm, DEVICE_NOT_FOUND );
	    return 0;
	}
}

u32 TimerIsRunning( TimerAPI *pa, u32 id )
{
	DeviceManager *pm = (DeviceManager *)pa->pDeviceManager;
	Timer *pt =         (Timer *)pm->GetDevicePtr( pm, id );
	if ( pt )
	{
		return( pt->IsRunning( pt ));
	}
	else
	{
		pm->SetDeviceErr( pm, DEVICE_NOT_FOUND );
	    return 0;
	}
}

void ResetTimer( TimerAPI *pa, u32 id )
{
	DeviceManager *pm = (DeviceManager *)pa->pDeviceManager;
	Timer *pt =         (Timer *)pm->GetDevicePtr( pm, id );
	if ( pt )
	{
		pt->Reset( pt );
	}
	else
	{
		pm->SetDeviceErr( pm, DEVICE_NOT_FOUND );
	}
}

void StartTimer( TimerAPI *pa, u32 id )
{
	DeviceManager *pm = (DeviceManager *)pa->pDeviceManager;
	Timer *pt =         (Timer *)pm->GetDevicePtr( pm, id );
	if ( pt )
	{
		pt->Start( pt );
	}
	else
	{
		pm->SetDeviceErr( pm, DEVICE_NOT_FOUND );
	}
}

void StopTimer( TimerAPI *pa, u32 id )
{
	DeviceManager *pm = (DeviceManager *)pa->pDeviceManager;
	Timer *pt =         (Timer *)pm->GetDevicePtr( pm, id );
	if ( pt )
	{
		pt->Stop( pt );
	}
	else
	{
		pm->SetDeviceErr( pm, DEVICE_NOT_FOUND );
	}
}

u32 ReadTimerHi( TimerAPI *pa, u32 id )
{
	DeviceManager *pm = (DeviceManager *)pa->pDeviceManager;
	Timer *pt =         (Timer *)pm->GetDevicePtr( pm, id );
	if ( pt )
	{
		return ( pt->ReadHiWord( pt ));
	}
	else
	{
		pm->SetDeviceErr( pm, DEVICE_NOT_FOUND );
	    return 0;
	}
}

u32 ReadTimerLo( TimerAPI *pa, u32 id )
{
	DeviceManager *pm = (DeviceManager *)pa->pDeviceManager;
	Timer *pt =         (Timer *)pm->GetDevicePtr( pm, id );
	if ( pt )
	{
		return ( pt->ReadLoWord( pt ));
	}
	else
	{
		pm->SetDeviceErr( pm, DEVICE_NOT_FOUND );
	    return 0;
	}
}

void ClearTimer( TimerAPI *pa, u32 id )
{
	DeviceManager *pm = (DeviceManager *)pa->pDeviceManager;
	Timer *pt =         (Timer *)pm->GetDevicePtr( pm, id );
	if ( pt )
	{
		pt->Clear( pt );
	}
	else
	{
		pm->SetDeviceErr( pm, DEVICE_NOT_FOUND );
	}
}

double ReadTimerLong( TimerAPI *pa, u32 id )
{
	DeviceManager *pm = (DeviceManager *)pa->pDeviceManager;
	Timer *pt =         (Timer *)pm->GetDevicePtr( pm, id );
	if ( pt )
	{
		return ( pt->ReadLong( pt ));
	}
	else
	{
		pm->SetDeviceErr( pm, DEVICE_NOT_FOUND );
		return( 0.0 );
	}
}

double ReadTimerLongSec( TimerAPI *pa, u32 id )
{
	return ( pa->ReadTimerLong( pa, id ) / PL_FABRIC_CLOCK_FREQ );
}

u32 ReadClearTimer( TimerAPI *pa, u32 id )
{
	DeviceManager *pm = (DeviceManager *)pa->pDeviceManager;
	Timer *pt =         (Timer *)pm->GetDevicePtr( pm, id );
	if ( pt )
	{
		return ( pt->ReadClear( pt ));
	}
	else
	{
		pm->SetDeviceErr( pm, DEVICE_NOT_FOUND );
	    return 0;
	}
}

static void ReadTimerSlvRegs( TimerAPI *pa, u32 id, u32 *pData )
{
	DeviceManager *pm = (DeviceManager *)pa->pDeviceManager;
	Timer *pt =         (Timer *)pm->GetDevicePtr( pm, id );
	if ( pt )
	{
		pt->ReadSlvRegs( pt, pData );
	}
	else
	{
		pm->SetDeviceErr( pm, DEVICE_NOT_FOUND );
	}
}

static u32 ReadTimerSlvReg(  TimerAPI *pa, u32 id, u32 reg )
{
	DeviceManager *pm = (DeviceManager *)pa->pDeviceManager;
	Timer *pt =         (Timer *)pm->GetDevicePtr( pm, id );
	if ( pt )
	{
		return ( pt->ReadSlvReg( pt, reg ));
	}
	else
	{
		pm->SetDeviceErr( pm, DEVICE_NOT_FOUND );
	    return 0;
	}
}

static void WriteTimerSlvReg(  TimerAPI *pa, u32 id, u32 reg, u32 data )
{
	DeviceManager *pm = (DeviceManager *)pa->pDeviceManager;
	Timer *pt =         (Timer *)pm->GetDevicePtr( pm, id );
	if ( pt )
	{
		pt->WriteSlvReg( pt, reg, data );
	}
	else
	{
		pm->SetDeviceErr( pm, DEVICE_NOT_FOUND );
	}
}

static u32 GetTimerNSlvReg(  TimerAPI *pa, u32 id  )
{
	DeviceManager *pm = (DeviceManager *)pa->pDeviceManager;
	Timer *pt =         (Timer *)pm->GetDevicePtr( pm, id );
	if ( pt )
	{
		return ( pt->GetNSlvReg( pt ));
	}
	else
	{
		pm->SetDeviceErr( pm, DEVICE_NOT_FOUND );
	    return 0;
	}
}

static u32 GetTimerDeviceId(  TimerAPI *pa, u32 id )
{
	DeviceManager *pm = (DeviceManager *)pa->pDeviceManager;
	Timer *pt =         (Timer *)pm->GetDevicePtr( pm, id );
	if ( pt )
	{
		return ( pt->GetDeviceId( pt ));
	}
	else
	{
		pm->SetDeviceErr( pm, DEVICE_NOT_FOUND );
	    return 0;
	}
}


// Access allowed only from this module.
static u32 IsRunning( Timer *p )
{
	return( p->ReadLoWord( p ) != p->ReadLoWord( p ) );
}

static void Reset( Timer *p )
{
   hwWriteCtrlReg( p->BaseAddr( p ), p->OfsCtrlReg( p ), 1 << p->BitCtrlReset( p ));
}

static void IntAck( Timer *p )
{
   hwWriteCtrlReg( p->BaseAddr( p ), p->OfsCtrlReg( p ), 1 << p->BitCtrlIntAck( p ));
}

static void Start( Timer *p )
{
   hwWriteCtrlReg( p-> BaseAddr( p ), p->OfsCtrlReg( p ), 1 << p->BitCtrlStart( p ));
}

static void Stop( Timer *p )
{
   hwWriteCtrlReg( p->BaseAddr( p ), p->OfsCtrlReg( p ), 1 << p->BitCtrlStop( p ));
}

static void Clear( Timer *p )
{
   hwWriteCtrlReg( p->BaseAddr( p ), p->OfsCtrlReg( p ), 1 << p->BitCtrlClear( p ));
}

static u32 ReadClear( Timer *p )
{
   hwWriteCtrlReg( p->BaseAddr( p ), p->OfsCtrlReg( p ), 1 << p->BitCtrlReadClear( p ));
   return hwReadReg( p->BaseAddr( p ), p->OfsCtrLoReg( p ));
}

static u32 ReadHiWord( Timer *p )
{
   hwWriteCtrlReg( p->BaseAddr( p ), p->OfsCtrlReg( p ), 1 << p->BitCtrlRead( p ));
   return hwReadReg( p->BaseAddr( p ), p->OfsCtrHiReg( p ));
}

static u32 ReadLoWord( Timer *p )
{
   hwWriteCtrlReg( p->BaseAddr( p ), p->OfsCtrlReg( p ), 1 << p->BitCtrlRead( p ));
   return hwReadReg( p->BaseAddr( p ), p->OfsCtrLoReg( p ));
}

static double ReadLong( Timer *p )
{
   hwWriteCtrlReg( p->BaseAddr( p ), p->OfsCtrlReg( p ), 1 << p->BitCtrlRead( p ));
   return (double)hwReadReg( p->BaseAddr( p ), p->OfsCtrHiReg( p )) * 4294967296.0
	    + (double)hwReadReg( p->BaseAddr( p ), p->OfsCtrLoReg( p ));
}

static void ReadSlvRegs( Timer *p, u32 *pData )
{
	u32 i;
	for (i = 0; i < p->GetNSlvReg( p ); i++ )
	{
		pData[i] = hwReadReg( p->BaseAddr( p ), i << 2 );
	}
}

static u32 ReadSlvReg( Timer *p, u32 reg )
{
	return hwReadReg( p->BaseAddr( p ), reg << 2);
}

static void WriteSlvReg( Timer *p, u32 reg, u32 data )
{
	return hwWriteReg( p->BaseAddr( p ), reg << 2, data );
}

static u32 GetNSlvReg( Timer * p )
{
    TimerInfo info;
    info.slvReg = hwReadReg( p->BaseAddr( p ),  p->OfsInfoReg( p ));
    return( ( u32 ) info.nSlvReg );
}

static u32 GetDeviceId( Timer * p )
{
    return (hwReadReg( p->BaseAddr( p ),  p->OfsInfoReg( p )));
}


