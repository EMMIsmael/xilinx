/*
 * SysDevice.c

 *
 *  Created on: 10 Aug 2018
 *      Author: JULIAN MORTIMER
 */
#include "include.h"

static SysDvceInfo SysDvceRegistry[ SYS_MAX_DEVICES ] =
        {
        { 0, NULL },
        { 0, NULL },
        { 0, NULL },
        { 0, NULL },
        { 0, NULL },
        { 0, NULL },
        { 0, NULL },
        { 0, NULL },
        { 0, NULL },
        { 0, NULL },
        { 0, NULL },
        { 0, NULL },
        { 0, NULL },
        { 0, NULL },
        { 0, NULL },
        { 0, NULL },

        { 0, NULL },
        { 0, NULL },
        { 0, NULL },
        { 0, NULL },
        { 0, NULL },
        { 0, NULL },
        { 0, NULL },
        { 0, NULL },
        { 0, NULL },
        { 0, NULL },
        { 0, NULL },
        { 0, NULL },
        { 0, NULL },
        { 0, NULL },
        { 0, NULL },
        { 0, NULL },

        { 0, NULL },
        { 0, NULL },
        { 0, NULL },
        { 0, NULL },
        { 0, NULL },
        { 0, NULL },
        { 0, NULL },
        { 0, NULL },
        { 0, NULL },
        { 0, NULL },
        { 0, NULL },
        { 0, NULL },
        { 0, NULL },
        { 0, NULL },
        { 0, NULL },
        { 0, NULL },

        { 0, NULL },
        { 0, NULL },
        { 0, NULL },
        { 0, NULL },
        { 0, NULL },
        { 0, NULL },
        { 0, NULL },
        { 0, NULL },
        { 0, NULL },
        { 0, NULL },
        { 0, NULL },
        { 0, NULL },
        { 0, NULL },
        { 0, NULL },
        { 0, NULL },
        { 0, NULL }
};

static u32 SysDvceRegPtr = 0;

static SysDvce *SysDvceRegistrySrch( u32 baseAddr );
static _Bool RegisterSysDvce( SysDvce *pd );
u32 SysDvceSimHwReg( u32 baseAddr, u32 reg, u32 data, _Bool bWrite );

static void SyncRdSlv( SysDvce *pd );

void SysDvceInit(
		SysDvce  *pd,
		const char *Name,
		SysCtrl  *pSysCtrl,
		SysDesc  *pSysDesc,
		SysPtrs  *pSysPtrs,
		SysData  *pSysData,
		u32      ( *SimHwReg ) ( u32, u32, u32, _Bool ),
		u32      baseAddr,
		u32      HwInfoReg,
		u32      CtrlReg,
		u32      InfoReg,
		u32      StatReg,
		u32      ResetBit,
		u32      IntAckBit,
		u32      RdSyncBit,
	    u32      RxDataCtrReg,
	    u32      RxSymbCtrReg,
	    u32      RxFrameCtrReg,
	    u32      RxBlockCtrReg,
	    u32      TxDataCtrReg,
	    u32      TxSymbCtrReg,
	    u32      TxFrameCtrReg,
	    u32      TxBlockCtrReg,
	    u32      RxTVCtrReg,
	    u32      TxTVCtrReg
		);

static inline void Reset         ( SysDvce *pd ) { hwWriteCtrlReg( pd->BaseAddr, pd->CtrlReg << 2, pd->ResetBit ); }
static inline u32  ReadStatReg   ( SysDvce *pd ) { return( hwRdReg( pd->BaseAddr, pd->StatReg << 2 ));       }
static inline u32  ReadSlvReg    ( SysDvce *pd, u32 reg ) { SyncRdSlv( pd ); return( hwRdReg   ( pd->BaseAddr, reg << 2 ));   }
static inline u64  ReadSlvReg64  ( SysDvce *pd, u32 reg ) { SyncRdSlv( pd ); return( hwRdReg64 ( pd->BaseAddr, reg << 2 ));   }
static inline void WriteSlvReg   ( SysDvce *pd, u32 reg, u32 u )  { hwWrReg   ( pd->BaseAddr, reg << 2, u ); }
static inline void WriteSlvReg64 ( SysDvce *pd, u32 reg, u64 u )  { hwWrReg64 ( pd->BaseAddr, reg << 2, u ); }
static inline u32  NSlvReg       ( SysDvce *pd ) { return( pd->ReadSlvReg( pd, pd->InfoReg ) & PL_DEVICE_N_SLV_REG_MASK ); }
static inline u32  DeviceId      ( SysDvce *pd ) { return( pd->ReadSlvReg( pd, pd->InfoReg ) >> PL_DEVICE_ID_LSB ); }
static inline void IntAck        ( SysDvce *pd ) { hwWriteCtrlReg( pd->BaseAddr, pd->CtrlReg << 2, pd->IntAckBit ); }
static void        Dump          ( SysDvce *pd );
static void        DumpSlvReg    ( SysDvce *pd, _Bool bBin );
static inline void SyncRdSlv     ( SysDvce *pd ) { hwWriteCtrlReg( pd->BaseAddr, pd->CtrlReg << 2, pd->RdSyncBit ); }


SysDvce SysDvceInitStr = {
        "",   // char                  Name[ SYS_MAX_OBJ_NAME_LEN + 1 ];
        NULL, // const SysCtrl         *pSysCtrl;
        NULL, // const SysDesc         *pSysDesc;
        NULL, // const SysPtrs         *pSysPtrs;
        NULL, // SysData               *pSysData;
        0,    // u32                   BaseAddr;
        0,    // u32                   HwInfoReg;
        0,    // u32                   CtrlReg;
        0,    // u32                   InfoReg;
        0,    // u32                   StatReg;
        0,    // u32                   ResetBit;
        0,    // u32                   IntAckBit;
        0,    // u32                   RdSyncBit;
        0,    // u32                   RxDataCtrReg;
        0,    // u32                   RxSymbCtrReg;
        0,    // u32                   RxFrameCtrReg;
        0,    // u32                   RxBlockCtrReg;
        0,    // u32                   TxDataCtrReg;
        0,    // u32                   TxSymbCtrReg;
        0,    // u32                   TxFrameCtrReg;
        0,    // u32                   TxBlockCtrReg;
        0,    // u32                   RxTVCtrReg;
        0,    // u32                   TxTVCtrReg;
        {
                0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0,     // u32                   SimSlvRegs[ SYS_N_SLV_REG_MAX ];
                0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0
        },
        Reset         ,  // void (* Reset         ) ( SysDvce *pd );
        ReadStatReg   ,  // u32  (* ReadStatReg   ) ( SysDvce *pd );
        ReadSlvReg    ,  // u32  (* ReadSlvReg    ) ( SysDvce *pd, u32 reg );
        ReadSlvReg64  ,  // u64  (* ReadSlvReg64  ) ( SysDvce *pd, u32 reglo, u32 reghi );
        WriteSlvReg   ,  // void (* WriteSlvReg   ) ( SysDvce *pd, u32 reg, u32 data );
        WriteSlvReg64 ,  // void (* WriteSlvReg64 ) ( SysDvce *pd, u32 reg, u64 data );
        NSlvReg       ,  // u32  (* NSlvReg       ) ( SysDvce *pd );
        DeviceId      ,  // u32  (* DeviceId      ) ( SysDvce *pd );
        IntAck        ,  // void (* IntAck        ) ( SysDvce *pd );
        Dump          ,  // void (* Dump          ) ( SysDvce *pd );
        DumpSlvReg    ,  // void (* DumpSlvReg    ) ( SysDvce *pd, _Bool bBin );
        NULL             // u32  (* SimHwReg      ) ( u32 baseAddr, u32 reg, u32 data, _Bool bWrite );
};


void SysDvceInit(
		SysDvce  *pd,
		const char *Name,
		SysCtrl  *pSysCtrl,
		SysDesc  *pSysDesc,
		SysPtrs  *pSysPtrs,
		SysData  *pSysData,
        u32      ( *SimHwReg ) ( u32, u32, u32, _Bool ),
		u32      BaseAddr,
		u32      HwInfoReg,
		u32      CtrlReg,
		u32      InfoReg,
		u32      StatReg,
		u32      ResetBit,
		u32      IntAckBit,
        u32      RdSyncBit,
	    u32      RxDataCtrReg,
	    u32      RxSymbCtrReg,
	    u32      RxFrameCtrReg,
	    u32      RxBlockCtrReg,
	    u32      TxDataCtrReg,
	    u32      TxSymbCtrReg,
	    u32      TxFrameCtrReg,
	    u32      TxBlockCtrReg,
        u32      RxTVCtrReg,
        u32      TxTVCtrReg
		)
{
    *pd             = SysDvceInitStr;
	strcpy( pd->Name, Name );
	pd->BaseAddr    = BaseAddr;
	pd->HwInfoReg   = HwInfoReg;
	pd->CtrlReg     = CtrlReg;
	pd->InfoReg     = InfoReg;
	pd->StatReg     = StatReg;
	pd->ResetBit    = ResetBit;
	pd->IntAckBit   = IntAckBit;
    pd->RdSyncBit   = RdSyncBit;
	pd->pSysCtrl    = pSysCtrl;
	pd->pSysDesc    = pSysDesc;
	pd->pSysPtrs    = pSysPtrs;
	pd->pSysData    = pSysData;
	pd->SimHwReg    = SimHwReg;
	SysDataInit(
		pd,
		RxDataCtrReg,
		RxSymbCtrReg,
		RxFrameCtrReg,
		RxBlockCtrReg,
		TxDataCtrReg,
		TxSymbCtrReg,
		TxFrameCtrReg,
		TxBlockCtrReg,
        RxTVCtrReg,
        TxTVCtrReg
	);
	RegisterSysDvce( pd );
}

static void DumpSlvReg( SysDvce *pd, _Bool bBin )
{
	char Name[ SYS_MAX_OBJ_NAME_LEN + 1 ]; strcpy( Name, pd->Name );
	SyncRdSlv( pd );
	if ( bBin )
		utilScreenDumpData( 0, DUMP_FMT_BIN32  , DUMP_TYP_SLV_REG, pd->BaseAddr, pd->NSlvReg( pd ), 0, 0, strcat( Name, " Slave Registers:" ) );
	else
		utilScreenDumpData( 0, DUMP_FMT_HEX_DEC, DUMP_TYP_SLV_REG, pd->BaseAddr, pd->NSlvReg( pd ), 0, 0, strcat( Name, " Slave Registers:" ) );

}

static void Dump( SysDvce *pd )
{
	printf("\nName             : %s"    , pd->Name                       ); fflush( stdout );
	printf("\nBase address     : 0x%08X", toUint( pd->BaseAddr )         ); fflush( stdout );
	printf("\nDevice ID        : 0x%06X", toUint( pd->DeviceId( pd ) )   ); fflush( stdout );
	printf("\nSlave regs       : %u"    , toUint( pd->NSlvReg( pd ) )    ); fflush( stdout );
	printf("\n"); fflush( stdout );
}



static _Bool RegisterSysDvce( SysDvce *pd )
{
    _Bool rc = false;
    if ( SysDvceRegPtr < SYS_MAX_DEVICES )
    {
        if ( SysDvceRegistrySrch( pd->BaseAddr ) == NULL )
        {
            SysDvceRegistry [ SysDvceRegPtr ].BaseAddr = pd->BaseAddr;
            SysDvceRegistry [ SysDvceRegPtr ].pSysDvce = pd;
            SysDvceRegPtr++;
            rc = true;
        }
    }
    return ( rc );
}


SysDvce *SysDvceRegistrySrch( u32 baseAddr )
{
    SysDvce *ps = NULL;
    for ( int i = 0; i < SYS_MAX_DEVICES; i++ )
    {
        if ( SysDvceRegistry[ i ].BaseAddr == baseAddr )
        {
            ps = SysDvceRegistry[ i ].pSysDvce;
            break;
        }
    }
    return ( ps );
}

// Simulate HW reg - for testing when device not in HW spec
// To enable simulation, assign a unique base address below
// 0x40000000
// If a SimHwReg address is supplied by the parent device object
// it is used, otherwise the system defaults to the code shown here
u32 SysDvceSimHwReg( u32 baseAddr, u32 reg, u32 data, _Bool bWrite )
{
    u32 rp = 0;
    SysDvce *pd = SysDvceRegistrySrch( baseAddr );
    if ( pd == NULL )
    {
        rp = 0x0DDBA115;
    }
    else if ( pd->SimHwReg == NULL )
    {
        if ( bWrite )
        {
            pd->SimSlvRegs[ reg ] = data;
            rp = data;
        }
        else
        {
            if ( reg == pd->InfoReg )
            {
                rp = pd->HwInfoReg;
            }
            else if ( reg == pd->StatReg )
            {
                rp = -1;
            }
            else
            {
                rp = pd->SimSlvRegs [ reg ];
            }
        }
    }
    else
    {
        rp = pd->SimHwReg( baseAddr, reg, data, bWrite );
    }
    return ( rp );
}
