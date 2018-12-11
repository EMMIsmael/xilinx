/*
 * TREmulSw.c
 *
 *  Created on: 2 Aug 2018
 *      Author: JULIAN MORTIMER
 */
#include "include.h"

// =================================
// Boilerplate Function Declarations
// =================================
void TREmulSwInit(
    TREmulSw *pd,
    const char      *Name,
    u32             BaseAddr,
    SysCtrl         *pSysCtrl,
    SysDvce         *pSysDvce,
    SysDesc         *pSysDesc,
    SysPtrs         *pSysPtrs,
    SysData         *pSysData,
    TREmulSwDesc *pDvcDesc,
    TREmulSwData *pDvcData
    );
static void SetParams    ( TREmulSw *pd );
static void SysDvceDump  ( TREmulSw *pd );
static void SysDescDump  ( TREmulSw *pd );
static void SysDataDump  ( TREmulSw *pd );
static void DvcDataInit  ( TREmulSw *pd );
static void DvcDataUpd   ( TREmulSw *pd );
static void DvcDescDump  ( TREmulSw *pd );
static void DvcDataDump  ( TREmulSw *pd );

// =================================================
// Device Structure Member Functions (public access)
// =================================================
#define BSET( u, mask, b ) do { if( b ) u |= mask; else u &= ~mask; } while( 0 )
static inline void SelChan           ( TREmulSw *pd, u32 SelEmul ) { DVCPAR( ConfigReg ) = SelEmul; SYSWRSLVREG( pd, TREMULSW_CONFIG_REG, DVCPAR( ConfigReg )); }
static void        Dump              ( TREmulSw *pd );

// Declare init structure here with member functions
TREmulSw TREmulSwInitStr = {
	NULL,
	NULL,
	NULL,
	SelChan,
	Dump
};

//===================
// Initialize Device
//===================
void TREmulSwInit( TREmulSw *pd, const char *Name, u32 BaseAddr, SysCtrl *pSysCtrl, SysDvce *pSysDvce, SysDesc *pSysDesc, SysPtrs *pSysPtrs, SysData *pSysData, TREmulSwDesc *pDvcDesc, TREmulSwData *pDvcData )
{
    *pd             = TREmulSwInitStr;
    pd->pSysDvce    = pSysDvce;
    SysDvceInit(
        pd->pSysDvce,
        Name,
        pSysCtrl,
        pSysDesc,
        pSysPtrs,
        pSysData,
        NULL,
        BaseAddr,
        TREMULSW_HW_INFO_REG,
        TREMULSW_CTRL_REG,
        TREMULSW_INFO_REG,
        TREMULSW_STAT_REG,
        TREMULSW_RESET,
        0,
        0,
        TREMULSW_RX_DATA_CTR_REG,
        TREMULSW_RX_SYMBOL_CTR_REG,
        TREMULSW_RX_FRAME_CTR_REG,
        TREMULSW_RX_BLOCK_CTR_REG,
        TREMULSW_TX_DATA_CTR_REG,
        TREMULSW_TX_SYMBOL_CTR_REG,
        TREMULSW_TX_FRAME_CTR_REG,
        TREMULSW_TX_BLOCK_CTR_REG,
        0,
        0
        );
    pd->pDvcDesc  = pDvcDesc;
    pd->pDvcData  = pDvcData;
    DvcDataInit( pd );
    SYSRESET( pd );
    SetParams( pd );
}

//========================
// Dump Device to Console
//========================
static void Dump( TREmulSw *pd ) {
    cprintf( "\n" );
	SysDvceDump( pd );
	SysDescDump( pd );
	DvcDescDump( pd );
	SysDataDump( pd );
	DvcDataDump( pd );
	cprintf( "\n" );
}

// End of boilerplate code
//========================
// Device-specific code
//========================
static void SetParams( TREmulSw *pd )
{
	SYSWRSLVREG( pd,  TREMUL_N_FFT_REG,         SYSPAR( nFFT        ));
    SYSWRSLVREG( pd,  TREMUL_N_CPRE_REG,        SYSPAR( nCPre       ));
	SYSWRSLVREG( pd,  TREMUL_N_SYMB_FRAME_REG,  SYSPAR( nSymbAGC ) + SYSPAR( nSymbNSD) + SYSPAR( nSymbPRS ) + SYSPAR( nSymbData ));
	SYSWRSLVREG( pd,  TREMUL_N_FRAME_BLOCK_REG, SYSPAR( nFrameBlock ));
	SYSWRSLVREG( pd,  TREMULSW_CONFIG_REG,      DVCPAR( ConfigReg   ));
}

static void DvcDataInit( TREmulSw *pd )
{

}

static void DvcDataUpd( TREmulSw *pd )
{
    u32 statReg = SYSRDSLVREG( pd, TREMULSW_STAT_REG );
    DVCDAT( sel_emul                 ) = ( statReg & TREMULSW_MSK_SEL_EMUL ) >> TREMULSW_BIT_SEL_EMUL;
    DVCDAT( fifo_clk_lock            ) = ( statReg & TREMULSW_FIFO_CLK_LOCK             ? true : false );
    DVCDAT( s_axis_in_tvalid         ) = ( statReg & TREMULSW_S_AXIS_IN_TVALID          ? true : false );
    DVCDAT( s_axis_in_tready         ) = ( statReg & TREMULSW_S_AXIS_IN_TREADY          ? true : false );
    DVCDAT( m_axis_out_tvalid        ) = ( statReg & TREMULSW_M_AXIS_OUT_TVALID         ? true : false );
    DVCDAT( m_axis_out_tready        ) = ( statReg & TREMULSW_M_AXIS_OUT_TREADY         ? true : false );
    DVCDAT( s_axis_device_out_tvalid ) = ( statReg & TREMULSW_S_AXIS_DEVICE_OUT_TVALID  ? true : false );
    DVCDAT( s_axis_device_out_tready ) = ( statReg & TREMULSW_S_AXIS_DEVICE_OUT_TREADY  ? true : false );
    DVCDAT( m_axis_device_in_tvalid  ) = ( statReg & TREMULSW_M_AXIS_DEVICE_IN_TVALID   ? true : false );
    DVCDAT( m_axis_device_in_tready  ) = ( statReg & TREMULSW_M_AXIS_DEVICE_IN_TREADY   ? true : false );
    DVCDAT( s_axis_emul_out_tvalid   ) = ( statReg & TREMULSW_S_AXIS_EMUL_OUT_TVALID    ? true : false );
    DVCDAT( s_axis_emul_out_tready   ) = ( statReg & TREMULSW_S_AXIS_EMUL_OUT_TREADY    ? true : false );
    DVCDAT( m_axis_emul_in_tvalid    ) = ( statReg & TREMULSW_M_AXIS_EMUL_IN_TVALID     ? true : false );
    DVCDAT( m_axis_emul_in_tready    ) = ( statReg & TREMULSW_M_AXIS_EMUL_IN_TREADY     ? true : false );
}

static void DvcDataDump( TREmulSw *pd )
{
	DvcDataUpd( pd );
    cprintf( "\nfifo_clk_lock              : %s", TRUEFALSE( DVCDAT( fifo_clk_lock      ) ) );
	cprintf( "\nsel_emul                   : %u", toUint( DVCDAT( sel_emul ) ) );
    if ( DVCDAT( sel_emul ))
    {
        cprintf( "\ns_axis_device_out_tvalid   : %s", OFFON( DVCDAT( s_axis_device_out_tvalid ) ) );
        cprintf( "\ns_axis_device_out_tready   : %s", OFFON( DVCDAT( s_axis_device_out_tready ) ) );
        cprintf( "\nm_axis_device_in_tvalid    : %s", OFFON( DVCDAT( m_axis_device_in_tvalid  ) ) );
        cprintf( "\nm_axis_device_in_tready    : %s", OFFON( DVCDAT( m_axis_device_in_tready  ) ) );
        cprintf( "\n" );
        cprintf( "\ns_axis_in_tvalid           : %s", OFFON( DVCDAT( s_axis_in_tvalid         ) ) );
        cprintf( "\nm_axis_emul_in_tvalid      : %s", OFFON( DVCDAT( m_axis_emul_in_tvalid    ) ) );
        cprintf( "\ns_axis_in_tready           : %s", OFFON( DVCDAT( s_axis_in_tready         ) ) );
        cprintf( "\nm_axis_emul_in_tready      : %s", OFFON( DVCDAT( m_axis_emul_in_tready    ) ) );
        cprintf( "\nm_axis_out_tvalid          : %s", OFFON( DVCDAT( m_axis_out_tvalid        ) ) );
        cprintf( "\ns_axis_emul_out_tvalid     : %s", OFFON( DVCDAT( s_axis_emul_out_tvalid   ) ) );
        cprintf( "\nm_axis_out_tready          : %s", OFFON( DVCDAT( m_axis_out_tready        ) ) );
        cprintf( "\ns_axis_emul_out_tready     : %s", OFFON( DVCDAT( s_axis_emul_out_tready   ) ) );
    }
    else
    {
        cprintf( "\ns_axis_emul_out_tvalid     : %s", OFFON( DVCDAT( s_axis_emul_out_tvalid   ) ) );
        cprintf( "\ns_axis_emul_out_tready     : %s", OFFON( DVCDAT( s_axis_emul_out_tready   ) ) );
        cprintf( "\nm_axis_emul_in_tvalid      : %s", OFFON( DVCDAT( m_axis_emul_in_tvalid    ) ) );
        cprintf( "\nm_axis_emul_in_tready      : %s", OFFON( DVCDAT( m_axis_emul_in_tready    ) ) );
        cprintf( "\n" );
        cprintf( "\ns_axis_in_tvalid           : %s", OFFON( DVCDAT( s_axis_in_tvalid         ) ) );
        cprintf( "\nm_axis_device_in_tvalid    : %s", OFFON( DVCDAT( m_axis_device_in_tvalid  ) ) );
        cprintf( "\ns_axis_in_tready           : %s", OFFON( DVCDAT( s_axis_in_tready         ) ) );
        cprintf( "\nm_axis_device_in_tready    : %s", OFFON( DVCDAT( m_axis_device_in_tready  ) ) );
        cprintf( "\nm_axis_out_tvalid          : %s", OFFON( DVCDAT( m_axis_out_tvalid        ) ) );
        cprintf( "\ns_axis_device_out_tvalid   : %s", OFFON( DVCDAT( s_axis_device_out_tvalid ) ) );
        cprintf( "\nm_axis_out_tready          : %s", OFFON( DVCDAT( m_axis_out_tready        ) ) );
        cprintf( "\ns_axis_device_out_tready   : %s", OFFON( DVCDAT( s_axis_device_out_tready ) ) );
    }
}

static const char s_SEL_RF             [ ] = "RF Channel";     // 0
static const char s_SEL_EMUL           [ ] = "Emulator";       // 1
static const char s_SEL_WIRE2          [ ] = "Wire";           // 2
static const char s_SEL_WIRE3          [ ] = "Wire";           // 2
static const char *SelEmulStr [ 4 ] = {
        s_SEL_RF,
        s_SEL_EMUL,
        s_SEL_WIRE2,
        s_SEL_WIRE3
};


static void DvcDescDump( TREmulSw *pd )
{
    cprintf( "Selected device  : %s\n",        SelEmulStr[ ( DVCPAR( ConfigReg ) & TREMULSW_CFG_MSK_SEL_EMUL ) >> TREMULSW_CFG_BIT_SEL_EMUL ]);
    cprintf( "\n" );
}

static void SysDvceDump( TREmulSw *pd )
{
	SYSDVCEDUMP( pd );
}

static void SysDataDump( TREmulSw *pd )
{
	SYSDATADUMP( pd );
}

static void SysDescDump( TREmulSw *pd )
{
	cprintf("\nnFFT             : %u",        toUint( SYSPAR( nFFT        )));
    cprintf("\nnCPre            : %u",        toUint( SYSPAR( nCPre       )));
	cprintf("\nnSymbFrame       : %u",        toUint( SYSPAR( nSymbAGC ) + SYSPAR( nSymbNSD ) + SYSPAR( nSymbPRS ) + SYSPAR( nSymbData )));
	cprintf("\nnFrameBlock      : %u",        toUint( SYSPAR( nFrameBlock )));
    cprintf( "\n" );
}

