/*
 * TxRx_Emulator.c
 *
 *  Created on: 1 May 2018
 *      Author: JULIAN MORTIMER
 */
#include "include.h"


// =================================
// Boilerplate Function Declarations
// =================================
void TREmulInit(
    TREmul *pd,
    const char      *Name,
    u32             BaseAddr,
    SysCtrl         *pSysCtrl,
    SysDvce         *pSysDvce,
    SysDesc         *pSysDesc,
    SysPtrs         *pSysPtrs,
    SysData         *pSysData,
    TREmulDesc *pDvcDesc,
    TREmulData *pDvcData
    );
static void SetParams    ( TREmul *pd );
static void SysDvceDump  ( TREmul *pd );
static void SysDescDump  ( TREmul *pd );
static void SysDataDump  ( TREmul *pd );
static void DvcDataInit  ( TREmul *pd );
static void DvcDataUpd   ( TREmul *pd );
static void DvcDescDump  ( TREmul *pd );
static void DvcDataDump  ( TREmul *pd );

// =================================================
// Device Structure Member Functions (public access)
// =================================================
static inline void Start             ( TREmul *pd ) { hwWriteCtrlReg( SYSDVC( BaseAddr ), TREMUL_CTRL_REG << 2, TREMUL_START   ); }
static inline void Stop              ( TREmul *pd ) { hwWriteCtrlReg( SYSDVC( BaseAddr ), TREMUL_CTRL_REG << 2, TREMUL_STOP    ); }
static u32         TestBram          ( TREmul *pd );
static u32         ReadPsLUT         ( TREmul *pd, u32 addr );
static void        WritePsLUT        ( TREmul *pd, u32 addr, u32 data );
static void        InitPsLUT         ( TREmul *pd, u32 addr, u32 count, u32 data );
static void        LoadPsLUT         ( TREmul *pd );
static void        DumpPsLUT         ( TREmul *pd, u32 addr, u32 nBits );
static u32         ReadHwLUT         ( TREmul *pd, u32 addr );
static void        WriteHwLUT        ( TREmul *pd, u32 addr, u32 data );
static void        InitHwLUT         ( TREmul *pd, u32 addr, u32 count, u32 data );
static void        LoadHwLUT         ( TREmul *pd );
static void        DumpHwLUT         ( TREmul *pd, u32 addr, u32 nBits );
static void        Dump              ( TREmul *pd );


//static _Bool       SelfTest          ( _Bool verbose );

// Declare init structure here with member functions
TREmul TREmulInitStr = {
	NULL,
	NULL,
	NULL,
	Start,
	Stop,
	TestBram,
	ReadPsLUT,
	WritePsLUT,
	InitPsLUT,
	LoadPsLUT,
	DumpPsLUT,
	ReadHwLUT,
	WriteHwLUT,
	InitHwLUT,
	LoadHwLUT,
	DumpHwLUT,
	Dump
};


//===================
// Initialize Device
//===================
void TREmulInit( TREmul *pd, const char *Name, u32 BaseAddr, SysCtrl *pSysCtrl, SysDvce *pSysDvce, SysDesc *pSysDesc, SysPtrs *pSysPtrs, SysData *pSysData, TREmulDesc *pDvcDesc, TREmulData *pDvcData )
{
    *pd             = TREmulInitStr;
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
        TREMUL_HW_INFO_REG,
        TREMUL_CTRL_REG,
        TREMUL_INFO_REG,
        TREMUL_STAT_REG,
        TREMUL_RESET,
        TREMUL_INT_ACK,
        0,
        TREMUL_RX_DATA_CTR_REG,
        TREMUL_RX_SYMBOL_CTR_REG,
        TREMUL_RX_FRAME_CTR_REG,
        TREMUL_RX_BLOCK_CTR_REG,
        TREMUL_TX_DATA_CTR_REG,
        TREMUL_TX_SYMBOL_CTR_REG,
        TREMUL_TX_FRAME_CTR_REG,
        TREMUL_TX_BLOCK_CTR_REG,
        0,
        0
        );
    pd->pDvcDesc  = pDvcDesc;
    pd->pDvcData  = pDvcData;
    DvcDataInit( pd );
    SYSRESET( pd );
  	( void ) pd->TestBram( pd );
    pd->LoadPsLUT( pd );
    pd->LoadHwLUT( pd );
    SetParams( pd );
    //    return SelfTest( pde, true );
}

//========================
// Dump Device to Console
//========================
static void Dump( TREmul *pd ) {
    cprintf("\n");
	SysDvceDump( pd );
	SysDescDump( pd );
	DvcDescDump( pd );
	SysDataDump( pd );
	DvcDataDump( pd );
	cprintf("\n");
}
// End of boilerplate code
//========================
// Device-specific code
//========================
static void SetParams( TREmul *pd )
{
	SYSWRSLVREG( pd,  TREMUL_N_FFT_REG,         SYSPAR( nFFT ));
    SYSWRSLVREG( pd,  TREMUL_N_CPRE_REG,        SYSPAR( nCPre ));
	SYSWRSLVREG( pd,  TREMUL_N_SYMB_FRAME_REG,  SYSPAR( nSymbAGC ) + SYSPAR( nSymbNSD ) + SYSPAR( nSymbPRS ) + SYSPAR( nSymbData ));
	SYSWRSLVREG( pd,  TREMUL_N_FRAME_BLOCK_REG, SYSPAR( nFrameBlock ));
}

static void DvcDataInit( TREmul *pd )
{
// Initialise DvcData members here
// Example: DVCDAT( Bob ) = 0;
}

static void DvcDataUpd( TREmul *pd )
{
// Update DvcData members here
// Example:	DVCDAT( Bob ) = SYSRDSLVREG( pd, TREMUL_BOB_REG );
}

static void DvcDataDump( TREmul *pd )
{
	DvcDataUpd( pd );
// Insert device-specific SysDesc members here
// Example: cprintf("\nBob  : %u", toUint( DVCDAT( Bob )) );
}

static void DvcDescDump( TREmul *pd )
{
// Insert device-specific TREmulDesc members here
// Example: cprintf("\nAlice: %u", toUint( DVCPAR( Alice )) );
}

static void SysDvceDump( TREmul *pd )
{
	SYSDVCEDUMP( pd );
}

static void SysDataDump( TREmul *pd )
{
	SYSDATADUMP( pd );
}

static void SysDescDump( TREmul *pd )
{
    cprintf( "\nnFFT             : %u",   toUint( SYSPAR( nFFT         ) ));
    cprintf( "\nnCPre            : %u",   toUint( SYSPAR( nCPre        ) ));
    cprintf( "\nnSymbframe       : %u",   toUint( SYSPAR( nSymbAGC ) + SYSPAR( nSymbNSD ) + SYSPAR( nSymbPRS ) + SYSPAR( nSymbData )));
    cprintf( "\nnframeblock      : %u",   toUint( SYSPAR( nFrameBlock  ) ));
    cprintf( "\n" );
}


// Insert device-specific code here
//======================================
// Device-specific non member functions
//======================================

static inline void EnterAdmin        ( TREmul *pd ) { hwWriteCtrlReg( SYSDVC( BaseAddr ), TREMUL_CTRL_REG << 2, TREMUL_ENTER_ADMIN ); }
static inline void ExitAdmin         ( TREmul *pd ) { hwWriteCtrlReg( SYSDVC( BaseAddr ), TREMUL_CTRL_REG << 2, TREMUL_EXIT_ADMIN  ); }
static void DumpLUT( TREmul *pd, u32 addr, u32 nbits, dmpSelector dsel );

//=================================
// Device-specific member functions
//=================================

static void DumpLUT( TREmul *pd, u32 addr, u32 nbits, dmpSelector dsel )
{
	char s[ 2 * SYS_MAX_OBJ_NAME_LEN + 1 ]; strcpy( s, SYSDVC( Name ) );
	switch( dsel ) {
	  case DUMP_SEL_TRE_HW_LUT :
		  strcat( s, " Emulator BRAM Buffer:" );
		  break;
	  case DUMP_SEL_TRE_PS_LUT :
	  default:
		  strcat( s, " Local PC Buffer" );
		  break;
	};

	if ( nbits == 0 )
		utilScreenDumpData( dsel, DUMP_FMT_HEX,       DUMP_TYP_MEMORY, addr, 4,     0, ( u32 )pd, s );
	else if(( nbits >= 8 ) && ( nbits <= 16 ))
     	utilScreenDumpData( dsel, DUMP_FMT_IQ_SIGNED, DUMP_TYP_MEMORY, addr, nbits, 0, ( u32 )pd, s );
	else
		cprintf( "Format must be 0 (hex) or 8-16 (decimal Q:I)\n" );
}


static u32 ReadPsLUT( TREmul *pd, u32 addr )
{
    return( DVCPAR( pLUT[ addr ] ));
}

static void WritePsLUT( TREmul *pd, u32 addr, u32 data )
{
    DVCPAR( pLUT[ addr ] ) = data;
}

static void InitPsLUT( TREmul *pd, u32 addr, u32 count, u32 val )
{
	for( u32 i = 0; i < count; i++ )
	{
       DVCPAR( pLUT[ i + addr ] ) = val;
	}
}

static void LoadPsLUT( TREmul *pd )
{
	u32 EmulBitMask = DVCPAR( PsLUTSize ) - 1;
	for( u32 i = 0; i < DVCPAR( PsLUTSize ); i++ )
	{
	   u32 u = i & EmulBitMask;
       DVCPAR( pLUT[ i ] ) = u + ( u << 12 );
	}
}

static void DumpPsLUT( TREmul *pd, u32 addr, u32 nbits )
{
	DumpLUT( pd, addr, nbits, DUMP_SEL_TRE_PS_LUT );
}

static u32 ReadHwLUT( TREmul *pd, u32 addr )
{
    u32 u = 0;
    EnterAdmin( pd );
    u = hwReadData( SYSDVC( BaseAddr ),
        TREMUL_CTRL_REG << 2, TREMUL_AD_RD_LUT, TREMUL_AD_RD_LUT_ACK,
        TREMUL_STAT_REG << 2, TREMUL_AD_LUT_RD_RDY,
        TREMUL_AD_RD_WR_ADDR_REG << 2, addr,
        TREMUL_AD_RD_DATA_REG << 2 );
    ExitAdmin( pd );
    return( u );
}

static void WriteHwLUT( TREmul *pd, u32 addr, u32 data )
{
	EnterAdmin( pd );
    hwWriteData( SYSDVC( BaseAddr ),
    	TREMUL_CTRL_REG << 2, TREMUL_AD_WR_LUT,
    	TREMUL_AD_RD_WR_ADDR_REG << 2, addr,
    	TREMUL_AD_WR_DATA_REG << 2, data );
    ExitAdmin( pd );
}

static void InitHwLUT( TREmul *pd, u32 addr, u32 count, u32 val )
{
    EnterAdmin( pd );
	for( u32 i = 0; i < count; i++ )
	{
	    hwWriteData( SYSDVC( BaseAddr ),
	        TREMUL_CTRL_REG << 2, TREMUL_AD_WR_LUT,
	        TREMUL_AD_RD_WR_ADDR_REG << 2, addr + i,
	        TREMUL_AD_WR_DATA_REG << 2, val );
	}
    ExitAdmin( pd );
}

static void LoadHwLUT( TREmul *pd )
{
    EnterAdmin( pd );
	for( u32 i = 0; i < SYS_MAX_CARRIERS; i++ )
	{
	    hwWriteData( SYSDVC( BaseAddr ),
	        TREMUL_CTRL_REG << 2, TREMUL_AD_WR_LUT,
	        TREMUL_AD_RD_WR_ADDR_REG << 2, i,
	        TREMUL_AD_WR_DATA_REG << 2, DVCPAR( pLUT[ i ] ));
	}
    ExitAdmin( pd );
}

static void DumpHwLUT( TREmul *pd, u32 addr, u32 nbits )
{
	DumpLUT( pd, addr, nbits, DUMP_SEL_TRE_HW_LUT );
}


static u32 ReadLUT( u32 addr )
{
	static TREmul *pd = NULL;
	if ( addr & 0x80000000 ) { pd = ( TREmul * ) ( addr & 0x7FFFFFFF ); }
	else if ( pd != NULL ) { return ( pd->ReadHwLUT( pd, ( u32 )addr )); }
	return( 0 );
}

static void WriteLUT( u32 addr, u32 data )
{
	static TREmul *pd = NULL;
	if ( addr & 0x80000000 ) { pd = ( TREmul * )( addr & 0x7FFFFFFF ); }
	else if ( pd != NULL ) { pd->WriteHwLUT( pd, ( u32 )addr, data ); }
}

static u32 TestLUT( TREmul *pd )
{
	char Name[ 2 * SYS_MAX_OBJ_NAME_LEN + 1]; strcpy( Name, SYSDVC( Name ) );
	return ( UtilTestBram( ReadLUT, WriteLUT, DVCPAR( HwLUTSize ), 0x00000FFF, strcat( Name, " emulator LUT" ) ));
}

static u32 TestBram( TREmul *pd )
{
    ( void ) ReadLUT( ( u32 ) pd | 0x80000000 );
    WriteLUT( ( u32 ) pd | 0x80000000, 0 );
    ( void ) TestLUT( pd );
	return( 0 );
}


