/*
 * PRS_CIR.c
 *
 *  Created on: 12 Sep 2018
 *      Author: JULIAN MORTIMER
 */
#include "include.h"

// Simulated CIR
static u32  nPass = 0;

static cplx TxPRSFd           [ SYS_MAX_CARRIERS ];
static cplx RxPRSFd           [ SYS_MAX_CARRIERS ];
static cplx RxPRSTd           [ SYS_MAX_CARRIERS ];
static cplx RxPRSTdCorDeltaT  [ SYS_MAX_CARRIERS ];
static cplx RxPRSTdCordeltafp [ SYS_MAX_CARRIERS ];
static cplx RxPRSFdCordeltafp [ SYS_MAX_CARRIERS ];

cplx RxPRSFdCorFT      [ SYS_MAX_CARRIERS ];
cplx RxPRSTdCorFT      [ SYS_MAX_CARRIERS ];

static s32  TimeCorrSampleIdx [ SYS_MAX_CARRIERS ];

static cplx FFT00 [ SYS_MAX_CARRIERS ];
static cplx FFT01 [ SYS_MAX_CARRIERS ];
static cplx FFT02 [ SYS_MAX_CARRIERS ];
static cplx FFT03 [ SYS_MAX_CARRIERS ];
static cplx FFT04 [ SYS_MAX_CARRIERS ];
static cplx FFT05 [ SYS_MAX_CARRIERS ];
static cplx FFT06 [ SYS_MAX_CARRIERS ];
static cplx FFT07 [ SYS_MAX_CARRIERS ];
static cplx FFT08 [ SYS_MAX_CARRIERS ];
static cplx FFT09 [ SYS_MAX_CARRIERS ];
static cplx FFT10 [ SYS_MAX_CARRIERS ];
static cplx FFT11 [ SYS_MAX_CARRIERS ];
static cplx FFT12 [ SYS_MAX_CARRIERS ];
static cplx FFT13 [ SYS_MAX_CARRIERS ];
static cplx FFT14 [ SYS_MAX_CARRIERS ];
static cplx FFT15 [ SYS_MAX_CARRIERS ];
static cplx FFT16 [ SYS_MAX_CARRIERS ];
static cplx FFT17 [ SYS_MAX_CARRIERS ];
static cplx FFT18 [ SYS_MAX_CARRIERS ];
static cplx FFT19 [ SYS_MAX_CARRIERS ];
static cplx FFT20 [ SYS_MAX_CARRIERS ];

static cplx *FFTPtrs [ SYS_N_CHAN_PRSCIR ] = {
    FFT00,
    FFT01,
    FFT02,
    FFT03,
    FFT04,
    FFT05,
    FFT06,
    FFT07,
    FFT08,
    FFT09,
    FFT10,
    FFT11,
    FFT12,
    FFT13,
    FFT14,
    FFT15,
    FFT16,
    FFT17,
    FFT18,
    FFT19,
    FFT20
};

static double Abs00 [ SYS_MAX_CARRIERS ];
static double Abs01 [ SYS_MAX_CARRIERS ];
static double Abs02 [ SYS_MAX_CARRIERS ];
static double Abs03 [ SYS_MAX_CARRIERS ];
static double Abs04 [ SYS_MAX_CARRIERS ];
static double Abs05 [ SYS_MAX_CARRIERS ];
static double Abs06 [ SYS_MAX_CARRIERS ];
static double Abs07 [ SYS_MAX_CARRIERS ];
static double Abs08 [ SYS_MAX_CARRIERS ];
static double Abs09 [ SYS_MAX_CARRIERS ];
static double Abs10 [ SYS_MAX_CARRIERS ];
static double Abs11 [ SYS_MAX_CARRIERS ];
static double Abs12 [ SYS_MAX_CARRIERS ];
static double Abs13 [ SYS_MAX_CARRIERS ];
static double Abs14 [ SYS_MAX_CARRIERS ];
static double Abs15 [ SYS_MAX_CARRIERS ];
static double Abs16 [ SYS_MAX_CARRIERS ];
static double Abs17 [ SYS_MAX_CARRIERS ];
static double Abs18 [ SYS_MAX_CARRIERS ];
static double Abs19 [ SYS_MAX_CARRIERS ];
static double Abs20 [ SYS_MAX_CARRIERS ];

static double *AbsPtrs [ SYS_N_CHAN_PRSCIR ] = {
    Abs00,
    Abs01,
    Abs02,
    Abs03,
    Abs04,
    Abs05,
    Abs06,
    Abs07,
    Abs08,
    Abs09,
    Abs10,
    Abs11,
    Abs12,
    Abs13,
    Abs14,
    Abs15,
    Abs16,
    Abs17,
    Abs18,
    Abs19,
    Abs20
};

static u32 ColIdx;
static u32 RowIdx;

static double DeltaTAcc = 0.0;
static double DeltaFAcc = 0.0;


// =================================
// Boilerplate Function Declarations
// =================================
void PrsCIRInit(
    PrsCIR *pd,
    const char *Name,
    u32        BaseAddr,
    SysCtrl    *pSysCtrl,
    SysDvce    *pSysDvce,
    SysDesc    *pSysDesc,
    SysPtrs    *pSysPtrs,
    SysData    *pSysData,
    PrsDesc    *pPrsCIRDesc,
    PrsCIRData *pPrsCIRData
    );


static void SetParams    ( PrsCIR *pd );
static void SysDvceDump  ( PrsCIR *pd );
static void SysDescDump  ( PrsCIR *pd );
static void SysDataDump  ( PrsCIR *pd );
static void DvcDataInit  ( PrsCIR *pd );
static void DvcDataUpd   ( PrsCIR *pd );
static void DvcDescDump  ( PrsCIR *pd );
static void DvcDataDump  ( PrsCIR *pd );



// =================================================
// Device Structure Member Functions (public access)
// =================================================
// Insert device-specific member function prototypes here
// Example: static inline void Start         ( PrsCIR *pd ) { hwWriteCtrlReg( SYSDVC( BaseAddr ), PRSCIR_CTRL_REG << 2, PRSCIR_START ); }
// Example: static void        WriteBuf      ( PrsCIR *pd, u32 addr, u32 data );
static u32  TestBram                 ( PrsCIR *pd );
static u64  ReadBram                 ( PrsCIR *pd, u32 addr );
static void WriteBram                ( PrsCIR *pd, u32 addr, u64 data );
static u32  ReadBramFFT              ( PrsCIR *pd, u32 addr );
static void WriteBramFFT             ( PrsCIR *pd, u32 addr, u32 data );
static void InitBramFFT              ( PrsCIR *pd, u32 addr, u32 count, u32 data );
static void DumpBramFFT              ( PrsCIR *pd, u32 addr, u32 nBits );
static u32  ReadBramABS              ( PrsCIR *pd, u32 addr );
static void WriteBramABS             ( PrsCIR *pd, u32 addr, u32 data );
static void InitBramABS              ( PrsCIR *pd, u32 addr, u32 count, u32 data );
static void DumpBramABS              ( PrsCIR *pd, u32 addr, u32 nBits );
static void Dump                     ( PrsCIR *pd );
static void CalcDeltas               ( PrsCIR *pd );
static void DumpCIRTable             ( PrsCIR *pd );
static void DumpRxPRSTd              ( PrsCIR *pd );
static void DumpRxPRSFd              ( PrsCIR *pd );
static void DumpTxPRSFd              ( PrsCIR *pd );
static void DumpRxPRSTdCorDeltaT     ( PrsCIR *pd );
static void DumpRxPRSTdCordeltafp    ( PrsCIR *pd );
static void DumpRxPRSFdCordeltafp    ( PrsCIR *pd );
static void DumpRxPRSFdCorFT         ( PrsCIR *pd );
static void DumpRxPRSTdCorFT         ( PrsCIR *pd );




// =================================================
// Device Non Member Functions
// =================================================
static void SimReset   ( void );
static u32  SimHwReg   ( u32 baseAddr, u32 reg, u32 data, _Bool bWrite );
static u32  SimHwRdReg ( u32 reg );
static void SimHwWrReg ( u32 reg, u32 data );


static void ReadPRSData  ( PrsCIR *pd );
static void CopyPRS      ( PrsCIR *pd );
static void CalcFFT      ( PrsCIR *pd );
static void CalcAbs      ( PrsCIR *pd );
static void CalcDelta    ( PrsCIR *pd );
static void Calc         ( PrsCIR *pd );
static void Recalc       ( PrsCIR *pd );
static void CorrectRxPRS ( PrsCIR *pd );


// Declare init structure here with member functions
PrsCIR PrsCIRInitStr = {
    NULL,
    NULL,
    NULL,
    TestBram              ,
    ReadBram              ,
    WriteBram             ,
    ReadBramFFT           ,
    WriteBramFFT          ,
    InitBramFFT           ,
    DumpBramFFT           ,
    ReadBramABS           ,
    WriteBramABS          ,
    InitBramABS           ,
    DumpBramABS           ,
    Dump                  ,
    Calc                  ,
    Recalc                ,
    DumpCIRTable          ,
    DumpRxPRSTd           ,
    DumpRxPRSFd           ,
    DumpTxPRSFd           ,
    DumpRxPRSTdCorDeltaT  ,
    DumpRxPRSTdCordeltafp ,
    DumpRxPRSFdCordeltafp ,
    DumpRxPRSFdCorFT      ,
    DumpRxPRSTdCorFT
};

//===================
// Initialize Device
//===================
void PrsCIRInit( PrsCIR *pd, const char *Name, u32 BaseAddr, SysCtrl *pSysCtrl, SysDvce *pSysDvce, SysDesc *pSysDesc, SysPtrs *pSysPtrs, SysData *pSysData, PrsDesc *pDvcDesc, PrsCIRData *pDvcData )
{
    *pd             = PrsCIRInitStr;
    pd->pSysDvce    = pSysDvce;
    SysDvceInit(
        pd->pSysDvce,
        Name,
        pSysCtrl,
        pSysDesc,
        pSysPtrs,
        pSysData,
        SimHwReg,
        BaseAddr,
        CIR_HW_INFO_REG,
        CIR_CTRL_REG,
        CIR_INFO_REG,
        CIR_STAT_REG,
        CIR_RESET,
        CIR_INT_ACK,
        0,
        CIR_PRS_RX_DATA_CTR_REG,
        CIR_PRS_RX_SYMB_CTR_REG,
        CIR_PRS_RX_FRAME_CTR_REG,
        CIR_PRS_RX_BLOCK_CTR_REG,
        CIR_PRS_TX_DATA_CTR_REG,
        CIR_PRS_TX_SYMB_CTR_REG,
        CIR_PRS_TX_FRAME_CTR_REG,
        CIR_PRS_TX_BLOCK_CTR_REG,
        0,
        0
        );
    pd->pDvcDesc  = pDvcDesc;
    pd->pDvcData  = pDvcData;
    DvcDataInit( pd );
    SimReset();
    ( void ) pd->TestBram( pd );    SYSRESET( pd );
    SetParams( pd );
}

//========================
// Dump Device to Console
//========================
static void Dump( PrsCIR *pd ) {
    cprintf( "\n" );
    SysDvceDump( pd );
    SysDescDump( pd );
    DvcDescDump( pd );
//    SysDataDump( pd );
//    DvcDataDump( pd );
    cprintf( "\n" );
}

// End of boilerplate code
//========================
// Device-specific code
//========================
static void SetParams( PrsCIR *pd )
{
// Set PL device slave regs here
// Example SYSWRSLVREG( pd, PRSCIR_BOB_REG,    SYSPAR( Bob   ));
// Example SYSWRSLVREG( pd, PRSCIR_ALICE_REG , SYSPAR( Alice ));
}

static void DvcDataInit( PrsCIR *pd )
{

}

static void DvcDataUpd( PrsCIR *pd )
{
// Update DvcData members here
// Example: DVCDAT( Bob ) = SYSRDSLVREG( pd, PRSCIR_BOB_REG );
}

static void DvcDataDump( PrsCIR *pd )
{
    DvcDataUpd( pd );
// Insert device-specific SysDesc members here
// Example: cprintf("Bob  : %u", toUint( DVCDAT( Bob )) );   fflush( stdout );
}

static void DvcDescDump( PrsCIR *pd )
{
    cprintf("DeltaT           : %i\n", toInt( DVCPAR( DeltaT )));
    cprintf("DeltaF           : %i\n", toInt( DVCPAR( DeltaF )));
    cprintf("deltat'          : %.6g\n", DVCPAR( deltatp ));
    cprintf("deltaf'          : %.6g\n", DVCPAR( deltafp ));
    cprintf("Combined delta f : %.2f kHz\n", fabs( ( float ) (( DVCPAR( DeltaF ) + DVCPAR( deltafp )) * AD9364RxSampleFreq() / SYSPAR( nFFT ) ) / 1000.00 ) );
    cprintf("Passes           : %u\n", toUint( nPass ));
    cprintf("Accumulated delta T : %.6f sams\n", DeltaTAcc );
    cprintf("Accumulated delta F : %.6f kHz\n", DeltaFAcc * ( AD9364RxSampleFreq() / SYSPAR( nFFT ) ) / 1000.00 );
    cprintf("\n");
}

static void SysDvceDump( PrsCIR *pd )
{
    SYSDVCEDUMP( pd );
}

static void SysDataDump( PrsCIR *pd )
{
    SYSDATADUMP( pd );
}

static void SysDescDump( PrsCIR *pd )
{
// Insert device-specific PrsCIRDesc members here
// Example: cprintf("Alice: %u", toUint( SYSPAR( Alice )) ); fflush( stdout );
}


// Insert device-specific code here
//======================================
// Device-specific non member functions
//======================================
static inline void EnterAdmin        ( PrsCIR *pd ) { hwWriteCtrlReg( SYSDVC( BaseAddr ), CIR_CTRL_REG << 2, CIR_ENTER_ADMIN ); }
static inline void ExitAdmin         ( PrsCIR *pd ) { hwWriteCtrlReg( SYSDVC( BaseAddr ), CIR_CTRL_REG << 2, CIR_EXIT_ADMIN  ); }
// Example: static inline void EnterAdmin    ( PrsCIR *pd ) { hwWriteCtrlReg( SYSDVC( BaseAddr ), PRSCIR_CTRL_REG << 2, PRSCIR_ENTER_ADMIN ); }


//=================================
// Device-specific member functions
//=================================

// Example:
// static u32 ReadBuf( PrsCIR *pd, u32 addr )
// {
//     return ( SYSPTR( pBuf[ addr ] ));
// }

// BRAM Access from PS
static void DumpCIRTable( PrsCIR *pd )
{
    u32 nFFT     = SYSPAR( nFFT );
    cprintf( "\n" );
    for( int i = 0; i < nFFT; i++ ) {
        cprintf( "%4u", toUint( i ) );
        for( int j = 0; j < SYS_N_CHAN_PRSCIR; j++ ) {
             double *pAbs = AbsPtrs [ j ];
             cprintf(" %.3e", pAbs[ i ] );
        }
        cprintf( "\n" );
        if ( ( ( i + 1 ) % N_LINES_DUMP ) == 0 ) {
            char c = 0;
            if (( c = GetChar() ) == CC_ESC ) {
                break;
            }
        }
    }
}


static void DumpRxPRSTd( PrsCIR *pd )
{
    DumpCplx( RxPRSTd, SYSPAR( nFFT ) );
}

static void DumpRxPRSTdCorDeltaT( PrsCIR *pd )
{
    DumpCplx( SYSPTR( pRxPRSTdCorDeltaT ), SYSPAR( nFFT ) );
}

static void DumpRxPRSTdCordeltafp( PrsCIR *pd )
{
    DumpCplx( RxPRSTdCordeltafp, SYSPAR( nFFT ) );
}

static void DumpRxPRSFdCordeltafp( PrsCIR *pd )
{
    DumpCplx( RxPRSFdCordeltafp, SYSPAR( nFFT ) );
}

static void DumpRxPRSFdCorFT( PrsCIR *pd )
{
    DumpCplx( SYSPTR( pRxPRSFdCorFT ), SYSPAR( nFFT ) );
}

static void DumpRxPRSFd( PrsCIR *pd )
{
    DumpCplx( RxPRSFd, SYSPAR( nFFT ) );
}

static void DumpTxPRSFd( PrsCIR *pd )
{
    DumpCplx( TxPRSFd, SYSPAR( nFFT ) );
}

static void DumpRxPRSTdCorFT( PrsCIR *pd )
{
    DumpCplx( SYSPTR( pRxPRSTdCorFT ), SYSPAR( nFFT ) );
}


static void DumpBramFFT( PrsCIR *pd, u32 addr, u32 nbits )
{
    char s[ 2 * SYS_MAX_OBJ_NAME_LEN + 1 ]; strcpy( s, SYSDVC( Name ) );
    strcat( s, " Bram FFT" );
    dmpSelector dsel = DUMP_SEL_CIR_BRAM_FFT;

    if ( nbits == 0 )
        utilScreenDumpData( dsel, DUMP_FMT_HEX,       DUMP_TYP_MEMORY, addr, 4,     0, ( u32 )pd, s );
    else if(( nbits >= 8 ) && ( nbits <= 16 ))
        utilScreenDumpData( dsel, DUMP_FMT_IQ_SIGNED, DUMP_TYP_MEMORY, addr, nbits, 0, ( u32 )pd, s );
    else
        cprintf( "Format must be 0 (hex) or 8-16 (decimal Q:I)\n" );
}

static void DumpBramABS( PrsCIR *pd, u32 addr, u32 nbits )
{
    char s[ 2 * SYS_MAX_OBJ_NAME_LEN + 1 ]; strcpy( s, SYSDVC( Name ) );
    strcat( s, " Bram ABS" );
    dmpSelector dsel = DUMP_SEL_CIR_BRAM_ABS;

    if ( nbits == 0 )
        utilScreenDumpData( dsel, DUMP_FMT_HEX,       DUMP_TYP_MEMORY, addr, 4,     0, ( u32 )pd, s );
    else if(( nbits >= 8 ) && ( nbits <= 16 ))
        utilScreenDumpData( dsel, DUMP_FMT_IQ_SIGNED, DUMP_TYP_MEMORY, addr, nbits, 0, ( u32 )pd, s );
    else
        cprintf( "Format must be 0 (hex) or 8-16 (decimal Q:I)\n" );
}

static u64 ReadBram( PrsCIR *pd, u32 addr )
{
    EnterAdmin( pd );
    u64 u = hwReadData64( SYSDVC( BaseAddr ),
                   CIR_CTRL_REG << 2, CIR_AD_RD_BRAM, CIR_AD_RD_ACK,
                   CIR_STAT_REG << 2, CIR_AD_BRAM_RD_RDY,
                   CIR_AD_RD_WR_ADDR_REG << 2, addr,
                   CIR_AD_RD_DATA_REG << 2
                   );
    ExitAdmin( pd );
    return( u );
}

static u32 ReadBramFFT( PrsCIR *pd, u32 addr )
{
    return( ReadBram( pd, addr ) & 0x00000000FFFFFFFF );
}

static u32 ReadBramABS( PrsCIR *pd, u32 addr )
{
    return( ( ReadBram( pd, addr ) >> 32 ) & 0x00000000FFFFFFFF );
}

static void WriteBram( PrsCIR *pd, u32 addr, u64 data )
{
    EnterAdmin( pd );
    hwWriteData64( SYSDVC( BaseAddr ),
                CIR_CTRL_REG << 2, CIR_AD_WR_BRAM,
                CIR_AD_RD_WR_ADDR_REG << 2, addr,
                CIR_AD_WR_DATA_REG << 2, data
                );
    ExitAdmin( pd );
}

static void WriteBramFFT( PrsCIR *pd, u32 addr, u32 data )
{
    u64 u = ReadBram( pd, addr );
    WriteBram( pd, addr, ( u & 0xFFFFFFFF00000000 ) | ( u64 ) data );
}

static void WriteBramABS( PrsCIR *pd, u32 addr, u32 data )
{
    u64 u = ReadBram( pd, addr );
    WriteBram( pd, addr, ( u & 0x00000000FFFFFFFF ) | ( ( u64 ) data << 32 ));
}


static void InitBramFFT( PrsCIR *pd, u32 addr, u32 count, u32 val )
{
    for( int i = 0; i < count; i++ )
    {
        WriteBramFFT( pd, addr + i, val );
    }
}

static void InitBramABS( PrsCIR *pd, u32 addr, u32 count, u32 val )
{
    for( int i = 0; i < count; i++ )
    {
        WriteBramABS( pd, addr + i, val );
    }
}



// BRAM Test Code
static u64 ReadBR( u32 addr )
{
    static PrsCIR *pd = NULL;
    if ( addr & 0x80000000 ) { pd = ( PrsCIR * ) ( addr & 0x7FFFFFFF ); }
    else if ( pd != NULL ) { return ( pd->ReadBram( pd, ( u32 )addr )); }
    return( 0 );
}

static void WriteBR( u32 addr, u64 data )
{
    static PrsCIR *pd = NULL;
    if ( addr & 0x80000000 ) { pd = ( PrsCIR * )( addr & 0x7FFFFFFF ); }
    else if ( pd != NULL ) { pd->WriteBram( pd, ( u32 )addr, data ); }
}


static u32 TestBram( PrsCIR *pd )
{
    char Name[ 2 * SYS_MAX_OBJ_NAME_LEN + 1]; strcpy( Name, SYSDVC( Name ) );
    ( void ) ReadBR( ( u32 ) pd | 0x80000000 );
    WriteBR( ( u32 ) pd | 0x80000000, 0 );
    return ( UtilTestBram64( ReadBR, WriteBR, DVCPAR( CIRBramSize ), 0x0000FFFFFFFFFFFF, Name ));
}


// Simulator
static u32 SlvRegs[ CIR_N_SLV_REGS ];
static u64 BRAM[ PRS_CIR_BRAM_SIZE ];

static void SimReset( void )
{
    memset( SlvRegs, 0, CIR_N_SLV_REGS * sizeof( u32 ));
    memset( BRAM, 0, PRS_CIR_BRAM_SIZE * sizeof( u64 ));
    SlvRegs[ CIR_INFO_REG        ] = CIR_HW_INFO_REG;
    SlvRegs[ CIR_STATE_REG       ] = CIR_STATE_REG_IDLE;
    SlvRegs[ CIR_AM_STATE1_REG   ] = CIR_AM_STATE1_REG_IDLE;
    SlvRegs[ CIR_AM_STATE2_REG   ] = CIR_AM_STATE2_REG_IDLE;
    SlvRegs[ CIR_AM_STATE3_REG   ] = CIR_AM_STATE3_REG_IDLE;
    SlvRegs[ CIR_DO_STATE1_REG   ] = CIR_DO_STATE1_REG_IDLE;
    SlvRegs[ CIR_DO_STATE2_REG   ] = CIR_DO_STATE2_REG_IDLE;
    SlvRegs[ CIR_DO_STATE3_REG   ] = CIR_DO_STATE3_REG_IDLE;
    SlvRegs[ CIR_STRM_STATE1_REG ] = CIR_STRM_STATE1_IDLE;
    SlvRegs[ CIR_STRM_STATE2_REG ] = CIR_STRM_STATE2_IDLE;
}


static u32 SimHwReg ( u32 baseAddr, u32 reg, u32 data, _Bool bWrite )
{
    u32 rv = 0;
    if ( bWrite )
    {
        SimHwWrReg( reg, data );
        rv = 1;
    }
    else
    {
        rv = SimHwRdReg( reg );
    }
    return ( rv );
}

static u32 SimHwRdReg( u32 reg )
{
    u32 rv = 0;
    switch( reg )
    {
    case CIR_CTRL_REG:
        break;
    case CIR_N_CPRE_REG :
    case CIR_N_FFT_REG :
    case CIR_N_SYMB_FRAME_REG :
    case CIR_N_FRAME_BLOCK_REG :
    case CIR_AD_RD_WR_ADDR_REG :
    case CIR_AD_WR_DATA_REG :
    case CIR_AD_WR_DATA_REG_H :
    case CIR_LOG_N_FFT_REG :
    case CIR_KFF_REG :
    case CIR_KTF_REG :
    case CIR_KTT_REG :
    case CIR_KFT_REG :
    case CIR_BRAM_CHAN_SEL_REG :
    case CIR_INFO_REG :
    case CIR_STAT_REG :
    case CIR_AD_RD_DATA_REG :
    case CIR_AD_RD_DATA_REG_H :
    case CIR_HNDSHK_REG :
    case CIR_AXIS_CIR_TVALID_REG :
    case CIR_AXIS_CIR_TREADY_REG :
    case CIR_AXIS_CIR_ROW_TVALID_REG:
    case CIR_AXIS_CIR_ROW_TREADY_REG:
    case CIR_DELTA_F_REG :
    case CIR_DELTA_T_REG :
    case CIR_DELTAF_RAW_LO_REG :
    case CIR_DELTAF_RAW_HI_REG :
    case CIR_DELTAT_RAW_LO_REG :
    case CIR_DELTAT_RAW_HI_REG :
    case CIR_STATE_REG :
    case CIR_AM_STATE1_REG :
    case CIR_AM_STATE2_REG :
    case CIR_AM_STATE3_REG :
    case CIR_DO_STATE1_REG :
    case CIR_DO_STATE2_REG :
    case CIR_DO_STATE3_REG :
    case CIR_STRM_STATE1_REG :
    case CIR_STRM_STATE2_REG :
    case CIR_DELTAF_REG :
    case CIR_DELTAF_C_REG :
    case CIR_DELTAT_REG :
    case CIR_DELTAT_C_REG :
    case CIR_AXIS_INSTR_REG :
    case CIR_PRS_RX_DATA_CTR_REG :
    case CIR_PRS_RX_SYMB_CTR_REG :
    case CIR_PRS_RX_FRAME_CTR_REG :
    case CIR_PRS_RX_BLOCK_CTR_REG :
    case CIR_PRS_TX_DATA_CTR_REG :
    case CIR_PRS_TX_SYMB_CTR_REG :
    case CIR_PRS_TX_FRAME_CTR_REG :
    case CIR_PRS_TX_BLOCK_CTR_REG :
    case CIR_BPB_RX_DATA_CTR_REG :
    case CIR_BPB_RX_SYMB_CTR_REG :
    case CIR_BPB_RX_FRAME_CTR_REG :
    case CIR_BPB_RX_BLOCK_CTR_REG :
    case CIR_BPB_TX_DATA_CTR_REG :
    case CIR_BPB_TX_SYMB_CTR_REG :
    case CIR_BPB_TX_FRAME_CTR_REG :
    case CIR_BPB_TX_BLOCK_CTR_REG :
    case CIR_BPA_RX_DATA_CTR_REG :
    case CIR_BPA_RX_SYMB_CTR_REG :
    case CIR_BPA_RX_FRAME_CTR_REG :
    case CIR_BPA_RX_BLOCK_CTR_REG :
    case CIR_BPA_TX_DATA_CTR_REG :
    case CIR_BPA_TX_SYMB_CTR_REG :
    case CIR_BPA_TX_FRAME_CTR_REG :
    case CIR_BPA_TX_BLOCK_CTR_REG :
    case CIR_DTA_RX_DATA_CTR_REG :
    case CIR_DTA_RX_SYMB_CTR_REG :
    case CIR_DTA_RX_FRAME_CTR_REG :
    case CIR_DTA_RX_BLOCK_CTR_REG :
    case CIR_DTA_TX_DATA_CTR_REG :
    case CIR_DTA_TX_SYMB_CTR_REG :
    case CIR_DTA_TX_FRAME_CTR_REG :
    case CIR_DTA_TX_BLOCK_CTR_REG :
    default :
        rv = SlvRegs[ reg ];
        break;
    }
    return( rv );
}


static void SimHwWrReg( u32 reg, u32 data )
{
    if ( reg < CIR_STAT_REG )
    {
        switch( reg )
        {
        case CIR_CTRL_REG :
            switch( data )
            {
            case CIR_RESET       :
            case CIR_SYS_RESET   :
                SimReset();
                break;
            case CIR_ENTER_ADMIN:
                SlvRegs[ CIR_STATE_REG ] &= 0xFFFFFFF0;
                SlvRegs[ CIR_STATE_REG ] |= CIR_EX_ADMIN;
                break;
            case CIR_EXIT_ADMIN :
                SlvRegs[ CIR_STATE_REG ] &= 0xFFFFFFF0;
                SlvRegs[ CIR_STATE_REG ] |= CIR_EX_IDLE;
                break;
            case CIR_AD_RD_BRAM :
                if ( ( SlvRegs [ CIR_STATE_REG ] & 0xF ) == CIR_EX_ADMIN )
                {
                    u64 u = BRAM[ SlvRegs[ CIR_AD_RD_WR_ADDR_REG ] ];
                    u32 ul = u & 0xFFFFFFFF;
                    u32 uh = ( u >> 32 ) & 0xFFFFFFFF;
                    SlvRegs[ CIR_AD_RD_DATA_REG   ] = ul;
                    SlvRegs[ CIR_AD_RD_DATA_REG_H ] = uh;
                    SlvRegs[ CIR_STAT_REG         ]  |= CIR_AD_BRAM_RD_RDY;
                }
                break;
            case CIR_AD_RD_ACK  :
                if ( ( SlvRegs[ CIR_STATE_REG ] & 0xF ) == CIR_EX_ADMIN )
                {
                    SlvRegs [ CIR_STAT_REG ]  &= ~CIR_AD_BRAM_RD_RDY;
                    break;
                }
            case CIR_AD_WR_BRAM :
                if ( ( SlvRegs [ CIR_STATE_REG ] & 0xF ) == CIR_EX_ADMIN )
                {
                    u64 u = SlvRegs[ CIR_AD_WR_DATA_REG_H ];
                    u = ( u << 32 ) | SlvRegs[ CIR_AD_WR_DATA_REG ];
                    BRAM[ SlvRegs [ CIR_AD_RD_WR_ADDR_REG ] ] = u;
                    break;
                }
            case CIR_INT_ACK    :
            case CIR_START      :
            default :
                break;
            }
            break;
            case CIR_INFO_REG :
                break;
            case CIR_N_CPRE_REG :
            case CIR_N_FFT_REG :
            case CIR_N_SYMB_FRAME_REG :
            case CIR_N_FRAME_BLOCK_REG :
            case CIR_AD_RD_WR_ADDR_REG :
            case CIR_AD_WR_DATA_REG :
            case CIR_AD_WR_DATA_REG_H :
            case CIR_LOG_N_FFT_REG :
            case CIR_KFF_REG :
            case CIR_KTF_REG :
            case CIR_KTT_REG :
            case CIR_KFT_REG :
            case CIR_BRAM_CHAN_SEL_REG :
            default:
                SlvRegs[ reg ] = data;
                break;
        }
    }
    else
    {

    }
}



//Simulated CIR

static void ReadPRSData( PrsCIR *pd )
{
    u32 nFFT = SYSPAR( nFFT );
    Cplx32 s;
    AXISink *ps = &RxPrsAXISink;
    for ( u32 i = 0; i < nFFT; i++ )
    {
        u32 u = U32Sext( ps->ReadBram( ps, i ), 12 );
        s = U32toCplx32( u );
        RxPRSTd [ i ] = Cplx32toCplx( s );
    }
}

static void CopyPRS( PrsCIR *pd )
{
    u32 nFFT = SYSPAR( nFFT );
    for ( int i = 0; i < nFFT; i++ ) TxPRSFd [ i ] = conj( Cplx32toCplx( U32toCplx32( SYSPTR( pTxSourceBuf[ i ] ) ) ) );
}

static void CalcFFT( PrsCIR *pd )
{
    u32 nFFT = SYSPAR( nFFT );
    u32 addrMask = nFFT - 1;
    s32 k = -( SYS_N_CHAN_PRSCIR - 1 ) / 2;

    fft( RxPRSFd, RxPRSTd, ( int ) nFFT );
    for( int i = 0; i < SYS_N_CHAN_PRSCIR; i++ )
    {
        cplx *pFFT = FFTPtrs [ i ];
        for ( s32 j = 0; j < nFFT; j++ )
        {
            pFFT [ j ] = TxPRSFd [ ( j + k ) & addrMask ] * RxPRSFd [ j ];
        }
        ifft( pFFT, pFFT, nFFT );
        k++;
    }
}

static void CalcAbs( PrsCIR *pd )
{
    u32 nFFT       = SYSPAR( nFFT );
    u32 absMaxRidx = 0;
    u32 absMaxCidx = 0;
    double absMax  = 0.0;

    for( int ridx = 0; ridx < SYS_N_CHAN_PRSCIR; ridx++ )
    {
        cplx   *pFFT   = FFTPtrs [ ridx ];
        double *pAbs   = AbsPtrs [ ridx ];
        double absm       = 0.0;
        u32    cidx       = 0;

        for ( int j = 0; j < nFFT; j++ ) {
            double d = cabs( pFFT [ j ] );
            pAbs [ j ] = d;
            if ( d > absm )
            {
                absm  = d;
                cidx  = j;
            }
        }

        if ( absm > absMax ) {
            absMax     = absm;
            absMaxRidx = ridx;
            absMaxCidx = cidx;
        }
    }

    RowIdx = absMaxRidx;
    ColIdx = absMaxCidx;
}

static void CalcDelta( PrsCIR *pd )
{
    u32 nFFT     = SYSPAR( nFFT );
    u32 addrMask = nFFT - 1;
    cplx *pRm    = FFTPtrs [ MAX( RowIdx - 1, 0 ) ];
    cplx *pR0    = FFTPtrs [ RowIdx ];
    cplx *pRp    = FFTPtrs [ MIN( RowIdx + 1, SYS_N_CHAN_PRSCIR - 1 ) ];
    u32   pCm    = ( ColIdx - 1 ) & addrMask;
    u32   pC0    = ColIdx;
    u32   pCp    = ( ColIdx + 1 ) & addrMask;

    cplx w1 = pRp [ pC0 ];
    cplx w2 = pRm [ pC0 ];
    cplx z1 = pR0 [ pC0 ] - pRp [ pC0 ];
    cplx z2 = pRm [ pC0 ] - pR0 [ pC0 ];



    cplx w  = pR0 [ pCm ] - pR0 [ pCp ];
    cplx z  = pR0 [ pCm ] + pR0 [ pCp ] + ( M_SQRT2 + 0.0 * I ) * pR0 [ pC0 ];

    u32 fftIdx = FFTIndex( SYSPAR( nFFT) );
    DVCPAR( DeltaF  ) = RowIdx - ( SYS_N_CHAN_PRSCIR - 1 ) / 2;
    double deltaf     = -creal( w1 * conj( z1 ) + w2 * conj( z2 ) ) / ( cabs( z1 ) * cabs( z1 ) + cabs( z2 ) * cabs( z2 ) );
    DVCPAR( DeltaT  ) = ( pC0 < nFFT / 2 - 1 ? pC0 : pC0 - nFFT );
    double deltat     = -creal( w * conj( z ) ) / ( cabs( z ) * cabs( z ) );
    DVCPAR( deltafp ) = deltaf * DVCPAR( KFF [ fftIdx ] ) + deltat * DVCPAR( KTF [ fftIdx ] );
    DVCPAR( deltatp ) = deltat * DVCPAR( KTT [ fftIdx ] ) + deltaf * DVCPAR( KFT [ fftIdx ] );
}



static void Calc( PrsCIR * pd )
{
    nPass = 0;
    DeltaTAcc = 0.0;
    DeltaFAcc = 0.0;
    ReadPRSData( pd );
    CopyPRS( pd );
    CalcDeltas( pd );
}

static void Recalc( PrsCIR *pd )
{
    u32 nFFT     = SYSPAR( nFFT );
    for ( int i = 0; i < nFFT; i++ )
    {
        RxPRSTd  [ i ] = RxPRSTdCorFT [ i ];
    }
    CalcDeltas( pd );
}

static void CalcDeltas( PrsCIR *pd )
{
    CalcFFT( pd );
    CalcAbs( pd );
    CalcDelta( pd );
    CorrectRxPRS( pd );
    nPass = nPass + 1;
    DeltaTAcc += ( DVCPAR( DeltaT ) + DVCPAR( deltatp ) );
    DeltaFAcc += ( DVCPAR( DeltaF ) + DVCPAR( deltafp ) );
}


static void CorrectRxPRS( PrsCIR *pd )
{
    u32 nFFT = SYSPAR( nFFT );
    u32 idxMask = nFFT - 1;
    for ( int i = 0; i < nFFT; i++ )
    {
         RxPRSTdCorDeltaT [ i ]  = RxPRSTd [ ( i + DVCPAR( DeltaT ) ) & idxMask ];

         SYSPTR( pRxPRSTdCorDeltaT [ i ] ) = RxPRSTdCorDeltaT [ i ];
         RxPRSTdCordeltafp [ i ] = RxPRSTdCorDeltaT [ i ];
    }


    //cplx cordeltafp = cexp( I * 2 * M_PI * ( ( double ) DVCPAR( DeltaF ) + DVCPAR( deltafp ) ) / ( double ) ( nFFT ) );
    for ( int j = 0; j < nFFT; j++ )
    {
          cplx cordeltafp = cexp( I * 2 * M_PI * ( ( double ) DVCPAR( DeltaF ) + DVCPAR( deltafp ) ) * ( double ) ( j + 1 ) / ( double ) ( nFFT ) );
          RxPRSTdCordeltafp [ j ] = RxPRSTdCordeltafp [ j ] * cordeltafp;
//        for ( int k = j; k < nFFT; k++ )
//        {
//            RxPRSTdCordeltafp [ k ] = RxPRSTdCordeltafp [ k ] * cordeltafp;
//        }
    }
//    for ( int i = 0; i < nFFT; i++ )
//    {
//         cplx cordeltafp = cexp( I * 2 * M_PI * ( ( double ) DVCPAR( DeltaF ) + DVCPAR( deltafp ) ) * ( double ) ( i + 1 ) / ( double ) ( nFFT ) );
//
//         RxPRSTdCorDeltaT [ i ]  = RxPRSTd [ ( i + DVCPAR( DeltaT ) ) & idxMask ];
//
//         SYSPTR( pRxPRSTdCorDeltaT [ i ] ) = RxPRSTdCorDeltaT [ i ];
//         RxPRSTdCordeltafp [ i ] = RxPRSTdCorDeltaT [ i ] * cordeltafp;
////         RxPRSTdCorFT [ i ] = RxPRSTdCordeltafp [ i ];
//    }
    fft( RxPRSFdCordeltafp, RxPRSTdCordeltafp, nFFT );
    for ( int i = 0; i < nFFT; i++ )
    {
         TimeCorrSampleIdx [ i ] = ( i < nFFT / 2 ? i : i - nFFT );
         cplx cordeltatp = cexp( I * 2 * M_PI * ( ( double ) DVCPAR( deltatp ) * ( double ) TimeCorrSampleIdx [ i ] ) / ( double ) ( nFFT ));
         RxPRSFdCorFT [ i ]            = RxPRSFdCordeltafp [ i ] * cordeltatp;
         SYSPTR( pRxPRSFdCorFT [ i ] ) = RxPRSFdCorFT [ i ];
    }
    ifft( RxPRSTdCorFT, RxPRSFdCorFT, nFFT );
    for ( int i = 0; i < nFFT; i++ )
    {
        SYSPTR( pRxPRSTdCorFT[ i ] ) = RxPRSTdCorFT [ i ];
    }

}

