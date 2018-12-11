/* * AXI4_Stream_Data_Source.c
 *
 *  Created on: 20 May 2017
 *      Author: JULIAN MORTIMER
 */
#include "include.h"


// =================================
// Boilerplate Function Declarations
// =================================
void AXISourceInit(
    AXISource *pd,
    const char      *Name,
    u32             BaseAddr,
    SysCtrl         *pSysCtrl,
    SysDvce         *pSysDvce,
    SysDesc         *pSysDesc,
    SysPtrs         *pSysPtrs,
    SysData         *pSysData,
    AXISourceDesc *pDvcDesc,
    AXISourceData *pDvcData
    );
static void SetParams    ( AXISource *pd );
static void SysDvceDump  ( AXISource *pd );
static void SysDescDump  ( AXISource *pd );
static void SysDataDump  ( AXISource *pd );
static void DvcDataInit  ( AXISource *pd );
static void DvcDataUpd   ( AXISource *pd );
static void DvcDescDump  ( AXISource *pd );
static void DvcDataDump  ( AXISource *pd );

// =================================================
// Device Structure Member Functions (public access)
// =================================================
static inline void Start1Shot    ( AXISource *pd ) { hwWriteCtrlReg( SYSDVC( BaseAddr ), AXISOURCE_CTRL_REG << 2, AXISOURCE_START_1_SHOT ); }
static inline void StartContin   ( AXISource *pd ) { hwWriteCtrlReg( SYSDVC( BaseAddr ), AXISOURCE_CTRL_REG << 2, AXISOURCE_START_CONTIN ); }
static inline void Stop          ( AXISource *pd ) { hwWriteCtrlReg( SYSDVC( BaseAddr ), AXISOURCE_CTRL_REG << 2, AXISOURCE_STOP );         }
static inline void SetTxAddr     ( AXISource *pd, u32 addr ) { DVCPAR( TxAddr    ) = addr; SYSWRSLVREG( pd, AXISOURCE_TX_ADDR_REG   , addr ); }
static inline void SetBlockSize  ( AXISource *pd, u32 size ) { DVCPAR( BlockSize ) = size; SYSWRSLVREG( pd, AXISOURCE_BLOCK_SIZE_REG, size ); }
static void        TxStream      ( AXISource *pd, u32 RxAddr, u32 BlockSize );
static u32         ReadBram      ( AXISource *pd, u32 addr );
static void        WriteBram     ( AXISource *pd, u32 addr, u32 data );
static void        InitBram      ( AXISource *pd, u32 addr, u32 count, u32 data );
static void        LoadBram      ( AXISource *pd );
static void        DumpBram      ( AXISource *pd, u32 addr, u32 nBits );
static inline u32  HwBramSize    ( AXISource *pd ) { return( hwReadReg( SYSDVC( BaseAddr ), AXISOURCE_BRAM_SIZE_REG << 2 ) ); }
static void        Dump          ( AXISource *pd );

// Declare init structure here with member functions
AXISource AXISourceInitStr = {
    NULL,
    NULL,
    NULL,
    Start1Shot,
    StartContin,
    Stop,
    SetTxAddr,
	SetBlockSize,
    TxStream,
    ReadBram,
    WriteBram,
    InitBram,
    LoadBram,
    DumpBram,
    HwBramSize,
    Dump
};

//===================
// Initialize Device
//===================
void AXISourceInit( AXISource *pd, const char *Name, u32 BaseAddr, SysCtrl *pSysCtrl, SysDvce *pSysDvce, SysDesc *pSysDesc, SysPtrs *pSysPtrs, SysData *pSysData, AXISourceDesc *pDvcDesc, AXISourceData *pDvcData )
{
	*pd             = AXISourceInitStr;
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
        AXISOURCE_HW_INFO_REG,
        AXISOURCE_CTRL_REG,
        AXISOURCE_INFO_REG,
        AXISOURCE_STAT_REG,
        AXISOURCE_RESET,
        AXISOURCE_INT_ACK,
        0,
        0,
        0,
        0,
        0,
        AXISOURCE_TX_DATA_CTR_REG,
        AXISOURCE_TX_SYMBOL_CTR_REG,
        AXISOURCE_TX_FRAME_CTR_REG,
        AXISOURCE_TX_BLOCK_CTR_REG,
        0,
        0
        );
    pd->pDvcDesc  = pDvcDesc;
    pd->pDvcData  = pDvcData;
    DvcDataInit( pd );
    pd->InitBram( pd, 0, DVCPAR( BramSize ), DVCPAR( BramInitVal ));
    SYSRESET( pd );
    SetParams( pd );
}

//========================
// Dump Device to Console
//========================
static void Dump( AXISource *pd ) {
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
static const char s_TF_UNK0             [ ] = "UNKNOWN";
static const char s_TF_UNK1             [ ] = "UNKNOWN";
static const char s_TF_UNK2             [ ] = "UNKNOWN";
static const char s_TF_UNK3             [ ] = "UNKNOWN";
static const char s_TF_UNK4             [ ] = "UNKNOWN";
static const char s_TF_UNK5             [ ] = "UNKNOWN";
static const char s_TF_UNK6             [ ] = "UNKNOWN";
static const char s_TF_UNK7             [ ] = "UNKNOWN";
static const char s_TF_UNK8             [ ] = "UNKNOWN";
static const char s_TF_IDLE             [ ] = "TF_IDLE";
static const char s_TF_WAIT_1_SHOT      [ ] = "TF_WAIT_1_SHOT";
static const char s_TF_WAIT_CONTIN      [ ] = "TF_WAIT_CONTIN";
static const char s_TF_UNKC             [ ] = "UNKNOWN";
static const char s_TF_UNKD             [ ] = "UNKNOWN";
static const char s_TF_TX_1_SHOT        [ ] = "TF_TX_1_SHOT";
static const char s_TF_TX_CONTIN        [ ] = "TF_TX_CONTIN";

static const char *TfStateStr [ 16 ] = {
         s_TF_UNK0,
         s_TF_UNK1,
         s_TF_UNK2,
         s_TF_UNK3,
         s_TF_UNK4,
         s_TF_UNK5,
         s_TF_UNK6,
         s_TF_UNK7,
         s_TF_UNK8,
         s_TF_IDLE,
         s_TF_WAIT_1_SHOT,
         s_TF_WAIT_CONTIN,
         s_TF_UNKC,
         s_TF_UNKD,
         s_TF_TX_1_SHOT,
         s_TF_TX_CONTIN
};

static const char s_TM_UNK0             [ ] = "UNKNOWN";
static const char s_TM_UNK1             [ ] = "UNKNOWN";
static const char s_TM_UNK2             [ ] = "UNKNOWN";
static const char s_TM_UNK3             [ ] = "UNKNOWN";
static const char s_TM_UNK4             [ ] = "UNKNOWN";
static const char s_TM_UNK5             [ ] = "UNKNOWN";
static const char s_TM_UNK6             [ ] = "UNKNOWN";
static const char s_TM_UNK7             [ ] = "UNKNOWN";
static const char s_TM_UNK8             [ ] = "UNKNOWN";
static const char s_TM_IDLE             [ ] = "TM_IDLE";
static const char s_TM_UNKA             [ ] = "UNKNOWN";
static const char s_TM_WAIT_BRAM        [ ] = "TM_WAIT_BRAM";
static const char s_TM_UNKC             [ ] = "UNKNOWN";
static const char s_TM_RUNNING          [ ] = "TM_RUNNING";
static const char s_TM_UNKE             [ ] = "UNKNOWN";
static const char s_TM_UNKF             [ ] = "UNKNOWN";

static const char *TmStateStr [ 16 ] = {
         s_TM_UNK0,
         s_TM_UNK1,
         s_TM_UNK2,
         s_TM_UNK3,
         s_TM_UNK4,
         s_TM_UNK5,
         s_TM_UNK6,
         s_TM_UNK7,
         s_TM_UNK8,
         s_TM_IDLE,
         s_TM_UNKA,
         s_TM_WAIT_BRAM,
         s_TM_UNKC,
         s_TM_RUNNING,
         s_TM_UNKE,
         s_TM_UNKF
};


static void SetParams( AXISource *pd )
{
	SYSWRSLVREG( pd,  AXISOURCE_N_FFT_REG,         SYSPAR( nFFT        ));
	SYSWRSLVREG( pd,  AXISOURCE_N_SYMB_FRAME_REG,  SYSPAR( nSymbData   ) + SYSPAR( nSymbPRS ));
	SYSWRSLVREG( pd,  AXISOURCE_N_FRAME_BLOCK_REG, SYSPAR( nFrameBlock ));
	SYSWRSLVREG( pd,  AXISOURCE_TX_ADDR_REG,       DVCPAR( TxAddr      ));
	SYSWRSLVREG( pd,  AXISOURCE_BLOCK_SIZE_REG,    DVCPAR( BlockSize   ));
}

static void DvcDataInit( AXISource *pd )
{
   DVCDAT( AxsTxFrameDone          ) = false;
   DVCDAT( AxsMAxisTready          ) = false;
   DVCDAT( AxsMAxisTvalid          ) = false;
   DVCDAT( AxsTfState              ) = 0;
   DVCDAT( AxsTmState              ) = 0;
}

static void DvcDataUpd( AXISource *pd )
{
    u32 statreg = SYSRDSLVREG( pd, AXISOURCE_STAT_REG );
    DVCDAT( AxsTxFrameDone ) = ( statreg & AXISOURCE_TX_FRAME_DONE_MASK ) ? true : false;
    DVCDAT( AxsMAxisTready ) = ( statreg & AXISOURCE_M_AXIS_TREADY_MASK ) ? true : false;
    DVCDAT( AxsMAxisTvalid ) = ( statreg & AXISOURCE_M_AXIS_TVALID_MASK ) ? true : false;
    DVCDAT( AxsTfState     ) = ( statreg & AXISOURCE_TF_STATE_MASK      ) >> AXISOURCE_TF_STATE_SHR;
    DVCDAT( AxsTmState     ) = ( statreg & AXISOURCE_TM_STATE_MASK      ) >> AXISOURCE_TM_STATE_SHR;
}

static void DvcDataDump( AXISource *pd )
{
	DvcDataUpd( pd );
    cprintf("\ntx_frame_done      : %s",   OFFON( DVCDAT( AxsTxFrameDone  ) ) );
	cprintf("\nm_axis_data_tready : %s",   OFFON( DVCDAT( AxsMAxisTready  ) ) );
    cprintf("\nm_axis_data_tvalid : %s",   OFFON( DVCDAT( AxsMAxisTvalid  ) ) );
    cprintf("\ntf_state           : %s",   TfStateStr[ DVCDAT( AxsTfState ) ] );
    cprintf("\ntm_state           : %s",   TmStateStr[ DVCDAT( AxsTmState ) ] );
    cprintf("\n" );
}

static void DvcDescDump( AXISource *pd )
{
	cprintf("\nBRAM Size          : %u",   toUint( DVCPAR( BramSize    )) );
	cprintf("\nHW BRAM Size       : %u",   toUint( pd->HwBramSize( pd  )) );
	cprintf("\nBRAM Init Val      : %08X", toUint( DVCPAR( BramInitVal )) );
	cprintf("\nTx address         : %u",   toUint( DVCPAR( TxAddr      )) );
	cprintf("\nBlock Size         : %u",   toUint( DVCPAR( BlockSize   )) );
	cprintf("\n");
}

static void SysDvceDump( AXISource *pd )
{
	SYSDVCEDUMP( pd );
}



static void SysDataDump( AXISource *pd )
{
	SYSDATADUMP( pd );
}

static void SysDescDump( AXISource *pd )
{

	cprintf("\nnFFT               : %u",   toUint( SYSPAR( nFFT        ) ) );
    cprintf("\nnCPre              : %u",   toUint( SYSPAR( nCPre       ) ) );
	cprintf("\nnSymbFrame         : %u",   toUint( SYSPAR( nSymbPRS ) + SYSPAR( nSymbData ) ) );
	cprintf("\nnFrameBlock        : %u",   toUint( SYSPAR( nFrameBlock ) ) );
	cprintf("\n");
}

// Insert device-specific code here
//======================================
// Device-specific non member functions
//======================================
// Example: static inline void EnterAdmin    ( AXISource *pd ) { hwWriteCtrlReg( SYSDVC( BaseAddr ), AXISOURCE_CTRL_REG << 2, AXISOURCE_ENTER_ADMIN ); }

//=================================
// Device-specific member functions
//=================================

static u32 ReadBram( AXISource *pd, u32 addr )
{
	return( DVCPAR( pBramC )->ReadBram( DVCPAR( pBramC ), addr ));
}

static void WriteBram( AXISource *pd, u32 addr, u32 data )
{
	DVCPAR( pBramC )->WriteBram( DVCPAR( pBramC ), addr, data );
}

static void InitBram( AXISource *pd, u32 addr, u32 count, u32 data )
{
	DVCPAR( pBramC )->InitBram( DVCPAR( pBramC ), addr, count, data );
}

static void LoadBram( AXISource *pd )
{
	DVCPAR( pBramC )->LoadBram( DVCPAR( pBramC ) );
}

static void DumpBram( AXISource *pd, u32 addr, u32 nbits)
{
	char Name[ SYS_MAX_OBJ_NAME_LEN + 1 ]; strcpy( Name, SYSDVC( Name ) );
	if ( nbits == 0 )
		utilScreenDumpData( DUMP_SEL_SOURCE, DUMP_FMT_HEX,       DUMP_TYP_MEMORY, addr, 4,     0, ( dmpArg )pd, strcat( Name, "Tx Buffer:" ));
	else if(( nbits >= 8 ) && ( nbits <= 16 ))
     	utilScreenDumpData( DUMP_SEL_SOURCE, DUMP_FMT_IQ_SIGNED, DUMP_TYP_MEMORY, addr, nbits, 0, ( dmpArg )pd, strcat( Name, "Tx Buffer:" ));
	else
		cprintf( "Format must be 0 (hex) or 8-16 (decimal Q:I)\n" );
}


static void TxStream( AXISource *pd, u32 TxAddr, u32 BlockSize )
{
	//SYSRESET( pd );
	DVCPAR( TxAddr )    = TxAddr;
	DVCPAR( BlockSize ) = BlockSize;
	SYSWRSLVREG( pd,  AXISOURCE_TX_ADDR_REG,       DVCPAR( TxAddr    ));
	SYSWRSLVREG( pd,  AXISOURCE_BLOCK_SIZE_REG,    DVCPAR( BlockSize ));
    pd->Start1Shot( pd );
}
