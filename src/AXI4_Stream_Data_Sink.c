/*
 * AXI4_Stream_Data_Sink.c

 *
 *  Created on: 16 May 2017
 *      Author: JULIAN MORTIMER
 */
#include "include.h"


// =================================
// Boilerplate Function Declarations
// =================================
void AXISinkInit(
		AXISink *pd,
		const char *Name,
		u32        BaseAddr,
		SysCtrl    *pSysCtrl,
		SysDvce    *pSysDvce,
		SysDesc    *pSysDesc,
		SysPtrs    *pSysPtrs,
		SysData    *pSysData,
		AXISinkDesc *pDvcDesc,
		AXISinkData *pDvcData
);
static void SetParams    ( AXISink *pd );
static void SysDvceDump  ( AXISink *pd );
static void SysDescDump  ( AXISink *pd );
static void SysDataDump  ( AXISink *pd );
static void DvcDataInit  ( AXISink *pd );
static void DvcDataUpd   ( AXISink *pd );
static void DvcDescDump  ( AXISink *pd );
static void DvcDataDump  ( AXISink *pd );

// =================================================
// Device Structure Member Functions (public access)
// =================================================
static inline void Start1Shot      ( AXISink *pd ) { hwWriteCtrlReg( SYSDVC( BaseAddr ), AXISINK_CTRL_REG << 2, AXISINK_START_RX ); }
static inline void StartContin     ( AXISink *pd ) { hwWriteCtrlReg( SYSDVC( BaseAddr ), AXISINK_CTRL_REG << 2, AXISINK_START );    }
static inline void Stop            ( AXISink *pd ) { hwWriteCtrlReg( SYSDVC( BaseAddr ), AXISINK_CTRL_REG << 2, AXISINK_STOP );     }
static inline void SetRxAddr       ( AXISink *pd, u32 addr ) { DVCPAR( RxAddr )    = addr; SYSWRSLVREG( pd, AXISINK_RX_ADDR_REG    , addr ); }
static inline void SetBlockSize    ( AXISink *pd, u32 size ) { DVCPAR( BlockSize ) = size; SYSWRSLVREG( pd, AXISINK_BLOCK_SIZE_REG , size ); }
static u32         TestBram        ( AXISink *pd );
static u32         ReadBram        ( AXISink *pd, u32 addr );
static void        WriteBram       ( AXISink *pd, u32 addr, u32 data );
static void        InitBram        ( AXISink *pd, u32 addr, u32 count, u32 data );
static void        DumpBram        ( AXISink *pd, u32 addr, u32 nBits );
static inline u32  HwBramSize      ( AXISink *pd ) { return( hwReadReg( SYSDVC( BaseAddr ), AXISINK_BRAM_SIZE_REG << 2 ) ); }
static inline void WriteInhibit    ( AXISink *pd ) { SYSWRSLVREG( pd, AXISINK_CFG_REG, SYSRDSLVREG( pd, AXISINK_CFG_REG ) |  AXISINK_WRITE_INHIBIT ); }
static inline void WriteEnable     ( AXISink *pd ) { SYSWRSLVREG( pd, AXISINK_CFG_REG, SYSRDSLVREG( pd, AXISINK_CFG_REG ) & ~AXISINK_WRITE_INHIBIT ); }
static void        RxStream        ( AXISink *pd, u32 RxAddr, u32 BlockSize );
static void        Dump            ( AXISink *pd );

// Declare init structure here with member functions
AXISink AXISinkInitStr = {
		NULL,
		NULL,
		NULL,
		Start1Shot      ,
		StartContin     ,
		Stop            ,
		SetRxAddr       ,
		SetBlockSize    ,
		WriteInhibit    ,
		WriteEnable     ,
		RxStream        ,
		TestBram        ,
		ReadBram        ,
		WriteBram       ,
		InitBram        ,
		DumpBram        ,
		HwBramSize      ,
		Dump
};

//===================
// Initialize Device
//===================
void AXISinkInit( AXISink *pd, const char *Name, u32 BaseAddr, SysCtrl *pSysCtrl, SysDvce *pSysDvce, SysDesc *pSysDesc, SysPtrs *pSysPtrs, SysData *pSysData, AXISinkDesc *pDvcDesc, AXISinkData *pDvcData )
{
	*pd             = AXISinkInitStr;
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
			AXISINK_HW_INFO_REG,
			AXISINK_CTRL_REG,
			AXISINK_INFO_REG,
			AXISINK_STAT_REG,
			AXISINK_RESET,
			AXISINK_INT_ACK,
			0,
			AXISINK_RX_DATA_CTR_REG,
			AXISINK_RX_SYMB_CTR_REG,
			AXISINK_RX_FRAME_CTR_REG,
			AXISINK_RX_BLOCK_CTR_REG,
			0,
			0,
			0,
			0,
			0,
			0
	);
	pd->pDvcDesc  = pDvcDesc;
	pd->pDvcData  = pDvcData;
	DvcDataInit( pd );
	SYSRESET( pd );
	( void ) pd->TestBram( pd );
	pd->InitBram( pd, 0, pd->HwBramSize( pd ) / 4, DVCPAR( BramInitVal ) );
	SetParams( pd );
}

//========================
// Dump Device to Console
//========================
static void Dump( AXISink *pd ) {
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
static void SetParams( AXISink *pd )
{
	SYSWRSLVREG( pd,  AXISINK_N_FFT_REG,         SYSPAR( nFFT ));
	SYSWRSLVREG( pd,  AXISINK_N_SYMB_FRAME_REG,  SYSPAR( nSymbData ));
	SYSWRSLVREG( pd,  AXISINK_N_FRAME_BLOCK_REG, SYSPAR( nFrameBlock ));
	SYSWRSLVREG( pd,  AXISINK_RX_ADDR_REG,       DVCPAR( RxAddr ));
	SYSWRSLVREG( pd,  AXISINK_BLOCK_SIZE_REG,    DVCPAR( BlockSize ));
	SYSWRSLVREG( pd,  AXISINK_CFG_REG,           0 );
}

static void DvcDataInit( AXISink *pd )
{
	DVCDAT( WrDataCtr  ) = 0;
	DVCDAT( WrSymbCtr  ) = 0;
	DVCDAT( WrFrameCtr ) = 0;
	DVCDAT( WrBlockCtr ) = 0;
}

static void DvcDataUpd( AXISink *pd )
{
    u32 statReg = SYSRDSLVREG( pd, AXISINK_STAT_REG );
	DVCDAT( WrDataCtr  ) = SYSRDSLVREG( pd, AXISINK_WR_DATA_CTR_REG );
	DVCDAT( WrSymbCtr  ) = SYSRDSLVREG( pd, AXISINK_WR_SYMB_CTR_REG );
	DVCDAT( WrFrameCtr ) = SYSRDSLVREG( pd, AXISINK_WR_FRAME_CTR_REG );
	DVCDAT( WrBlockCtr ) = SYSRDSLVREG( pd, AXISINK_WR_BLOCK_CTR_REG );
    DVCDAT( rd_data_rdy )        = ( statReg & AXISINK_RD_DATA_RDY    ? true : false );
    DVCDAT( rx_stream_done )     = ( statReg & AXISINK_RX_STREAM_DONE ? true : false );
	DVCDAT( s_axis_data_tvalid ) = ( statReg & AXISINK_S_AXIS_DATA_TVALID ? true : false );
    DVCDAT( s_axis_data_tready ) = ( statReg & AXISINK_S_AXIS_DATA_TREADY ? true : false );
    DVCDAT( ex_state )           = ( statReg & AXISINK_EX_STATE_MASK ) >> AXISINK_EX_STATE_SHR;
    DVCDAT( rx_state )           = ( statReg & AXISINK_RX_STATE_MASK ) >> AXISINK_RX_STATE_SHR;
    DVCDAT( rd_data_state )      = ( statReg & AXISINK_RD_DATA_STATE_MASK ) >> AXISINK_RD_DATA_STATE_SHR;
}

static const char s_EX_UNK0           [ ] = "UNKNOWN";     // 0
static const char s_EX_UNK1           [ ] = "UNKNOWN";     // 1
static const char s_EX_UNK2           [ ] = "UNKNOWN";     // 2
static const char s_EX_UNK3           [ ] = "UNKNOWN";     // 3
static const char s_EX_UNK4           [ ] = "UNKNOWN";     // 4
static const char s_EX_UNK5           [ ] = "UNKNOWN";     // 5
static const char s_EX_UNK6           [ ] = "UNKNOWN";     // 6
static const char s_EX_UNK7           [ ] = "UNKNOWN";     // 7
static const char s_EX_UNK8           [ ] = "UNKNOWN";     // 8
static const char s_EX_IDLE           [ ] = "EX_IDLE";     // 9
static const char s_EX_UNKA           [ ] = "UNKNOWN";     // A
static const char s_EX_SOFT_LOAD_DATA [ ] = "EX_SOFT_LOAD_DATA"; // B
static const char s_EX_UNKC           [ ] = "UNKNOWN";     // C
static const char s_EX_UNKD           [ ] = "UNKNOWN";     // D
static const char s_EX_UNKE           [ ] = "UNKNOWN";     // E
static const char s_EX_UNKF           [ ] = "UNKNOWN";     // F
static const char *EXStateStr [ 16 ] = {
        s_EX_UNK0,
        s_EX_UNK1,
        s_EX_UNK2,
        s_EX_UNK3,
        s_EX_UNK4,
        s_EX_UNK5,
        s_EX_UNK6,
        s_EX_UNK7,
        s_EX_UNK8,
        s_EX_IDLE,
        s_EX_UNKA,
        s_EX_SOFT_LOAD_DATA,
        s_EX_UNKC,
        s_EX_UNKD,
        s_EX_UNKE,
        s_EX_UNKF
};

static const char s_RX_UNK0           [ ] = "UNKNOWN";     // 0
static const char s_RX_UNK1           [ ] = "UNKNOWN";     // 1
static const char s_RX_UNK2           [ ] = "UNKNOWN";     // 2
static const char s_RX_UNK3           [ ] = "UNKNOWN";     // 3
static const char s_RX_UNK4           [ ] = "UNKNOWN";     // 4
static const char s_RX_UNK5           [ ] = "UNKNOWN";     // 5
static const char s_RX_UNK6           [ ] = "UNKNOWN";     // 6
static const char s_RX_UNK7           [ ] = "UNKNOWN";     // 7
static const char s_RX_UNK8           [ ] = "UNKNOWN";     // 8
static const char s_RX_IDLE           [ ] = "RX_IDLE";     // 9
static const char s_RX_UNKA           [ ] = "UNKNOWN";     // A
static const char s_RX_BLCK           [ ] = "RX_BLCK";     // B
static const char s_RX_UNKC           [ ] = "UNKNOWN";     // C
static const char s_RX_CONT           [ ] = "RX_CONT";     // D
static const char s_RX_UNKE           [ ] = "UNKNOWN";     // E
static const char s_RX_UNKF           [ ] = "UNKNOWN";     // F
static const char *RXStateStr [ 16 ] = {
        s_RX_UNK0,
        s_RX_UNK1,
        s_RX_UNK2,
        s_RX_UNK3,
        s_RX_UNK4,
        s_RX_UNK5,
        s_RX_UNK6,
        s_RX_UNK7,
        s_RX_UNK8,
        s_RX_IDLE,
        s_RX_UNKA,
        s_RX_BLCK,
        s_RX_UNKC,
        s_RX_CONT,
        s_RX_UNKE,
        s_RX_UNKF
};


static const char s_RD_UNK0           [ ] = "UNKNOWN";     // 0
static const char s_RD_UNK1           [ ] = "UNKNOWN";     // 1
static const char s_RD_UNK2           [ ] = "UNKNOWN";     // 2
static const char s_RD_UNK3           [ ] = "UNKNOWN";     // 3
static const char s_RD_UNK4           [ ] = "UNKNOWN";     // 4
static const char s_RD_UNK5           [ ] = "UNKNOWN";     // 5
static const char s_RD_UNK6           [ ] = "UNKNOWN";     // 6
static const char s_RD_UNK7           [ ] = "UNKNOWN";     // 7
static const char s_RD_UNK8           [ ] = "UNKNOWN";     // 8
static const char s_RD_IDLE           [ ] = "RD_IDLE";     // 9
static const char s_RD_WAITING        [ ] = "RD_WAITING";  // A
static const char s_RD_UNKB           [ ] = "UNKNOWN";     // B
static const char s_RD_UNKC           [ ] = "UNKNOWN";     // C
static const char s_RD_READY          [ ] = "RD_READY";    // D
static const char s_RD_UNKE           [ ] = "UNKNOWN";     // E
static const char s_RD_UNKF           [ ] = "UNKNOWN";     // F
static const char *RdDataStateStr [ 16 ] = {
        s_RD_UNK0,
        s_RD_UNK1,
        s_RD_UNK2,
        s_RD_UNK3,
        s_RD_UNK4,
        s_RD_UNK5,
        s_RD_UNK6,
        s_RD_UNK7,
        s_RD_UNK8,
        s_RD_IDLE,
        s_RD_WAITING,
        s_RD_UNKB,
        s_RD_UNKC,
        s_RD_READY,
        s_RD_UNKE,
        s_RD_UNKF
};

static void DvcDataDump( AXISink *pd )
{
	DvcDataUpd( pd );
	cprintf( "\nWr data ctr        : %u", toUint( DVCDAT( WrDataCtr  )));
	cprintf( "\nWr symbol ctr      : %u", toUint( DVCDAT( WrSymbCtr  )));
	cprintf( "\nWr frame ctr       : %u", toUint( DVCDAT( WrFrameCtr )));
	cprintf( "\nWr block ctr       : %u", toUint( DVCDAT( WrBlockCtr )));
    cprintf( "\nrd_data_rdy        : %s", OFFON( DVCDAT( rd_data_rdy ) ) );
    cprintf( "\nrx_stream_done     : %s", OFFON( DVCDAT( rx_stream_done ) ) );
	cprintf( "\ns_axis_data_tvalid : %s", OFFON( DVCDAT( s_axis_data_tvalid ) ) );
    cprintf( "\ns_axis_data_tready : %s", OFFON( DVCDAT( s_axis_data_tready ) ) );
    cprintf( "\nex_state           : %s", EXStateStr     [ DVCDAT( ex_state  ) ] );
    cprintf( "\nrx_state           : %s", RXStateStr     [ DVCDAT( rx_state ) ] );
    cprintf( "\nrd_data_state      : %s", RdDataStateStr [ DVCDAT( rd_data_state ) ] );
	cprintf( "\n" );
}

static void DvcDescDump( AXISink *pd )
{
	cprintf("\nBRAM size           : %u dwords", toUint( DVCPAR( BramSize    )) );
	cprintf("\nHw BRAM size        : %u bytes",  toUint( pd->HwBramSize( pd  )) );
	cprintf("\nBRAM init val       : %08X",      toUint( DVCPAR( BramInitVal )) );
	cprintf("\nRx address          : %u",        toUint( DVCPAR( RxAddr      )) );
	cprintf("\nBlock size          : %u",        toUint( DVCPAR( BlockSize   )) );
	cprintf("\n");
}

static void SysDvceDump( AXISink *pd )
{
	SYSDVCEDUMP( pd );
}

static void SysDataDump( AXISink *pd )
{
	SYSDATADUMP( pd );
}

static void SysDescDump( AXISink *pd )
{
	cprintf("\nFFT length          : %u",   toUint( SYSPAR( nFFT        )) );
	cprintf("\nCyclic prefix       : %u",   toUint( SYSPAR( nCPre       )) );
	cprintf("\nSymbols/frame       : %u",   toUint( SYSPAR( nSymbData   )) );
	cprintf("\nFrames/block        : %u",   toUint( SYSPAR( nFrameBlock )) );
	cprintf("\n");
}



// Insert device-specific code here
//======================================
// Device-specific non member functions
//======================================
static inline void LoadData        ( AXISink *pd ) { hwWriteCtrlReg( SYSDVC( BaseAddr ), AXISINK_CTRL_REG << 2, AXISINK_LOAD_DATA);    }
static inline void EndLoad         ( AXISink *pd ) { hwWriteCtrlReg( SYSDVC( BaseAddr ), AXISINK_CTRL_REG << 2, AXISINK_END_LOAD );    }

//=================================
// Device-specific member functions
//=================================

static void WriteBram( AXISink *pd, u32 addr, u32 data )
{
	LoadData( pd );
	hwWriteData( SYSDVC( BaseAddr ),
			AXISINK_CTRL_REG << 2, AXISINK_WR_DATA,
			AXISINK_ADDR_REG << 2, addr,
			AXISINK_WR_DATA_REG << 2, data
	);
	EndLoad( pd );
}

static u32 ReadBram( AXISink *pd, u32 addr )
{
	return hwReadData( SYSDVC( BaseAddr ),
			AXISINK_CTRL_REG << 2, AXISINK_RD_DATA, AXISINK_RD_ACK,
			AXISINK_STAT_REG << 2, AXISINK_RD_DATA_RDY,
			AXISINK_ADDR_REG << 2, addr,
			AXISINK_RD_DATA_REG << 2
	);
}

static void InitBram( AXISink *pd, u32 addr, u32 count, u32 val )
{
	LoadData( pd );
	for ( u32 i = 0; i < count; i++ )
	{
		hwWriteData( SYSDVC( BaseAddr ),
				AXISINK_CTRL_REG << 2, AXISINK_WR_DATA,
				AXISINK_ADDR_REG << 2, i + addr,
				AXISINK_WR_DATA_REG << 2, val
		);

	}
	EndLoad( pd );
}

static void DumpBram( AXISink *pd, u32 addr, u32 nbits)
{
	char Name[ SYS_MAX_OBJ_NAME_LEN + 1 ]; strcpy( Name, SYSDVC( Name ));
	if ( nbits == 0 )
		utilScreenDumpData( DUMP_SEL_SINK, DUMP_FMT_HEX,       DUMP_TYP_MEMORY, addr, 4,     0, ( dmpArg ) pd, strcat( Name , "Rx Buffer:" ));
	else if(( nbits >= 8 ) && ( nbits <= 16 ))
		utilScreenDumpData( DUMP_SEL_SINK, DUMP_FMT_IQ_SIGNED, DUMP_TYP_MEMORY, addr, nbits, 0, ( dmpArg ) pd, strcat( Name , "Rx Buffer:" ));
	else
		cprintf( "Format must be 0 (hex) or 8-16 (decimal Q:I)\n" );
}

// Receive stream, set block size to 0 to run continuously
static void RxStream( AXISink *pd, u32 RxAddr, u32 BlockSize )
{
	DVCPAR( RxAddr    )  = RxAddr;
	DVCPAR( BlockSize )  = BlockSize;
	SYSWRSLVREG( pd,  AXISINK_RX_ADDR_REG,         RxAddr );
	SYSWRSLVREG( pd,  AXISINK_BLOCK_SIZE_REG,   BlockSize );
	pd->WriteEnable( pd );
	if ( BlockSize == 0 ) {
		pd->StartContin( pd );
	} else {
		pd->Start1Shot( pd );
	}
}

static u32 ReadBR( u32 addr )
{
	static AXISink *pd = NULL;
	if ( addr & 0x80000000 ) { pd = ( AXISink * ) ( addr & 0x7FFFFFFF ); }
	else if ( pd != NULL ) { return ( pd->ReadBram( pd, ( u32 )addr )); }
	return( 0 );
}

static void WriteBR( u32 addr, u32 data )
{
	static AXISink *pd = NULL;
	if ( addr & 0x80000000 ) { pd = ( AXISink * )( addr & 0x7FFFFFFF ); }
	else if ( pd != NULL ) { pd->WriteBram( pd, ( u32 )addr, data ); }
}

static u32 TestBR( AXISink *pd )
{
	char Name[ 2 * SYS_MAX_OBJ_NAME_LEN + 1]; strcpy( Name, SYSDVC( Name ) );
	return ( UtilTestBram( ReadBR, WriteBR, pd->HwBramSize( pd ) / 4, 0xFFFFFFFF, Name ));
}

static u32 TestBram( AXISink *pd )
{
    ( void ) ReadBR( ( u32 ) pd | 0x80000000 );
    WriteBR( ( u32 ) pd | 0x80000000, 0 );
    ( void ) TestBR( pd );
	return( 0 );
}

