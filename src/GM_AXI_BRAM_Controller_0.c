/*
 * GM_AXI_BRAM_Controller_0.c
 *
 *  Created on: 20 May 2017
 *      Author: JULIAN MORTIMER
 */
#include "include.h"


// =================================
// Boilerplate Function Declarations
// =================================
void AXIBramCInit(
    AXIBramC *pd,
    const char      *Name,
    u32             BaseAddr,
    SysCtrl         *pSysCtrl,
    SysDvce         *pSysDvce,
    SysDesc         *pSysDesc,
    SysPtrs         *pSysPtrs,
    SysData         *pSysData,
    AXIBramCDesc *pDvcDesc,
    AXIBramCData *pDvcData
    );
static void SetParams    ( AXIBramC *pd );
static void SysDvceDump  ( AXIBramC *pd );
static void SysDescDump  ( AXIBramC *pd );
static void SysDataDump  ( AXIBramC *pd );
static void DvcDataInit  ( AXIBramC *pd );
static void DvcDataUpd   ( AXIBramC *pd );
static void DvcDescDump  ( AXIBramC *pd );
static void DvcDataDump  ( AXIBramC *pd );

// =================================================
// Device Structure Member Functions (public access)
// =================================================
static u32         TestBram      ( AXIBramC *pd );
static u32         ReadBram      ( AXIBramC *pd, u32 addr );
static void        WriteBram     ( AXIBramC *pd, u32 addr, u32 data );
static void        InitBram      ( AXIBramC *pd, u32 addr, u32 count, u32 data );
static void        LoadBram      ( AXIBramC *pd );
static inline u32  HwBramSize    ( AXIBramC *pd ) { return( SYSRDSLVREG( pd, AXIBRAMC_BRAM_SIZE_REG )); }
static void        Dump          ( AXIBramC *pd );

// Declare init structure here with member functions
AXIBramC AXIBramCInitStr = {
		NULL,
		NULL,
		NULL,
		TestBram,
		ReadBram,
		WriteBram,
		InitBram,
		LoadBram,
		HwBramSize,
		Dump
};

//===================
// Initialize Device
//===================
void AXIBramCInit( AXIBramC *pd, const char *Name, u32 BaseAddr, SysCtrl *pSysCtrl, SysDvce *pSysDvce, SysDesc *pSysDesc, SysPtrs *pSysPtrs, SysData *pSysData, AXIBramCDesc *pDvcDesc, AXIBramCData *pDvcData )
{
    if ( BaseAddr == 0 ) return;
    *pd             = AXIBramCInitStr;
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
        AXIBRAMC_HW_INFO_REG,
        AXIBRAMC_CTRL_REG,
        AXIBRAMC_INFO_REG,
        AXIBRAMC_STAT_REG,
        AXIBRAMC_RESET,
        0,
        0,
        0,
        0,
        0,
        0,
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
    SetParams( pd );
}

//========================
// Dump Device to Console
//========================
static void Dump( AXIBramC *pd ) {
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
static void SetParams( AXIBramC *pd )
{
// Set PL devive slave regs here
// Example SYSWRSLVREG( pd, AXIBRAMC_BOB_REG,    SYSPAR( Bob   ));
// Example SYSWRSLVREG( pd, AXIBRAMC_ALICE_REG , SYSPAR( Alice ));
}

static void DvcDataInit( AXIBramC *pd )
{
// Initialise DvcData members here
// Example: DVCDAT( Bob ) = 0;
}

static void DvcDataUpd( AXIBramC *pd )
{
// Update DvcData members here
// Example:	DVCDAT( Bob ) = SYSRDSLVREG( pd, AXIBRAMC_BOB_REG );
}

static void DvcDataDump( AXIBramC *pd )
{
	DvcDataUpd( pd );
// Insert device-specific SysDesc members here
// Example: cprintf("\nBob  : %u", toUint( DVCDAT( Bob )) );
}

static void DvcDescDump( AXIBramC *pd )
{
// Insert device-specific AXIBramCDesc members here
// Example: cprintf("\nAlice: %u", toUint( DVCPAR( Alice )) );
}

static void SysDvceDump( AXIBramC *pd )
{
	SYSDVCEDUMP( pd );
}

static void SysDataDump( AXIBramC *pd )
{
	SYSDATADUMP( pd );
}

static void SysDescDump( AXIBramC *pd )
{
// Insert device-specific SysDesc members here
// Example: cprintf("\nBob  : %u", toUint( SYSPAR( Bob )) );
}


// Insert device-specific code here
//======================================
// Device-specific non member functions
//======================================
static void WrData ( AXIBramC *pd, u32 addr );
static void EndWr  ( AXIBramC *pd );
static void RdData ( AXIBramC *pd, u32 addr );
static void EndRd  ( AXIBramC *pd );

//=================================
// Device-specific member functions
//=================================

static void WrData( AXIBramC *pd, u32 addr )
{
   hwWrReg( SYSDVC( BaseAddr ), AXIBRAMC_WR_ADDR_REG << 2, addr );
   hwWriteCtrlReg( SYSDVC( BaseAddr ), AXIBRAMC_CTRL_REG << 2, AXIBRAMC_WR_DATA );
}

static void EndWr( AXIBramC *pd )
{
   hwWriteCtrlReg( SYSDVC( BaseAddr ), AXIBRAMC_CTRL_REG << 2, AXIBRAMC_END_WR );
}

static void RdData( AXIBramC *pd, u32 addr )
{
   hwWrReg( SYSDVC( BaseAddr ), AXIBRAMC_RD_ADDR_REG << 2, addr );
   hwWriteCtrlReg( SYSDVC( BaseAddr ), AXIBRAMC_CTRL_REG << 2, AXIBRAMC_RD_DATA );
}

static void EndRd( AXIBramC *pd )
{
   hwWriteCtrlReg( SYSDVC( BaseAddr ), AXIBRAMC_CTRL_REG << 2, AXIBRAMC_END_RD );
}

static u32 ReadBram( AXIBramC *pd, u32 addr )
{
   RdData( pd, addr );
   u32 u = hwRdReg( SYSDVC( BaseAddr ), AXIBRAMC_RD_DATA_REG << 2 );
   EndRd( pd );
   return u;
}

static void WriteBram( AXIBramC *pd, u32 addr, u32 data )
{
	WrData( pd, addr );
	hwWrReg( SYSDVC( BaseAddr ), AXIBRAMC_WR_DATA_REG << 2, data );
	EndWr( pd );
}

static void InitBram( AXIBramC *pd, u32 addr, u32 count, u32 data )
{
	WrData( pd, addr );
	for ( int i = 0; i < count; i++ )
	{
      hwWrReg( SYSDVC( BaseAddr ), AXIBRAMC_WR_DATA_REG << 2, data );
	}
	EndWr( pd );
}

static void LoadBram( AXIBramC *pd )
{
	WrData( pd, 0 );
	u32 *pBuf = SYSPTR( pTxSourceBuf );
	for ( int i = 0; i < DVCPAR( BramSize ); i++ )
	{
      hwWrReg( SYSDVC( BaseAddr ), AXIBRAMC_WR_DATA_REG << 2, *pBuf++ );
	}
	EndWr( pd );
}

static u32 ReadBR( u32 addr )
{
	static AXIBramC *pd = NULL;
	if ( addr & 0x80000000 ) { pd = ( AXIBramC * ) ( addr & 0x7FFFFFFF ); }
	else if ( pd != NULL ) { return ( pd->ReadBram( pd, ( u32 )addr )); }
	return( 0 );
}

static void WriteBR( u32 addr, u32 data )
{
	static AXIBramC *pd = NULL;
	if ( addr & 0x80000000 ) { pd = ( AXIBramC * )( addr & 0x7FFFFFFF ); }
	else if ( pd != NULL ) { pd->WriteBram( pd, ( u32 )addr, data ); }
}

static u32 TestBR( AXIBramC *pd )
{
	char Name[ 2 * SYS_MAX_OBJ_NAME_LEN + 1]; strcpy( Name, SYSDVC( Name ) );
	return ( UtilTestBram( ReadBR, WriteBR, pd->HwBramSize( pd ) / 4, 0xFFFFFFFF, Name ));
}

static u32 TestBram( AXIBramC *pd )
{
    ( void ) ReadBR( ( u32 ) pd | 0x80000000 );
    WriteBR( ( u32 ) pd | 0x80000000, 0 );
    ( void ) TestBR( pd );
	return( 0 );
}







