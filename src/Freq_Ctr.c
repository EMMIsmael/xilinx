/*
 * Freq_Ctr.c
 *
 *  Created on: 13 Jun 2017
 *      Author: JULIAN MORTIMER
 */
#include "include.h"


// =================================
// Boilerplate Function Declarations
// =================================
void FreqCtrInit(
    FreqCtr *pd,
    const char      *Name,
    u32             BaseAddr,
    SysCtrl         *pSysCtrl,
    SysDvce         *pSysDvce,
    SysDesc         *pSysDesc,
    SysPtrs         *pSysPtrs,
    SysData         *pSysData,
    FreqCtrDesc *pDvcDesc,
    FreqCtrData *pDvcData
    );
static void SetParams    ( FreqCtr *pd );
static void SysDvceDump  ( FreqCtr *pd );
static void SysDescDump  ( FreqCtr *pd );
static void SysDataDump  ( FreqCtr *pd );
static void DvcDataInit  ( FreqCtr *pd );
static void DvcDataUpd   ( FreqCtr *pd );
static void DvcDescDump  ( FreqCtr *pd );
static void DvcDataDump  ( FreqCtr *pd );

// =================================================
// Device Structure Member Functions (public access)
// =================================================
static inline u32  ReadCtrReg        ( FreqCtr *pd ) { return ( SYSRDSLVREG( pd, FREQCTR_CTR_REG )); }
static inline u32  ReadTmrReg        ( FreqCtr *pd ) { return ( SYSRDSLVREG( pd, FREQCTR_TMR_REG )); }
static double      ReadFreqCtr       ( FreqCtr *pd );
static void        Dump              ( FreqCtr *pd );

// Declare init structure here with member functions
FreqCtr FreqCtrInitStr = {
    NULL,
	NULL,
	NULL,
	ReadCtrReg,
	ReadTmrReg,
	ReadFreqCtr,
	Dump
};

//===================
// Initialize Device
//===================
void FreqCtrInit( FreqCtr *pd, const char *Name, u32 BaseAddr, SysCtrl *pSysCtrl, SysDvce *pSysDvce, SysDesc *pSysDesc, SysPtrs *pSysPtrs, SysData *pSysData, FreqCtrDesc *pDvcDesc, FreqCtrData *pDvcData )
{
	*pd             = FreqCtrInitStr;
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
        FREQCTR_HW_INFO_REG,
        FREQCTR_CTRL_REG,
        FREQCTR_INFO_REG,
        FREQCTR_STAT_REG,
        FREQCTR_RESET,
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
    SetParams( pd );
}

//========================
// Dump Device to Console
//========================
static void Dump( FreqCtr *pd ) {
    printf("\n"); fflush( stdout );
	SysDvceDump( pd );
	SysDescDump( pd );
	DvcDescDump( pd );
	SysDataDump( pd );
	DvcDataDump( pd );
	printf("\n"); fflush( stdout );
}

// End of boilerplate code
//========================
// Device-specific code
//========================
static void SetParams( FreqCtr *pd )
{
// Set PL devive slave regs here
// Example SYSWRSLVREG( pd, FREQCTR_BOB_REG,    SYSPAR( Bob   ));
// Example SYSWRSLVREG( pd, FREQCTR_ALICE_REG , SYSPAR( Alice ));
}

static void DvcDataInit( FreqCtr *pd )
{
// Initialise DvcData members here
// Example: DVCDAT( Bob ) = 0;
}

static void DvcDataUpd( FreqCtr *pd )
{
// Update DvcData members here
// Example:	DVCDAT( Bob ) = SYSRDSLVREG( pd, FREQCTR_BOB_REG );
}

static void DvcDataDump( FreqCtr *pd )
{
	DvcDataUpd( pd );
// Insert device-specific SysDesc members here
// Example: printf("\nBob  : %u", toUint( DVCDAT( Bob )) );   fflush( stdout );
}

static void DvcDescDump( FreqCtr *pd )
{
// Insert device-specific FreqCtrDesc members here
// Example: printf("\nAlice: %u", toUint( DVCPAR( Alice )) ); fflush( stdout );
}

static void SysDvceDump( FreqCtr *pd )
{
	SYSDVCEDUMP( pd );
}

static void SysDataDump( FreqCtr *pd )
{
	SYSDATADUMP( pd );
}

static void SysDescDump( FreqCtr *pd )
{
// Insert device-specific SysDesc members here
// Example: printf("\nBob  : %u", toUint( SYSPAR( Bob )) );   fflush( stdout );
}


// Insert device-specific code here
//======================================
// Device-specific non member functions
//======================================
// Example: static inline void EnterAdmin    ( FreqCtr *pd ) { hwWriteCtrlReg( SYSDVC( BaseAddr ), FREQCTR_CTRL_REG << 2, FREQCTR_ENTER_ADMIN ); }

//=================================
// Device-specific member functions
//=================================

static double ReadFreqCtr( FreqCtr *pd )
{
    u32 u = ReadCtrReg( pd );
    return( ( double ) u / 1000000.0 );
}

