/*
 * FFT_Controller.c
 *
 *  Created on: 27 Jun 2017
 *      Author: JULIAN MORTIMER
 */
#include "include.h"


// =================================
// Boilerplate Function Declarations
// =================================
void FFTCInit(
    FFTC *pd,
    const char      *Name,
    u32             BaseAddr,
    SysCtrl         *pSysCtrl,
    SysDvce         *pSysDvce,
    SysDesc         *pSysDesc,
    SysPtrs         *pSysPtrs,
    SysData         *pSysData,
    FFTCDesc *pDvcDesc,
    FFTCData *pDvcData
    );
static void SetParams    ( FFTC *pd );
static void SysDvceDump  ( FFTC *pd );
static void SysDescDump  ( FFTC *pd );
static void SysDataDump  ( FFTC *pd );
static void DvcDataInit  ( FFTC *pd );
static void DvcDataUpd   ( FFTC *pd );
static void DvcDescDump  ( FFTC *pd );
static void DvcDataDump  ( FFTC *pd );

static void LocDvceDump ( FFTC *pd );

// =================================================
// Device Structure Member Functions (public access)
// =================================================
static inline u32  HwCfgRegLen   ( FFTC *pd ) { return ( SYSRDSLVREG( pd, FFTC_CFG_REG_LEN_REG ) ); }
static void        CfgFFT        ( FFTC *pd );
static void        Dump          ( FFTC *pd );

static void CreateCfgRegs( FFTC *pd );

// Declare init structure here with member functions
FFTC FFTCInitStr = {
			NULL,
			NULL,
			NULL,
			HwCfgRegLen,
			CfgFFT,
			Dump
};

//===================
// Initialize Device
//===================
void FFTCInit( FFTC *pd, const char *Name, u32 BaseAddr, SysCtrl *pSysCtrl, SysDvce *pSysDvce, SysDesc *pSysDesc, SysPtrs *pSysPtrs, SysData *pSysData, FFTCDesc *pDvcDesc, FFTCData *pDvcData )
{
	*pd             = FFTCInitStr;
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
        FFTC_HW_INFO_REG,
        FFTC_CTRL_REG,
        FFTC_INFO_REG,
        FFTC_STAT_REG,
        FFTC_RESET,
        FFTC_INT_ACK,
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
    CreateCfgRegs( pd );
    SetParams( pd );
}

//========================
// Dump Device to Console
//========================
static void Dump( FFTC *pd ) {
    cprintf("\n");
	SysDvceDump( pd );
	SysDescDump( pd );
	DvcDescDump( pd );
	SysDataDump( pd );
	DvcDataDump( pd );
	LocDvceDump( pd );
	cprintf("\n");
}



// End of boilerplate code
//========================
// Device-specific code
//========================

static void SetParams( FFTC *pd )
{
	SYSWRSLVREG( pd, FFTC_CFG_FFT_LO_REG, DVCPAR( cfgFFTLoReg ));
	SYSWRSLVREG( pd, FFTC_CFG_FFT_HI_REG, DVCPAR( cfgFFTHiReg ));
    pd->CfgFFT( pd );
}

static void DvcDataInit( FFTC *pd )
{
// Initialise DvcData members here
// Example: DVCDAT( Bob ) = 0;
}

static void DvcDataUpd( FFTC *pd )
{
// Update DvcData members here
// Example:	DVCDAT( Bob ) = SYSRDSLVREG( pd, FFTC_BOB_REG );
}
static void DvcDataDump( FFTC *pd )
{
	DvcDataUpd( pd );
// Insert device-specific SysDesc members here
// Example: cprintf("\nBob  : %u", toUint( DVCDAT( Bob )) );
}

static void DvcDescDump( FFTC *pd )
{
// Insert device-specific FFTCDesc members here
// Example: cprintf("\nAlice: %u", toUint( DVCPAR( Alice )) );
}

static void SysDvceDump( FFTC *pd )
{
	SYSDVCEDUMP( pd );
}

static void SysDataDump( FFTC *pd )
{
	SYSDATADUMP( pd );
}

static void SysDescDump( FFTC *pd )
{
// Insert device-specific SysDesc members here
// Example: cprintf("\nBob  : %u", toUint( SYSPAR( Bob )) );
}



// Insert device-specific code here
//======================================
// Device-specific non member functions
//======================================
// Example: static inline void EnterAdmin    ( FFTC *pd ) { hwWriteCtrlReg( SYSDVC( BaseAddr ), FFTC_CTRL_REG << 2, FFTC_ENTER_ADMIN ); }

//=================================
// Device-specific member functions
//=================================
static const char scale0[] = "1";
static const char scale1[] = "1/2";
static const char scale2[] = "1/4";
static const char scale3[] = "1/8";

static const char *pscale[ 4 ] = { scale0, scale1, scale2, scale3 };


static void LocDvceDump( FFTC *pd )
{
	char sh[ 48 ];
	char sl[ 48 ];
	for ( int i = 0; i < sizeof( sh ); i++ ) sh[ i ] = '\0';
	for ( int i = 0; i < sizeof( sh ); i++ ) sl[ i ] = '\0';
	cprintf( "\n" );
	if (! DVCPAR( nFFTFixed )) {
		cprintf( "FFT length           :  %u (%u - %u)\n", toUint( DVCPAR( nFFT )), toUint( DVCPAR( nFFTMin )), toUint( DVCPAR( nFFTMax )));
	}
	else {
		cprintf( "FFT length           :  %u (fixed)\n", toUint( DVCPAR( nFFT )));
	}
	if ( DVCPAR( hasCP )) {
		cprintf( "Cyclic prefix        :  %u\n", toUint( DVCPAR( nCPre )));
	}
	else {
		cprintf( "Cyclic prefix        :  No insertion\n" );
	}
	if ( DVCPAR( fwdInv ) == FFTC_FWD ) {
		cprintf( "Type                 :  Forward\n" );
	} else {
		cprintf( "Type                 :  Inverse\n" );
	}
	if ( DVCPAR( isScaled ))
	{
		u32 i;
	    cprintf( "Scaling\n" );
        for ( i = 0; i < ( log2u32( DVCPAR( nFFT )) - 1 ) / 2; i ++ )
        {
    	    cprintf( "Stages %2u -%2u        : %s\n", 2 * toUint( i ), 1 + 2 * toUint( i ), pscale[ DVCPAR( scaleSch[ i ] ) & 0x3 ] );
        }
        if ( IS_ODD( log2u32( DVCPAR( nFFT ))))
        {
    	    cprintf( "Stage  %2u            : %s\n", 2 * toUint( i ), pscale[ DVCPAR( scaleSch[ i ] ) & 0x3 ] );
        }
        else
        {
    	    cprintf( "Stages %2u -%2u        : %s\n", 2 * toUint( i ), 1 + 2 * toUint( i ), pscale[ DVCPAR( scaleSch[ i ] ) & 0x3 ] );
        }
    }
	else
	{
		cprintf( "Scaling              :  No scaling\n" );
	}
	u32bin2str( sh, DVCPAR( cfgFFTHiReg ), pd->HwCfgRegLen( pd ) / 8 - 4, 4 );
	u32bin2str( sl, DVCPAR( cfgFFTLoReg ), MIN( pd->HwCfgRegLen( pd ) / 8, 4 ), 4 );
	cprintf( "\n" );
	cprintf( "Cfg out (length %-2u)  : %s%s\n", toUint( pd->HwCfgRegLen( pd )), sh, sl );
	cprintf( "Cfg reg hi           : " );
	if ( strlen( sh ) == 0 ) {
		for ( int i = 0; i < pd->HwCfgRegLen( pd ) / 8; i++ ) cprintf( "0000 0000 " );
		cprintf( "\n" );
	}
	else
		cprintf( "%s\n", sh );

	cprintf( "Cfg reg lo           : " );
	if ( strlen( sl ) == 0 )
		cprintf( "0000 0000 0000 0000 0000 0000 0000 0000\n");
	else
		cprintf( "%s\n", sl );
	cprintf( "Cfg reg lgth from hw : %u\n", toUint( DVCPAR( hwCfgRegLen )));
	cprintf( "\n" );
}


static void CreateCfgRegs( FFTC *pd )
{
	u32 cfglo  = 0;
	u32 cfghi  = 0;
	u32 cbits  = 0;
//	u64 cfgrg  = 0;
	u32 i      = 0;

	u32 bptr   = 0;
	if (! DVCPAR( nFFTFixed )) {
		cfglo = log2u32( DVCPAR( nFFT ) );
		bptr += 8;
	}

	if ( DVCPAR( hasCP )) {
	   cbits = DVCPAR( nCPre ) << (( log2u32( DVCPAR( nFFTMax ))) - log2u32( DVCPAR( nFFT )));
       cfglo += cbits << 8;
       bptr  += 8 * ( log2u32( DVCPAR( nFFTMax )) > 8 ? 2 : 1 );
	}
	cfglo |= DVCPAR( fwdInv ) << bptr;
	bptr ++ ;

    if ( DVCPAR( isScaled )) {
    	for ( i = 1; i < ( log2u32( DVCPAR( nFFT ) - 1 ) / 2 ); i++ ) {
            DVCPAR( scaleSch[ i ] ) = 0x02;
        }
        DVCPAR( scaleSch[ 0 ]     ) = 0x03;
        DVCPAR( scaleSch[ i - 1 ] ) = 0x01;
        if ( DVCPAR( nFFT ) == 16 )
        {
            DVCPAR( scaleSch[ 0 ]     ) = 0x02;
            DVCPAR( scaleSch[ 1 ]     ) = 0x03;
        }
        else
        {
            DVCPAR( scaleSch[ i ]     ) = ( IS_ODD( log2u32( DVCPAR( nFFT ) )) ? 0x00 : 0x01 /* 0x01 : 0x02 */ );
        }
        cbits = 0;
	    for ( i = 0; i < ( log2u32( DVCPAR( nFFT ) ) + 1 ) / 2; i ++ )
    	    cbits |= DVCPAR( scaleSch[ i ] ) << ( 2 * i );

	    if ( bptr == 25 ) cfghi = cbits >> 7;
        cfglo |= cbits << bptr;
        bptr += ( log2u32( DVCPAR( nFFTMax ) ) + 1 ) & 0xFFFFFFFE;
    } else {
    	for ( int i = 0; i < FFTC_SCALE_CHAN; i++ ) {
    		DVCPAR( scaleSch[ i ] ) = 0;
    	}
    }
    DVCPAR( cfgRegLen )   = bptr;
    DVCPAR( cfgFFTLoReg ) = cfglo;
    DVCPAR( cfgFFTHiReg ) = cfghi;
    DVCPAR( hwCfgRegLen ) = pd->HwCfgRegLen( pd );
}


static void CfgFFT( FFTC *pd )
{
    hwWriteCtrlReg( SYSDVC( BaseAddr ), FFTC_CTRL_REG << 2, FFTC_TX_CFG_FFT );
	double t = SYS_TIMER( ReadTimerLongSec );
	while (( !( SYSRDSTATREG( pd ) & FFTC_TX_CFG_FFT_DONE )) && ( SYS_TIMER( ReadTimerLongSec ) - t < 0.0001 ))
		;
	if ( SYS_TIMER( ReadTimerLongSec ) - t > 0.0001 )
	{
		cprintf("\nTimed out while configuring %s\n", SYSDVC( Name ) );
	}
}


