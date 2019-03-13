/*
 *  Created on: 16 May 2017
 *      Author: JULIAN MORTIMER
 */

#ifndef EXTERN_H_
#define EXTERN_H_
//===============
// String tables
//===============
extern const char *RFFstrs [ 3 ];
extern const char *Modstrs [ 4 ];

//==============
// System Timer
//==============
extern TimerAPI *SysTimerAPI;
extern TimerAPI *SysTimer;
extern u32      SysTimerId;

extern void ChkNSDRxFreq( void );
extern void ChkTSITxFreq( void );

extern void WaitSec( double tsec );


//============
// Modulation
//============
extern Modem      *GetModem( ModType m );
extern const char *GetModemName( ModType m );
extern Modem UtilModem;
extern Modem TxModem;
extern Modem RxModem;

extern double TxRxFreqs[ 3 ];

extern void       DDUtilCreateDemod( ModType modType );
extern void       DDUtilDumpLUT( u32 addr, u32 fmt );
extern u32        *DDUtilGetLUTPtr( void );
extern DemodData  DDUtilDemodulate( Cplx32 din );
extern u32        DDUtilRawLUT  ( u32 addr );
extern u32        DDUtilDmodLUT ( u32 addr );
//=====
// PRS
//=====

extern u32 *PRSymbPtr         ( u32 nFFT );
extern u32 *CopyPRSymb        ( u32 nFFT, u32 *pBuf, double sf );
extern u32 *ToneConstSymbPtr  ( u32 nFFT );
extern u32 *CopyToneConstSymb ( u32 nFFT, u32 *pBuf, double sf );

//=========
// Numeric
//=========
extern u32 u32bit( u32 u, u32 i );
extern void fft  (cplx dest[], cplx src[], int n );
extern void ifft (cplx dest[], cplx src[], int n );
extern u32 prng32s( PRNG32 *, u32 );
extern u32 prng32( void );
extern u64 prng64s( PRNG64 *, u32 );
extern u64 prng64( void );
extern u32 clog2u32( u32 );
extern u32 log2u32( u32 );
extern u32 HammingWeight( u32 w );
extern u32 HammingDist( u32 w1, u32 w2 );

extern void *memsetd( void *p, double v, u32 n );
extern void *memsetc( void *p, cplx   v, u32 n );


extern double MeanXd    ( double *px, u32 n );
extern double SigmaXd   ( double *px, u32 n );
extern double SigmaX2d  ( double *px, u32 n );
extern double SigmaXYd  ( double *px, double *py, u32 n );
extern double MeanXs    ( s32    *px, u32 n );
extern double SigmaXs   ( s32    *px, u32 n );
extern double SigmaX2s  ( s32    *px, u32 n );
extern double SigmaXYs  ( s32    *px, s32 *py, u32 n );
extern double MeanXu    ( u32    *px, u32 n );
extern double SigmaXu   ( u32    *px, u32 n );
extern double SigmaX2u  ( u32    *px, u32 n );
extern double SigmaXYu  ( u32    *px, u32 *py, u32 n );

extern u32  RMeantoRNoise( double rmean, u32 nBits );
extern s16  Sext16( s16 i, u16 nbits );
extern u32  RMeantoRNoise( double rmean, u32 nBits );
extern s16  Sext16( s16 i, u16 nbits );

//==============================================
// Cplx32 - complex held in 32 bit unsigned int
//==============================================
extern _Bool  Cplx32Eq     ( Cplx32 z1, Cplx32 z2 );
extern Cplx32 ScaleCplx32  ( Cplx32 z, double sf );
extern Cplx32 CplxtoCplx32 ( cplx z );
extern cplx   Cplx32toCplx ( Cplx32 z );
extern Cplx32 U32toCplx32  ( u32 u );
extern u32    Cplx32toU32  ( Cplx32 z );
extern u32    U32Sext      ( u32 z, u32 nbits );

//=========
// Console
//=========
// console flags
extern _Bool verboseMode;
extern _Bool escSeqDet;
// Console I/O
extern char GetChar( void );
extern char cCode( char *s, const CtrlChar cctbl[] );
extern char *cDisp( char c, const CtrlChar cctbl[] );
extern void testKB( void );
extern _Bool scandbl ( char buf [ ], double *res );
extern _Bool scanuns ( char buf [ ], unsigned int *res );
extern _Bool scanint ( char buf [ ], signed int *res );
extern _Bool scanstr ( char buf [ ], int treatAsSigned, unsigned int *res );
extern unsigned int satMulAddu ( unsigned int p1, unsigned int p2, unsigned int s );
extern signed   int satMulAdds ( signed int p1, signed int p2, signed int s );
extern char *strucase ( char *dest, char *srce );
extern char *strlcase ( char *dest, char *srce );
extern char pause( char msg[] );
extern void dispSysTimerRunning( void );
extern int  sysTimerRunning( void );
extern void utilScreenDumpData( dmpSelector sel, dmpFormat fmt, dmpType typ, dmpAddr addr, dmpArg arg1, dmpArg arg2, dmpArg arg3, dmpTitle title );
extern u32  utilGetDumpData( dmpSelector sel, u32 addr, dmpArg arg, dmpArg ptr );
extern void u32bin2str( char s[], u32 u, s32 nbytes, u32 spccnt );
extern void DumpCplx( cplx *pz, u32 nData );



//===========
// HW Access
//===========
extern u32  hwRdReg   ( u32 addr, u32 offs );
extern void hwWrReg   ( u32 addr, u32 offs, u32 data );
extern u64  hwRdReg64 ( u32 addr, u32 offs );
extern void hwWrReg64 ( u32 addr, u32 offs, u64 data );
extern u32  SysDvceSimHwReg ( u32 baseAddr, u32 reg, u32 data, _Bool bWrite );

extern void hwWriteCtrlReg( INTPTR, u32, u32 );

extern void hwWriteData( INTPTR ,
		u32, u32,
		u32, u32,
		u32, u32 );

extern void hwWriteData64( INTPTR ,
        u32, u32,
        u32, u32,
        u32, u64 );

extern u32 hwReadData( INTPTR ,
		u32, u32, u32,
		u32, u32,
		u32, u32,
		u32 );

extern u64 hwReadData64( INTPTR ,
        u32, u32, u32,
        u32, u32,
        u32, u32,
        u32 );


//================
// Analog Devices
//================
extern double TxRxFreqs[];

extern void ad9361_console( void );
extern int  ad9361_main( void );
extern int  ad9364_Init_EMM( enum dds_data_select, enum samp_freq_select, uint32_t );
extern u32  AD9364RxSampleFreq( void );


extern void clkWizReset    ( tClkWiz * );
extern void clkWizWaitLock ( tClkWiz * );
extern void clkWizConfig   ( tClkWiz * );
extern int  clkWizSetFreq  ( tClkWiz *, u32 );
extern int  clkWizInit     ( tClkWiz *, u32 , u32, u32, const char [] );

//==========
// Emulator
//==========
extern tClkWiz clkWizDac;

//======================================
// System PL Devices
//======================================
extern AXIBramC        TxSourceAXIBramC;
extern AXISource       TxAXISource;
extern FFTC            TxFFTC;
extern FFTC            RxFFTC;
extern FFTC            RxNsdFFTC;
extern FFTC            RxPrsFFTC;
extern FFTC            RxCirFFTC;
extern FFTC            RxBpaFFTC;
extern TSInj           TxTSInj;
extern TREmul          RFChanTREmul;
extern TREmulSw        RFChanTREmulSw;
extern FreqCtr         RFChanTREmulFreqCtr;
extern FreqCtr         AD9364ClkFreqCtr;
extern FreqCtr         AD9364DACFreqCtr;
extern FreqCtr         AD9364ADCFreqCtr;
extern NSDet           RxNSDet;
extern AXISink         RxTddAXISink;
extern AXISink         RxAsqAXISink;
extern AXISink         RxPrsAXISink;
extern AXISink         RxDtaAXISink;
extern PrsCIR          RxPrsCIR;
extern PrsBPA          RxPrsBPA;

//========================================
// PL System Interface
//========================================
extern u32 FFTIndex( u32 nFFT );
//================
// System Control
//================
extern SysCtrl EMMSysCtrl;
extern SysDesc EMMSysDesc;
extern SysPtrs EMMSysPtrs;

extern SysCtrl *pEMMSysCtrl;
extern SysDesc *pEMMSysDesc;
extern SysPtrs *pEMMSysPtrs;


extern cplx PrsBPATimeDom    [ SYS_MAX_CARRIERS ];
extern cplx PrsBPATimeDomCorr[ SYS_MAX_CARRIERS ];
extern cplx PrsBPAFreqDom    [ SYS_MAX_CARRIERS ];
extern cplx PrsBPAFreqDomCorr[ SYS_MAX_CARRIERS ];
extern double PrsBPAPRSPhase [ SYS_MAX_CARRIERS ];
extern double PrsBPADtaPhase [ SYS_MAX_CARRIERS ];



extern void SysSetSysParams( _Bool ShowHelp, _Bool LoadDefault, SysCtrl *ps,
        u32      nCPre,
        u32      nFFT,
        u32      nSymbFrame,
        u32      nFrameBlock,
        u32      nSymbPRS,
        u32      nSymbAGC,
        u32      nSymbNSD,
        u32      TxRxFreq,
        u32      BWPercent,
        ModType  ModType,
        SoftBits SoftBits
);

//==================================
// General Device Testing
//==================================
extern void SysSetTestParams( _Bool ShowHelp, _Bool LoadDefault, SysCtrl *ps,
        u32     TxSourceMode,
        u32     SelEmul,
        u32     NsdTstMode,
        u32     RxTddSinkMode,
        u32     RxAsqSinkMode,
        u32     RxAsqSinkBlkSize,
        u32     RxPrsSinkMode
);


extern void    InitPLSystem( SysCtrl *ps );
extern void    StopPLSystem( SysCtrl *ps );
extern void    SysRunTest( SysCtrl *ps, u32 selEmulChan, double rmean, u32 srcmode );

extern u32     UtilTestBram( u32 ( * Read )( u32 addr ), void ( * Write )( u32 addr, u32 data ), u32 Size, u32 BitMask, char *Name );
extern u32     UtilTestBram64( u64 ( * Read )( u32 addr ), void ( * Write )( u32 addr, u64 data ), u32 Size, u64 BitMask, char *Name );


// =================================
// General Device Init
// =================================
extern void    SysDvceInit( SysDvce *, const char *, SysCtrl *, SysDesc *, SysPtrs *, SysData *,
        u32 ( * )( u32, u32, u32, _Bool ),
        u32, u32, u32, u32,
        u32, u32, u32, u32,
        u32, u32, u32, u32,
        u32, u32, u32,
        u32, u32, u32	);

extern void    SysDataInit( SysDvce *, u32, u32, u32, u32, u32, u32, u32, u32, u32, u32 );

// =================================
// AXIBramC Device Init
// =================================
extern void    AXIBramCInit( AXIBramC *, const char *, u32, SysCtrl *, SysDvce *, SysDesc *, SysPtrs *, SysData *, AXIBramCDesc *, AXIBramCData * );

// =================================
// AXISource Device Init
// =================================
extern void    AXISourceInit( AXISource *, const char *, u32, SysCtrl *, SysDvce *, SysDesc *, SysPtrs *, SysData *, AXISourceDesc *, AXISourceData * );

// =================================
// FFTC Device Init
// =================================
extern void    FFTCInit( FFTC *, const char *, u32, SysCtrl *, SysDvce *, SysDesc *, SysPtrs *, SysData *, FFTCDesc *, FFTCData * );

// =================================
// TSInj Device Init
// =================================
extern void    TSInjInit( TSInj *, const char *, u32, SysCtrl *, SysDvce *, SysDesc *, SysPtrs *, SysData *, TSInjDesc *, TSInjData * );

// =================================
// TREmulSw Device Init
// =================================
extern void    TREmulSwInit( TREmulSw *, const char *, u32, SysCtrl *, SysDvce *, SysDesc *, SysPtrs *, SysData *, TREmulSwDesc *, TREmulSwData * );

// =================================
// TREmul Device Init
// =================================
extern void    TREmulInit( TREmul *, const char *, u32, SysCtrl *, SysDvce *, SysDesc *, SysPtrs *, SysData *, TREmulDesc *, TREmulData * );

// =================================
// FreqCtr Device Init
// =================================
extern void    FreqCtrInit( FreqCtr *, const char *, u32, SysCtrl *, SysDvce *, SysDesc *, SysPtrs *, SysData *, FreqCtrDesc *, FreqCtrData * );

// =================================
// NSDet Device Init
// =================================
extern void    NSDetInit( NSDet *, const char *, u32, SysCtrl *, SysDvce *, SysDesc *, SysPtrs *, SysData *, NSDetDesc *, NSDetData * );

// =================================
// AXISink Device Init
// =================================
extern void    AXISinkInit( AXISink *, const char *, u32, SysCtrl *, SysDvce *, SysDesc *, SysPtrs *, SysData *, AXISinkDesc *, AXISinkData * );

// =================================
// PrsCIR Device Init
// =================================
extern void    PrsCIRInit( PrsCIR *, const char *, u32, SysCtrl *, SysDvce *, SysDesc *, SysPtrs *, SysData *, PrsDesc *, PrsCIRData * );

// =================================
// PrsBPA Device Init
// =================================
extern void    PrsBPAInit( PrsBPA *, const char *, u32, SysCtrl *, SysDvce *, SysDesc *, SysPtrs *, SysData *, PrsDesc *, PrsBPAData * );

#endif /* EXTERN_H_ */
