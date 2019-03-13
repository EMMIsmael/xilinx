/*
 * SystemControl.c
 *
 *  Created on: 30 Jul 2018
 *      Author: JULIAN MORTIMER
 */

#include "include.h"

double TxRxFreqs[ 3 ] = { 26000000.0, 30720000.0, 61440000.0 };

static void HelpSysParams   ( SysDesc *pd );
static void ShowSysParams   ( SysDesc *pd );
static void CalcTxSourceData ( SysDesc *pdsc, AXISource *psrc, u32 *pbuf );
static void AXISourceStart ( AXISource *pd, u32 TxMode, u32 addr, u32 blksize );
static void TSInjStart ( TSInj *pd );
static void TREmulSwStart ( TREmulSw *pd, u32 SelEmul );
static void TREmulStart ( TREmul *pd );
static void NSDetStart ( NSDet *pd, u32 tstmode, double RMean );
static void AXISinkStart ( AXISink *pd, u32 rxmode, u32 addr, u32 blksize );
static void PrsCIRStart ( PrsCIR *pd );
static void PrsBPAStart ( PrsBPA *pd );

void InitPLSystem ( SysCtrl *ps );
u32 FFTIndex( u32 nFFT );


u32 FFTIndex( u32 nFFT )
{
    u32 idx  = 0;
    switch ( nFFT )
    {
    case 256:
        idx = 0;
        break;
    case 512:
        idx = 1;
        break;
    case 1024:
        idx = 2;
        break;
    case 2048:
        idx = 3;
        break;
    default:
        assert( false );
        break;
    }
    return ( idx );
}




//--------------------------------------------------------
// The following protects against erroneous buffer sizing
// Which can destroy device structures when buffers are
// initialised
//--------------------------------------------------------
typedef union protPLSysBuffers PLSysBuf;
typedef struct sPLSysBuffers PLSysBufProt;
struct sPLSysBuffers {
    u32  TxSourceBuf      [ TX_AXISOURCE_BRAM_SIZE ];
    u32  TxTSInjTsBuf     [ 4 * TSINJ_BUF_SIZE ];
    u32  RFChanTREmulLUT  [ TREMUL_LUT_SIZE ];
    cplx RxPRSTdCorDeltaT [ SYS_MAX_CARRIERS ];
    cplx RxPRSFdCorFT     [ SYS_MAX_CARRIERS ];
    cplx RxPRSTdCorFT     [ SYS_MAX_CARRIERS ];
};

union protPLSysBuffers {
    PLSysBufProt Bufs;
    u32 Space [ ONE_MEGAWORD ];
};

PLSysBuf SysBuf;

// End of buffer protection
//--------------------------------------------------------

AXIBramC      TxSourceAXIBramC;
AXISource     TxAXISource;
FFTC          TxFFTC;
FFTC          RxFFTC;
FFTC          RxNsdFFTC;
FFTC          RxPrsFFTC;
FFTC          RxCirFFTC;
FFTC          RxBpaFFTC;
TSInj         TxTSInj;
TREmul        RFChanTREmul;
TREmulSw      RFChanTREmulSw;
FreqCtr       RFChanTREmulFreqCtr;
FreqCtr       AD9364ClkFreqCtr;
FreqCtr       AD9364DACFreqCtr;
FreqCtr       AD9364ADCFreqCtr;
NSDet         RxNSDet;
AXISink       RxTddAXISink;
AXISink       RxPrsAXISink;
AXISink       RxDtaAXISink;
AXISink       RxNsdB00AXISink;
AXISink       RxNsdB01AXISink;
AXISink       RxNsdB15AXISink;
AXISink       RxNsdFFTAXISink;
AXISink       RxNsdCordicAXISink;
AXISink       RxNsdA00AXISink;
AXISink       RxNsdA01AXISink;
AXISink       RxNsdA15AXISink;
PrsCIR        RxPrsCIR;
PrsBPA        RxPrsBPA;

SysDvce       TxSourceAXIBramCSysDvce;
SysDvce       TxAXISourceSysDvce;
SysDvce       TxFFTCSysDvce;
SysDvce       RxFFTCSysDvce;
SysDvce       RxNsdFFTCSysDvce;
SysDvce       RxPrsFFTCSysDvce;
SysDvce       RxCirFFTCSysDvce;
SysDvce       RxBpaFFTCSysDvce;
SysDvce       TxTSInjSysDvce;
SysDvce       RFChanTREmulSysDvce;
SysDvce       RFChanTREmulSwSysDvce;
SysDvce       RFChanTREmulFreqCtrSysDvce;
SysDvce       AD9364ClkFreqCtrSysDvce;
SysDvce       AD9364DACFreqCtrSysDvce;
SysDvce       AD9364ADCFreqCtrSysDvce;
SysDvce       RxNSDetSysDvce;
SysDvce       RxTddAXISinkSysDvce;
SysDvce       RxPrsAXISinkSysDvce;
SysDvce       RxDtaAXISinkSysDvce;
SysDvce       RxNsdB00AXISinkSysDvce;
SysDvce       RxNsdB01AXISinkSysDvce;
SysDvce       RxNsdB15AXISinkSysDvce;
SysDvce       RxNsdFFTAXISinkSysDvce;
SysDvce       RxNsdCordicAXISinkSysDvce;
SysDvce       RxNsdA00AXISinkSysDvce;
SysDvce       RxNsdA01AXISinkSysDvce;
SysDvce       RxNsdA15AXISinkSysDvce;
SysDvce       RxPrsCIRSysDvce;
SysDvce       RxPrsBPASysDvce;

SysData       TxSourceAXIBramCSysData;
SysData       TxAXISourceSysData;
SysData       TxFFTCSysData;
SysData       RxFFTCSysData;
SysData       RxNsdFFTCSysData;
SysData       RxPrsFFTCSysData;
SysData       RxCirFFTCSysData;
SysData       RxBpaFFTCSysData;
SysData       TxTSInjSysData;
SysData       RFChanTREmulSysData;
SysData       RFChanTREmulSwSysData;
SysData       RFChanTREmulFreqCtrSysData;
SysData       AD9364ClkFreqCtrSysData;
SysData       AD9364DACFreqCtrSysData;
SysData       AD9364ADCFreqCtrSysData;
SysData       RxNSDetSysData;
SysData       RxTddAXISinkSysData;
SysData       RxPrsAXISinkSysData;
SysData       RxDtaAXISinkSysData;
SysData       RxNsdB00AXISinkSysData;
SysData       RxNsdB01AXISinkSysData;
SysData       RxNsdB15AXISinkSysData;
SysData       RxNsdFFTAXISinkSysData;
SysData       RxNsdCordicAXISinkSysData;
SysData       RxNsdA00AXISinkSysData;
SysData       RxNsdA01AXISinkSysData;
SysData       RxNsdA15AXISinkSysData;
SysData       RxPrsCIRSysData;
SysData       RxPrsBPASysData;

AXIBramCData  TxSourceAXIBramCData;
AXISourceData TxAXISourceData;
FFTCData      TxFFTCData;
FFTCData      RxFFTCData;
FFTCData      RxNsdFFTCData;
FFTCData      RxPrsFFTCData;
FFTCData      RxCirFFTCData;
FFTCData      RxBpaFFTCData;
TSInjData     TxTSInjData;
TREmulData    RFChanTREmulData;
TREmulSwData  RFChanTREmulSwData;
FreqCtrData   RFChanTREmulFreqCtrData;
FreqCtrData   AD9364ClkFreqCtrData;
FreqCtrData   AD9364DACFreqCtrData;
FreqCtrData   AD9364ADCFreqCtrData;
NSDetData     RxNSDetData;
AXISinkData   RxTddAXISinkData;
AXISinkData   RxPrsAXISinkData;
AXISinkData   RxDtaAXISinkData;
AXISinkData   RxNsdA00AXISinkData;
AXISinkData   RxNsdA01AXISinkData;
AXISinkData   RxNsdA15AXISinkData;
AXISinkData   RxNsdFFTAXISinkData;
AXISinkData   RxNsdCordicAXISinkData;
AXISinkData   RxNsdB00AXISinkData;
AXISinkData   RxNsdB01AXISinkData;
AXISinkData   RxNsdB15AXISinkData;
PrsCIRData    RxPrsCIRData;
PrsBPAData    RxPrsBPAData;


SysDesc RFTstSysDesc = {
        0,             // u32      nCPre;
        1024,          // u32      nFFT;
        0,             // u32      nSymbAGC;
        0,             // u32      nSymbNSD;
        0,             // u32      nSymbPRS;
        16,            // u32      nSymbData;
        10,            // u32      nFrameBlock;
        2,             // u32      TxRxFreq;
        75,            // u32      BWPercent;
        MOD_4QAM,      // ModType  ModType;
        FOUR_SOFT_BITS // SoftBits SoftBits;
        };

SysDesc DefSysDesc = {
        64,            // u32      nCPre;
        1024,           // u32      nFFT;
        3,             // u32      nSymbAGC;
        1,             // u32      nSymbNSD;
        1,             // u32      nSymbPRS;
        0,             // u32      nSymbData;
        10,            // u32      nFrameBlock;
        2,             // u32      TxRxFreq;
        75,            // u32      BWPercent;
        MOD_4QAM,      // ModType  ModType;
        FOUR_SOFT_BITS // SoftBits SoftBits;
        };

SysDesc EMMSysDesc = {
        64,            // u32      nCPre;
        1024,          // u32      nFFT;
        3,             // u32      nSymbAGC;
        1,             // u32      nSymbNSD;
        1,             // u32      nSymbPRS;
        50,             // u32      nSymbData;
        10,            // u32      nFrameBlock;
        2,             // u32      TxRxFreq;
        75,            // u32      BWPercent;
        MOD_4QAM,      // ModType  ModType;
        FOUR_SOFT_BITS // SoftBits SoftBits;
        };

SysPtrs EMMSysPtrs = {
        SysBuf.Bufs.TxSourceBuf,      // u32     *pTxSourceBuf;
        SysBuf.Bufs.TxTSInjTsBuf,     // u32     *pTsBuf;
        SysBuf.Bufs.RFChanTREmulLUT,  // u32     *pRFChanTREmulLUT;
        SysBuf.Bufs.RxPRSTdCorDeltaT,
        SysBuf.Bufs.RxPRSFdCorFT,     // cplx    *pRxPRSTdCorDeltaT;
        SysBuf.Bufs.RxPRSTdCorFT      // cplx    *pRxPRSTdCorDeltaT;
};



SysCtrl EMMSysCtrl;


SysDesc *pEMMSysDesc = &EMMSysDesc;
SysPtrs *pEMMSysPtrs = &EMMSysPtrs;
SysCtrl *pEMMSysCtrl = &EMMSysCtrl;


AXIBramCDesc TxSourceAXIBramCDesc = {
        TX_AXISOURCE_BRAM_SIZE
};

AXISourceDesc TxAXISourceDesc = {
        0,
        0,
        TX_AXISOURCE_BRAM_SIZE,
        0xDEADBEEF,
        &TxSourceAXIBramC
};

//AXIBramCDesc RxPrsAXIBramCDesc = {
//        RX_PRS_BRAM_SIZE
//};


FFTCDesc TxFFTCDesc = {
        false,                         // _Bool    nFFTFixed;
        0,                             // u32      nFFT;
        SYS_MIN_CARRIERS,              // u32      nFFTMin;
        SYS_MAX_CARRIERS,              // u32      nFFTMax;
        false,                         // _Bool    hasCP;
        0,                             // u32      nCPre;
        FFTC_INV,                      // u32      fwdInv;
        true,                          // _Bool    isScaled;
        { 0, 0, 0, 0, 0, 0, 0, 0 },    // u32      scaleSch[ FFTC_SCALE_CHAN ];
        21,                            // u32      cfgRegLen;
        24,                            // u32      hwCfgRegLen;
        0,                             // u32      cfgFFTLoReg;
        0,                             // u32      cfgFFTHiReg;
        true,                          // _Bool    idxBitReversed;
        true                           // _Bool    hasXkIdx;
        };

FFTCDesc RxFFTCDesc = {
        false,                         // _Bool    nFFTFixed;
        0,                             // u32      nFFT;
        SYS_MIN_CARRIERS,              // u32      nFFTMin;
        SYS_MAX_CARRIERS,              // u32      nFFTMax;
        false,                         // _Bool    hasCP;
        0,                             // u32      nCPre;
        FFTC_FWD,                      // u32      fwdInv;
        true,                          // _Bool    isScaled;
        { 0, 0, 0, 0, 0, 0, 0, 0 },    // u32      scaleSch[ FFTC_SCALE_CHAN ];
        21,                            // u32      cfgRegLen;
        24,                            // u32      hwCfgRegLen;
        0,                             // u32      cfgFFTLoReg;
        0,                             // u32      cfgFFTHiReg;
        true,                          // _Bool    idxBitReversed;
        true                           // _Bool    hasXkIdx;
        };

FFTCDesc RxNsdFFTCDesc = {
        true,                          // _Bool    nFFTFixed;
        NSDET_N_FFT,                   // u32      nFFT;
        0,                             // u32      nFFTMin;
        0,                             // u32      nFFTMax;
        false,                         // _Bool    hasCP;
        0,                             // u32      nCPre;
        FFTC_FWD,                      // u32      fwdInv;
        true,                          // _Bool    isScaled;
        { 0, 0, 0, 0, 0, 0, 0, 0 },    // u32      scaleSch[ FFTC_SCALE_CHAN ];
        5,                             // u32      cfgRegLen;
        8,                             // u32      hwCfgRegLen;
        0,                             // u32      cfgFFTLoReg;
        0,                             // u32      cfgFFTHiReg;
        true,                          // _Bool    idxBitReversed;
        false                          // _Bool    hasXkIdx;
        };

FFTCDesc RxPrsFFTCDesc = {
        false,                         // _Bool    nFFTFixed;
        0,                             // u32      nFFT;
        SYS_MIN_CARRIERS,              // u32      nFFTMin;
        SYS_MAX_CARRIERS,              // u32      nFFTMax;
        false,                         // _Bool    hasCP;
        0,                             // u32      nCPre;
        FFTC_FWD,                      // u32      fwdInv;
        true,                          // _Bool    isScaled;
        { 0, 0, 0, 0, 0, 0, 0, 0 },    // u32      scaleSch[ FFTC_SCALE_CHAN ];
        21,                            // u32      cfgRegLen;
        24,                            // u32      hwCfgRegLen;
        0,                             // u32      cfgFFTLoReg;
        0,                             // u32      cfgFFTHiReg;
        false,                         // _Bool    idxBitReversed;
        false                          // _Bool    hasXkIdx;
        };

FFTCDesc RxCirFFTCDesc = {
        false,                         // _Bool    nFFTFixed;
        0,                             // u32      nFFT;
        SYS_MIN_CARRIERS,              // u32      nFFTMin;
        SYS_MAX_CARRIERS,              // u32      nFFTMax;
        false,                         // _Bool    hasCP;
        0,                             // u32      nCPre;
        FFTC_INV,                      // u32      fwdInv;
        true,                          // _Bool    isScaled;
        { 0, 0, 0, 0, 0, 0, 0, 0 },    // u32      scaleSch[ FFTC_SCALE_CHAN ];
        21,                            // u32      cfgRegLen;
        24,                            // u32      hwCfgRegLen;
        0,                             // u32      cfgFFTLoReg;
        0,                             // u32      cfgFFTHiReg;
        false,                         // _Bool    idxBitReversed;
        false                          // _Bool    hasXkIdx;
        };

FFTCDesc RxBpaFFTCDesc = {
        false,                         // _Bool    nFFTFixed;
        0,                             // u32      nFFT;
        SYS_MIN_CARRIERS,              // u32      nFFTMin;
        SYS_MAX_CARRIERS,              // u32      nFFTMax;
        false,                         // _Bool    hasCP;
        0,                             // u32      nCPre;
        FFTC_FWD,                      // u32      fwdInv;
        true,                          // _Bool    isScaled;
        { 0, 0, 0, 0, 0, 0, 0, 0 },    // u32      scaleSch[ FFTC_SCALE_CHAN ];
        21,                            // u32      cfgRegLen;
        24,                            // u32      hwCfgRegLen;
        0,                             // u32      cfgFFTLoReg;
        0,                             // u32      cfgFFTHiReg;
        true,                          // _Bool    idxBitReversed;
        true                           // _Bool    hasXkIdx;
        };

TSInjDesc TxTSInjDesc = {
        TSINJ_TSBUF_SIZE, // u32 TsBufSize;
        TSINJ_RXBUF_SIZE, // u32 TxBufSize;
        0xBADBEEF5,       // u32 TsBufInitVal;
        0xDEADFACE,       // u32 RxBufInitVal;
        TSINJ_DOUT_LSB    // u32 TSInjDoutLsb;
        };

TREmulDesc RFChanTREmulDesc = {
        SysBuf.Bufs.RFChanTREmulLUT,
        TREMUL_LUT_SIZE,
        TREMUL_LUT_SIZE
};

TREmulSwDesc RFChanTREmulSwDesc = {
        0
};

FreqCtrDesc RFChanTREmulFreqCtrDesc = {
        0
};

FreqCtrDesc AD9364ClkFreqCtrDesc = {
        0
};

FreqCtrDesc AD9364DACFreqCtrDesc = {
        0
};

FreqCtrDesc AD9364ADCFreqCtrDesc = {
        0
};

NSDetDesc RxNSDetDesc = {
		0,                          // u32   RNoise;
		0,                          // u32   FrqIdxRef;
		0,                          // u32   NsdDetThr;
		0,                          // u32   ChanDataSel;
		NSDET_TM_OFF,               // u32   NsdTstMode;
		{
				4, 4, 1, 1,      // u32   NsdFrqBlk    [ NSDET_N_BLK_DET   ];
				4, 1, 4, 1,
				1, 4, 1, 4,
				1, 1, 4, 4
		},
		{                           // u32   SampleMapTbl [ NSDET_MAX_FREQ    ];
				0b0100000000000001,
				0b0010000000000010,
				0b0001000000000100,
				0b0000100000001000,
				0b0000010000010000,
				0b0000001000100000,
				0b0000000101000000
		},
		{ 0, 0, 0, 0, 0, 0, 0 },    // u32   SampleMapIdx [ NSDET_MAX_FREQ  ];
		{ 0, 0, 0, 0 },             // u32   SampleMap    [ NSDET_N_FREQ      ];
		{ 0, 0, 0, 0 },             // u32   SampleMapFrq [ NSDET_N_FREQ    ];
		{ 2, 2, 2, 2 }              // u32   TsFreq       [ TSINJ_N_TSYMB_MAX ];
};

PrsDesc RxPrsDesc = {
    PRS_CIR_BRAM_SIZE,           // u32    CIRBramSize;
    0,                           // s32    DeltaT;
    0,                           // s32    DeltaF;
    0.0,                         // double deltafp;
    0.0,                         // double deltatp;
    0.0,                         // double deltafpp;
    0.0,                         // double deltatpp;
    0.0,                         // double alpha;
    0.0,                         // double beta;
    0.0,                         // double theta0;
    0.0,                         // double lastTheta;
    0.0,                         // double dataAlpha;
    0.0,                         // double dataBeta;
    0.0,                         // double fineDeltaf;
    0.0,                         // double fineDeltat;
//    { 1.29, 1.34, 1.32, 1.29 },  // double KFF[ SYS_N_CARRIER_CFG ];
    { 1.00, 1.00, 1.00, 1.00 },  // double KFF[ SYS_N_CARRIER_CFG ];
    { 0.00, 0.00, 0.00, 0.00 },  // double KTF[ SYS_N_CARRIER_CFG ];
    { 1.00, 1.00, 1.00, 1.00 },  // double KTT[ SYS_N_CARRIER_CFG ];
    { 0.00, 0.00, 0.00, 0.00 }   // double KFT[ SYS_N_CARRIER_CFG ];
};


AXISinkDesc RxTddAXISinkDesc = {
        0,
        0,
        RX_AXISINK_TDD_BRAM_SIZE,
        0xDEADBA11
};

AXISinkDesc RxPrsAXISinkDesc = {
        0,
        0,
        RX_AXISINK_PRS_BRAM_SIZE,
        0xDEADBA11
};

AXISinkDesc RxDtaAXISinkDesc = {
        0,
        0,
        RX_AXISINK_DTA_BRAM_SIZE,
        0xBAD1DEA5
};

AXISinkDesc RxNsdA00AXISinkDesc = {
        0,
        0,
        RX_AXISINK_ABS00_BRAM_SIZE,
        0xAB551000
};

AXISinkDesc RxNsdA01AXISinkDesc = {
        0,
        0,
        RX_AXISINK_ABS01_BRAM_SIZE,
        0xAB551001
};

AXISinkDesc RxNsdA15AXISinkDesc = {
        0,
        0,
        RX_AXISINK_ABS15_BRAM_SIZE,
        0xAB551015
};

AXISinkDesc RxNsdFFTAXISinkDesc = {
        0,
        0,
        RX_AXISINK_NSD_FFT_BRAM_SIZE,
        0xBAFF1ED
};

AXISinkDesc RxNsdCordicAXISinkDesc = {
        0,
        0,
        RX_AXISINK_NSD_CORDIC_BRAM_SIZE,
        0x2C0FFEE5
};

AXISinkDesc RxNsdB00AXISinkDesc = {
        0,
        0,
        RX_AXISINK_NSD_BCAST00_BRAM_SIZE,
        0xBCA57000
};

AXISinkDesc RxNsdB01AXISinkDesc = {
        0,
        0,
        RX_AXISINK_NSD_BCAST01_BRAM_SIZE,
        0xBCA57001
};

AXISinkDesc RxNsdB15AXISinkDesc = {
        0,
        0,
        RX_AXISINK_NSD_BCAST15_BRAM_SIZE,
        0xBCA57015
};

const char TxSourceAXIBramCName    [ SYS_MAX_OBJ_NAME_LEN + 1 ] = "Tx Data Source BRAM Controller";
const char TxAXISourceName         [ SYS_MAX_OBJ_NAME_LEN + 1 ] = "Tx Data Source";
const char TxFFTCName              [ SYS_MAX_OBJ_NAME_LEN + 1 ] = "Tx FFT Controller";
const char RxFFTCName              [ SYS_MAX_OBJ_NAME_LEN + 1 ] = "Rx FFT Controller";
const char RxNsdFFTCName           [ SYS_MAX_OBJ_NAME_LEN + 1 ] = "Rx NSD FFT Controller";
const char RxPrsFFTCName           [ SYS_MAX_OBJ_NAME_LEN + 1 ] = "Rx PRS FFT Controller";
const char RxCirFFTCName           [ SYS_MAX_OBJ_NAME_LEN + 1 ] = "Rx CIR FFT Controller";
const char RxBpaFFTCName           [ SYS_MAX_OBJ_NAME_LEN + 1 ] = "Rx BPA FFT Controller";
const char TxTSInjName             [ SYS_MAX_OBJ_NAME_LEN + 1 ] = "Tx Training Symbol Injector";
const char RFChanTREmulName        [ SYS_MAX_OBJ_NAME_LEN + 1 ] = "Tx Rx Emulator";
const char RFChanTREmulSwName      [ SYS_MAX_OBJ_NAME_LEN + 1 ] = "Tx Rx Emulator Switch";
const char RFChanTREmulFreqCtrName [ SYS_MAX_OBJ_NAME_LEN + 1 ] = "Tx Rx Emulator Frequency Counter";
const char AD9364ClkFreqCtrName    [ SYS_MAX_OBJ_NAME_LEN + 1 ] = "AD9364 CLK Frequency Counter";
const char AD9364DACFreqCtrName    [ SYS_MAX_OBJ_NAME_LEN + 1 ] = "AD9364 DAC Frequency Counter";
const char AD9364ADCFreqCtrName    [ SYS_MAX_OBJ_NAME_LEN + 1 ] = "AD9364 ADc Frequency Counter";
const char RxNSDetName             [ SYS_MAX_OBJ_NAME_LEN + 1 ] = "Rx Null Symbol Detector";
const char RxTddAXISinkName        [ SYS_MAX_OBJ_NAME_LEN + 1 ] = "Rx Time Domain Data Sink";
const char RxPrsAXISinkName        [ SYS_MAX_OBJ_NAME_LEN + 1 ] = "Rx PRS Data Sink";
const char RxDtaAXISinkName        [ SYS_MAX_OBJ_NAME_LEN + 1 ] = "Rx DTA Data Sink";
const char RxNsdB00AXISinkName     [ SYS_MAX_OBJ_NAME_LEN + 1 ] = "Rx NSD Bcast 00 Data Sink";
const char RxNsdB01AXISinkName     [ SYS_MAX_OBJ_NAME_LEN + 1 ] = "Rx NSD Bcast 01 Data Sink";
const char RxNsdB15AXISinkName     [ SYS_MAX_OBJ_NAME_LEN + 1 ] = "Rx NSD Bcast 15 Data Sink";
const char RxNsdFFTAXISinkName     [ SYS_MAX_OBJ_NAME_LEN + 1 ] = "Rx NSD FFT Data Sink";
const char RxNsdCordicAXISinkName  [ SYS_MAX_OBJ_NAME_LEN + 1 ] = "Rx NSD CORDIC Data Sink";
const char RxNsdA00AXISinkName     [ SYS_MAX_OBJ_NAME_LEN + 1 ] = "Rx NSD Abs 00 Data Sink";
const char RxNsdA01AXISinkName     [ SYS_MAX_OBJ_NAME_LEN + 1 ] = "Rx NSD Abs 01 Data Sink";
const char RxNsdA15AXISinkName     [ SYS_MAX_OBJ_NAME_LEN + 1 ] = "Rx NSD Abs 15 Data Sink";
const char RxPrsCIRName            [ SYS_MAX_OBJ_NAME_LEN + 1 ] = "Rx PRS Channel Impulse Response Calculator";
const char RxPrsBPAName            [ SYS_MAX_OBJ_NAME_LEN + 1 ] = "Rx PRS Brutal PRS Attach Mechanism";

SysCtrl SysCtrlInitStr = {
        &EMMSysDesc,           // SysDesc        *pSysDesc;
        &EMMSysPtrs,           // SysPtrs        *pSysPtrs;
        &TxSourceAXIBramC,     // AXIBramC       *pTxSourceAXIBramC;
        &TxAXISource,          // AXISource      *pTxAXISource;
        &TxFFTC,               // FFTC           *pTxFFTC;
        &RxFFTC,               // FFTC           *pRxFFTC;
        &RxNsdFFTC,            // FFTC           *pRxNsdFFTC;
        &RxPrsFFTC,            // FFTC           *pRxPrsFFTC;
        &RxCirFFTC,            // FFTC           *pRxCirFFTC;
        &RxBpaFFTC,            // FFTC           *pRxBpaFFTC;
        &TxTSInj,              // TSInj          *pTxTSInj;
        &RFChanTREmul,         // TREmul         *pRFChanTREmul;
        &RFChanTREmulSw,       // TREmulSw       *pRFChanTREmulSw;
        &RFChanTREmulFreqCtr,  // FreqCtr        *pRFChanTREmulFreqCtr;
        &AD9364ClkFreqCtr,     // FreqCtr        *pAD9364ClkFreqCtr;
        &AD9364DACFreqCtr,     // FreqCtr        *pAD9364DACFreqCtr;
        &AD9364ADCFreqCtr,     // FreqCtr        *pAD9364ADCFreqCtr;
        &RxNSDet,              // NSDet          *pRxNSDet;
        &RxTddAXISink,         // AXISink        *pRxTddAXISink;
        &RxPrsAXISink,         // AXISink        *pRxPrsAXISink;
        &RxDtaAXISink,         // AXISink        *pRxPrsAXISink;
        &RxNsdB00AXISink,      // AXISink        *pRxNsdB00AXISink;
        &RxNsdB01AXISink,      // AXISink        *pRxNsdB01AXISink;
        &RxNsdB15AXISink,      // AXISink        *pRxNsdB15AXISink;
        &RxNsdFFTAXISink,      // AXISink        *pRxNsdFFTAXISink;
        &RxNsdCordicAXISink,   // AXISink        *pRxNsdCordicAXISink;
        &RxNsdA00AXISink,      // AXISink        *pRxNsdA00AXISink;
        &RxNsdA01AXISink,      // AXISink        *pRxNsdA01AXISink;
        &RxNsdA15AXISink,      // AXISink        *pRxNsdA15AXISink;
        &RxPrsCIR,             // PrsCIR         *pRxPrsCIR;
        &RxPrsBPA,             // PrsBPA         *pRxPrsBPA;
        };

void InitPLSystem ( SysCtrl *ps )
{
    *ps = SysCtrlInitStr;
    SysDesc *pd = ps->pSysDesc;
    SysPtrs *pp = ps->pSysPtrs;
//    SysTest *pt = ps->pSysTest;
//    *pd = DefSysDesc;
//    *pt = DefSysTest;

    TxFFTCDesc.nFFT    = pd->nFFT;
    RxFFTCDesc.nFFT    = pd->nFFT;
    RxPrsFFTCDesc.nFFT = pd->nFFT;
    RxCirFFTCDesc.nFFT = pd->nFFT;
    RxBpaFFTCDesc.nFFT = pd->nFFT;

//                           Device Id   Device Name              Device Descriptor Name    Name Used in Macros
//                           ---------   -----------              ----------------------    -------------------
    DEVICE_INIT( ps, pd, pp, AXIBramC  , TxSourceAXIBramC       , TxSourceAXIBramC        , TX_SOURCE_AXIBRAMC     );
    DEVICE_INIT( ps, pd, pp, AXISource , TxAXISource            , TxAXISource             , TX_AXISOURCE           );
    DEVICE_INIT( ps, pd, pp, FFTC      , TxFFTC                 , TxFFTC                  , TX_FFTC                );
    DEVICE_INIT( ps, pd, pp, FFTC      , RxFFTC                 , RxFFTC                  , RX_FFTC                );
    DEVICE_INIT( ps, pd, pp, FFTC      , RxPrsFFTC              , RxPrsFFTC               , RX_PRS_FFTC            );
    DEVICE_INIT( ps, pd, pp, FFTC      , RxCirFFTC              , RxCirFFTC               , RX_CIR_FFTC            );
    DEVICE_INIT( ps, pd, pp, FFTC      , RxBpaFFTC              , RxBpaFFTC               , RX_BPA_FFTC            );
    DEVICE_INIT( ps, pd, pp, NSDet     , RxNSDet                , RxNSDet                 , RX_NSDET               );
    DEVICE_INIT( ps, pd, pp, TSInj     , TxTSInj                , TxTSInj                 , TX_TSINJ               );
    DEVICE_INIT( ps, pd, pp, TREmulSw  , RFChanTREmulSw         , RFChanTREmulSw          , RFCHAN_TREMULSW        );
    DEVICE_INIT( ps, pd, pp, TREmul    , RFChanTREmul           , RFChanTREmul            , RFCHAN_TREMUL          );
    DEVICE_INIT( ps, pd, pp, FreqCtr   , RFChanTREmulFreqCtr    , RFChanTREmulFreqCtr     , RFCHAN_TREMUL_FREQCTR  );
    DEVICE_INIT( ps, pd, pp, FreqCtr   , AD9364ClkFreqCtr       , AD9364ClkFreqCtr        , AD9364_CLK_FREQCTR     );
    DEVICE_INIT( ps, pd, pp, FreqCtr   , AD9364DACFreqCtr       , AD9364DACFreqCtr        , AD9364_DAC_FREQCTR     );
    DEVICE_INIT( ps, pd, pp, FreqCtr   , AD9364ADCFreqCtr       , AD9364ADCFreqCtr        , AD9364_ADC_FREQCTR     );
    DEVICE_INIT( ps, pd, pp, FFTC      , RxNsdFFTC              , RxNsdFFTC               , RX_NSD_FFTC            );
    DEVICE_INIT( ps, pd, pp, AXISink   , RxTddAXISink           , RxTddAXISink            , RX_TDD_AXISINK         );
    DEVICE_INIT( ps, pd, pp, AXISink   , RxPrsAXISink           , RxPrsAXISink            , RX_PRS_AXISINK         );
    DEVICE_INIT( ps, pd, pp, AXISink   , RxDtaAXISink           , RxDtaAXISink            , RX_DTA_AXISINK         );
    DEVICE_INIT( ps, pd, pp, AXISink   , RxNsdB00AXISink        , RxNsdB00AXISink         , RX_NSD_B00_AXISINK     );
    DEVICE_INIT( ps, pd, pp, AXISink   , RxNsdB01AXISink        , RxNsdB01AXISink         , RX_NSD_B01_AXISINK     );
    DEVICE_INIT( ps, pd, pp, AXISink   , RxNsdB15AXISink        , RxNsdB15AXISink         , RX_NSD_B15_AXISINK     );
    DEVICE_INIT( ps, pd, pp, AXISink   , RxNsdFFTAXISink        , RxNsdFFTAXISink         , RX_NSD_FFT_AXISINK     );
    DEVICE_INIT( ps, pd, pp, AXISink   , RxNsdCordicAXISink     , RxNsdCordicAXISink      , RX_NSD_CORDIC_AXISINK  );
    DEVICE_INIT( ps, pd, pp, AXISink   , RxNsdA00AXISink        , RxNsdA00AXISink         , RX_NSD_A00_AXISINK     );
    DEVICE_INIT( ps, pd, pp, AXISink   , RxNsdA01AXISink        , RxNsdA01AXISink         , RX_NSD_A01_AXISINK     );
    DEVICE_INIT( ps, pd, pp, AXISink   , RxNsdA15AXISink        , RxNsdA15AXISink         , RX_NSD_A15_AXISINK     );
    DEVICE_INIT( ps, pd, pp, PrsCIR    , RxPrsCIR               , RxPrs                   , RX_PRS_CIR             );
    DEVICE_INIT( ps, pd, pp, PrsBPA    , RxPrsBPA               , RxPrs                   , RX_PRS_BPA             );
#ifdef verbose
    cprintf ( "\n" );
#endif
}



const char RF2600 [ ] = "26.00 MHz";
const char RF3072 [ ] = "30.00 MHz";
const char RF6144 [ ] = "61.44 MHz";
const char *RFFstrs [ 3 ] = { RF2600, RF3072, RF6144 };

const char QPSKMod [ ] = "QPSK";
const char QAM4Mod [ ] = "4-QAM";
const char QM16Mod [ ] = "16-QAM";
const char QM64Mod [ ] = "64-QAM";
const char *Modstrs [ 4 ] = { QPSKMod, QAM4Mod, QM16Mod, QM64Mod };


// create a frame of data
static void CalcTxSourceData ( SysDesc *pdsc, AXISource *psrc, u32 *pbuf )
{
    u32 carrierMargin = pdsc->nFFT / 2 * ( 100 - pdsc->BWPercent ) / 100;
    u32 cLimLo = pdsc->nFFT / 2 - carrierMargin;
    u32 cLimHi = pdsc->nFFT / 2 + carrierMargin;

    u32 A_RND = 22695477;
    u32 C_RND = 1;
    u32 ACC = 0;
    u32 mod = 4;
    Modem *pMod = GetModem ( pdsc->ModType );
    if (pdsc->ModType == MOD_16QAM) {
        mod = 16;
    } else if (pdsc->ModType == MOD_64QAM) {
        mod = 64;
    }
    for ( int nPrsSymb = 0; nPrsSymb < pdsc->nSymbPRS; nPrsSymb++ ) {
        pbuf = CopyPRSymb ( pdsc->nFFT, pbuf, 1.0 / 4.0 );
    }
    for ( u32 nSymb = 0; nSymb < pdsc->nSymbData; nSymb++ ) {
       ACC = nSymb*nSymb;
       for ( u32 nCarrier = 0; nCarrier < pdsc->nFFT; nCarrier++ ) {
            if ( ( nCarrier <= cLimLo || nCarrier >= cLimHi ) && nCarrier > 0 ) {
                 ACC = ACC*A_RND + C_RND;
                 Cplx32 z = pMod->Modulate ( pMod, ( ACC>>16 ) % mod, 1 );
                 *pbuf++ = Cplx32toU32 ( z );
            }
            else
            {
                *pbuf++ = 0;
            }
        }
    }
}


static void AXISourceStart ( AXISource *pd, u32 TxMode, u32 addr, u32 blksize )
{
    pd->SetTxAddr ( pd, addr );
    pd->SetBlockSize ( pd, blksize );
    pd->LoadBram( pd );
    if ( TxMode == AXISOURCE_CONTIN )
        pd->StartContin ( pd );
    else
        pd->Start1Shot ( pd );
#ifdef verbose
    cprintf( "\n%s started.", SYSDVC( Name ) );
#endif
    }

static void TSInjStart ( TSInj *pd )
{
    pd->LoadTsBuf( pd );
    pd->Start ( pd );
#ifdef verbose
    cprintf( "\n%s started.", SYSDVC( Name ) );
#endif
}

static void TREmulSwStart ( TREmulSw *pd, u32 SelEmul )
{
    pd->SelChan ( pd, SelEmul );
#ifdef verbose
    cprintf( "\n%s started.", SYSDVC( Name ) );
#endif
}

static void TREmulStart ( TREmul *pd )
{
    pd->Start ( pd );
#ifdef verbose
    cprintf( "\n%s started.", SYSDVC( Name ) );
#endif
}

static void NSDetStart ( NSDet *pd, u32 tstmode, double rmean )
{
    pd->SetTstMode ( pd, tstmode );
    pd->SetRNoise( pd, rmean );
    pd->SetDbgTready( pd, true );
    pd->Start( pd );
#ifdef verbose
    cprintf( "\n%s started.", SYSDVC( Name ) );
#endif
}

static void AXISinkStart ( AXISink *pd, u32 rxmode, u32 addr, u32 blksize )
{
    pd->SetRxAddr ( pd, addr );
    pd->SetBlockSize ( pd, blksize );
    if ( rxmode == AXISINK_CONTIN )
        pd->StartContin ( pd );
    else
    pd->Start1Shot ( pd );
#ifdef verbose
    cprintf( "\n%s started.", SYSDVC( Name ) );
#endif
}

static void PrsCIRStart( PrsCIR *pd )
{

}

static void PrsBPAStart( PrsBPA *pd )
{

}


static void HelpSysParams ( SysDesc *pd )
{
    cprintf( "\nSysParam < nCPre nFFT nSymbData nFrameBlock nSymbPRS nSymbAGC nSymbNSD TxRxFreq BWPercent ModType >" );
    cprintf( "\n" );
    cprintf( "\nnCPre       : 0-( NFFT - 1 )" );
    cprintf( "\nnFFT        : { 256, 512, 1024, 2048 }" );
    cprintf( "\nnSymbData   : 1-MAXUINT" );
    cprintf( "\nnFrameBlock : 1-MAXUINT" );
    cprintf( "\n" );
    cprintf( "\nnSymbPRS    : 1-MAXUINT" );
    cprintf( "\nnSymbAGC    : 1-4" );
    cprintf( "\nnSymbNSD    : 1-4" );
    cprintf( "\nNote: nSymbAGC + nSymbNSD <= 4" );
    cprintf( "\n" );
    cprintf( "\nRF Freq     : 0 = 26.00 MHz, 1 = 30.72 MHz, 2 = 61.44 MHz" );
    cprintf( "\nBWPercent   : 75 = 75%% BW" );
    cprintf( "\nModType     : 0 = QPSK, 1 = 4QAM, 2 = 16QAM 3 = 64QAM" );
    cprintf( "\nSoftBits    : 0 = 3, 1 = 4" );
    cprintf( "\n" );
 }

void ShowSysParams ( SysDesc *pd )
{
    cprintf( "\nCurrent system params are:" );
    cprintf( "\n--------------------------" );
    cprintf( "\nCyclic Prefix            : %u", toUint( pd->nCPre ) );
    cprintf( "\nNFFT                     : %u", toUint( pd->nFFT ) );
    cprintf( "\nPRS Symbols              : %u", toUint( pd->nSymbPRS ) );
    cprintf( "\nAGC Symbols              : %u", toUint( pd->nSymbAGC ) );
    cprintf( "\nNSD Symbols              : %u", toUint( pd->nSymbNSD ) );
    cprintf( "\nData Symbols             : %u", toUint( pd->nSymbData ) );
    cprintf( "\nTrain Symbols / Frame    : %u", toUint( pd->nSymbAGC + pd->nSymbNSD  ) );
    cprintf( "\nPRS+Data Symbols / Frame : %u", toUint( pd->nSymbPRS + pd->nSymbData ) );
    cprintf( "\nFrames / Block           : %u", toUint( pd->nFrameBlock ) );
    cprintf( "\nRF Frequency             : %s", RFFstrs [ pd->TxRxFreq ] );
    cprintf( "\nUseful BW                : %u%%", toUint( pd->BWPercent ) );
    cprintf( "\nModulation               : %s", Modstrs [ pd->ModType ] );
    cprintf( "\nSoft bits (fixed)        : %u", pd->SoftBits );
    cprintf( "\n" );
}

void SysSetSysParams ( _Bool ShowHelp, _Bool LoadDefault, SysCtrl *ps,
        u32 nCPre,
        u32 nFFT,
        u32 nSymbData,
        u32 nFrameBlock,
        u32 nSymbPRS,
        u32 nSymbAGC,
        u32 nSymbNSD,
        u32 TxRxFreq,
        u32 BWPercent,
        ModType ModType,
        SoftBits SoftBits
        )
{
    SysDesc *pd = ps->pSysDesc;
    if ( ShowHelp )
    {
        HelpSysParams ( pd );
    }
    else if ( LoadDefault )
    {
        *pd = DefSysDesc;
    }
    else
    {
        pd->nFFT = ( ( nFFT == 2048 ) || ( nFFT == 1024 ) || ( nFFT = 512 ) || ( nFFT == 256 ) ) ? nFFT : 1024;
        pd->nCPre = min( nCPre, nFFT - 1 );
        pd->nSymbData = max( min( nSymbData, 1024 ), 1 );
        pd->nFrameBlock = nFrameBlock;
        pd->nSymbPRS = nSymbPRS;
        pd->nSymbAGC = min( nSymbAGC, SYS_N_TSYMB_MAX );
        pd->nSymbNSD = min( nSymbNSD, SYS_N_TSYMB_MAX - nSymbAGC );
        pd->TxRxFreq = ( TxRxFreq <= 2 ? TxRxFreq : 2 );
        pd->BWPercent = max( min( BWPercent, 100 ), 50 );
        pd->ModType = ( ModType <= 3 ? ModType : 1 );
        pd->SoftBits = SYS_N_SOFT_BITS;
    }
    ShowSysParams ( pd );
}

void SysRunTest( SysCtrl *ps, u32 selEmulChan, double rmean, u32 srcmode )
{
    SysDesc *pd = ps->pSysDesc;
    SysPtrs *pp = ps->pSysPtrs;
    u32 sel = ( selEmulChan == SEL_EMUL_WIRE ) ? SEL_EMUL_WIRE : SEL_EMUL_RF_CHAN;
#ifdef verbose
    if ( selEmulChan == SEL_EMUL_WIRE ) {
        cprintf( "\nWire channel selected." );
    } else {
        cprintf( "\nRF channel selected." );
    }
#endif
    //( void ) ad9364_Init_EMM( DATA_SEL_DMA, pd->TxRxFreq, 1 );
    CalcTxSourceData ( pd, ps->pTxAXISource, pp->pTxSourceBuf );
    InitPLSystem   ( ps );
    TREmulSwStart  ( ps->pRFChanTREmulSw, sel );
    TREmulStart    ( ps->pRFChanTREmul );
    NSDetStart     ( ps->pRxNSDet, NSDET_TM_OFF, rmean );
    WaitSec( 1.0 );
    AXISinkStart   ( ps->pRxTddAXISink,       AXISINK_CONTIN, 0, 0 );
    AXISinkStart   ( ps->pRxPrsAXISink,       AXISINK_CONTIN, 0, 0 );
    AXISinkStart   ( ps->pRxDtaAXISink,       AXISINK_CONTIN, 0, 0 );
    AXISinkStart   ( ps->pRxNsdB00AXISink,    AXISINK_CONTIN, 0, 0 );
    AXISinkStart   ( ps->pRxNsdB01AXISink,    AXISINK_CONTIN, 0, 0 );
    AXISinkStart   ( ps->pRxNsdB15AXISink,    AXISINK_CONTIN, 0, 0 );
    AXISinkStart   ( ps->pRxNsdCordicAXISink, AXISINK_CONTIN, 0, 0 );
    AXISinkStart   ( ps->pRxNsdFFTAXISink,    AXISINK_CONTIN, 0, 0 );
    AXISinkStart   ( ps->pRxNsdA00AXISink,    AXISINK_CONTIN, 0, 0 );
    AXISinkStart   ( ps->pRxNsdA01AXISink,    AXISINK_CONTIN, 0, 0 );
    AXISinkStart   ( ps->pRxNsdA15AXISink,    AXISINK_CONTIN, 0, 0 );
    PrsCIRStart    ( ps->pRxPrsCIR );
    PrsBPAStart    ( ps->pRxPrsBPA );
    if ( srcmode ) {
        AXISourceStart ( ps->pTxAXISource, AXISOURCE_CONTIN, 0, pd->nFFT * ( pd->nSymbData + pd->nSymbPRS ) );
    } else {
        AXISourceStart ( ps->pTxAXISource, AXISOURCE_ONESHOT, 0, pd->nFFT * ( pd->nSymbData + pd->nSymbPRS ) );
    }
    TSInjStart     ( ps->pTxTSInj );
#ifdef verbose
    cprintf( "\n\n" );
#endif
}


void StopPLSystem( SysCtrl *ps )
{
    ps->pTxAXISource->Stop( ps->pTxAXISource );
    WaitSec( 1.0 );
    ps->pTxTSInj->Stop            ( ps->pTxTSInj );
    ps->pRFChanTREmul->Stop       ( ps->pRFChanTREmul );
    ps->pRxNSDet->Stop            ( ps->pRxNSDet );
    ps->pRxTddAXISink->Stop       ( ps->pRxTddAXISink );
    ps->pRxPrsAXISink->Stop       ( ps->pRxPrsAXISink );
    ps->pRxDtaAXISink->Stop       ( ps->pRxDtaAXISink );
    ps->pRxNsdB00AXISink->Stop    ( ps->pRxNsdB00AXISink );
    ps->pRxNsdB01AXISink->Stop    ( ps->pRxNsdB01AXISink );
    ps->pRxNsdB15AXISink->Stop    ( ps->pRxNsdB15AXISink );
    ps->pRxNsdFFTAXISink->Stop    ( ps->pRxNsdFFTAXISink );
    ps->pRxNsdCordicAXISink->Stop ( ps->pRxNsdCordicAXISink );
    ps->pRxNsdA00AXISink->Stop    ( ps->pRxNsdA00AXISink );
    ps->pRxNsdA01AXISink->Stop    ( ps->pRxNsdA01AXISink );
    ps->pRxNsdA15AXISink->Stop    ( ps->pRxNsdA15AXISink );
#ifdef verbose
    cprintf( "Test halted.\n" );
#endif
}

void ChkTSITxFreq( void )
{
    TSInj *pd = &TxTSInj;
    u32 u0 = SYSRDSLVREG( pd, TSINJ_TX_SYMBOL_CTR_REG );
    double t0 = SysTimer->ReadTimerLongSec( SysTimer, SysTimerId );
    double t1 = t0;
    while ( t1 < t0 + 1.0 ) {
        t1 = SysTimer->ReadTimerLongSec( SysTimer, SysTimerId );
    }
    u32 u1 = SYSRDSLVREG( pd, TSINJ_TX_SYMBOL_CTR_REG );
    double r = ( double ) ( u1 - u0 ) / ( t1 - t0 );
    cprintf("\nData Rate = %.4g Mt/s", 1088.0 * r / 1e6 );
    cprintf("\n");

}

void ChkNSDRxFreq( void )
{
    NSDet *pd = &RxNSDet;
    u32 u0 = SYSRDSLVREG( pd, NSDET_RX_SYMB_CTR_REG );
    double t0 = SysTimer->ReadTimerLongSec( SysTimer, SysTimerId );
    double t1 = t0;
    while ( t1 < t0 + 1.0 ) {
        t1 = SysTimer->ReadTimerLongSec( SysTimer, SysTimerId );
    }
    u32 u1 = SYSRDSLVREG( pd, NSDET_RX_SYMB_CTR_REG );
    double r = ( double ) ( u1 - u0 ) / ( t1 - t0 );
    cprintf("\nData Rate = %.4g Mt/s", r * 1088.0 / 1e6 );
    cprintf("\n");

}
