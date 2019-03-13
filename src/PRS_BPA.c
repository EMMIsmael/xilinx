/*
 * PRS_BPA.c
 *
 *  Created on: 12 Sep 2018
 *      Author: JULIAN MORTIMER
 */
#include "include.h"

shell PrsBPAShell;

PrsPhi txPrsPhi;

extern cplx RxPRSTdCorFT [ SYS_MAX_CARRIERS ];
extern cplx RxPRSFdCordeltatp [ SYS_MAX_CARRIERS ];

static void ClearData ( PrsBPA *pd );
static double PrsPhaseError ( PrsBPA *pd, cplx *p, u32 i );
static PrsPhase PrsPhiPhase ( PrsBPA *pd, u32 n );
static void PrsPhiInit ( PrsBPA *pd );

void CMatDump ( cplx *m, u32 nCol, u32 nRow, char name [ ] );
void CRowDump ( cplx *m, u32 ncol, char name [ ] );
void CColDump ( cplx *m, u32 nRow, char name [ ] );
void DMatDump ( double *m, u32 nCol, u32 nRow, char name [ ] );
void DRowDump ( double *m, u32 ncol, char name [ ] );
void DColDump ( double *m, u32 nRow, char name [ ] );
void SMatDump ( s32 *m, u32 nCol, u32 nRow, char name [ ] );
void SRowDump ( s32 *m, u32 ncol, char name [ ] );
void SColDump ( s32 *m, u32 nRow, char name [ ] );
void UMatDump ( u32 *m, u32 nCol, u32 nRow, char name [ ] );
void URowDump ( u32 *m, u32 ncol, char name [ ] );
void UColDump ( u32 *m, u32 nRow, char name [ ] );

double SigmaYFinal = 0.0;
double SigmaXYFinal = 0.0;
cplx RxPrsTDCorrDfpp [ SYS_MAX_CARRIERS ];
cplx RxPrsFDCorrDfpp [ SYS_MAX_CARRIERS ];
double RxPrsFDFinalCorrPhaseErr [ SYS_MAX_CARRIERS ];

const double FrqCorrFactor [ SYS_N_CHAN_PRSBPA ] = { -0.20, -0.15, -0.10, -0.05, 0.00, 0.05, 0.10, 0.15, 0.20 };

static s32 BPARowIndex = 0;
static s32 ridx = 0;

static double SigmaErrPhiCIR = 0.0;
static double Theta0CIR = 0.0;
//static double ErrPhiCIR          [ SYS_MAX_CARRIERS  ];
//static cplx   RxPrsFDCIR         [ SYS_MAX_CARRIERS  ];

static cplx PRSBuf [ SYS_MAX_CARRIERS ];

static cplx TxPrs [ SYS_MAX_CARRIERS ];
static double TxPrsPhase [ SYS_MAX_CARRIERS ];
static cplx RxPrsTD [ SYS_MAX_CARRIERS ];
static cplx RxPrsTDCorr [ SYS_N_CHAN_PRSBPA + 1 ] [ SYS_MAX_CARRIERS ];
static cplx RxPrsFD [ SYS_N_CHAN_PRSBPA + 1 ] [ SYS_MAX_CARRIERS ];
static cplx RxPrsFDCorr [ SYS_N_CHAN_PRSBPA + 1 ] [ SYS_MAX_CARRIERS ];
static cplx RxPrsFDFinalCorr [ SYS_MAX_CARRIERS ];
static cplx RxPrsTDFinalCorr [ SYS_MAX_CARRIERS ];
static double RxPrsPhase [ SYS_N_CHAN_PRSBPA + 1 ] [ SYS_MAX_CARRIERS ];
static double RxPrsErrPhase [ SYS_N_CHAN_PRSBPA + 1 ] [ SYS_MAX_CARRIERS ];
static s32 sampleIdx1 [ SYS_MAX_CARRIERS ];
static s32 sampleIdx2 [ SYS_MAX_CARRIERS ];
static u32 nSams = 0;

static cplx complexGain [ SYS_MAX_CARRIERS ];

static cplx RcvDataTd [ SYS_MAX_CARRIERS ];
static cplx RcvDataTdp1 [ SYS_MAX_CARRIERS ];

static cplx RcvDataFd [ SYS_MAX_CARRIERS ];
static cplx RcvDataFdp1 [ SYS_MAX_CARRIERS ];
static cplx TxDataFd [ SYS_MAX_CARRIERS ];
static cplx TxDataFdEq [ SYS_MAX_CARRIERS ];
static cplx RxDataTDFinalCorr [ SYS_MAX_CARRIERS ];
static cplx RxDataFDFinalCorra [ SYS_MAX_CARRIERS ];
static cplx RxDataFDFinalCorrb [ SYS_MAX_CARRIERS ];
static cplx RxDataFDFinalCorrc [ SYS_MAX_CARRIERS ];
static cplx RxDataFDFinalCorrd [ SYS_MAX_CARRIERS ];
static cplx RxDataTDFinalCorrp1 [ SYS_MAX_CARRIERS ];
static cplx RxDataFDFinalCorrap1 [ SYS_MAX_CARRIERS ];
static cplx RxDataFDFinalCorrbp1 [ SYS_MAX_CARRIERS ];
static cplx RxDataFDFinalCorrcp1 [ SYS_MAX_CARRIERS ];
static cplx RxDataFDFinalCorrdp1 [ SYS_MAX_CARRIERS ];
static DemodData DmodRxData4QAM [ SYS_MAX_CARRIERS ];
static DemodData DmodRxData4QAMp1 [ SYS_MAX_CARRIERS ];
static DemodData DmodTxData4QAM [ SYS_MAX_CARRIERS ];
static u32 DModData [ SYS_MAX_CARRIERS ];
static u32 DModRxData [ SYS_MAX_CARRIERS ];

static double YA [ SYS_N_CHAN_PRSBPA + 1 ] [ SYS_MAX_CARRIERS ];
static double MeanYA [ SYS_N_CHAN_PRSBPA + 1 ] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
static double NormYA [ SYS_N_CHAN_PRSBPA + 1 ] [ SYS_MAX_CARRIERS ];
static double SigmaYA [ SYS_N_CHAN_PRSBPA + 1 ] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
static double SigmaY2A [ SYS_N_CHAN_PRSBPA + 1 ] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
static double SigmaXYA [ SYS_N_CHAN_PRSBPA + 1 ] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

static double YB [ SYS_N_CHAN_PRSBPA + 1 ] [ SYS_MAX_CARRIERS ];
static double MeanYB [ SYS_N_CHAN_PRSBPA + 1 ] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
static double NormYB [ SYS_N_CHAN_PRSBPA + 1 ] [ SYS_MAX_CARRIERS ];
static double SigmaYB [ SYS_N_CHAN_PRSBPA + 1 ] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
static double SigmaY2B [ SYS_N_CHAN_PRSBPA + 1 ] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
static double SigmaXYB [ SYS_N_CHAN_PRSBPA + 1 ] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

static double Y [ SYS_N_CHAN_PRSBPA + 1 ] [ SYS_MAX_CARRIERS ];
static double MeanY [ SYS_N_CHAN_PRSBPA + 1 ] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
static double NormY [ SYS_N_CHAN_PRSBPA + 1 ] [ SYS_MAX_CARRIERS ];
static double SigmaY [ SYS_N_CHAN_PRSBPA + 1 ] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
static double SigmaY2 [ SYS_N_CHAN_PRSBPA + 1 ] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
static double SigmaXY [ SYS_N_CHAN_PRSBPA + 1 ] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

static double X [ SYS_MAX_CARRIERS ];
static double MeanX = 0.0;
static double NormX [ SYS_MAX_CARRIERS ];
static double SigmaX = 0.0;
static double SigmaX2 = 0.0;

static cplx PRSBufRot [ SYS_MAX_CARRIERS ];
static cplx RxPrsBufRot [ SYS_MAX_CARRIERS ];

static double LUTPE [ 256 ] [ 256 ]; // ToDo: Probably Change
double BerQPSK = 0.0;
double Ber4QAM = 0.0;
double Ber4QAMp1 = 0.0;
double Ber16QAM = 0.0;
double Ber64QAM = 0.0;
double Qfr4QAM = 0.0;
double Qfr16QAM = 0.0;
double Qfr64QAM = 0.0;


double alpha,beta,alphap1,betap1;


// =================================
// Boilerplate Function Declarations
// =================================
void PrsBPAInit (
        PrsBPA *pd,
        const char *Name,
        u32 BaseAddr,
        SysCtrl *pSysCtrl,
        SysDvce *pSysDvce,
        SysDesc *pSysDesc,
        SysPtrs *pSysPtrs,
        SysData *pSysData,
        PrsDesc *pPrsBPADesc,
        PrsBPAData *pPrsBPAData
        );

static void SetParams ( PrsBPA *pd );
static void SysDvceDump ( PrsBPA *pd );
static void SysDescDump ( PrsBPA *pd );
static void SysDataDump ( PrsBPA *pd );
static void DvcDataInit ( PrsBPA *pd );
static void DvcDataUpd ( PrsBPA *pd );
static void DvcDescDump ( PrsBPA *pd );
static void DvcDataDump ( PrsBPA *pd );

// static non-member functions
static void ReadPRSData ( PrsBPA *pd );
static void CopyPRSPhase ( PrsBPA *pd );
static void CalcDelta ( PrsBPA *pd, s32 ridx );
static void CalcPhaseError ( PrsBPA *pd );
static void FinalCorr ( PrsBPA *pd );
static void CalcLinReg ( PrsBPA *pd );
static void CalcCorr ( PrsBPA *pd, u32 cidx, double freqCorrFact, double timeCorrFact );
static void ReadRxData ( PrsBPA *pd , int nSymbol);
static void CopyTxData ( PrsBPA *pd , int nSymbol);
static void DemodRxData ( PrsBPA *pd );
static void CreatePhaseErrorLUT ();
static void FreqDataCorr ( PrsBPA *pd , int nSymbol);
static void TimeDataCorr ( PrsBPA *pd , int nSymbol);
static void NormalizeDataVal ( PrsBPA *pd );
static void FineDataCorr ( PrsBPA *pd );
//static void CalcBER( PrsBPA *pd, s32 ridx );

// =================================================
// Device Structure Member Functions (public access)
// =================================================
// Insert device-specific member function prototypes here
// Example: static inline void Start         ( PrsBPA *pd ) { hwWriteCtrlReg( SYSDVC( BaseAddr ), PRSBPA_CTRL_REG << 2, PRSBPA_START ); }
// Example: static void        WriteBuf      ( PrsBPA *pd, u32 addr, u32 data );
static void Dump ( PrsBPA *pd );
static void DumpCorr ( PrsBPA *pd );
static void DumpCorrSymbFd ( PrsBPA *pd );
static void DumpCorrSymbTd ( PrsBPA *pd );
static void CalcDeltas ( PrsBPA *pd );
static void DumpCorrDataSymbFd ( PrsBPA *pd );
static void DumpCorrDataSymbTd ( PrsBPA *pd );
static void RxDataSymbolDumpRawFd ( PrsBPA *pd );
static void RxDataSymbolDumpFCFd ( PrsBPA *pd );
static void RxDataSymbolDumpFTCFd ( PrsBPA *pd );

// =================================================
// Device Shell Commands
// =================================================

static inline void ShCalcDeltas ( shell *psh ) {
    PrsBPA *pd = &RxPrsBPA;
    CalcDeltas ( pd );
}
static inline void DumpRxPrsTD ( shell *psh ) {
    PrsBPA *pd = &RxPrsBPA;
    CMatDump ( &RxPrsTD [ 0 ], 0, SYSPAR( nFFT ), "Rx PRS time domain" );
}
static inline void DumpRxPrsTDCorr ( shell *psh ) {
    PrsBPA *pd = &RxPrsBPA;
    CMatDump ( &RxPrsTDCorr [ 0 ] [ 0 ], SYS_N_CHAN_PRSBPA, SYSPAR( nFFT ), "Rx PRS time domain corrected" );
}
static inline void DumpRxPrsFD ( shell *psh ) {
    PrsBPA *pd = &RxPrsBPA;
    CMatDump ( &RxPrsFD [ 0 ] [ 0 ], SYS_N_CHAN_PRSBPA, SYSPAR( nFFT ), "Rx PRS frequency domain" );
}
static inline void DumpRxPrsFDCorr ( shell *psh ) {
    PrsBPA *pd = &RxPrsBPA;
    CMatDump ( &RxPrsFDCorr [ 0 ] [ 0 ], SYS_N_CHAN_PRSBPA, SYSPAR( nFFT ), "Rx PRS frequency domain corrected" );
}
static inline void DumpRxPrsPhase ( shell *psh ) {
    PrsBPA *pd = &RxPrsBPA;
    DMatDump ( &RxPrsPhase [ 0 ] [ 0 ], SYS_N_CHAN_PRSBPA, SYSPAR( nFFT ), "Rx PRS phase" );
}
static inline void DumpRxPrsErrPhase ( shell *psh ) {
    PrsBPA *pd = &RxPrsBPA;
    DMatDump ( &RxPrsErrPhase [ 0 ] [ 0 ], SYS_N_CHAN_PRSBPA, SYSPAR( nFFT ), "Rx PRS phase error" );
}
static inline void DumpRxPrsFDFinalCorr ( shell *psh ) {
    PrsBPA *pd = &RxPrsBPA;
    CMatDump ( &RxPrsFDFinalCorr [ 0 ], 0, SYSPAR( nFFT ), "Rx PRS freq domain final Corr" );
}
static inline void DumpRxPrsTDFinalCorr ( shell *psh ) {
    PrsBPA *pd = &RxPrsBPA;
    CMatDump ( &RxPrsTDFinalCorr [ 0 ], 0, SYSPAR( nFFT ), "Rx PRS time domain final Corr" );
}
static inline void DumpTxPrs ( shell *psh ) {
    PrsBPA *pd = &RxPrsBPA;
    CMatDump ( &TxPrs [ 0 ], 0, SYSPAR( nFFT ), "Tx PRS" );
}
static inline void DumpTxPrsPhase ( shell *psh ) {
    PrsBPA *pd = &RxPrsBPA;
    DMatDump ( &TxPrsPhase [ 0 ], 0, SYSPAR( nFFT ), "Tx PRS phase" );
}
static inline void DumpTxData ( shell *psh ) {
    PrsBPA *pd = &RxPrsBPA;
    CMatDump ( &TxDataFd [ 0 ], 0, SYSPAR( nFFT ), "Tx Mod Data frequency domain" );
}
static inline void DumpRxTDData ( shell *psh ) {
    PrsBPA *pd = &RxPrsBPA;
    CMatDump ( &RcvDataTd [ 0 ], 0, SYSPAR( nFFT ), "Rx time domain" );
}
static inline void DumpRxFDData ( shell *psh ) {
    PrsBPA *pd = &RxPrsBPA;
    CMatDump ( &RcvDataFd [ 0 ], 0, SYSPAR( nFFT ), "Rx freq domain" );
}
static inline void DumpRxTDDataFinalCorr ( shell *psh ) {
    PrsBPA *pd = &RxPrsBPA;
    CMatDump ( &RxDataTDFinalCorr [ 0 ], 0, SYSPAR( nFFT ), "Rx time domain final corr" );
}
//static inline void DumpRxFDDataFinalCorr   ( shell *psh ) { PrsBPA *pd = &RxPrsBPA; CMatDump( &RxDataFDFinalCorr     [ 0 ], 0                , SYSPAR( nFFT ), "Rx freq domain final corr"          ); }
static inline void DumpX ( shell *psh ) {
    PrsBPA *pd = &RxPrsBPA;
    DMatDump ( &X [ 0 ], 0, SYSPAR( nFFT ), "Linear interpolation X" );
}
static inline void DumpNormX ( shell *psh ) {
    PrsBPA *pd = &RxPrsBPA;
    DMatDump ( &NormX [ 0 ], 0, SYSPAR( nFFT ), "Linear interpolation X - mean( X )" );
}
static inline void DumpSampleIdx1 ( shell *psh ) {
    PrsBPA *pd = &RxPrsBPA;
    SMatDump ( &sampleIdx1 [ 0 ], 0, SYSPAR( nFFT ), "Sample index 1" );
}
static inline void DumpSampleIdx2 ( shell *psh ) {
    PrsBPA *pd = &RxPrsBPA;
    SMatDump ( &sampleIdx2 [ 0 ], 0, SYSPAR( nFFT ), "Sample index 2" );
}

static inline void DumpY ( shell *psh ) {
    PrsBPA *pd = &RxPrsBPA;
    DMatDump ( &Y [ 0 ] [ 0 ], SYS_N_CHAN_PRSBPA, SYSPAR( nFFT ), "Y" );
}
static inline void DumpNormY ( shell *psh ) {
    PrsBPA *pd = &RxPrsBPA;
    DMatDump ( &NormY [ 0 ] [ 0 ], SYS_N_CHAN_PRSBPA, SYSPAR( nFFT ), "NormY" );
}
static inline void DumpMeanY ( shell *psh ) {
    DMatDump ( &MeanY [ 0 ], SYS_N_CHAN_PRSBPA, 0, "Linear interpolation Y" );
}
static inline void DumpSigmaY ( shell *psh ) {
    DMatDump ( &SigmaY [ 0 ], SYS_N_CHAN_PRSBPA, 0, "Linear interpolation sigma Y" );
}
static inline void DumpSigmaY2 ( shell *psh ) {
    DMatDump ( &SigmaY2 [ 0 ], SYS_N_CHAN_PRSBPA, 0, "Linear interpolation sigma Y^2" );
}
static inline void DumpSigmaXY ( shell *psh ) {
    DMatDump ( &SigmaXY [ 0 ], SYS_N_CHAN_PRSBPA, 0, "Linear interpolation Sigma XY" );
}
static void DumpChanRes ( shell *psh );

shellCommandDef PrsBPACmdDef [ ] = {
        { "all", "", 0, "", NULL },
        { "all", "q", 0, "                                        Exit this shell", shellConsoleExit },
        { "all", "", 0, "", NULL },
        { "all", "Calc", 0, "                                        Calculate deltas", ShCalcDeltas },
        { "all", "RxPTD", 0, "                                        Dump Rx PRS Time Domain", DumpRxPrsTD },
        { "all", "RxPTDC", 0, "                                        Dump Rx PRS Time Domain Corrected",
                DumpRxPrsTDCorr },
        { "all", "RxPFD", 0, "                                        Dump Rx PRS Freq Domain", DumpRxPrsFD },
        { "all", "RxPFDC", 0, "                                        Dump Rx PRS Freq Domain Corrected",
                DumpRxPrsFDCorr },
        { "all", "RxPPh", 0, "                                        Dump Rx PRS Phase", DumpRxPrsPhase },

        { "all", "", 0, "", NULL },
        { "all", "ErrPh", 0, "                                        Dump Phase Error", DumpRxPrsErrPhase },
        { "all", "", 0, "", NULL },
        { "all", "RxFDFinC", 0, "                                        Dump Rx PRS Final Corr Freq Domain",
                DumpRxPrsFDFinalCorr },
        { "all", "RxTDFinC", 0, "                                        Dump Rx PRS Final Corr Time Domain",
                DumpRxPrsTDFinalCorr },
        { "all", "TxP", 0, "                                        Dump Tx PRS", DumpTxPrs },
        { "all", "TxPPh", 0, "                                        Dump Tx PRS Phase", DumpTxPrsPhase },
        { "all", "", 0, "", NULL },
        { "all", "TxD", 0, "                                        Dump Tx Data (Modulated Freq Domain)", DumpTxData },
        { "all", "RxTDD", 0, "                                        Dump Rx Data Time Domain            ",
                DumpRxTDData },
        { "all", "RxFDD", 0, "                                        Dump Rx Data Freq Domain            ",
                DumpRxFDData },
        { "all", "RxTDC", 0, "                                        Dump Rx Data Time Domain Corrected  ",
                DumpRxTDDataFinalCorr },
        //{ "all"  , "RxFDC"      ,  0, "                                        Dump Rx Data Freq Domain Corrected  "   , DumpRxFDDataFinalCorr },
        { "all", "", 0, "", NULL },
        { "all", "X", 0, "                                        Dump X", DumpX },
        { "all", "NormX", 0, "                                        Dump X - mean( x )", DumpNormX },
        { "all", "", 0, "", NULL },
        { "all", "Y", 0, "                                        Dump Y", DumpY },
        { "all", "MeanY", 0, "                                        Dump Mean( Y )", DumpMeanY },
        { "all", "NormY", 0, "                                        Dump Y  - mean( Y )", DumpNormY },
        { "all", "SigY", 0, "                                        Dump Sigma Y", DumpSigmaY },
        { "all", "SigXY", 0, "                                        Dump Sigma XY", DumpSigmaXY },
        { "all", "SigY2", 0, "                                        Dump Sigma Y^2b", DumpSigmaY2 },
        { "all", "", 0, "", NULL },
        { "all", "SI1", 0, "                                        Dump Sample Index 1", DumpSampleIdx1 },
        { "all", "SI2", 0, "                                        Dump Sample Index 2", DumpSampleIdx2 },
        { "all", "CRes", 0, "                                        Dump Chan Results", DumpChanRes },
        { "\0" }
};

// Declare init structure here with member functions
PrsBPA PrsBPAInitStr = {
NULL,              // SysDvce        *pSysDvce;
        NULL,              // PrsDesc        *pDvcDesc;
        NULL,              // PrsBPAData     *pDvcData;
        Dump,              // void (* Dump           ) ( PrsBPA *pd );
        DumpCorr,          // void (* DumpCorr       ) ( PrsBPA *pd );
        DumpCorrSymbFd,    // void (* DumpCorrSymbFd ) ( PrsBPA *pd, u32 offs );
        DumpCorrSymbTd,    // void (* DumpCorrSymbTd ) ( PrsBPA *pd, u32 offs );
        DumpCorrDataSymbFd,    // void (* DumpCorrDataSymbFd ) ( PrsBPA *pd );
        DumpCorrDataSymbTd,    // void (* DumpCorrDataSymbTd ) ( PrsBPA *pd );
        CalcDeltas,         // void (* CalcDeltas     ) ( PrsBPA *pd );
        RxDataSymbolDumpRawFd,
        RxDataSymbolDumpFCFd,
        RxDataSymbolDumpFTCFd
};

//===================
// Initialize Device
//===================
void PrsBPAInit ( PrsBPA *pd, const char *Name, u32 BaseAddr, SysCtrl *pSysCtrl, SysDvce *pSysDvce, SysDesc *pSysDesc,
        SysPtrs *pSysPtrs, SysData *pSysData, PrsDesc *pDvcDesc, PrsBPAData *pDvcData )
{
    *pd = PrsBPAInitStr;
    pd->pSysDvce = pSysDvce;
    SysDvceInit (
            pd->pSysDvce,
            Name,
            pSysCtrl,
            pSysDesc,
            pSysPtrs,
            pSysData,
            NULL,
            BaseAddr,
            BPA_HW_INFO_REG,
            BPA_CTRL_REG,
            BPA_INFO_REG,
            BPA_STAT_REG,
            BPA_RESET,
            BPA_INT_ACK,
            0,
            0,
            0,
            0,
            0,
            BPA_TX_DATA_CTR_REG,
            BPA_TX_SYMB_CTR_REG,
            BPA_TX_FRAME_CTR_REG,
            BPA_TX_BLOCK_CTR_REG,
            0,
            0
            );
    pd->pDvcDesc = pDvcDesc;
    pd->pDvcData = pDvcData;
    DvcDataInit ( pd );
    SYSRESET( pd );
    SetParams ( pd );
    txPrsPhi.init = false;
    txPrsPhi.Initialise = PrsPhiInit;
    txPrsPhi.Phase = PrsPhiPhase;
}

//========================
// Dump Device to Console
//========================
static void Dump ( PrsBPA *pd ) {
    cprintf( "\n" );
    SysDvceDump ( pd );
    SysDescDump ( pd );
    DvcDescDump ( pd );
    SysDataDump ( pd );
    DvcDataDump ( pd );
    cprintf( "\n" );
}

// End of boilerplate code
//========================
// Device-specific code
//========================
static void SetParams ( PrsBPA *pd )
{
// Set PL device slave regs here
// Example SYSWRSLVREG( pd, PRSBPA_BOB_REG,    SYSPAR( Bob   ));
// Example SYSWRSLVREG( pd, PRSBPA_ALICE_REG , SYSPAR( Alice ));
}

static void DvcDataInit ( PrsBPA *pd )
{
// Initialise DvcData members here
// Example: DVCDAT( Bob ) = 0;
}

static void DvcDataUpd ( PrsBPA *pd )
{
// Update DvcData members here
// Example: DVCDAT( Bob ) = SYSRDSLVREG( pd, PRSBPA_BOB_REG );
}

static void DvcDataDump ( PrsBPA *pd )
{
    DvcDataUpd ( pd );
// Insert device-specific SysDesc members here
// Example: cprintf("\nBob  : %u", toUint( DVCDAT( Bob )) );
}

static void DvcDescDump ( PrsBPA *pd )
{
// Insert device-specific PrsBPADesc members here
// Example: cprintf("\nAlice: %u", toUint( SYSPAR( Alice )) );
}

static void SysDvceDump ( PrsBPA *pd )
{
    SYSDVCEDUMP( pd );
}

static void SysDataDump ( PrsBPA *pd )
{
    SYSDATADUMP( pd );
}

static void SysDescDump ( PrsBPA *pd )
{
// Insert device-specific PrsBPADesc members here
// Example: cprintf("\nAlice: %u", toUint( SYSPAR( Alice )) );
}

// Insert device-specific code here
//======================================
// Device-specific non member functions
//======================================
// Example: static inline void EnterAdmin    ( PrsBPA *pd ) { hwWriteCtrlReg( SYSDVC( BaseAddr ), PRSBPA_CTRL_REG << 2, PRSBPA_ENTER_ADMIN ); }

//=================================
// Device-specific member functions
//=================================

// Example:
// static u32 ReadBuf( PrsBPA *pd, u32 addr )
// {
//     return ( SYSPTR( pBuf[ addr ] ));
// }

static void ClearData ( PrsBPA * pd )
{
    BPARowIndex = 0;
    ridx = 0;
    nSams = 0;

    SigmaErrPhiCIR = 0.0;
    Theta0CIR = 0.0;

//    ( void ) memsetd( ( void * ) ErrPhiCIR        , 0.0, SIZE( ErrPhiCIR        , double ) );
//    ( void ) memsetc( ( void * ) RxPrsFDCIR       , 0.0, SIZE( RxPrsFDCIR       , cplx   ) );
    ( void ) memsetc ( ( void * ) PRSBuf, 0.0, SIZE( PRSBuf, cplx ) );
    ( void ) memsetc ( ( void * ) PRSBufRot, 0.0, SIZE( PRSBufRot, cplx ) );
    ( void ) memsetc ( ( void * ) RxPrsBufRot, 0.0, SIZE( RxPrsBufRot, cplx ) );
    ( void ) memsetc ( ( void * ) TxPrs, 0.0, SIZE( TxPrs, cplx ) );
    ( void ) memsetd ( ( void * ) TxPrsPhase, 0.0, SIZE( TxPrsPhase, double ) );
    ( void ) memsetc ( ( void * ) RxPrsTD, 0.0, SIZE( RxPrsTD, cplx ) );
    ( void ) memsetc ( ( void * ) RxPrsTDCorr, 0.0, SIZE( RxPrsTDCorr, cplx ) );
    ( void ) memsetc ( ( void * ) RxPrsFD, 0.0, SIZE( RxPrsFD, cplx ) );
    ( void ) memsetc ( ( void * ) RxPrsFDCorr, 0.0, SIZE( RxPrsFDCorr, cplx ) );
    ( void ) memsetc ( ( void * ) RxPrsFDFinalCorr, 0.0, SIZE( RxPrsFDFinalCorr, cplx ) );
    ( void ) memsetc ( ( void * ) complexGain, 1.0, SIZE( complexGain, cplx ) );
    ( void ) memsetc ( ( void * ) RxPrsTDFinalCorr, 0.0, SIZE( RxPrsTDFinalCorr, cplx ) );
    ( void ) memsetd ( ( void * ) RxPrsPhase, 0.0, SIZE( RxPrsPhase, double ) );
    ( void ) memsetd ( ( void * ) RxPrsErrPhase, 0.0, SIZE( RxPrsErrPhase, double ) );
    ( void ) memsetd ( ( void * ) NormX, 0.0, SIZE( NormX, double ) );

    ( void ) memsetc ( ( void * ) RcvDataTd, 0.0, SIZE( RcvDataTd, cplx ) );
    ( void ) memsetc ( ( void * ) RcvDataFd, 0.0, SIZE( RcvDataFd, cplx ) );
    ( void ) memsetc ( ( void * ) TxDataFd, 0.0, SIZE( TxDataFd, cplx ) );
    ( void ) memsetc ( ( void * ) TxDataFdEq, 0.0, SIZE( TxDataFdEq, cplx ) );
    ( void ) memsetc ( ( void * ) RxDataTDFinalCorr, 0.0, SIZE( RxDataTDFinalCorr, cplx ) );
    ( void ) memsetc ( ( void * ) RxDataFDFinalCorra, 0.0, SIZE( RxDataFDFinalCorra, cplx ) );
    ( void ) memsetc ( ( void * ) RxDataFDFinalCorrb, 0.0, SIZE( RxDataFDFinalCorrb, cplx ) );
    ( void ) memsetc ( ( void * ) RxDataFDFinalCorrc, 0.0, SIZE( RxDataFDFinalCorrc, cplx ) );
    ( void ) memsetc ( ( void * ) RxDataFDFinalCorrd, 0.0, SIZE( RxDataFDFinalCorrd, cplx ) );
    ( void ) memsetd ( ( void * ) DModData, 0.0, SIZE( DModData, double ) );
    ( void ) memsetd ( ( void * ) DModRxData, 0.0, SIZE( DModRxData, double ) );

    ( void ) memset ( ( void * ) sampleIdx1, 0, sizeof ( sampleIdx1 ) );
    ( void ) memset ( ( void * ) sampleIdx2, 0, sizeof ( sampleIdx2 ) );

    ( void ) memsetd ( ( void * ) &X [ 0 ], 0.0, SIZE( X, double ) );
    MeanX = 0.0;
    SigmaX = 0.0;
    SigmaX2 = 0.0;
    ( void ) memsetd ( ( void * ) &NormX [ 0 ], 0.0, SIZE( NormX, double ) );
    ( void ) memsetd ( ( void * ) &YA [ 0 ], 0.0, SIZE( YA, double ) );
    ( void ) memsetd ( ( void * ) &MeanYA [ 0 ], 0.0, SIZE( MeanYA, double ) );
    ( void ) memsetd ( ( void * ) &NormYA [ 0 ], 0.0, SIZE( NormYA, double ) );
    ( void ) memsetd ( ( void * ) &SigmaYA [ 0 ], 0.0, SIZE( SigmaYA, double ) );
    ( void ) memsetd ( ( void * ) &SigmaY2A [ 0 ], 0.0, SIZE( SigmaY2A, double ) );
    ( void ) memsetd ( ( void * ) &SigmaXYA [ 0 ], 0.0, SIZE( SigmaXYA, double ) );
    ( void ) memsetd ( ( void * ) &YB [ 0 ], 0.0, SIZE( YB, double ) );
    ( void ) memsetd ( ( void * ) &MeanYB [ 0 ], 0.0, SIZE( MeanYB, double ) );
    ( void ) memsetd ( ( void * ) &NormYB [ 0 ], 0.0, SIZE( NormYB, double ) );
    ( void ) memsetd ( ( void * ) &SigmaYB [ 0 ], 0.0, SIZE( SigmaYB, double ) );
    ( void ) memsetd ( ( void * ) &SigmaY2B [ 0 ], 0.0, SIZE( SigmaY2B, double ) );
    ( void ) memsetd ( ( void * ) &SigmaXYB [ 0 ], 0.0, SIZE( SigmaXYB, double ) );
    ( void ) memsetd ( ( void * ) &Y [ 0 ], 0.0, SIZE( Y, double ) );
    ( void ) memsetd ( ( void * ) &MeanY [ 0 ], 0.0, SIZE( MeanY, double ) );
    ( void ) memsetd ( ( void * ) &NormY [ 0 ], 0.0, SIZE( NormY, double ) );
    ( void ) memsetd ( ( void * ) &SigmaY [ 0 ], 0.0, SIZE( SigmaY, double ) );
    ( void ) memsetd ( ( void * ) &SigmaY2 [ 0 ], 0.0, SIZE( SigmaY2, double ) );
    ( void ) memsetd ( ( void * ) &SigmaXY [ 0 ], 0.0, SIZE( SigmaXY, double ) );

    BerQPSK = 0.0;
    Ber4QAM = 0.0;
    Ber16QAM = 0.0;
    Ber64QAM = 0.0;
    Qfr4QAM = 0.0;
    Qfr16QAM = 0.0;
    Qfr64QAM = 0.0;
}

static void FreqDataCorr ( PrsBPA *pd , int nSymbol)
{
    /*
     * Frequency Correction
     *  */
    double lastTheta = DVCPAR( lastTheta );
    u32 nFFT = SYSPAR( nFFT );
    double nFFTd = ( double ) nFFT;
    double G = (double) SYSPAR( nCPre );
    double fineDeltaf = DVCPAR( fineDeltaf )-DVCPAR( dataAlpha )/(2.0 * M_PI )*1.0/4.0;
    double deltafpp = DVCPAR( deltafpp );
    double deltafp = DVCPAR( deltafp );
    double deltaf = DVCPAR( DeltaF );
    double deltaftot = ( deltaf + deltafp + deltafpp + fineDeltaf);


    double thetaN = lastTheta ;//+ ( deltaftot * 2.0 * M_PI ) / nFFTd*((nFFTd+G)*(nSymbol+1));


    for ( int n = 0; n < nFFT; n++ ) {
        thetaN += ( deltaftot * 2.0 * M_PI ) / nFFTd;
        RxDataTDFinalCorr [ n ] = RcvDataTd [ n ] * cexp ( I * thetaN );
        RxDataTDFinalCorrp1 [ n ] = RcvDataTdp1 [ n ] * cexp ( I * thetaN );
    }
    fft ( RxDataFDFinalCorra, RxDataTDFinalCorr, nFFT );
    fft ( RxDataFDFinalCorrap1, RxDataTDFinalCorrp1, nFFT );

    DVCPAR( lastTheta ) = lastTheta + ( deltaftot * 2.0 * M_PI ) / nFFTd*(nFFTd+G);
    DVCPAR( fineDeltaf ) = fineDeltaf;
}

static void TimeDataCorr ( PrsBPA *pd , int nSymbol)
{
    double deltatime = DVCPAR( DeltaT ) - (double)((int)(1+(double)(0*nSymbol)));
    double deltatp = DVCPAR( deltatp );
    double deltatpp = DVCPAR( deltatpp );
    u32 nFFT = SYSPAR( nFFT );
    double fineDeltat = DVCPAR( fineDeltat )-DVCPAR( dataBeta )/(2.0 * M_PI )*1.0*nFFT;
    double deltattot =deltatime+ deltatp + deltatpp +fineDeltat;

    /*
     * Time Correction
     *  */
    double deltaTN = 0;//( deltattot * 2.0 * M_PI ) / nFFT;
    for ( int n = 0; n < nFFT / 2; n++ ) {
        if ( ( sampleIdx1 [ n ] >= 0 ) ) {
            deltaTN += ( deltattot * 2.0 * M_PI ) / nFFT;
            RxDataFDFinalCorrb [ n ] = RxDataFDFinalCorra [ n ] * cexp ( I * deltaTN );
            RxDataFDFinalCorrbp1 [ n ] = RxDataFDFinalCorrap1 [ n ] * cexp ( I * deltaTN );
        }
    }
    deltaTN = 0;//- ( deltattot * 2.0 * M_PI ) / nFFT;
    for ( int n = nFFT - 1; n > nFFT / 2; n-- ) {
        if ( ( sampleIdx1 [ n ] >= 0 ) ) {
            deltaTN -= ( deltattot * 2.0 * M_PI ) / nFFT;
            RxDataFDFinalCorrb [ n ] = RxDataFDFinalCorra [ n ] * cexp ( I * deltaTN );
            RxDataFDFinalCorrbp1 [ n ] = RxDataFDFinalCorrap1 [ n ] * cexp ( I * deltaTN );
        }
    }
    DVCPAR( fineDeltat ) = fineDeltat;
}

static void NormalizeDataVal ( PrsBPA *pd )
{
    /* Normalize Corrected Values */
    cplx DemodScaleFactor = 1448.0 + 0.0 * I;
    ModType Modulation = SYSPAR( ModType );
    u32 nFFT = SYSPAR( nFFT );
    double absMeanrx = 0.0;
    double absMeanrxp1 = 0.0;
    double absMeantx = 0.0;
    double nMean = 0.0;

    /*Create Demodulator and Normalize Corrected Values */
    if ( Modulation == MOD_4QAM ) {
        DemodScaleFactor = 1448.0 + 0.0 * I;
    } else if ( Modulation == MOD_16QAM ) {
        DemodScaleFactor = 1236.0 + 0.0 * I;
    } else if ( Modulation == MOD_64QAM ) {
        DemodScaleFactor = 1558.0 + 0.0 * I;
    }
    for ( int i = 0; i < nFFT; i++ ) {
        if ( ( sampleIdx1 [ i ] >= 0 ) )
        {
            absMeanrx += cabs ( RxDataFDFinalCorrb [ i ] );
            absMeanrxp1 += cabs ( RxDataFDFinalCorrbp1 [ i ] );
            absMeantx += cabs ( TxDataFd [ i ] );
            nMean += 1;
        }
    }
    absMeanrx /= nMean;
    absMeanrxp1 /= nMean;
    absMeantx /= nMean;

    for ( u32 i = 1; i < nFFT; i++ ) {
        RxDataFDFinalCorrc [ i ] = RxDataFDFinalCorrb [ i ] * DemodScaleFactor / ( absMeanrx + 0.0 * I )
                * complexGain [ i ];
        RxDataFDFinalCorrcp1 [ i ] = RxDataFDFinalCorrbp1 [ i ] * DemodScaleFactor / ( absMeanrxp1 + 0.0 * I )
                * complexGain [ i ];
        TxDataFdEq [ i ] = TxDataFd [ i ] * DemodScaleFactor / ( absMeantx + 0.0 * I );
    }
}
static double YD [ 768 ];
static double XD [ 768 ];
static double NormYD [ 768 ];

static void FineDataCorr ( PrsBPA *pd )
{
    /*
     * Fine Correction
     * */
    s32 nFFT = SYSPAR( nFFT );
    u8 modType = SYSPAR( ModType );

    double th = 0;

    if ( modType == 1 )
        th = 0;
    else if ( modType == 2 )
        th = 2048 / 16 * 2 / sqrt ( 10 ) / 2;
//        th = 2048 / 2/32;
    else if ( modType == 3 )
        th = 2048 / 16 * 6 / sqrt ( 42 ) / 2;
//        th = 3*2048 / 4/32;

    /* Create Phase Error LUT */
    CreatePhaseErrorLUT ();

    /* Calculate alpha and Beta for Data Values */
    s8 RxDataFDFinalCorrI = 0;
    s8 RxDataFDFinalCorrQ = 0;
    Cplx32 RxDataFDFinalCorrcx;
    double ErrDataPhase = 0.0;
    double MeanYD = 0.0;
    double SigmaX2D = 0.0;
    double SigmaXYD = 0.0;
    u32 k = 0;

    for ( int i = 0; i < nFFT; i++ ) {
        if ( sampleIdx1 [ i ] >= 0 ) {
            RxDataFDFinalCorrcx = CplxtoCplx32 ( RxDataFDFinalCorrc [ i ] );
            RxDataFDFinalCorrI = ( RxDataFDFinalCorrcx.real ) / 16;
            RxDataFDFinalCorrQ = ( RxDataFDFinalCorrcx.imag ) / 16;
            if ( ( abs ( RxDataFDFinalCorrI ) > th ) && ( abs ( RxDataFDFinalCorrQ ) > th ) )
                    {
                ErrDataPhase = LUTPE [ ( s16 ) RxDataFDFinalCorrI + 128 ] [ ( s16 ) RxDataFDFinalCorrQ + 128 ];
                YD [ k ] = ErrDataPhase;
                MeanYD += ErrDataPhase;
                XD [ k ] = ( double ) sampleIdx2 [ i ];
                SigmaX2D += XD [ k ] * XD [ k ];
                k++;
            }
        }
    }
    MeanYD = MeanYD / k;
    for ( int i = 0; i < k; i++ ) {
        NormYD [ i ] = YD [ i ] - MeanYD;
        SigmaXYD += NormYD [ i ] * XD [ i ];
    }
    beta = -SigmaXYD / SigmaX2D;
    alpha = -MeanYD;

    MeanYD = 0.0;
    SigmaX2D = 0.0;
    SigmaXYD = 0.0;
    k = 0;

    for ( int i = 0; i < nFFT; i++ ) {
        if ( sampleIdx1 [ i ] >= 0 ) {
            RxDataFDFinalCorrcx = CplxtoCplx32 ( RxDataFDFinalCorrcp1 [ i ] );
            RxDataFDFinalCorrI = ( RxDataFDFinalCorrcx.real ) / 16;
            RxDataFDFinalCorrQ = ( RxDataFDFinalCorrcx.imag ) / 16;
            if ( ( abs ( RxDataFDFinalCorrI ) > th ) && ( abs ( RxDataFDFinalCorrQ ) > th ) )
                    {
                ErrDataPhase = LUTPE [ ( s16 ) RxDataFDFinalCorrI + 128 ] [ ( s16 ) RxDataFDFinalCorrQ + 128 ];
                YD [ k ] = ErrDataPhase;
                MeanYD += ErrDataPhase;
                XD [ k ] = ( double ) sampleIdx2 [ i ];
                SigmaX2D += XD [ k ] * XD [ k ];
                k++;
            }
        }
    }
    MeanYD = MeanYD / k;
    for ( int i = 0; i < k; i++ ) {
        NormYD [ i ] = YD [ i ] - MeanYD;
        SigmaXYD += NormYD [ i ] * XD [ i ];
    }
    betap1 = -SigmaXYD / SigmaX2D;
    alphap1 = -MeanYD;

    double deltaFine = 0.0;

    //correction=exp(1j*(alpha + beta*[0:N/2-1 -N/2:-1]))';
    /* Apply Fine Correction */
    for ( int n = 0; n < nFFT / 2; n++ ) {
        if ( sampleIdx1 [ n ] >= 0 ) {
            deltaFine = - ( ( n * beta ) + alpha );         // * 2.0 * M_PI; //ToDo: Changed n for deltaFine
            RxDataFDFinalCorrd [ n ] = RxDataFDFinalCorrc [ n ] * cexp ( I * deltaFine );
            deltaFine = - ( ( n * betap1 ) + alphap1 );         // * 2.0 * M_PI; //ToDo: Changed n for deltaFine
            RxDataFDFinalCorrdp1 [ n ] = RxDataFDFinalCorrcp1 [ n ] * cexp ( I * deltaFine );
        }
    }
    //deltaFine = -((deltaFine*beta)/nFFT + alpha) * 2.0 * M_PI;
    for ( int n = nFFT - 1; n > nFFT / 2; n-- ) {
        if ( sampleIdx1 [ n ] >= 0 ) {
            deltaFine = - ( ( ( n - nFFT ) * beta ) + alpha );         // * 2.0 * M_PI;
            RxDataFDFinalCorrd [ n ] = RxDataFDFinalCorrc [ n ] * cexp ( I * deltaFine );
            deltaFine = - ( ( ( n - nFFT ) * betap1 ) + alphap1 );         // * 2.0 * M_PI;
            RxDataFDFinalCorrdp1 [ n ] = RxDataFDFinalCorrcp1 [ n ] * cexp ( I * deltaFine );
        }
    }


    /*
              symbolFineFCorrection = b(1)/2/pi/(N+G);
%             lastPhase=b(1);
%             deltaPhaseCorrection = exp(1j*(lastPhase+symbolFineFCorrection*(N+G)));
%             lastPhaseCorrection= lastPhaseCorrection*deltaPhaseCorrection;

     */


}

static void DemodRxData ( PrsBPA *pd )
{
    ModType Modulation = SYSPAR( ModType );
    u32 nFFT = SYSPAR( nFFT );
    Qfr4QAM = 0.0;
    u32 errCtr = 0;
    u32 errCtrp1 = 0;
    double nMean = 0.0;
    /*Create Demodulator */
    if ( Modulation == MOD_4QAM ) {
        DDUtilCreateDemod ( MOD_4QAM );
    } else if ( Modulation == MOD_16QAM ) {
        DDUtilCreateDemod ( MOD_16QAM );
    } else if ( Modulation == MOD_64QAM ) {
        DDUtilCreateDemod ( MOD_64QAM );
    }


    for ( u32 i = 1; i < nFFT; i++ ) {
        DmodRxData4QAM [ i ] = DDUtilDemodulate ( CplxtoCplx32 ( RxDataFDFinalCorrd [ i ] ) );
        DmodRxData4QAMp1 [ i ] = DDUtilDemodulate ( CplxtoCplx32 ( RxDataFDFinalCorrdp1 [ i ] ) );
        DmodTxData4QAM [ i ] = DDUtilDemodulate ( CplxtoCplx32 ( TxDataFdEq [ i ] ) );

        if ( sampleIdx1 [ i ] >= 0 ) {
            nMean++;
            if ( DmodRxData4QAM [ i ].word != DmodTxData4QAM [ i ].word )
                errCtr += HammingDist ( DmodRxData4QAM [ i ].word, DmodTxData4QAM [ i ].word );
            errCtrp1 += HammingDist ( DmodRxData4QAMp1 [ i ].word, DmodTxData4QAM [ i ].word );
        }
    }
    for ( u32 i = 0; i < nFFT; i++ ) {
        SYSPTR( pRxDataFdCorFT[ i ] )= RxDataFDFinalCorrd [ i ];
        SYSPTR( pRxDataTdCorFT[ i ] ) = RxDataTDFinalCorr [ i ];
    }
    Ber4QAM = ( double ) errCtr / nMean;
    Ber4QAMp1 = ( double ) errCtrp1 / nMean;
}

static void CreatePhaseErrorLUT ()
{
    for ( int i = -128; i < 0; i++ ) {
        for ( int j = -128; j < 0; j++ ) {
            LUTPE [ i + 128 ] [ j + 128 ] = -atan ( ( double ) j / ( double ) i ) + M_PI / 4.0;
        }
        LUTPE [ i + 128 ] [ 128 ] = 0;
        for ( int j = 1; j < 128; j++ ) {
            LUTPE [ i + 128 ] [ j + 128 ] = -atan ( ( double ) j / ( double ) i ) - M_PI / 4.0;
        }
    }
    for ( int j = -128; j < 128; j++ ) {
        LUTPE [ 128 ] [ j + 128 ] = 0;
    }
    for ( int i = 1; i < 128; i++ ) {
        for ( int j = -128; j < 0; j++ ) {
            LUTPE [ i + 128 ] [ j + 128 ] = -atan ( ( double ) j / ( double ) i ) - M_PI / 4.0;
        }
        LUTPE [ i + 128 ] [ 128 ] = 0;
        for ( int j = 1; j < 128; j++ ) {
            LUTPE [ i + 128 ] [ j + 128 ] = -atan ( ( double ) j / ( double ) i ) + M_PI / 4.0;
        }
    }
}

static void ReadRxData ( PrsBPA *pd , int nSymbol)
{
    u32 nFFT = SYSPAR( nFFT );
    s32 deltatime = 0;//( s32 ) DVCPAR( DeltaT ) - 2;
    Cplx32 s;
    Cplx32 sp1;
    AXISink *ps = &RxDtaAXISink;
    u32 idxMask = nFFT - 1;
    u32 corSam = 0;
    u32 samp1 = 0;
    u32 G =  SYSPAR( nCPre );

    for ( u32 i = 0; i < nFFT; i++ ) // ToDo: nSymb
            {
        corSam = U32Sext ( ps->ReadBram ( ps, (( i + deltatime ) & idxMask )+nSymbol*(nFFT)), 12 );
        samp1 = U32Sext ( ps->ReadBram ( ps, (( i + deltatime - 1 ) & idxMask )+nSymbol*(nFFT)), 12 );

        s = U32toCplx32 ( corSam );
        sp1 = U32toCplx32 ( samp1 );

        RcvDataTd [ i ] = Cplx32toCplx ( s );
        RcvDataTdp1 [ i ] = Cplx32toCplx ( sp1 );

    }
    fft ( RcvDataFd, RcvDataTd, nFFT );
    fft ( RcvDataFdp1, RcvDataTdp1, nFFT );
}

static void CopyTxData ( PrsBPA *pd , int nSymbol)
{
    u32 nFFT = SYSPAR( nFFT );
    u32 nSymbPRS = SYSPAR( nSymbPRS );
    u32 Offset = nSymbPRS * nFFT+nSymbol*nFFT;

    for ( u32 i = 0; i < nFFT; i++ ) {
        TxDataFd [ i ] = ( Cplx32toCplx ( U32toCplx32 ( SYSPTR( pTxSourceBuf[ Offset + i ] ) ) ) );
    }
}

static void ReadPRSData ( PrsBPA *pd )
{
    u32 nFFT = SYSPAR( nFFT );
    for ( int i = 0; i < nFFT; i++ )
            {
        RxPrsTD [ i ] = SYSPTR( pRxPRSTdCorFT [ i ] );
    }
}

static void CopyPRSPhase ( PrsBPA *pd )
{
    s32 j = 0;
    u32 nFFT = SYSPAR( nFFT );
    for ( u32 i = 0; i < nFFT; i++ )
            {
        u32 u = SYSPTR( pTxSourceBuf[ i ] );
        if ( u ) {
            sampleIdx1 [ i ] = i;
            sampleIdx2 [ i ] = ( i < nFFT / 2 ? i : i - nFFT );
            Cplx32 zu = U32toCplx32 ( u );
            cplx z = Cplx32toCplx ( zu );
            TxPrs [ i ] = z;
            TxPrsPhase [ i ] = carg ( z );
            j++;
        }
        else
        {
            sampleIdx1 [ i ] = -1;
            sampleIdx2 [ i ] = 0;
            TxPrs [ i ] = 0.0;
            TxPrsPhase [ i ] = 0.0;
        }
    }
    nSams = j;
}

static void CalcPhaseError ( PrsBPA *pd )
{
    u32 nFFT = SYSPAR( nFFT );
    double nFFTd = ( double ) nFFT;
    double nSamsd = ( double ) nSams;
    u32 j = 0;

    MeanX = 0.0;
    SigmaX = 0.0;
    SigmaX2 = 0.0;
    for ( u32 i = 0; i < nFFT; i++ )
            {
        if ( sampleIdx1 [ i ] >= 0 )
                {
            X [ j ] = ( double ) sampleIdx2 [ i ]; // ToDo: Check if correct
            SigmaX += X [ j ];
            j++;
        }
    }
    MeanX = SigmaX / nSamsd;

    j = 0;
    for ( u32 i = 0; i < nFFT; i++ )
            {
        if ( sampleIdx1 [ i ] >= 0 )
                {
            NormX [ j ] = X [ j ] - MeanX;
            SigmaX += NormX [ j ];
            SigmaX2 += NormX [ j ] * NormX [ j ];
            j++;
        }
    }

    ridx = 0;
    for ( u32 i = 0; i < SYS_N_CHAN_PRSBPA; i++ )
            {
        double freqCorrFact = 2.0 * M_PI * ( double ) FrqCorrFactor [ i ] / nFFTd;
        CalcCorr ( pd, i, freqCorrFact, 0.0 );
        if ( SigmaY2 [ i ] < SigmaY2 [ ridx ] ) {
            ridx = i;
        }
    }
    BPARowIndex = ridx - ( SYS_N_CHAN_PRSBPA - 1 ) / 2;

}
#define FilterSpan 15
u16 nextValue ( u16 index, u16 offset, u32 nFFT )
{
    if ( index <= ( nFFT * 3 / 8 - offset ) )
        return index + offset;
    else if ( index <= ( nFFT * 3 / 8 ) )
        return index + offset - FilterSpan*2-1;
    else if ( index < ( nFFT - offset ) )
        return index + offset;
    else return offset + index - nFFT + 1;
}
u16 prevValue ( u16 index, u16 offset, u32 nFFT )
{
    if ( index <= offset )
        return nFFT + index - offset - 1;
    else if ( index <= ( nFFT * 3 / 8 ) )
        return index - offset;
    else if ( index < ( nFFT * 5 / 8 + offset  ) )
        return index - offset + FilterSpan*2+1;
    else return index - offset;
}

static void FinalCorr ( PrsBPA *pd )
{
    u32 nFFT = SYSPAR( nFFT );
    double nFFTd = ( double ) nFFT;
    double G = (double) SYSPAR( nCPre );
    double deltafpp = DVCPAR( deltafpp );
    double deltafp = DVCPAR( deltafp );
    double deltaf = DVCPAR( DeltaF );
    double deltaftot = ( deltaf + deltafp + deltafpp );



    double freqCorrFact = 2.0 * M_PI * deltafpp / nFFTd;
    CalcCorr ( pd, SYS_N_CHAN_PRSBPA, freqCorrFact, 0.0 );
    fft ( RxPrsFDCorrDfpp, RxPrsTDCorr [ SYS_N_CHAN_PRSBPA ], nFFT );

    CalcLinReg ( pd );

    DVCPAR( theta0 ) = -DVCPAR( alpha );
    double theta0 = DVCPAR( theta0 );
    DVCPAR( lastTheta ) = theta0 + ( deltaftot * 2.0 * M_PI ) / nFFTd*(nFFTd+G);


    cplx phaseCorr = cexp ( I * theta0 );
    for ( int i = 0; i < nFFT; i++ )
            {
        RxPrsFDFinalCorr [ i ] = RxPrsFDCorrDfpp [ i ] * phaseCorr; //creal(RxPrsFDCorrDfpp[ i ])
    }
    double SigmaA = 0.0;
    double GainCorr = 0.0;
    for ( int i = 0; i < nFFT; i++ )
            {
        if ( sampleIdx1 [ i ] >= 0 ) {
            SigmaA += cabs ( RxPrsFDFinalCorr [ i ] );
        }
    }
    GainCorr = ( double ) nSams / SigmaA;
    for ( int i = 0; i < nFFT; i++ )
            {
        RxPrsFDFinalCorr [ i ] *= GainCorr; /// nFFTd;
    }
    ifft ( RxPrsTDFinalCorr, RxPrsFDFinalCorr, nFFT );

    u32 *pPRS = PRSymbPtr ( nFFT );
    cplx unfilteredComplexGain [ SYS_MAX_CARRIERS ];
    for ( int i = 0; i < nFFT; i++ ) {
        cplx z = Cplx32toCplx ( U32toCplx32 ( *pPRS++ ) ) / 32767;
        if ( sampleIdx1 [ i ] >= 0 )
            unfilteredComplexGain [ i ] = z / RxPrsFDFinalCorr [ i ]; /// nFFTd;
    }

    for ( int i = 1; i < nFFT; i++ )
            {
        if ( sampleIdx1 [ i ] >= 0 )
                {
            complexGain [ i ] = unfilteredComplexGain [ i ] / ( 2 * FilterSpan - 1 );
            for ( int j = 1; j < FilterSpan; j++ )
                    {
                complexGain [ i ] += unfilteredComplexGain [ nextValue ( i, j, nFFT ) ] / ( 2 * FilterSpan - 1 );
                complexGain [ i ] += unfilteredComplexGain [ prevValue ( i, j, nFFT ) ] / ( 2 * FilterSpan - 1 );
            }
        }
    }

}

// The PrsPHI Object is not currently used.
PrsPhase PrsPhiPhase ( PrsBPA *pd, u32 n ) {
    PrsPhi *ps = pd->pPrsPhi;
    if ( !ps->init ) {
        ps->Initialise ( pd );
    }
    return ( ps->phase [ n ] );
}

void PrsPhiInit ( PrsBPA *pd )
{
    u32 nFFT = SYSPAR( nFFT );
    PrsPhi *ps = pd->pPrsPhi;
    if ( ps->init == false )
    {
        for ( int n = 0; n < nFFT; n++ )
                {
            u32 u = SYSPTR( pTxSourceBuf [ n ] );
            if ( u == 0 ) {
                ps->phase [ n ] = PHI_UNDEF;
            } else if ( u & 0x00008000 ) {
                ps->phase [ n ] = PHI_PI;
            } else if ( u & 0x00007FFF ) {
                ps->phase [ n ] = PHI_ZERO;
            } else if ( u & 0x80000000 ) {
                ps->phase [ n ] = PHI_3_PI_BY_2;
            } else if ( u & 0x7FFF0000 ) {
                ps->phase [ n ] = PHI_PI_BY_2;
            }
        }
        ps->init = true;
    }

}

static double PrsPhaseError ( PrsBPA *pd, cplx *p, u32 i )
{
    cplx z = 0.0;
    txPrsPhi.Initialise ( pd );
    switch ( txPrsPhi.Phase ( pd, i ) )
    {
    case PHI_ZERO :
        z = p [ i ];
        break;
    case PHI_PI_BY_2 :
        z = cimag ( p [ i ] );
        break;
    case PHI_PI :
        z = I * cimag ( p [ i ] ) - creal ( p [ i ] );
        break;
    case PHI_3_PI_BY_2 :
        z = -cimag ( p [ i ] );
        break;
    default :
        return ( 0 );
    }
    return ( carg ( z ) > 0 ? carg ( z ) : 2 * M_PI + carg ( z ) );
}
// End of PrsPhi object

static void CalcCorr ( PrsBPA *pd, u32 cidx, double freqCorrFact, double timeCorrFact )
{
    u32 nFFT = SYSPAR( nFFT );
    s32 deltatime = ( s32 ) DVCPAR( DeltaT ) ;
    u32 idxMask = nFFT - 1;

    for ( u32 i = 0; i < nFFT; i++ )
            {
        cplx freqCorr = cexp ( I * freqCorrFact * ( double ) (  ( i - deltatime ) & idxMask) );
        RxPrsTDCorr [ cidx ] [ i ] = freqCorr * RxPrsTD [ i ];
    }
    fft ( RxPrsFD [ cidx ], RxPrsTDCorr [ cidx ], nFFT );

    double SYA = 0.0;
    double SYB = 0.0;
    double ErrPhase = 0.0;

    u32 j = 0;
    for ( int i = 0; i < nFFT; i++ ) {
        if ( sampleIdx1 [ i ] >= 0 ) {
            cplx timeCorr = cexp ( I * timeCorrFact * ( double ) sampleIdx2 [ i ] );
            RxPrsFDCorr [ cidx ] [ i ] = timeCorr * RxPrsFD [ cidx ] [ i ];
            RxPrsPhase [ cidx ] [ i ] = carg ( RxPrsFDCorr [ cidx ] [ i ] );
            ErrPhase = RxPrsPhase [ cidx ] [ i ] - TxPrsPhase [ i ];
            YA [ cidx ] [ j ] = ErrPhase;
            if ( YA [ cidx ] [ j ] > M_PI ) {
                YA [ cidx ] [ j ] = YA [ cidx ] [ j ] - 2.0 * M_PI;
            }
            if ( YA [ cidx ] [ j ] < -M_PI ) {
                YA [ cidx ] [ j ] = YA [ cidx ] [ j ] + 2.0 * M_PI;
            }
            SYA += YA [ cidx ] [ j ];

            YB [ cidx ] [ j ] = ErrPhase;
            if ( YB [ cidx ] [ j ] < 0 ) {
                YB [ cidx ] [ j ] = YB [ cidx ] [ j ] + 2.0 * M_PI;
            }
            if ( YB [ cidx ] [ j ] > ( 2 * M_PI ) ) {
                YB [ cidx ] [ j ] = YB [ cidx ] [ j ] - 2.0 * M_PI;
            }
            SYB += YB [ cidx ] [ j ];

            j++;
        }
    }

    MeanYA [ cidx ] = SYA / ( double ) nSams;
    MeanYB [ cidx ] = SYB / ( double ) nSams;
    for ( u32 j = 0; j < nSams; j++ ) {
        NormYA [ cidx ] [ j ] = YA [ cidx ] [ j ] - MeanYA [ cidx ];
        SigmaYA [ cidx ] += YA [ cidx ] [ j ];
        SigmaY2A [ cidx ] += NormYA [ cidx ] [ j ] * NormYA [ cidx ] [ j ];
        SigmaXYA [ cidx ] += NormYA [ cidx ] [ j ] * NormX [ j ];
        NormYB [ cidx ] [ j ] = YB [ cidx ] [ j ] - MeanYB [ cidx ];
        SigmaYB [ cidx ] += YB [ cidx ] [ j ];
        SigmaY2B [ cidx ] += NormYB [ cidx ] [ j ] * NormYB [ cidx ] [ j ];
        SigmaXYB [ cidx ] += NormYB [ cidx ] [ j ] * NormX [ j ];
    }

    for ( u32 j = 0; j < nSams; j++ ) {
        if ( SigmaY2A [ cidx ] < SigmaY2B [ cidx ] ) {
            Y [ cidx ] [ j ] = YA [ cidx ] [ j ];
            MeanY [ cidx ] = MeanYA [ cidx ];
            NormY [ cidx ] [ j ] = NormYA [ cidx ] [ j ];
        } else {
            Y [ cidx ] [ j ] = YB [ cidx ] [ j ];
            MeanY [ cidx ] = MeanYB [ cidx ];
            NormY [ cidx ] [ j ] = NormYB [ cidx ] [ j ];
        }
    }

    for ( u32 j = 0; j < nSams; j++ ) {
        SigmaY [ cidx ] += Y [ cidx ] [ j ];
        SigmaY2 [ cidx ] += NormY [ cidx ] [ j ] * NormY [ cidx ] [ j ];
        SigmaXY [ cidx ] += NormY [ cidx ] [ j ] * NormX [ j ];
    }
}

static void CalcDelta ( PrsBPA *pd, s32 ridx )
{
    u32 nFFT = SYSPAR( nFFT );
    DVCPAR( fineDeltaf )=0;
    DVCPAR( fineDeltat )=0;
    DVCPAR( dataAlpha )=0;
    DVCPAR( dataBeta )=0;

    if ( ridx == 0 )
            {
        DVCPAR( deltafpp ) = FrqCorrFactor [ 0 ];
        DVCPAR( theta0 ) = MeanY [ 0 ];
    }
    else if ( ridx == SYS_N_CHAN_PRSBPA - 1 )
            {
        DVCPAR( deltafpp ) = FrqCorrFactor [ SYS_N_CHAN_PRSBPA - 1 ];
        DVCPAR( theta0 ) = MeanY [ SYS_N_CHAN_PRSBPA - 1 ];
    }
    else
    {
        double A = SigmaY2 [ ridx - 1 ];
        double B = SigmaY2 [ ridx ];
        double C = SigmaY2 [ ridx + 1 ];
        double R = ( A - C ) / ( 2 * ( A + C - 2 * B ) );
        DVCPAR( deltafpp ) = ( FrqCorrFactor [ ridx ] + R * ( FrqCorrFactor [ ridx + 1 ] - FrqCorrFactor [ ridx ] ) )
                / ( double ) nFFT;
        Theta0CIR = -MeanY [ ( SYS_N_CHAN_PRSBPA - 1 ) / 2 ];
    }
}

static void CalcLinReg ( PrsBPA *pd )
{
    double beta = SigmaXY [ SYS_N_CHAN_PRSBPA ] / SigmaX2;
    double alpha = MeanY [ SYS_N_CHAN_PRSBPA ] - beta * SigmaX;
    DVCPAR( beta ) = beta;
    DVCPAR( deltatpp ) = beta;
    DVCPAR( alpha ) = alpha;
}

void CalcDeltas ( PrsBPA *pd )
{
    int nSymbols= SYSPAR( nSymbData   ) ;

    ClearData ( pd );
    ReadPRSData ( pd );
    CopyPRSPhase ( pd );
    CalcPhaseError ( pd );
    CalcDelta ( pd, ridx );
    FinalCorr ( pd );
    for (int i=0; i< nSymbols; i++)
    {
        ReadRxData ( pd ,i);
        CopyTxData ( pd ,i);
        FreqDataCorr ( pd ,i);
        TimeDataCorr ( pd ,i);
        NormalizeDataVal ( pd );
        FineDataCorr ( pd );
        DemodRxData ( pd );
        DumpCorr ( pd );
    }
}

void DumpCorr ( PrsBPA *pd )
{

    double BERToPrint=0;
    if (Ber4QAM< Ber4QAMp1)
    {
        BERToPrint=Ber4QAM;
        DVCPAR( dataAlpha ) = alpha;
        DVCPAR( dataBeta )  = beta;
    } else {
        BERToPrint=Ber4QAMp1;
        DVCPAR( dataAlpha ) = alphap1;
        DVCPAR( dataBeta )  = betap1;
    }
    DVCPAR( lastTheta ) = DVCPAR( lastTheta )-DVCPAR( dataAlpha );

#ifdef verbose
    cprintf( "\n" );
    cprintf( "________________________________________________________\n" );
    cprintf( "for developer use:\n" );
    //cprintf( "Sigma X            = %g\n", SigmaX  );
    //cprintf( "Mean  X            = %g\n", MeanX   );
    //cprintf( "Sigma X^2          = %g\n", SigmaX2 );
    //cprintf( "Sigma Y            = %g\n", SigmaY  [ SYS_N_CHAN_PRSBPA ] );
    cprintf( "Mean  Y            = %g\n", MeanY [ SYS_N_CHAN_PRSBPA ] );
    cprintf( "Sigma Y^2          = %g\n", SigmaY2 [ SYS_N_CHAN_PRSBPA ] );
    cprintf( "Sigma XY           = %g\n", SigmaXY [ SYS_N_CHAN_PRSBPA ] );
    cprintf( "BPA Row Index      = %i\n", toInt( BPARowIndex ) );
    cprintf( "________________________________________________________\n" );
    cprintf( "\n" );
    cprintf( "========================\n" );
    cprintf( "PRS CIR/BPA Synchronizer\n" );
    cprintf( "========================\n" );
    cprintf( "NFFT               = %u\n", toUint( SYSPAR( nFFT ) ) );
    cprintf( "NCPre              = %u\n", toUint( SYSPAR( nCPre ) ) );
    cprintf( "Frequency          = %s\n", RFFstrs [ SYSPAR( TxRxFreq ) ] );
    cprintf( "Modulation         = %s\n", Modstrs [ SYSPAR( ModType ) ] );
    cprintf( "Bandwidth          = %u%%\n", toUint( SYSPAR( BWPercent ) ) );
    cprintf( "Soft decision bits = %u\n", SYS_N_SOFT_BITS );
    cprintf( "\n" );
    cprintf( "DeltaT             = %i\n", toInt( DVCPAR( DeltaT ) ) );
    cprintf( "deltat'            = %g\n", DVCPAR( deltatp ) );
    cprintf( "deltat''           = %g\n", DVCPAR( deltatpp ) );
    cprintf( "DeltaF             = %i\n", toInt( DVCPAR( DeltaF ) ) );
    cprintf( "deltaf'            = %g\n", DVCPAR( deltafp ) );
    cprintf( "deltaf''           = %g\n", DVCPAR( deltafpp ) );
    cprintf( "theta0             = %g deg\n", DVCPAR( theta0 ) * 180.0 / M_PI );
    cprintf( "theta0CIR          = %g deg\n", Theta0CIR * 180.0 / M_PI );
    cprintf( "alpha              = %g\n", DVCPAR( alpha ) );
    cprintf( "beta               = %g\n", DVCPAR( beta ) );
    cprintf( "Combined delta f   = %g\n", DVCPAR( deltafp ) + DVCPAR( deltafpp ));
    cprintf( "Sigma Error^2      = %g\n", SigmaY2 [ SYS_N_CHAN_PRSBPA ] );
    cprintf( "data alpha         = %g\n", DVCPAR( dataAlpha ) );
    cprintf( "data beta          = %g\n", DVCPAR( dataBeta ) );
    cprintf( "BER ( 4-QAM )      = %g, %g, %g\n", Ber4QAM, Ber4QAMp1,BERToPrint);
    cprintf( "\n" );

#else
    cprintf( ",%u", toUint( SYSPAR( nFFT ) ) );
    cprintf( ",%u", toUint( SYSPAR( nCPre ) ) );
    cprintf( ",%s", RFFstrs [ SYSPAR( TxRxFreq ) ] );
    cprintf( ",%s", Modstrs [ SYSPAR( ModType ) ] );
    cprintf( ",%u%%", toUint( SYSPAR( BWPercent ) ) );
    cprintf( ",%u", SYS_N_SOFT_BITS );
    cprintf( ",%i", toInt( DVCPAR( DeltaT ) ) );
    cprintf( ",%g", DVCPAR( deltatp ) );
    cprintf( ",%g", DVCPAR( deltatpp ) );
    cprintf( ",%i", toInt( DVCPAR( DeltaF ) ) );
    cprintf( ",%g", DVCPAR( deltafp ) );
    cprintf( ",%g", DVCPAR( deltafpp ) );
    cprintf( ",%g deg", DVCPAR( theta0 ) * 180.0 / M_PI );
    cprintf( ",%g deg", Theta0CIR * 180.0 / M_PI );
    cprintf( ",%g", DVCPAR( alpha ) );
    cprintf( ",%g", DVCPAR( beta ) );
    cprintf( ",%g", DVCPAR( fineDeltaf ) );
    cprintf( ",%g", SigmaY2 [ SYS_N_CHAN_PRSBPA ] );
    cprintf( ",%g", DVCPAR( dataAlpha ) );
    cprintf( ",%g", DVCPAR( dataBeta ) );
    cprintf( ",%g,%g,%g", Ber4QAM, Ber4QAMp1, BERToPrint );
    cprintf( "\n" );
#endif
}

static void DumpCorrSymbFd ( PrsBPA *pd )
{
    DumpCplx ( RxPrsFDFinalCorr, SYSPAR( nFFT ) );
}

static void DumpCorrSymbTd ( PrsBPA *pd )
{
    DumpCplx ( RxPrsTDFinalCorr, SYSPAR( nFFT ) );
}

static void DumpCorrDataSymbFd ( PrsBPA *pd )
{
    DumpCplx ( SYSPTR( pRxDataFdCorFT ), SYSPAR( nFFT ) );
}

static void DumpCorrDataSymbTd ( PrsBPA *pd )
{
    DumpCplx ( SYSPTR( pRxDataTdCorFT ), SYSPAR( nFFT ) );
}

static void RxDataSymbolDumpRawFd ( PrsBPA *pd )
{
    DumpCplx ( RxDataFDFinalCorra, SYSPAR( nFFT ) );
}
static void RxDataSymbolDumpFCFd ( PrsBPA *pd )
{
    DumpCplx ( RxDataFDFinalCorrc, SYSPAR( nFFT ) );
}
static void RxDataSymbolDumpFTCFd ( PrsBPA *pd )
{
    DumpCplx ( RxDataFDFinalCorrd, SYSPAR( nFFT ) );
}

void DMatDump ( double *m, u32 nCol, u32 nRow, char name [ ] )
{
    if ( !nCol ) {
        DColDump ( m, nRow, name );
    } else if ( !nRow ) {
        DRowDump ( m, nCol, name );
    } else {
        cprintf( "\n" );
        if ( strlen ( name ) ) {
            cprintf( "\n%s\n", name );
        }
        double *p = m;
        for ( int i = 0; i < nRow; i++ ) {
            cprintf( "%-5i ", toInt( i ) );
            for ( int j = 0; j < nCol; j++ ) {
                cprintf( " %8.5g", *p++ );
            }
            cprintf( "\n" );
            if ( ( ( i + 1 ) % 32 == 0 ) && ( pause ( "" ) == CC_ESC ) ) {
                break;
            }
        }
    }
    cprintf( "\n" );
}

void DRowDump ( double *m, u32 nCol, char name [ ] )
{
    cprintf( "\n" );
    if ( strlen ( name ) ) {
        cprintf( "%s  ", name );
    }
    for ( int i = 0; i < nCol; i++ ) {
        cprintf( "%8.5g  ", m [ i ] );
    }
    cprintf( "\n" );
}

void DColDump ( double *m, u32 nRow, char name [ ] )
{
    cprintf( "\n" );
    if ( strlen ( name ) ) {
        cprintf( "\n%s\n", name );
    }
    for ( int i = 0; i < nRow; i++ )
            {
        cprintf( "%-5i  %8.5g\n", toInt( i ), m [ i ] );
        if ( ( ( i + 1 ) % 32 == 0 ) && ( pause ( "" ) == CC_ESC ) ) {
            break;
        }
    }
    cprintf( "\n" );
}

void CMatDump ( cplx *m, u32 nCol, u32 nRow, char name [ ] )
{
    if ( !nCol ) {
        CColDump ( m, nRow, name );
    } else if ( !nRow ) {
        CRowDump ( m, nCol, name );
    } else {
        cprintf( "\n" );
        if ( strlen ( name ) ) {
            cprintf( "\n%s\n", name );
        }
        cplx *p = m;
        for ( int i = 0; i < nRow; i++ ) {
            cprintf( "%-5i", toInt( i ) );
            for ( int j = 0; j < nCol; j++ ) {
                cprintf( "  %8.5g  %8.5g", creal ( *p ), cimag ( *p ) );
                p++;
            }
            cprintf( "\n" );
            if ( ( ( i + 1 ) % 32 == 0 ) && ( pause ( "" ) == CC_ESC ) ) {
                break;
            }
        }
    }
    cprintf( "\n" );
}

void CRowDump ( cplx *m, u32 nCol, char name [ ] )
{
    cprintf( "\n" );
    if ( strlen ( name ) ) {
        cprintf( "%s\n", name );
    }
    for ( int i = 0; i < nCol; i++ ) {
        cprintf( "%8.5g %8.5g ", creal ( m [ i ] ), cimag ( m [ i ] ) );
    }
    cprintf( "\n" );
}

void CColDump ( cplx *m, u32 nRow, char name [ ] )
{
    cprintf( "\n" );
    if ( strlen ( name ) ) {
        cprintf( "%s\n", name );
    }
    for ( int i = 0; i < nRow; i++ )
            {
        cprintf( "%-5i  %8.5g  %8.5g\n", toInt( i ), creal ( m [ i ] ), cimag ( m [ i ] ) );
        if ( ( ( i + 1 ) % 32 == 0 ) && ( pause ( "" ) == CC_ESC ) ) {
            break;
        }
    }
    cprintf( "\n" );
}

void SMatDump ( s32 *m, u32 nCol, u32 nRow, char name [ ] )
{
    if ( !nCol ) {
        SColDump ( m, nRow, name );
    } else if ( !nRow ) {
        SRowDump ( m, nCol, name );
    } else {
        cprintf( "\n" );
        if ( strlen ( name ) ) {
            cprintf( "%s\n", name );
        }
        s32 *p = m;
        for ( int i = 0; i < nRow; i++ ) {
            cprintf( "%-5i", toInt( i ) );
            for ( int j = 0; j < nCol; j++ ) {
                cprintf( "  %-10i", toInt( *p ) );
            }
            cprintf( "\n" );
            if ( ( ( i + 1 ) % 32 == 0 ) && ( pause ( "" ) == CC_ESC ) ) {
                break;
            }
        }
    }
    cprintf( "\n" );
}

void SRowDump ( s32 *m, u32 nCol, char name [ ] )
{
    cprintf( "\n" );
    if ( strlen ( name ) ) {
        cprintf( "%s ", name );
    }
    for ( int i = 0; i < nCol; i++ ) {
        cprintf( "%-10i ", toInt( m[ i ] ) );
    }
    cprintf( "\n" );
}

void SColDump ( s32 *m, u32 nRow, char name [ ] )
{
    if ( strlen ( name ) ) {
        cprintf( "\n%s", name );
    }
    for ( int i = 0; i < nRow; i++ )
            {
        cprintf( "\n%-5i %-10i", toInt( i ), toInt( m[ i ] ) );
        if ( ( ( i + 1 ) % 32 == 0 ) && ( pause ( "" ) == CC_ESC ) ) {
            break;
        }
    }
    cprintf( "\n" );
}

void UMatDump ( u32 *m, u32 nCol, u32 nRow, char name [ ] )
{
    if ( !nCol ) {
        UColDump ( m, nRow, name );
    } else if ( !nRow ) {
        URowDump ( m, nCol, name );
    } else {
        if ( strlen ( name ) ) {
            cprintf( "\n%s", name );
        }
        u32 *p = m;
        for ( int i = 0; i < nRow; i++ ) {
            cprintf( "\n%-5i ", toInt( i ) );
            for ( int j = 0; j < nCol; j++ ) {
                cprintf( " %-10u", toUint( *p++ ) );
            }
            if ( ( ( i + 1 ) % 32 == 0 ) && ( pause ( "" ) == CC_ESC ) ) {
                break;
            }
        }
    }
    cprintf( "\n" );
}

void URowDump ( u32 *m, u32 nCol, char name [ ] )
{
    if ( strlen ( name ) ) {
        cprintf( "\n%s", name );
    }
    for ( int i = 0; i < nCol; i++ ) {
        cprintf( "%-10u ", toUint( m[ i ] ) );
    }
    cprintf( "\n" );
}

void UColDump ( u32 *m, u32 nRow, char name [ ] )
{
    if ( strlen ( name ) ) {
        cprintf( "\n%s", name );
    }
    for ( int i = 0; i < nRow; i++ )
            {
        cprintf( "\n%-5u %-10u", toInt( i ), toUint( m[ i ] ) );
        if ( ( ( i + 1 ) % 32 == 0 ) && ( pause ( "" ) == CC_ESC ) ) {
            break;
        }
    }
    cprintf( "\n" );
}

void DumpChanRes ( shell *psh )
{
    cprintf( "\n" );
    DRowDump ( MeanY, SYS_N_CHAN_PRSBPA, "Mean Y   " );
    DRowDump ( SigmaY, SYS_N_CHAN_PRSBPA, "Sigma Y  " );
    DRowDump ( SigmaY2, SYS_N_CHAN_PRSBPA, "Sigma Y2 " );
    DRowDump ( SigmaXY, SYS_N_CHAN_PRSBPA, "Sigma XY " );
    cprintf( "\nMean  X      : %g", MeanX );
    cprintf( "\nSigma X      : %g", SigmaX );
    cprintf( "\nSigma X^2    : %g", SigmaX2 );
    cprintf( "\nRow Idx      : %u", toUint( BPARowIndex ) );
    cprintf( "\n" );
}

