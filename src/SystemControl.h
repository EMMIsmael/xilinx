/*
 * SystemControl.h
 *
 *  Created on: 30 Jul 2018
 *      Author: JULIAN MORTIMER
 */

#ifndef SYSTEMCONTROL_H_
#define SYSTEMCONTROL_H_

#define COM_STAT_NAME_LENGTH ( 8 )

#define TX_AXISOURCE_BRAM_SIZE            ( 65536 )
#define RX_AXISINK_TDD_BRAM_SIZE          ( 65536 )
#define RX_AXISINK_PRS_BRAM_SIZE          (  8192 )
#define RX_AXISINK_DTA_BRAM_SIZE          ( 65536 )
#define RX_AXISINK_ABS00_BRAM_SIZE        (  8192 )
#define RX_AXISINK_ABS01_BRAM_SIZE        (  8192 )
#define RX_AXISINK_ABS15_BRAM_SIZE        (  8192 )
#define RX_AXISINK_NSD_FFT_BRAM_SIZE      ( 65536 )
#define RX_AXISINK_NSD_CORDIC_BRAM_SIZE   ( 65536 )
#define RX_AXISINK_NSD_BCAST00_BRAM_SIZE  (  8192 )
#define RX_AXISINK_NSD_BCAST01_BRAM_SIZE  (  8192 )
#define RX_AXISINK_NSD_BCAST15_BRAM_SIZE  (  8192 )
#define PRS_CIR_BRAM_SIZE                 (  2048 )
#define PRS_BPA_BRAM_SIZE                 (  2048 )


extern u32 TsSourceBuf[];
extern u32 TxTSInjTsBuf[];

typedef struct sSysCtrl   SysCtrl;

//================================
//================================
// IP Primitives Common Structures
//================================
//================================
typedef struct sSysDesc   SysDesc;
typedef struct sSysPtrs   SysPtrs;
typedef struct sSysData   SysData;
typedef struct sSysDvce   SysDvce;



//=====================
// System Descriptor
//=====================
struct sSysDesc {
	u32      nCPre;
	u32      nFFT;
	u32      nSymbAGC;
	u32      nSymbNSD;
    u32      nSymbPRS;
	u32      nSymbData;
    u32      nFrameBlock;
	u32      TxRxFreq;
	u32      BWPercent;
	ModType  ModType;
	SoftBits SoftBits;
};

//=====================
// System Pointers
//=====================
struct sSysPtrs {
	u32     *pTxSourceBuf;
	u32     *pTsBuf;
	u32     *pRFChanTREmulLUT;
	cplx    *pRxPRSTdCorDeltaT;
	cplx    *pRxPRSFdCorFT;
    cplx    *pRxPRSTdCorFT;
    cplx    *pRxDataTdCorFT;
    cplx    *pRxDataFdCorFT;
};


//==========================
// System Tx / Rx Stats etc.
//==========================
struct sSysData {
	void (* Clear  ) ( SysData *ps );
	void (* Dump   ) ( SysData *ps );
	void (* Update ) ( SysData *ps );
	SysDvce *pSysDvce;

	u32 RxDataCtrReg;
	u32 RxSymbCtrReg;
	u32 RxFrameCtrReg;
	u32 RxBlockCtrReg;
	u32 TxDataCtrReg;
	u32 TxSymbCtrReg;
	u32 TxFrameCtrReg;
	u32 TxBlockCtrReg;
    u32 RxTVCtrReg;
    u32 TxTVCtrReg;

	u32 RxDataCtr;
	u32 RxSymbCtr;
	u32 RxFrameCtr;
	u32 RxBlockCtr;
	u32 TxDataCtr;
	u32 TxSymbCtr;
	u32 TxFrameCtr;
	u32 TxBlockCtr;
	u32 RxTVCtr;
    u32 TxTVCtr;
};

//================================
//================================
// IP Primitives Device Structures
//================================
//================================
//================
// FFT Controller
//================

typedef struct sFFTCDesc FFTCDesc;
struct sFFTCDesc {
	_Bool    nFFTFixed;
	u32      nFFT;
	u32      nFFTMin;
	u32      nFFTMax;
	_Bool    hasCP;
	u32      nCPre;
	u32      fwdInv;
	_Bool    isScaled;
	u32      scaleSch[ FFTC_SCALE_CHAN ];
	u32      cfgRegLen;
	u32      hwCfgRegLen;
	u32      cfgFFTLoReg;
	u32      cfgFFTHiReg;
	_Bool    idxBitReversed;
	_Bool    hasXkIdx;
};

typedef struct sFFTCData FFTCData;
struct sFFTCData {
    u32 dummy;
};

typedef struct sFFTC  FFTC;
struct sFFTC {
	SysDvce   *pSysDvce;
	FFTCDesc  *pDvcDesc;
	FFTCData  *pDvcData;
	u32  (* HwCfgRegLen   ) ( FFTC *pd );
	void (* CfgFFT        ) ( FFTC *pd );
	void (* Dump          ) ( FFTC *pd );
};

//================================
// RF Channel Emulator and Switch
//================================
typedef struct sTREmulDesc TREmulDesc;
struct sTREmulDesc {
	u32 *pLUTf;
    u32 PsLUTSize;
    u32 HwLUTSize;
};

typedef struct sTREmulData TREmulData;
struct sTREmulData {
    u32 dummy;
};

typedef struct sTREmul  TREmul;
struct sTREmul {
	SysDvce     *pSysDvce;
	TREmulDesc  *pDvcDesc;
	TREmulData  *pDvcData;
	void (* Start           ) ( TREmul *pd );
	void (* Stop            ) ( TREmul *pd );
	u32  (* TestBram        ) ( TREmul *pd );
	u32  (* ReadPsLUT       ) ( TREmul *pd, u32 addr );
	void (* WritePsLUT      ) ( TREmul *pd, u32 addr, u32 data );
	void (* InitPsLUT       ) ( TREmul *pd, u32 addr, u32 count, u32 data );
	void (* LoadPsLUT       ) ( TREmul *pd );
	void (* DumpPsLUT       ) ( TREmul *pd, u32 addr, u32 nBits );
	u32  (* ReadHwLUT       ) ( TREmul *pd, u32 addr );
	void (* WriteHwLUT      ) ( TREmul *pd, u32 addr, u32 data );
	void (* InitHwLUT       ) ( TREmul *pd, u32 addr, u32 count, u32 data );
	void (* LoadHwLUT       ) ( TREmul *pd );
	void (* DumpHwLUT       ) ( TREmul *pd, u32 addr, u32 nBits );
	void (* Dump            ) ( TREmul *pd );
};

typedef struct sTREmulSwDesc TREmulSwDesc;
struct sTREmulSwDesc {
	u32 ConfigReg;
};

typedef struct sTREmulSwData TREmulSwData;
struct sTREmulSwData {
    u32   sel_emul;
    _Bool fifo_clk_lock;
    _Bool s_axis_in_tvalid;
    _Bool s_axis_in_tready;
    _Bool m_axis_out_tvalid;
    _Bool m_axis_out_tready;
    _Bool s_axis_device_out_tvalid;
    _Bool s_axis_device_out_tready;
    _Bool m_axis_device_in_tvalid;
    _Bool m_axis_device_in_tready;
    _Bool s_axis_emul_out_tvalid;
    _Bool s_axis_emul_out_tready;
    _Bool m_axis_emul_in_tvalid;
    _Bool m_axis_emul_in_tready;
};



typedef struct sTREmulSw  TREmulSw;
struct sTREmulSw {
	SysDvce       *pSysDvce;
	TREmulSwDesc  *pDvcDesc;
	TREmulSwData  *pDvcData;
	void (* SelChan          ) ( TREmulSw *pd, u32 SelEmul );
	void (* Dump             ) ( TREmulSw *pd );
};

//======================
// AXI4 BRAM Controller
//======================
typedef struct sAXIBramCDesc AXIBramCDesc;
struct sAXIBramCDesc {
	u32 BramSize;
};

typedef struct sAXIBramCData AXIBramCData;
struct sAXIBramCData {
    u32 dummy;
};


typedef struct sAXIBramC  AXIBramC;
struct sAXIBramC {
	SysDvce       *pSysDvce;
	AXIBramCDesc  *pDvcDesc;
	AXIBramCData  *pDvcData;
	u32  (* TestBram         ) ( AXIBramC *pd );
	u32  (* ReadBram         ) ( AXIBramC *pd, u32 addr );
	void (* WriteBram        ) ( AXIBramC *pd, u32 addr, u32 data );
	void (* InitBram         ) ( AXIBramC *pd, u32 addr, u32 count, u32 data );
	void (* LoadBram         ) ( AXIBramC *pd );
	u32  (* HwBramSize       ) ( AXIBramC *pd );
	void (* Dump             ) ( AXIBramC *pd );
};

//==================
// AXI4 Data Source
//==================
typedef struct sAXISourceDesc AXISourceDesc;
struct sAXISourceDesc {
	u32             TxAddr;
	u32             BlockSize;
	u32             BramSize;
	u32             BramInitVal;
	AXIBramC        *pBramC;
};

typedef struct sAXISourceData AXISourceData;
struct sAXISourceData {
   _Bool            AxsTxFrameDone;
   _Bool            AxsMAxisTready;
   _Bool            AxsMAxisTvalid;
   u32              AxsTfState;
   u32              AxsTmState;
};

typedef struct sAXISource  AXISource;
struct sAXISource {
	SysDvce        *pSysDvce;
	AXISourceDesc  *pDvcDesc;
	AXISourceData  *pDvcData;
	void (* Start1Shot        ) ( AXISource *pd );
	void (* StartContin       ) ( AXISource *pd );
	void (* Stop              ) ( AXISource *pd );
    void (* SetTxAddr         ) ( AXISource *pd, u32 addr );
    void (* SetBlockSize      ) ( AXISource *pd, u32 size );
	void (* TxStream          ) ( AXISource *pd, u32 addr, u32 nwords );
	u32  (* ReadBram          ) ( AXISource *pd, u32 addr );
	void (* WriteBram         ) ( AXISource *pd, u32 addr, u32 data );
	void (* InitBram          ) ( AXISource *pd, u32 addr, u32 count, u32 data );
	void (* LoadBram          ) ( AXISource *pd );
	void (* DumpBram          ) ( AXISource *pd, u32 addr, u32 nBits );
	u32  (* HwBramSize        ) ( AXISource *pd );
	void (* Dump              ) ( AXISource *pd );
};

//==================
// Frequency Counter
//==================
typedef struct sFreqCtrDesc FreqCtrDesc;
struct sFreqCtrDesc {
    u32 dummy;
};

typedef struct sFreqCtrData FreqCtrData;
struct sFreqCtrData {
    u32 dummy;
};


typedef struct sFreqCtr  FreqCtr;
struct sFreqCtr {
	SysDvce      *pSysDvce;
	FreqCtrDesc  *pDvcDesc;
	FreqCtrData  *pDvcData;
	u32    (* ReadCtrReg    ) ( FreqCtr *pd );
	u32    (* ReadTmrReg    ) ( FreqCtr *pd );
	double (* ReadFreqCtr   ) ( FreqCtr *pd );
	void   (* Dump          ) ( FreqCtr *pd );
};

//=========================
// Training Symbol Injector
//=========================
typedef struct sTSInjDesc TSInjDesc;
struct sTSInjDesc {
    u32 TsBufSize;
    u32 TxBufSize;
    u32 TsBufInitVal;
    u32 RxBufInitVal;
    u32 TSInjDoutLsb;
};

typedef struct sTSInjData TSInjData;
struct sTSInjData {
	u32 RxDsAccum;
	u32 TxDsAccum;
	u32 TxTsAccum;
	u32 RxFrameCount;
    u32 TxFrameCount;
	u32 NWdSymbRx;
    u32 NWdSymbTx;
    u32 NSymbFrameRx;
    u32 NSymbFrameTx;
    u32 NSymbTrain;
	_Bool TsiRxSymbolDone;
	_Bool TsiTxSymbolDone;
	_Bool TsiSAxisDataTready;
	_Bool TsiSAxisDataTvalid;
	_Bool TsiMAxisDataTready;
	_Bool TsiMAxisDataTvalid;
	u32 TsiExState;
    u32 TsiPfState;
    u32 TsiRmState;
    u32 TsiTmState;
    u32 TsiAdRdTsState;
    u32 TsiAdRdTxState;
};


typedef struct sTSInj TSInj;
struct sTSInj {
	SysDvce    *pSysDvce;
	TSInjDesc  *pDvcDesc;
	TSInjData  *pDvcData;
	void (* Start           ) ( TSInj *pd );
	void (* Stop            ) ( TSInj *pd );
	u32  (* TestBram        ) ( TSInj *pd );
	u32  (* ReadTsBuf       ) ( TSInj *pd, u32 addr );
	void (* WriteTsBuf      ) ( TSInj *pd, u32 addr, u32 data );
	void (* InitTsBuf       ) ( TSInj *pd, u32 addr, u32 count, u32 val );
	void (* LoadTsBuf       ) ( TSInj *pd );
	void (* DumpTsBuf       ) ( TSInj *pd, u32 addr, u32 nBits );
	u32  (* ReadTxBuf       ) ( TSInj *pd, u32 addr );
	void (* WriteRxBuf      ) ( TSInj *pd, u32 addr, u32 data );
	void (* InitRxBuf       ) ( TSInj *pd, u32 addr, u32 count, u32 val );
	void (* DumpTxBuf       ) ( TSInj *pd, u32 addr, u32 nBits );
	u32  (* ReadPsBuf       ) ( TSInj *pd, u32 addr );
	void (* DumpPsBuf       ) ( TSInj *pd, u32 addr, u32 nBits );
	void (* Dump            ) ( TSInj *pd );
};

//======================
// NULL Symbol Detector
//======================
typedef struct sNSDetDesc NSDetDesc;
struct sNSDetDesc {
	u32   RNoise;
	u32   FrqIdxRef;
	u32   NsdDetThr;
	u32   ChanDataSel;
	u32   NsdTstMode;
	u32   NsdFrqBlk    [ NSDET_N_BLK_DET ];
	u32   SampleMapTbl [ NSDET_MAX_FREQ  ];
	u32   SampleMapIdx [ NSDET_MAX_FREQ  ];
	u32   SampleMap    [ NSDET_N_FREQ    ];
	u32   SampleMapFrq [ NSDET_N_FREQ    ];
	u32   TsFreq       [ SYS_N_TSYMB_MAX ];
};


typedef struct sNSDetData NSDetData;
struct sNSDetData {
	u64   NsdDtaRxCtr;
	u64   NsdBlkIdx;
	u32   NsdBlkRxCtr;
	u32   NsdSumNoiseMin;
	u32   NsdSumPeaksMax;
	u32   NsdDetCtr;
	u32   NsdChnDetMap;
    u32   NsdSymbDetCtr;
	u64   NsdChnFrqIdx   [ NSDET_N_CHAN ];
	u32   NsdChnSumNoise [ NSDET_N_CHAN ];
	u32   NsdChnSumPeaks [ NSDET_N_CHAN ];
	u32   NsdChnDetHist  [ NSDET_N_CHAN ];
	u32   NsdChnDb0      [ NSDET_N_CHAN ];
	u32   NsdChnDb1      [ NSDET_N_CHAN ];
	u32   NsdChnDb2      [ NSDET_N_CHAN ];
	u32   NsdChnDb3      [ NSDET_N_CHAN ];
    u32   NsdChnDb4      [ NSDET_N_CHAN ];
    u32   NsdChnDb5      [ NSDET_N_CHAN ];
    u32   NsdChnDb6      [ NSDET_N_CHAN ];
    u32   NsdChnDb7      [ NSDET_N_CHAN ];
    u32   NsdChnDb8      [ NSDET_N_CHAN ];
    u32   NsdChnDb9      [ NSDET_N_CHAN ];
    u32   NsdChnDb10     [ NSDET_N_CHAN ];
    u32   NsdChnDb11     [ NSDET_N_CHAN ];
    u32   NsdSdState     [ NSDET_N_CHAN ];
    u32   NsdTdState     [ NSDET_N_CHAN ];
	_Bool NsdFFTCfgDone;
    _Bool NsdDetDone;
    _Bool NsdTxdTxDone;
    _Bool NsdSAxisRxdTvalid;
    _Bool NsdSAxisRxdTready;
    _Bool NsdMAxisFFTTvalid;
    _Bool NsdMAXisFFTTready;
    _Bool NsdMAxisPrsTvalid;
    _Bool NsdMAXisPrsTready;
    _Bool NsdMAxisDtaTvalid;
    _Bool NsdMAXisDtaTready;
    _Bool NsdBufAbsTvalidF;
    _Bool NsdSysResetn;
    _Bool NsdFifoResetn;
    _Bool NsdFFTResetn;
    _Bool NsdFFTCfgRqst;
    u32   NsdExState;
    u32   NsdNsdState;
    u32   NsdTxdState;
	u32   NsdSAxisAbsTready;
	u32   NsdSAxisAbsTvalid;
    u32   NsdBufAbsTready;
    u32   NsdBufAbsTvalid;
    u32   NsdPsSampleMap0;
    u32   NsdPsSampleMap1;
    u32   NsdPsSampleMap2;
    u32   NsdPsSampleMap3;
    u32   NsdExTimer;

    u64   NsdSymbDtaIdx;
    u32   NsdChanDetIdx;
    u32   NsdSumNoiseIdx;
    u32   NsdSumPeaksIdx;
    u64   NsdRxRxdCtr;
    u64   NsdRxFFTCtr;
    u64   NsdTxFFTCtr;
    u32   NsdTxPrsCtr;
    u32   NsdTxDtaCtr;
    u32   NsdTxSymCtr;
    u64   NsdTxPrsStart;
    u64   NsdTxDtaStart;
};


typedef struct sNSDet NSDet;
struct sNSDet {
	SysDvce    *pSysDvce;
	NSDetDesc  *pDvcDesc;
	NSDetData  *pDvcData;
	void (* Start           ) ( NSDet *pd );
	void (* Stop            ) ( NSDet *pd );
	void (* SetTstMode      ) ( NSDet *pd, NSDetTstMode m );
	void (* SetRNoise       ) ( NSDet *pd, double rmean );
    void (* SyncDataRead    ) ( NSDet *pd );
    void (* SetDbgTready    ) ( NSDet *pd, _Bool b );
	void (* Dstate          ) ( NSDet *pd );
	void (* Dump            ) ( NSDet *pd );
};

//=======================
// PRS-CIR Calculation
//=======================
//============================
// Common Descriptor Structure
//============================
typedef struct sPrsDesc PrsDesc;
struct sPrsDesc {
    u32    CIRBramSize;
    s32    DeltaT;
    s32    DeltaF;
    double deltafp;
    double deltatp;
    double deltafpp;
    double deltatpp;
    double alpha;
    double beta;
    double theta0;
    double lastTheta;
    double dataAlpha;
    double dataBeta;
    double fineDeltaf;
    double fineDeltat;
    double KFF[ SYS_N_CARRIER_CFG ];
    double KTF[ SYS_N_CARRIER_CFG ];
    double KTT[ SYS_N_CARRIER_CFG ];
    double KFT[ SYS_N_CARRIER_CFG ];
};

//=================
// PRS - CIR Table:
//=================
typedef struct sPrsCIRData PrsCIRData;
struct sPrsCIRData {
    u32 RxDataCtr;
    u32 RxSymbCtr;
    u32 RxFrameCtr;
    u32 RxBlockCtr;
    u32 TxDataCtr;
    u32 TxSymbCtr;
    u32 TxFrameCtr;
    u32 TxBlockCtr;
};

typedef struct sPrsCIR PrsCIR;
struct sPrsCIR {
    SysDvce        *pSysDvce;
    PrsDesc        *pDvcDesc;
    PrsCIRData     *pDvcData;
    u32  (* TestBram                 ) ( PrsCIR *pd );
    u64  (* ReadBram                 ) ( PrsCIR *pd, u32 addr );
    void (* WriteBram                ) ( PrsCIR *pd, u32 addr, u64 data );
    u32  (* ReadBramFFT              ) ( PrsCIR *pd, u32 addr );
    void (* WriteBramFFT             ) ( PrsCIR *pd, u32 addr, u32 data );
    void (* InitBramFFT              ) ( PrsCIR *pd, u32 addr, u32 count, u32 data );
    void (* DumpBramFFT              ) ( PrsCIR *pd, u32 addr, u32 nBits );
    u32  (* ReadBramABS              ) ( PrsCIR *pd, u32 addr );
    void (* WriteBramABS             ) ( PrsCIR *pd, u32 addr, u32 data );
    void (* InitBramABS              ) ( PrsCIR *pd, u32 addr, u32 count, u32 data );
    void (* DumpBramABS              ) ( PrsCIR *pd, u32 addr, u32 nBits );
    void (* Dump                     ) ( PrsCIR *pd );
    void (* Calc                     ) ( PrsCIR *pd );
    void (* Recalc                   ) ( PrsCIR *pd );
    void (* DumpCIRTable             ) ( PrsCIR *pd );
    void (* DumpRxPRSTd              ) ( PrsCIR *pd );
    void (* DumpRxPRSFd              ) ( PrsCIR *pd );
    void (* DumpTxPRSFd              ) ( PrsCIR *pd );
    void (* DumpRxPRSTdCorDeltaT     ) ( PrsCIR *pd );
    void (* DumpRxPRSTdCordeltafp    ) ( PrsCIR *pd );
    void (* DumpRxPRSFdCordeltafp    ) ( PrsCIR *pd );
    void (* DumpRxPRSFdCordeltatp    ) ( PrsCIR *pd );
    void (* DumpRxPRSTdCorFT         ) ( PrsCIR *pd );

};

//=======================
// PRS-BPA Calculation
//=======================
typedef struct sPrsBPAData PrsBPAData;
struct sPrsBPAData {
    u32 TxDataCtr;
    u32 TxSymbCtr;
    u32 TxFrameCtr;
    u32 TxBlockCtr;
};

typedef struct sprsphi PrsPhi;

typedef struct sPrsBPA PrsBPA;
struct sPrsBPA {
    SysDvce        *pSysDvce;
    PrsDesc        *pDvcDesc;
    PrsBPAData     *pDvcData;
    void (* Dump           )    ( PrsBPA *pd );
    void (* DumpCorr       )    ( PrsBPA *pd );
    void (* DumpCorrSymbFd )    ( PrsBPA *pd );
    void (* DumpCorrSymbTd )    ( PrsBPA *pd );
    void (* DumpCorrDataSymbFd) ( PrsBPA *pd );
    void (* DumpCorrDataSymbTd) ( PrsBPA *pd );
    void (* CalcDeltas     ) ( PrsBPA *pd );
    void (* RxDataSymbolDumpRawFd     ) ( PrsBPA *pd );
    void (* RxDataSymbolDumpFCFd     ) ( PrsBPA *pd );
    void (* RxDataSymbolDumpFTCFd      ) ( PrsBPA *pd );
    PrsPhi         *pPrsPhi;
};


//=======================
// AXI4-Stream Data Sink
//=======================
typedef struct sAXISinkDesc AXISinkDesc;
struct sAXISinkDesc {
	u32 RxAddr;
	u32 BlockSize;
	u32 BramSize;
	u32 BramInitVal;
};

typedef struct sAXISinkData AXISinkData;
struct sAXISinkData {
	u32 WrDataCtr;
	u32 WrSymbCtr;
	u32 WrFrameCtr;
	u32 WrBlockCtr;
	_Bool rd_data_rdy;
	_Bool rx_stream_done;
	_Bool s_axis_data_tready;
	_Bool s_axis_data_tvalid;
	u32 ex_state;
	u32 rx_state;
	u32 rd_data_state;
};


typedef struct sAXISink  AXISink;
struct sAXISink {
	SysDvce      *pSysDvce;
	AXISinkDesc  *pDvcDesc;
	AXISinkData  *pDvcData;
	void (* Start1Shot      ) ( AXISink *pd );
	void (* StartContin     ) ( AXISink *pd );
	void (* Stop            ) ( AXISink *pd );
	void (* SetRxAddr       ) ( AXISink *pd, u32 addr );
	void (* SetBlockSize    ) ( AXISink *pd, u32 size );
	void (* WriteInhibit    ) ( AXISink *pd );
	void (* WriteEnable     ) ( AXISink *pd );
	void (* RxStream        ) ( AXISink *pd, u32 addr, u32 nwords );
	u32  (* TestBram        ) ( AXISink *pd );
	u32  (* ReadBram        ) ( AXISink *pd, u32 addr );
	void (* WriteBram       ) ( AXISink *pd, u32 addr, u32 data );
	void (* InitBram        ) ( AXISink *pd, u32 addr, u32 count, u32 data );
	void (* DumpBram        ) ( AXISink *pd, u32 addr, u32 nBits );
	u32  (* HwBramSize      ) ( AXISink *pd );
	void (* Dump            ) ( AXISink *pd );
};

struct sSysCtrl {
	SysDesc        *pSysDesc;
    SysPtrs        *pSysPtrs;
	AXIBramC       *pTxSourceAXIBramC;
	AXISource      *pTxAXISource;
	FFTC           *pTxFFTC;
	FFTC           *pRxFFTC;
	FFTC           *pRxNsdFFTC;
	FFTC           *pRxPrsFFTC;
	FFTC           *pRxCirFFTC;
	FFTC           *pRxBpaFFTC;
	TSInj          *pTxTSInj;
	TREmul         *pRFChanTREmul;
	TREmulSw       *pRFChanTREmulSw;
	FreqCtr        *pRFChanTREmulFreqCtr;
	FreqCtr        *pAD9364ClkFreqCtr;
	FreqCtr        *pAD9364DACFreqCtr;
	FreqCtr        *pAD9364ADCFreqCtr;
	NSDet          *pRxNSDet;
	AXISink        *pRxTddAXISink;
    AXISink        *pRxPrsAXISink;
    AXISink        *pRxDtaAXISink;
    AXISink        *pRxNsdB00AXISink;
    AXISink        *pRxNsdB01AXISink;
    AXISink        *pRxNsdB15AXISink;
    AXISink        *pRxNsdFFTAXISink;
    AXISink        *pRxNsdCordicAXISink;
    AXISink        *pRxNsdA00AXISink;
    AXISink        *pRxNsdA01AXISink;
    AXISink        *pRxNsdA15AXISink;
	PrsCIR         *pRxPrsCIR;
    PrsBPA         *pRxPrsBPA;
};


extern FFTCDesc TxFFTDesc;
extern FFTCDesc NsdFFTDesc;
extern FFTCDesc RxFFTDesc;
extern FFTCDesc PrsFFTDesc;
extern FFTCDesc CirFFTDesc;
extern FFTCDesc BpaFFTDesc;

#endif /* SYSTEMCONTROL_H_ */
