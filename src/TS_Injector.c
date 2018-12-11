/*
 * TS_Injector.c
 *
 *  Created on: 6 Jul 2017
 *      Author: JULIAN MORTIMER
 */
#include "include.h"


// =================================
// Boilerplate Function Declarations
// =================================
void TSInjInit(
    TSInj *pd,
    const char      *Name,
    u32             BaseAddr,
    SysCtrl         *pSysCtrl,
    SysDvce         *pSysDvce,
    SysDesc         *pSysDesc,
    SysPtrs         *pSysPtrs,
    SysData         *pSysData,
    TSInjDesc *pDvcDesc,
    TSInjData *pDvcData
    );
static void SetParams    ( TSInj *pd );
static void SysDvceDump  ( TSInj *pd );
static void SysDescDump  ( TSInj *pd );
static void SysDataDump  ( TSInj *pd );
static void DvcDataInit  ( TSInj *pd );
static void DvcDataUpd   ( TSInj *pd );
static void DvcDescDump  ( TSInj *pd );
static void DvcDataDump  ( TSInj *pd );

// =================================================
// Device Structure Member Functions (public access)
// =================================================
static inline void Start        ( TSInj *pd ) { hwWriteCtrlReg( SYSDVC( BaseAddr ), TSINJ_CTRL_REG << 2, TSINJ_START ); }
static inline void Stop         ( TSInj *pd ) { hwWriteCtrlReg( SYSDVC( BaseAddr ), TSINJ_CTRL_REG << 2, TSINJ_STOP  ); }
static        u32  TestBram     ( TSInj *pd );
static        u32  ReadTsBuf    ( TSInj *pd, u32 addr );
static        void WriteTsBuf   ( TSInj *pd, u32 addr, u32 data );
static        void InitTsBuf    ( TSInj *pd, u32 addr, u32 count, u32 val );
static        void LoadTsBuf    ( TSInj *pd );
static        void DumpTsBuf    ( TSInj *pd, u32 addr, u32 nBits );
static        u32  ReadTxBuf    ( TSInj *pd, u32 addr );
static        void WriteRxBuf   ( TSInj *pd, u32 addr, u32 data );
static        void InitRxBuf    ( TSInj *pd, u32 addr, u32 count, u32 val );
static        void DumpTxBuf    ( TSInj *pd, u32 addr, u32 nBits );
static        u32  ReadPsBuf    ( TSInj *pd, u32 addr );
static        void DumpPsBuf    ( TSInj *pd, u32 addr, u32 nBits );
static        void Dump         ( TSInj *pd );


static void DumpBram( TSInj *pd, u32 addr, u32 nbits, dmpSelector dsel );


// Declare init structure here with member functions
TSInj TSInjInitStr = {
    NULL,
    NULL,
    NULL,
    Start,
    Stop,
    TestBram,
    ReadTsBuf,
    WriteTsBuf,
    InitTsBuf,
    LoadTsBuf,
    DumpTsBuf,
    ReadTxBuf,
    WriteRxBuf,
    InitRxBuf,
    DumpTxBuf,
    ReadPsBuf,
    DumpPsBuf,
    Dump
};

//===================
// Initialize Device
//===================
void TSInjInit( TSInj *pd, const char *Name, u32 BaseAddr, SysCtrl *pSysCtrl, SysDvce *pSysDvce, SysDesc *pSysDesc, SysPtrs *pSysPtrs, SysData *pSysData, TSInjDesc *pDvcDesc, TSInjData *pDvcData )
{
    *pd             = TSInjInitStr;
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
        TSINJ_HW_INFO_REG,
        TSINJ_CTRL_REG,
        TSINJ_INFO_REG,
        TSINJ_STAT_REG,
        TSINJ_RESET,
        TSINJ_INT_ACK,
        TSINJ_RD_SYNC_F,
        TSINJ_RX_DATA_CTR_REG,
        TSINJ_RX_SYMBOL_CTR_REG,
        TSINJ_RX_FRAME_CTR_REG,
        TSINJ_RX_BLOCK_CTR_REG,
        TSINJ_TX_DATA_CTR_REG,
        TSINJ_TX_SYMBOL_CTR_REG,
        TSINJ_TX_FRAME_CTR_REG,
        TSINJ_TX_BLOCK_CTR_REG,
        TSINJ_RX_TV_CTR_REG,
        TSINJ_TX_TV_CTR_REG
        );
    pd->pDvcDesc  = pDvcDesc;
    pd->pDvcData  = pDvcData;
    DvcDataInit( pd );
    SYSRESET( pd );
   	( void ) pd->TestBram( pd );
	pd->InitTsBuf( pd, 0, DVCPAR( TsBufSize ), DVCPAR( TsBufInitVal ) );
	pd->InitRxBuf( pd, 0, DVCPAR( TxBufSize ), DVCPAR( RxBufInitVal ) );
    SetParams( pd );
}

//========================
// Dump Device to Console
//========================
static void Dump( TSInj *pd ) {
    cprintf("\n"); fflush( stdout );
	SysDvceDump( pd );
	SysDescDump( pd );
	DvcDescDump( pd );
	SysDataDump( pd );
	DvcDataDump( pd );
	cprintf("\n"); fflush( stdout );
}

// End of boilerplate code
//========================
// Device-specific code
//========================
static void SetParams( TSInj *pd )
{
    SYSWRSLVREG( pd, TSINJ_N_CPRE_REG        , SYSPAR( nCPre        ));
    SYSWRSLVREG( pd, TSINJ_N_FFT_REG         , SYSPAR( nFFT         ));
    SYSWRSLVREG( pd, TSINJ_N_SYMB_AGC_REG    , SYSPAR( nSymbAGC     ));
    SYSWRSLVREG( pd, TSINJ_N_SYMB_NSD_REG    , SYSPAR( nSymbNSD     ));
    SYSWRSLVREG( pd, TSINJ_N_SYMB_PRS_REG    , SYSPAR( nSymbPRS     ));
    SYSWRSLVREG( pd, TSINJ_N_SYMB_DTA_REG    , SYSPAR( nSymbData    ));
    SYSWRSLVREG( pd, TSINJ_N_FRAME_BLOCK_REG , SYSPAR( nFrameBlock  ));
}

static void DvcDataInit( TSInj *pd )
{
	DVCDAT( RxDsAccum ) = 0;
	DVCDAT( TxDsAccum ) = 0;
	DVCDAT( TxTsAccum ) = 0;
	DVCDAT( NWdSymbTx ) = 0;
}

static void DvcDataUpd( TSInj *pd )
{
	DVCDAT( RxDsAccum          ) = SYSRDSLVREG( pd, TSINJ_RX_DS_ACCUM_REG     );
	DVCDAT( TxDsAccum          ) = SYSRDSLVREG( pd, TSINJ_TX_DS_ACCUM_REG     );
	DVCDAT( TxTsAccum          ) = SYSRDSLVREG( pd, TSINJ_TX_TS_ACCUM_REG     );
    DVCDAT( RxFrameCount       ) = SYSRDSLVREG( pd, TSINJ_RX_FRAME_COUNT_REG  );
    DVCDAT( TxFrameCount       ) = SYSRDSLVREG( pd, TSINJ_TX_FRAME_COUNT_REG  );
	DVCDAT( NWdSymbRx          ) = SYSRDSLVREG( pd, TSINJ_N_WD_SYMB_RX_REG    );
    DVCDAT( NWdSymbTx          ) = SYSRDSLVREG( pd, TSINJ_N_WD_SYMB_TX_REG    );
    DVCDAT( NSymbFrameRx       ) = SYSRDSLVREG( pd, TSINJ_N_SYMB_FRAME_RX_REG );
    DVCDAT( NSymbFrameTx       ) = SYSRDSLVREG( pd, TSINJ_N_SYMB_FRAME_TX_REG );
    DVCDAT( NSymbTrain         ) = SYSRDSLVREG( pd, TSINJ_N_SYMB_TRAIN_REG    );
	u32 statreg                  = SYSRDSLVREG( pd, TSINJ_STAT_REG            );
    DVCDAT( TsiRxSymbolDone    ) = ( statreg & TSINJ_RX_SYMBOL_DONE_MASK      ) ? true : false;
    DVCDAT( TsiTxSymbolDone    ) = ( statreg & TSINJ_TX_SYMBOL_DONE_MASK      ) ? true : false;
	DVCDAT( TsiSAxisDataTready ) = ( statreg & TSINJ_S_AXIS_DATA_TREADY_MASK  ) ? true : false;
    DVCDAT( TsiSAxisDataTvalid ) = ( statreg & TSINJ_S_AXIS_DATA_TVALID_MASK  ) ? true : false;
    DVCDAT( TsiMAxisDataTready ) = ( statreg & TSINJ_M_AXIS_DATA_TREADY_MASK  ) ? true : false;
    DVCDAT( TsiMAxisDataTvalid ) = ( statreg & TSINJ_M_AXIS_DATA_TVALID_MASK  ) ? true : false;
    u32 statereg                 = SYSRDSLVREG( pd, TSINJ_STATE_REG           );
    DVCDAT( TsiExState         ) = ( statereg & TSINJ_EX_STATE_MASK           ) >> TSINJ_EX_STATE_SHR;
    DVCDAT( TsiPfState         ) = ( statereg & TSINJ_PF_STATE_MASK           ) >> TSINJ_PF_STATE_SHR;
    DVCDAT( TsiRmState         ) = ( statereg & TSINJ_RM_STATE_MASK           ) >> TSINJ_RM_STATE_SHR;
    DVCDAT( TsiTmState         ) = ( statereg & TSINJ_TM_STATE_MASK           ) >> TSINJ_TM_STATE_SHR;
    DVCDAT( TsiAdRdTsState     ) = ( statereg & TSINJ_AD_RD_TS_STATE_MASK     ) >> TSINJ_AD_RD_TS_STATE_SHR;
    DVCDAT( TsiAdRdTxState     ) = ( statereg & TSINJ_AD_RD_TX_STATE_MASK     ) >> TSINJ_AD_RD_TX_STATE_SHR;
}

static const char s_EX_UNK0             [ ] = "UNKNOWN";
static const char s_EX_UNK1             [ ] = "UNKNOWN";
static const char s_EX_UNK2             [ ] = "UNKNOWN";
static const char s_EX_UNK3             [ ] = "UNKNOWN";
static const char s_EX_UNK4             [ ] = "UNKNOWN";
static const char s_EX_UNK5             [ ] = "UNKNOWN";
static const char s_EX_UNK6             [ ] = "UNKNOWN";
static const char s_EX_UNK7             [ ] = "UNKNOWN";
static const char s_EX_UNK8             [ ] = "UNKNOWN";
static const char s_EX_IDLE             [ ] = "EX_IDLE";
static const char s_EX_ADMIN            [ ] = "EX_ADMIN";
static const char s_EX_UNKB             [ ] = "UNKNOWN";
static const char s_EX_UNKC             [ ] = "UNKNOWN";
static const char s_EX_RUNNING          [ ] = "EX_RUNNING";
static const char s_EX_UNKE             [ ] = "UNKNOWN";
static const char s_EX_UNKF             [ ] = "UNKNOWN";

static const char *ExStateStr [ 16 ] = {
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
         s_EX_ADMIN,
         s_EX_UNKB,
         s_EX_UNKC,
         s_EX_RUNNING,
         s_EX_UNKE,
         s_EX_UNKF
};

static const char s_PF_UNK0             [ ] = "UNKNOWN";             // 0
static const char s_PF_UNK1             [ ] = "UNKNOWN";             // 1
static const char s_PF_UNK2             [ ] = "UNKNOWN";             // 2
static const char s_PF_START_FRAME      [ ] = "PF_START_FRAME";      // 3
static const char s_PF_UNK4             [ ] = "UNKNOWN";             // 4
static const char s_PF_WAIT_RX_DS       [ ] = "PF_WAIT_RX_DS";       // 5
static const char s_PF_UNK6             [ ] = "UNKNOWN";             // 6
static const char s_PF_UNK7             [ ] = "UNKNOWN";             // 7
static const char s_PF_UNK8             [ ] = "UNKNOWN";             // 8
static const char s_PF_IDLE             [ ] = "PF_IDLE";             // 9
static const char s_PF_UNKA             [ ] = "UNKNOWN";             // A
static const char s_PF_UNKB             [ ] = "UNKNOWN";             // B
static const char s_PF_WAIT_TX_TS_START [ ] = "PF_WAIT_TX_TS_START"; // C
static const char s_PF_WAIT_TX_TS       [ ] = "PF_WAIT_TX_TS";       // D
static const char s_PF_WAIT_TX_DS_START [ ] = "PF_WAIT_TX_DS_START"; // E
static const char s_PF_WAIT_TX_DS       [ ] = "PF_WAIT_TX_DS";       // F

static const char *PfStateStr [ 16 ] = {
         s_PF_UNK0,
         s_PF_UNK1,
         s_PF_UNK2,
         s_PF_START_FRAME,
         s_PF_UNK4,
         s_PF_WAIT_RX_DS,
         s_PF_UNK6,
         s_PF_UNK7,
         s_PF_UNK8,
         s_PF_IDLE,
         s_PF_UNKA,
         s_PF_UNKB,
         s_PF_WAIT_TX_TS_START,
         s_PF_WAIT_TX_TS,
         s_PF_WAIT_TX_DS_START,
         s_PF_WAIT_TX_DS
};

static const char s_RM_UNK0             [ ] = "UNKNOWN";
static const char s_RM_UNK1             [ ] = "UNKNOWN";
static const char s_RM_UNK2             [ ] = "UNKNOWN";
static const char s_RM_UNK3             [ ] = "UNKNOWN";
static const char s_RM_UNK4             [ ] = "UNKNOWN";
static const char s_RM_UNK5             [ ] = "UNKNOWN";
static const char s_RM_UNK6             [ ] = "UNKNOWN";
static const char s_RM_UNK7             [ ] = "UNKNOWN";
static const char s_RM_UNK8             [ ] = "UNKNOWN";
static const char s_RM_IDLE             [ ] = "RM_IDLE";
static const char s_RM_UNKA             [ ] = "UNKNOWN";
static const char s_RM_RUNNING          [ ] = "RM_RUNNING";
static const char s_RM_UNKC             [ ] = "UNKNOWN";
static const char s_RM_UNKD             [ ] = "UNKNOWN";
static const char s_RM_UNKE             [ ] = "UNKNOWN";
static const char s_RM_UNKF             [ ] = "UNKNOWN";

static const char *RmStateStr [ 16 ] = {
         s_RM_UNK0,
         s_RM_UNK1,
         s_RM_UNK2,
         s_RM_UNK3,
         s_RM_UNK4,
         s_RM_UNK5,
         s_RM_UNK6,
         s_RM_UNK7,
         s_RM_UNK8,
         s_RM_IDLE,
         s_RM_UNKA,
         s_RM_RUNNING,
         s_RM_UNKC,
         s_RM_UNKD,
         s_RM_UNKE,
         s_RM_UNKF
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

static const char s_AD_RD_TS_UNK0       [ ] = "UNKNOWN";
static const char s_AD_RD_TS_UNK1       [ ] = "UNKNOWN";
static const char s_AD_RD_TS_UNK2       [ ] = "UNKNOWN";
static const char s_AD_RD_TS_UNK3       [ ] = "UNKNOWN";
static const char s_AD_RD_TS_UNK4       [ ] = "UNKNOWN";
static const char s_AD_RD_TS_UNK5       [ ] = "UNKNOWN";
static const char s_AD_RD_TS_UNK6       [ ] = "UNKNOWN";
static const char s_AD_RD_TS_UNK7       [ ] = "UNKNOWN";
static const char s_AD_RD_TS_UNK8       [ ] = "UNKNOWN";
static const char s_AD_RD_TS_IDLE       [ ] = "AD_RD_TS_IDLE";
static const char s_AD_RD_TS_UNKA       [ ] = "UNKNOWN";
static const char s_AD_RD_TS_WAIT_BRAM  [ ] = "AD_RD_TS_WAIT_BRAM";
static const char s_AD_RD_TS_UNKC       [ ] = "UNKNOWN";
static const char s_AD_RD_TS_UNKD       [ ] = "UNKNOWN";
static const char s_AD_RD_TS_UNKE       [ ] = "UNKNOWN";
static const char s_AD_RD_TS_RDY        [ ] = "AD_RD_TS_RDY";

static const char *AdRdTsStateStr [ 16 ] = {
         s_AD_RD_TS_UNK0,
         s_AD_RD_TS_UNK1,
         s_AD_RD_TS_UNK2,
         s_AD_RD_TS_UNK3,
         s_AD_RD_TS_UNK4,
         s_AD_RD_TS_UNK5,
         s_AD_RD_TS_UNK6,
         s_AD_RD_TS_UNK7,
         s_AD_RD_TS_UNK8,
         s_AD_RD_TS_IDLE,
         s_AD_RD_TS_UNKA,
         s_AD_RD_TS_WAIT_BRAM,
         s_AD_RD_TS_UNKC,
         s_AD_RD_TS_UNKD,
         s_AD_RD_TS_UNKE,
         s_AD_RD_TS_RDY
};

static const char s_AD_RD_TX_UNK0       [ ] = "UNKNOWN";
static const char s_AD_RD_TX_UNK1       [ ] = "UNKNOWN";
static const char s_AD_RD_TX_UNK2       [ ] = "UNKNOWN";
static const char s_AD_RD_TX_UNK3       [ ] = "UNKNOWN";
static const char s_AD_RD_TX_UNK4       [ ] = "UNKNOWN";
static const char s_AD_RD_TX_UNK5       [ ] = "UNKNOWN";
static const char s_AD_RD_TX_UNK6       [ ] = "UNKNOWN";
static const char s_AD_RD_TX_UNK7       [ ] = "UNKNOWN";
static const char s_AD_RD_TX_UNK8       [ ] = "UNKNOWN";
static const char s_AD_RD_TX_IDLE       [ ] = "AD_RD_TX_IDLE";
static const char s_AD_RD_TX_UNKA       [ ] = "UNKNOWN";
static const char s_AD_RD_TX_WAIT_BRAM  [ ] = "AD_RD_TX_WAIT_BRAM";
static const char s_AD_RD_TX_UNKC       [ ] = "UNKNOWN";
static const char s_AD_RD_TX_UNKD       [ ] = "UNKNOWN";
static const char s_AD_RD_TX_UNKE       [ ] = "UNKNOWN";
static const char s_AD_RD_TX_RDY        [ ] = "AD_RD_TX_RDY";

static const char *AdRdTxStateStr [ 16 ] = {
         s_AD_RD_TX_UNK0,
         s_AD_RD_TX_UNK1,
         s_AD_RD_TX_UNK2,
         s_AD_RD_TX_UNK3,
         s_AD_RD_TX_UNK4,
         s_AD_RD_TX_UNK5,
         s_AD_RD_TX_UNK6,
         s_AD_RD_TX_UNK7,
         s_AD_RD_TX_UNK8,
         s_AD_RD_TX_IDLE,
         s_AD_RD_TX_UNKA,
         s_AD_RD_TX_WAIT_BRAM,
         s_AD_RD_TX_UNKC,
         s_AD_RD_TX_UNKD,
         s_AD_RD_TX_UNKE,
         s_AD_RD_TX_RDY
};


static void DvcDataDump( TSInj *pd )
{
	DvcDataUpd( pd );
	cprintf("\nRxDsAccum           : %u", toUint( DVCDAT( RxDsAccum               ) ) );
	cprintf("\nTxDsAccum           : %u", toUint( DVCDAT( TxDsAccum               ) ) );
	cprintf("\nTxTsAccum           : %u", toUint( DVCDAT( TxTsAccum               ) ) );
    cprintf("\nRxFrameCount        : %u", toUint( DVCDAT( RxFrameCount            ) ) );
    cprintf("\nTxFrameCount        : %u", toUint( DVCDAT( TxFrameCount            ) ) );
	cprintf("\nnWdSymbRx           : %u", toUint( DVCDAT( NWdSymbRx               ) ) );
    cprintf("\nnWdSymbTx           : %u", toUint( DVCDAT( NWdSymbTx               ) ) );
    cprintf("\nnSymbFrameRx        : %u", toUint( DVCDAT( NSymbFrameRx            ) ) );
    cprintf("\nnSymbFrameTx        : %u", toUint( DVCDAT( NSymbFrameTx            ) ) );
    cprintf("\nnSymbTrain          : %u", toUint( DVCDAT( NSymbTrain              ) ) );
    cprintf("\n");
    cprintf("\nrx_symbol_done      : %s", OFFON(  DVCDAT( TsiRxSymbolDone         ) ) );
    cprintf("\ntx_symbol_done      : %s", OFFON(  DVCDAT( TsiTxSymbolDone         ) ) );
    cprintf("\n");
    cprintf("\ns_axis_data_tready  : %s", OFFON(  DVCDAT( TsiSAxisDataTready      ) ) );
    cprintf("\ns_axis_data_tvalid  : %s", OFFON(  DVCDAT( TsiSAxisDataTvalid      ) ) );
    cprintf("\nm_axis_data_tready  : %s", OFFON(  DVCDAT( TsiMAxisDataTready      ) ) );
    cprintf("\nm_axis_data_tvalid  : %s", OFFON(  DVCDAT( TsiMAxisDataTvalid      ) ) );
    cprintf("\n");
	cprintf("\nex_state            : %s", ExStateStr     [ DVCDAT( TsiExState     ) ] );
    cprintf("\npf_state            : %s", PfStateStr     [ DVCDAT( TsiPfState     ) ] );
    cprintf("\nrm_state            : %s", RmStateStr     [ DVCDAT( TsiRmState     ) ] );
    cprintf("\ntm_state            : %s", TmStateStr     [ DVCDAT( TsiTmState     ) ] );
    cprintf("\nad_rd_ts_state      : %s", AdRdTsStateStr [ DVCDAT( TsiAdRdTsState ) ] );
    cprintf("\nad_rd_tx_state      : %s", AdRdTxStateStr [ DVCDAT( TsiAdRdTxState ) ] );
	cprintf("\n");
}

static void DvcDescDump( TSInj *pd )
{
    cprintf("\nTsBufSize        : %u"  , toUint( DVCPAR( TsBufSize    ) ) );
    cprintf("\nTxBufSize        : %u"  , toUint( DVCPAR( TxBufSize    ) ) );
    cprintf("\nTsBufInitVal     : %08X", toUint( DVCPAR( TsBufInitVal ) ) );
    cprintf("\nRxBufInitVal     : %08X", toUint( DVCPAR( RxBufInitVal ) ) );
    cprintf("\nTSInjDoutLsb     : %u"  , toUint( DVCPAR( TSInjDoutLsb ) ) );
    cprintf("\n");
}

static void SysDvceDump( TSInj *pd )
{
	SYSDVCEDUMP( pd );
}

static void SysDataDump( TSInj *pd )
{
	SYSDATADUMP( pd );
}

static void SysDescDump( TSInj *pd )
{
	cprintf("\nnFFT             : %u",   toUint( SYSPAR( nFFT        )) );
	cprintf("\nnCPre            : %u",   toUint( SYSPAR( nCPre       )) );
    cprintf("\nnWdSymbFrameRx   : %u",   toUint( SYSPAR( nFFT        )) );
    cprintf("\nnWdSymbFrameTx   : %u",   toUint( SYSPAR( nFFT ) + SYSPAR( nCPre )) );
	cprintf("\nnSymbFrameRx     : %u",   toUint( SYSPAR( nSymbPRS ) + SYSPAR( nSymbData )) );
    cprintf("\nnSymbFrameTx     : %u",   toUint( SYSPAR( nSymbAGC ) + SYSPAR( nSymbNSD ) + SYSPAR( nSymbPRS ) + SYSPAR( nSymbData )) );
	cprintf("\nnFrameBlock      : %u",   toUint( SYSPAR( nFrameBlock )) );
	cprintf("\nnSymbAGC         : %u",   toUint( SYSPAR( nSymbAGC    )) );
	cprintf("\nnSymbNSD         : %u",   toUint( SYSPAR( nSymbNSD    )) );
    cprintf("\nnSymbPRS         : %u",   toUint( SYSPAR( nSymbPRS    )) );
    cprintf("\nnSymbDTA         : %u",   toUint( SYSPAR( nSymbData   )) );
	cprintf("\n");
}



// Insert device-specific code here
//======================================
// Device-specific non member functions
//======================================
static inline void EnterAdmin    ( TSInj *pd ) { hwWriteCtrlReg( SYSDVC( BaseAddr ), TSINJ_CTRL_REG << 2, TSINJ_ENTER_ADMIN ); }
static inline void ExitAdmin     ( TSInj *pd ) { hwWriteCtrlReg( SYSDVC( BaseAddr ), TSINJ_CTRL_REG << 2, TSINJ_EXIT_ADMIN  ); }

//=================================
// Device-specific member functions
//=================================
static void DumpTsBuf( TSInj *pd, u32 addr, u32 nbits )
{
	DumpBram( pd, addr, nbits, DUMP_SEL_TSINJ_TSBUF );
}

static void DumpTxBuf( TSInj *pd, u32 addr, u32 nbits )
{
	DumpBram( pd, addr, nbits, DUMP_SEL_TSINJ_TXBUF );
}
static void DumpPsBuf( TSInj *pd, u32 addr, u32 nbits )
{
	DumpBram( pd, addr, nbits, DUMP_SEL_TSINJ_PSBUF );
}

static void DumpBram( TSInj *pd, u32 addr, u32 nbits, dmpSelector dsel )
{
	char s[ SYS_MAX_OBJ_NAME_LEN + 1  ]; strcpy( s, SYSDVC( Name ) );
	switch( dsel ) {
	  case DUMP_SEL_TSINJ_PSBUF :
		  strcat( s, " Ps Buffer" );
		  break;
	  case DUMP_SEL_TSINJ_TXBUF :
		  strcat( s, " Tx Buffer" );
		  break;
	  case DUMP_SEL_TSINJ_TSBUF :
	  default:
		  dsel = DUMP_SEL_TSINJ_TSBUF;
		  strcat( s, " Ts Buffer" );
		  break;
	};

	if ( nbits == 0 )
		utilScreenDumpData( dsel, DUMP_FMT_HEX,       DUMP_TYP_MEMORY, addr, 4,     0, ( u32 )pd, s );
	else if(( nbits >= 8 ) && ( nbits <= 16 ))
     	utilScreenDumpData( dsel, DUMP_FMT_IQ_SIGNED, DUMP_TYP_MEMORY, addr, nbits, 0, ( u32 )pd, s );
	else
		cprintf( "Format must be 0 (hex) or 8-16 (decimal Q:I)\n" );
}


static u32 ReadTsBuf( TSInj *pd, u32 addr )
{
	EnterAdmin( pd );
    u32 u = hwReadData( SYSDVC( BaseAddr ),
                   TSINJ_CTRL_REG << 2, TSINJ_AD_RD_TS_BUF, TSINJ_AD_RD_ACK,
                   TSINJ_STAT_REG << 2, TSINJ_AD_RD_RDY,
                   TSINJ_AD_RD_WR_ADDR_REG << 2, addr,
                   TSINJ_AD_RD_TS_BUF_REG << 2
                   );
    ExitAdmin( pd );
    return ( u );
}

static void WriteTsBuf( TSInj *pd, u32 addr, u32 data )
{
	EnterAdmin( pd );
    hwWriteData( SYSDVC( BaseAddr ),
    		    TSINJ_CTRL_REG << 2, TSINJ_AD_WR_TS_BUF,
    		    TSINJ_AD_RD_WR_ADDR_REG << 2, addr,
    		    TSINJ_AD_WR_DATA_REG << 2, data
    		    );
    ExitAdmin( pd );
}

static void InitTsBuf( TSInj *pd, u32 addr, u32 count, u32 val )
{
    EnterAdmin( pd );
    for( int i = 0; i < count; i++ )
	{
	    hwWriteData( SYSDVC( BaseAddr ),
	    		    TSINJ_CTRL_REG << 2, TSINJ_AD_WR_TS_BUF,
	    		    TSINJ_AD_RD_WR_ADDR_REG << 2, addr + i,
	    		    TSINJ_AD_WR_DATA_REG << 2, val
	    		    );
	}
    ExitAdmin( pd );
}


static void LoadTsBuf( TSInj *pd )
{
	EnterAdmin( pd );
	for ( int i = 0; i < SYS_MAX_CARRIERS * SYS_N_TSYMB_MAX; i++ )
	{
		hwWriteData( SYSDVC( BaseAddr ),
				TSINJ_CTRL_REG << 2, TSINJ_AD_WR_TS_BUF,
				TSINJ_AD_RD_WR_ADDR_REG << 2, i,
				TSINJ_AD_WR_DATA_REG << 2, SYSPTR( pTsBuf[ i ] )
		);
	}
	ExitAdmin( pd );
}

static u32 ReadTxBuf( TSInj *pd, u32 addr )
{
    EnterAdmin( pd );
    u32 u = hwReadData( SYSDVC( BaseAddr ),
                   TSINJ_CTRL_REG << 2, TSINJ_AD_RD_TX_BUF, TSINJ_AD_RD_ACK,
                   TSINJ_STAT_REG << 2, TSINJ_AD_RD_RDY,
                   TSINJ_AD_RD_WR_ADDR_REG << 2, addr,
                   TSINJ_AD_RD_TX_BUF_REG << 2
                   );
    ExitAdmin( pd );
    return( u );
}

static void WriteRxBuf( TSInj *pd, u32 addr, u32 data )
{
	EnterAdmin( pd );
    hwWriteData( SYSDVC( BaseAddr ),
    		    TSINJ_CTRL_REG << 2, TSINJ_AD_WR_RX_BUF,
    		    TSINJ_AD_RD_WR_ADDR_REG << 2, addr,
    		    TSINJ_AD_WR_DATA_REG << 2, data
    		    );
    ExitAdmin( pd );
}

static void InitRxBuf( TSInj *pd, u32 addr, u32 count, u32 val )
{
    EnterAdmin( pd );
	for( int i = 0; i < count; i++ )
	{
	    hwWriteData( SYSDVC( BaseAddr ),
	    		    TSINJ_CTRL_REG << 2, TSINJ_AD_WR_RX_BUF,
	    		    TSINJ_AD_RD_WR_ADDR_REG << 2, addr + i,
	    		    TSINJ_AD_WR_DATA_REG << 2, val
	    		    );
	}
    ExitAdmin( pd );
}

static u32 ReadPsBuf( TSInj *pd, u32 addr )
{
    return ( SYSPTR( pTsBuf[ addr ] ));
}


static u32 ReadTxB( u32 addr )
{
	static TSInj *pd = NULL;
	if ( addr & 0x80000000 ) { pd = ( TSInj * ) ( addr & 0x7FFFFFFF ); }
	else if ( pd != NULL ) { return ( pd->ReadTxBuf( pd, ( u32 )addr )); }
	return( 0 );
}

static void WriteTxB( u32 addr, u32 data )
{
	static TSInj *pd = NULL;
	if ( addr & 0x80000000 ) { pd = ( TSInj * )( addr & 0x7FFFFFFF ); }
	else if ( pd != NULL ) { pd->WriteRxBuf( pd, ( u32 )addr, data ); }
}

static u32 ReadTsB( u32 addr )
{
	static TSInj *pd = NULL;
	if ( addr & 0x80000000 ) { pd = ( TSInj * ) ( addr & 0x7FFFFFFF ); }
	else if ( pd != NULL ) { return ( pd->ReadTsBuf( pd, ( u32 )addr )); }
	return( 0 );
}

static void WriteTsB( u32 addr, u32 data )
{
	static TSInj *pd = NULL;
	if ( addr & 0x80000000 ) { pd = ( TSInj * )( addr & 0x7FFFFFFF ); }
	else if ( pd != NULL ) { pd->WriteTsBuf( pd, ( u32 )addr, data ); }
}


static u32 TestTxBuf( TSInj *pd )
{
	char Name[ 2 * SYS_MAX_OBJ_NAME_LEN + 1]; strcpy( Name, SYSDVC( Name ) );
//	u32 Read( u32 addr ) { return( pd->ReadTxBuf( pd, addr )); }
//	void Write( u32 addr, u32 data ) { pd->WriteRxBuf( pd, addr, data   ); }
	return ( UtilTestBram( ReadTxB, WriteTxB, DVCPAR( TxBufSize ), 0xFFFFFFFF, strcat( Name, " tx/rx buffer" )));
}

static u32 TestTsBuf( TSInj *pd )
{
	char Name[ 2 * SYS_MAX_OBJ_NAME_LEN + 1]; strcpy( Name, SYSDVC( Name ) );
//	u32 Read( u32 addr ) { return( pd->ReadTsBuf( pd, addr )); }
//	void Write( u32 addr, u32 data ) { pd->WriteTsBuf( pd, addr, data   ); }
	return ( UtilTestBram( ReadTsB, WriteTsB, DVCPAR( TsBufSize ), 0xFFFFFFFF, strcat( Name, " ts buffer" ) ));
}

static u32 TestBram( TSInj *pd )
{
    ( void ) ReadTsB( ( u32 ) pd | 0x80000000 );
    ( void ) ReadTxB( ( u32 ) pd | 0x80000000 );
    WriteTsB( ( u32 ) pd | 0x80000000, 0 );
    WriteTxB( ( u32 ) pd | 0x80000000, 0 );
    ( void ) TestTsBuf( pd );
	( void ) TestTxBuf( pd );
	return( 0 );
}


