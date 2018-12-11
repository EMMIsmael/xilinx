/*
 * Null_Symbol_Detector.c
 *
 *  Created on: 20 Jun 2018
 *      Author: JULIAN MORTIMER
 */
#include "include.h"




// =================================
// Boilerplate Function Declarations
// =================================
void NSDetInit(
		NSDet *pd,
		const char      *Name,
		u32             BaseAddr,
		SysCtrl         *pSysCtrl,
		SysDvce         *pSysDvce,
		SysDesc         *pSysDesc,
		SysPtrs         *pSysPtrs,
		SysData         *pSysData,
		NSDetDesc *pDvcDesc,
		NSDetData *pDvcData
);
static void SetParams    ( NSDet *pd );
static void SysDvceDump  ( NSDet *pd );
static void SysDescDump  ( NSDet *pd );
static void SysDataDump  ( NSDet *pd );
static void DvcDataInit  ( NSDet *pd );
static void DvcDataUpd   ( NSDet *pd );
static void DvcDescDump  ( NSDet *pd );
static void DvcDataDump  ( NSDet *pd );

// =================================================
// Device Structure Member Functions (public access)
// =================================================
static inline void Start         ( NSDet *pd )                 { hwWriteCtrlReg( SYSDVC( BaseAddr ), NSDET_CTRL_REG << 2, NSDET_START );   }
static inline void Stop          ( NSDet *pd )                 { hwWriteCtrlReg( SYSDVC( BaseAddr ), NSDET_CTRL_REG << 2, NSDET_STOP );    }
static inline void SetTstMode    ( NSDet *pd, NSDetTstMode m ) { DVCPAR( NsdTstMode ) = m; SYSWRSLVREG( pd, NSDET_TST_MODE_REG, m ); }
static inline void SetRNoise     ( NSDet *pd, double rmean )   { DVCPAR( RNoise ) = RMeantoRNoise( rmean, NSDET_N_BITS_SCALE_SUM_BLOCK ); SYSWRSLVREG( pd, NSDET_R_NOISE_REG, DVCPAR( RNoise ) ); }
static inline void SyncDataRead  ( NSDet *pd )                 { hwWriteCtrlReg( SYSDVC( BaseAddr ), NSDET_CTRL_REG << 2, NSDET_SYNC );    }
static inline void SetDbgTready  ( NSDet *pd, _Bool b )        { u32 r = SYSRDSLVREG( pd, NSDET_GPO_REG ); r &= ~NSDET_S_AXIS_DATA_TREADY_G; r |= ( b ? NSDET_S_AXIS_DATA_TREADY_G : 0 ); SYSWRSLVREG( pd, NSDET_GPO_REG, r ); }
static void        Dstate        ( NSDet *pd );
static void        Dump          ( NSDet *pd );

// Declare init structure here with member functions
NSDet NSDetInitStr = {
		NULL,
		NULL,
		NULL,
		Start,
		Stop,
		SetTstMode,
		SetRNoise,
		SyncDataRead,
		SetDbgTready,
		Dstate,
		Dump
};


static inline void SelChanData( NSDet *pd, u32 n ) { if( n < NSDET_N_CHAN ) { DVCPAR( ChanDataSel ) = n; SYSWRSLVREG( pd, NSDET_CHAN_DATA_SEL_REG, DVCPAR( ChanDataSel )); SyncDataRead( pd ); }}
static void MkTSBuf( NSDet *pd );
static void MkFrqRef(NSDet *pd);


//===================
// Initialize Device
//===================
void NSDetInit( NSDet *pd, const char *Name, u32 BaseAddr, SysCtrl *pSysCtrl, SysDvce *pSysDvce, SysDesc *pSysDesc, SysPtrs *pSysPtrs, SysData *pSysData, NSDetDesc *pDvcDesc, NSDetData *pDvcData )
{
	*pd             = NSDetInitStr;
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
			NSDET_HW_INFO_REG,
			NSDET_CTRL_REG,
			NSDET_INFO_REG,
			NSDET_STAT_REG,
			NSDET_RESET,
			NSDET_INT_ACK,
            0,
			NSDET_RX_DATA_CTR_REG,
			NSDET_RX_SYMB_CTR_REG,
			NSDET_RX_FRAME_CTR_REG,
			NSDET_RX_BLOCK_CTR_REG,
			NSDET_TX_DATA_CTR_REG,
			NSDET_TX_SYMB_CTR_REG,
			NSDET_TX_FRAME_CTR_REG,
			NSDET_TX_BLOCK_CTR_REG,
			0,
			0
	);
	pd->pDvcDesc  = pDvcDesc;
	pd->pDvcData  = pDvcData;
	DvcDataInit( pd );
	SYSRESET( pd );
    MkTSBuf( pd );
    MkFrqRef( pd );
	SetParams( pd );
}

//========================
// Dump Device to Console
//========================
static void Dump( NSDet *pd ) {
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
static void SetParams( NSDet *pd )
{
	SYSWRSLVREG( pd, NSDET_SAMPLE_MAP_0_REG,        DVCPAR( SampleMap[ 0 ] ));
	SYSWRSLVREG( pd, NSDET_SAMPLE_MAP_1_REG,        DVCPAR( SampleMap[ 1 ] ));
	SYSWRSLVREG( pd, NSDET_SAMPLE_MAP_2_REG,        DVCPAR( SampleMap[ 2 ] ));
	SYSWRSLVREG( pd, NSDET_SAMPLE_MAP_3_REG,        DVCPAR( SampleMap[ 3 ] ));
	SYSWRSLVREG( pd, NSDET_N_FFT_REG,               SYSPAR( nFFT ));
	SYSWRSLVREG( pd, NSDET_N_SYMB_AGC_REG,          SYSPAR( nSymbAGC ));
    SYSWRSLVREG( pd, NSDET_N_SYMB_NSD_REG,          SYSPAR( nSymbNSD ));
    SYSWRSLVREG( pd, NSDET_N_SYMB_PRS_REG,          SYSPAR( nSymbPRS ));
    SYSWRSLVREG( pd, NSDET_N_SYMB_DTA_REG,          SYSPAR( nSymbData ));
	SYSWRSLVREG( pd, NSDET_N_FRAME_BLOCK_REG,       SYSPAR( nFrameBlock ));
	SYSWRSLVREG( pd, NSDET_N_CPRE_REG,              SYSPAR( nCPre ));
	SYSWRSLVREG( pd, NSDET_R_NOISE_REG,             DVCPAR( RNoise ));
	SYSWRSLVREG( pd, NSDET_FRQ_IDX_REF_REG,         DVCPAR( FrqIdxRef ));
	SYSWRSLVREG( pd, NSDET_DET_THR_REG,             DVCPAR( NsdDetThr ));
	SYSWRSLVREG( pd, NSDET_CHAN_DATA_SEL_REG,       DVCPAR( ChanDataSel ));
	SYSWRSLVREG( pd, NSDET_TST_MODE_REG,            DVCPAR( NsdTstMode ));
}

static void DvcDataInit( NSDet *pd )
{
	SYSWRSLVREG ( pd, NSDET_CHAN_DATA_SEL_REG, 0 );
	DVCDAT( NsdDtaRxCtr         ) = 0;
	DVCDAT( NsdBlkIdx           ) = 0;
	DVCDAT( NsdBlkRxCtr         ) = 0;
	DVCDAT( NsdSumNoiseMin      ) = 0;
	DVCDAT( NsdSumPeaksMax      ) = 0;
	DVCDAT( NsdDetCtr           ) = 0;
	DVCDAT( NsdSymbDetCtr       ) = 0;
	DVCDAT( NsdChnDetMap        ) = 0;
	for( int i = 0; i < NSDET_N_CHAN; i++ ) {
		DVCDAT( NsdChnFrqIdx  [ i ] ) = 0;
		DVCDAT( NsdChnSumNoise[ i ] ) = 0;
		DVCDAT( NsdChnSumPeaks[ i ] ) = 0;
		DVCDAT( NsdChnDetHist [ i ] ) = 0;
		DVCDAT( NsdChnDb0     [ i ] ) = 0;
		DVCDAT( NsdChnDb1     [ i ] ) = 0;
		DVCDAT( NsdChnDb2     [ i ] ) = 0;
		DVCDAT( NsdChnDb3     [ i ] ) = 0;
        DVCDAT( NsdChnDb4     [ i ] ) = 0;
        DVCDAT( NsdChnDb5     [ i ] ) = 0;
        DVCDAT( NsdChnDb6     [ i ] ) = 0;
        DVCDAT( NsdChnDb7     [ i ] ) = 0;
        DVCDAT( NsdChnDb8     [ i ] ) = 0;
        DVCDAT( NsdChnDb9     [ i ] ) = 0;
        DVCDAT( NsdChnDb10    [ i ] ) = 0;
        DVCDAT( NsdChnDb11    [ i ] ) = 0;
        DVCDAT( NsdSdState    [ i ] ) = 0;
        DVCDAT( NsdTdState    [ i ] ) = 0;
	}
    DVCDAT( NsdFFTCfgDone       ) = false;
    DVCDAT( NsdDetDone          ) = false;
    DVCDAT( NsdTxdTxDone        ) = false;
    DVCDAT( NsdSAxisRxdTvalid   ) = false;
    DVCDAT( NsdSAxisRxdTready   ) = false;
    DVCDAT( NsdMAxisFFTTvalid   ) = false;
    DVCDAT( NsdMAXisFFTTready   ) = false;
    DVCDAT( NsdMAxisPrsTvalid   ) = false;
    DVCDAT( NsdMAXisPrsTready   ) = false;
    DVCDAT( NsdMAxisDtaTvalid   ) = false;
    DVCDAT( NsdMAXisDtaTready   ) = false;
    DVCDAT( NsdBufAbsTvalidF    ) = false;
    DVCDAT( NsdSysResetn        ) = false;
    DVCDAT( NsdFifoResetn       ) = false;
    DVCDAT( NsdFFTResetn        ) = false;
    DVCDAT( NsdFFTCfgRqst       ) = false;
    DVCDAT( NsdExState          ) = 0;
    DVCDAT( NsdNsdState         ) = 0;
    DVCDAT( NsdTxdState         ) = 0;
    DVCDAT( NsdSAxisAbsTready   ) = 0;
    DVCDAT( NsdSAxisAbsTvalid   ) = 0;
    DVCDAT( NsdBufAbsTready     ) = 0;
    DVCDAT( NsdBufAbsTvalid     ) = 0;
    DVCDAT( NsdPsSampleMap0     ) = 0;
    DVCDAT( NsdPsSampleMap1     ) = 0;
    DVCDAT( NsdPsSampleMap2     ) = 0;
    DVCDAT( NsdPsSampleMap3     ) = 0;
    DVCDAT( NsdExTimer          ) = 0;
    DVCDAT( NsdSymbDtaIdx       ) = 0;
    DVCDAT( NsdChanDetIdx       ) = 0;
    DVCDAT( NsdSumNoiseIdx      ) = 0;
    DVCDAT( NsdSumPeaksIdx      ) = 0;
    DVCDAT( NsdRxRxdCtr         ) = 0;
    DVCDAT( NsdRxFFTCtr         ) = 0;
    DVCDAT( NsdTxFFTCtr         ) = 0;
    DVCDAT( NsdTxPrsCtr         ) = 0;
    DVCDAT( NsdTxDtaCtr         ) = 0;
    DVCDAT( NsdTxSymCtr         ) = 0;
    DVCDAT( NsdTxPrsStart       ) = 0;
    DVCDAT( NsdTxDtaStart       ) = 0;
}


void u32toFidx( char s[], u32 u )
{
    u32 uu = u;
    char c = 0;
    s[ 16 ] = '\0';
    for ( int i = 15; i >= 0; i-- )
    {
        switch ( uu & 3 )
        {
        case 0 :
            c = '0';
            break;
        case 1 :
            c = '1';
            break;
        case 2:
            c = '2';
            break;
        case 3:
            c = '3';
            break;
        default:
            c = '?';
            break;
        }
        s[ i ] = c;
        uu >>= 2;
    }
}


static void DvcDataUpd( NSDet *pd )
{
    SyncDataRead( pd );

    u32 StatReg                       = SYSRDSLVREG   ( pd, NSDET_STAT_REG );
	DVCDAT( NsdBlkIdx           )     = SYSRDSLVREG64 ( pd, NSDET_NSD_BLK_IDX_LO_REG );
	DVCDAT( NsdBlkRxCtr         )     = SYSRDSLVREG   ( pd, NSDET_NSD_BLK_RX_CTR_REG    );
	DVCDAT( NsdSumNoiseMin      )     = SYSRDSLVREG   ( pd, NSDET_NSD_SUM_NOISE_MIN_REG );
	DVCDAT( NsdSumPeaksMax      )     = SYSRDSLVREG   ( pd, NSDET_NSD_SUM_PEAKS_MAX_REG );
	DVCDAT( NsdDetCtr           )     = SYSRDSLVREG   ( pd, NSDET_NSD_DET_CTR_REG       );
	DVCDAT( NsdSymbDetCtr       )     = SYSRDSLVREG   ( pd, NSDET_NSD_SYMB_DET_CTR_REG  );
	DVCDAT( NsdChnDetMap        )     = SYSRDSLVREG   ( pd, NSDET_NSD_CHN_DET_REG       );

	for( int i = 0; i < NSDET_N_CHAN; i++ ) {
		SelChanData( pd, i );
		DVCDAT( NsdChnFrqIdx  [ i ] ) = SYSRDSLVREG64 ( pd, NSDET_NSD_CHN_FRQ_IDX_LO_REG );
		DVCDAT( NsdChnSumNoise[ i ] ) = SYSRDSLVREG   ( pd, NSDET_NSD_SUM_NOISE_REG      );
		DVCDAT( NsdChnSumPeaks[ i ] ) = SYSRDSLVREG   ( pd, NSDET_NSD_SUM_PEAKS_REG      );
		DVCDAT( NsdChnDetHist [ i ] ) = SYSRDSLVREG   ( pd, NSDET_NSD_CHN_DET_HIST_REG   );
		DVCDAT( NsdChnDb0     [ i ] ) = SYSRDSLVREG   ( pd, NSDET_NSD_DEBUG_0_REG        );
		DVCDAT( NsdChnDb1     [ i ] ) = SYSRDSLVREG   ( pd, NSDET_NSD_DEBUG_1_REG        );
		DVCDAT( NsdChnDb2     [ i ] ) = SYSRDSLVREG   ( pd, NSDET_NSD_DEBUG_2_REG        );
		DVCDAT( NsdChnDb3     [ i ] ) = SYSRDSLVREG   ( pd, NSDET_NSD_DEBUG_3_REG        );
        DVCDAT( NsdChnDb4     [ i ] ) = SYSRDSLVREG   ( pd, NSDET_NSD_DEBUG_4_REG        );
        DVCDAT( NsdChnDb5     [ i ] ) = SYSRDSLVREG   ( pd, NSDET_NSD_DEBUG_5_REG        );
        DVCDAT( NsdChnDb6     [ i ] ) = SYSRDSLVREG   ( pd, NSDET_NSD_DEBUG_6_REG        );
        DVCDAT( NsdChnDb7     [ i ] ) = SYSRDSLVREG   ( pd, NSDET_NSD_DEBUG_7_REG        );
        DVCDAT( NsdChnDb8     [ i ] ) = SYSRDSLVREG   ( pd, NSDET_NSD_DEBUG_8_REG        );
        DVCDAT( NsdChnDb9     [ i ] ) = SYSRDSLVREG   ( pd, NSDET_NSD_DEBUG_9_REG        );
        DVCDAT( NsdChnDb10    [ i ] ) = SYSRDSLVREG   ( pd, NSDET_NSD_DEBUG_10_REG       );
        DVCDAT( NsdChnDb11    [ i ] ) = SYSRDSLVREG   ( pd, NSDET_NSD_DEBUG_11_REG       );
        DVCDAT( NsdSdState    [ i ] ) = SYSRDSLVREG   ( pd, NSDET_NSD_SD_STATE_REG       );
        DVCDAT( NsdTdState    [ i ] ) = SYSRDSLVREG   ( pd, NSDET_NSD_TD_STATE_REG       );
	}

    DVCDAT( NsdSAxisAbsTvalid   )     = SYSRDSLVREG   ( pd, NSDET_S_AXIS_ABS_TVALID_REG  );
    DVCDAT( NsdSAxisAbsTready   )     = SYSRDSLVREG   ( pd, NSDET_S_AXIS_ABS_TREADY_REG  );
    DVCDAT( NsdBufAbsTvalid     )     = SYSRDSLVREG   ( pd, NSDET_BUF_ABS_TVALID_REG     );
    DVCDAT( NsdBufAbsTready     )     = SYSRDSLVREG   ( pd, NSDET_BUF_ABS_TREADY_REG     );
    DVCDAT( NsdPsSampleMap0     )     = SYSRDSLVREG   ( pd, NSDET_PS_SAMPLE_MAP_0_REG    );
    DVCDAT( NsdPsSampleMap1     )     = SYSRDSLVREG   ( pd, NSDET_PS_SAMPLE_MAP_1_REG    );
    DVCDAT( NsdPsSampleMap2     )     = SYSRDSLVREG   ( pd, NSDET_PS_SAMPLE_MAP_2_REG    );
    DVCDAT( NsdPsSampleMap3     )     = SYSRDSLVREG   ( pd, NSDET_PS_SAMPLE_MAP_3_REG    );
    DVCDAT( NsdExTimer          )     = SYSRDSLVREG   ( pd, NSDET_EX_TIMER_REG           );

    DVCDAT( NsdSymbDtaIdx       )     = SYSRDSLVREG64 ( pd, NSDET_SYMB_DTA_IDX_REG       );
    DVCDAT( NsdChanDetIdx       )     = SYSRDSLVREG   ( pd, NSDET_CHAN_DET_IDX_REG       );
    DVCDAT( NsdSumNoiseIdx      )     = SYSRDSLVREG   ( pd, NSDET_SUM_NOISE_IDX_REG      );
    DVCDAT( NsdSumPeaksIdx      )     = SYSRDSLVREG   ( pd, NSDET_SUM_PEAKS_IDX_REG      );
    DVCDAT( NsdRxRxdCtr         )     = SYSRDSLVREG64 ( pd, NSDET_RX_RXD_CTR_REG         );
    DVCDAT( NsdRxFFTCtr         )     = SYSRDSLVREG64 ( pd, NSDET_RX_FFT_CTR_REG         );
    DVCDAT( NsdTxFFTCtr         )     = SYSRDSLVREG64 ( pd, NSDET_TX_FFT_CTR_REG         );
    DVCDAT( NsdTxPrsCtr         )     = SYSRDSLVREG   ( pd, NSDET_TX_PRS_CTR_REG         );
    DVCDAT( NsdTxDtaCtr         )     = SYSRDSLVREG   ( pd, NSDET_TX_DTA_CTR_REG         );
    DVCDAT( NsdTxSymCtr         )     = SYSRDSLVREG   ( pd, NSDET_TX_SYM_CTR_REG         );
    DVCDAT( NsdTxPrsStart       )     = SYSRDSLVREG64 ( pd, NSDET_TX_PRS_START_REG       );
    DVCDAT( NsdTxDtaStart       )     = SYSRDSLVREG64 ( pd, NSDET_TX_DTA_START_REG       );

	DVCDAT( NsdFFTCfgDone       )     = ( StatReg & NSDET_FFT_CFG_DONE_MASK        ? true : false );
	DVCDAT( NsdDetDone          )     = ( StatReg & NSDET_NSD_DET_DONE_MASK        ? true : false );
    DVCDAT( NsdTxdTxDone        )     = ( StatReg & NSDET_TXD_TX_DONE_MASK         ? true : false );
    DVCDAT( NsdSAxisRxdTvalid   )     = ( StatReg & NSDET_S_AXIS_RXD_TVALID_MASK   ? true : false );
    DVCDAT( NsdSAxisRxdTready   )     = ( StatReg & NSDET_S_AXIS_RXD_TREADY_MASK   ? true : false );
    DVCDAT( NsdMAxisFFTTvalid   )     = ( StatReg & NSDET_M_AXIS_FFT_TVALID_MASK   ? true : false );
    DVCDAT( NsdMAXisFFTTready   )     = ( StatReg & NSDET_M_AXIS_FFT_TREADY_MASK   ? true : false );
    DVCDAT( NsdMAxisPrsTvalid   )     = ( StatReg & NSDET_M_AXIS_PRS_TVALID_MASK   ? true : false );
    DVCDAT( NsdMAXisPrsTready   )     = ( StatReg & NSDET_M_AXIS_PRS_TREADY_MASK   ? true : false );
    DVCDAT( NsdMAxisDtaTvalid   )     = ( StatReg & NSDET_M_AXIS_DTA_TVALID_MASK   ? true : false );
    DVCDAT( NsdMAXisDtaTready   )     = ( StatReg & NSDET_M_AXIS_DTA_TREADY_MASK   ? true : false );
    DVCDAT( NsdBufAbsTvalidF    )     = ( StatReg & NSDET_BUF_ABS_TVALID_F_MASK    ? true : false );
    DVCDAT( NsdExState          )     = ( StatReg & NSDET_EX_STATE_MASK  ) >> NSDET_EX_STATE_SHR;
    DVCDAT( NsdNsdState         )     = ( StatReg & NSDET_NSD_STATE_MASK ) >> NSDET_NSD_STATE_SHR;
    DVCDAT( NsdTxdState         )     = ( StatReg & NSDET_TX_STATE_MASK  ) >> NSDET_TX_STATE_SHR;
    DVCDAT( NsdSysResetn        )     = ( StatReg & NSDET_SYS_RESETN_MASK   ? false : true  );
    DVCDAT( NsdFFTResetn        )     = ( StatReg & NSDET_FFT_RESETN_MASK   ? false : true  );
    DVCDAT( NsdFifoResetn       )     = ( StatReg & NSDET_FIFO_RESETN_MASK  ? false : true  );
    DVCDAT( NsdFFTCfgRqst       )     = ( StatReg & NSDET_FFT_CFG_RQST_MASK ? true  : false );
    SelChanData( pd, 0 );
}

static const char s_EX_UNK0             [ ] = "UNKNOWN";             // 0
static const char s_EX_UNK1             [ ] = "UNKNOWN";             // 1
static const char s_EX_UNK2             [ ] = "UNKNOWN";             // 2
static const char s_EX_UNK3             [ ] = "UNKNOWN";             // 3
static const char s_EX_UNK4             [ ] = "UNKNOWN";             // 4
static const char s_EX_UNK5             [ ] = "UNKNOWN";             // 5
static const char s_EX_UNK6             [ ] = "UNKNOWN";             // 6
static const char s_EX_UNK7             [ ] = "UNKNOWN";             // 7
static const char s_EX_TX_DATA          [ ] = "EX_TX_DATA";          // 8
static const char s_EX_IDLE             [ ] = "EX_IDLE";             // 9
static const char s_EX_UNKA             [ ] = "UNKNOWN";             // A
static const char s_EX_WAIT_NSD         [ ] = "EX_WAIT_NSD";         // B
static const char s_EX_TX_PRS           [ ] = "EX_TX_PRS";           // C
static const char s_EX_WAIT_FFT_CFG_ACK [ ] = "EX_WAIT_FFT_CFG_ACK"; // D
static const char s_EX_WAIT_FFT_RST     [ ] = "EX_WAIT_FFT_RST";     // E
static const char s_EX_WAIT_FFT_CFG     [ ] = "EX_WAIT_FFT_CFG";     // F

static const char *ExStateStr [ 16 ] = {
        s_EX_UNK0,
        s_EX_UNK1,
        s_EX_UNK2,
        s_EX_UNK3,
        s_EX_UNK4,
        s_EX_UNK5,
        s_EX_UNK6,
        s_EX_UNK7,
        s_EX_TX_DATA,
        s_EX_IDLE,
        s_EX_UNKA,
        s_EX_WAIT_NSD,
        s_EX_TX_PRS,
        s_EX_WAIT_FFT_CFG_ACK,
        s_EX_WAIT_FFT_RST,
        s_EX_WAIT_FFT_CFG
};

static const char s_NSD_UNK0           [ ] = "UNKNOWN";            // 0
static const char s_NSD_UNK1           [ ] = "UNKNOWN";            // 1
static const char s_NSD_UNK2           [ ] = "UNKNOWN";            // 2
static const char s_NSD_UNK3           [ ] = "UNKNOWN";            // 3
static const char s_NSD_UNK4           [ ] = "UNKNOWN";            // 4
static const char s_NSD_UNK5           [ ] = "UNKNOWN";            // 5
static const char s_NSD_UNK6           [ ] = "UNKNOWN";            // 6
static const char s_NSD_UNK7           [ ] = "UNKNOWN";            // 7
static const char s_NSD_UNK8           [ ] = "UNKNOWN";            // 8
static const char s_NSD_IDLE           [ ] = "NSD_IDLE";           // 9
static const char s_NSD_UNKA           [ ] = "UNKNOWN";            // A
static const char s_NSD_WAIT_PASS_2    [ ] = "NSD_WAIT_PASS_2";    // B
static const char s_NSD_UNKC           [ ] = "UNKNOWN";            // C
static const char s_NSD_COMP_DATA      [ ] = "NSD_CALC_FRAME_IDX"; // D
static const char s_NSD_UNKE           [ ] = "UNKNOWN";            // E
static const char s_NSD_CALC_FRAME_IDX [ ] = "NSD_COMP_DATA";      // F

static const char *NsdStateStr [ 16 ] = {
        s_NSD_UNK0,
        s_NSD_UNK1,
        s_NSD_UNK2,
        s_NSD_UNK3,
        s_NSD_UNK4,
        s_NSD_UNK5,
        s_NSD_UNK6,
        s_NSD_UNK7,
        s_NSD_UNK8,
        s_NSD_IDLE,
        s_NSD_UNKA,
        s_NSD_WAIT_PASS_2,
        s_NSD_UNKC,
        s_NSD_COMP_DATA,
        s_NSD_UNKE,
        s_NSD_CALC_FRAME_IDX
};

static const char s_NSD_TX_UNK0           [ ] = "UNKNOWN";            // 0
static const char s_NSD_TX_UNK1           [ ] = "UNKNOWN";            // 1
static const char s_NSD_TX_UNK2           [ ] = "UNKNOWN";            // 2
static const char s_NSD_TX_UNK3           [ ] = "UNKNOWN";            // 3
static const char s_NSD_TX_UNK4           [ ] = "UNKNOWN";            // 4
static const char s_NSD_TX_UNK5           [ ] = "UNKNOWN";            // 5
static const char s_NSD_TX_UNK6           [ ] = "UNKNOWN";            // 6
static const char s_NSD_TX_UNK7           [ ] = "UNKNOWN";            // 7
static const char s_NSD_TX_UNK8           [ ] = "UNKNOWN";            // 8
static const char s_NSD_TX_IDLE           [ ] = "TX_IDLE";            // 9
static const char s_NSD_TX_FFT            [ ] = "TX_FFT";             // A
static const char s_NSD_TX_UNKB           [ ] = "UNKNOWN";            // B
static const char s_NSD_TX_WAIT_PRS       [ ] = "TX_WAIT_PRS";        // C
static const char s_NSD_TX_PRS            [ ] = "TX_PRS";             // D
static const char s_NSD_TX_WAIT_DTA       [ ] = "TX_WAIT_DTA";        // E
static const char s_NSD_TX_DTA            [ ] = "TX_DTA";             // F

static const char *TxStateStr [ 16 ] = {
        s_NSD_TX_UNK0,
        s_NSD_TX_UNK1,
        s_NSD_TX_UNK2,
        s_NSD_TX_UNK3,
        s_NSD_TX_UNK4,
        s_NSD_TX_UNK5,
        s_NSD_TX_UNK6,
        s_NSD_TX_UNK7,
        s_NSD_TX_UNK8,
        s_NSD_TX_IDLE,
        s_NSD_TX_FFT,
        s_NSD_TX_UNKB,
        s_NSD_TX_WAIT_PRS,
        s_NSD_TX_PRS,
        s_NSD_TX_WAIT_DTA,
        s_NSD_TX_DTA
};

static const char s_NSD_TD_UNK0            [ ] = "UNKNOWN";            // 0
static const char s_NSD_TD_UNK1            [ ] = "UNKNOWN";            // 1
static const char s_NSD_TD_UNK2            [ ] = "UNKNOWN";            // 2
static const char s_NSD_TD_UNK3            [ ] = "UNKNOWN";            // 3
static const char s_NSD_TD_UNK4            [ ] = "UNKNOWN";            // 4
static const char s_NSD_TD_UNK5            [ ] = "UNKNOWN";            // 5
static const char s_NSD_TD_UNK6            [ ] = "UNKNOWN";            // 6
static const char s_NSD_TD_UNK7            [ ] = "UNKNOWN";            // 7
static const char s_NSD_TD_UNK8            [ ] = "UNKNOWN";            // 8
static const char s_NSD_TD_IDLE            [ ] = "TD_IDLE";            // 9
static const char s_NSD_TD_UNKA            [ ] = "UNKNOWN";            // A
static const char s_NSD_TD_CALC_SCALED_MAX [ ] = "TD_CALC_SCALED_MAX"; // B
static const char s_NSD_TD_UNKC            [ ] = "UNKNOWN";            // C
static const char s_NSD_TD_UNKD            [ ] = "UNKNOWN";            // D
static const char s_NSD_TD_UNKE            [ ] = "UNKNOWN";            // E
static const char s_NSD_TD_TD_OUTPUT       [ ] = "TD_OUTPUT";          // F

static const char *NsdTdStateStr [ 16 ] = {
        s_NSD_TD_UNK0,
        s_NSD_TD_UNK1,
        s_NSD_TD_UNK2,
        s_NSD_TD_UNK3,
        s_NSD_TD_UNK4,
        s_NSD_TD_UNK5,
        s_NSD_TD_UNK6,
        s_NSD_TD_UNK7,
        s_NSD_TD_UNK8,
        s_NSD_TD_IDLE,
        s_NSD_TD_UNKA,
        s_NSD_TD_CALC_SCALED_MAX,
        s_NSD_TD_UNKC,
        s_NSD_TD_UNKD,
        s_NSD_TD_UNKE,
        s_NSD_TD_TD_OUTPUT
};

static const char s_NSD_SD_UNK0      [ ] = "UNKNOWN";            // 0
static const char s_NSD_SD_UNK1      [ ] = "UNKNOWN";            // 1
static const char s_NSD_SD_UNK2      [ ] = "UNKNOWN";            // 2
static const char s_NSD_SD_UNK3      [ ] = "UNKNOWN";            // 3
static const char s_NSD_SD_UNK4      [ ] = "UNKNOWN";            // 4
static const char s_NSD_SD_UNK5      [ ] = "UNKNOWN";            // 5
static const char s_NSD_SD_UNK6      [ ] = "UNKNOWN";            // 6
static const char s_NSD_SD_UNK7      [ ] = "UNKNOWN";            // 7
static const char s_NSD_SD_UNK8      [ ] = "UNKNOWN";            // 8
static const char s_NSD_SD_IDLE      [ ] = "SD_IDLE";            // 9
static const char s_NSD_SD_UNKA      [ ] = "UNKNOWN";            // A
static const char s_NSD_SD_BLOCK_END [ ] = "SD_BLOCK_END";       // B
static const char s_NSD_SD_UNKC      [ ] = "UNKNOWN";            // C
static const char s_NSD_SD_UNKD      [ ] = "UNKNOWN";            // D
static const char s_NSD_SD_UNKE      [ ] = "UNKNOWN";            // E
static const char s_NSD_SD_UNKF      [ ] = "UNKNOWN";            // F

static const char *NsdSdStateStr [ 16 ] = {
        s_NSD_SD_UNK0,
        s_NSD_SD_UNK1,
        s_NSD_SD_UNK2,
        s_NSD_SD_UNK3,
        s_NSD_SD_UNK4,
        s_NSD_SD_UNK5,
        s_NSD_SD_UNK6,
        s_NSD_SD_UNK7,
        s_NSD_SD_UNK8,
        s_NSD_SD_IDLE,
        s_NSD_SD_UNKA,
        s_NSD_SD_BLOCK_END,
        s_NSD_SD_UNKC,
        s_NSD_SD_UNKD,
        s_NSD_SD_UNKE,
        s_NSD_SD_UNKF
};

static void DvcDataDump( NSDet *pd )
{
    char sbin [ 2 * NSDET_N_TONE + 1 ];
    char s1   [ NSDET_N_TONE + 1 ];
    char s2   [ NSDET_N_TONE + 1 ];
	DvcDataUpd( pd );
	cprintf( "\nCurrent block index         : %lu", toULint( DVCDAT( NsdBlkIdx           ) ) );
	cprintf( "\nData Blocks Received        : %u",  toUint ( DVCDAT( NsdBlkRxCtr         ) ) );
	cprintf( "\n" );
    cprintf( "\nFFT CFG done                : %s", TRUEFALSE( DVCDAT( NsdFFTCfgDone      ) ) );
    cprintf( "\nDet done                    : %s", TRUEFALSE( DVCDAT( NsdDetDone         ) ) );
    cprintf( "\nTxd tx done                 : %s", TRUEFALSE( DVCDAT( NsdTxdTxDone       ) ) );
    cprintf( "\ns_axis_rxd_tvalid           : %s", OFFON( DVCDAT( NsdSAxisRxdTvalid      ) ) );
    cprintf( "\ns_axis_rxd_tready           : %s", OFFON( DVCDAT( NsdSAxisRxdTready      ) ) );
    cprintf( "\nm_axis_fft_tvalid           : %s", OFFON( DVCDAT( NsdMAxisFFTTvalid      ) ) );
    cprintf( "\nm_axis_fft_tready           : %s", OFFON( DVCDAT( NsdMAXisFFTTready      ) ) );
    cprintf( "\nm_axis_prs_tvalid           : %s", OFFON( DVCDAT( NsdMAxisPrsTvalid      ) ) );
    cprintf( "\nm_axis_prs_tready           : %s", OFFON( DVCDAT( NsdMAXisPrsTready      ) ) );
    cprintf( "\nm_axis_dta_tvalid           : %s", OFFON( DVCDAT( NsdMAxisDtaTvalid      ) ) );
    cprintf( "\nm_axis_dta_tready           : %s", OFFON( DVCDAT( NsdMAXisDtaTready      ) ) );
    cprintf( "\nbuf_abs_tvalid_f            : %s", OFFON( DVCDAT( NsdBufAbsTvalidF       ) ) );
    cprintf( "\nsys_resetn                  : %s", OFFON( DVCDAT( NsdSysResetn           ) ) );
    cprintf( "\nfft_resetn                  : %s", OFFON( DVCDAT( NsdFFTResetn           ) ) );
    cprintf( "\ndsp_resetn                  : %s", OFFON( DVCDAT( NsdFifoResetn          ) ) );
    cprintf( "\nfft_cfg_rqst                : %s", OFFON( DVCDAT( NsdFFTCfgRqst          ) ) );
    u32bin2str( sbin, DVCDAT( NsdPsSampleMap0   ), NSDET_N_TONE / 8, 0 );
    cprintf( "\nps_sample_map_0             : %s", sbin );
    u32bin2str( sbin, DVCDAT( NsdPsSampleMap1   ), NSDET_N_TONE / 8, 0 );
    cprintf( "\nps_sample_map_1             : %s", sbin );
    u32bin2str( sbin, DVCDAT( NsdPsSampleMap2   ), NSDET_N_TONE / 8, 0 );
    cprintf( "\nps_sample_map_2             : %s", sbin );
    u32bin2str( sbin, DVCDAT( NsdPsSampleMap3   ), NSDET_N_TONE / 8, 0 );
    cprintf( "\nps_sample_map_3             : %s", sbin );
    u32bin2str( sbin, DVCDAT( NsdSAxisAbsTvalid ), NSDET_N_TONE / 8, 0 );
    cprintf( "\ns_axis_abs_tvalid           : %s", sbin );
    u32bin2str( sbin, DVCDAT( NsdSAxisAbsTready ), NSDET_N_TONE / 8, 0 );
    cprintf( "\ns_axis_abs_tready           : %s", sbin );
    u32bin2str( sbin, DVCDAT( NsdBufAbsTvalid   ), NSDET_N_TONE / 8, 0 );
    cprintf( "\nbuf_abs_tvalid              : %s", sbin );
    u32bin2str( sbin, DVCDAT( NsdBufAbsTready   ), NSDET_N_TONE / 8, 0 );
    cprintf( "\nbuf_abs_tready              : %s", sbin );
//    cprintf( "\nex_timer                    : %u", toUint( DVCDAT( NsdExTimer ) ) );

    cprintf( "\nChan Data:" );
    cprintf( "\nnch sd_state           td_state        CDHist     SumNse     SumPks     FrqIdx           sbk    sbrn   sp0    sp1    sp2    sp3    i symbsumpk  symbsumns  symbblkct" );
    cprintf( "\n=== ========           ========        ======     ======     ======     ======           ===    ====   ===    ===    ===    ===    = =========  =========  =========" );
//              00: BULL_SHIT          BALLSACE        0123456789 0123456789 0123456789 4321432143214321 012345 012345 012345 012345 012345 012345 0 0123456789 0123456789 0123456789
    for( int i = 0; i < NSDET_N_CHAN; i++ ) {
        u32toFidx( s1, DVCDAT( NsdChnFrqIdx [ i ] ) );
        u32toFidx( s2, DVCDAT( NsdChnDb4 [ i ] ) );
        cprintf("\n%02u: %-16s%-19s%-11u%-11u%-11u%-17s%-7u%-7u%-7u%-7u%-7u%-7u%-2u%-11u%-11u%-11u",
                i,
                NsdSdStateStr [ DVCDAT( NsdSdState[ i ] ) ],
                NsdTdStateStr [ DVCDAT( NsdTdState[ i ] ) ],
                toUint( DVCDAT( NsdChnDetHist  [ i ] ) ),
                toUint( DVCDAT( NsdChnSumNoise [ i ] ) ),
                toUint( DVCDAT( NsdChnSumPeaks [ i ] ) ),
                s1,
                toUint( DVCDAT( NsdChnDb0      [ i ] ) ),
                toUint( DVCDAT( NsdChnDb1      [ i ] ) ),
                toUint( DVCDAT( NsdChnDb2      [ i ] ) ),
                toUint( DVCDAT( NsdChnDb3      [ i ] ) ),
                toUint( DVCDAT( NsdChnDb4      [ i ] ) ),
                toUint( DVCDAT( NsdChnDb5      [ i ] ) ),
                toUint( DVCDAT( NsdChnDb6      [ i ] ) ),
                toUint( DVCDAT( NsdChnDb9      [ i ] ) ),
                toUint( DVCDAT( NsdChnDb8      [ i ] ) ),
                toUint( DVCDAT( NsdChnDb11     [ i ] ) )
                );
    }
    cprintf( "\nex_state             : %s", ExStateStr  [ DVCDAT( NsdExState  ) ] );
    cprintf( "\nnsd_state            : %s", NsdStateStr [ DVCDAT( NsdNsdState ) ] );
    cprintf( "\ntxd_state            : %s", TxStateStr  [ DVCDAT( NsdTxdState ) ] );
    cprintf( "\n" );
    cprintf( "\nSum Noise Min (last block ) : %u",  toUint ( DVCDAT( NsdSumNoiseMin      ) ) );
    cprintf( "\nSum Peaks Max (last block ) : %u",  toUint ( DVCDAT( NsdSumPeaksMax      ) ) );
    cprintf( "\nChans detected              : %u",  toUint ( DVCDAT( NsdDetCtr           ) ) );
    u32bin2str( sbin, DVCDAT( NsdChnDetMap  ), NSDET_N_TONE / 8, 0 );
    cprintf( "\nLast channel detect map     : %s",  sbin );
    cprintf( "\nTotal Symbols detected      : %u",  toUint ( DVCDAT( NsdSymbDetCtr       ) ) );
    cprintf( "\nnsd_symb_dta_idx     : %lu", toULint ( DVCDAT( NsdSymbDtaIdx  ) ) );
    cprintf( "\nnsd_chan_det_idx     : %u",  toUint  ( DVCDAT( NsdChanDetIdx  ) ) );
    cprintf( "\nnsd_sum_noise_idx    : %u",  toUint  ( DVCDAT( NsdSumNoiseIdx ) ) );
    cprintf( "\nnsd_sum_peaks_idx    : %u",  toUint  ( DVCDAT( NsdSumPeaksIdx ) ) );
    cprintf( "\nrx_rxd_ctr           : %lu", toULint ( DVCDAT( NsdRxRxdCtr    ) ) );
    cprintf( "\nrx_fft_ctr           : %lu", toULint ( DVCDAT( NsdRxFFTCtr    ) ) );
    cprintf( "\ntx_fft_ctr           : %lu", toULint ( DVCDAT( NsdTxFFTCtr    ) ) );
    cprintf( "\ntx_prs_ctr           : %u",  toUint  ( DVCDAT( NsdTxPrsCtr    ) ) );
    cprintf( "\ntx_dta_ctr           : %u",  toUint  ( DVCDAT( NsdTxDtaCtr    ) ) );
    cprintf( "\ntx_sym_ctr           : %u",  toUint  ( DVCDAT( NsdTxSymCtr    ) ) );
    cprintf( "\ntx_prs_start         : %lu", toULint ( DVCDAT( NsdTxPrsStart  ) ) );
    cprintf( "\ntx_dta_start         : %lu", toULint ( DVCDAT( NsdTxDtaStart  ) ) );
}

void FrqIdxStr( NSDet *pd, char *sidx, char * sfrq, u64 idxref )
{
	int j = 2 * NSDET_N_BLK_DET - 1;
	for (int i = 0; i < NSDET_N_BLK_DET; i++)
	{
		if (j < 2 * NSDET_N_BLK_DET - 1)
		{
			sidx[ j ] = ' ';
			sfrq[ j ] = ' ';
		}
		else
		{
			sidx[ j ] = '\0';
			sfrq[ j ] = '\0';
		}
		j--;
		sidx[ j ] = ( idxref & 3 ) + '0';
		sfrq[ j ] = DVCPAR( SampleMapFrq[ idxref & 3 ] ) + '0';
		idxref >>= 2;
		j--;
	}
}

static void DvcDescDump( NSDet *pd )
{
    char sidx[ 2 * NSDET_N_BLK_DET ];
    char sfrq[ 2 * NSDET_N_BLK_DET ];
	char sbin[ 2 * NSDET_N_TONE + 1 ];
    u32  FrqIdxRef = DVCPAR( FrqIdxRef );

	cprintf("\nRNoise           : %u",   toUint( DVCPAR( RNoise ))          );
	u32bin2str( sbin, DVCPAR( FrqIdxRef ), NSDET_N_BLK_DET / 4, 2 );
	cprintf("\nFrqIdxRef        : 0b%s", sbin );
    FrqIdxStr( pd, sidx, sfrq, FrqIdxRef );
	cprintf("\nRef idx          : %s",   sidx );
	cprintf("\nRef frq          : %s",   sfrq );
	cprintf("\nRef frq (MHz)    : ");
	for ( int i = 0; i < NSDET_N_BLK_DET; i++ ) {
		cprintf( "%-3.6g ", ( double ) 1.0e-6 * DVCPAR( NsdFrqBlk[ i ] ) * TxRxFreqs [ SYSPAR( TxRxFreq ) ] / NSDET_N_TONE );
	}
	cprintf("\n");
	cprintf("\nChan data sel    : %u",   toUint( DVCPAR( ChanDataSel ))     );
	cprintf("\nNSD test mode    : %u",   toUint( DVCPAR( NsdTstMode ))      );
	u32bin2str( sbin, DVCPAR( SampleMap[ 0 ] ), NSDET_N_TONE / 8, 0 );
	cprintf("\nSampleMap[ 0 ]   : %s",   sbin );
	u32bin2str( sbin, DVCPAR( SampleMap[ 1 ] ), NSDET_N_TONE / 8, 0 );
	cprintf("\nSampleMap[ 1 ]   : %s",   sbin );
	u32bin2str( sbin, DVCPAR( SampleMap[ 2 ] ), NSDET_N_TONE / 8, 0 );
	cprintf("\nSampleMap[ 2 ]   : %s",   sbin );
	u32bin2str( sbin, DVCPAR( SampleMap[ 3 ] ), NSDET_N_TONE / 8, 0 );
	cprintf("\nSampleMap[ 3 ]   : %s",   sbin );
	cprintf("\n");
}

static void SysDvceDump( NSDet *pd )
{
	SYSDVCEDUMP( pd );
}

static void SysDataDump( NSDet *pd )
{
	SYSDATADUMP( pd );
}

static void SysDescDump( NSDet *pd )
{
	cprintf("\nnFFT             : %u",   toUint( SYSPAR( nFFT        )) );
	cprintf("\nnCPre            : %u",   toUint( SYSPAR( nCPre       )) );
	cprintf("\nnSymbFrame       : %u",   toUint( SYSPAR( nSymbAGC ) + SYSPAR( nSymbNSD ) + SYSPAR ( nSymbPRS ) + SYSPAR( nSymbData )) );
    cprintf("\nnDataSymbFrame   : %u",   toUint( SYSPAR( nSymbData   )) );
	cprintf("\nnFrameBlock      : %u",   toUint( SYSPAR( nFrameBlock )) );
	cprintf("\nnSymbAGC         : %u",   toUint( SYSPAR( nSymbAGC    )) );
	cprintf("\nnSymbNSD         : %u",   toUint( SYSPAR( nSymbNSD    )) );
	cprintf("\n");
}


// Insert device-specific code here
//======================================
// Device-specific non member functions
//======================================
// Example: static inline void EnterAdmin    ( NSDet *pd ) { hwWriteCtrlReg( SYSDVC( BaseAddr ), NSDET_CTRL_REG << 2, NSDET_ENTER_ADMIN ); }

//=================================
// Device-specific member functions
//=================================
static const char *ExState( u32 StatReg )
{
	return( ExStateStr [ ( StatReg & NSDET_EX_STATE_MASK ) >> NSDET_EX_STATE_SHR ] );
}

static const char *NsdState( u32 StatReg )
{
	return( NsdStateStr [ ( StatReg & NSDET_NSD_STATE_MASK ) >> NSDET_NSD_STATE_SHR ] );
}

static const char *TxState( u32 StatReg )
{
	return ( TxStateStr [ ( StatReg & NSDET_EX_STATE_MASK ) >> NSDET_EX_STATE_SHR  ] );
}

static void Dstate( NSDet *pd )
{
	cprintf( "%s States:\n", SYSDVC( Name ) );
	cprintf( "[ex_state ]  : %s\n", ExState  ( SYSRDSLVREG( pd, NSDET_STAT_REG )) );
	cprintf( "[nsd_state]  : %s\n", NsdState ( SYSRDSLVREG( pd, NSDET_STAT_REG )) );
	cprintf( "[tx_state]   : %s\n", TxState  ( SYSRDSLVREG( pd, NSDET_STAT_REG )) );
    cprintf( "\n" );
}

static void MkFrqRef( NSDet *pd )
{
	DVCPAR( FrqIdxRef ) = 0;

	_Bool SampleMapSel[ NSDET_MAX_FREQ ];
	u32   FrqIdxRef = 0;

	for ( int i = 0; i < NSDET_MAX_FREQ; i ++ ) {
		SampleMapSel[ i ] = false;
		DVCPAR( SampleMapIdx[ i ] ) = 0xFFFFFFFF;
	}
	for ( int i = 0; i < NSDET_N_BLK_DET; i++ ) {
		SampleMapSel[ DVCPAR( NsdFrqBlk[ i ] ) - 1 ] = true;
	}
	for( int j = 0; j < NSDET_N_FREQ; j++ ) {
		DVCPAR( SampleMap[ j ]    ) = 0;
		DVCPAR( SampleMapFrq[ j ] ) = 0;
	}
	for ( int i = 0, j = 0; i < NSDET_MAX_FREQ; i++ ) {
		if ( SampleMapSel[ i ] )
		{
			DVCPAR( SampleMap[ j ] )    = DVCPAR( SampleMapTbl[ i ] );
			DVCPAR( SampleMapFrq[ j ] ) = i + 1;
			DVCPAR( SampleMapIdx[ i ] ) = j++;
		}
	}
	for ( int i = 0; i < NSDET_N_BLK_DET; i++ ) {
		FrqIdxRef = FrqIdxRef << 2 | DVCPAR( SampleMapIdx[ DVCPAR( NsdFrqBlk[ i ] )) - 1 ];
	}
    DVCPAR( FrqIdxRef ) = FrqIdxRef;


	SYSWRSLVREG( pd, NSDET_FRQ_IDX_REF_REG,   DVCPAR( FrqIdxRef ));
	SYSWRSLVREG( pd, NSDET_SAMPLE_MAP_0_REG,  DVCPAR( SampleMap[ 0 ] ));
	SYSWRSLVREG( pd, NSDET_SAMPLE_MAP_1_REG,  DVCPAR( SampleMap[ 1 ] ));
	SYSWRSLVREG( pd, NSDET_SAMPLE_MAP_2_REG,  DVCPAR( SampleMap[ 2 ] ));
	SYSWRSLVREG( pd, NSDET_SAMPLE_MAP_3_REG,  DVCPAR( SampleMap[ 3 ] ));
}

static void MkTSBuf( NSDet *pd )
{

	for ( int ntsymb = 0; ntsymb < SYS_N_TSYMB_MAX; ntsymb++ ) {
		for( int i = 0; i < SYS_MAX_CARRIERS; i++ ) {
			cplx s =  ccosh( I * 2 * M_PI * ( double ) i * DVCPAR( TsFreq[ ntsymb ] ) / NSDET_N_TONE );
			s32  re = ( s32 ) ( creal( s ) * SYS_RF_DAC_MAX );
			s32  im = ( s32 ) ( cimag( s ) * SYS_RF_DAC_MAX );
			SYSPTR( pTsBuf[ ntsymb * SYS_MAX_CARRIERS + i ] ) = ( *( u32 *) &re & 0x0000ffff ) | ( *( u32 *) &im << 16 );
		}
	}

    for( int nsymb = 0; nsymb < SYSPAR( nSymbNSD ); nsymb++ ) {
	    int sgn = 1.0;
	    u32 SymbOfs = ( SYSPAR( nSymbAGC ) + nsymb ) * SYS_MAX_CARRIERS;
	    for( int nblk = 0; nblk < NSDET_N_BLK_DET; nblk++ ) {
	        for ( int j = 0; j < NSDET_N_TONE; j ++ ) {
	            cplx s =  sgn * ccosh( I * 2 * M_PI * ( double ) j * DVCPAR( NsdFrqBlk[ nblk ] ) / NSDET_N_TONE );
	            s32  re = ( s32 ) ( creal( s ) * SYS_RF_DAC_MAX / 1.0 );
	            s32  im = ( s32 ) ( cimag( s ) * SYS_RF_DAC_MAX / 1.0 );
	            SYSPTR( pTsBuf[ ( nblk * NSDET_N_TONE + SymbOfs + j ) ] ) = ( *( u32 *) &re & 0x0000ffff ) | ( *( u32 *) &im << 16 );
	        }
	        sgn = -sgn;
	    }
  }

}

