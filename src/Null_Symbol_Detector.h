/*
 * Null_Symbol_Detector.h
 *
 *  Created on: 20 Jun 2018
 *      Author: JULIAN MORTIMER
 */

#ifndef NULL_SYMBOL_DETECTOR_H_
#define NULL_SYMBOL_DETECTOR_H_

#define NSDPAR( Param ) DEVICE_LOC_PAR( pd, NSD, Param )
#define NSDDAT( Param ) DEVICE_LOC_DAT( pd, NSD, Param )


#define NSDET_N_SLAVE_REGS              ( 64 )
#define NSDET_N_CHAN                    ( 16 )
#define NSDET_N_TONE                    ( 16 )
#define NSDET_N_BLK_DET                 ( 16 )
#define NSDET_N_FREQ                    ( 4 )
#define NSDET_N_FFT                     ( NSDET_N_TONE )
#define NSDET_MAX_FREQ                  ( 7 )
#define NSDET_MAX_BLKS_SYMB             ( SYS_MAX_CARRIERS / NSDET_N_TONE )
#define NSDET_FRQ_IDX_WIDTH             ( 2 )
#define NSDET_N_FRQ_REF                 (( NSDET_N_BLK * NSDET_FRQ_IDX_WIDTH + 32 - 1 ) / 32 )
#define NSDET_DET_THR                   ( 0 )
#define NSDET_N_BITS_SCALE_SUM_BLOCK    ( 12 )
#define NSDET_CTRL_REG                  ( 0 )

#define NSDET_GPO_REG                   ( 2 )
#define NSDET_N_FFT_REG                 ( 3 )
#define NSDET_N_CPRE_REG                ( 4 )
#define NSDET_N_SYMB_AGC_REG            ( 5 )
#define NSDET_N_SYMB_NSD_REG            ( 6 )
#define NSDET_N_SYMB_PRS_REG            ( 7 )
#define NSDET_N_SYMB_DTA_REG            ( 8 )
#define NSDET_N_FRAME_BLOCK_REG         ( 9 )
#define NSDET_R_NOISE_REG               ( 10 )
#define NSDET_FRQ_IDX_REF_REG           ( 11 )
#define NSDET_DET_THR_REG               ( 12 )
#define NSDET_CHAN_DATA_SEL_REG         ( 13 )
#define NSDET_TST_MODE_REG              ( 14 )
#define NSDET_SAMPLE_MAP_0_REG          ( 15 )
#define NSDET_SAMPLE_MAP_1_REG          ( 16 )
#define NSDET_SAMPLE_MAP_2_REG          ( 17 )
#define NSDET_SAMPLE_MAP_3_REG          ( 18 )

#define NSDET_SYS_RESET_B               ( 0x00000001 )
#define NSDET_FIFO_RESET_B              ( 0x00000002 )
#define NSDET_FFT_RESET_B               ( 0x00000004 )
#define NSDET_FFT_CONFIG_B              ( 0x00000008 )
#define NSDET_S_AXIS_DATA_TREADY_G      ( 0x00000010 )


#define NSDET_INFO_REG                  ( 1 )  // slv_reg_01 <= info_reg;
#define NSDET_STAT_REG                 ( 32 )  // slv_reg_32 <= status_reg;


#define NSDET_NSD_BLK_IDX_LO_REG       ( 39 )  // slv_reg_39 <= to_slv_ul( nsd_blk_idx         , C_PSL );
#define NSDET_NSD_BLK_IDX_HI_REG       ( 40 )  // slv_reg_40 <= to_slv_uh( nsd_blk_idx         , C_PSL );

#define NSDET_NSD_CHN_FRQ_IDX_LO_REG   ( 41 )  // slv_reg_41  <= to_slv_vl( nsd_chan_frq_idx_p  ( PSCDS ), C_PSL );
#define NSDET_NSD_CHN_FRQ_IDX_HI_REG   ( 42 )  // slv_reg_42  <= ( others=> '0' );
#define NSDET_NSD_SUM_NOISE_REG        ( 43 )  // slv_reg_43  <= to_slv_u( nsd_chan_sum_noise_p ( PSCDS ), C_PSL );
#define NSDET_NSD_SUM_PEAKS_REG        ( 44 )  // slv_reg_44  <= to_slv_u( nsd_chan_sum_peaks_p ( PSCDS ), C_PSL );
#define NSDET_NSD_CHN_DET_HIST_REG     ( 45 )  // slv_reg_45  <= to_slv_u( nsd_chan_det_hist_p  ( PSCDS ), C_PSL );
#define NSDET_NSD_DEBUG_0_REG          ( 46 )  // slv_reg_46  <= to_slv_v( ps_in_db_0           ( PSCDS ), C_PSL );
#define NSDET_NSD_DEBUG_1_REG          ( 47 )  // slv_reg_47  <= to_slv_v( ps_in_db_1           ( PSCDS ), C_PSL );
#define NSDET_NSD_DEBUG_2_REG          ( 48 )  // slv_reg_48  <= to_slv_v( ps_in_db_2           ( PSCDS ), C_PSL );
#define NSDET_NSD_DEBUG_3_REG          ( 49 )  // slv_reg_49  <= to_slv_v( ps_in_db_3           ( PSCDS ), C_PSL );
#define NSDET_NSD_DEBUG_4_REG          ( 75 )  // slv_reg_75  <= to_slv_v( ps_in_db_4           ( PSCDS ), C_PSL );
#define NSDET_NSD_DEBUG_5_REG          ( 76 )  // slv_reg_76  <= to_slv_v( ps_in_db_5           ( PSCDS ), C_PSL );
#define NSDET_NSD_DEBUG_6_REG          ( 77 )  // slv_reg_77  <= to_slv_v( ps_in_db_6           ( PSCDS ), C_PSL );
#define NSDET_NSD_DEBUG_7_REG          ( 78 )  // slv_reg_78  <= to_slv_v( ps_in_db_7           ( PSCDS ), C_PSL );
#define NSDET_NSD_DEBUG_8_REG          ( 79 )  // slv_reg_79  <= to_slv_v( ps_in_db_8           ( PSCDS ), C_PSL );
#define NSDET_NSD_DEBUG_9_REG          ( 80 )  // slv_reg_80  <= to_slv_v( ps_in_db_9           ( PSCDS ), C_PSL );
#define NSDET_NSD_DEBUG_10_REG         ( 81 )  // slv_reg_81  <= to_slv_v( ps_in_db_10          ( PSCDS ), C_PSL );
#define NSDET_NSD_DEBUG_11_REG         ( 82 )  // slv_reg_82  <= to_slv_v( ps_in_db_11          ( PSCDS ), C_PSL );
#define NSDET_NSD_SD_STATE_REG         ( 73 )  // slv_reg_73  <= to_slv_v( nsd_sd_state         ( PSCDS ), C_PSL );
#define NSDET_NSD_TD_STATE_REG         ( 74 )  // slv_reg_74  <= to_slv_v( nsd_td_state         ( PSCDS ), C_PSL );


#define NSDET_NSD_BLK_RX_CTR_REG       ( 50 )  // slv_reg_50 <= to_slv_u( nsd_blk_rx_count      , C_PSL );
#define NSDET_NSD_SUM_PEAKS_MAX_REG    ( 51 )  // slv_reg_51 <= to_slv_u( nsd_sum_peaks_max     , C_PSL );
#define NSDET_NSD_SUM_NOISE_MIN_REG    ( 52 )  // slv_reg_52 <= to_slv_u( nsd_sum_noise_min     , C_PSL );
#define NSDET_NSD_DET_CTR_REG          ( 53 )  // slv_reg_53 <= to_slv_u( nsd_det_ctr           , C_PSL );
#define NSDET_NSD_SYMB_DET_CTR_REG     ( 54 )  // slv_reg_54 <= to_slv_u( nsd_symb_det_count    , C_PSL );
#define NSDET_NSD_CHN_DET_REG          ( 55 )  // slv_reg_55 <= to_slv_v( nsd_chan_det_p        , C_PSL );
#define NSDET_RX_DATA_CTR_REG          ( 56 )  // slv_reg_56 <= to_slv_u( rx_data_ctr_p         , C_PSL );
#define NSDET_RX_SYMB_CTR_REG          ( 57 )  // slv_reg_57 <= to_slv_u( rx_symbol_ctr_p       , C_PSL );
#define NSDET_RX_FRAME_CTR_REG         ( 58 )  // slv_reg_58 <= to_slv_u( rx_frame_ctr_p        , C_PSL );
#define NSDET_RX_BLOCK_CTR_REG         ( 59 )  // slv_reg_59 <= to_slv_u( rx_block_ctr_p        , C_PSL );
#define NSDET_TX_DATA_CTR_REG          ( 60 )  // slv_reg_60 <= to_slv_u( tx_data_ctr_p         , C_PSL );
#define NSDET_TX_SYMB_CTR_REG          ( 61 )  // slv_reg_61 <= to_slv_u( tx_symbol_ctr_p       , C_PSL );
#define NSDET_TX_FRAME_CTR_REG         ( 62 )  // slv_reg_62 <= to_slv_u( tx_frame_ctr_p        , C_PSL );
#define NSDET_TX_BLOCK_CTR_REG         ( 63 )  // slv_reg_63 <= to_slv_u( tx_block_ctr_p        , C_PSL );


#define NSDET_S_AXIS_ABS_TVALID_REG    ( 64 )  // slv_reg_64  <= to_slv_v( s_axis_abs_tvalid    , C_PSL );
#define NSDET_S_AXIS_ABS_TREADY_REG    ( 65 )  // slv_reg_65  <= to_slv_v( s_axis_abs_tready    , C_PSL );
#define NSDET_BUF_ABS_TVALID_REG       ( 66 )  // slv_reg_66  <= to_slv_v( buf_abs_tvalid       , C_PSL );
#define NSDET_BUF_ABS_TREADY_REG       ( 67 )  // slv_reg_67  <= to_slv_v( buf_abs_tready       , C_PSL );
#define NSDET_PS_SAMPLE_MAP_0_REG      ( 68 )  // slv_reg_68  <= to_slv_v( ps_sample_map_0      , C_PSL );
#define NSDET_PS_SAMPLE_MAP_1_REG      ( 69 )  // slv_reg_69  <= to_slv_v( ps_sample_map_1      , C_PSL );
#define NSDET_PS_SAMPLE_MAP_2_REG      ( 70 )  // slv_reg_70  <= to_slv_v( ps_sample_map_2      , C_PSL );
#define NSDET_PS_SAMPLE_MAP_3_REG      ( 71 )  // slv_reg_71  <= to_slv_v( ps_sample_map_3      , C_PSL );
#define NSDET_EX_TIMER_REG             ( 72 )  // slv_reg_72  <= to_slv_u( ex_timer             , C_PSL );

#define NSDET_SYMB_DTA_IDX_REG         ( 85 )  // slv_reg_85  <= to_slv_ul( nsd_symb_dta_idx    , C_PSL );
                                               // slv_reg_86  <= to_slv_uh( nsd_symb_dta_idx    , C_PSL );
#define NSDET_CHAN_DET_IDX_REG         ( 87 )  // slv_reg_87  <= to_slv_u( nsd_chan_det_idx     , C_PSL );
#define NSDET_SUM_NOISE_IDX_REG        ( 88 )  // slv_reg_88  <= to_slv_u( nsd_sum_noise_idx    , C_PSL );
#define NSDET_SUM_PEAKS_IDX_REG        ( 89 )  // slv_reg_89  <= to_slv_u( nsd_sum_peaks_idx    , C_PSL );
#define NSDET_RX_RXD_CTR_REG           ( 90 )  // slv_reg_90  <= to_slv_ul( rx_rxd_ctr          , C_PSL );
                                               // slv_reg_91  <= to_slv_uh( rx_rxd_ctr          , C_PSL );
#define NSDET_RX_FFT_CTR_REG           ( 92 )  // slv_reg_92  <= to_slv_ul( rx_fft_ctr          , C_PSL );
                                               // slv_reg_93  <= to_slv_uh( rx_fft_ctr          , C_PSL );
#define NSDET_TX_FFT_CTR_REG           ( 94 )  // slv_reg_94  <= to_slv_ul( tx_fft_ctr          , C_PSL );
                                               // slv_reg_95  <= to_slv_uh( tx_fft_ctr          , C_PSL );
#define NSDET_TX_PRS_CTR_REG           ( 96 )  // slv_reg_96  <= to_slv_u( tx_prs_ctr           , C_PSL );
#define NSDET_TX_DTA_CTR_REG           ( 97 )  // slv_reg_97  <= to_slv_u( tx_dta_ctr           , C_PSL );
#define NSDET_TX_SYM_CTR_REG           ( 98 )  // slv_reg_98  <= to_slv_u( tx_sym_ctr           , C_PSL );
#define NSDET_TX_PRS_START_REG         ( 99 )  // slv_reg_99  <= to_slv_ul( tx_prs_start        , C_PSL );
                                               // slv_reg_100 <= to_slv_uh( tx_prs_start        , C_PSL );
#define NSDET_TX_DTA_START_REG         ( 101 ) // slv_reg_101 <= to_slv_ul( tx_dta_start        , C_PSL );
                                               // slv_reg_102 <= to_slv_uh( tx_dta_start        , C_PSL );


#define NSDET_TST_TM_OFF                 (  0 )
#define NSDET_TST_TM_DET_COUNT           (  1 )

#define NSDET_FFT_CFG_DONE_MASK          ( 0x00000004 )
#define NSDET_NSD_DET_DONE_MASK          ( 0x00000010 )
#define NSDET_TXD_TX_DONE_MASK           ( 0x00000020 )
#define NSDET_BUF_ABS_TVALID_F_MASK      ( 0x00000080 )

#define NSDET_S_AXIS_RXD_TVALID_MASK     ( 0x00000100 )
#define NSDET_S_AXIS_RXD_TREADY_MASK     ( 0x00000200 )
#define NSDET_M_AXIS_FFT_TVALID_MASK     ( 0x00000400 )
#define NSDET_M_AXIS_FFT_TREADY_MASK     ( 0x00000800 )
#define NSDET_M_AXIS_PRS_TVALID_MASK     ( 0x00001000 )
#define NSDET_M_AXIS_PRS_TREADY_MASK     ( 0x00002000 )
#define NSDET_M_AXIS_DTA_TVALID_MASK     ( 0x00004000 )
#define NSDET_M_AXIS_DTA_TREADY_MASK     ( 0x00008000 )

#define NSDET_EX_STATE_MASK              ( 0x000F0000 )
#define NSDET_EX_STATE_SHR               ( 16 )
#define NSDET_NSD_STATE_MASK             ( 0x00F00000 )
#define NSDET_NSD_STATE_SHR              ( 20 )
#define NSDET_TX_STATE_MASK              ( 0x0F000000 )
#define NSDET_TX_STATE_SHR               ( 24 )
#define NSDET_SYS_RESETN_MASK            ( 0x10000000 )
#define NSDET_FIFO_RESETN_MASK           ( 0x20000000 )
#define NSDET_FFT_RESETN_MASK            ( 0x40000000 )
#define NSDET_FFT_CFG_RQST_MASK          ( 0x80000000 )



#define NSDET_EX_IDLE                  ( 0x9 )
#define NSDET_EX_WAIT_NSD              ( 0xB )
#define NSDET_EX_WAIT_FFT_RST          ( 0xE )
#define NSDET_EX_WAIT_FFT_CFG          ( 0xF )
#define NSDET_EX_WAIT_FFT_CFG_ACK      ( 0xD )
#define NSDET_EX_TX_PRS                ( 0xC )
#define NSDET_EX_TX_DATA               ( 0x8 )

#define NSDET_NSD_IDLE                 ( 0x9 )
#define NSDET_NSD_WAIT_PASS_2          ( 0xB )
#define NSDET_NSD_COMP_DATA            ( 0xF )
#define NSDET_NSD_CALC_FRAME_IDX       ( 0xD )

#define NSDET_TX_IDLE                  ( 0x9 )
#define NSDET_TX_FFT                   ( 0xA )
#define NSDET_TX_WAIT_PRS              ( 0xC )
#define NSDET_TX_PRS                   ( 0xD )
#define NSDET_TX_WAIT_DTA              ( 0xE )
#define NSDET_TX_DTA                   ( 0xF )

#define NSDET_RESET                    ( 0x00000001 )
#define NSDET_INT_ACK                  ( 0x00000002 )
#define NSDET_START                    ( 0x00000004 )
#define NSDET_STOP                     ( 0x00000008 )
#define NSDET_SYNC                     ( 0x00000010 )

#define NSDET_HW_INFO_REG              ( 0x05D05D80 )

//typedef enum e_nullSymbDetErr  { NSDET_NAME_TOO_LONG,
//	                             NSDET_DEVICE_ID_FAIL,
//	                             NSDET_FFT_INIT_FAIL,
//                                 NSDET_SUCCESS
//                               } nullSymbDetErr;

typedef enum __nsdtstmode {
	NSDET_TM_OFF = 0,
	NSDET_TM_DET_COUNT = 1,
	NSDET_TM_ONE_SHOT = 2
} NSDetTstMode;

#endif /* NULL_SYMBOL_DETECTOR_H_ */
