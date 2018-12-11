/*
 * TS_Injector_h.h
 *
 *  Created on: 6 Jul 2017
 *      Author: JULIAN MORTIMER
 */

#ifndef TS_INJECTOR_H_H_
#define TS_INJECTOR_H_H_

#define TSIPAR( Param ) DEVICE_LOC_PAR( pd, TSInj, Param )
#define TSIDAT( Param ) DEVICE_LOC_DAT( pd, TSInj, Param )


#define TSINJ_BUF_WORD_WIDTH                ( 32 )
#define TSINJ_BUF_SIZE                      ( SYS_MAX_CARRIERS )
#define TSINJ_TSBUF_SIZE                    ( SYS_MAX_CARRIERS * 4 )
#define TSINJ_RXBUF_SIZE                    ( SYS_MAX_CARRIERS * 2 )
#define TSINJ_DOUT_LSB                      ( 0 )
#define TSINJ_N_SLAVE_REGS                  ( 32 )

#define TSINJ_CTRL_REG                      ( 0 )
#define TSINJ_INFO_REG                      ( 1 )
#define TSINJ_N_CPRE_REG                    ( 2 )
#define TSINJ_N_FFT_REG                     ( 3 )
#define TSINJ_N_SYMB_AGC_REG                ( 4 )
#define TSINJ_N_SYMB_NSD_REG                ( 5 )
#define TSINJ_N_SYMB_PRS_REG                ( 6 )
#define TSINJ_N_SYMB_DTA_REG                ( 7 )
#define TSINJ_N_FRAME_BLOCK_REG             ( 8 )
#define TSINJ_AD_RD_WR_ADDR_REG             ( 9 )
#define TSINJ_AD_WR_DATA_REG                ( 10 )


#define TSINJ_SOFT_RESET           ( 0x00000001 )
#define TSINJ_INT_ACK              ( 0x00000004 )
#define TSINJ_START                ( 0x00000008 )
#define TSINJ_STOP                 ( 0x00000010 )
#define TSINJ_ENTER_ADMIN          ( 0x00000020 )
#define TSINJ_EXIT_ADMIN           ( 0x00000040 )
#define TSINJ_AD_WR_TS_BUF         ( 0x00000080 )
#define TSINJ_AD_WR_RX_BUF         ( 0x00000100 )
#define TSINJ_AD_RD_TS_BUF         ( 0x00000200 )
#define TSINJ_AD_RD_TX_BUF         ( 0x00000400 )
#define TSINJ_AD_RD_ACK            ( 0x00000800 )
#define TSINJ_RD_SYNC_F            ( 0x00001000 )

#define TSINJ_STAT_REG                     ( 16 )
#define TSINJ_STATE_REG                    ( 17 )
#define TSINJ_AD_RD_TS_BUF_REG             ( 18 )
#define TSINJ_AD_RD_TX_BUF_REG             ( 19 )
#define TSINJ_RX_DATA_CTR_REG              ( 20 )
#define TSINJ_RX_SYMBOL_CTR_REG            ( 21 )
#define TSINJ_RX_FRAME_CTR_REG             ( 22 )
#define TSINJ_RX_BLOCK_CTR_REG             ( 23 )
#define TSINJ_TX_DATA_CTR_REG              ( 24 )
#define TSINJ_TX_SYMBOL_CTR_REG            ( 25 )
#define TSINJ_TX_FRAME_CTR_REG             ( 26 )
#define TSINJ_TX_BLOCK_CTR_REG             ( 27 )
#define TSINJ_RX_TV_CTR_REG                ( 28 )
#define TSINJ_TX_TV_CTR_REG                ( 29 )
#define TSINJ_RX_DS_ACCUM_REG              ( 32 )
#define TSINJ_TX_TS_ACCUM_REG              ( 33 )
#define TSINJ_TX_DS_ACCUM_REG              ( 34 )
#define TSINJ_RX_FRAME_COUNT_REG           ( 35 )
#define TSINJ_TX_FRAME_COUNT_REG           ( 36 )
#define TSINJ_N_WD_SYMB_RX_REG             ( 37 )
#define TSINJ_N_WD_SYMB_TX_REG             ( 38 )
#define TSINJ_N_SYMB_FRAME_RX_REG          ( 39 )
#define TSINJ_N_SYMB_FRAME_TX_REG          ( 40 )
#define TSINJ_N_SYMB_TRAIN_REG             ( 41 )

#define TSINJ_RX_SYMBOL_DONE_MASK         ( 0x00000008 )
#define TSINJ_TX_SYMBOL_DONE_MASK         ( 0x00000010 )

#define TSINJ_S_AXIS_DATA_TREADY_MASK     ( 0x00000100 )
#define TSINJ_S_AXIS_DATA_TVALID_MASK     ( 0x00000200 )
#define TSINJ_M_AXIS_DATA_TREADY_MASK     ( 0x00000400 )
#define TSINJ_M_AXIS_DATA_TVALID_MASK     ( 0x00000800 )

#define TSINJ_RESET                ( 0x00000001 )
#define TSINJ_IRQ                  ( 0x00000004 )
#define TSINJ_RX_FRAME_DONE        ( 0x00000008 )
#define TSINJ_TX_FRAME_DONE        ( 0x00000010 )
#define TSINJ_AD_RD_RDY            ( 0x00000020 )

#define TSINJ_EX_STATE_MASK              ( 0x0000000F )
#define TSINJ_EX_STATE_SHR               (  0 )
#define TSINJ_PF_STATE_MASK              ( 0x000000F0 )
#define TSINJ_PF_STATE_SHR               (  4 )
#define TSINJ_RM_STATE_MASK              ( 0x00000F00 )
#define TSINJ_RM_STATE_SHR               (  8 )
#define TSINJ_TM_STATE_MASK              ( 0x0000F000 )
#define TSINJ_TM_STATE_SHR               ( 12 )
#define TSINJ_AD_RD_TS_STATE_MASK        ( 0x000F0000 )
#define TSINJ_AD_RD_TS_STATE_SHR         ( 16 )
#define TSINJ_AD_RD_TX_STATE_MASK        ( 0x00F00000 )
#define TSINJ_AD_RD_TX_STATE_SHR         ( 20 )

#define TSINJ_HW_INFO_REG                ( 0x75175120 )

//typedef enum e_tsInjectorError { TSINJ_NAME_TOO_LONG,
//                                 TSINJ_N_FFT_ERROR,
//                                 TSINJ_N_CPRE_ERROR,
//                                 TSINJ_N_DSYMB_ERROR,
//                                 TSINJ_N_TSYMB_ERROR,
//                                 TSINJ_SUCCESS
//                               } TSInjError;

#endif /* TS_INJECTOR_H_H_ */
