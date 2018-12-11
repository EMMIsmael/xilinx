/*
 * PRS_CIR.h
 *
 *  Created on: 12 Sep 2018
 *      Author: JULIAN MORTIMER
 */

#ifndef PRS_CIR_H_
#define PRS_CIR_H_

#define CIRPAR( Param ) DEVICE_LOC_PAR( pd, CIR, Param )
#define CIRDAT( Param ) DEVICE_LOC_DAT( pd, CIR, Param )

#define CIR_N_SLV_REGS                 ( 128 )

#define CIR_CTRL_REG                   (  0 )
#define CIR_N_CPRE_REG                 (  2 )
#define CIR_N_FFT_REG                  (  3 )
#define CIR_N_SYMB_FRAME_REG           (  4 )
#define CIR_N_FRAME_BLOCK_REG          (  5 )
#define CIR_AD_RD_WR_ADDR_REG          (  6 )
#define CIR_AD_WR_DATA_REG             (  7 ) // fft
#define CIR_AD_WR_DATA_REG_H           (  8 ) // abs
#define CIR_LOG_N_FFT_REG              ( 16 )
#define CIR_KFF_REG                    ( 17 )
#define CIR_KTF_REG                    ( 18 )
#define CIR_KTT_REG                    ( 19 )
#define CIR_KFT_REG                    ( 20 )
#define CIR_BRAM_CHAN_SEL_REG          ( 21 )



#define CIR_INFO_REG                   (  1 )
#define CIR_STAT_REG                   ( 64 )
#define CIR_AD_RD_DATA_REG             ( 65 ) // fft
#define CIR_AD_RD_DATA_REG_H           ( 66 ) // abs
#define CIR_HNDSHK_REG                 ( 67 )
#define CIR_AXIS_CIR_TVALID_REG        ( 68 )
#define CIR_AXIS_CIR_TREADY_REG        ( 60 )
#define CIR_AXIS_CIR_ROW_TVALID_REG    ( 70 )
#define CIR_AXIS_CIR_ROW_TREADY_REG    ( 71 )
#define CIR_DELTA_F_REG                ( 72 )
#define CIR_DELTA_T_REG                ( 73 )
#define CIR_DELTAF_RAW_LO_REG          ( 76 )
#define CIR_DELTAF_RAW_HI_REG          ( 77 )
#define CIR_DELTAT_RAW_LO_REG          ( 78 )
#define CIR_DELTAT_RAW_HI_REG          ( 79 )
#define CIR_STATE_REG                  ( 80 )
#define CIR_AM_STATE1_REG              ( 81 )
#define CIR_AM_STATE2_REG              ( 82 )
#define CIR_AM_STATE3_REG              ( 83 )
#define CIR_DO_STATE1_REG              ( 84 )
#define CIR_DO_STATE2_REG              ( 85 )
#define CIR_DO_STATE3_REG              ( 86 )
#define CIR_STRM_STATE1_REG            ( 87 )
#define CIR_STRM_STATE2_REG            ( 88 )
#define CIR_DELTAF_REG                 ( 89 )
#define CIR_DELTAF_C_REG               ( 90 )
#define CIR_DELTAT_REG                 ( 91 )
#define CIR_DELTAT_C_REG               ( 92 )
#define CIR_AXIS_INSTR_REG             ( 95 )
#define CIR_PRS_RX_DATA_CTR_REG        ( 96  )
#define CIR_PRS_RX_SYMB_CTR_REG        ( 97  )
#define CIR_PRS_RX_FRAME_CTR_REG       ( 98  )
#define CIR_PRS_RX_BLOCK_CTR_REG       ( 99  )
#define CIR_PRS_TX_DATA_CTR_REG        ( 100 )
#define CIR_PRS_TX_SYMB_CTR_REG        ( 101 )
#define CIR_PRS_TX_FRAME_CTR_REG       ( 102 )
#define CIR_PRS_TX_BLOCK_CTR_REG       ( 103 )
#define CIR_BPB_RX_DATA_CTR_REG        ( 104 )
#define CIR_BPB_RX_SYMB_CTR_REG        ( 105 )
#define CIR_BPB_RX_FRAME_CTR_REG       ( 106 )
#define CIR_BPB_RX_BLOCK_CTR_REG       ( 107 )
#define CIR_BPB_TX_DATA_CTR_REG        ( 108 )
#define CIR_BPB_TX_SYMB_CTR_REG        ( 109 )
#define CIR_BPB_TX_FRAME_CTR_REG       ( 110 )
#define CIR_BPB_TX_BLOCK_CTR_REG       ( 111 )
#define CIR_BPA_RX_DATA_CTR_REG        ( 112 )
#define CIR_BPA_RX_SYMB_CTR_REG        ( 113 )
#define CIR_BPA_RX_FRAME_CTR_REG       ( 114 )
#define CIR_BPA_RX_BLOCK_CTR_REG       ( 115 )
#define CIR_BPA_TX_DATA_CTR_REG        ( 116 )
#define CIR_BPA_TX_SYMB_CTR_REG        ( 117 )
#define CIR_BPA_TX_FRAME_CTR_REG       ( 118 )
#define CIR_BPA_TX_BLOCK_CTR_REG       ( 119 )
#define CIR_DTA_RX_DATA_CTR_REG        ( 120 )
#define CIR_DTA_RX_SYMB_CTR_REG        ( 121 )
#define CIR_DTA_RX_FRAME_CTR_REG       ( 122 )
#define CIR_DTA_RX_BLOCK_CTR_REG       ( 123 )
#define CIR_DTA_TX_DATA_CTR_REG        ( 124 )
#define CIR_DTA_TX_SYMB_CTR_REG        ( 125 )
#define CIR_DTA_TX_FRAME_CTR_REG       ( 126 )
#define CIR_DTA_TX_BLOCK_CTR_REG       ( 127 )

// Control register bits
#define CIR_RESET                      ( 0x00000001 )
#define CIR_SYS_RESET                  ( 0x00000002 )
#define CIR_INT_ACK                    ( 0x00000004 )
#define CIR_START                      ( 0x00000008 )
#define CIR_ENTER_ADMIN                ( 0x00000010 )
#define CIR_EXIT_ADMIN                 ( 0x00000020 )
#define CIR_AD_RD_BRAM                 ( 0x00000040 )
#define CIR_AD_RD_ACK                  ( 0x00000080 )
#define CIR_AD_WR_BRAM                 ( 0x00000100 )

// Status register bits
#define CIR_AD_BRAM_RD_RDY             ( 0x00000002 )
#define CIR_EXT_RESET_DONE             ( 0x00000008 )
#define CIR_SYS_RESETN_I               ( 0x00000010 )
#define CIR_SYS_RESETP_I               ( 0x00000020 )
#define CIR_DSP_RESETN_I               ( 0x00000040 )
#define CIR_DSP_RESETP_I               ( 0x00000080 )
#define CIR_SYS_RESET_S                ( 0x00000100 )
#define CIR_EX_SYS_RESET               ( 0x00000200 )
#define CIR_PRS_BRAM_PRS_DONE          ( 0x00001000 )
#define CIR_PRS_BRAM_RXD_IN_DONE       ( 0x00002000 )
#define CIR_PRS_BRAM_RXD_OUT_DONE      ( 0x00004000 )
#define CIR_NSD_DONE                   ( 0x00008000 )


#define CIR_HW_INFO_REG                ( 0xCDE17A80 )

#define CIR_STATE_REG_IDLE             ( 0x00999999 )
#define CIR_AM_STATE1_REG_IDLE         ( 0x00099999 )
#define CIR_AM_STATE2_REG_IDLE         ( 0x99999999 )
#define CIR_AM_STATE3_REG_IDLE         ( 0x99999999 )
#define CIR_DO_STATE1_REG_IDLE         ( 0x00099999 )
#define CIR_DO_STATE2_REG_IDLE         ( 0x99999999 )
#define CIR_DO_STATE3_REG_IDLE         ( 0x99999999 )
#define CIR_STRM_STATE1_IDLE           ( 0x00009999 )
#define CIR_STRM_STATE2_IDLE           ( 0x99999999 )

#define CIR_EX_IDLE                    ( 0x00000009 )
#define CIR_EX_ADMIN_WAIT_RST          ( 0x0000000B )
#define CIR_EX_ADMIN                   ( 0x0000000D )
#define CIR_EX_RUN_WAIT_RST            ( 0x0000000E )
#define CIR_EX_RUN                     ( 0x0000000F )

#endif /* PRS_CIR_H_ */
