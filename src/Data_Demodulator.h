/*
 * Data_Demodulator.h
 *
 *  Created on: 14 Apr 2018
 *      Author: JULIAN MORTIMER
 */

#ifndef DATA_DEMODULATOR_H_
#define DATA_DEMODULATOR_H_

#define DATA_DEMOD_ADR_BITS           ( 12 )
#define DATA_DEMOD_LUT_BITS           ( 12 )


typedef struct sDemodData DemodData;
struct sDemodData {
    u32 raw;
    u32 word;
    s32 bit0;
    s32 bit1;
    s32 bit2;
    s32 bit3;
    s32 bit4;
    s32 bit5;
};


//
//#define DATA_DEMOD_BASEADDR         ( 0x43CB0000 )
//#define DATA_DEMOD_DEVICE_ID        ( 0x000c0c20 )
//#define DATA_DEMOD_MODE             ( DEMOD_DATA )
//#define DATA_DEMOD_N_FFT            ( 1024 )
//#define DATA_DEMOD_N_SYMBOLS_FRAME  ( 50652 )
//#define DATA_DEMOD_N_FRAMES_BLOCK   ( 30 )
//#define DATA_DEMOD_N_SLV_REG        ( 32 )
//#define DATA_DEMOD_LUT_BITS         ( 12 )
//#define DATA_DEMOD_FFT_RES          ( 12 )
//
//
//#define DATA_DEMOD_TEST_VEC_SIZE    ( 65536 )
//


//#define DDM_CTRL_REG              ( 0 )
//#define DDM_INFO_REG              ( 1 )
//#define DDM_N_FFT_REG             ( 2 )
//#define DDM_N_SYMBOLS_FRAME_REG   ( 3 )
//#define DDM_N_FRAMES_BLOCK_REG    ( 4 )
//#define DDM_AD_RD_WR_ADDR_REG     ( 5 )
//#define DDM_AD_WR_DATA_REG        ( 6 )
//
//#define DDM_STATUS_REG            ( 16 )
//#define DDM_AD_RD_DATA_REG        ( 17 )
//
//#define DDM_RX_DATA_CTR_REG       ( 20 )
//#define DDM_RX_SYMBOL_CTR_REG     ( 21 )
//#define DDM_RX_FRAME_CTR_REG      ( 22 )
//#define DDM_RX_BLOCK_CTR_REG      ( 23 )
//#define DDM_TX_DATA_CTR_REG       ( 24 )
//#define DDM_TX_SYMBOL_CTR_REG     ( 25 )
//#define DDM_TX_FRAME_CTR_REG      ( 26 )
//#define DDM_TX_BLOCK_CTR_REG      ( 27 )
//
//// CTRL_REG Bits
//#define DDM_BIT_SOFT_RESET        ( 0 )
//#define DDM_BIT_INT_ACK           ( 1 )
//#define DDM_BIT_START             ( 2 )
//#define DDM_BIT_STOP              ( 3 )
//#define DDM_BIT_ENTER_ADMIN       ( 4 )
//#define DDM_BIT_EXIT_ADMIN        ( 5 )
//#define DDM_BIT_AD_WR_LUT         ( 6 )
//#define DDM_BIT_AD_RD_LUT         ( 7 )
//#define DDM_BIT_AD_RD_LUT_ACK     ( 8 )
//
//// STATUS_REG Bits
//#define DDM_BIT_PS_RESET          ( 0 )
//#define DDM_BIT_IRQ               ( 1 )
//#define DDM_BIT_DDM_MODE          ( 2 )
//#define DDM_BIT_AD_LUT_RD_RDY     ( 3 )
//
//// Bit masks
//#define DDM_SOFT_RESET            ( 1 << DDM_BIT_SOFT_RESET )
//#define DDM_INT_ACK               ( 1 << DDM_BIT_INT_ACK )
//#define DDM_START                 ( 1 << DDM_BIT_START )
//#define DDM_STOP                  ( 1 << DDM_BIT_STOP        )
//#define DDM_ENTER_ADMIN           ( 1 << DDM_BIT_ENTER_ADMIN )
//#define DDM_EXIT_ADMIN            ( 1 << DDM_BIT_EXIT_ADMIN  )
//#define DDM_AD_WR_LUT             ( 1 << DDM_BIT_AD_WR_LUT      )
//#define DDM_AD_RD_LUT             ( 1 << DDM_BIT_AD_RD_LUT      )
//#define DDM_AD_RD_LUT_ACK         ( 1 << DDM_BIT_AD_RD_LUT_ACK  )
//
//#define DDM_PS_RESET              ( 1 << DDM_BIT_PS_RESET )
//#define DDM_AD_LUT_RD_RDY         ( 1 << DDM_BIT_AD_LUT_RD_RDY )
//
//
//
//#define DDM_BIT_MASK_EXEC_STATE   ( 0x00030000 )
//#define DDM_BITS_EXEC_STATE_IDLE  ( 0x00000000 )
//#define DDM_BITS_EXEC_STATE_ADMIN ( 0x00010000 )
//#define DDM_BITS_EXEC_STATE_RUN   ( 0x00020000 )
//



#endif /* DATA_DEMODULATOR_H_ */
