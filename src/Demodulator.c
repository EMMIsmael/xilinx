/*
 * Demodulator.c
 *
 *  Created on: 13 Apr 2018
 *      Author: JULIAN MORTIMER
 */

#include "include.h"
#include "Device_Manager.h"


//#define DEMOD_CTRL_REG              ( 0 )
//#define DEMOD_INFO_REG              ( 1 )
//#define DEMOD_N_FFT_REG             ( 2 )
//#define DEMOD_N_SYMBOLS_FRAME_REG   ( 3 )
//#define DEMOD_N_FRAMES_BLOCK_REG    ( 4 )
//#define DEMOD_AD_RD_WR_ADDR_REG     ( 5 )
//#define DEMOD_AD_WR_DATA_REG        ( 6 )
//
//#define DEMOD_STAT_REG              ( 16 )
//#define DEMOD_AD_RD_DATA_REG        ( 17 )
//
//#define DEMOD_RX_DATA_CTR_REG       ( 20 )
//#define DEMOD_RX_SYMBOL_CTR_REG     ( 21 )
//#define DEMOD_RX_FRAME_CTR_REG      ( 22 )
//#define DEMOD_RX_BLOCK_CTR_REG      ( 23 )
//#define DEMOD_TX_DATA_CTR_REG       ( 24 )
//#define DEMOD_TX_SYMBOL_CTR_REG     ( 25 )
//#define DEMOD_TX_FRAME_CTR_REG      ( 26 )
//#define DEMOD_TX_BLOCK_CTR_REG      ( 27 )
//
//// CTRL_REG Bits
//#define DEMOD_BIT_RESET             ( 0 )
//#define DEMOD_BIT_INT_ACK           ( 1 )
//#define DEMOD_BIT_START             ( 2 )
//#define DEMOD_BIT_STOP              ( 3 )
//#define DEMOD_BIT_ENTER_ADMIN       ( 4 )
//#define DEMOD_BIT_EXIT_ADMIN        ( 5 )
//#define DEMOD_BIT_WR_DEMOD          ( 6 )
//#define DEMOD_BIT_RD_DEMOD          ( 7 )
//#define DEMOD_BIT_RD_DEMOD_ACK      ( 8 )
//
//// STATUS_REG Bits
//#define DEMOD_BIT_RESET_ASSERTED    ( 0 )
//#define DEMOD_BIT_IRQ               ( 1 )
//#define DEMOD_BIT_DEMOD_MODE        ( 2 )
//#define DEMOD_BIT_AD_DEMOD_RD_RDY   ( 3 )
//
//#define DEMOD_BIT_MASK_EXEC_STATE   ( 0x00030000 )
//#define DEMOD_BITS_EXEC_STATE_IDLE  ( 0x00000000 )
//#define DEMOD_BITS_EXEC_STATE_ADMIN ( 0x00010000 )
//#define DEMOD_BITS_EXEC_STATE_RUN   ( 0x00020000 )
//
