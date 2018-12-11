/*
 * TxRx_Emulator.h
 *
 *  Created on: 1 May 2018
 *      Author: JULIAN MORTIMER
 */

#ifndef TR_EMULATOR_H_
#define TR_EMULATOR_H_

#include "TxRxEmulsw.h"

#define TREMUL_DEVICE_ID        ( 0x000c0c20 )
#define TREMUL_SPAN_MAX         ( 4096 )
#define TREMUL_LUT_SIZE         ( 4096 )
#define TREMUL_BIT_MASK         ( TREMUL_LUT_SIZE - 1 )
#define TREMUL_LUT_BITS         ( 12 )
#define TREMUL_FFT_RES          ( 12 )


#define TREMUL_CTRL_REG              ( 0 )
#define TREMUL_INFO_REG              ( 1 )
#define TREMUL_N_FFT_REG             ( 2 )
#define TREMUL_N_CPRE_REG            ( 3 )
#define TREMUL_N_SYMB_FRAME_REG      ( 4 )
#define TREMUL_N_FRAME_BLOCK_REG     ( 5 )
#define TREMUL_AD_RD_WR_ADDR_REG     ( 6 )
#define TREMUL_AD_WR_DATA_REG        ( 7 )

#define TREMUL_STAT_REG              ( 16 )
#define TREMUL_AD_RD_DATA_REG        ( 17 )

#define TREMUL_RX_DATA_CTR_REG       ( 20 )
#define TREMUL_RX_SYMBOL_CTR_REG     ( 21 )
#define TREMUL_RX_FRAME_CTR_REG      ( 22 )
#define TREMUL_RX_BLOCK_CTR_REG      ( 23 )
#define TREMUL_TX_DATA_CTR_REG       ( 24 )
#define TREMUL_TX_SYMBOL_CTR_REG     ( 25 )
#define TREMUL_TX_FRAME_CTR_REG      ( 26 )
#define TREMUL_TX_BLOCK_CTR_REG      ( 27 )

// CTRL_REG Bits
#define TREMUL_BIT_SOFT_RESET        ( 0 )
#define TREMUL_BIT_INT_ACK           ( 1 )
#define TREMUL_BIT_START             ( 2 )
#define TREMUL_BIT_STOP              ( 3 )
#define TREMUL_BIT_ENTER_ADMIN       ( 4 )
#define TREMUL_BIT_EXIT_ADMIN        ( 5 )
#define TREMUL_BIT_AD_WR_LUT         ( 6 )
#define TREMUL_BIT_AD_RD_LUT         ( 7 )
#define TREMUL_BIT_AD_RD_LUT_ACK     ( 8 )

// STATUS_REG Bits
#define TREMUL_BIT_PS_RESET          ( 0 )
#define TREMUL_BIT_IRQ               ( 1 )
#define TREMUL_BIT_TREMUL_MODE       ( 2 )
#define TREMUL_BIT_AD_LUT_RD_RDY     ( 3 )

// Bit masks
#define TREMUL_RESET                 ( 1 << TREMUL_BIT_SOFT_RESET )
#define TREMUL_INT_ACK               ( 1 << TREMUL_BIT_INT_ACK )
#define TREMUL_START                 ( 1 << TREMUL_BIT_START )
#define TREMUL_STOP                  ( 1 << TREMUL_BIT_STOP )
#define TREMUL_ENTER_ADMIN           ( 1 << TREMUL_BIT_ENTER_ADMIN )
#define TREMUL_EXIT_ADMIN            ( 1 << TREMUL_BIT_EXIT_ADMIN )
#define TREMUL_AD_WR_LUT             ( 1 << TREMUL_BIT_AD_WR_LUT )
#define TREMUL_AD_RD_LUT             ( 1 << TREMUL_BIT_AD_RD_LUT )
#define TREMUL_AD_RD_LUT_ACK         ( 1 << TREMUL_BIT_AD_RD_LUT_ACK )

#define TREMUL_PS_RESET              ( 1 << TREMUL_BIT_PS_RESET )
#define TREMUL_AD_LUT_RD_RDY         ( 1 << TREMUL_BIT_AD_LUT_RD_RDY )

#define TREMUL_BIT_MASK_EXEC_STATE   ( 0x00030000 )
#define TREMUL_BITS_EXEC_STATE_IDLE  ( 0x00000000 )
#define TREMUL_BITS_EXEC_STATE_ADMIN ( 0x00010000 )
#define TREMUL_BITS_EXEC_STATE_RUN   ( 0x00020000 )

#define TREMUL_HW_INFO_REG           ( 0x000C0C20 )


#endif /* TXRX_EMULATOR_H_ */
