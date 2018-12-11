/*
 * TxRxEmulSw.h
 *
 *  Created on: 2 Aug 2018
 *      Author: JULIAN MORTIMER
 */

#ifndef TXRXEMULSW_H_
#define TXRXEMULSW_H_

#define TR_EMUL_SWITCH_DEVICE_ID ( 0xB19A5520 )
#define TR_EMUL_SWITCH_N_SLV_REG ( 32 )

#define TREMULSW_CTRL_REG              ( 0 )
#define TREMULSW_INFO_REG              ( 1 )
#define TREMULSW_N_FFT_REG             ( 2 )
#define TREMULSW_N_CPRE_REG            ( 3 )
#define TREMULSW_N_SYMB_FRAME_REG      ( 4 )
#define TREMULSW_N_FRAME_BLOCK_REG     ( 5 )
#define TREMULSW_CONFIG_REG            ( 7 )

#define TREMULSW_STAT_REG              ( 16 )

#define TREMULSW_RX_DATA_CTR_REG       ( 20 )
#define TREMULSW_RX_SYMBOL_CTR_REG     ( 21 )
#define TREMULSW_RX_FRAME_CTR_REG      ( 22 )
#define TREMULSW_RX_BLOCK_CTR_REG      ( 23 )
#define TREMULSW_TX_DATA_CTR_REG       ( 24 )
#define TREMULSW_TX_SYMBOL_CTR_REG     ( 25 )
#define TREMULSW_TX_FRAME_CTR_REG      ( 26 )
#define TREMULSW_TX_BLOCK_CTR_REG      ( 27 )

// CTRL_REG Bits
#define TREMULSW_BIT_RESET             ( 0 )

// STATUS_REG Bits
#define TREMULSW_BIT_RESET                     (  0 )
#define TREMULSW_BIT_SEL_EMUL                  (  4 )
#define TREMULSW_BIT_FIFO_CLK_LOCK             (  3 )
#define TREMUL_SW_BIT_S_AXIS_IN_TVALID         (  8 )
#define TREMUL_SW_BIT_S_AXIS_IN_TREADY         (  9 )
#define TREMUL_SW_BIT_M_AXIS_OUT_TVALID        ( 10 )
#define TREMUL_SW_BIT_M_AXIS_OUT_TREADY        ( 11 )
#define TREMUL_SW_BIT_S_AXIS_DEVICE_OUT_TVALID ( 12 )
#define TREMUL_SW_BIT_S_AXIS_DEVICE_OUT_TREADY ( 13 )
#define TREMUL_SW_BIT_M_AXIS_DEVICE_IN_TVALID  ( 14 )
#define TREMUL_SW_BIT_M_AXIS_DEVICE_IN_TREADY  ( 15 )
#define TREMUL_SW_BIT_S_AXIS_EMUL_OUT_TVALID   ( 16 )
#define TREMUL_SW_BIT_S_AXIS_EMUL_OUT_TREADY   ( 17 )
#define TREMUL_SW_BIT_M_AXIS_EMUL_IN_TVALID    ( 18 )
#define TREMUL_SW_BIT_M_AXIS_EMUL_IN_TREADY    ( 19 )


// STATUS_REG Bit masks
#define TREMULSW_RESET                      ( 1 << TREMULSW_BIT_RESET )
#define TREMULSW_MSK_SEL_EMUL               ( 3 << TREMULSW_BIT_SEL_EMUL )
#define TREMULSW_FIFO_CLK_LOCK              ( 1 << TREMULSW_BIT_FIFO_CLK_LOCK )
#define TREMULSW_S_AXIS_IN_TVALID           ( 1 << TREMUL_SW_BIT_S_AXIS_IN_TVALID         )
#define TREMULSW_S_AXIS_IN_TREADY           ( 1 << TREMUL_SW_BIT_S_AXIS_IN_TREADY         )
#define TREMULSW_M_AXIS_OUT_TVALID          ( 1 << TREMUL_SW_BIT_M_AXIS_OUT_TVALID        )
#define TREMULSW_M_AXIS_OUT_TREADY          ( 1 << TREMUL_SW_BIT_M_AXIS_OUT_TREADY        )
#define TREMULSW_S_AXIS_DEVICE_OUT_TVALID   ( 1 << TREMUL_SW_BIT_S_AXIS_DEVICE_OUT_TVALID )
#define TREMULSW_S_AXIS_DEVICE_OUT_TREADY   ( 1 << TREMUL_SW_BIT_S_AXIS_DEVICE_OUT_TREADY )
#define TREMULSW_M_AXIS_DEVICE_IN_TVALID    ( 1 << TREMUL_SW_BIT_M_AXIS_DEVICE_IN_TVALID  )
#define TREMULSW_M_AXIS_DEVICE_IN_TREADY    ( 1 << TREMUL_SW_BIT_M_AXIS_DEVICE_IN_TREADY  )
#define TREMULSW_S_AXIS_EMUL_OUT_TVALID     ( 1 << TREMUL_SW_BIT_S_AXIS_EMUL_OUT_TVALID   )
#define TREMULSW_S_AXIS_EMUL_OUT_TREADY     ( 1 << TREMUL_SW_BIT_S_AXIS_EMUL_OUT_TREADY   )
#define TREMULSW_M_AXIS_EMUL_IN_TVALID      ( 1 << TREMUL_SW_BIT_M_AXIS_EMUL_IN_TVALID    )
#define TREMULSW_M_AXIS_EMUL_IN_TREADY      ( 1 << TREMUL_SW_BIT_M_AXIS_EMUL_IN_TREADY    )

#define TREMULSW_PS_RESET              ( 1 << TREMULSW_BIT_PS_RESET )

// CONFIG_REG Bits
//#define TREMULSW_BIT_CFG_SEL_EMUL      ( 0 )
//#define TREMULSW_CFG_SEL_EMUL          ( 1 << TREMULSW_BIT_CFG_SEL_EMUL )
#define TREMULSW_CFG_MSK_SEL_EMUL        ( 0x00000003 )
#define TREMULSW_CFG_BIT_SEL_EMUL        ( 0x00000000 )

#define TREMULSW_HW_INFO_REG             ( 0xB19A5520 )

typedef enum eselemulchan {
    SEL_EMUL_RF_CHAN = 0,
    SEL_EMUL_EMUL    = 1,
    SEL_EMUL_WIRE    = 2
} SelEmulChan;

#endif /* TXRXEMULSW_H_ */
