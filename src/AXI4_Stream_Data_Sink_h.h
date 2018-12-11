/*
 * AXI4_Stream_Data_Sink_h.h
 *
 *  Created on: 16 May 2017
 *      Author: JULIAN MORTIMER
 */

#ifndef AXI4_STREAM_DATA_SINK_H_H_
#define AXI4_STREAM_DATA_SINK_H_H_

#define AXISINK_CTRL_REG                  ( 0 )
#define AXISINK_INFO_REG                  ( 1 )
#define AXISINK_CFG_REG                   ( 2 )
#define AXISINK_ADDR_REG                  ( 3 )
#define AXISINK_WR_DATA_REG               ( 4 )
#define AXISINK_BLOCK_SIZE_REG            ( 5 )
#define AXISINK_RX_ADDR_REG               ( 6 )

#define AXISINK_N_FFT_REG                 ( 8 )
#define AXISINK_N_SYMB_FRAME_REG          ( 9 )
#define AXISINK_N_FRAME_BLOCK_REG         ( 10 )

#define AXISINK_STAT_REG                  ( 16 )
#define AXISINK_RD_DATA_REG               ( 17 )
#define AXISINK_BRAM_SIZE_REG             ( 18 )

#define AXISINK_RX_DATA_CTR_REG           ( 20 )
#define AXISINK_RX_SYMB_CTR_REG           ( 21 )
#define AXISINK_RX_FRAME_CTR_REG          ( 22 )
#define AXISINK_RX_BLOCK_CTR_REG          ( 23 )
#define AXISINK_WR_DATA_CTR_REG           ( 24 )
#define AXISINK_WR_SYMB_CTR_REG           ( 25 )
#define AXISINK_WR_FRAME_CTR_REG          ( 26 )
#define AXISINK_WR_BLOCK_CTR_REG          ( 27 )

#define AXISINK_RESET             ( 0x00000001 )
#define AXISINK_INT_ACK           ( 0x00000002 )
#define AXISINK_START_RX          ( 0x00000004 )
#define AXISINK_START             ( 0x00000008 )
#define AXISINK_STOP              ( 0x00000010 )
#define AXISINK_RD_DATA           ( 0x00000020 )
#define AXISINK_RD_ACK            ( 0x00000040 )
#define AXISINK_LOAD_DATA         ( 0x00000080 )
#define AXISINK_WR_DATA           ( 0x00000100 )
#define AXISINK_END_LOAD          ( 0x00000200 )

#define AXISINK_RESET_P           ( 0x00000001 )
#define AXISINK_IRQ               ( 0x00000002 )
#define AXISINK_RD_DATA_RDY       ( 0x00000004 )
#define AXISINK_RX_STREAM_DONE    ( 0x00000008 )

#define AXISINK_WRITE_INHIBIT     ( 0x00000001 )
#define AXISINK_WRITE_ENABLE      ( ~AXISINK_WRITE_INHIBIT )

#define AXISINK_HW_INFO_REG       ( 0x1BEFED20 )

#define AXISINK_S_AXIS_DATA_TREADY ( 0x00000010 )
#define AXISINK_S_AXIS_DATA_TVALID ( 0x00000020 )
#define AXISINK_EX_STATE_MASK      ( 0x00000F00 )
#define AXISINK_EX_STATE_SHR       ( 8 )
#define AXISINK_RX_STATE_MASK      ( 0x0000F000 )
#define AXISINK_RX_STATE_SHR       ( 12 )
#define AXISINK_RD_DATA_STATE_MASK ( 0x000F0000 )
#define AXISINK_RD_DATA_STATE_SHR  ( 16 )

typedef enum eaxisinkmode {
    AXISINK_ONESHOT = 0,
    AXISINK_CONTIN  = 1
} AxiSinkMode;


typedef enum e_axi4sinkerr        { A4SK_NAME_TOO_LONG,
                                    A4SK_SUCCESS
                                  } axi4SinkErr;


#endif /* AXI4_STREAM_DATA_SINK_H_H_ */
