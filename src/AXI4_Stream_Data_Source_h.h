/*
 * AXI4_Stream_Data_Source_h.h
 *
 *  Created on: 20 May 2017
 *      Author: JULIAN MORTIMER
 */

#ifndef AXI4_STREAM_DATA_SOURCE_H_H_
#define AXI4_STREAM_DATA_SOURCE_H_H_


#define AXISOURCE_CTRL_REG                  ( 0 )
#define AXISOURCE_INFO_REG                  ( 1 )
#define AXISOURCE_BLOCK_SIZE_REG            ( 2 )
#define AXISOURCE_TX_ADDR_REG               ( 3 )

#define AXISOURCE_N_FFT_REG                 ( 5 )
#define AXISOURCE_N_SYMB_FRAME_REG          ( 6 )
#define AXISOURCE_N_FRAME_BLOCK_REG         ( 7 )

#define AXISOURCE_STAT_REG                  ( 8 )
#define AXISOURCE_BRAM_SIZE_REG             ( 9 )

#define AXISOURCE_TX_DATA_CTR_REG          ( 12 )
#define AXISOURCE_TX_SYMBOL_CTR_REG        ( 13 )
#define AXISOURCE_TX_FRAME_CTR_REG         ( 14 )
#define AXISOURCE_TX_BLOCK_CTR_REG         ( 15 )

#define AXISOURCE_RESET                    ( 0x00000001 )
#define AXISOURCE_INT_ACK                  ( 0x00000002 )
#define AXISOURCE_START_1_SHOT             ( 0x00000004 )
#define AXISOURCE_START_CONTIN             ( 0x00000008 )
#define AXISOURCE_STOP                     ( 0x00000010 )

#define AXISOURCE_TX_STREAM_DONE           ( 0x00000001 )

#define AXISOURCE_N_SLV_REG_MASK           ( 0x000000FF )

#define AXISOURCE_HW_INFO_REG              ( 0xD0FEED10 )

#define AXISOURCE_TF_STATE_MASK            ( 0x00000F00 )
#define AXISOURCE_TF_STATE_SHR             (  8 )
#define AXISOURCE_TM_STATE_MASK            ( 0x0000F000 )
#define AXISOURCE_TM_STATE_SHR             ( 12 )

#define AXISOURCE_TX_FRAME_DONE_MASK       ( 0x00000001 )
#define AXISOURCE_M_AXIS_TVALID_MASK       ( 0x00000010 )
#define AXISOURCE_M_AXIS_TREADY_MASK       ( 0x00000020 )


typedef enum eaxisourcemode {
    AXISOURCE_ONESHOT = 0,
    AXISOURCE_CONTIN  = 1
} AxiSourceMode;


typedef enum e_axi4sourceerr   {
	A4SO_NAME_TOO_LONG,
	A4SO_BRAMC,
	A4SO_SUCCESS
} axi4SourceErr;


#endif /* AXI4_STREAM_DATA_SOURCE_H_H_ */
