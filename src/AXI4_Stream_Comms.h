/*
 * AXI4_Stream_Comms.h
 *
 *  Created on: 15 Jul 2017
 *      Author: JULIAN MORTIMER
 */

#ifndef AXI4_STREAM_COMMS_H_
#define AXI4_STREAM_COMMS_H_

#define AXI4_CLOCK_FREQ        ( 100000000.0 )
#define AXI4_TIMEOUT_SEC   ( 5.0 )
#define AXI4_TIMEOUT_TICKS ( (double) (AXI4_CLOCK_FREQ * AXI4_TIMEOUT_SEC) )

typedef enum axi4StreamError {
	AXI4_TIMEOUT = 0,
	AXI4_SUCCESS
}  axi4StreamError;


#endif /* AXI4_STREAM_COMMS_H_ */
