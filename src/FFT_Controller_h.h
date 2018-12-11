/*
 * FFT_Controller_h.h
 *
 *  Created on: 27 Jun 2017
 *      Author: JULIAN MORTIMER
 */

#ifndef FFT_CONTROLLER_H_H_
#define FFT_CONTROLLER_H_H_


#define FFTC_NFFT_OPTS            ( 14 )
#define FFTC_NFFT_MAX             ( 65536 )

#define FFTC_NFFT_FIXED           ( true  )
#define FFTC_NFFT_PGMBL           ( false )
#define FFTC_SCALE_CHAN           (     8 )
#define FFTC_FWD                  (   0x1 )
#define FFTC_INV                  (   0x0 )
#define FFTC_IS_SCALED            (  true )
#define FFTC_NOT_SCALED           ( false )
#define FFTC_HAS_CP               (  true )
#define FFTC_NO_CP                ( false )

#define FFTC_N_SLAVE_REGS         ( 16 )
#define FFTC_DEVICE_ID            ( 0xFF7CC1 )

#define FFTC_CTRL_REG             ( 0 )
#define FFTC_INFO_REG             ( 1 )
#define FFTC_CFG_FFT_LO_REG       ( 2 )
#define FFTC_CFG_FFT_HI_REG       ( 3 )

#define FFTC_STAT_REG             ( 8 )
#define FFTC_EVENT_REG            ( 9 )
#define FFTC_CFG_REG_LEN_REG      ( 10 )

#define FFTC_RESET                ( 0x00000001 )
#define FFTC_INT_ACK              ( 0x00000002 )
#define FFTC_TX_CFG_FFT           ( 0x00000004 )

#define FFTC_PS_RESET             ( 0x00000001 )
#define FFTC_IRQ                  ( 0x00000002 )
#define FFTC_TX_CFG_FFT_DONE      ( 0x00000004 )

#define FFTC_HW_INFO_REG          ( 0xFF7CC110 )


typedef enum e_fftCtrlError { FFT_NAME_TOO_LONG,
                              FFT_MIN_NOT_SUPPORTED,
                              FFT_MAX_NOT_SUPPORTED,
                              FFT_LENGTH_NOT_SUPPORTED,
                              FFT_LENGTH_MIN_EXCEEDS_MAX,
                              FFT_LENGTH_OUT_OF_LIMITS,
                              FFT_HAS_CP_ERROR,
                              FFT_CP_TOO_LONG,
                              FFT_CP_NOT_SPECIFIED,
                              FFT_FWD_INV_ERROR,
                              FFT_IS_SCALED_ERROR,
                              FFT_CFG_REG_MISMATCH,
                              FFT_MIN_MUST_EQUAL_MAX,
                              FFT_SUCCESS
                            } fftCError;


#endif /* FFT_CONTROLLER_H_H_ */
