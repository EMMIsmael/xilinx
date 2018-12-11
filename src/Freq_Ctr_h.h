/*
 * Freq_Ctr_h.h
 *
 *  Created on: 13 Jun 2017
 *      Author: JULIAN MORTIMER
 */

#ifndef FREQ_CTR_H_H_
#define FREQ_CTR_H_H_

#define FREQCTR_CTRL_REG     ( 0 )
#define FREQCTR_INFO_REG     ( 1 )

#define FREQCTR_STAT_REG     ( 4 )
#define FREQCTR_CTR_REG      ( 5 )
#define FREQCTR_TMR_REG      ( 6 )

#define FREQCTR_RESET        ( 0x00000001 )

#define FREQCTR_PS_RESET     ( 0x00000001 )
#define FREQCTR_MMCM_LOCK    ( 0x00000002 )
#define FREQCTR_RST_IDLE     ( 0x00000004 )

#define FREQCTR_HW_INFO_REG  ( 0xFC0FC008 )
#endif /* FREQ_CTR_H_H_ */
