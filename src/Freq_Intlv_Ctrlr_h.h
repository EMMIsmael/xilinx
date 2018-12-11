/*
 * Freq_Intlv_Ctrlr_h.h
 *
 *  Created on: 21 Jun 2017
 *      Author: JULIAN MORTIMER
 */

#ifndef FREQ_INTLV_CTRLR_H_H_
#define FREQ_INTLV_CTRLR_H_H_
//
//
//#define FREQ_INTLV_CTRLR_RX_BUF_WORD_WIDTH    ( 32 )
//#define FREQ_INTLV_CTRLR_RX_BUF_SIZE          ( 512 )
//#define FREQ_INTLV_CTRLR_TX_BUF_SIZE          ( 16384 )
//#define FREQ_INTLV_CTRLR_INTLVR_SIZE          ( SYS_MAX_CARRIERS )
//#define FREQ_INTLV_CTRLR_MOD_SIZE             ( 1 << SYS_MAX_BITS_PER_CARRIER )
//#define FREQ_INTLVR_N_BITS_PER_SYMBOL         ( SYS_MAX_BITS_PER_CARRIER )
//
//#define FREQ_INTLV_CTRLR_N_SLAVE_REGS         ( 32 )
//
//#define FREQ_INTLV_CTRLR_BASE_ADDR            ( 0x00000000 )
//
//#define FREQ_INTLV_CTRLR_CTRL_REG             ( 0 )
//#define FREQ_INTLV_CTRLR_N_DATA_WORDS_RX_REG  ( 1 )
//#define FREQ_INTLV_CTRLR_N_DATA_WORDS_TX_REG  ( 2 )
//#define FREQ_INTLV_CTRLR_RX_START_ADDR_REG    ( 3 )
//#define FREQ_INTLV_CTRLR_TX_START_ADDR_REG    ( 4 )
//#define FREQ_INTLV_CTRLR_RX_BITS_PER_BYTE_REG ( 5 )
//#define FREQ_INTLV_CTRLR_AD_RD_WR_ADDR_REG    ( 6 )
//#define FREQ_INTLV_CTRLR_AD_WR_DATA_REG       ( 7 )
//#define FREQ_INTLV_CTRLR_AD_BIT_SEL_REG       ( 8 )
//
//#define FREQ_INTLV_CTRLR_STATUS_REG           ( 16 )
//#define FREQ_INTLV_CTRLR_AD_INTLVR_RD_REG_0   ( 17 )
//#define FREQ_INTLV_CTRLR_AD_INTLVR_RD_REG_1   ( 18 )
//#define FREQ_INTLV_CTRLR_AD_INTLVR_RD_REG_2   ( 19 )
//#define FREQ_INTLV_CTRLR_AD_INTLVR_RD_REG_3   ( 20 )
//#define FREQ_INTLV_CTRLR_AD_INTLVR_RD_REG_4   ( 21 )
//#define FREQ_INTLV_CTRLR_AD_INTLVR_RD_REG_5   ( 22 )
//#define FREQ_INTLV_CTRLR_AD_TX_BUF_RD_REG     ( 23 )
//#define FREQ_INTLV_CTRLR_AD_MOD_RD_REG        ( 24 )
//
//#define FREQ_INTLV_CTRLR_SOFT_RESET           ( 0x00000001 )
//#define FREQ_INTLV_CTRLR_INT_ACK              ( 0x00000002 )
//#define FREQ_INTLV_CTRLR_START                ( 0x00000004 )
//#define FREQ_INTLV_CTRLR_ENTER_ADMIN          ( 0x00000008 )
//#define FREQ_INTLV_CTRLR_EXIT_ADMIN           ( 0x00000010 )
//#define FREQ_INTLV_CTRLR_AD_WR_INTLVR         ( 0x00000020 )
//#define FREQ_INTLV_CTRLR_AD_WR_RX_BUF         ( 0x00000040 )
//#define FREQ_INTLV_CTRLR_AD_WR_MOD            ( 0x00000080 )
//#define FREQ_INTLV_CTRLR_AD_RD_INTLVR         ( 0x00000100 )
//#define FREQ_INTLV_CTRLR_AD_RD_TX_BUF         ( 0x00000200 )
//#define FREQ_INTLV_CTRLR_AD_RD_MOD            ( 0x00000400 )
//#define FREQ_INTLV_CTRLR_AD_RD_INTLVR_ACK     ( 0x00000800 )
//#define FREQ_INTLV_CTRLR_AD_RD_TX_BUF_ACK     ( 0x00001000 )
//#define FREQ_INTLV_CTRLR_AD_RD_MOD_ACK        ( 0x00002000 )
//
//
//#define FREQ_INTLV_CTRLR_PS_RESET             ( 0x00000001 )
//#define FREQ_INTLV_CTRLR_IRQ                  ( 0x00000002 )
//#define FREQ_INTLV_CTRLR_RX_SYMBOL_DONE       ( 0x00000004 )
//#define FREQ_INTLV_CTRLR_TX_SYMBOL_DONE       ( 0x00000008 )
//#define FREQ_INTLV_CTRLR_AD_INTLVR_RD_RDY     ( 0x00000010 )
//#define FREQ_INTLV_CTRLR_AD_TX_BUF_RD_RDY     ( 0x00000020 )
//#define FREQ_INTLV_CTRLR_AD_MOD_RD_RDY        ( 0x00000040 )
//
//
//
//typedef struct sFreqIntlvCtrlr tFreqIntlvCtrlr;
//
//struct sFreqIntlvCtrlr {
//	SysCtrl *pSysCtrl;
//	char name[SYS_MAX_OBJ_NAME_LEN+1];
//	u32 baseAddr;
//	u32 id;
//	u32 nChanSymb;
//	u32 nChanUsed;
//	u32 nChanBits;
//	u32 nSymbBits;
//	u32 nRxWordBits;
//	u32 nRxWordsSymb;
//	u32 nRxBitsUsed;
//	u32 nRxWordPack;
//	u32 nRxBitsPack;
//	const u32 *txChanMask;
//	u32 *modulator;
//	u32 (*intlvr)[FREQ_INTLV_CTRLR_INTLVR_SIZE];
//};
//
//

#endif /* FREQ_INTLV_CTRLR_H_H_ */
