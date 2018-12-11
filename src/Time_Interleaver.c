/*
 * Time_Interleaver.c
 *
 *  Created on: 31 May 2017
 *      Author: JULIAN MORTIMER
 */
#include "include.h"

//u32 timeInterleaverMap[TIME_INTERLEAVER_BUF_SIZE];
//
//void timeInterleaverReset(void);
//void timeInterleaverIntAck(void);
//void timeInterleaverStart(void);
//void timeInterleaverSetNWordsRx(u32 nwords);
//void timeInterleaverSetNWordsTx(u32 nwords);
//void timeInterleaverSetRxAddr(u32 nwords);
//void timeInterleaverSetTxAddr(u32 nwords);
//
//void timeInterleaverEnterAdmin(void);
//void timeInterleaverExitAdmin(void);
//void timeInterleaverAdWriteIntlvr(u32 addr, u32 data);
//void timeInterleaverAdWriteDataBuf(u32 addr, u32 data);
//u32 timeInterleaverAdReadIntlvr(u32 addr);
//u32 timeInterleaverAdReadDataBuf(u32 addr);
//void timeInterleaverAdInitIntlvr(u32 addr, u32 count, u32 val);
//void timeInterleaverAdLoadIntlvr(u32 addr, u32 count, u32* data);
//void timeInterleaverAdInitDataBuf(u32 addr, u32 count, u32 val);
//void timeInterleaverAdLoadDataBuf(u32 addr, u32 count, u32* data);
//
//void timeInterleaverWriteSlaveReg(u32 reg, u32 data);
//u32 timeInterleaverReadSlaveReg(u32 reg);
//u32 timeInterleaverReadStatusReg(void);
//
//void timeInterleaverInitMap(void);
//
//void timeInterleaverInitMap(void)
//{
//	u32 i = 0;
//	for (i = 0; i < TIME_INTERLEAVER_BUF_SIZE; i++)
//		timeInterleaverMap[i] = i;
//}
//
//void timeInterleaverReset(void)
//{
//    hwWriteCtrlReg(TIME_INTERLEAVER_BASE_ADDR, TIME_INTERLEAVER_CTRL_REG << 2, TIME_INTERLEAVER_SOFT_RESET);
//    while((timeInterleaverReadStatusReg() & TIME_INTERLEAVER_PS_RESET) != 0)
//	   ;
//}
//
//void timeInterleaverIntAck(void)
//{
//    hwWriteCtrlReg(TIME_INTERLEAVER_BASE_ADDR, TIME_INTERLEAVER_CTRL_REG << 2, TIME_INTERLEAVER_INT_ACK);
//}
//
//void timeInterleaverStart(void)
//{
//    hwWriteCtrlReg(TIME_INTERLEAVER_BASE_ADDR, TIME_INTERLEAVER_CTRL_REG << 2, TIME_INTERLEAVER_START);
//}
//
//void timeInterleaverEnterAdmin(void)
//{
//    hwWriteCtrlReg(TIME_INTERLEAVER_BASE_ADDR, TIME_INTERLEAVER_CTRL_REG << 2, TIME_INTERLEAVER_ENTER_ADMIN);
//}
//
//void timeInterleaverExitAdmin(void)
//{
//    hwWriteCtrlReg(TIME_INTERLEAVER_BASE_ADDR, TIME_INTERLEAVER_CTRL_REG << 2, TIME_INTERLEAVER_EXIT_ADMIN);
//}
//
//void timeInterleaverAdWriteIntlvr(u32 addr, u32 data)
//{
//	timeInterleaverEnterAdmin();
//    hwWriteData(TIME_INTERLEAVER_BASE_ADDR,
//    		    TIME_INTERLEAVER_CTRL_REG << 2, TIME_INTERLEAVER_AD_WR_INTLVR,
//    		    TIME_INTERLEAVER_AD_RD_WR_ADDR_REG << 2, addr,
//    		    TIME_INTERLEAVER_AD_WR_DATA_REG << 2, data
//    		    );
//    timeInterleaverExitAdmin();
//}
//
//void timeInterleaverAdWriteDataBuf(u32 addr, u32 data)
//{
//	timeInterleaverEnterAdmin();
//    hwWriteData(TIME_INTERLEAVER_BASE_ADDR,
//    		    TIME_INTERLEAVER_CTRL_REG << 2, TIME_INTERLEAVER_AD_WR_RX_BUF,
//    		    TIME_INTERLEAVER_AD_RD_WR_ADDR_REG << 2, addr,
//    		    TIME_INTERLEAVER_AD_WR_DATA_REG << 2, data
//    		    );
//    timeInterleaverExitAdmin();
//}
//
//u32 timeInterleaverAdReadIntlvr(u32 addr)
//{
//    u32 u = 0;
//    timeInterleaverEnterAdmin();
//    u = hwReadData(TIME_INTERLEAVER_BASE_ADDR,
//                   TIME_INTERLEAVER_CTRL_REG << 2, TIME_INTERLEAVER_AD_RD_INTLVR, TIME_INTERLEAVER_AD_RD_INTLVR_ACK,
//                   TIME_INTERLEAVER_STATUS_REG << 2, TIME_INTERLEAVER_AD_INTLVR_RD_RDY,
//                   TIME_INTERLEAVER_AD_RD_WR_ADDR_REG << 2, addr,
//                   TIME_INTERLEAVER_AD_INTLVR_RD_REG << 2
//                   );
//    timeInterleaverExitAdmin();
//    return u;
//}
//
//u32 timeInterleaverAdReadDataBuf(u32 addr)
//{
//    u32 u = 0;
//    timeInterleaverEnterAdmin();
//    u = hwReadData(TIME_INTERLEAVER_BASE_ADDR,
//                   TIME_INTERLEAVER_CTRL_REG << 2, TIME_INTERLEAVER_AD_RD_TX_BUF, TIME_INTERLEAVER_AD_RD_TX_BUF_ACK,
//                   TIME_INTERLEAVER_STATUS_REG << 2, TIME_INTERLEAVER_AD_TX_BUF_RD_RDY,
//                   TIME_INTERLEAVER_AD_RD_WR_ADDR_REG << 2, addr,
//                   TIME_INTERLEAVER_AD_TX_BUF_RD_REG << 2
//                   );
//    timeInterleaverExitAdmin();
//    return u;
//}
//
//void timeInterleaverAdInitIntlvr(u32 addr, u32 count, u32 val)
//{
//    u32 i = 0;
//    timeInterleaverEnterAdmin();
//	for(i = 0; i < count; i++)
//	{
//	    hwWriteData(TIME_INTERLEAVER_BASE_ADDR,
//	    		    TIME_INTERLEAVER_CTRL_REG << 2, TIME_INTERLEAVER_AD_WR_INTLVR,
//	    		    TIME_INTERLEAVER_AD_RD_WR_ADDR_REG << 2, addr + i,
//	    		    TIME_INTERLEAVER_AD_WR_DATA_REG << 2, val
//	    		    );
//	}
//    timeInterleaverExitAdmin();
//}
//
//void timeInterleaverAdLoadIntlvr(u32 addr, u32 count, u32* data)
//{
//    u32 i = 0;
//    timeInterleaverEnterAdmin();
//	for(i = 0; i < count; i++)
//	{
//	    hwWriteData(TIME_INTERLEAVER_BASE_ADDR,
//	    		    TIME_INTERLEAVER_CTRL_REG << 2, TIME_INTERLEAVER_AD_WR_INTLVR,
//	    		    TIME_INTERLEAVER_AD_RD_WR_ADDR_REG << 2, addr + i,
//	    		    TIME_INTERLEAVER_AD_WR_DATA_REG << 2, data[i]
//	    		    );
//	}
//    timeInterleaverExitAdmin();
//}
//
//void timeInterleaverAdInitDataBuf(u32 addr, u32 count, u32 val)
//{
//    u32 i = 0;
//    timeInterleaverEnterAdmin();
//	for(i = 0; i < count; i++)
//	{
//	    hwWriteData(TIME_INTERLEAVER_BASE_ADDR,
//	    		    TIME_INTERLEAVER_CTRL_REG << 2, TIME_INTERLEAVER_AD_WR_RX_BUF,
//	    		    TIME_INTERLEAVER_AD_RD_WR_ADDR_REG << 2, addr + i,
//	    		    TIME_INTERLEAVER_AD_WR_DATA_REG << 2, val
//	    		    );
//	}
//    timeInterleaverExitAdmin();
//}
//
//void timeInterleaverAdLoadDataBuf(u32 addr, u32 count, u32* data)
//{
//    u32 i = 0;
//    timeInterleaverEnterAdmin();
//	for(i = 0; i < count; i++)
//	{
//	    hwWriteData(TIME_INTERLEAVER_BASE_ADDR,
//	    		    TIME_INTERLEAVER_CTRL_REG << 2, TIME_INTERLEAVER_AD_WR_RX_BUF,
//	    		    TIME_INTERLEAVER_AD_RD_WR_ADDR_REG << 2, addr + i,
//	    		    TIME_INTERLEAVER_AD_WR_DATA_REG << 2, data[i]
//	    		    );
//
//	}
//    timeInterleaverExitAdmin();
//}
//
//void timeInterleaverWriteSlaveReg(u32 reg, u32 data)
//{
//	hwWriteReg(TIME_INTERLEAVER_BASE_ADDR, reg << 2, data);
//}
//
//u32 timeInterleaverReadSlaveReg(u32 reg)
//{
//	return hwReadReg(TIME_INTERLEAVER_BASE_ADDR, reg << 2);
//}
//
//u32 timeInterleaverReadStatusReg(void)
//{
//	return timeInterleaverReadSlaveReg(TIME_INTERLEAVER_STATUS_REG);
//}
//
//void timeInterleaverSetNWordsRx(u32 nwords)
//{
//	timeInterleaverWriteSlaveReg(TIME_INTERLEAVER_N_DATA_WORDS_RX_REG, nwords);
//}
//
//void timeInterleaverSetNWordsTx(u32 nwords)
//{
//	timeInterleaverWriteSlaveReg(TIME_INTERLEAVER_N_DATA_WORDS_TX_REG, nwords);
//}
//
//void timeInterleaverSetRxAddr(u32 addr)
//{
//	timeInterleaverWriteSlaveReg(TIME_INTERLEAVER_RX_START_ADDR_REG, addr);
//}
//
//void timeInterleaverSetTxAddr(u32 addr)
//{
//	timeInterleaverWriteSlaveReg(TIME_INTERLEAVER_TX_START_ADDR_REG, addr);
//}
//
