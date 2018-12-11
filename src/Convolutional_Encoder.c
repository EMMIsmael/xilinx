/*
 * Convolutional_Encoder.c
 *
 *  Created on: 9 Jun 2017
 *      Author: JULIAN MORTIMER
 */

#include "include.h"

//u32 convEncdrMap[CONV_ENCDR_MAP_SIZE];
//u32 convEncdrTstMap[CONV_ENCDR_MAP_SIZE];
//u32 convEncdrBuf[CONV_ENCDR_BUF_SIZE];
//
//void convEncdrReset(void);
//void convEncdrIntAck(void);
//void convEncdrStart(void);
//void convEncdrSetNWordsRx(u32 nwords);
//void convEncdrSetNWordsTx(u32 nwords);
//void convEncdrSetRxAddr(u32 addr);
//void convEncdrSetTxAddr(u32 addr);
//
//void convEncdrEnterAdmin(void);
//void convEncdrExitAdmin(void);
//void convEncdrAdWriteEncdr(u32 addr, u32 data);
//void convEncdrAdWriteDataBuf(u32 addr, u32 data);
//u32 convEncdrAdReadEncdr(u32 addr);
//u32 convEncdrAdReadDataBuf(u32 addr);
//void convEncdrAdInitEncdr(u32 addr, u32 count, u32 val);
//void convEncdrAdLoadEncdr(u32 addr, u32 count, u32* data);
//void convEncdrAdInitDataBuf(u32 addr, u32 count, u32 val);
//void convEncdrAdLoadDataBuf(u32 addr, u32 count, u32* data);
//
//void convEncdrWriteSlaveReg(u32 reg, u32 data);
//u32 convEncdrReadSlaveReg(u32 reg);
//u32 convEncdrReadStatusReg(void);
//
//void convEncdrInitMap(void);
//void convEncdrInitTstMap(u32 *pdata);
//void convEncdrInitBuf(void);
//
//void convEncdrInitMap(void)
//{
//	u32 i = 0;
//	for (i = 0; i < CONV_ENCDR_MAP_SIZE; i++)
//		convEncdrMap[i] = i >> 6; //(((i & 0x3fff) << 2) & 0xfff0) + (i & 0x3f); // ((i << 8) & 0x0000ff00) + (i & 0x3f); // i>>6;
//}
//
////Convolutional encoder map address is input byte << 6 plus last input byte >> 2
////This loads map with data required by following block
//void convEncdrInitTstMap(u32 *pdata)
//{
//	u32 i = 0;
//	u32 *p = pdata;
//	for (i = 0; i < CONV_ENCDR_MAP_SIZE; i++)
//	{
//		convEncdrTstMap[i] = *p;
//		if ((i + 1) % 64 == 0) p++;
//	}
//}
//
//void convEncdrInitBuf(void)
//{
//	u32 i = 0;
//	for (i = 0; i < CONV_ENCDR_BUF_SIZE; i++)
//		convEncdrBuf[i] = i;
//}
//
//void convEncdrReset(void)
//{
//    hwWriteCtrlReg(CONV_ENCDR_BASE_ADDR, CONV_ENCDR_CTRL_REG << 2, CONV_ENCDR_SOFT_RESET);
//    while((convEncdrReadStatusReg() & CONV_ENCDR_PS_RESET) != 0)
//	   ;
//}
//
//void convEncdrIntAck(void)
//{
//    hwWriteCtrlReg(CONV_ENCDR_BASE_ADDR, CONV_ENCDR_CTRL_REG << 2, CONV_ENCDR_INT_ACK);
//}
//
//void convEncdrStart(void)
//{
//    hwWriteCtrlReg(CONV_ENCDR_BASE_ADDR, CONV_ENCDR_CTRL_REG << 2, CONV_ENCDR_START);
//}
//
//void convEncdrEnterAdmin(void)
//{
//    hwWriteCtrlReg(CONV_ENCDR_BASE_ADDR, CONV_ENCDR_CTRL_REG << 2, CONV_ENCDR_ENTER_ADMIN);
//}
//
//void convEncdrExitAdmin(void)
//{
//    hwWriteCtrlReg(CONV_ENCDR_BASE_ADDR, CONV_ENCDR_CTRL_REG << 2, CONV_ENCDR_EXIT_ADMIN);
//}
//
//void convEncdrAdWriteEncdr(u32 addr, u32 data)
//{
//	convEncdrEnterAdmin();
//    hwWriteData(CONV_ENCDR_BASE_ADDR,
//    		    CONV_ENCDR_CTRL_REG << 2, CONV_ENCDR_AD_WR_ENCDR,
//    		    CONV_ENCDR_AD_RD_WR_ADDR_REG << 2, addr,
//    		    CONV_ENCDR_AD_WR_DATA_REG << 2, data
//    		    );
//    convEncdrExitAdmin();
//}
//
//void convEncdrAdWriteDataBuf(u32 addr, u32 data)
//{
//	convEncdrEnterAdmin();
//    hwWriteData(CONV_ENCDR_BASE_ADDR,
//    		    CONV_ENCDR_CTRL_REG << 2, CONV_ENCDR_AD_WR_TX_BUF,
//    		    CONV_ENCDR_AD_RD_WR_ADDR_REG << 2, addr,
//    		    CONV_ENCDR_AD_WR_DATA_REG << 2, data
//    		    );
//    convEncdrExitAdmin();
//}
//
//u32 convEncdrAdReadEncdr(u32 addr)
//{
//    u32 u = 0;
//    convEncdrEnterAdmin();
//    u = hwReadData(CONV_ENCDR_BASE_ADDR,
//                   CONV_ENCDR_CTRL_REG << 2, CONV_ENCDR_AD_RD_ENCDR, CONV_ENCDR_AD_RD_ENCDR_ACK,
//                   CONV_ENCDR_STATUS_REG << 2, CONV_ENCDR_AD_ENCDR_RD_RDY,
//                   CONV_ENCDR_AD_RD_WR_ADDR_REG << 2, addr,
//                   CONV_ENCDR_AD_ENCDR_RD_REG << 2
//                   );
//    convEncdrExitAdmin();
//    return u;
//}
//
//u32 convEncdrAdReadDataBuf(u32 addr)
//{
//    u32 u = 0;
//    convEncdrEnterAdmin();
//    u = hwReadData(CONV_ENCDR_BASE_ADDR,
//                   CONV_ENCDR_CTRL_REG << 2, CONV_ENCDR_AD_RD_TX_BUF, CONV_ENCDR_AD_RD_TX_BUF_ACK,
//                   CONV_ENCDR_STATUS_REG << 2, CONV_ENCDR_AD_TX_BUF_RD_RDY,
//                   CONV_ENCDR_AD_RD_WR_ADDR_REG << 2, addr,
//                   CONV_ENCDR_AD_TX_BUF_RD_REG << 2
//                   );
//    convEncdrExitAdmin();
//    return u;
//}
//
//void convEncdrAdInitEncdr(u32 addr, u32 count, u32 val)
//{
//    u32 i = 0;
//    convEncdrEnterAdmin();
//	for(i = 0; i < count; i++)
//	{
//	    hwWriteData(CONV_ENCDR_BASE_ADDR,
//	    		    CONV_ENCDR_CTRL_REG << 2, CONV_ENCDR_AD_WR_ENCDR,
//	    		    CONV_ENCDR_AD_RD_WR_ADDR_REG << 2, addr + i,
//	    		    CONV_ENCDR_AD_WR_DATA_REG << 2, val
//	    		    );
//	}
//    convEncdrExitAdmin();
//}
//
//void convEncdrAdLoadEncdr(u32 addr, u32 count, u32* data)
//{
//    u32 i = 0;
//    convEncdrEnterAdmin();
//	for(i = 0; i < count; i++)
//	{
//	    hwWriteData(CONV_ENCDR_BASE_ADDR,
//	    		    CONV_ENCDR_CTRL_REG << 2, CONV_ENCDR_AD_WR_ENCDR,
//	    		    CONV_ENCDR_AD_RD_WR_ADDR_REG << 2, addr + i,
//	    		    CONV_ENCDR_AD_WR_DATA_REG << 2, data[i]
//	    		    );
//	}
//    convEncdrExitAdmin();
//}
//
//void convEncdrAdInitDataBuf(u32 addr, u32 count, u32 val)
//{
//    u32 i = 0;
//    convEncdrEnterAdmin();
//	for(i = 0; i < count; i++)
//	{
//	    hwWriteData(CONV_ENCDR_BASE_ADDR,
//	    		    CONV_ENCDR_CTRL_REG << 2, CONV_ENCDR_AD_WR_TX_BUF,
//	    		    CONV_ENCDR_AD_RD_WR_ADDR_REG << 2, addr + i,
//	    		    CONV_ENCDR_AD_WR_DATA_REG << 2, val
//	    		    );
//	}
//    convEncdrExitAdmin();
//}
//
//void convEncdrAdLoadDataBuf(u32 addr, u32 count, u32* data)
//{
//    u32 i = 0;
//    convEncdrEnterAdmin();
//	for(i = 0; i < count; i++)
//	{
//	    hwWriteData(CONV_ENCDR_BASE_ADDR,
//	    		    CONV_ENCDR_CTRL_REG << 2, CONV_ENCDR_AD_WR_TX_BUF,
//	    		    CONV_ENCDR_AD_RD_WR_ADDR_REG << 2, addr + i,
//	    		    CONV_ENCDR_AD_WR_DATA_REG << 2, data[i]
//	    		    );
//	}
//    convEncdrExitAdmin();
//}
//
//void convEncdrWriteSlaveReg(u32 reg, u32 data)
//{
//	hwWriteReg(CONV_ENCDR_BASE_ADDR, reg << 2, data);
//}
//
//u32 convEncdrReadSlaveReg(u32 reg)
//{
//	return hwReadReg(CONV_ENCDR_BASE_ADDR, reg << 2);
//}
//
//u32 convEncdrReadStatusReg(void)
//{
//	return convEncdrReadSlaveReg(CONV_ENCDR_STATUS_REG);
//}
//
//void convEncdrSetNWordsRx(u32 nwords)
//{
//	convEncdrWriteSlaveReg(CONV_ENCDR_N_DATA_WORDS_RX_REG, nwords);
//}
//
//void convEncdrSetNWordsTx(u32 nwords)
//{
//	convEncdrWriteSlaveReg(CONV_ENCDR_N_DATA_WORDS_TX_REG, nwords);
//}
//
//void convEncdrSetRxAddr(u32 addr)
//{
//	convEncdrWriteSlaveReg(CONV_ENCDR_RX_START_ADDR_REG, addr);
//}
//
//void convEncdrSetTxAddr(u32 addr)
//{
//	convEncdrWriteSlaveReg(CONV_ENCDR_TX_START_ADDR_REG, addr);
//}



