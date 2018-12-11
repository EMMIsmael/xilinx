/*
 * Freq_Intlv_Ctrlr.c *
 *  Created on: 21 Jun 2017
 *      Author: JULIAN MORTIMER
 */

#include "include.h"
//
//
//const u32 freqIntlvrAdIntlvrReadReg[FREQ_INTLVR_N_BITS_PER_SYMBOL] = {
//    FREQ_INTLV_CTRLR_AD_INTLVR_RD_REG_0,
//    FREQ_INTLV_CTRLR_AD_INTLVR_RD_REG_1,
//    FREQ_INTLV_CTRLR_AD_INTLVR_RD_REG_2,
//    FREQ_INTLV_CTRLR_AD_INTLVR_RD_REG_3,
//    FREQ_INTLV_CTRLR_AD_INTLVR_RD_REG_4,
//    FREQ_INTLV_CTRLR_AD_INTLVR_RD_REG_5
//};
//
//u32 freqIntlvCtrlrTstBuf[A4SOURCE_BRAM_SIZE];
//u32 freqIntlvCtrlrRndBuf[A4SOURCE_BRAM_SIZE];
//u32 freqIntlvCtrlrIntlvr[FREQ_INTLVR_N_BITS_PER_SYMBOL][FREQ_INTLV_CTRLR_INTLVR_SIZE];
//u32 freqIntlvCtrlrRxBuf[FREQ_INTLV_CTRLR_RX_BUF_SIZE];
//u32 freqIntlvCtrlrModulator[FREQ_INTLV_CTRLR_MOD_SIZE];
//
//u32 txBits[384];
//
//const u32 txChanMask2048[64] = { 0x0000ffff, 0xffffffff, 0xffffffff, 0xffffffff,
//		                         0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff,
//		                         0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff,
//		                         0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff,
//		                         0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff,
//		                         0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff,
//		                         0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff,
//		                         0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff,
//		                         0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff,
//		                         0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff,
//		                         0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff,
//		                         0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff,
//		                         0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff,
//		                         0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff,
//		                         0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff,
//		                         0xffffffff, 0xffffffff, 0xffffffff, 0xffff0000 };
//
//const u32 txChanMask1024[32] = { 0x0000ffff, 0xffffffff, 0xffffffff, 0xffffffff,
//		                         0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff,
//		                         0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff,
//		                         0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff,
//		                         0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff,
//		                         0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff,
//		                         0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff,
//		                         0xffffffff, 0xffffffff, 0xffffffff, 0xffff0000 };
//
//const u32 txChanMask512[16] =  { 0x0000ffff, 0xffffffff, 0xffffffff, 0xffffffff,
//		                         0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff,
//		                         0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff,
//		                         0xffffffff, 0xffffffff, 0xffffffff, 0xffff0000 };
//
//const u32 txChanMask256[8] =   { 0x0000ffff, 0xffffffff, 0xffffffff, 0xffffffff,
//                                 0xffffffff, 0xffffffff, 0xffffffff, 0xffff0000 };
//
//const u32 txChanMask128[4] =   { 0x0000ffff, 0xffffffff, 0xffffffff, 0xffff0000 };
//
//const u32 txChanMask64[2]  =   { 0xffffffff, 0xffffffff };
//
//
//
//
//void freqIntlvCtrlrInit(tFreqIntlvCtrlr *ps, SysCtrl *pSysCtrl, u32 baseAddr, u32 id, const char name[], u32 nChanSymb, u32 nBitsRxWord, u32 nBitsChan,
//		                const u32 *txChanMask, u32 *modulator, u32 intlvr[][FREQ_INTLV_CTRLR_INTLVR_SIZE]);
//void freqIntlvCtrlrDump(tFreqIntlvCtrlr * ps);
//
//void freqIntlvCtrlrReset(tFreqIntlvCtrlr *ps);
//void freqIntlvCtrlrIntAck(tFreqIntlvCtrlr *ps);
//void freqIntlvCtrlrStart(tFreqIntlvCtrlr *ps);
//void freqIntlvCtrlrSetNWordsRx(tFreqIntlvCtrlr *ps, u32 nwords);
//void freqIntlvCtrlrSetNWordsTx(tFreqIntlvCtrlr *ps, u32 nwords);
//void freqIntlvCtrlrSetRxAddr(tFreqIntlvCtrlr *ps, u32 addr);
//void freqIntlvCtrlrSetTxAddr(tFreqIntlvCtrlr *ps, u32 addr);
//void freqIntlvCtrlrSetBitsPerByte(tFreqIntlvCtrlr *ps, u32 nbits);
//void freqIntlvCtrlrSetRegs(tFreqIntlvCtrlr *ps);
//
//void freqIntlvCtrlrEnterAdmin(tFreqIntlvCtrlr *ps);
//void freqIntlvCtrlrExitAdmin(tFreqIntlvCtrlr *ps);
//
//u32 freqIntlvCtrlrAdReadIntlvr(tFreqIntlvCtrlr *ps, u32 addr, u32 bitSel);
//void freqIntlvCtrlrAdWriteIntlvr(tFreqIntlvCtrlr *ps, u32 addr, u32 bitSel, u32 data);
//void freqIntlvCtrlrAdInitIntlvr(tFreqIntlvCtrlr *ps, u32 addr, u32 count, u32 bitSel, u32 val);
//void freqIntlvCtrlrAdLoadIntlvr(tFreqIntlvCtrlr *ps);
//
//u32 freqIntlvCtrlrAdReadDataBuf(tFreqIntlvCtrlr *ps, u32 addr);
//void freqIntlvCtrlrAdWriteDataBuf(tFreqIntlvCtrlr *ps, u32 addr, u32 data);
//void freqIntlvCtrlrAdInitDataBuf(tFreqIntlvCtrlr *ps, u32 addr, u32 count, u32 val);
//void freqIntlvCtrlrAdLoadDataBuf(tFreqIntlvCtrlr *ps, u32 addr, u32 count, u32* data);
//
//u32 freqIntlvCtrlrAdReadMod(tFreqIntlvCtrlr *ps, u32 addr);
//void freqIntlvCtrlrAdWriteMod(tFreqIntlvCtrlr *ps, u32 addr, u32 data);
//void freqIntlvCtrlrAdInitMod(tFreqIntlvCtrlr *ps, u32 addr, u32 count, u32 val);
//void freqIntlvCtrlrAdLoadMod(tFreqIntlvCtrlr *ps);
//
//void freqIntlvCtrlrWriteSlaveReg(tFreqIntlvCtrlr *ps, u32 reg, u32 data);
//u32 freqIntlvCtrlrReadSlaveReg(tFreqIntlvCtrlr *ps, u32 reg);
//u32 freqIntlvCtrlrReadStatusReg(tFreqIntlvCtrlr *ps);
//
//void freqIntlvCtrlrInitTstBuf(tFreqIntlvCtrlr *ps);
//void freqIntlvCtrlrInitRndBuf(tFreqIntlvCtrlr *ps);
//void freqIntlvCtrlrInitIntlvr(tFreqIntlvCtrlr *ps);
//void freqIntlvCtrlrInitRxBuf(void);
//void freqIntlvCtrlrInitModulator(tFreqIntlvCtrlr *ps);
//
//u32 lcm(u32 a, u32 b)
//{
//	u32 l = a;
//	while ((l % b) != 0) l += a;
//	return(l);
//}
//
//u32 freqIntlvCtrlrFitRxWords(tFreqIntlvCtrlr *ps)
//{
//	u32 u = lcm(ps->nRxBitsPack, ps->nChanBits);
//	return (ps->nSymbBits / u) * (u / ps->nRxWordBits);
//}
//
//void freqIntlvCtrlrInit(tFreqIntlvCtrlr *ps, SysCtrl *pSysCtrl, u32 baseAddr, u32 id, const char name[], u32 nChanSymb, u32 nBitsChan, u32 nBitsRxWord,
//		                const u32 *txChanMask, u32 *modulator, u32 intlvr[][FREQ_INTLV_CTRLR_INTLVR_SIZE])
//{
//	u32 i = 0;
//	u32 j = 0;
//	u32 n = 0;
//	u32 m = 0;
//	strcpy(ps->name, name);
//	ps->SysCtrl      = pSysCtrl,
//	ps->baseAddr     = baseAddr;
//	ps->id           = id;
//	ps->nChanSymb    = nChanSymb;
//	ps->nRxWordBits  = nBitsRxWord;
//	ps->nChanBits    = nBitsChan;
//
//	ps->txChanMask   = txChanMask;
//	for (i = 0; i < (ps->nChanSymb >> 5); i++)
//	{
//		m = txChanMask[i];
//		for ( j = 0; j < 32; j++ )
//		{
//			if (m & 0x80000000) n++;
//			m <<= 1;
//		}
//	}
//	ps->nChanUsed    = n;
//
//	ps->nRxWordPack  = FREQ_INTLV_CTRLR_RX_BUF_WORD_WIDTH / ps->nRxWordBits;
//	ps->nRxBitsPack  = ps->nRxWordPack * ps->nRxWordBits;
//	ps->nSymbBits    = ps->nChanBits * ps->nChanUsed;
//
//	ps->nRxWordsSymb = freqIntlvCtrlrFitRxWords(ps);
//	ps->nRxBitsUsed  = ps->nRxWordsSymb * ps->nRxWordBits;
//	ps->modulator    = modulator;
//	ps->intlvr       = intlvr;
//	freqIntlvCtrlrInitModulator(ps);
//	freqIntlvCtrlrInitIntlvr(ps);
//}
//
//void freqIntlvCtrlrDump(tFreqIntlvCtrlr * ps)
//{
//	printf("\n%s ", ps->name );
//	printf("System Specification:\n");
//	printf("\n");
//    printf("Base address %08x id 0x%02x\n", (unsigned int)ps->baseAddr, (unsigned int)ps->id);
//	printf("\n");
//	printf("Number of Carriers (N) :  %u\n",   (unsigned int)ps->nChanSymb);
//	printf("Used Carriers          :  %u\n",   (unsigned int)ps->nChanUsed);
//	printf("Carrier Bits           :  %u\n",   (unsigned int)ps->nChanBits);
//	printf("Symbol Bits            :  %u\n",   (unsigned int)ps->nSymbBits);
//	printf("\n");
//	printf("Rx word Bits           :  %u\n",   (unsigned int)ps->nRxWordBits);
//	printf("Rx Symbol Words        :  %u\n",   (unsigned int)ps->nRxWordsSymb);
//	printf("Rx Symbol Bits Used    :  %u\n",   (unsigned int)ps->nRxBitsUsed);
//	printf("\n");
//	printf("Rx Buffer Words Packed :  %u\n",   (unsigned int)ps->nRxWordPack);
//	printf("Rx Buffer Bits Packed  :  %u\n",   (unsigned int)ps->nRxBitsPack);
//	printf("\n");
//}
//
//
//// Interleaver addresses packed data bits,
//// unused channels set to point to last location in TX buffer memory, which should contain zero.
//void freqIntlvCtrlrInitIntlvr(tFreqIntlvCtrlr *ps)
//{
//	u32 i = 0;
//	u32 j = 0;
//    u32 bit_addr = 0;
//    u32 m = 0;
//    u32 addr_zero = FREQ_INTLV_CTRLR_TX_BUF_SIZE-1;
//    u32 rxBitCnt = 0;
//
//    for (i = 0; i < FREQ_INTLV_CTRLR_INTLVR_SIZE; i++)
//    {
//    	for (j = 0; j < FREQ_INTLVR_N_BITS_PER_SYMBOL; j++)
//    	{
//    		ps->intlvr[j][i] = addr_zero;
//    	}
//    }
//
//	for (i = 0; i < ps->nChanSymb; i++)
//	{
//		if((i % 32) == 0) m = ps->txChanMask[i / 32];
//
//		for (j = 0; j < ps->nChanBits; j++)
//		{
//			if (m & 0x80000000)
//			{
//				if (rxBitCnt < ps->nRxBitsUsed)
//				{
//					rxBitCnt++;
//		            ps->intlvr[j][i] = bit_addr++;
//		            if ((bit_addr % 32) == ps->nRxBitsPack) bit_addr += (32 - ps->nRxBitsPack);
//				}
//		        else
//		        {
//				    ps->intlvr[j][i] = addr_zero;
//		        }
//			}
//			else
//			{
//			    ps->intlvr[j][i] = addr_zero;
//			}
//		}
//		m <<= 1;
//	}
//}
//
//void freqIntlvCtrlrInitRxBuf()
//{
//	u32 i = 0;
//	for (i = 0; i < FREQ_INTLV_CTRLR_RX_BUF_SIZE; i++)
//		freqIntlvCtrlrRxBuf[i] = i;
//}
//
//void freqIntlvCtrlrInitModulator(tFreqIntlvCtrlr *ps)
//{
//	u32 i = 0;
//	u32 a = 0;
//	u32 amask = 0;
//	const u32 *modtbl;
//	switch(ps->nChanBits)
//	{
//	  case 2:
//		  modtbl = modulators[1]->constellation;
//		  break;
//	  case 4:
//		  modtbl = modulators[2]->constellation;
//		  break;
//	  case 6:
//		  modtbl = modulators[3]->constellation;
//		  break;
//	  default:
//		  assert(FALSE);
//		  break;
//	}
//	amask = (1 << ps->nChanBits) - 1;
//	for (i = 0; i < FREQ_INTLV_CTRLR_MOD_SIZE /*(1 << ps->nChanBits)*/; i++)
//	{
//		ps->modulator[i] = modtbl[a++ & amask];
//	}
//}
//
//void freqIntlvCtrlrInitTstBuf(tFreqIntlvCtrlr *ps)
//{
//	u32 i = 0;
//	u32 j = 0;
//	u32 u = 0;
//
//	u32 bitptr = 0;
//	u32 nSymbs = (1 << ps->nChanBits);
//	char s[FREQ_INTLVR_N_BITS_PER_SYMBOL+1];
//
//	for (i = 0; i < nSymbs; i++)
//	{
//		for (j = 0; j < ps->nChanBits; j++)
//		{
//			txBits[bitptr++] = (i >> j /*(ps->nChanBits - 1 - j)*/) & 0x1;
//		}
//	}
//
//
//	bitptr = 0;
//	for (i = 0; i < A4SOURCE_BRAM_SIZE; i++)
//	{
//		u = 0;
//		for (j = 0; j < ps->nRxWordBits; j++)
//		{
//    		u = (u >> 1) + (txBits[bitptr] << (ps->nRxWordBits-1));
//			bitptr = ((bitptr < (nSymbs * ps->nChanBits - 1)) ? bitptr + 1 : 0);
//		}
//		freqIntlvCtrlrTstBuf[i] = u;
//	}
//
//	if(verboseMode)
//	{
//		printf("\n");
//		bitptr = 0;
//		for (i = 0; i < 3 * nSymbs; i++)
//		{
//			for (j = 0; j < ps->nChanBits; j++ )
//			{
//				if (((freqIntlvCtrlrTstBuf[bitptr / 32] >> (bitptr % 32)) & 0x1) == 0x1)
//					s[ps->nChanBits-1-j] = '1';
//				else
//					s[ps->nChanBits-1-j] = '0';
//
//				bitptr = (((bitptr % 32) < (ps->nRxWordBits - 1)) ? bitptr + 1 : ((bitptr / 32) << 5) + 32);
//			}
//			s[j] ='\0';
//			if ((i + 1) % 4 == 0)
//				printf("%s\n", s);
//			else
//				printf("%s ", s);
//		}
//	}
//}
//
//void freqIntlvCtrlrInitRndBuf(tFreqIntlvCtrlr *ps)
//{
//	u32 i = 0;
//	for(i = 0; i < A4SOURCE_BRAM_SIZE; i++)
//	    freqIntlvCtrlrRndBuf[i] = prng32();
//}
//
//void freqIntlvCtrlrReset(tFreqIntlvCtrlr *ps)
//{
//    hwWriteCtrlReg(ps->baseAddr, FREQ_INTLV_CTRLR_CTRL_REG << 2, FREQ_INTLV_CTRLR_SOFT_RESET);
//    while((freqIntlvCtrlrReadStatusReg(ps) & FREQ_INTLV_CTRLR_PS_RESET) != 0)
//	   ;
//}
//
//void freqIntlvCtrlrIntAck(tFreqIntlvCtrlr *ps)
//{
//    hwWriteCtrlReg(ps->baseAddr, FREQ_INTLV_CTRLR_CTRL_REG << 2, FREQ_INTLV_CTRLR_INT_ACK);
//}
//
//void freqIntlvCtrlrStart(tFreqIntlvCtrlr *ps)
//{
//    hwWriteCtrlReg(ps->baseAddr, FREQ_INTLV_CTRLR_CTRL_REG << 2, FREQ_INTLV_CTRLR_START);
//}
//
//void freqIntlvCtrlrEnterAdmin(tFreqIntlvCtrlr *ps)
//{
//    hwWriteCtrlReg(ps->baseAddr, FREQ_INTLV_CTRLR_CTRL_REG << 2, FREQ_INTLV_CTRLR_ENTER_ADMIN);
//}
//
//void freqIntlvCtrlrExitAdmin(tFreqIntlvCtrlr *ps)
//{
//    hwWriteCtrlReg(ps->baseAddr, FREQ_INTLV_CTRLR_CTRL_REG << 2, FREQ_INTLV_CTRLR_EXIT_ADMIN);
//}
//
//
//u32 freqIntlvCtrlrAdReadIntlvr(tFreqIntlvCtrlr *ps, u32 addr, u32 bitSel)
//{
//    u32 u = 0;
//    if (bitSel >= FREQ_INTLVR_N_BITS_PER_SYMBOL) bitSel = 0;
//    freqIntlvCtrlrEnterAdmin(ps);
//    u = hwReadData(ps->baseAddr,
//                   FREQ_INTLV_CTRLR_CTRL_REG << 2, FREQ_INTLV_CTRLR_AD_RD_INTLVR, FREQ_INTLV_CTRLR_AD_RD_INTLVR_ACK,
//                   FREQ_INTLV_CTRLR_STATUS_REG << 2, FREQ_INTLV_CTRLR_AD_INTLVR_RD_RDY,
//                   FREQ_INTLV_CTRLR_AD_RD_WR_ADDR_REG << 2, addr,
//                   freqIntlvrAdIntlvrReadReg[bitSel] << 2
//                   );
//    freqIntlvCtrlrExitAdmin(ps);
//    return u;
//}
//
//void freqIntlvCtrlrAdWriteIntlvr(tFreqIntlvCtrlr *ps, u32 addr, u32 bitSel, u32 data)
//{
//	freqIntlvCtrlrEnterAdmin(ps);
//    if (bitSel >= FREQ_INTLVR_N_BITS_PER_SYMBOL) bitSel = 0;
//    freqIntlvCtrlrWriteSlaveReg(ps, FREQ_INTLV_CTRLR_AD_BIT_SEL_REG, bitSel);
//
//    hwWriteData(ps->baseAddr,
//    		    FREQ_INTLV_CTRLR_CTRL_REG << 2, FREQ_INTLV_CTRLR_AD_WR_INTLVR,
//    		    FREQ_INTLV_CTRLR_AD_RD_WR_ADDR_REG << 2, addr,
//    		    FREQ_INTLV_CTRLR_AD_WR_DATA_REG << 2, data
//    		    );
//    freqIntlvCtrlrExitAdmin(ps);
//}
//
//void freqIntlvCtrlrAdInitIntlvr(tFreqIntlvCtrlr *ps, u32 addr, u32 count, u32 bitSel, u32 val)
//{
//    u32 i = 0;
//    freqIntlvCtrlrEnterAdmin(ps);
//    if (bitSel >= FREQ_INTLVR_N_BITS_PER_SYMBOL) bitSel = 0;
//    freqIntlvCtrlrWriteSlaveReg(ps, FREQ_INTLV_CTRLR_AD_BIT_SEL_REG, bitSel);
//
//    for(i = 0; i < count; i++)
//	{
//	    hwWriteData(ps->baseAddr,
//	    		    FREQ_INTLV_CTRLR_CTRL_REG << 2, FREQ_INTLV_CTRLR_AD_WR_INTLVR,
//	    		    FREQ_INTLV_CTRLR_AD_RD_WR_ADDR_REG << 2, addr + i,
//	    		    FREQ_INTLV_CTRLR_AD_WR_DATA_REG << 2, val
//	    		    );
//	}
//    freqIntlvCtrlrExitAdmin(ps);
//}
//
//void freqIntlvCtrlrAdLoadIntlvr(tFreqIntlvCtrlr *ps)
//{
//    u32 i = 0;
//    u32 j = 0;
//    freqIntlvCtrlrEnterAdmin(ps);
//
//    for(j = 0; j < FREQ_INTLVR_N_BITS_PER_SYMBOL; j++)
//    {
//    	freqIntlvCtrlrWriteSlaveReg(ps, FREQ_INTLV_CTRLR_AD_BIT_SEL_REG, j);
//    	for (i = 0; i < FREQ_INTLV_CTRLR_INTLVR_SIZE; i++)
//	    {
//	        hwWriteData(ps->baseAddr,
//	        		    FREQ_INTLV_CTRLR_CTRL_REG << 2, FREQ_INTLV_CTRLR_AD_WR_INTLVR,
//	        		    FREQ_INTLV_CTRLR_AD_RD_WR_ADDR_REG << 2, i,
//	        		    FREQ_INTLV_CTRLR_AD_WR_DATA_REG << 2, ps->intlvr[j][i]
//	        		    );
//	    }
//    }
//    freqIntlvCtrlrExitAdmin(ps);
//}
//
//u32 freqIntlvCtrlrAdReadDataBuf(tFreqIntlvCtrlr *ps, u32 addr)
//{
//    u32 u = 0;
//    freqIntlvCtrlrEnterAdmin(ps);
//    u = hwReadData(ps->baseAddr,
//                   FREQ_INTLV_CTRLR_CTRL_REG << 2, FREQ_INTLV_CTRLR_AD_RD_TX_BUF, FREQ_INTLV_CTRLR_AD_RD_TX_BUF_ACK,
//                   FREQ_INTLV_CTRLR_STATUS_REG << 2, FREQ_INTLV_CTRLR_AD_TX_BUF_RD_RDY,
//                   FREQ_INTLV_CTRLR_AD_RD_WR_ADDR_REG << 2, addr,
//                   FREQ_INTLV_CTRLR_AD_TX_BUF_RD_REG << 2
//                   );
//    freqIntlvCtrlrExitAdmin(ps);
//    return u;
//}
//
//void freqIntlvCtrlrAdWriteDataBuf(tFreqIntlvCtrlr *ps, u32 addr, u32 data)
//{
//	freqIntlvCtrlrEnterAdmin(ps);
//    hwWriteData(ps->baseAddr,
//    		    FREQ_INTLV_CTRLR_CTRL_REG << 2, FREQ_INTLV_CTRLR_AD_WR_RX_BUF,
//    		    FREQ_INTLV_CTRLR_AD_RD_WR_ADDR_REG << 2, addr,
//    		    FREQ_INTLV_CTRLR_AD_WR_DATA_REG << 2, data
//    		    );
//    freqIntlvCtrlrExitAdmin(ps);
//}
//
//void freqIntlvCtrlrAdInitDataBuf(tFreqIntlvCtrlr *ps, u32 addr, u32 count, u32 val)
//{
//    u32 i = 0;
//    freqIntlvCtrlrEnterAdmin(ps);
//	for(i = 0; i < count; i++)
//	{
//	    hwWriteData(ps->baseAddr,
//	    		    FREQ_INTLV_CTRLR_CTRL_REG << 2, FREQ_INTLV_CTRLR_AD_WR_RX_BUF,
//	    		    FREQ_INTLV_CTRLR_AD_RD_WR_ADDR_REG << 2, addr + i,
//	    		    FREQ_INTLV_CTRLR_AD_WR_DATA_REG << 2, val
//	    		    );
//	}
//    freqIntlvCtrlrExitAdmin(ps);
//}
//
//void freqIntlvCtrlrAdLoadDataBuf(tFreqIntlvCtrlr *ps, u32 addr, u32 count, u32* data)
//{
//    u32 i = 0;
//    freqIntlvCtrlrEnterAdmin(ps);
//	for(i = 0; i < count; i++)
//	{
//	    hwWriteData(ps->baseAddr,
//	    		    FREQ_INTLV_CTRLR_CTRL_REG << 2, FREQ_INTLV_CTRLR_AD_WR_RX_BUF,
//	    		    FREQ_INTLV_CTRLR_AD_RD_WR_ADDR_REG << 2, addr + i,
//	    		    FREQ_INTLV_CTRLR_AD_WR_DATA_REG << 2, data[i]
//	    		    );
//	}
//    freqIntlvCtrlrExitAdmin(ps);
//}
//
//u32 freqIntlvCtrlrAdReadMod(tFreqIntlvCtrlr *ps, u32 addr)
//{
//    u32 u = 0;
//    freqIntlvCtrlrEnterAdmin(ps);
//    u = hwReadData(ps->baseAddr,
//                   FREQ_INTLV_CTRLR_CTRL_REG << 2, FREQ_INTLV_CTRLR_AD_RD_MOD, FREQ_INTLV_CTRLR_AD_RD_MOD_ACK,
//                   FREQ_INTLV_CTRLR_STATUS_REG << 2, FREQ_INTLV_CTRLR_AD_MOD_RD_RDY,
//                   FREQ_INTLV_CTRLR_AD_RD_WR_ADDR_REG << 2, addr,
//                   FREQ_INTLV_CTRLR_AD_MOD_RD_REG << 2
//                   );
//    freqIntlvCtrlrExitAdmin(ps);
//    return u;
//}
//
//void freqIntlvCtrlrAdWriteMod(tFreqIntlvCtrlr *ps, u32 addr, u32 data)
//{
//	freqIntlvCtrlrEnterAdmin(ps);
//    hwWriteData(ps->baseAddr,
//    		    FREQ_INTLV_CTRLR_CTRL_REG << 2, FREQ_INTLV_CTRLR_AD_WR_MOD,
//    		    FREQ_INTLV_CTRLR_AD_RD_WR_ADDR_REG << 2, addr,
//    		    FREQ_INTLV_CTRLR_AD_WR_DATA_REG << 2, data
//    		    );
//    freqIntlvCtrlrExitAdmin(ps);
//}
//
//void freqIntlvCtrlrAdInitMod(tFreqIntlvCtrlr *ps, u32 addr, u32 count, u32 val)
//{
//    u32 i = 0;
//    freqIntlvCtrlrEnterAdmin(ps);
//	for(i = 0; i < count; i++)
//	{
//	    hwWriteData(ps->baseAddr,
//	    		    FREQ_INTLV_CTRLR_CTRL_REG << 2, FREQ_INTLV_CTRLR_AD_WR_MOD,
//	    		    FREQ_INTLV_CTRLR_AD_RD_WR_ADDR_REG << 2, addr + i,
//	    		    FREQ_INTLV_CTRLR_AD_WR_DATA_REG << 2, val
//	    		    );
//	}
//    freqIntlvCtrlrExitAdmin(ps);
//}
//
//void freqIntlvCtrlrAdLoadMod(tFreqIntlvCtrlr *ps)
//{
//    u32 i = 0;
//    freqIntlvCtrlrEnterAdmin(ps);
//	for(i = 0; i < FREQ_INTLV_CTRLR_MOD_SIZE; i++)
//	{
//	    hwWriteData(ps->baseAddr,
//	    		    FREQ_INTLV_CTRLR_CTRL_REG << 2, FREQ_INTLV_CTRLR_AD_WR_MOD,
//	    		    FREQ_INTLV_CTRLR_AD_RD_WR_ADDR_REG << 2, i,
//	    		    FREQ_INTLV_CTRLR_AD_WR_DATA_REG << 2, ps->modulator[i]
//	    		    );
//	}
//    freqIntlvCtrlrExitAdmin(ps);
//}
//
//void freqIntlvCtrlrWriteSlaveReg(tFreqIntlvCtrlr *ps, u32 reg, u32 data)
//{
//	hwWriteReg(ps->baseAddr, reg << 2, data);
//}
//
//u32 freqIntlvCtrlrReadSlaveReg(tFreqIntlvCtrlr *ps, u32 reg)
//{
//	return hwReadReg(ps->baseAddr, reg << 2);
//}
//
//u32 freqIntlvCtrlrReadStatusReg(tFreqIntlvCtrlr *ps)
//{
//	return freqIntlvCtrlrReadSlaveReg(ps, FREQ_INTLV_CTRLR_STATUS_REG);
//}
//
//void freqIntlvCtrlrSetNWordsRx(tFreqIntlvCtrlr *ps, u32 nwords)
//{
//	freqIntlvCtrlrWriteSlaveReg(ps, FREQ_INTLV_CTRLR_N_DATA_WORDS_RX_REG, nwords);
//}
//
//void freqIntlvCtrlrSetNWordsTx(tFreqIntlvCtrlr *ps, u32 nwords)
//{
//	freqIntlvCtrlrWriteSlaveReg(ps, FREQ_INTLV_CTRLR_N_DATA_WORDS_TX_REG, nwords);
//}
//
//void freqIntlvCtrlrSetRxAddr(tFreqIntlvCtrlr *ps, u32 addr)
//{
//	freqIntlvCtrlrWriteSlaveReg(ps, FREQ_INTLV_CTRLR_RX_START_ADDR_REG, addr);
//}
//
//void freqIntlvCtrlrSetTxAddr(tFreqIntlvCtrlr *ps, u32 addr)
//{
//	freqIntlvCtrlrWriteSlaveReg(ps, FREQ_INTLV_CTRLR_TX_START_ADDR_REG, addr);
//}
//
//void freqIntlvCtrlrSetBitsPerByte(tFreqIntlvCtrlr *ps, u32 nbits)
//{
//	freqIntlvCtrlrWriteSlaveReg(ps, FREQ_INTLV_CTRLR_RX_BITS_PER_BYTE_REG, nbits);
//}
//
//void freqIntlvCtrlrSetRegs(tFreqIntlvCtrlr *ps)
//{
//	freqIntlvCtrlrSetNWordsRx(ps, ps->nRxWordsSymb);
//	freqIntlvCtrlrSetNWordsTx(ps, ps->nChanSymb);
//	freqIntlvCtrlrSetRxAddr(ps, 0);
//	freqIntlvCtrlrSetTxAddr(ps, 0);
//	freqIntlvCtrlrSetBitsPerByte(ps, ps->nRxWordBits);
//
//}
//
