/*
 * AXI4_Stream_Comms.c
 *
 *  Created on: 15 Jul 2017
 *      Author: JULIAN MORTIMER
 */
#include "include.h"

axi4StreamError axi4StreamTfrNFrames(u32 nFrames, u32 nFrmWordsTx, u32 nFrmWordsRx);
axi4StreamError axi4StreamTfrBlock(u32 txAddr, u32 rxAddr, u32 nwordsTx, u32 nwordsRx);

axi4StreamError axi4StreamTxNFrames(u32 nFrames, u32 nFrmWordsTx);


//axi4StreamError axi4StreamTfrBlock(u32 txAddr, u32 rxAddr, u32 nwordsTx, u32 nwordsRx)
//{
//    u32 statusRegTx = 0;
//    u32 statusRegRx = 0;
//    double f = 0.0;
//	axi4StreamDataSourceTxStream(txAddr, nwordsTx);
//	axi4StreamDataSinkRxStream( &AXI4Sink0, rxAddr, nwordsRx);
//    f = SysTimer->ReadTimerLongSec(SysTimer, SysTimerId);
//	while (((statusRegRx & A4SINK_RX_STREAM_DONE) == 0 ) || ((statusRegTx & A4SOURCE_TX_STREAM_DONE) == 0))
//    {
//        statusRegRx = axi4StreamDataSinkReadSlave( &AXI4Sink0, A4SINK_STAT_REG);
//	    statusRegTx = axi4StreamDataSourceReadSlave(A4SOURCE_STAT_REG);
//		if ((SysTimer->ReadTimerLongSec(SysTimer, SysTimerId) - f) > AXI4_TIMEOUT_SEC)
//		{
//            return AXI4_TIMEOUT;
//		}
//    }
//	return AXI4_SUCCESS;
//}
//
//axi4StreamError axi4StreamTfrNFrames(u32 nFrames, u32 nFrmWordsTx, u32 nFrmWordsRx)
//{
//	u32 nfTx = 0;
//	u32 nfRx = 0;
//	float f = SysTimer->ReadTimerLongSec(SysTimer, SysTimerId);
//    float g = 0.0;
//
//	axi4StreamDataSourceTxStream(0, nFrmWordsTx); while((axi4StreamDataSourceReadSlave(A4SOURCE_STAT_REG) & A4SOURCE_TX_STREAM_DONE) == A4SOURCE_TX_STREAM_DONE);
//	axi4StreamDataSinkStart( &AXI4Sink0 );
//	while ((nfTx < nFrames) || (nfRx < nFrames))
//	{
//		if (nfTx < nFrames)
//		{
//			if ((axi4StreamDataSourceReadSlave(A4SOURCE_STAT_REG) & A4SOURCE_TX_STREAM_DONE) == A4SOURCE_TX_STREAM_DONE )
//			{
//				if(++nfTx < nFrames)
//				{
//					axi4StreamDataSourceTxStream(0, nFrmWordsTx);
//					while ((axi4StreamDataSourceReadSlave(A4SOURCE_STAT_REG) & A4SOURCE_TX_STREAM_DONE) == A4SOURCE_TX_STREAM_DONE )
//						;
//				}
//			}
//		}
//		nfRx = axi4StreamDataSinkReadSlave( &AXI4Sink0, A4SINK_RX_FRAME_CTR_REG );
//		g = SysTimer->ReadTimerLongSec(SysTimer, SysTimerId);
//		if ((g - f) > AXI4_TIMEOUT_SEC)
//		{
//			axi4StreamDataSinkStop( &AXI4Sink0 );
//            return AXI4_TIMEOUT;
//		}
//	}
//	axi4StreamDataSinkStop( &AXI4Sink0 );
//	return AXI4_SUCCESS;
//}
//
//axi4StreamError axi4StreamTxNFrames(u32 nFrames, u32 nFrmWordsTx)
//{
//	u32 nfTx = 0;
//	double f = SysTimer->ReadTimerLong(SysTimer, SysTimerId);
//
//	axi4StreamDataSourceTxStream(0, nFrmWordsTx); while((axi4StreamDataSourceReadSlave(A4SOURCE_STAT_REG) & A4SOURCE_TX_STREAM_DONE) == A4SOURCE_TX_STREAM_DONE);
//	while (nfTx < nFrames)
//	{
//		if ((axi4StreamDataSourceReadSlave(A4SOURCE_STAT_REG) & A4SOURCE_TX_STREAM_DONE) == A4SOURCE_TX_STREAM_DONE )
//		{
//			if(++nfTx < nFrames)
//			{
//				axi4StreamDataSourceTxStream(0, nFrmWordsTx);
//				while ((axi4StreamDataSourceReadSlave(A4SOURCE_STAT_REG) & A4SOURCE_TX_STREAM_DONE) == A4SOURCE_TX_STREAM_DONE )
//					;
//			}
//		}
//		if ((SysTimer->ReadTimerLong(SysTimer, SysTimerId) - f) > AXI4_TIMEOUT_TICKS)
//		{
//            return AXI4_TIMEOUT;
//		}
//	}
//	return AXI4_SUCCESS;
//}
