/*
 * Test_TxRx.c
 *
 *  Created on: 2 May 2018
 *      Author: JULIAN MORTIMER
 */
#include "include.h"









//void testTxRxNFrames( ModType modType, u32 softBits, u32 nFFT, u32 nSymbolsFrame, u32 nFramesBlock, u32 nBlocks, u32 sf );
//
//
//void testTxRxNFrames( ModType modType, u32 softBits, u32 nFFT, u32 nSymbolsFrame, u32 nFramesBlock, u32 nBlocks, u32 sf )
//{
//
//	float f = SysTimer->ReadTimerLongSec( SysTimer, SysTimerId );
//    float g = 0.0;
//
//    u32 nWordsTx  = nBlocks * nFramesBlock * nSymbolsFrame * nFFT;
//    u32 nFramesTx = nBlocks * nFramesBlock;
//
//    printf("\nTx Rx Test:\n"); fflush( stdout );
//
//    ( void ) TxRxEmulInit( nFFT, nSymbolsFrame, nFramesBlock );
//
//    ( void ) DataDemodInit( modType, softBits, nFFT, nSymbolsFrame, nFramesBlock );
//    DataDemodCalcTstVec( sf );
//    axi4StreamDataSourceInit( ( u32* )dataDemodTstVec, DATA_DEMOD_TEST_VEC_SIZE, 0, nFFT * nSymbolsFrame, nFFT, nSymbolsFrame, nFramesBlock );
//
//    ( void ) fftCtrlrInit( &txFFTCtrlr, TX_FFT_CTRLR_BASE_ADDR, FFT_CTRLR_DEVICE_ID, "Tx FFT Controller",
//   		      nFFT, TX_FFT_CTRLR_NFFT_MIN, TX_FFT_CTRLR_NFFT_MAX,
//   		      FFT_CTRLR_INV, 0,
//   		      FFT_CTRLR_HAS_CP, FFT_CTRLR_IS_SCALED );
//    ( void ) fftCtrlrInit( &rxFFTCtrlr, RX_FFT_CTRLR_BASE_ADDR, FFT_CTRLR_DEVICE_ID, "Rx FFT Controller",
//    		  nFFT, RX_FFT_CTRLR_NFFT_MIN, RX_FFT_CTRLR_NFFT_MAX,
//    		  FFT_CTRLR_FWD, 0,
//    		  FFT_CTRLR_NO_CP,  FFT_CTRLR_IS_SCALED );
//
//    ( void ) tsInjInit( &txTsInj, TSI_0_BASE_ADDR, 0, "Tx Training Symbol Injector", 1024, 0, 1, 0, 1, 1, txTsInjTsBuf );
//	tsInjAdInitRxBuf  ( &txTsInj, 0, TSI_BUF_SIZE, 0xDEADBEEF );
//	tsInjAdInitRxBuf  ( &txTsInj, TSI_BUF_SIZE, TSI_BUF_SIZE, 0xBEEFDEAD );
//
//    ( void ) axi4StreamDataSinkInit( &AXI4Sink0, "AXI4 Stream Data Sink 0", A4SINK_0_BASE_ADDR, 0xDEADBEEF, 0, 0, nFFT, nSymbolsFrame, nFramesBlock );
//    ( void ) axi4StreamDataSinkInit( &AXI4Sink1, "AXI4 Stream Data Sink 1", A4SINK_1_BASE_ADDR, 0xDEADBEEF, 0, 0, nFFT, nSymbolsFrame, nFramesBlock );
//
//    f = SysTimer->ReadTimerLongSec( SysTimer, SysTimerId ) - f;
//    printf("\nHardware Initialisation time: %0.9fs\n", f); fflush( stdout );
//
//    printf("Transferring %d blocks of %d frames of %d Symbols, of %d words. Total %d words.\n",
//              (unsigned int)nBlocks,
//              (unsigned int)nFramesBlock,
//              (unsigned int)nSymbolsFrame,
//              (unsigned int)nFFT,
//              (unsigned int)nWordsTx
//              ); fflush( stdout );
//
//    axi4StreamDataSourceStartTx( );
//	tsInjStart( &txTsInj );
//	DataDemodStart( );
//	axi4StreamDataSinkWriteInhibit( &AXI4Sink0 );
//	axi4StreamDataSinkWriteInhibit( &AXI4Sink1 );
//    axi4StreamDataSinkStart( &AXI4Sink0 );
//	axi4StreamDataSinkStart( &AXI4Sink1 );
//
//    u32 nf = 0;
//    f = SysTimer->ReadTimerLongSec( SysTimer, SysTimerId );
//	for ( nf = 0; nf < nFramesTx; nf++  )
//	{
//		g = SysTimer->ReadTimerLongSec( SysTimer, SysTimerId );
//		axi4StreamDataSourceTxStream( 0,  nFFT * nSymbolsFrame  );
//	    while ((( axi4StreamDataSourceReadSlave( A4SOURCE_STAT_REG ) & A4SOURCE_TX_STREAM_DONE ) == A4SOURCE_TX_STREAM_DONE )
//	        &&  ( SysTimer->ReadTimerLongSec( SysTimer, SysTimerId ) - g < 1.0 ))
//	    	;
//	    while ((( axi4StreamDataSourceReadSlave( A4SOURCE_STAT_REG ) & A4SOURCE_TX_STREAM_DONE ) != A4SOURCE_TX_STREAM_DONE )
//	        &&  ( SysTimer->ReadTimerLongSec( SysTimer, SysTimerId ) - g < 1.0 ))
//	    	;
//
//	    if ( nf == nFramesTx / 4 )
//	    {
//	    	axi4StreamDataSinkWriteEnable( &AXI4Sink0 );
//	    	axi4StreamDataSinkWriteEnable( &AXI4Sink1 );
//	    }
//	    else if ( nf == 3 * nFramesTx / 4 )
//	    {
//	    	axi4StreamDataSinkWriteInhibit( &AXI4Sink0 );
//	    	axi4StreamDataSinkWriteInhibit( &AXI4Sink1 );
//	    }
//
//    	if ( SysTimer->ReadTimerLongSec( SysTimer, SysTimerId ) - g >= 1.0 )
//    	       { printf( "\nTimed out waiting Tx from AXI4 Source after %d tx frames\n", ( unsigned int ) nf ); fflush( stdout ); break; }
//	}
//
//
//	g = SysTimer->ReadTimerLongSec( SysTimer, SysTimerId );
//	while((( nf = TxRxEmulSwitchRxFrameCtr() ) < nFramesTx ) && ( SysTimer->ReadTimerLongSec( SysTimer, SysTimerId ) - g < 1.0 ))
//        ;
//	if ( SysTimer->ReadTimerLongSec( SysTimer, SysTimerId ) - g >= 1.0 )
//	       { printf( "\nTimed out waiting Rx at TxRx Emulator Switch after %d rx frames\n", ( unsigned int ) nf ); fflush( stdout ); }
//
//	tsInjStop( &txTsInj );
//	DataDemodStop( );
//	axi4StreamDataSinkStop( &AXI4Sink0 );
//	axi4StreamDataSinkStop( &AXI4Sink1 );
//	axi4StreamDataSinkWriteEnable( &AXI4Sink0 );
//	axi4StreamDataSinkWriteEnable( &AXI4Sink1 );
//
//	f = SysTimer->ReadTimerLongSec( SysTimer, SysTimerId ) - f;
//    printf("\n%d words transferred in %0.9fs, data rate %0.0f words/s\n", (unsigned int) nWordsTx, f, nWordsTx / f );
//}
//
