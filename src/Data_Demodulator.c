/*
 * Data_Demodulator.c
 *
 *  Created on: 14 Apr 2018
 *      Author: JULIAN MORTIMER
 */
#include "include.h"
#define DDM_BASE_ADDR             ( DATA_DEMOD_BASEADDR )

static s32      TransFunc   [ DATA_DEMOD_LUT_SIZE ];
u32             utilDmodLUT [ DATA_DEMOD_LUT_SIZE ];
u32             utilRawLUT  [ DATA_DEMOD_LUT_SIZE ];

void      DDUtilCreateDemod ( ModType modType );
void      DDUtilDumpDmodLUT ( u32 addr, u32 fmt );
void      DDUtilDumpRawLUT  ( u32 addr, u32 fmt );
DemodData DDUtilDemodulate  ( Cplx32 din );
u32       *DDUtilGetLUTPtr  ( void );

static void CalcTransFunc   ( Modem *pd, u32 bit );
static void CalcLUTBit      ( Modem *pd, u32 bit );
static void MakeLUT         ( Modem *pd );


u32 *DDUtilGetLUTPtr( void )
{
    return ( UtilModem.pLUT );
}


void DDUtilCreateDemod( ModType modType )
{
    Modem *pMod = &UtilModem;
    pMod->InitModem (
            pMod,
            modType,
            utilDmodLUT,
            NULL,
            NULL,
            CalcTransFunc,
            CalcLUTBit,
            MakeLUT
    );

    pMod->MakeDmodLUT( pMod );
}

inline u32 DDUtilRawLUT  ( u32 addr ) { return( utilRawLUT  [ addr ] ); }
inline u32 DDUtilDmodLUT ( u32 addr ) { return( utilDmodLUT [ addr ] ); }

void DDUtilDumpDmodLUT( u32 addr, u32 fmt )
{
    Modem *pMod = &UtilModem;
    u32 arg1 = 0, arg2 = 0;
      if ( fmt == 0 )
    {
        fmt = DUMP_FMT_HEX;
        arg1 = 4;
        arg2 = 0;
    }
    else
    {
        fmt  = DUMP_FMT_SOFT_BITS;
        arg1 = pMod->SoftBits;
        arg2 = pMod->DataBits;
    }
    utilScreenDumpData( DUMP_SEL_DDM_UTIL_DMOD_LUT, fmt, DUMP_TYP_MEMORY, addr, arg1, arg2, 0, "Data Demodulator Util Dmod LUT:" );
}

void DDUtilDumpRawLUT( u32 addr, u32 fmt )
{
    Modem *pMod = &UtilModem;
    u32 arg1 = 0, arg2 = 0;
      if ( fmt == 0 )
    {
        fmt = DUMP_FMT_HEX;
        arg1 = 4;
        arg2 = 0;
    }
    else
    {
        fmt  = DUMP_FMT_SOFT_BITS;
        arg1 = pMod->SoftBits;
        arg2 = pMod->DataBits;
    }
    utilScreenDumpData( DUMP_SEL_DDM_UTIL_RAW_LUT, fmt, DUMP_TYP_MEMORY, addr, arg1, arg2, 0, "Data Demodulator Util Raw LUT:" );
}

DemodData DDUtilDemodulate( Cplx32 din )
{
    Modem    *pd      = &UtilModem;
    u32      *pLUT    = pd->pLUT;
    u32      wordMask = ( 1 << pd->DataBits ) - 1;
    DemodData data = { 0, 0, 0, 0, 0, 0, 0, 0 };

    u32 softBits = pd->SoftBits;
    u32 bitMask  = ( 1 << softBits ) - 1;
    u32 signBit  = 1 << ( softBits - 1 );
    u32 ibits = pLUT [ din.real & wordMask ];
    u32 qbits = pLUT [ wordMask - ( din.imag & wordMask ) ];

    data.raw = Cplx32toU32( din );

    data.bit0  = ( ( ( ibits >> ( 0 * softBits ) ) & bitMask ) ^ signBit ) - signBit;
    data.bit1  = ( ( ( ibits >> ( 1 * softBits ) ) & bitMask ) ^ signBit ) - signBit;
    data.bit2  = ( ( ( ibits >> ( 2 * softBits ) ) & bitMask ) ^ signBit ) - signBit;
    data.bit3  = ( ( ( qbits >> ( 0 * softBits ) ) & bitMask ) ^ signBit ) - signBit;
    data.bit4  = ( ( ( qbits >> ( 1 * softBits ) ) & bitMask ) ^ signBit ) - signBit;
    data.bit5  = ( ( ( qbits >> ( 2 * softBits ) ) & bitMask ) ^ signBit ) - signBit;

    data.word |= ( data.bit0 < 0 ? 1 << 0 : 0 );
    data.word |= ( data.bit1 < 0 ? 1 << 1 : 0 );
    data.word |= ( data.bit2 < 0 ? 1 << 2 : 0 );
    data.word |= ( data.bit3 < 0 ? 1 << 3 : 0 );
    data.word |= ( data.bit4 < 0 ? 1 << 4 : 0 );
    data.word |= ( data.bit5 < 0 ? 1 << 5 : 0 );

    return ( data );
}


// Construct the function which gives the soft codes for bits between points on the
// constellation where they are 1 and points where they are zero.
// The Xilinx Viterbi decoder requires codes between
// 0 and 2^( softBits - 1 ) - 1, with the msb denoting the value of the bit
// where n is the number of soft decision bits.
//
// At the moment, we are using cos( n/2 * pi/n )^0.5, negated
// and shifted to give a transition between +1 and -1
//
const u32 softMapOffsetBinary [ 16 ] = {
     0x0,  // Strongest 0
     0x1,
     0x2,
     0x3,
     0x4,
     0x5,
     0x6,
     0x7,  // Weakest 0
     0x8,  // Weakest 1
     0x9,
     0xA,
     0xB,
     0xC,
     0xD,
     0xE,
     0xF   // Strongest 1
};

const u32 softMapSignedMagnitude [ 16 ] = {
     0x7,  // Strongest 0
     0x6,
     0x5,
     0x4,
     0x3,
     0x2,
     0x1,
     0x0,  // Weakest 0
     0x8,  // Weakest 1
     0x9,
     0xA,
     0xB,
     0xC,
     0xD,
     0xE,
     0xF   // Strongest 1
};

const u32 *softBitMap = softMapOffsetBinary;
const double   TFGrad = 4;

// Make a transfer function to offset binary spec:
// Both low-to-high ( offset 0 )
// and  high-to-low ( offset npts )
// In the case of 4 soft bits:
// 0  = strongest zero
// 7  = weakest zero
// 8  = weakest  1
// 15 = strongest 1
static void CalcTransFunc( Modem *pm, u32 bit )
{
    const DmodTFDesc *pd = pm [ bit ].pTFDesc;
    u32 tfl        = pd->TFLength;
    s32 ofs        = 1 << ( pm->SoftBits - 1 );  // 8
    s32 max        = ofs - 1;                    // 7

	for ( u32 i = 0; i < DATA_DEMOD_LUT_SIZE; i++ ) TransFunc [ i ] = 0;

	for ( u32 i = tfl / 2; i < tfl; i++ )
	{
	    double k  = ( double ) i;
	    double d  = ( round( tanh( k * TFGrad / tfl ) ) * ( double ) max );
		s32 TF = ( s32 ) d + ofs;
		TransFunc [ i ] = TF;
	}
    for ( u32 i = 0; i < tfl / 2; i++ )
    {
        TransFunc [ i ] = max - ( TransFunc [ tfl - 1 - i ] - ofs );
    }
    for ( u32 i = 0; i < tfl; i++ )
    {
        TransFunc [ tfl + i ] = TransFunc [ tfl - 1 - i ];
    }
}


// Use a table to indicate when bit value transitions between its two possible states.
// During the transition, generate a function to input to the Viterbi soft decision decoder
//
static void CalcLUTBit( Modem *pm, u32 bit )
{
    const DmodTFDesc *pd = pm->pTFDesc;
    u32  b           = pd [ bit ].InitVal;
    u32  lenTF       = pd [ bit ].TFLength;
    const u32 *pTbl  = pd [ bit ].TFTable;
    u32 inTF         = false;                      // true when using transition function table to transition between bit values
    s32 max          = ( 1 << pm->SoftBits ) - 1;  // 15
    s32 softBit      = ( b ? max : 0 );
    u32 iTF          = 0;
    u32 cTF          = 0;
    u32 k            = 1 << ( pm->AddrBits - 1 );

    for( u32 i = 0; i < ( 1 << pm->AddrBits ); i++ )
    {
        if ( inTF )
        {
            softBit = TransFunc [ iTF++ ];
            if( cTF >= ( lenTF - 1 ) )
            {
                inTF = false;
            }
            else
            {
                cTF++;
            }
        }
        else if ( i == *pTbl )
        {
            pTbl++;
            iTF     = ( ( b == 1 ) ? lenTF : 0 );
            cTF     = 0;
            b       = b ^ 1;

            softBit = TransFunc [ iTF++ ];
            cTF++;
            inTF      = true;
        }
        else
        {
            softBit = ( ( b == 1 ) ? max : 0 );
        }

        utilRawLUT [ k ] |= softBit << ( bit * pm->SoftBits );
        if ( k >= ( ( 1 << pm->AddrBits ) - 1 ) )
        {
            k = 0;
        }
        else
        {
            k++;
        }
    }
}

static void MakeLUT( Modem *pm )
{
    for( u32 i = 0; i < ( 1 << pm->AddrBits ); i++ )
    {
        pm->pLUT [ i ] = 0;
        utilRawLUT [ i ] = 0;
        utilDmodLUT [ i ] = 0;
    }
    for( u32 i = 0; i < ( pm->SymbBits / 2 ); i++ )
    {
        pm->CalcDmodTF     ( pm, i );
        pm->CalcDmodLUTBit ( pm, i );
    }
    for( u32 i = 0; i < ( 1 << pm->AddrBits ); i++ )
    {
        pm->pLUT [ i ] = softBitMap [ utilRawLUT [ i ] ];
    }
}


//static DemodLUT DataDemodLUT;
//static s32      dataDemodTransFunc [ DATA_DEMOD_LUT_SIZE / 2 + 1 ];
//static u32      dataDemodLUT [ DATA_DEMOD_LUT_SIZE ];
//static u32      dataDemodLUTRdBuf [ DATA_DEMOD_LUT_SIZE ];
//static u32      dataDemodCountLUT [ DATA_DEMOD_LUT_SIZE ];

//static u32      dataDemodId;
//_Bool DataDemodInit( ModType modType, SoftBits softBits, u32 nFFT, u32 nSymbolsFrame, u32 nFramesBlock );
//void DataDemodResetDDM( void );
//u32  DataDemodRdLocalLUTData( u32 addr );
//void DataDemodInitHWLUT( u32 addr, u32 count, u32 data );
//void DataDemodLoadHWLUT( void );
//void DataDemodLoadHWLUTCount( void );
//void DataDemodWrHWLUTData( u32 addr, u32 data );
//u32  DataDemodRdHWLUTData( u32 addr );
//void DataDemodWrSlvReg( u32 reg, u32 data );
//u32  DataDemodRdSlvReg( u32 reg );
//_Bool DataDemodRdHWLUTBuf( void );
//u32  DataDemodRdHWLUTBufData( u32 addr );
//u32  DataDemodPRNG( void );
//void DataDemodCalcLUT( void );
//void DataDemodCalcTstVec( u32 sf );
//void DataDemodLoadTstVec( void );
//void DataDemodDumpSlvRegs( u32 dbin );
//void DataDemodDumpLocalLUT( u32 addr, u32 fmt );
//void DataDemodDumpHWLUT( u32 addr, u32 fmt );
//u32  DataDemodRdTstVec( u32 addr );
//void DataDemodDumpTstVec( u32 addr, u32 nbits );
//u32  DataDemodRdTFData( u32 addr );
//void DataDemodDumpTF( u32 addr, u32 fmt );
//void DataDemodSetNFFT( u32 n );
//void DataDemodSetNSymbolsFrame( u32 n );
//void DataDemodSetNFramesBlock( u32 n );
//void DataDemodTest( ModType modType, SoftBits softBits, u32 nFFT, u32 nSymbolsFrame, u32 nFramesBlock, u32 nBlocks );
//_Bool DataDemodSelfTest( _Bool verbose );
//void DataDemodDump( void );
//u32  DataDemodModType( void );
//void DataDemodStart( void );
//void DataDemodStop( void );
//static void DDMReset( void );
////static void DDMIntAck( void );
//static void DDMStart( void );
//static void DDMStop( void );
//static void DDMEnterAdmin( void );
//static void DDMExitAdmin( void );
//static void DDMWriteSlaveReg( u32 reg, u32 data );
//static u32  DDMReadSlaveReg( u32 reg );
//static u32  DDMReadStatusReg( void );
//static u32  DDMAdReadLUT( u32 addr );
//static void DDMAdWriteLUT( u32 addr, u32 data );
//static void DDMAdInitLUT( u32 addr, u32 count, u32 val );
//static void DDMAdLoadLUT( u32 addr, u32 count, u32* data );
//static _Bool DDMAdReadLUTBuf( u32 *pData, u32 nWords );
//static void DDMSetNFFT( u32 n );
//static void DDMSetNSymbolsFrame( u32 n );
//static void DDMSetNFramesBlock( u32 n );
//static u32 DDMNFFT( void );
//static u32 DDMNSymbolsFrame( void );
//static u32 DDMNFramesBlock( void );
//static u32 DDMDeviceId( void );
//static u32 DDMNSlvReg( void );
//static u32 DDMNWordsRx( void );
//static u32 DDMNSymbolsRx( void );
//static u32 DDMNFramesRx( void );
//static u32 DDMNBlocksRx( void );
//static u32 DDMNWordsTx( void );
//static u32 DDMNSymbolsTx( void );
//static u32 DDMNFramesTx( void );
//static u32 DDMNBlocksTx( void );

//_Bool DataDemodInit( ModType modType, SoftBits softBits, u32 nFFT, u32 nSymbolsFrame, u32 nFramesBlock )
//{
//	DataDemodLUT = DEMODLUT_TEMPLATE;
////	DemodAPI *pAPI = DemodulatorAPI;
//	DemodLUT *pLUT = &DataDemodLUT;
//    Modulator *pMod = GetModulator( modType );
//	pLUT->InitDemodLUT(
//		pLUT,
//		DMOD_MODE_DATA,
//		modType,
//		DATA_DEMOD_LUT_SIZE,
//		softBits,
//        dataDemodLUT,
//        dataDemodTransFunc,
//        pMod->demod_ttbl,
//        NBitsPerSymbol,
//        LUTBits,
//        CalcTransFunc,
//        CalcLUTBit,
//	    CalcLUT
//	  );
////    dataDemodId = pAPI->InitDemod(
////    	pAPI,
////    	DATA_DEMOD_AXI_BASEADDR,
////    	ID_NEW_DEVICE,
////    	"Rx Data Demodulator",
////    	DATA_DEMOD_N_FFT,
////    	DATA_DEMOD_N_SYMBOLS_FRAME,
////    	DATA_DEMOD_N_FRAMES_BLOCK,
////    	DATA_DEMOD_DEVICE_ID
////    );
////	  DataDemodResetDDM( );
////      DataDemodSetNFFT( nFFT );
////      DataDemodSetNSymbolsFrame( nSymbolsFrame );
////      DataDemodSetNFramesBlock( nFramesBlock );
////      DataDemodLoadHWLUT( );
////    DataDemodLoadTstVec( );
////      return DataDemodSelfTest( true );
//      return ( true );
//}
//

//void DataDemodResetDDM( void )
//{
//    DDMReset( );
//}
//
//u32 DataDemodModType( void )
//{
//	return DataDemodLUT.modType;
//}

//void DataDemodCalcLUT( void )
//{
//	DemodLUT *pLUT = &DataDemodLUT;
//	pLUT->CalcDemodLUT( pLUT );
//}
//
//void DataDemodLoadHWLUTCount( void )
//{
//	for( int i = 0; i < DATA_DEMOD_LUT_SIZE; i++ )
//	{
//		dataDemodCountLUT [ i ] = i;
//	}
//	DDMAdLoadLUT( 0, DATA_DEMOD_LUT_SIZE, dataDemodCountLUT );
//}
//
//
//void DataDemodInitHWLUT( u32 addr, u32 count, u32 data )
//{
//	DDMAdInitLUT( addr, count, data );
//}
//
//void DataDemodLoadHWLUT( void )
//{
//	DemodLUT *pLUT = &DataDemodLUT;
//	DDMAdLoadLUT( 0                 , pLUT->LUTSize >> 1, pLUT->pLUT + ( pLUT->LUTSize >> 1 ) );
//	DDMAdLoadLUT( pLUT->LUTSize >> 1, pLUT->LUTSize >> 1, pLUT->pLUT                        );
//}

//u32 DataDemodRdLocalLUTData( u32 addr )
//{
//	DemodLUT *pLUT = &DataDemodLUT;
//	return pLUT->pLUT [ addr ];
//}

//void DataDemodWrHWLUTData( u32 addr, u32 data )
//{
////	DemodAPI *pAPI = DemodulatorAPI;
////    pAPI->WriteDemodLUTData( pAPI, dataDemodId, addr, data );
//	DDMAdWriteLUT( addr, data );
//
//}

//u32 DataDemodRdHWLUTData( u32 addr )
//{
////	DemodAPI *pAPI = DemodulatorAPI;
////    return pAPI->ReadDemodLUTData( pAPI, dataDemodId, addr );;
//	return ( DDMAdReadLUT( addr ) );
//}

//void DataDemodWrSlvReg( u32 reg, u32 data )
//{
////	DemodAPI *pAPI = DemodulatorAPI;
////    pAPI->WriteDemodSlvReg( pAPI, dataDemodId, reg, data );
//    DDMWriteSlaveReg( reg, data );
//}
//
//u32 DataDemodRdSlvReg( u32 reg )
//{
////	DemodAPI *pAPI = DemodulatorAPI;
////    return pAPI->ReadDemodSlvReg( pAPI, dataDemodId, reg );
//	return ( DDMReadSlaveReg( reg ) );
//}
//
//_Bool DataDemodRdHWLUTBuf( void )
//{
////	DemodAPI *pAPI = DemodulatorAPI;
////	pAPI->ReadDemodLUT( pAPI, dataDemodId, 0, DATA_DEMOD_LUT_SIZE, dataDemodLUTRdBuf );
//	return ( DDMAdReadLUTBuf( dataDemodLUTRdBuf, DATA_DEMOD_LUT_SIZE ) );
//}
//
//u32 DataDemodRdHWLUTBufData( u32 addr )
//{
//	return dataDemodLUTRdBuf [ addr ];
//}

//u32 DataDemodPRNG( void )
//{
//	static PRNG32 prng = { 0x12345678, 0x12345678, prng32s };
//	return prng.prng( &prng, false );
//}

//void DataDemodCalcTstVec( u32 sf )
//{
////    DemodLUT *pLUT =  &DataDemodLUT;
////    ModType modType = pLUT->modType;
////    Modulator *mod  = modulators [ modType ];
//    u32 nwds = DATA_DEMOD_TEST_VEC_SIZE;
//    u32 nfft = DDMNFFT( );
////    for( int i = 0; i < nwds; i++ )
////	{
////
////        u32 u  = DataDemodPRNG( );
////        dataDemodTstVec [ i ] = mod->Modulate( mod, u, sf );
////		dataDemodTstVec [ i ] =  0;
////	}
////	for ( int i = 1; i < nwds; i++ )
////	{
////		if( ( i % nfft ) == 0 )
////		{
////			dataDemodTstVec [ i - 128 ] = 32767;
////		    dataDemodTstVec [ i + 128 ] = 32767;
////			dataDemodTstVec [ i - 64 ] = 32767;
////		    dataDemodTstVec [ i + 64 ] = 32767;
////			dataDemodTstVec [ i - 32 ] = 32767;
////		    dataDemodTstVec [ i + 32 ] = 32767;
////		}
////	}
////	dataDemodTstVec [ 128 ]        = 32767;
////	dataDemodTstVec [ nwds - 128 ] = 32767;
////	dataDemodTstVec [ 64 ]         = 32767;
////	dataDemodTstVec [ nwds - 64 ] = 32767;
////	dataDemodTstVec [ 32 ]         = 32767;
////	dataDemodTstVec [ nwds - 32 ] = 32767;
//    srcCalcPRS( nfft, dataDemodTstVec, nwds );
//
////	for( s32 i = 0; i < nWords / 8; i ++ )
////	{
////		s32 s = ( s32 ) DataDemodPRNG( ) & 0x0000001F;
////		s = ( s - 16 ) & 0xfff;
////		dataDemodTstVec [ i ] += s + ( s << 16 );
////	}
//}
//
//void DataDemodLoadTstVec( void )
//{
//	u32 nWords = DATA_DEMOD_TEST_VEC_SIZE;
//	DataDemodCalcTstVec( DATA_MOD_SCALE );
//	double f = SYS_TIMER( ReadTimerLongSec );
//	gmAxiBramController0InitBramDDR3( 0, nWords, ( u32 * )dataDemodTstVec );
//	f = SYS_TIMER( ReadTimerLongSec ) - f;
//	printf( " Load time for %d bytes = %.8fs, %.0f bytes/sec\n", ( unsigned int ) ( 4 * nWords ), f, 4.0 * nWords / f );
//}
//
//void DataDemodDumpSlvRegs( u32 dbin )
//{
//	u32 dfmt = ( dbin ? DUMP_FMT_BIN32 : DUMP_FMT_HEX_DEC );
//    utilScreenDumpData( 0, dfmt, DUMP_TYP_SLV_REG, DATA_DEMOD_BASEADDR, DATA_DEMOD_N_SLV_REG, 0, 0, "Data Demodulator Slave Registers:" );
//}

//void DataDemodDumpLocalLUT( u32 addr, u32 fmt )
//{
//	DemodLUT *pLUT = &DataDemodLUT;
//	u32 arg1 = 0, arg2 = 0;
//      if ( fmt == 0 )
//    {
//    	fmt = DUMP_FMT_HEX;
//        arg1 = 4;
//        arg2 = 0;
//    }
//    else
//    {
//    	fmt  = DUMP_FMT_SOFT_BITS;
//        arg1 = pLUT->softBits;
//        arg2 = pLUT->LUTSize;
//    }
//	utilScreenDumpData( DUMP_SEL_DDM_LOCAL_LUT, fmt, DUMP_TYP_MEMORY, addr, arg1, arg2, 0, "Data Demodulator Local LUT:" );
//}




//void DataDemodDumpHWLUT( u32 addr, u32 fmt )
//{
//	u32 arg1 = 0, arg2 = 0;
//    if ( !DataDemodRdHWLUTBuf( ) )
//    {
//        printf ( "\nTimed out reading Data Demodulator HW LUT\n" ); fflush( stdout );
//    }
//    else
//    {
//        if ( fmt == 0 )
//        {
//        	fmt = DUMP_FMT_HEX;
//            arg1 = 4;
//            arg2 = 0;
//        }
//        else
//        {
//        	DemodLUT *pLUT = &DataDemodLUT;
//        	fmt  = DUMP_FMT_SOFT_BITS_HW;
//            arg1 = pLUT->softBits;
//            arg2 = pLUT->LUTSize;
//        }
//	    utilScreenDumpData( DUMP_SEL_DDM_HW_LUT, fmt, DUMP_TYP_MEMORY, addr, arg1, arg2, 0, "Data Demodulator HW LUT:" );
//    }
//}
//
//u32  DataDemodRdTstVec( u32 addr )
//{
//	return( dataDemodTstVec [ addr ] );
//}
//
//void DataDemodDumpTstVec( u32 addr, u32 nbits )
//{
//	u32 arg1 = 0, fmt = 0;
//	if ( nbits == 0 )
//	{
//		fmt = DUMP_FMT_HEX;
//        arg1 = 4;
//	}
//	else
//	{
//		fmt = DUMP_FMT_IQ_SIGNED;
//		arg1 = nbits;
//	}
//	utilScreenDumpData( DUMP_SEL_DDM_TST_VEC, fmt, DUMP_TYP_MEMORY, addr, arg1, 0, 0, "Data Demodulator Test Vectors:" );
//}
//
//u32  DataDemodRdTFData( u32 addr )
//{
//	DemodLUT *pLUT = &DataDemodLUT;
//	return ( *( u32* )&pLUT->pTransFunc [ addr ] );
//}
//
//void DataDemodDumpTF( u32 addr, u32 fmt )
//{
//   	if ( fmt == 0 )
//	{
//		fmt = DUMP_FMT_HEX;
//	}
//	else if ( fmt == 1 )
//	{
//		fmt = DUMP_FMT_DEC_SIGNED;
//	}
//	else
//	{
//		fmt = DUMP_FMT_GRAPH;
//
//	}
//
//	utilScreenDumpData( DUMP_SEL_DDM_TF, fmt, DUMP_TYP_MEMORY, addr, 3, 56, 4, "Data Demodulator Transition Function:" );
//}
//
//void DataDemodSetNFFT( u32 n )
//{
//	DDMSetNFFT( n );
//}
//
//void DataDemodSetNSymbolsFrame( u32 n )
//{
//	DDMSetNSymbolsFrame( n );
//}
//
//void DataDemodSetNFramesBlock( u32 n )
//{
//	DDMSetNFramesBlock( n );
//}
//
//void DataDemodStart( void )
//{
//	DDMStart( );
//}
//
//void DataDemodStop( void )
//{
//	DDMStop( );
//}
//
//void DataDemodTest( ModType modType, SoftBits softBits, u32 nFFT, u32 nSymbolsFrame, u32 nFramesBlock, u32 nBlocks )
//{
//	DataDemodInit( modType, softBits, nFFT, nSymbolsFrame, nFramesBlock );
//	axi4StreamDataSourceInit( ( u32* )dataDemodTstVec,  DATA_DEMOD_TEST_VEC_SIZE, 0, 0, nFFT, nSymbolsFrame, nFramesBlock );
//	axi4StreamDataSinkInit(  &AXI4Sink0, "AXI4 Stream Data Sink 0", A4SINK_0_BASE_ADDR, 0, nFFT * nSymbolsFrame, 0, nFFT, nSymbolsFrame, nFramesBlock );
//
//	u32             nWordsFrame = nFFT         * nSymbolsFrame;
//	u32             nWordsBlock = nFramesBlock * nWordsFrame;
//	axi4StreamError ecode;
//	double          t = SysTimer->ReadTimerLongSec( SysTimer, SysTimerId );
//    int             i = 0;
//
//	DDMStart( );
//	for ( i = 0; i < nBlocks; i ++ )
//	{
//	  ecode = axi4StreamTfrNFrames( nFramesBlock, nWordsFrame, nWordsFrame );
//	  if ( ecode == AXI4_TIMEOUT ) break;
//    }
//	DDMStop( );
//	if ( ecode == AXI4_TIMEOUT )
//	{
//		printf( "\nBlock transfer timed out after %d block%s.\n", i, SNMBR( i ) );
//	}
//	else
//	{
//		t = SysTimer->ReadTimerLongSec( SysTimer, SysTimerId ) - t;
//		printf( "\nTime to transfer %d blocks of %d frames of %d words: %0.9fs.  Data rate = %0.0f Mbytes/sec\n",
//				( unsigned int ) nBlocks, ( unsigned int ) nFramesBlock, ( unsigned int ) nWordsFrame, t, ( double ) 4.0 * nBlocks * nWordsBlock / ( 1.0e6 * t ) );
//	}
//
//}
//
//_Bool DataDemodSelfTest( _Bool verbose )
//{
//	DemodLUT *pLUT = &DataDemodLUT;
//	u32 deviceId = DDMDeviceId( );
//	if ( deviceId != DATA_DEMOD_DEVICE_ID )
//    {
//        if ( verbose )
//        {
//            printf ( "\nData Demodulator Device Id Incorrect: is %06X, should be %06X\n", ( unsigned int )deviceId, DATA_DEMOD_DEVICE_ID ); fflush( stdout );
//        }
//        return false;
//    }
//    DataDemodLoadHWLUT( );
//    if ( !DataDemodRdHWLUTBuf( ) )
//    {
//        if ( verbose )
//        {
//            printf ( "\nTimed out reading Data Demodulator HW LUT\n" ); fflush( stdout );
//        }
//        return false;
//    }
//    for ( u32 i = 0; i < pLUT->LUTSize / 2; i++ )
//    {
//    	u32 idxl = ( pLUT->LUTSize / 2 + i );
//    	u32 idxh = ( pLUT->LUTSize / 2 - i );
//    	u32 u    = dataDemodLUTRdBuf [ i ];
//    	u32 v    = pLUT->pLUT [ idxl ] + ( u32 )( pLUT->pLUT [ idxh ] << 16 );
//    	if ( v != u )
//        {
//            if ( verbose )
//            {
//                printf ( "\nData Demodulator HW LUT data error at address %04X, reads %08X, should be %08X\n", ( unsigned int ) i, ( unsigned int ) v, ( unsigned int ) u ); fflush( stdout );
//            }
//            return false;
//        }
//    }
//    if ( verbose )
//    {
//        printf ( "\nData Demodulator self test passed\n" ); fflush( stdout );
//    }
//    return true;
//}
//
//void DataDemodDump( void )
//{
//    printf( "\n\n" ); fflush( stdout );
//	printf( "Data Demodulator:\n" ); fflush( stdout );
//	printf( "=================\n" ); fflush( stdout );
//    printf( "\n" ); fflush( stdout );
//    printf( "Device Id        : %08X\n", ( unsigned int )DDMDeviceId( ) ); fflush( stdout );
//    printf( "slave registers  : %d\n"  , ( unsigned int )DDMNSlvReg( ) ); fflush( stdout );
//    printf( "\n" ); fflush( stdout );
//    printf( "NFFT             : %d\n"  , ( unsigned int )DDMNFFT( ) ); fflush( stdout );
//    printf( "symbols / frame  : %d\n"  , ( unsigned int )DDMNSymbolsFrame( ) ); fflush( stdout );
//    printf( "frames / block   : %d\n"  , ( unsigned int )DDMNFramesBlock( ) ); fflush( stdout );
//    printf( "\n" ); fflush( stdout );
//    printf( "Rx Words         : %d\n"  , ( unsigned int )DDMNWordsRx( ) ); fflush( stdout );
//    printf( "Rx Symbols       : %d\n"  , ( unsigned int )DDMNSymbolsRx( ) ); fflush( stdout );
//    printf( "Rx Frames        : %d\n"  , ( unsigned int )DDMNFramesRx( ) ); fflush( stdout );
//    printf( "Rx Blocks        : %d\n"  , ( unsigned int )DDMNBlocksRx( ) ); fflush( stdout );
//    printf( "Tx Words         : %d\n"  , ( unsigned int )DDMNWordsTx( ) ); fflush( stdout );
//    printf( "Tx Symbols       : %d\n"  , ( unsigned int )DDMNSymbolsTx( ) ); fflush( stdout );
//    printf( "Tx Frames        : %d\n"  , ( unsigned int )DDMNFramesTx( ) ); fflush( stdout );
//    printf( "Tx Blocks        : %d\n"  , ( unsigned int )DDMNBlocksTx( ) ); fflush( stdout );
//    printf( "\n" ); fflush( stdout );
//    DataDemodSelfTest( true );
//}
//
//
////====================================================================
////Temporary code
////====================================================================
//static void DDMReset( void )
//{
//    hwWriteCtrlReg( DDM_BASE_ADDR, DDM_CTRL_REG << 2, DDM_SOFT_RESET );
//    while( ( DDMReadStatusReg( ) & DDM_PS_RESET ) != 0 )
//	   ;
//}
//
////static void DDMIntAck( void )
////{
////    hwWriteCtrlReg( DDM_BASE_ADDR, DDM_CTRL_REG << 2, DDM_INT_ACK );
////}
//
//static void DDMStart( void )
//{
//    hwWriteCtrlReg( DDM_BASE_ADDR, DDM_CTRL_REG << 2, DDM_START );
//}
//
//static void DDMStop( void )
//{
//    hwWriteCtrlReg( DDM_BASE_ADDR, DDM_CTRL_REG << 2, DDM_STOP );
//}
//
//static void DDMEnterAdmin( void )
//{
//    hwWriteCtrlReg( DDM_BASE_ADDR, DDM_CTRL_REG << 2, DDM_ENTER_ADMIN );
//}
//
//static void DDMExitAdmin( void )
//{
//    hwWriteCtrlReg( DDM_BASE_ADDR, DDM_CTRL_REG << 2, DDM_EXIT_ADMIN );
//}
//
//static void DDMAdWriteLUT( u32 addr, u32 data )
//{
//	DDMEnterAdmin( );
//    hwWriteData( DDM_BASE_ADDR,
//    		    DDM_CTRL_REG << 2, DDM_AD_WR_LUT,
//    		    DDM_AD_RD_WR_ADDR_REG << 2, addr,
//    		    DDM_AD_WR_DATA_REG << 2, data
//    		    );
//    DDMExitAdmin( );
//}
//
//
//u32 DDMAdReadLUT( u32 addr )
//{
//    u32 u = 0;
//    DDMEnterAdmin( );
//    u = hwReadData( DDM_BASE_ADDR,
//                   DDM_CTRL_REG << 2, DDM_AD_RD_LUT, DDM_AD_RD_LUT_ACK,
//                   DDM_STATUS_REG << 2, DDM_AD_LUT_RD_RDY,
//                   DDM_AD_RD_WR_ADDR_REG << 2, addr,
//                   DDM_AD_RD_DATA_REG << 2
//                   );
//    DDMExitAdmin( );
//    return u;
//}
//
//_Bool DDMAdReadLUTBuf( u32 *pData, u32 nWords )
//{
//    _Bool res = true;
//    DDMEnterAdmin( );
//    for( int i = 0; i < nWords; i++ )
//    {
//        double f = SysTimer->ReadTimerLongSec( SysTimer, SysTimerId );
//    	pData [ i ] = hwReadData( DDM_BASE_ADDR,
//                   DDM_CTRL_REG << 2, DDM_AD_RD_LUT, DDM_AD_RD_LUT_ACK,
//                   DDM_STATUS_REG << 2, DDM_AD_LUT_RD_RDY,
//                   DDM_AD_RD_WR_ADDR_REG << 2, ( u32 )  i,
//                   DDM_AD_RD_DATA_REG << 2
//                   );
//    	if ( SysTimer->ReadTimerLongSec( SysTimer, SysTimerId ) - f > HW_READ_BRAM_TIMEOUT ) { res = false; break; }
//    }
//    DDMExitAdmin( );
//    return( res );
//}
//
//
//static void DDMAdInitLUT( u32 addr, u32 count, u32 val )
//{
//    u32 i = 0;
//    DDMEnterAdmin( );
//	for( i = 0; i < count; i++ )
//	{
//	    hwWriteData( DDM_BASE_ADDR,
//	    		    DDM_CTRL_REG << 2, DDM_AD_WR_LUT,
//	    		    DDM_AD_RD_WR_ADDR_REG << 2, addr + i,
//	    		    DDM_AD_WR_DATA_REG << 2, val
//	    		    );
//	}
//    DDMExitAdmin( );
//}
//
//static void DDMAdLoadLUT( u32 addr, u32 count, u32* data )
//{
//    u32 i = 0;
//    DDMEnterAdmin( );
//	for( i = 0; i < count; i++ )
//	{
//	    hwWriteData( DDM_BASE_ADDR,
//	    		    DDM_CTRL_REG << 2, DDM_AD_WR_LUT,
//	    		    DDM_AD_RD_WR_ADDR_REG << 2, addr + i,
//	    		    DDM_AD_WR_DATA_REG << 2, data [ i ]
//	    		    );
//	}
//    DDMExitAdmin( );
//}
//
//
//static void DDMWriteSlaveReg( u32 reg, u32 data )
//{
//	hwWriteReg( DDM_BASE_ADDR, reg << 2, data );
//}
//
//u32 DDMReadSlaveReg( u32 reg )
//{
//	return hwReadReg( DDM_BASE_ADDR, reg << 2 );
//}
//
//u32 DDMReadStatusReg( void )
//{
//	return DDMReadSlaveReg( DDM_STATUS_REG );
//}
//
//static void DDMSetNFFT( u32 n )
//{
//	DDMWriteSlaveReg( DDM_N_FFT_REG, n );
//}
//
//static void DDMSetNSymbolsFrame( u32 n )
//{
//	DDMWriteSlaveReg( DDM_N_SYMBOLS_FRAME_REG, n );
//}
//
//static void DDMSetNFramesBlock( u32 n )
//{
//	DDMWriteSlaveReg( DDM_N_FRAMES_BLOCK_REG, n );
//}
//
//static u32 DDMNFFT( void )
//{
//	return DDMReadSlaveReg( DDM_N_FFT_REG );
//}
//
//static u32 DDMNSymbolsFrame( void )
//{
//	return DDMReadSlaveReg( DDM_N_SYMBOLS_FRAME_REG );
//}
//
//static u32 DDMNFramesBlock( void )
//{
//	return DDMReadSlaveReg( DDM_N_FRAMES_BLOCK_REG );
//}
//
//static u32 DDMDeviceId( void )
//{
//	return DDMReadSlaveReg( DDM_INFO_REG );
//}
//
//static u32 DDMNSlvReg( void )
//{
//	return DDMReadSlaveReg( DDM_INFO_REG ) & 0x000000ff;
//}
//
//static u32 DDMNWordsRx( void )
//{
//	return DDMReadSlaveReg( DDM_RX_DATA_CTR_REG );
//}
//
//static u32 DDMNSymbolsRx( void )
//{
//	return DDMReadSlaveReg( DDM_RX_SYMBOL_CTR_REG );
//}
//
//static u32 DDMNFramesRx( void )
//{
//	return DDMReadSlaveReg( DDM_RX_FRAME_CTR_REG );
//}
//
//static u32 DDMNBlocksRx( void )
//{
//	return DDMReadSlaveReg( DDM_RX_BLOCK_CTR_REG );
//}
//
//static u32 DDMNWordsTx( void )
//{
//	return DDMReadSlaveReg( DDM_TX_DATA_CTR_REG );
//}
//
//static u32 DDMNSymbolsTx( void )
//{
//	return DDMReadSlaveReg( DDM_TX_SYMBOL_CTR_REG );
//}
//
//static u32 DDMNFramesTx( void )
//{
//	return DDMReadSlaveReg( DDM_TX_FRAME_CTR_REG );
//}
//
//static u32 DDMNBlocksTx( void )
//{
//	return DDMReadSlaveReg( DDM_TX_BLOCK_CTR_REG );
//}









