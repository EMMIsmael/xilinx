/*
 * Data_Demodulator.c
 *
 *  Created on: 14 Apr 2018
 *      Author: JULIAN MORTIMER
 */
#include "include.h"
#define DDM_BASE_ADDR             ( DATA_DEMOD_BASEADDR )

static s32      TransFuncf   [ DATA_DEMOD_LUT_SIZE ];
static s32      TransFuncg   [ DATA_DEMOD_LUT_SIZE ];
static s32      TransFunch   [ DATA_DEMOD_LUT_SIZE ];
u32             utilDmodLUTf [ DATA_DEMOD_LUT_SIZE ];
u32             utilDmodLUTg [ DATA_DEMOD_LUT_SIZE ];
u32             utilDmodLUTh [ DATA_DEMOD_LUT_SIZE ];
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
    return ( UtilModem.pLUTf );
}


void DDUtilCreateDemod( ModType modType )
{
    Modem *pMod = &UtilModem;
    pMod->InitModem (
            pMod,
            modType,
            utilDmodLUTf,
            utilDmodLUTg,
            utilDmodLUTh,
            NULL,
            NULL,
            CalcTransFunc,
            CalcLUTBit,
            MakeLUT
    );

    pMod->MakeDmodLUT( pMod );
}

inline u32 DDUtilRawLUT  ( u32 addr ) { return( utilRawLUT  [ addr ] ); }
inline u32 DDUtilDmodLUT ( u32 addr ) { return( utilDmodLUTf [ addr ] ); }

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
    u32      *pLUTf    = pd->pLUTf;
    u32      *pLUTg    = pd->pLUTg;
    u32      *pLUTh    = pd->pLUTh;
    u32      wordMask = ( 1 << pd->DataBits ) - 1;
    DemodData data = { 0, 0, 0, 0, 0, 0, 0, 0 };
    u32 tfl         = DATA_DEMOD_LUT_SIZE;
    u32 tlfOff      = tfl/2;
    u32 imag = din.imag + tlfOff;
    u32 real = din.real + tlfOff;

    u32 softBits = pd->SoftBits;

    if (pd->ModType == MOD_4QAM) {
        data.bit0 = (pLUTf[ imag] >> (softBits  - 1));
        data.bit1 = (pLUTf[ real] >> (softBits  - 1));
        data.bit2 = 0;
        data.bit3 = 0;
        data.bit4 = 0;
        data.bit5 = 0;
    }

    if (pd->ModType == MOD_16QAM) {
        data.bit0 = (pLUTg[ imag ] >> (softBits  - 1));
        data.bit1 = (pLUTf[ imag ] >> (softBits  - 1));
        data.bit2 = (pLUTg[ real ] >> (softBits  - 1));
        data.bit3 = (pLUTf[ real ] >> (softBits  - 1));
        data.bit4 = 0;
        data.bit5 = 0;
    }

    if (pd->ModType == MOD_64QAM) {
        data.bit0 = (pLUTh[ imag ] >> (softBits  - 1));
        data.bit1 = (pLUTg[ imag ] >> (softBits  - 1));
        data.bit2 = (pLUTf[ imag ] >> (softBits  - 1));
        data.bit3 = (pLUTh[ real ] >> (softBits  - 1));
        data.bit4 = (pLUTg[ real ] >> (softBits  - 1));
        data.bit5 = (pLUTf[ real ] >> (softBits  - 1));
    }
    data.raw = Cplx32toU32( din );


    data.word |= ( data.bit0 );
    data.word |= ( data.bit1 << 1 );
    data.word |= ( data.bit2 << 2 );
    data.word |= ( data.bit3 << 3 );
    data.word |= ( data.bit4 << 4 );
    data.word |= ( data.bit5 << 5 );

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
    //u32 tfl        = pd->TFLength;
    u32 tfl         = DATA_DEMOD_LUT_SIZE;
    s32 ofs        = 1 << ( pm->SoftBits - 1 );  // 8
    s32 max        = ofs - 1;                    // 7

    for ( u32 i = 0; i < DATA_DEMOD_LUT_SIZE; i++ ) {
	    TransFuncf [ i ] = 0;
	    TransFuncg [ i ] = 0;
	    TransFunch [ i ] = 0;
	}

/*	for ( u32 i = tfl / 2; i < tfl; i++ )
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
    }*/
    for ( u32 i = 0; i< tfl; i++ ) {
        if ( i < tfl/2 ) {
            TransFuncf [ i ] = 0xf;
            if ( i < tfl/4 ){
                TransFuncg [ i ] = 0xf;
                if ( i < tfl/8 ) {
                    TransFunch [ i ] = 0x0;
                } else {
                    TransFunch [ i ] = 0xf;
                }
            } else {
                TransFuncg [ i ] = 0x0;
                if ( i < (3*tfl/8) ) {
                    TransFunch [ i ] = 0xf;
                } else {
                    TransFunch [ i ] = 0x0;
                }
            }
        } else {
            TransFuncf [ i ] = 0x0;
            if ( i < ( 3*tfl/4 ) ){
                TransFuncg [ i ] = 0x0;
                if ( i < ( 5*tfl/8 ) ) {
                    TransFunch [ i ] = 0x0;
                } else {
                    TransFunch [ i ] = 0xf;
                }
            } else {
                TransFuncg [ i ] = 0xf;
                if ( i < (7*tfl/8) ) {
                    TransFunch [ i ] = 0xf;
                } else {
                    TransFunch [ i ] = 0x0;
                }
            }
        }

    }
/*    for ( u32 i = tfl / 2; i < tfl; i++ )
    {
        TransFuncf [ i ] = 0x0;
    }
    for ( u32 i = 0; i < tfl / 2; i++ )
    {
        TransFuncf [ i ] = 0xf;
    }*/
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
            softBit = TransFuncf [ iTF++ ];
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

            softBit = TransFuncf [ iTF++ ];
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
        pm->pLUTf [ i ] = 0;
        pm->pLUTg [ i ] = 0;
        pm->pLUTh [ i ] = 0;
        utilRawLUT [ i ] = 0;
        utilDmodLUTf [ i ] = 0;
    }
    for( u32 i = 0; i < ( pm->SymbBits / 2 ); i++ )
    {
        pm->CalcDmodTF     ( pm, i );
        //pm->CalcDmodLUTBit ( pm, i );
    }
    for( u32 i = 0; i < ( 1 << pm->AddrBits ); i++ )
    {

        //pm->pLUT [ i ] = softBitMap [ utilRawLUT [ i ] ];
        pm->pLUTf [ i ]  = TransFuncf [ i ];
        pm->pLUTg [ i ]  = TransFuncg [ i ];
        pm->pLUTh [ i ]  = TransFunch [ i ];

    }
}
