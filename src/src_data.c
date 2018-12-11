/*
 * src_data.c

 *
 *  Created on: 29 Nov 2017
 *      Author: JULIAN MORTIMER
 */

#include "include.h"


void srcCalcRandomSymbols(u32 nfft, Modem *pmod, u32 buf[], u32 bufSize);
void srcCalcTone(u32 nfft, u32 toneIdx, u32 ii, u32 qq, u32 buf[], u32 bufSize);
void srcCalcPRS(u32 nfft, u32 buf[], u32 bufSize);







void srcCalcRandomSymbols(u32 nfft, Modem *pmod, u32 buf[], u32 bufSize)
{
	u32 i = 0;
    u32 k = 0;

	for( i = 0; i < bufSize; i++ )
	{
		if ((k = i % nfft) && ((k <= ((nfft/2) - (nfft/8))) || (k >= ((nfft/2) + (nfft/8)))))
		{
          buf[i] = Cplx32toU32( pmod->Modulate( pmod, prng32(), DATA_MOD_SCALE ) );
		}
		else
		{
		  buf[i] = 0;
		}
	}
}

void srcCalcTone(u32 nfft, u32 toneIdx, u32 ii, u32 qq, u32 buf[], u32 bufSize)
{
	u32 i = 0;
	u32 k = 0;

	for( i = 0; i < bufSize; i++ )
	{
        if( k == toneIdx )
		  buf[ i ] = ii + ( qq << 16 );
        else
          buf[ i ] = 0;

        if( ++k == nfft )
        	k = 0;
	}
}

void srcCalcPRS(u32 nfft, u32 buf[], u32 bufSize)
{
	u32 i = 0, m = 0;
	u32 *psrc;

	switch(nfft)
	{
	    case 2048:
            psrc = PRS_2048;
	    	break;
	    case 1024:
            psrc = PRS_1024;
	    	break;
	    case 512:
            psrc = PRS_512;
	    	break;
	    case 256:
            psrc = PRS_256;
	    	break;
	    case 128:
            psrc = PRS_128;
	    	break;
	    case 64:
            psrc = PRS_64;
	    	break;
	    default:
            psrc = PRS_1024;
	    	break;
	}
	for(i=0; i < bufSize; i++)
	{
		buf[i] = Cplx32toU32( ScaleCplx32( U32toCplx32( psrc[ m++ ] ), 0.0625 ));
		if ( m == nfft )
			m = 0;
	}
}


