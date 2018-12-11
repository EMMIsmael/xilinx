/*
 * cplx32.c
 *
 *  Created on: 16 Aug 2018
 *      Author: JULIAN MORTIMER
 */
#include "include.h"

_Bool Cplx32Eq( Cplx32 z1, Cplx32 z2 );
Cplx32 ScaleCplx32( Cplx32 z, double sf );
Cplx32 CplxtoCplx32( cplx z );
cplx Cplx32toCplx( Cplx32 z );
Cplx32 U32toCplx32( u32 u );
u32 Cplx32toU32( Cplx32 z );
u32 U32Sext( u32 z, u32 nbits );

u32 U32Sext( u32 z, u32 nbits )
{

    u64 u = ( ~( ( u64 ) z & ( 0x00010001LL << ( nbits - 1 ) ) ) + 0x00010001LL ) & ~0x00010001LL;
    return ( z | ( u32 ) u );
}

// Scale a 32 bit complex number stored as imaginary:real by a power of 2.
// If sf is not a power of two, it will be rounded up to nearest.
_Bool Cplx32Eq( Cplx32 z1, Cplx32 z2 )
{
    return( ( ( u32 )( z1.real ) == ( u32 ) ( z2.real ) ) && ( ( u32 )( z1.imag ) == ( u32 ) ( z2.imag ) ) );
}

// Scale a 32 bit complex number stored as imaginary:real by a power of 2.
// If sf is not a power of two, it will be rounded up to nearest.
Cplx32 ScaleCplx32( Cplx32 z, double sf )
{
    cplx zz = Cplx32toCplx( z );
    zz = ( sf + 0.0 * I ) * zz;
    return ( CplxtoCplx32( zz ) );
}


// Convert cplx to cplx32,
// Binary point 0 - scaling done by caller
// range -32768:32767
Cplx32 CplxtoCplx32( cplx z )
{
    Cplx32 zz;
    zz.real = ( s16 ) creal( z );
    zz.imag = ( s16 ) cimag( z );
    return ( zz );
}

// Convert complex-encoded u32 to cplx
cplx Cplx32toCplx( Cplx32 z )
{
    return ( ( double ) ( z.real ) + I * ( double ) ( z.imag ) );
}


Cplx32 U32toCplx32( u32 u )
{
    u16 r = ( u16 ) ( u & 0xFFFF );
    u16 i = ( u16 ) ( ( u >> 16 ) & 0xFFFF );
    s16 real = *( s16 * ) &r;
    s16 imag = *( s16 * ) &i;
    Cplx32 z = { real, imag };
    return ( z );
}

u32 Cplx32toU32( Cplx32 z )
{
    u32 real = ( u32 ) z.real;
    u32 imag = ( u32 ) z.imag;
    return ( ( imag << 16 ) | ( real & 0xFFFF ) );
}













