/*
 * gm_util.c
 *
 *  Created on: 22 Mar 2016
 *      Author: RTL3
 */
#include "include.h"

static cplx fft_out[ MAX_FFT_SIZE ];

void fft  (cplx dest[], cplx src[], int n );
void ifft (cplx dest[], cplx src[], int n );
inline u32 u32bit( u32 u, u32 i ) { return( u & ( 1 << i ) ? 1 : 0 ); }
u32 prng32s( PRNG32 *p, u32 load );
u32 prng32( void );
u64 prng64s( PRNG64 *p, u32 load );
u64 prng64( void );
u32 clog2u32( u32 n );
u32 log2u32( u32 u );
u32  UtilTestBram   ( u32 ( * Read )( u32 addr ), void ( * Write )( u32 addr, u32 data ), u32 Size, u32 BitMask, char *Name );
u32  UtilTestBram64 ( u64 ( * Read )( u32 addr ), void ( * Write )( u32 addr, u64 data ), u32 Size, u64 BitMask, char *Name );
void WaitSec( double tsec );
u32  RMeantoRNoise( double rmean, u32 nBits );

void *memsetd( void *p, double v, u32 n );
void *memsetc( void *p, cplx   v, u32 n );

double MeanXd    ( double *px, u32 n );
double SigmaXd   ( double *px, u32 n );
double SigmaX2d  ( double *px, u32 n );
double SigmaXYd  ( double *px, double *py, u32 n );
double MeanXs    ( s32    *px, u32 n );
double SigmaXs   ( s32    *px, u32 n );
double SigmaX2s  ( s32    *px, u32 n );
double SigmaXYs  ( s32    *px, s32 *py, u32 n );
double MeanXu    ( u32    *px, u32 n );
double SigmaXu   ( u32    *px, u32 n );
double SigmaX2u  ( u32    *px, u32 n );
double SigmaXYu  ( u32    *px, u32 *py, u32 n );

void WaitSec( double tsec )
{
    double t = SysTimer->ReadTimerLongSec( SysTimer, SysTimerId );
    while ( SysTimer->ReadTimerLongSec( SysTimer, SysTimerId ) - t < tsec );
}

//============================================
// Numerical section
//============================================

// Generate Random Numbers of the range of 32-bit signed integers.
u32 prng32( void )
{
	static PRNG32 rng = { 0x12345678, 0x12345678, prng32s };
	PRNG32 *p = &rng;
	return( p->prng( p, false ));
}

u32 prng32s( PRNG32 *p, u32 load )
{
	if ( load ) p->state = p->seed;
	u32 x = p->state;
	x ^= x << 13;
	x ^= x >> 17;
	x ^= x << 5;
	p->state = x;
	return( x );
}

u64 prng64( void )
{
    static PRNG64 rng = { 0x1234567812345678, 0x1234567812345678, prng64s };
    PRNG64 *p = &rng;
    return( p->prng( p, false ));
}

u64 prng64s( PRNG64 *p, u32 load )
{
    if ( load ) p->state = p->seed;
    u64 x = p->state;
    x ^= x >> 12;
    x ^= x << 25;
    x ^= x >> 27;
    p->state = x;
    return( x );
}


u32 log2u32( u32 n )
{
	int i;
	for( i = -1; n > 0; n >>= 1, i++ );
	return (( u32 ) i );
}

u32 clog2u32( u32 n )
{
	int i;
	for ( i = 0; n > 0; n >>= 1, i++ );
	return (( u32 ) i );
}


//static void _fft( cplx buf[], cplx out[], int n, int step )
//{
//	if ( step < n )
//	{
//		_fft( out, buf, n, step * 2 );
//		_fft( out + step, buf + step, n, step * 2 );
//
//		for ( int i = 0; i < n; i += 2 * step ) {
//			cplx t = cexp( -I * M_PI * i / n ) * out[ i + step ];
//			buf[ i / 2 ]       = out[ i ] + t;
//			buf[ ( i + n )/2 ] = out[ i ] - t;
//		}
//	}
//}

//// FFT - recursive.
//void fft( cplx dest[], cplx src[], int n )
//{
//	for ( int i = 0; i < n; i++ ) {
//        dest [ i ] = src [ i ] / ( double ) n;
//	}
//    memcpy( fft_out, dest, n * sizeof( cplx ) );
//
//	_fft( dest, fft_out, n, 1 );
//}
//
//// Inverse FFT
//void ifft( cplx dest[], cplx src[], int n )
//{
//    for ( int i = 0; i < n; i++ ) {
//        dest [ i ] = conj( src [ i ] );
//    }
//    memcpy( fft_out, dest, n * sizeof( cplx ) );
//
//    _fft( dest, fft_out, n, 1 );
//    for ( int i = 0; i < n; i++ ) {
//        dest [ i ] = conj ( dest [ i ] );
//    }
//}



static void __fft( cplx buf[], cplx out[], int n, int step )
{
      if ( step < n )
      {
            __fft( out, buf, n, step * 2 );
            __fft( out + step, buf + step, n, step * 2 );

            for ( int i = 0; i < n; i += 2 * step ) {
                  cplx t = cexp( -I * M_PI * i / n ) * out[ i + step ];
                  buf[ i / 2 ]     = out[ i ] + t;
                  buf[ ( i + n )/2 ] = out[ i ] - t;
            }
      }
}

static void _fft( cplx buf[], int n )
{
      int i;
      for ( i = 0; i < n; i++ ) fft_out[ i ] = buf[ i ];

      __fft( buf, fft_out, n, 1 );
}

void fft( cplx dest[], cplx src[], int n )
{
    for ( int i = 0; i < n; i++ ) {
        dest[ i ] = src[ i ];
    }
    _fft( dest, n );
    for( int i = 0; i < n; i++ ) {
        dest[ i ] = dest[ i ] / ( double ) n;
    }
}

void ifft( cplx dest[], cplx src[], int n )
{
    for ( int i = 0; i < n; i++ ) {
        dest[ i ] = conj( src[ i ] );
    }
    _fft( dest, n );
    for( int i = 0; i < n; i++ ) {
        dest[ i ] = conj( dest[ i ] );
    }
}

//============================================
// Memory Utils
//============================================
void *memsetd( void *p, double v, u32 n )
{
    double *pp = ( double* ) p;
    for ( u32 i = 0; i < n; i++ ) {
        pp [ i ] = v;
    }
    return( p );
}

void *memsetc( void *p, cplx   v, u32 n )
{
    cplx *pp = ( cplx* ) p;
    for ( u32 i = 0; i < n; i++ ) {
        pp [ i ] = v;
    }
    return( p );
}


//============================================
// Stats
//============================================
double MeanXd( double *px, u32 n ) {
    double s = SigmaXd( px, n );
    return ( s / ( double ) n );
}

double SigmaXd( double *px, u32 n ) {
    double s = 0.0;
    for( u32 i = 0; i < n; i++ ) {
        s += px [ i ];
    }
    return( s );
}
double SigmaX2d( double *px, u32 n ) {
    double s = 0.0;
    for( u32 i = 0; i < n; i++ ) {
        s += px [ i ] * px [ i ];
    }
    return( s );
}

double SigmaXYd( double *px, double *py, u32 n ) {
    double s = 0.0;
    for( u32 i = 0; i < n; i++ ) {
        s += px [ i ] * py [ i ];
    }
    return( s );
}

double MeanXs( s32 *px, u32 n ) {
    double s = SigmaXs( px, n );
    return ( s / ( double ) n );
}

double SigmaXs( s32 *px, u32 n ) {
    double s = 0.0;
    for( u32 i = 0; i < n; i++ ) {
        s += ( double ) px [ i ];
    }
    return( s );
}
double SigmaX2s( s32 *px, u32 n ) {
    double s = 0.0;
    for( u32 i = 0; i < n; i++ ) {
        s += ( double ) px [ i ] * ( double ) px [ i ];
    }
    return( s );
}

double SigmaXYs( s32 *px, s32 *py, u32 n ) {
    double s = 0.0;
    for( u32 i = 0; i < n; i++ ) {
        s += ( double ) px [ i ] * ( double ) py [ i ];
    }
    return( s );
}

double MeanXu( u32 *px, u32 n ) {
    double s = SigmaXu( px, n );
    return ( s / ( double ) n );
}

double SigmaXu( u32 *px, u32 n ) {
    double s = 0.0;
    for( u32 i = 0; i < n; i++ ) {
        s += ( double ) px [ i ];
    }
    return( s );
}
double SigmaX2u( u32 *px, u32 n ) {
    double s = 0.0;
    for( u32 i = 0; i < n; i++ ) {
        s += ( double ) px [ i ] * ( double ) px [ i ];
    }
    return( s );
}

double SigmaXYu( u32 *px, u32 *py, u32 n ) {
    double s = 0.0;
    for( u32 i = 0; i < n; i++ ) {
        s += ( double ) px [ i ] * ( double ) py [ i ];
    }
    return( s );
}


//============================================
// PL Device Section
//============================================
u32 UtilTestBram( u32 ( * Read )( u32 addr ), void ( * Write )( u32 addr, u32 data ), u32 Size, u32 BitMask, char *Name )
{
   PRNG32 Prng = { 0x12345678, 0, prng32s };
   PRNG32 *p = &Prng;
   ( void ) p->prng( p, true ); // seed prng
   for ( int addr = 0; addr < Size; addr++ ) Write( addr, p->prng( p, false ) & BitMask );
   ( void ) p->prng( p, true ); // re-seed prng
   for ( int addr = 0; addr < Size; addr++ )
   {
	   u32 r = p->prng( p, false ) & BitMask;;
	   u32 v = Read( addr ) & BitMask;
       if ( v != r )
       {
    	   if( verboseMode ) { printf( "%s BRAM test failed at address %08X: expected: %08X read: %08X\n", Name, toUint( addr ), toUint( r ), toUint( v ) ); fflush( stdout ); }
    	   return( addr );
       }
   }
   if( verboseMode ) printf("\n%s BRAM test passed", Name ); fflush( stdout );
   return( 0xFFFFFFFF );
}

u32 UtilTestBram64( u64 ( * Read )( u32 addr ), void ( * Write )( u32 addr, u64 data ), u32 Size, u64 BitMask, char *Name )
{
   PRNG64 Prng = { 0x12345678, 0, prng64s };
   PRNG64 *p = &Prng;
   ( void ) p->prng( p, true ); // seed prng
   for ( int addr = 0; addr < Size; addr++ ) Write( addr, p->prng( p, false ) & BitMask );
   ( void ) p->prng( p, true ); // re-seed prng
   for ( int addr = 0; addr < Size; addr++ )
   {
       u64 r = p->prng( p, false ) & BitMask;;
       u64 v = Read( addr ) & BitMask;
       if ( v != r )
       {
           if( verboseMode ) { printf( "\n%s BRAM test failed at address %08X: expected: %08X read: %08X\n", Name, toUint( addr ), toUint( r ), toUint( v ) ); fflush( stdout ); }
           return( addr );
       }
   }
   if( verboseMode ) printf("\n%s BRAM test passed", Name ); fflush( stdout );
   return( 0xFFFFFFFF );
}

//========================================================================
// This, believe it or not, counts the number of set bits in a 32 bit word
// To the Code Review team:
// 1. It is cited in numerous sources
//    https://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetParallel
//    https://en.wikipedia.org/wiki/Hamming_weight
//    http://gurmeet.net/puzzles/fast-bit-counting-routines/
//    http://aggregate.ee.engr.uky.edu/MAGIC/#Population%20Count%20(Ones%20Count)
// 2. It requires no maintenance as long as nothing is changed.
//    Do not use it with signed integers, for instance, since
//    the >> operator is poorly implemented in that case.
// 3. To make it MISRA-compliant you have to substitute the magic numbers
//    in a MISRA-compliant way
//========================================================================
u32 HammingWeight( u32 w )
{
    return ( __builtin_popcount( w ) );
//    u32 d = w;
//    d = d - (( d >> 1 ) & 0x55555555 );
//    d = ( d & 0x33333333 ) + (( d >> 2 ) & 0x33333333 );
//    return ( ((( d + ( d >> 4 )) & 0x0F0F0F0F ) * 0x01010101 ) >> 24 );
}

u32 HammingDist( u32 w1, u32 w2 )
{
    return ( HammingWeight( w1 ^ w2 ) );
}

const double ra[ 32 ] = {
        1.0 / 1.0,
        1.0 / 2.0,
        1.0 / 4.0,
        1.0 / 8.0,

        1.0 / 16.0,
        1.0 / 32.0,
        1.0 / 64.0,
        1.0 / 128.0,

        1.0 / 256.0,
        1.0 / 512.0,
        1.0 / 1024.0,
        1.0 / 2048.0,

        1.0 / 4096.0,
        1.0 / 8192.0,
        1.0 / 16384.0,
        1.0 / 32768.0,

        1.0 / 65536.0,
        1.0 / 131072.0,
        1.0 / 262144.0,
        1.0 / 524288.0,

        1.0 / 1048576.0,
        1.0 / 2097152.0,
        1.0 / 4194304.0,
        1.0 / 8388608.0,

        1.0 / 16771216.0,
        1.0 / 33554432.0,
        1.0 / 67108864.0,
        1.0 / 134217728.0,
        1.0 / 268435456.0,
        1.0 / 536870912.0,
        1.0 / 1073741824.0,
        1.0 / 2147483648.0
};

u32 RMeantoRNoise( double rmean, u32 nBits )
{
    u32 rNoise = 0;
    double rrm = 1.0 / rmean;
    double rm  = 0.0;
    for( int i = 0; i < nBits; i++ )
    {
        rNoise <<= 1;
        if ( ( rm + ra[ i ] ) <= rrm )
        {
            rm += ra[ i ];
            rNoise |= 1;
        }

    }
    return ( rNoise );
}

s16 Sext16( s16 i, u16 nbits )
{
    s16 sgn = ( 1 << ( nbits - 1 ));
    return( ( i & sgn ) ? i | -sgn : i );
}
