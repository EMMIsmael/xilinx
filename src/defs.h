/*
 * defs.h
 *
 *  Created on: 16 May 2017
 *      Author: JULIAN MORTIMER
 */

#ifndef DEFS_H_
#define DEFS_H_

#include <complex.h>

#define UART_BASE    ( XPAR_PS7_UART_1_BASEADDR )
#define GETC         ( XUartPs_RecvByte( UART_BASE ) )
#define PUTC(C)      do { XUartPs_SendByte( UART_BASE, (C) ); } while (0)
#define KEYPRESSED   ( XUartPs_IsReceiveData( UART_BASE ) )
#define FLUSH        do { while( KEYPRESSED ){ (void) GETC; }} while (0)

#define islcase( c ) ( c >= 'a' && c <= 'z' )
#define isucase( c ) ( c >= 'A' && c <= 'Z' )
#define lcase( c )   ( isucase( c ) ? c - 'A' + 'a' : c )
#define ucase( c )   ( islcase( c ) ? c - 'a' + 'A' : c )
#define cprintf( ... ) do { printf( __VA_ARGS__ ); fflush( stdout ); } while( 0 )
#define cprintfv( ... ) do { if ( verboseMode ) { printf( __VA_ARGS__ ); fflush( stdout ); }} while( 0 )
#define toUint( u )  ( ( unsigned int )      ( u ) )
#define toULint( u ) ( ( unsigned long int ) ( u ) )
#define toInt( i )   ( ( int )               ( i ) )
#define toLInt ( i ) ( ( long int )          ( i ) )

#define ONE_MEGAWORD ( 0x100000 )

#define ESEQ_TIMEOUT     ( 0.001 )
#define ESEQ_RX_BUF_SIZE ( 10 )

#define CC_NUL       ( 0x00 )
#define CC_BEL       ( 0x07 )
#define CC_BS        ( 0x08 )
#define CC_TAB       ( 0x09 )
#define CC_LF        ( 0x0A )
#define CC_VTAB      ( 0x0B )
#define CC_FF        ( 0x0C )
#define CC_CR        ( 0x0D )
#define CC_ESC       ( 0x1B )
#define CC_UP        ( 0x80 )
#define CC_DOWN      ( 0x81 )
#define CC_RIGHT     ( 0x82 )
#define CC_LEFT      ( 0x83 )
#define CC_F1        ( 0x84 )
#define CC_F2        ( 0x85 )
#define CC_F3        ( 0x86 )
#define CC_F4        ( 0x87 )
#define CC_F5        ( 0x88 )
#define CC_F6        ( 0x89 )
#define CC_F7        ( 0x8A )
#define CC_F8        ( 0x8B )
#define CC_F9        ( 0x8C )
#define CC_F10       ( 0x8D )
#define CC_F11       ( 0x8E )
#define CC_F12       ( 0x8F )
#define CC_HOME      ( 0x90 )
#define CC_INS       ( 0x91 )
#define CC_DEL       ( 0x92 )
#define CC_END       ( 0x93 )
#define CC_PGUP      ( 0x94 )
#define CC_PGDN      ( 0x95 )
#define CC_NUML      ( 0x96 )
#define CC_NPSLASH   ( 0x97 )
#define CC_NPSTAR    ( 0x98 )
#define CC_NPMIN     ( 0x99 )
#define CC_NPPLUS    ( 0xA0 )

#define CURSOR_RIGHT( n ) do { for ( int __j = 0; __j < ( n ) ; __j++ ) { PUTC(0x1B); PUTC('[' ); PUTC('C' ); fflush(stdout); }} while( 0 )
#define CURSOR_LEFT( n )  do { for ( int __j = 0; __j < ( n ) ; __j++ ) { PUTC(0x1B); PUTC('[' ); PUTC('D' ); fflush(stdout); }} while( 0 )
#define ERASE_LEFT( n )   do { for ( int __j = 0; __j < ( n ) ; __j++ ) { PUTC('\b'); PUTC(' ' ); PUTC('\b'); fflush(stdout); }} while( 0 )

#define MAX( x, y )     ( ( ( x ) > ( y ) ) ? ( x ) : ( y ) )
#define MIN( x, y )     ( ( ( x ) < ( y ) ) ? ( x ) : ( y ) )
#define MAXI( x, i, j ) ( ( ( x [ i ] ) > ( x [ j ] ) ) ? ( i ) : ( j ) )
#define MINI( x, i, j ) ( ( ( x [ i ] ) < ( x [ j ] ) ) ? ( i ) : ( j ) )
#define MAXS( x0, x1 )  ( ( x0 > x1 ? 0 : 1 ) )
#define MINS( x0, x1 )  ( ( x0 < x1 ? 0 : 1 ) )

// The following must be used as RValues and CANNOT BE USED IN EXPRESSIONS
#define MINVI( v, n, idx )      0; do { for ( int __i = 0; __i < n; idx = MINI( v, __i, idx ), __i++ ); } while ( 0 )
#define MINV( v, n, vmin )      0; do { int idx = MINVI( v, n, idx ); vmin = v[ idx ]; } while ( 0 )
#define MAXVI( v, n, idx )      0; do { for ( int __i = 0; __i < n; idx = MAXI( v, __i, idx ), __i++ ); } while ( 0 )
#define MAXV( v, n, vmax )      0; do { int idx = MAXVI( v, n, idx ); vmax = v[ idx ]; } while ( 0 )
#define SUMV( v, n, sum )       0; do { for ( int __i = 0; __i < n; __i++ ) sum += v [ __i ]; } while ( 0 )
#define MEANV( v, n, mean )     SUMV( v, n, mean ); do { mean = mean / n; } while ( 0 )
#define NORMV( x, xn, n, mean ) MEANV( x, n, mean ); do { for ( int i = 0; i < n; i++ ) xn [ i ] = x [ i ] - mean; } while ( 0 )
#define SUMSQ( v, n, sum )      0; do { for ( int __i = 0; __i < n; __i++ ) sum += SQUARE( v [ __i ] ); } while ( 0 )
#define SUMXY( x, y, n, sum )   0; do { for ( int __i = 0; __i < n; __i++ ) sum += x [ __i ] * y [ __i ]; } while ( 0 )


#define IS_ODD( x )  (   ( x ) & 1 )
#define IS_EVEN( x ) ( ( ( x ) & 1 ) ^ 1 )

// WARNING these should be used with EXTREME CAUTION
// For instance, CUBE( x++ ) yields ( x++ ) * ( x++ ) * ( x++ )
#define SQUARE( x )  ( x * x )
#define CUBE( x )    ( x * SQUARE( x ) )

#define COPY( src, dest, n ) do { for ( int i = 0; i < n; i++ ) dest [ i ] = src [ i ]; } while ( 0 )


#define SIZE( X, type ) ( ( u32 ) sizeof( X ) / sizeof( type ) )

#define SNMBR(x) (( x ) == 1 ? "" : "s" )

#define UI( n ) (( unsigned int ) n )

#ifndef true
#define true  ( 1 )
#endif
#ifndef false
#define false ( 0 )
#endif

//Needed for -std=c11 but not -std=gnu11
#ifndef M_PI
#define M_PI    ( 3.14159265358979323846 )
#endif

#ifndef M_SQRT2
#define M_SQRT2 ( 1.41421356237309504880 )
#endif

#define PL_FABRIC_CLOCK_FREQ ( 100000000 )

#define SYS_TIMER( f ) ( SysTimer->f( SysTimer, SysTimerId ))

typedef double complex cplx;
typedef union uulng ulng;
struct uhl { u32 ul; u32 uh; };
union  uulng { u64 u; struct uhl uhl; };


#define TIME_BEG  { double __g = SysTimer->ReadTimerLongSec( SysTimer, SysTimerId );  do { } while ( 0 )
#define TIME_END    __g = SysTimer->ReadTimerLongSec( SysTimer, SysTimerId ) - __g; } do { } while ( 0 )

#define TIME( __s, s, su, k, doDisp ) do { double __g = SysTimer->ReadTimerLongSec( SysTimer, SysTimerId );\
    { __s; }\
    __g = SysTimer->ReadTimerLongSec( SysTimer, SysTimerId ) - __g;\
    if ( doDisp ) cprintf( "\n%s time = %0.8f %s\n", s, __g * k, su ); } while( 0 )

#define TIME_TFR( __s, nwd, nbw ) do { double g = SysTimer->ReadTimerLongSec( SysTimer, SysTimerId );\
	{ __s; }\
	g = SysTimer->ReadTimerLongSec( SysTimer, SysTimerId ) - g;\
	cprintf( " Load time for %d bytes = %.8fs, %.0f bytes/sec\n", ( unsigned int ) ( nbw * nwd ), g, nbw * nwd / g ); } while( 0 )

#define TRUEFALSE( b ) ( b ? "true" : "false" )
#define OFFON( b )     ( b ? "on"   : "off"   )

#endif /* DEFS_H_ */
