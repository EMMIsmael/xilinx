/*
 * console_util.c
 *
 *  Created on: 18 Aug 2018
 *      Author: JULIAN MORTIMER
 *
 *  Deal with console I/O
 *
 */
#include "include.h"

const CtrlChar CtrlChars[] = {
        { ""    , CC_NUL     ,"<NUL>"      },
        { ""    , CC_BEL     ,"<BEL>"      },
        { ""    , CC_BS      ,"<BS>"       },
        { ""    , CC_TAB     ,"<TAB>"      },
        { ""    , CC_LF      ,"<LF>"       },
        { ""    , CC_VTAB    ,"<VTAB>"     },
        { ""    , CC_FF      ,"<FF>"       },
        { ""    , CC_CR      ,"<CR>"       },
        { ""    , CC_ESC     ,"<ESC>"      },
        { "[A"  , CC_UP      ,"<UP>"       },
        { "[B"  , CC_DOWN    ,"<DOWN>"     },
        { "[C"  , CC_RIGHT   ,"<RIGHT>"    },
        { "[D"  , CC_LEFT    ,"<LEFT>"     },
        { "[11~", CC_F1      ,"<F1>"       },
        { "[12~", CC_F2      ,"<F2>"       },
        { "[13~", CC_F3      ,"<F3>"       },
        { "[14~", CC_F4      ,"<F4>"       },
        { "[15~", CC_F5      ,"<F5>"       },
        { "[17~", CC_F6      ,"<F6>"       },
        { "[18~", CC_F7      ,"<F7>"       },
        { "[19~", CC_F8      ,"<F8>"       },
        { "[20~", CC_F9      ,"<F9>"       },
        { "[21~", CC_F10     ,"<F10>"      },
        { "[23~", CC_F11     ,"<F11>"      },
        { "[24~", CC_F12     ,"<F12>"      },
        { "[1~" , CC_HOME    ,"<HOME>"     },
        { "[2~" , CC_INS     ,"<INS>"      },
        { "[3~" , CC_DEL     ,"<DEL>"      },
        { "[4~" , CC_END     ,"<END>"      },
        { "[5~" , CC_PGUP    ,"<PGUP>"     },
        { "[6~" , CC_PGDN    ,"<PGDN>"     },
        { "OP"  , CC_NUML    ,"<NUML>"     },
        { "OQ"  , CC_NPSLASH ,"<NP-SLASH>" },
        { "OR"  , CC_NPSTAR  ,"<NP-STAR>"  },
        { "OS"  , CC_NPMIN   ,"<NP-DASH>"  },
        { ""    , 0xFF       ,""           }
};

char GetChar( void );
char cCode( char *s, const CtrlChar cctbl[] );
char *cDisp( char c, const CtrlChar cctbl[] );
void testKB( void );

_Bool scandbl ( char buf [ ], double *res );
_Bool scanuns ( char buf [ ], unsigned int *res );
_Bool scanint ( char buf [ ], signed int *res );
_Bool scanstr ( char buf [ ], int treatAsSigned, unsigned int *res );
unsigned int satMulAddu ( unsigned int p1, unsigned int p2, unsigned int s );
signed   int satMulAdds ( signed int p1, signed int p2, signed int s );

char *strucase ( char *dest, char *srce );
char *strlcase ( char *dest, char *srce );
char pause( char msg[] );
void dispSysTimerRunning( void );
int  sysTimerRunning( void );
void utilScreenDumpData( dmpSelector sel, dmpFormat fmt, dmpType typ, dmpAddr addr, dmpArg arg1, dmpArg arg2, dmpArg arg3, dmpTitle title );
u32  utilGetDumpData( dmpSelector sel, u32 addr, dmpArg arg, dmpArg ptr );
void u32bin2str( char s[], u32 u, s32 nbytes, u32 spccnt );
void DumpCplx( cplx *pz, u32 nData );


void DumpCplx( cplx *pz, u32 nData )
{
    cprintf( "\n" );
    for( u32 i = 0; i < nData; i++ ) {
        cprintf( "%4u %10.4g %10.4g\n", toUint( i ), cimag( pz [ i ] ), creal( pz [ i ] ) );
        if ( ( ( i + 1 ) % N_LINES_DUMP ) == 0 ) {
            char c = 0;
            if (( c = GetChar() ) == CC_ESC ) {
                break;
            }
        }
    }
}


// Scan a buffer containing only non whitespace characters for an unsigned or signed int
_Bool scandbl ( char buf [ ], double *res )
{
    *res = 0.0;
    return sscanf ( buf, "%lg", res );
}
// Scan a buffer containing only non whitespace characters for an unsigned or signed int
_Bool scanuns ( char buf [ ], unsigned int *res )
{
    return ( scanstr ( buf, false, res ) );
}

_Bool scanint ( char buf [ ], signed int *res )
{
    unsigned int u = 0;
    _Bool b = scanstr ( buf, true, &u );
    *res = b ? *( signed int * ) &( u ) : 0;
    return( b );
}

_Bool scanstr ( char buf [ ], int treatAsSigned, unsigned int *res )
{
    int len = strlen ( buf );
    int i = 0;
    *res  = 0;

    if ( len == 0 ) return ( false );
    if ( ( strstr ( buf, "0x" ) == buf ) || ( strstr ( buf, "0X" ) == buf ) ) {
        unsigned int ui = 0;
        for ( i = 2; ( ( i < len ) && isxdigit( ( int )buf [ i ] ) ); i++ )
            if ( isdigit( ( int )buf [ i ] ) )
                ui = satMulAddu ( ui, 16, buf [ i ] - '0' );
            else
            ui = satMulAddu ( ui, 16, toupper(( int )buf[ i ] ) - 'A' + 10 );
        if ( i != len )
            return ( false );
        else {
            *res = ui;
            return ( true );
        }
    }
    else
    {
        if ( treatAsSigned )
        {
            int si = 0;
            if ( buf [ 0 ] == '-' )
                for ( i = 1; ( ( i < len ) && isdigit( ( int )buf [ i ] ) ); i++ )
                    si = satMulAdds ( si, 10, - ( buf [ i ] - '0' ) );
            else {
                if ( buf [ 0 ] == '+' ) i++;
                for ( ; ( ( i < len ) && isdigit( ( int )buf [ i ] ) ); i++ )
                    si = satMulAdds ( si, 10, buf [ i ] - '0' );
            }
            if ( i != len )
                return ( false );
            else {
                *res = *( unsigned int * ) &si;
                return ( true );
            }
        }
        else
        {
            unsigned int ui = 0;
            if ( ( buf [ 0 ] == '+' ) || ( buf [ 0 ] == '-' ) ) i++;
            for ( ; ( ( i < len ) && isdigit( ( int )buf [ i ] ) ); i++ )
                ui = satMulAddu ( ui, 10, buf [ i ] - '0' );
            if ( i != len )
                return ( false );
            else
            {
                //two's complement result if minus sign
                if ( buf [ 0 ] == '-' )
                    *res = ~ui + 1;
                else
                *res = ui;
                return ( true );
            }
        }
    }
}

// Saturate the expression p1 * p2 + s to UINT_MAX. p2 always non-zero
unsigned int satMulAddu ( unsigned int p1, unsigned int p2, unsigned int s )
{
    if ( ( UINT_MAX / p2 ) < p1 ) return ( UINT_MAX );
    if ( ( p1 * p2 + s ) < s ) return ( UINT_MAX );
    return ( p1 * p2 + s );
}

// Signed version : p2 always positive non-zero
signed int satMulAdds ( signed int p1, signed int p2, signed int s ) {
    if ( p1 < 0 ) {
        if ( ( INT_MIN / p2 ) > p1 ) return ( INT_MIN );
        if ( ( p1 * p2 + s ) > s ) return ( INT_MIN );
        return ( p1 * p2 + s );
    }
    else {
        if ( ( INT_MAX / p2 ) < p1 ) return ( INT_MAX );
        if ( ( p1 * p2 + s ) < s ) return ( INT_MAX );
        return ( p1 * p2 + s );
    }
}

// GNU warning-free case mapping
char *strlcase ( char *dest, char *srce ) {
    char *d = dest;
    for ( ; *srce; ++srce )
        *d++ = lcase( *srce );
    *d = '\0';
    return ( dest );
}

char *strucase ( char *dest, char *srce ) {
    char *d = dest;
    for ( ; *srce; ++srce )
        *d++ = ucase( *srce );
    *d = '\0';
    return ( dest );
}


static char cbuf[ 10 ];
char *cDisp( char c, const CtrlChar cctbl[] )
{
    int i = 0;
    if (( c < 0x80 ) && ( c >= 0x20 ))
    {
        cbuf[ 0 ] = c;
        cbuf[ 1 ] = '\0';
    }
    else
    {
        while ( cctbl[ i ].code != 0xFF )
        {
            if ( c == cctbl[ i ].code ) break;
            i++;
        }
        if (( cctbl[ i ].code ) == 0xFF )
        {
            if ( c < 0x20 )
                sprintf( cbuf, "<^%c>", c + '@' );
            else
                sprintf( cbuf, "<%d>", c );
        }
        else
        {
            sprintf( cbuf, cctbl[ i ].disp );
        }
    }
    return( cbuf );
}

// Input a character from the console
// returns:
//   printables from 0x20 to 0x7f
//   backspace
//   carriage return
//   naked escape
// If user types ESC tries to input a sequence.
//   If this times out, return ESC
//   If other characters are received it buffers them, then look for
//   a sequence used locally.
//   If found, it return local code.
//   Otherwise wait for more input.
char GetChar( void )
{
    char c = 0;
    u32 tr = sysTimerRunning( );
    while ( true )
    {
        if (( c = GETC ) == CC_ESC )
        {
            if ( escSeqDet && tr )
            {
                int  eidx = 0;
                char rxbuf[ ESEQ_RX_BUF_SIZE + 1 ];
                double t = SysTimer->ReadTimerLongSec( SysTimer, SysTimerId );
                do {
                    if( KEYPRESSED )
                    {
                        c = GETC;
                        if( eidx < ESEQ_RX_BUF_SIZE + 1 )
                        {
                            rxbuf[ eidx++ ] = c;
                            rxbuf[ eidx ] = '\0';
                        }
                    }
                }
                while (( SysTimer->ReadTimerLongSec( SysTimer, SysTimerId ) - t ) < ESEQ_TIMEOUT );
                if ( eidx == 0 ) return ( CC_ESC );
                if (( c = cCode( rxbuf, CtrlChars )) != 0xFF ) return ( c );
                continue;
            }
            else
            {
                return( c );
            }
        }
        return( c );
    }
    return( 0 ); // unreachable
}

// Convert escape sequence to locally-used code
char cCode( char *s, const CtrlChar cctbl[] )
{
    int i = 0;
    while ( cctbl[ i ].code != 0xFF )
    {
        if ( !strcmp( cctbl[ i ].eseq, s )) break;
        i++;
    }
    return ( cctbl[ i ].code );
}

// Test Keyboard
void testKB( void )
{
    char s[ 4 ] = "...";
    _Bool eSDSave = escSeqDet;
    if ( sysTimerRunning( )) escSeqDet = true;
    cprintf( "\n\nKeyboard Test. ""zzz"" to exit, ""e"" to enable/disable ESC sequence detection.\n" );
    cprintf( "\nESC sequence detection currently enabled.\n\n" );

    while ( strcmp( s, "zzz" ))
    {
        char c = GetChar( );
        if ( c == 'e' )
        {
            escSeqDet = !escSeqDet;
            if ( escSeqDet ) cprintf( " -ESC on- " ); else cprintf( " -ESC off- " );
        }
        cprintf( "%s", cDisp( c, CtrlChars ));
        for ( int i = 2; i > 0; i-- ) s[ i ] = s[ i-1 ];
        s[ 0 ] = c;
    }
    cprintf( "\n\n" );
    escSeqDet = eSDSave;
}

char pause( char msg[] )
{
    char prompt[] = "";
    cprintf( "%s", msg );
    char c = GetChar();
    ERASE_LEFT( strlen( msg ) + strlen( prompt ) );
    return( c );
}

int sysTimerRunning( void )
{
    u32 u = SysTimer->ReadTimerLo( SysTimer, SysTimerId );
    return ( SysTimer->ReadTimerLo( SysTimer, SysTimerId ) != u );
}

void dispSysTimerRunning( void )
{
    double f = SysTimer->ReadTimerLongSec( SysTimer, SysTimerId );
    double g = SysTimer->ReadTimerLongSec( SysTimer, SysTimerId );
    if ( g != f )
    {
        cprintf( "\nSystem Timer is running. Time for test %.3fus\n", ( g - f ) * 1.0e6 );
    }
    else
    {
        cprintf( "\nSystem Timer is stopped. Value = %.3fus\n", g * 1.0e6 );
    }
}


u32 utilGetDumpData( dmpSelector sel, u32 addr, dmpArg arg, dmpArg ptr )
{
    switch ( sel )
    {
    case DUMP_SEL_PS7_RAM:
        return( *( u32* ) addr );
        break;
    case DUMP_SEL_SOURCE:
    {
        AXISource *pd = ( AXISource * ) ptr;
        return ( RD_DVC_DATA( pd, Bram, addr ));
        break;
    }
    case DUMP_SEL_SINK:
    {
        AXISink *pd = ( AXISink * ) ptr;
        return ( RD_DVC_DATA( pd, Bram, addr ));
        break;
    }
        //    case DUMP_SEL_TI_BUF:
        //      return ( timeInterleaverAdReadDataBuf( addr ));
        //      break;
        //    case DUMP_SEL_TI_MAP:
        //      return ( timeInterleaverAdReadIntlvr( addr ));
        //      break;
        //    case DUMP_SEL_CE_BUF:
        //      return ( convEncdrAdReadDataBuf( addr ));
        //      break;
        //    case DUMP_SEL_CE_MAP:
        //      return ( convEncdrAdReadEncdr( addr ));
        //      break;
        //    case DUMP_SEL_FI_TX_BUF:
        //      freqIntlvCtrlrWriteSlaveReg(( tFreqIntlvCtrlr* )ptr, FREQ_INTLV_CTRLR_AD_BIT_SEL_REG, arg );
        //      return ( freqIntlvCtrlrAdReadDataBuf(( tFreqIntlvCtrlr* )ptr, addr ));
        //      break;
        //    case DUMP_SEL_FI_INTLVR:
        //      return ( freqIntlvCtrlrAdReadIntlvr(( tFreqIntlvCtrlr* )ptr, addr, arg ));
        //      break;
        //    case DUMP_SEL_FI_MOD:
        //      return ( freqIntlvCtrlrAdReadMod(( tFreqIntlvCtrlr* )ptr, addr ));
        //      break;
    case DUMP_SEL_TSINJ_PSBUF:
    {
        TSInj *pd = ( TSInj * ) ptr;
        return ( RD_DVC_DATA( pd, PsBuf, addr ));
        break;
    }
    case DUMP_SEL_TSINJ_TSBUF:
    {
        TSInj *pd = ( TSInj * ) ptr;
        return ( RD_DVC_DATA( pd, TsBuf, addr ));
        break;
    }
    case DUMP_SEL_TSINJ_TXBUF:
    {
        TSInj *pd = ( TSInj * ) ptr;
        return ( RD_DVC_DATA( pd, TxBuf, addr ));
        break;
    }
    case DUMP_SEL_DDM_UTIL_RAW_LUT:
        return ( DDUtilRawLUT( addr ));
        break;
    case DUMP_SEL_DDM_UTIL_DMOD_LUT:
        return ( DDUtilDmodLUT( addr ));
        break;
//    case DUMP_SEL_DDM_TST_VEC:
//        return ( DataDemodRdTstVec( addr ));
//        break;
//    case DUMP_SEL_DDM_TF:
//        return ( DataDemodRdTFData( addr ));
//        break;
    case DUMP_SEL_TRE_PS_LUT:
    {
        TREmul *pd = ( TREmul * ) ptr;
        return ( RD_DVC_DATA( pd, PsLUT, addr ));
        break;
    }
    case DUMP_SEL_TRE_HW_LUT:
    {
        TREmul *pd = ( TREmul * ) ptr;
        return ( RD_DVC_DATA( pd, HwLUT, addr ));
        break;
    }
    default:
        break;
    }
    return 0;
}

void utilScreenDumpData( dmpSelector sel, dmpFormat fmt, dmpType typ, dmpAddr addr, dmpArg arg1, dmpArg arg2, dmpArg arg3, dmpTitle title )
{
    u32 i = 0, j = 0, u = 0;

    char s[ 2 * 32 + 1 ];
    unsigned char c;
    cprintf( "\n%s\n\n", title );

    if ( typ == DUMP_TYP_SLV_REG )
    {
        u32 nregs = arg1;
        u32 cols  = ( nregs > 16 ? 2 : 1 );
        u32 rows  = ( nregs + 1 ) / cols;
        u32 rofs  = 0;
        for ( i = 0; i < rows; i ++ )
        {
            for ( j = 0; j < cols; j++ )
            {
                rofs = i + rows * j;
                if ( rofs < nregs )
                {
                    u = hwReadReg( addr, rofs * 4 );
                    switch( fmt )
                    {
                    case DUMP_FMT_HEX_DEC:
                        cprintf( "%03u: 0x%08x ( %10u )   ", ( unsigned int )rofs, ( unsigned int )u, ( unsigned int )u );
                        break;
                    case DUMP_FMT_BIN32:
                        u32bin2str( s, u, 4, 4 );
                        cprintf( "%03u: %s   ", ( unsigned int )rofs, s );
                        break;
                    case DUMP_FMT_HEX:
                    default:
                        cprintf( "%3u: 0x%08x   ", ( unsigned int )rofs, ( unsigned int )u );
                        break;
                    }
                }
            }
            cprintf( "\n" );
        }
        cprintf( "\n" );
    }
    else if ( typ == DUMP_TYP_MEMORY )
    {

        switch( fmt )
        {
        //            case DUMP_FMT_SOFT_BITS_HW:
        //            {
        //              u32 nbits    = MIN( MAX( arg1, 3 ), 5 );
        //              s16 s0 = 0, s1 = 0, s2 = 0, s3 = 0, s4 = 0, s5 = 0;
        //              u32 ibits = 0, qbits = 0;
        //              do {
        //                  for ( i = 0; i < 32; i++ )
        //                  {
        //                        u32 bitmsk =  ( 1 << nbits ) - 1;
        //                        u32 sgnbit =  ( 1 << ( nbits - 1 ));
        //                        u32 sxbits = -( 1 << nbits );
        //
        //                        u32 u = utilGetDumpData( sel, ( addr + i ), arg2, arg3 );
        //                      ibits = ( u & 0xffff );
        //                      qbits = ( u >> 16 ) & 0xffff;
        //                        s3 =  ibits                 & bitmsk; if( s3 & sgnbit ) s3 |= sxbits;
        //                      s4 = ( ibits >> nbits )       & bitmsk; if( s4 & sgnbit ) s4 |= sxbits;
        //                        s5 = ( ibits >> ( 2 * nbits )) & bitmsk; if( s5 & sgnbit ) s5 |= sxbits;
        //                        s0  = qbits                 & bitmsk; if( s0 & sgnbit ) s0 |= sxbits;
        //                      s1 = ( qbits >> nbits )       & bitmsk; if( s1 & sgnbit ) s1 |= sxbits;
        //                        s2 = ( qbits >> ( 2 * nbits )) & bitmsk; if( s2 & sgnbit ) s2 |= sxbits;
        //                      cprintf( "%6u  %3i  %3i  %3i  %3i  %3i  %3i  ", ( unsigned int )( addr + i ), ( int )s5, ( int )s4, ( int )s3, ( int )s2, ( int )s1, ( int )s0 );
        //                        { if( s5 > 0 ) cprintf( "0" ); else cprintf( "1" ); }
        //                        { if( s4 > 0 ) cprintf( "0" ); else cprintf( "1" ); }
        //                        { if( s3 > 0 ) cprintf( "0" ); else cprintf( "1" ); }
        //                        { if( s2 > 0 ) cprintf( "0" ); else cprintf( "1" ); }
        //                        { if( s1 > 0 ) cprintf( "0" ); else cprintf( "1" ); }
        //                        { if( s0 > 0 ) cprintf( "0" ); else cprintf( "1" ); }
        //                      cprintf( "\n" );
        //                  }
        //                  addr += 32;
        //              }  while (( c = GetChar( )) != CC_ESC );
        //                break;
        //
        //            }
        //
        //            case DUMP_FMT_HEX_DEMOD:
        //            {
        //                u32 modType = MIN( arg1, 3 );
        //                Modulator *mod  = modulators[ modType ];
        //              do {
        //                  for( u32 i = 0; i < 32; i++ )
        //                  {
        //                          u32 u = utilGetDumpData( sel, addr + i, arg2, arg3 );
        //                          u32 d = mod->Demodulate( mod, u, DATA_MOD_SCALE );
        //                          cprintf( "%6u  %08X  %2d " , ( unsigned int )( addr + i ), ( unsigned int )u, ( unsigned int )d );
        //                          for( int j = ( int )log2u32( mod->npoints ) - 1; j >= 0; j-- )
        //                          {
        //                              if( d & ( 1 << j )) cprintf( "1" ); else cprintf( "0" );
        //
        //                          }
        //                      cprintf( "\n" );
        //                  }
        //                  addr += 32;
        //              } while (( c = GetChar( )) != CC_ESC );
        //
        //              break;
        //
        //            }
        //
        //            case DUMP_FMT_SOFT_BITS:
        //            {
        //              u32 nbits    = MIN( MAX( arg1, 3 ), 5 );
        //              u32 addrMask = arg2 - 1;
        //              s16 s0 = 0, s1 = 0, s2 = 0, s3 = 0, s4 = 0, s5 = 0;
        //              u32 ibits = 0, qbits = 0;
        //              do {
        //                  for ( i = 0; i < 32; i++ )
        //                  {
        //                        u32 bitmsk =  ( 1 << nbits ) - 1;
        //                        u32 sgnbit =  ( 1 << ( nbits - 1 ));
        //                        u32 sxbits = -( 1 << nbits );
        //
        //                      ibits = utilGetDumpData( sel,   ( addr + i ),                   arg2, arg3 );
        //                      qbits = utilGetDumpData( sel, ((( addr + i ) ^ addrMask ) + 1 ), arg2, arg3 );
        //                        s0 =  ibits                 & bitmsk; if( s0 & sgnbit ) s0 |= sxbits;
        //                      s1 = ( ibits >> nbits )       & bitmsk; if( s1 & sgnbit ) s1 |= sxbits;
        //                        s2 = ( ibits >> ( 2 * nbits )) & bitmsk; if( s2 & sgnbit ) s2 |= sxbits;
        //                        s3  = qbits                 & bitmsk; if( s3 & sgnbit ) s3 |= sxbits;
        //                      s4 = ( qbits >> nbits )       & bitmsk; if( s4 & sgnbit ) s4 |= sxbits;
        //                        s5 = ( qbits >> ( 2 * nbits )) & bitmsk; if( s5 & sgnbit ) s5 |= sxbits;
        //                      cprintf( "%6u  ", ( unsigned int )( addr+i ));
        //                      cprintf( "%3i  %3i  %3i  %3i  %3i  %3i\n", ( int )s5, ( int )s4, ( int )s3, ( int )s2, ( int )s1, ( int )s0 );
        //                  }
        //                  addr += 32;
        //              }  while (( c = GetChar( )) != CC_ESC );
        //                break;
        //
        //            }
        case DUMP_FMT_IQ_SIGNED:
        {
            do {
                u32 nbits = MIN( MAX( arg1, 10 ), 16 );
                s16 si = 0, sq = 0;
                for ( i = 0; i < 32; i++ )
                {
                    s32 sgn = ( 1 << ( nbits - 1 ));
                    u = utilGetDumpData( sel, addr + i, arg2, arg3 );
                    si =  u        & 0xffff; if ( si & sgn ) si |= -sgn;
                    sq = ( u >> 16 ) & 0xffff; if ( sq & sgn ) sq |= -sgn;
                    cprintf( "%6u  %6i  %6i\n", ( unsigned int )( addr+i ), ( int )sq, ( int )si );
                }
                addr += 32;
            }  while (( c = GetChar( )) != CC_ESC );
            break;
        }

        case DUMP_FMT_CPLX:
        {
            do {
                cplx *array = ( cplx* ) &addr;
                u32 offs = 0;
                for ( int i = 0; i < 32; i++ )
                {
                    cprintf( "%6u  %15g + 1i * %15g\n", toUint( offs ), creal( array [ offs ] ), cimag( array [ offs ] ) );
                    offs++;
                }
            }  while (( c = GetChar( )) != CC_ESC );
        }

        case DUMP_FMT_DEC:
            do {
                for ( i = 0; i < 32; i++ )
                {
                    u = utilGetDumpData( sel, addr + i, arg2, arg3 );
                    cprintf( "%u\n", ( unsigned int )u );
                }
                addr += 32;
            }  while (( c = GetChar( )) != CC_ESC );
            break;

        case DUMP_FMT_DEC_SIGNED :
            do {
                u32 u;
                for ( i = 0; i < 32; i++ )
                {
                    u = utilGetDumpData( sel, addr + i, arg2, arg3 );
                    cprintf( "%6i\n", *( int* )&u );
                }
                addr += 32;
            } while (( c = GetChar( )) != CC_ESC );
            break;

        case DUMP_FMT_GRAPH:
            do {
                char column[ 115 + 1 ];
                int i1, i2, i3, i4;
                i2 = *( int* ) &arg1;
                i3 = *( int* ) &arg2;
                i4 = *( int* ) &arg3;

                u32 u, a;
                for ( int i = 0; i < 32; i++ )
                {
                    a = i << i4;
                    u = utilGetDumpData( sel, addr + a, arg2, arg3 );
                    i1 = *( int* ) &u;

                    for ( int j = 0; j < 115; j++ )
                    {
                        if (( j & 0x7 ) == 0 )
                            if (( i & 0x07 ) == 0 )
                                column[ j ] = '+';
                            else
                                column[ j ] = '|';

                        else if (( i & 0x07 ) == 0 )
                            column[ j ] = '-';
                        else
                            column[ j ] = ' ';
                    }

                    column[ i3 ] = '+';
                    column[ ( i1 << i2 ) + i3 ] = '@';

                    if (( i & 0x7 ) == 0 )
                        cprintf( "%8d-+-", ( int ) ( addr + a ));
                    else
                        cprintf( "         | " );

                    cprintf( "%s\n", column );

                }
                addr += ( 32 << i4 );
            } while (( c = GetChar( )) != CC_ESC );
            break;

        case DUMP_FMT_HEX:
        default:
        {
            char sfmt[ 20 ];
            u32 nbytesWord = arg1;
            switch( nbytesWord )
            {
            case 1:
                strcpy( sfmt, "%02X" );
                break;
            case 2:
                strcpy( sfmt, "%04X " );
                break;
            case 4:
            default:
                strcpy( sfmt, "%08X " );
                break;
            }
            do {
                for( i = 0; i < 16; i++ )
                {
                    cprintf( "%08X: ", ( unsigned int )addr );
                    for ( j = 0; j < ( 32 / nbytesWord ); j++ )
                    {
                        u = utilGetDumpData( sel, addr + j, arg2, arg3 );
                        cprintf( sfmt, ( unsigned int )u );

                        if (( nbytesWord == 1 ) && ((( j + 1 ) % 4 ) == 0 )) cprintf( " " );
                    }
                    cprintf( "\n" );
                    addr += 32 / nbytesWord;
                }
            } while (( c = GetChar( )) != CC_ESC );

            break;

        }
        }
    }
    else
    {
        cprintf( "\nScreen Dump : Selected Type Not Found - %i\n", ( int )typ );
    }
}

void u32bin2str( char s[], u32 u, s32 nbytes, u32 spccnt )
{
    u32 i  = 0;
    u32 p  = 0;
    u32 uu = u;
    if ( nbytes <= 0 )
    {
        s[ 0 ] = '\0';
    }
    else
    {
        uu <<= 8 * ( 4 - nbytes );
        for( i = 0; i < nbytes * 8; i++ )
        {
            s[ p++ ] = ((( uu << i ) & 0x80000000 ) == 0 ? '0' : '1' );
            if (( spccnt != 0 ) && ((( i + 1 ) % spccnt ) == 0 ) && ( i < nbytes * 8 - 1 ))
                s[ p++ ] = ' ';
        }
        s[ p ] = '\0';
    }
}

