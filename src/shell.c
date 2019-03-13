/*


 * shell.c
 *
 *  Created on: 22 Mar 2016
 *      Author: RTL3
 */
#include "include.h"

//static const char macro1 = "runtest 0 6 0; shellWait 0.2; nsddmp; cirdelta";
//static const char macro2 = "shellEcho ###==rxdata\n;b00dbram 0 12 $nFFT";
//static const char macro3 = "shellEcho ###==rxprs\n; prsdbram 0 12 $nFFT";
//static const char macro4 = "shellEcho ###==rxdta\n; dtadbram 0 12 $nFFT";
//static const char macro5 = "shellEcho ###==cirtable\n; cirdtable $nFFT";

const u32 nMacros = 5;
_Bool    verboseMode = true;



shellError shellConsole       ( shell *ps, shellCommandDef *pdef, const char name [ ], const char startupCmd [ ], int helpLinesPause );
void       shellConsoleExit   ( shell *ps );
void       shellShowHelp      ( shell *ps );
void       shellLoop          ( shell *ps );
void       shellEndLoop       ( shell *ps );
void       shellSetVerboseMode( shell *ps );

static shellError shellRegisterCommand ( shell *ps, shellCommandDef *pdef );
static void       shellGetCommand  ( shell *ps );
static shellError shellInitialise  ( shell *ps, shellCommandDef *pdef, const char name [ ], const char startupCmd [ ], int helpLinesPause );
static shellError shellExecCommand ( shell *ps );
static shellError shellDoCommand   ( shell *ps );

static shellCommandDef shellCmdDef[ ] = {
    { "all"  , "help"      ,  0, ""                                                                           , shellShowHelp                           },
    { "all"  , "help"      ,  1, ""                                                                           , shellShowHelp                           },
    { "all"  , "help"      ,  2, ""                                                                           , shellShowHelp                           },
    { "all"  , "h, help"   ,  0, "                                        Show help"                          , shellShowHelp                           },
    { "all"  , "h, help"   ,  1, "<context>                               Show help for given context"        , shellShowHelp                           },
    { "all"  , "h, help"   ,  2, "pause n                                 Set help pause lines"               , shellShowHelp                           },
    { "all"  , "wait"      ,  1, "<timesec>                               Wait for a given time"              , shellWait                               },
    { "all"  , "echo"      ,  1, "<string>                                Echo a string to console"           , shellEcho                               },
    { "all"  , "loop"      ,  1, "                                        Loop for n iterations"              , shellLoop                               },
    { "all"  , "endl"      ,  0, "                                        End innermost loop"                 , shellEndLoop                            },
    { "all"  , "v"         ,  1, "< 1 | 0 >                               Set Verbose Mode"                   , shellSetVerboseMode                     },
    { "\0" }
};

static shellError shellInitialise ( shell *ps, shellCommandDef *pdef, const char name [ ], const char startupCmd [ ], int helpLinesPause )
{
    shellError ecode = SHELL_ERROR;
    shellCommandDef *pd = shellCmdDef;
    if ( strlen ( name ) <= SHELL_NAME_LEN ) {
        ps->cmdCount = 0;
        strcpy ( ps->name, name );
        strcpy ( ps->cmdLine, startupCmd );
        ps->hlpLinesPause = helpLinesPause;
        ps->cmdBufIdx = 0;
        ps->cmdBufLen = 0;
        for ( int i = 0; i < SHELL_CMDLINE_BUF_SIZE; i++ )
            strcpy ( ps->cmdLineBuf [ i ], "" );
        for ( int i = 0; i < SHELL_MAX_ARGS; i++ )
            strcpy ( ps->hlpCtx [ i ], "" );
        strcpy ( ps->hlpCtx [ 0 ], "ALL" );
        ecode = SHELL_SUCCESS;
    }
    pd = shellCmdDef;
    do {
        shellRegisterCommand( ps, pd++ );
    } while ( *pd->hlpCtx );
    pd = pdef;
    do {
        shellRegisterCommand( ps, pd++ );
    } while ( *pd->hlpCtx );
    ps->bIsInitialised = true;
    return ( ecode );
}


static shellError shellRegisterCommand ( shell *ps, shellCommandDef *pdef )
{
    shellError ecode = SHELL_ERROR;
    u32 idx;
    shellCommandDef *pd = pdef;

    if ( ( idx = ps->cmdCount ) < SHELL_MAX_CMDS ) {
        if ( ( strlen ( pdef->hlpMsg ) <= SHELL_HLP_MSG_LEN )
                && ( strlen ( pd->hlpCtx ) <= SHELL_MAX_ARG_LEN )
                && ( strlen ( pd->cmd ) <= SHELL_CMD_LEN ) && ( pd->narg <= SHELL_MAX_ARGS ) ) {
            strcpy ( ps->commands [ idx ].hlpMsg, pd->hlpMsg );
            strucase ( ps->commands [ idx ].hlpCtx, ( char * ) pd->hlpCtx );
            strcpy ( ps->commands [ idx ].cmd, pd->cmd );
            ps->commands [ idx ].narg = pd->narg;
            ps->commands [ idx ].callbackFunction = pd->callbackFunction;
            ps->commands [ idx ].arg.narg = 0;
            for ( int i = 0; i < SHELL_MAX_ARGS; i++ ) {
                ps->commands [ idx ].arg.uarg [ i ] = 0;
                ps->commands [ idx ].arg.darg [ i ] = 0.0;
                strcpy ( ps->commands [ idx ].arg.sarg [ i ], "" );
            }
            ps->cmdCount++;
            ecode = SHELL_SUCCESS;
        }
    }
    return ( ecode );
}

static void shellGetCommand ( shell *ps ) {
    int i = 0, cp = 0, cps = 0, k = 0;
    int bufidx = ps->cmdBufIdx;
    char c = 0;

    if ( !sysTimerRunning () ) {
        printf ( "\nSYSTEM TIMER FAILURE\n" );
        fflush ( stdout );
    }
    FLUSH;
    printf ( "%s ", ps->name );
    fflush ( stdout );
    ps->cmdLine [ 0 ] = '\0';

    do {
        if ( k == 0 ) strcpy ( ps->cmdLineBuf [ ps->cmdBufIdx ], ps->cmdLine );

        switch ( c = GetChar () ) {
        case CC_UP :
            case CC_DOWN :
            switch ( c ) {
            case CC_UP :
                for ( int idx = bufidx, kk = k; kk < ps->cmdBufLen - 1; ) {
                    if ( ++idx == SHELL_CMDLINE_BUF_SIZE )
                        idx = 0;
                    kk++;
                    char clbi[ SHELL_CMD_LEN + 1 ];
                    char clbf[ SHELL_CMD_LEN + 1 ];
                    strlcase( clbi, ps->cmdLineBuf [ idx ] );
                    strlcase( clbf, ps->cmdLineBuf [ ps->cmdBufIdx ] );
                    if ( strstr ( clbi, clbf ) == clbi ) {
                        if ( k == 0 )
                            cps = cp;
                        bufidx = idx;
                        k = kk;
                        break;
                    }
                }
                break;
            case CC_DOWN :
                default :
                for ( int idx = bufidx, kk = k; kk > 0; ) {
                    if ( --idx < 0 ) idx = SHELL_CMDLINE_BUF_SIZE - 1;
                    kk--;
                    char clbi[ SHELL_CMD_LEN + 1 ];
                    char clbf[ SHELL_CMD_LEN + 1 ];
                    strlcase( clbi, ps->cmdLineBuf [ idx ] );
                    strlcase( clbf, ps->cmdLineBuf [ ps->cmdBufIdx ] );
                    if ( strstr ( clbi, clbf ) == clbi ) {
                        bufidx = idx;
                        k = kk;
                        break;
                    }
                }
                break;
            }
            strcpy ( ps->cmdLine, ps->cmdLineBuf [ bufidx ] );
            CURSOR_RIGHT( i - cp );
            ERASE_LEFT( i );
            cprintf ( "%s", ps->cmdLine );
            i = strlen ( ps->cmdLine );
            if ( k == 0 )
                cp = cps;
            else
            cp = i;
            CURSOR_LEFT( i - cp );
            break;
        case CC_LEFT :
            if ( cp > 0 ) {
                CURSOR_LEFT( 1 );
                cp--;
            }
            break;
        case CC_HOME :
            if ( cp > 0 ) {
                CURSOR_LEFT( cp );
                cp = 0;
            }
            break;
        case CC_RIGHT :
            if ( cp < i ) {
                CURSOR_RIGHT( 1 );
                cp++;
            }
            break;
        case CC_END :
            if ( cp < i ) {
                CURSOR_RIGHT( i - cp );
                cp = i;
            }
            break;
        case CC_DEL :
            if ( cp < i ) {
                CURSOR_RIGHT( 1 );
                cp++;
            } // now it can be dealt with as a backspace
        case CC_BS :
            if ( cp > 0 ) {
                for ( int j = cp; j <= i; j++ )
                    ps->cmdLine [ j - 1 ] = ps->cmdLine [ j ];
                CURSOR_RIGHT( i - cp );
                ERASE_LEFT( i );
                cprintf ( "%s", ps->cmdLine );
                CURSOR_LEFT( i - cp );
                i--;
                cp--;
            }
            break;
        case 0x20 ... 0x7F :  // gnu11 extension.
            if ( i < ( SHELL_CMD_LINE_LEN - 1 ) ) {
                for ( int j = i; j >= cp; j-- )
                    ps->cmdLine [ j + 1 ] = ps->cmdLine [ j ];
                ps->cmdLine [ cp ] = c;
                CURSOR_RIGHT( i - cp );
                ERASE_LEFT( i );
                cprintf ( "%s", ps->cmdLine );
                CURSOR_LEFT( i - cp );
                i++;
                cp++;
                k = 0;
                bufidx = ps->cmdBufIdx;
            }
            break;
        default :
            break;
        }
    } while ( c != '\r' );
    cprintf( "\n" );
}

static shellError shellExecCommand ( shell *ps ) {
    char *pci = ps->cmdLine;
    shellError errCode = SHELL_ERROR;
    ps->nLoop = 0;
    ps->nxtCmdAddr  = pci;
    while ( 1 ) {
        pci = ps->nxtCmdAddr;
        char *pco = ps->cmd;
        while ( ( *pci != ';' ) && ( *pci != '\0' ) ) {
            *pco++ = *pci++;
        };
        *pco = '\0';
        ps->nxtCmdAddr = ( *pci == '\0' ? pci : ++pci );
        errCode = shellDoCommand( ps );
        if ( *ps->nxtCmdAddr == '\0' ) {
            break;
        }
    }
    if ( strlen ( ps->cmdLine ) > 0 ) {
        strcpy ( ps->cmdLineBuf [ ps->cmdBufIdx ], ps->cmdLine );
        if ( ps->cmdBufLen < SHELL_CMDLINE_BUF_SIZE ) {
            ps->cmdBufLen++;
        }
        ps->cmdBufIdx = ps->cmdBufIdx ? ps->cmdBufIdx - 1 : SHELL_CMDLINE_BUF_SIZE - 1;
    }
    return( errCode );
}

void shellLoop( shell *ps )
{
    if ( ( ps->nLoop < SHELL_NESTED_LOOPS ) && ( ps->uarg [ 0 ] ) ) {
        ps->loopCounter [ ps->nLoop ] = 0;
        ps->loopCount   [ ps->nLoop ] = ps->uarg[ 0 ];
        ps->loopAddr    [ ps->nLoop ] = ps->nxtCmdAddr;
        if( verboseMode ) {
            cprintfv( "Repeating for %i iterations {\n", ps->loopCount[ ps->nLoop ] );
        }
        ps->nLoop++;
    }
}

void shellEndLoop( shell *ps )
{
    if ( ps->nLoop ) {
        if ( ++ps->loopCounter [ ps->nLoop - 1 ] < ps->loopCount[ ps->nLoop - 1 ] ) {
           ps->nxtCmdAddr = ps->loopAddr[ ps-> nLoop - 1 ];
        }
        else
        {
            cprintfv( "} end repeat\n" );
            ps->nLoop--;
        }
    }
}

static shellError shellDoCommand ( shell *ps ) {
    shellError errCode = SHELL_ERROR;
    int idx = 0, pc = 0, i = 0;
    char *pcl = ps->cmd;
    char argBuf [ SHELL_MAX_ARG_LEN + 1 ];
    char cmdbuf [ SHELL_CMD_LEN + 1 ];
    char cbflcs [ SHELL_CMD_LEN + 1 ];
    char cmdlcs [ SHELL_CMD_LEN + 1 ];
    char afmt   [ SHELL_MAX_ARG_LEN + 1 ];

    for ( ; *pcl != '\0' && isspace( ( int ) ( *pcl ) ); pcl++ )
        ;
    pc = 0;
    for ( ; pc < SHELL_CMD_LEN && *pcl > 0x20 && *pcl < 0x80; pc++, pcl++ )
        cmdbuf [ pc ] = *pcl;
    cmdbuf [ pc ] = '\0';

    if ( strlen( cmdbuf ) == 0 ) {
        return ( SHELL_SUCCESS );
    }

    sprintf( afmt, "%%%us%%n", SHELL_MAX_ARG_LEN );
    for ( i = 0; i < SHELL_MAX_ARGS; i++ ) {
        int pos = 0;
        if ( sscanf ( pcl, afmt, argBuf, &pos ) != 1 )
            break;
        pcl += pos;
        strcpy ( ps->sarg [ i ], argBuf );
        ( void ) scanuns ( argBuf, &ps->uarg [ i ] );
        ( void ) scanint ( argBuf, &ps->iarg [ i ] );
        ( void ) scandbl ( argBuf, &ps->darg [ i ] );
    }
    ps->narg = i;

    ( void ) strlcase ( cbflcs, cmdbuf );
    for ( ; idx < ps->cmdCount; idx++ ) {
        if ( ps->commands [ idx ].callbackFunction == NULL ) {
            continue;
        }
        ( void ) strlcase ( cmdlcs, ps->commands [ idx ].cmd );
        for ( int i = 0; i < strlen( cmdlcs ); i++ ) if ( cmdlcs[ i ] == ',' ) cmdlcs[ i ] = '\0';
        if ( ps->commands [ idx ].narg == ps->narg && !strcmp ( cmdlcs, cbflcs ) ) {
            ps->parg       = &ps->commands [ idx ].arg;
            ps->parg->narg = ps->narg;
            for ( int a = 0; a < ps->narg; a++ ) {
                ps->parg->uarg [ a ] = ps->uarg [ a ];
                ps->parg->iarg [ a ] = ps->iarg [ a ];
                ps->parg->darg [ a ] = ps->darg [ a ];
                strcpy ( ps->parg->sarg [ a ], ps->sarg [ a ] );
            }
            ps->commands [ idx ].callbackFunction ( ps );
            errCode = SHELL_SUCCESS;
            break;
        }
    }
    return ( errCode );
}

static int shellMaxCmdLen( shell *psh )
{
    int lmx = 0, l = 0;
    for ( int i = 0; i < psh->cmdCount; i++ )
    {
       if (( l = strlen( psh->commands[ i ].cmd )) > lmx ) lmx = l;
    }
    return lmx;
}

static int shellMaxHctLen( shell *psh )
{
    int lmx = 0, l = 0;
    for ( int i = 0; i < psh->cmdCount; i++ )
    {
       if (( l = strlen( psh->commands[ i ].hlpCtx )) > lmx ) lmx = l;
    }
    return lmx;
}

static int hctcmp( shell *ps, char *hc )
{
    for ( int i = 0; ; i++ ) {
        if ( ps->hlpCtx[ i ][ 0 ] == '\0' )
            break;
        if ( !strcmp( ps->hlpCtx[ i ], hc ) || !strcmp( ps->hlpCtx[ i ], "ALL" ) ) {
            return ( 0 );
        }
    }
    return ( 1 );
}


static void SetHelpContext( shell *ps )
{
    for ( int i = 0; i < SHELL_MAX_ARGS; i++ )
        strcpy ( ps->hlpCtx [ i ], "" );
    if ( ps->narg == 2 ) {
        strcpy ( ps->hlpCtx [ 0 ], "ALL" );
    } else {
        for ( int i = 2; i < ps->narg; i++ ) {
            strcpy ( ps->hlpCtx [ i - 2 ], ps->sarg [ i ] );
        }
    }
}

void shellShowHelp ( shell *ps ) {
    int idx = 0;
    char hc [ sizeof ( ps->commands [ 0 ].hlpCtx ) ];
    char fmt[ 20 ];
    int cmdlen = shellMaxCmdLen( ps );
    int hctlen = shellMaxHctLen( ps );
    sprintf( fmt, "%%-%us %%-%us %%s\n", hctlen + 3, cmdlen + 2 );
    cprintf( "\n" );
    if ( ps->cmdCount > 0 ) {
        switch ( ps->narg ) {
        case 2 ... SHELL_MAX_ARGS :
            if ( !strcmp ( ps->sarg [ 0 ], "pause" ) ) {
                int i = MAX( ( int )ps->uarg [ 1 ], 0 );
                ps->hlpLinesPause = i;
                if ( i == 0 ) {
                    cprintf( "help pause off\n" );
                }
                else if ( i == 1 ) {
                    cprintf( "Help will pause every line\n" );
                }
                else {
                    cprintf( "Help will pause every %i lines\n", i );
                }
                return;
            } else if ( !strcmp( ps->sarg [ 0 ], "ctx" ) ) {
                SetHelpContext( ps );
                cprintf( "\nHelp context set to [ " );
                for ( int i = 0; i < SHELL_MAX_ARGS; i++ ) {
                    if ( ps->hlpCtx [ i ][ 0 ] != '\0' ) {
                        cprintf( "%s ", ps->hlpCtx [ i ] );
                    }
                    cprintf( "]\n" );
                }
            }
            return;
        case 1 :
            ( void ) strucase ( hc, ps->sarg [ 0 ] );
            break;
        case 0 :
            strcpy ( hc, "ALL" );
            break;
        default :
            return;
        }
        int hcidx     = 0;
        _Bool showAll = !strcmp ( hc, "ALL" );
        for ( ; idx < ps->cmdCount; idx++ ) {
            if ( hcidx > 0 && hcidx % ps->hlpLinesPause == 0 && pause ( "..." ) == CC_ESC ) {
                break;
            }
            else if ( showAll || !strcmp ( ps->commands [ idx ].hlpCtx, "ALL" )
                    || !strcmp( ps->commands [ idx ].hlpCtx, hc ) ) {
                if ( !hctcmp ( ps, ps->commands[ idx ].hlpCtx ) && ( showAll || !hctcmp ( ps, hc ) ) ) {
                    hcidx++;
                    if ( ps->commands [ idx ].callbackFunction == NULL ) {
                        cprintf( "\n" );
                        continue;
                    }
                    if ( strlen ( ps->commands [ idx ].hlpMsg ) != 0 ) {
                        char shc [ sizeof ( ps->commands [ idx ].hlpCtx ) ];
                        sprintf ( shc, "%s:", ps->commands [ idx ].hlpCtx );
                        cprintf( fmt, shc, ps->commands [ idx ].cmd, ps->commands [ idx ].hlpMsg );
                    }
                }
            }
        }
    }
    cprintf( "\n" );
}

void shellConsoleExit( shell *ps )
{
   ps->bConsoleExit = true;
}

shellError shellConsole( shell *ps, shellCommandDef *pdef, const char name [ ], const char startupCmd [ ], int helpLinesPause )
{
    shellError res;
    if (!ps->bIsInitialised ) {
        shellInitialise( ps, pdef, name, startupCmd, helpLinesPause );
        shellExecCommand( ps ); // execute startup command
    }
    ps->bConsoleExit = false;
    while(!ps->bConsoleExit)
    {
        shellGetCommand( ps );
        switch ( res = shellExecCommand( ps ) )
        {
            case SHELL_ERROR:
            {
                printf( "Not found\n" );
                break;
            }
            default:
            {
                break;
            }
        }
    }
    return ( res );
}

void shellEcho( shell *ps )
{
    char *si = ps->sarg[ 0 ];
    char so [SHELL_CMD_LINE_LEN + 1];
    _Bool eos = false;
    int j = 0;
    int i = 0;
    do {
        char c = si[ i++ ];
        switch ( c )
        {
        case '\\' :
            c = si [ i++ ];
            switch ( c )
            {
            case '\0' :
                so [ j++ ] = c;
                eos = true;
                break;
            case '\\' :
                so [ j++ ] = '\\';
                break;
            case 'b' :
                so [ j++ ] = '\b';
                break;
            case 'n' :
                so [ j++ ] = '\n';
                break;
            case '_' :
                so [ j++ ] = ' ';
                break;
            case ';' :
                so [ j++ ] = ';';
                break;
            default :
                so [ j++ ] = '\\';
                so [ j++ ] = c;
            }
            break;
        case '%'  :
            c = si [ i++ ];
            switch ( c )
            {
            case '\0' :
                so [ j++ ] = c;
                eos = true;
                break;
            case 'l' :
            case 'L' :
                if ( ps->nLoop )
                {
                    int n = sprintf( so + j, "%i", ps->loopCounter [ ps->nLoop - 1 ] + 1 );
                    j += ( n > 0 ? n : 0 );
                }
                break;
            default :
                so [ j++ ] = '%';
                so [ j++ ] = c;
            }
            break;
        case '\0' :
            so[ j ] = c;
            eos = true;
            break;
        default :
            so [ j++ ] = c;
            break;
        }
    } while ( !eos );
    cprintf( "%s", so );
}

void shellWait( shell *ps )
{
    double d = ps->darg[ 0 ];
#ifdef verbose
    cprintf( "\nWaiting %0.3f s\n", d );
#endif
    WaitSec( d );
}

void shellSetVerboseMode( shell *ps )
{
    verboseMode = ( ps->uarg[ 0 ] ? true : false );
    if ( verboseMode ) {
        cprintf( "Verbose mode enabled\n" );
    }
    else {
        cprintf( "Verbose mode inhibited\n" );
    }
}


//void shellMacro( shell *ps )
//{
//    u32 n = ps->uarg[ 0 ];
//    if ( n <= ps->nMacros )
//    {
//        char *pm = ps->macros[ n ];
//        char *pc = ps->cmdLine;
//        do
//        {
//            *pc++ = *pm;
//        }
//        while ( *pm++ != '\0' );
//    }
//}

