/*
 * shell.h
 *
 *  Created on: 22 Mar 2016
 *      Author: RTL3
 */

#ifndef SHELL_H_
#define SHELL_H_

#define SHELL_CMD_LEN           ( 30 )
#define SHELL_HLP_MSG_LEN 		( 100 )
#define SHELL_MAX_ARGS          ( 20 )
#define SHELL_MAX_ARG_LEN       ( 50 )
#define SHELL_NESTED_LOOPS      ( 10 )
#define SHELL_CMD_LINE_LEN      ( 200 )
#define SHELL_NAME_LEN          ( 30 )
#define SHELL_MAX_CMDS        	( 1000 )
#define SHELL_CMDLINE_BUF_SIZE  ( 50 )

typedef struct shellCommandDef shellCommandDef;
typedef struct shellCommand    shellCommand;
typedef struct shellArgs       shellArgs;
typedef struct shell           shell;

typedef enum shellError { SHELL_ERROR, SHELL_SUCCESS, SHELL_EXIT } shellError;

struct shellCommandDef {
    const char hlpCtx[ SHELL_MAX_ARG_LEN + 1 ];
    const char cmd    [ SHELL_CMD_LEN + 1 ];
    const int  narg;
    const char hlpMsg [ SHELL_HLP_MSG_LEN + 1 ];
    void  ( *callbackFunction )( shell *p );
};

struct shellArgs {
    u32          narg;
    u32          uarg[ SHELL_MAX_ARGS ];
    s32          iarg[ SHELL_MAX_ARGS ];
    double       darg[ SHELL_MAX_ARGS ];
    char         sarg[ SHELL_MAX_ARGS ][ SHELL_MAX_ARG_LEN + 1 ];
};

struct shellCommand {
    char         cmd[ SHELL_CMD_LEN + 1 ];
    char         hlpCtx[ SHELL_MAX_ARG_LEN + 1 ];
    char         hlpMsg[ SHELL_HLP_MSG_LEN + 1 ];
    int          narg;
    shellArgs    arg;
    void         (* callbackFunction )( shell *p );
};

struct shell {
    char         name[ SHELL_NAME_LEN + 1];
    int          hlpLinesPause;
    char         hlpCtx[ SHELL_MAX_ARGS ][ SHELL_MAX_ARG_LEN + 1 ];
    int          cmdCount;
    char         cmdLine [SHELL_CMD_LINE_LEN + 1];
    char         cmd     [SHELL_CMD_LINE_LEN + 1];
    int          cmdBufIdx;
    int          cmdBufLen;
    char         cmdLineBuf[SHELL_CMDLINE_BUF_SIZE][SHELL_CMD_LINE_LEN + 1];
    unsigned int narg;
    unsigned int uarg [ SHELL_MAX_ARGS ];
    signed int   iarg [ SHELL_MAX_ARGS ];
    double       darg [ SHELL_MAX_ARGS ];
    char         sarg [ SHELL_MAX_ARGS ][ SHELL_MAX_ARG_LEN + 1 ];
    char         *nxtCmdAddr;
    int          loopCounter[ SHELL_NESTED_LOOPS ];
    int          loopCount[ SHELL_NESTED_LOOPS ];
    char         *loopAddr[ SHELL_NESTED_LOOPS ];
    int          nLoop;
    shellArgs    *parg;
    shellCommand commands[ SHELL_MAX_CMDS ];
    _Bool        bConsoleExit;
    _Bool        bIsInitialised;
};


extern shellCommandDef MainCmdDef[];
extern shellCommandDef PrsBPACmdDef[];
extern shellCommandDef AD9364CmdDef[];

extern shell mainShell;
extern shell AD9364Shell;
extern shell PrsBPAShell;

extern shellError shellConsole ( shell *ps, shellCommandDef *pdef, const char name [ ], const char startupCmd [ ], int helpLinesPause );
extern void shellShowHelp      ( shell *ps );
extern void shellEcho          ( shell *ps );
extern void shellWait          ( shell *ps );
extern void shellLoop          ( shell *ps );
extern void shellConsoleExit   ( shell *ps );
extern void shellLoop          ( shell *ps );
extern void shellEndLoop       ( shell *ps );




#endif /* SHELL_H_ */



