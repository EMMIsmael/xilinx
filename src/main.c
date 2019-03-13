/************************************************************
 User Console for Telemetry System Debug
 ( c ) G. Mortimer 2017
 Relevant Technologies Ltd
 ************************************************************/

#include "include.h"


_Bool    escSeqDet   = true;

u32      axiBramCtrlrWrBuf[ AXI_BRAM_CTRLR_BRAM_SIZE ];
u32      axiBramCtrlrRdBuf[ AXI_BRAM_CTRLR_BRAM_SIZE ];

shell    mainShell;

TimerAPI *SysTimer;
u32      SysTimerId;


// ********************************************************************
// **                       Shell Functions                          **
// ********************************************************************
void mainArgs( u32 *arg0, u32 *arg1, u32 *arg2, u32 *arg3, shell *psh );

void mainShowBanner( shell *psh );
void mainToggleEscSeqDet( shell *psh );
void mainTestKeyboard( shell *psh );
void mainMemDump( shell *psh );
void mainTestArg( shell *psh );
void mainResetTimer( shell *psh );
void mainStartTimer( shell *psh );
void mainStopTimer( shell *psh );
void mainClearTimer( shell *psh );
void mainReadTimer( shell *psh );
void mainReadClearTimer( shell *psh );
void mainReadTimerLong( shell *psh );
void mainDumpTimerSlvReg( shell *psh );
void mainSetDacFreq( shell *psh );
void mainSysInit( shell *psh );
void mainSysReset( shell *psh );

//===========================
// Standard Device Functions
//===========================
//===========================================
//// Script-generated code - do not edit ////
//===========================================
//=============
// TxAXISource
//=============
void mainReadTxAXISourceBram( shell *psh );
void mainWriteTxAXISourceBram( shell *psh );
void mainInitTxAXISourceBram( shell *psh );
void mainLoadTxAXISourceBram( shell *psh );
void mainDumpTxAXISourceBram( shell *psh );

void mainReadTxAXISourceSlvReg( shell *psh );
void mainWriteTxAXISourceSlvReg( shell *psh );
void mainDumpTxAXISourceSlvReg( shell *psh );

void mainDumpTxAXISource( shell *psh );
void mainResetTxAXISource( shell *psh );

//================
// TxSourceBramC
//================
void mainReadTxSourceAXIBramCBram( shell *psh );
void mainWriteTxSourceAXIBramCBram( shell *psh );
void mainInitTxSourceAXIBramCBram( shell *psh );
void mainLoadTxSourceAXIBramCBram( shell *psh );

void mainReadTxSourceAXIBramCSlvReg( shell *psh );
void mainWriteTxSourceAXIBramCSlvReg( shell *psh );
void mainDumpTxSourceAXIBramCSlvReg( shell *psh );

void mainDumpTxSourceAXIBramC( shell *psh );
void mainResetTxSourceAXIBramC( shell *psh );

//========
// TxFFTC
//========
void mainReadTxFFTCSlvReg( shell *psh );
void mainWriteTxFFTCSlvReg( shell *psh );
void mainDumpTxFFTCSlvReg( shell *psh );

void mainDumpTxFFTC( shell *psh );
void mainResetTxFFTC( shell *psh );

//===========
// RxNsdFFTC
//===========
void mainReadRxNsdFFTCSlvReg( shell *psh );
void mainWriteRxNsdFFTCSlvReg( shell *psh );
void mainDumpRxNsdFFTCSlvReg( shell *psh );

void mainDumpRxNsdFFTC( shell *psh );
void mainResetRxNsdFFTC( shell *psh );

//=========
// TxTSInj
//=========
void mainReadTxTSInjTsBuf( shell *psh );
void mainWriteTxTSInjTsBuf( shell *psh );
void mainInitTxTSInjTsBuf( shell *psh );
void mainLoadTxTSInjTsBuf( shell *psh );
void mainDumpTxTSInjTsBuf( shell *psh );

void mainReadTxTSInjTxBuf( shell *psh );
void mainDumpTxTSInjTxBuf( shell *psh );

void mainWriteTxTSInjRxBuf( shell *psh );
void mainInitTxTSInjRxBuf( shell *psh );

void mainReadTxTSInjPsBuf( shell *psh );
void mainDumpTxTSInjPsBuf( shell *psh );

void mainReadTxTSInjSlvReg( shell *psh );
void mainWriteTxTSInjSlvReg( shell *psh );
void mainDumpTxTSInjSlvReg( shell *psh );

void mainDumpTxTSInj( shell *psh );
void mainResetTxTSInj( shell *psh );

//==============
// RFChanTREmul
//==============
void mainReadRFChanTREmulPsLUT( shell *psh );
void mainWriteRFChanTREmulPsLUT( shell *psh );
void mainInitRFChanTREmulPsLUT( shell *psh );
void mainLoadRFChanTREmulPsLUT( shell *psh );
void mainDumpRFChanTREmulPsLUT( shell *psh );

void mainReadRFChanTREmulHwLUT( shell *psh );
void mainWriteRFChanTREmulHwLUT( shell *psh );
void mainInitRFChanTREmulHwLUT( shell *psh );
void mainLoadRFChanTREmulHwLUT( shell *psh );
void mainDumpRFChanTREmulHwLUT( shell *psh );

void mainReadRFChanTREmulSlvReg( shell *psh );
void mainWriteRFChanTREmulSlvReg( shell *psh );
void mainDumpRFChanTREmulSlvReg( shell *psh );

void mainDumpRFChanTREmul( shell *psh );
void mainResetRFChanTREmul( shell *psh );

//================
// RFChanTREmulSw
//================
void mainReadRFChanTREmulSwSlvReg( shell *psh );
void mainWriteRFChanTREmulSwSlvReg( shell *psh );
void mainDumpRFChanTREmulSwSlvReg( shell *psh );

void mainDumpRFChanTREmulSw( shell *psh );
void mainResetRFChanTREmulSw( shell *psh );

//=====================
// RFChanTREmulFreqCtr
//=====================
void mainReadRFChanTREmulFreqCtrSlvReg( shell *psh );
void mainWriteRFChanTREmulFreqCtrSlvReg( shell *psh );
void mainDumpRFChanTREmulFreqCtrSlvReg( shell *psh );

void mainDumpRFChanTREmulFreqCtr( shell *psh );
void mainResetRFChanTREmulFreqCtr( shell *psh );
void mainReadRFChanTREmulFreqCtr( shell *psh );

//==================
// AD9364ClkFreqCtr
//==================
void mainReadAD9364ClkFreqCtrSlvReg( shell *psh );
void mainWriteAD9364ClkFreqCtrSlvReg( shell *psh );
void mainDumpAD9364ClkFreqCtrSlvReg( shell *psh );

void mainDumpAD9364ClkFreqCtr( shell *psh );
void mainResetAD9364ClkFreqCtr( shell *psh );
void mainReadAD9364ClkFreqCtr( shell *psh );

//==================
// AD9364DACFreqCtr
//==================
void mainReadAD9364DACFreqCtrSlvReg( shell *psh );
void mainWriteAD9364DACFreqCtrSlvReg( shell *psh );
void mainDumpAD9364DACFreqCtrSlvReg( shell *psh );

void mainDumpAD9364DACFreqCtr( shell *psh );
void mainResetAD9364DACFreqCtr( shell *psh );
void mainReadAD9364DACFreqCtr( shell *psh );

//==================
// AD9364ADCFreqCtr
//==================
void mainReadAD9364ADCFreqCtrSlvReg( shell *psh );
void mainWriteAD9364ADCFreqCtrSlvReg( shell *psh );
void mainDumpAD9364ADCFreqCtrSlvReg( shell *psh );

void mainDumpAD9364ADCFreqCtr( shell *psh );
void mainResetAD9364ADCFreqCtr( shell *psh );
void mainReadAD9364ADCFreqCtr( shell *psh );

//=============
// RxNSDet
//=============
void mainReadRxNSDetSlvReg( shell *psh );
void mainWriteRxNSDetSlvReg( shell *psh );
void mainDumpRxNSDetSlvReg( shell *psh );

void mainDumpRxNSDet( shell *psh );
void mainResetRxNSDet( shell *psh );

//==============
// RxTddAXISink
//==============
void mainReadRxTddAXISinkBram( shell *psh );
void mainWriteRxTddAXISinkBram( shell *psh );
void mainInitRxTddAXISinkBram( shell *psh );
void mainDumpRxTddAXISinkBram( shell *psh );

void mainReadRxTddAXISinkSlvReg( shell *psh );
void mainWriteRxTddAXISinkSlvReg( shell *psh );
void mainDumpRxTddAXISinkSlvReg( shell *psh );

void mainDumpRxTddAXISink( shell *psh );
void mainResetRxTddAXISink( shell *psh );

//==============
// RxPrsAXISink
//==============
void mainReadRxPrsAXISinkBram( shell *psh );
void mainWriteRxPrsAXISinkBram( shell *psh );
void mainInitRxPrsAXISinkBram( shell *psh );
void mainDumpRxPrsAXISinkBram( shell *psh );

void mainReadRxPrsAXISinkSlvReg( shell *psh );
void mainWriteRxPrsAXISinkSlvReg( shell *psh );
void mainDumpRxPrsAXISinkSlvReg( shell *psh );

void mainDumpRxPrsAXISink( shell *psh );
void mainResetRxPrsAXISink( shell *psh );

//==============
// RxDtaAXISink
//==============
void mainReadRxDtaAXISinkBram( shell *psh );
void mainWriteRxDtaAXISinkBram( shell *psh );
void mainInitRxDtaAXISinkBram( shell *psh );
void mainDumpRxDtaAXISinkBram( shell *psh );

void mainReadRxDtaAXISinkSlvReg( shell *psh );
void mainWriteRxDtaAXISinkSlvReg( shell *psh );
void mainDumpRxDtaAXISinkSlvReg( shell *psh );

void mainDumpRxDtaAXISink( shell *psh );
void mainResetRxDtaAXISink( shell *psh );

//=================
// RxNsdA00AXISink
//=================
void mainReadRxNsdA00AXISinkBram( shell *psh );
void mainWriteRxNsdA00AXISinkBram( shell *psh );
void mainInitRxNsdA00AXISinkBram( shell *psh );
void mainDumpRxNsdA00AXISinkBram( shell *psh );

void mainReadRxNsdA00AXISinkSlvReg( shell *psh );
void mainWriteRxNsdA00AXISinkSlvReg( shell *psh );
void mainDumpRxNsdA00AXISinkSlvReg( shell *psh );

void mainDumpRxNsdA00AXISink( shell *psh );
void mainResetRxNsdA00AXISink( shell *psh );

//=================
// RxNsdA01AXISink
//=================
void mainReadRxNsdA01AXISinkBram( shell *psh );
void mainWriteRxNsdA01AXISinkBram( shell *psh );
void mainInitRxNsdA01AXISinkBram( shell *psh );
void mainDumpRxNsdA01AXISinkBram( shell *psh );

void mainReadRxNsdA01AXISinkSlvReg( shell *psh );
void mainWriteRxNsdA01AXISinkSlvReg( shell *psh );
void mainDumpRxNsdA01AXISinkSlvReg( shell *psh );

void mainDumpRxNsdA01AXISink( shell *psh );
void mainResetRxNsdA01AXISink( shell *psh );

//=================
// RxNsdA15AXISink
//=================
void mainReadRxNsdA15AXISinkBram( shell *psh );
void mainWriteRxNsdA15AXISinkBram( shell *psh );
void mainInitRxNsdA15AXISinkBram( shell *psh );
void mainDumpRxNsdA15AXISinkBram( shell *psh );

void mainReadRxNsdA15AXISinkSlvReg( shell *psh );
void mainWriteRxNsdA15AXISinkSlvReg( shell *psh );
void mainDumpRxNsdA15AXISinkSlvReg( shell *psh );

void mainDumpRxNsdA15AXISink( shell *psh );
void mainResetRxNsdA15AXISink( shell *psh );

//====================
// RxNsdCordicAXISink
//====================
void mainReadRxNsdCordicAXISinkBram( shell *psh );
void mainWriteRxNsdCordicAXISinkBram( shell *psh );
void mainInitRxNsdCordicAXISinkBram( shell *psh );
void mainDumpRxNsdCordicAXISinkBram( shell *psh );

void mainReadRxNsdCordicAXISinkSlvReg( shell *psh );
void mainWriteRxNsdCordicAXISinkSlvReg( shell *psh );
void mainDumpRxNsdCordicAXISinkSlvReg( shell *psh );

void mainDumpRxNsdCordicAXISink( shell *psh );
void mainResetRxNsdCordicAXISink( shell *psh );

//=================
// RxNsdFFTAXISink
//=================
void mainReadRxNsdFFTAXISinkBram( shell *psh );
void mainWriteRxNsdFFTAXISinkBram( shell *psh );
void mainInitRxNsdFFTAXISinkBram( shell *psh );
void mainDumpRxNsdFFTAXISinkBram( shell *psh );

void mainReadRxNsdFFTAXISinkSlvReg( shell *psh );
void mainWriteRxNsdFFTAXISinkSlvReg( shell *psh );
void mainDumpRxNsdFFTAXISinkSlvReg( shell *psh );

void mainDumpRxNsdFFTAXISink( shell *psh );
void mainResetRxNsdFFTAXISink( shell *psh );

//====================
// RxNsdBcast00AXISink
//====================
void mainReadRxNsdB00AXISinkBram( shell *psh );
void mainWriteRxNsdB00AXISinkBram( shell *psh );
void mainInitRxNsdB00AXISinkBram( shell *psh );
void mainDumpRxNsdB00AXISinkBram( shell *psh );

void mainReadRxNsdB00AXISinkSlvReg( shell *psh );
void mainWriteRxNsdB00AXISinkSlvReg( shell *psh );
void mainDumpRxNsdB00AXISinkSlvReg( shell *psh );

void mainDumpRxNsdB00AXISink( shell *psh );
void mainResetRxNsdB00AXISink( shell *psh );

//====================
// RxNsdBcast01AXISink
//====================
void mainReadRxNsdB01AXISinkBram( shell *psh );
void mainWriteRxNsdB01AXISinkBram( shell *psh );
void mainInitRxNsdB01AXISinkBram( shell *psh );
void mainDumpRxNsdB01AXISinkBram( shell *psh );

void mainReadRxNsdB01AXISinkSlvReg( shell *psh );
void mainWriteRxNsdB01AXISinkSlvReg( shell *psh );
void mainDumpRxNsdB01AXISinkSlvReg( shell *psh );

void mainDumpRxNsdB01AXISink( shell *psh );
void mainResetRxNsdB01AXISink( shell *psh );

//====================
// RxNsdBcast15AXISink
//====================
void mainReadRxNsdB15AXISinkBram( shell *psh );
void mainWriteRxNsdB15AXISinkBram( shell *psh );
void mainInitRxNsdB15AXISinkBram( shell *psh );
void mainDumpRxNsdB15AXISinkBram( shell *psh );

void mainReadRxNsdB15AXISinkSlvReg( shell *psh );
void mainWriteRxNsdB15AXISinkSlvReg( shell *psh );
void mainDumpRxNsdB15AXISinkSlvReg( shell *psh );

void mainDumpRxNsdB15AXISink( shell *psh );
void mainResetRxNsdB15AXISink( shell *psh );

//=======
// PrsCIR
//=======
void mainReadRxPrsCIRBramFFT( shell *psh );
void mainWriteRxPrsCIRBramFFT( shell *psh );
void mainInitRxPrsCIRBramFFT( shell *psh );
void mainDumpRxPrsCIRBramFFT( shell *psh );

void mainReadRxPrsCIRBramABS( shell *psh );
void mainWriteRxPrsCIRBramABS( shell *psh );
void mainInitRxPrsCIRBramABS( shell *psh );
void mainDumpRxPrsCIRBramABS( shell *psh );

void mainReadRxPrsCIRSlvReg( shell *psh );
void mainWriteRxPrsCIRSlvReg( shell *psh );
void mainDumpRxPrsCIRSlvReg( shell *psh );

void mainDumpRxPrsCIR( shell *psh );
void mainResetRxPrsCIR( shell *psh );
void mainRxPrsCIRCalc( shell *psh );
void mainRxPrsCIRRecalc( shell *psh );
void mainRxPrsCIRDumpCIRTable( shell *psh );
void mainRxPrsCIRDumpTxPRSFd( shell *psh );
void mainRxPrsCIRDumpRxPRSTd( shell *psh );
void mainRxPrsCIRDumpRxPRSFd( shell *psh );
void mainRxPrsCIRDumpRxPRSTdCorDeltaT( shell *psh );
void mainRxPrsCIRDumpRxPRSTdCordeltafp( shell *psh );
void mainRxPrsCIRDumpRxPRSFdCordeltafp( shell *psh );
void mainRxPrsCIRDumpRxPRSFdCordeltatp( shell *psh );
void mainRxPrsCIRDumpRxPRSTdCorFT( shell *psh );


//=======
// PrsBPA
//=======
void mainReadRxPrsBPASlvReg( shell *psh );
void mainWriteRxPrsBPASlvReg( shell *psh );
void mainDumpRxPrsBPASlvReg( shell *psh );

void mainDumpRxPrsBPA( shell *psh );
void mainResetRxPrsBPA( shell *psh );
void mainRxPrsBPACalcDelta( shell *psh );
void mainRxPrsBPACalcDelta2( shell *psh );

void mainRxPrsBPADumpCorSymbFd( shell *psh );
void mainRxPrsBPADumpCorSymbTd( shell *psh );

void mainRxDataBPADumpCorSymbFd( shell *psh );
void mainRxDataBPADumpCorSymbTd( shell *psh );


void mainRxDataSymbolDumpRawFd( shell *psh );
void mainRxDataSymbolDumpFCFd( shell *psh );
void mainRxDataSymbolDumpFTCFd( shell *psh );


//===========================================
//////// End of script-generated code ///////
//===========================================

//*********************
// User Shell Functions
//*********************
void mainAD9364Shell( shell *psh );
void mainPrsBPAShell( shell *psh );

void mainInitAD9364( shell *psh );
void mainSelChRFChanEmulSw( shell *psh );
void mainTestRxNSDet( shell *psh );
void mainSetTmodeRxNSDet( shell *psh );
void mainStartRxNSDet( shell *psh );
void mainStopTmodeRxNSDet( shell *psh );
void mainDstateRxNSDet( shell *psh );
void mainSetSysParams( shell *psh );
void mainRunTest( shell *psh );
void mainRunTestBPA( shell *psh );
void mainStopTest( shell * psh );
void mainChkTSITxFreq( shell * psh );
void mainChkNSDRxFreq( shell * psh );
// Shell functions end

// Command format is:
// [ help context ] [ command ] [ n args ] [ help message ] [ callback function ]
// [ n args ] = number of arguments.
//              These are collected as unsigned integers (shell.uargs[]), doubles(shell.dargs[]) and strings(shell.sargs[]).
//              however, it is possible to input them as signed integers, and recast in the callback function
//              For floating point input, use the string argument and convert in the callback function
//
// For different contextx, with different numbers of arguments, register a command for each.
//
// [ help message ] = [ arg description        explanation ]  to display
//                  | [ "" ] to hide display of different command contexts
//
// Example: Help, two commands, 3 different contexts
//    ( void ) shellRegisterCommand( &mainShell, "all"  , "help"      ,  0, ""                                                                           , mainShowHelp                            );
//    ( void ) shellRegisterCommand( &mainShell, "all"  , "help"      ,  1, ""                                                                           , mainShowHelp                            );
//    ( void ) shellRegisterCommand( &mainShell, "all"  , "help"      ,  2, ""                                                                           , mainShowHelp                            );
//    ( void ) shellRegisterCommand( &mainShell, "all"  , "h"         ,  0, "                                        Show help"                          , mainShowHelp                            );
//    ( void ) shellRegisterCommand( &mainShell, "all"  , "h"         ,  1, "<context>                               Show help for given context"        , mainShowHelp                            );
//    ( void ) shellRegisterCommand( &mainShell, "all"  , "h"         ,  2, "pause n                                 Set help pause lines"               , mainShowHelp                            );
//    Shell command interpreter is NOT case-sensitive, for either help context or command
//    In case a command can be accepted with no arguments, but needs them, the last values given will be reused
//    A command line buffer is kept: use up/down arrows to access
//    - this is searched on the basis of the first characters typed

shellCommandDef mainCmdDef[ ] = {
    { "all"  , "b"         ,  0, "                                        Show banner"                        , mainShowBanner                          },
    { "all"  , "e"         ,  0, "                                        Toggle ESC sequence detection"      , mainToggleEscSeqDet                     },
    { "all"  , "k"         ,  0, "                                        Test Keyboard"                      , mainTestKeyboard                        },
    { "all"  , ""          ,  0, "",                                                                                                               NULL },
    { "all"  , "ta"        ,  1, "<arg>                                   Test argument conversion"           , mainTestArg                             },
    { "all"  , ""          ,  0, "",                                                                                                               NULL },
    { "all"  , "d"         ,  0, ""                                                                           , mainMemDump                             },
    { "all"  , "d"         ,  1, ""                                                                           , mainMemDump                             },
    { "all"  , "d"         ,  2, ""                                                                           , mainMemDump                             },
    { "all"  , "d"         ,  3, ""                                                                           , mainMemDump                             },
    { "all"  , "d"         ,  4, "<addr>[ <fmt> ][ arg1 ][ arg2 ]         Memory Dump"                        , mainMemDump                             },
    { "all"  , ""          ,  0, "",                                                                                                               NULL },
    { "all"  , "rsttmr"    ,  0, "                                        Reset timer"                        , mainResetTimer                          },
    { "all"  , "stmr"      ,  0, "                                        Start timer"                        , mainStartTimer                          },
    { "all"  , "ptmr"      ,  0, "                                        Stop timer"                         , mainStopTimer                           },
    { "all"  , "ctmr"      ,  0, "                                        Clear timer"                        , mainClearTimer                          },
    { "all"  , "rtmr"      ,  0, "                                        Read timer"                         , mainReadTimer                           },
    { "all"  , "rctmr"     ,  0, "                                        Read/clear timer"                   , mainReadClearTimer                      },
    { "all"  , "rtmrl"     ,  0, "                                        Read timer long"                    , mainReadTimerLong                       },
    { "all"  , "dtmr"      ,  0, ""                                                                           , mainDumpTimerSlvReg                     },
    { "all"  , "dtmr"      ,  1, "[<fmt(0=hex 1=bin)>]                    Dump Timer Slave Regs"              , mainDumpTimerSlvReg                     },
    { "all"  , ""          ,  0, "",                                                                                                               NULL },
//===========================================
//// Script-generated code - do not edit ////
//===========================================
    { "AXS"  , "AxsRBram"  ,  1, "<adr>                                   Read  TxAXISource Bram"                , mainReadTxAXISourceBram              },
    { "AXS"  , "AxsWBram"  ,  2, "<adr val>                               Write TxAXISource Bram"                , mainWriteTxAXISourceBram             },
    { "AXS"  , "AxsIBram"  ,  3, "<adr nwd val>                           Init  TxAXISource Bram"                , mainInitTxAXISourceBram              },
    { "AXS"  , "AxsLBram"  ,  0, "                                        Load  TxAXISource Bram"                , mainLoadTxAXISourceBram              },
    { "AXS"  , "AxsDBram"  ,  0, ""                                                                              , mainDumpTxAXISourceBram              },
    { "AXS"  , "AxsDBram"  ,  1, ""                                                                              , mainDumpTxAXISourceBram              },
    { "AXS"  , "AxsDBram"  ,  2, "[<adr nbit(0, 8-16)>]                   Dump  TxAXISource Bram"                , mainDumpTxAXISourceBram              },
    { "AXS"  , "AxsRslv"   ,  1, "<reg>                                   Read  TxAXISource Slv Reg"             , mainReadTxAXISourceSlvReg            },
    { "AXS"  , "AxsWslv"   ,  2, "<reg val>                               Write TxAXISource Slv Reg"             , mainWriteTxAXISourceSlvReg           },
    { "AXS"  , "AxsDslv"   ,  0, ""                                                                              , mainDumpTxAXISourceSlvReg            },
    { "AXS"  , "AxsDslv"   ,  1, "[<fmt(0=hex 1=bin)>]                    Dump  TxAXISource Slv Reg"             , mainDumpTxAXISourceSlvReg            },
    { "AXS"  , "AxsDmp"    ,  0, "                                        Dump  TxAXISource State"               , mainDumpTxAXISource                  },
    { "AXS"  , "AxsRst"    ,  0, "                                        Reset TxAXISource"                     , mainResetTxAXISource                 },
    { "AXS"  , ""          ,  0, "",                                                                                                               NULL },
    { "BRC"  , "BrcRBram"  ,  1, "<adr>                                   Read  TxSourceAXIBramC Bram"           , mainReadTxSourceAXIBramCBram         },
    { "BRC"  , "BrcWBram"  ,  2, "<adr val>                               Write TxSourceAXIBramC Bram"           , mainWriteTxSourceAXIBramCBram        },
    { "BRC"  , "BrcIBram"  ,  3, "<adr nwd val>                           Init  TxSourceAXIBramC Bram"           , mainInitTxSourceAXIBramCBram         },
    { "BRC"  , "BrcLBram"  ,  0, "                                        Load  TxSourceAXIBramC Bram"           , mainLoadTxSourceAXIBramCBram         },
    { "BRC"  , "BrcRslv"   ,  1, "<reg>                                   Read  TxSourceAXIBramC Slv Reg"        , mainReadTxSourceAXIBramCSlvReg       },
    { "BRC"  , "BrcWslv"   ,  2, "<reg val>                               Write TxSourceAXIBramC Slv Reg"        , mainWriteTxSourceAXIBramCSlvReg      },
    { "BRC"  , "BrcDslv"   ,  0, ""                                                                              , mainDumpTxSourceAXIBramCSlvReg       },
    { "BRC"  , "BrcDslv"   ,  1, "[<fmt(0=hex 1=bin)>]                    Dump  TxSourceAXIBramC Slv Reg"        , mainDumpTxSourceAXIBramCSlvReg       },
    { "BRC"  , "BrcDmp"    ,  0, "                                        Dump  TxSourceAXIBramC State"          , mainDumpTxSourceAXIBramC             },
    { "BRC"  , "BrcRst"    ,  0, "                                        Reset TxSourceAXIBramC"                , mainResetTxSourceAXIBramC            },
    { "BRC"  , ""          ,  0, "",                                                                                                               NULL },
    { "TXF"  , "TxfRslv"   ,  1, "<reg>                                   Read  TxFFTC Slv Reg"                  , mainReadTxFFTCSlvReg                 },
    { "TXF"  , "TxfWslv"   ,  2, "<reg val>                               Write TxFFTC Slv Reg"                  , mainWriteTxFFTCSlvReg                },
    { "TXF"  , "TxfDslv"   ,  0, ""                                                                              , mainDumpTxFFTCSlvReg                 },
    { "TXF"  , "TxfDslv"   ,  1, "[<fmt(0=hex 1=bin)>]                    Dump  TxFFTC Slv Reg"                  , mainDumpTxFFTCSlvReg                 },
    { "TXF"  , "TxfDmp"    ,  0, "                                        Dump  TxFFTC State"                    , mainDumpTxFFTC                       },
    { "TXF"  , "TxfRst"    ,  0, "                                        Reset TxFFTC"                          , mainResetTxFFTC                      },
    { "TXF"  , ""          ,  0, "",                                                                                                               NULL },
    { "NSF"  , "NsfRslv"   ,  1, "<reg>                                   Read  RxNsdFFTC Slv Reg"               , mainReadRxNsdFFTCSlvReg              },
    { "NSF"  , "NsfWslv"   ,  2, "<reg val>                               Write RxNsdFFTC Slv Reg"               , mainWriteRxNsdFFTCSlvReg             },
    { "NSF"  , "NsfDslv"   ,  0, ""                                                                              , mainDumpRxNsdFFTCSlvReg              },
    { "NSF"  , "NsfDslv"   ,  1, "[<fmt(0=hex 1=bin)>]                    Dump  RxNsdFFTC Slv Reg"               , mainDumpRxNsdFFTCSlvReg              },
    { "NSF"  , "NsfDmp"    ,  0, "                                        Dump  RxNsdFFTC State"                 , mainDumpRxNsdFFTC                    },
    { "NSF"  , "NsfRst"    ,  0, "                                        Reset RxNsdFFTC"                       , mainResetRxNsdFFTC                   },
    { "NSF"  , ""          ,  0, "",                                                                                                               NULL },
    { "TSI"  , "TsiRTsBuf" ,  1, "<adr>                                   Read  TxTSInj TsBuf"                   , mainReadTxTSInjTsBuf                 },
    { "TSI"  , "TsiWTsBuf" ,  2, "<adr val>                               Write TxTSInj TsBuf"                   , mainWriteTxTSInjTsBuf                },
    { "TSI"  , "TsiITsBuf" ,  3, "<adr nwd val>                           Init  TxTSInj TsBuf"                   , mainInitTxTSInjTsBuf                 },
    { "TSI"  , "TsiLTsBuf" ,  0, "                                        Load  TxTSInj TsBuf"                   , mainLoadTxTSInjTsBuf                 },
    { "TSI"  , "TsiDTsBuf" ,  0, ""                                                                              , mainDumpTxTSInjTsBuf                 },
    { "TSI"  , "TsiDTsBuf" ,  1, ""                                                                              , mainDumpTxTSInjTsBuf                 },
    { "TSI"  , "TsiDTsBuf" ,  2, "[<adr nbit(0, 8-16)>]                   Dump  TxTSInj TsBuf"                   , mainDumpTxTSInjTsBuf                 },
    { "TSI"  , "TsiRTxBuf" ,  1, "<adr>                                   Read  TxTSInj TxBuf"                   , mainReadTxTSInjTxBuf                 },
    { "TSI"  , "TsiDTxBuf" ,  0, ""                                                                              , mainDumpTxTSInjTxBuf                 },
    { "TSI"  , "TsiDTxBuf" ,  1, ""                                                                              , mainDumpTxTSInjTxBuf                 },
    { "TSI"  , "TsiDTxBuf" ,  2, "[<adr nbit(0, 8-16)>]                   Dump  TxTSInj TxBuf"                   , mainDumpTxTSInjTxBuf                 },
    { "TSI"  , "TsiWRxBuf" ,  2, "<adr val>                               Write TxTSInj RxBuf"                   , mainWriteTxTSInjRxBuf                },
    { "TSI"  , "TsiIRxBuf" ,  3, "<adr nwd val>                           Init  TxTSInj RxBuf"                   , mainInitTxTSInjRxBuf                 },
    { "TSI"  , "TsiRPsBuf" ,  1, "<adr>                                   Read  TxTSInj PsBuf"                   , mainReadTxTSInjPsBuf                 },
    { "TSI"  , "TsiDPsBuf" ,  0, ""                                                                              , mainDumpTxTSInjPsBuf                 },
    { "TSI"  , "TsiDPsBuf" ,  1, ""                                                                              , mainDumpTxTSInjPsBuf                 },
    { "TSI"  , "TsiDPsBuf" ,  2, "[<adr nbit(0, 8-16)>]                   Dump  TxTSInj PsBuf"                   , mainDumpTxTSInjPsBuf                 },
    { "TSI"  , "TsiRslv"   ,  1, "<reg>                                   Read  TxTSInj Slv Reg"                 , mainReadTxTSInjSlvReg                },
    { "TSI"  , "TsiWslv"   ,  2, "<reg val>                               Write TxTSInj Slv Reg"                 , mainWriteTxTSInjSlvReg               },
    { "TSI"  , "TsiDslv"   ,  0, ""                                                                              , mainDumpTxTSInjSlvReg                },
    { "TSI"  , "TsiDslv"   ,  1, "[<fmt(0=hex 1=bin)>]                    Dump  TxTSInj Slv Reg"                 , mainDumpTxTSInjSlvReg                },
    { "TSI"  , "TsiDmp"    ,  0, "                                        Dump  TxTSInj State"                   , mainDumpTxTSInj                      },
    { "TSI"  , "TsiRst"    ,  0, "                                        Reset TxTSInj"                         , mainResetTxTSInj                     },
    { "TSI"  , ""          ,  0, "",                                                                                                               NULL },
    { "EMU"  , "EmuRPsLUT" ,  1, "<adr>                                   Read  RFChanTREmul PsLUT"              , mainReadRFChanTREmulPsLUT            },
    { "EMU"  , "EmuWPsLUT" ,  2, "<adr val>                               Write RFChanTREmul PsLUT"              , mainWriteRFChanTREmulPsLUT           },
    { "EMU"  , "EmuIPsLUT" ,  3, "<adr nwd val>                           Init  RFChanTREmul PsLUT"              , mainInitRFChanTREmulPsLUT            },
    { "EMU"  , "EmuLPsLUT" ,  0, "                                        Load  RFChanTREmul PsLUT"              , mainLoadRFChanTREmulPsLUT            },
    { "EMU"  , "EmuDPsLUT" ,  0, ""                                                                              , mainDumpRFChanTREmulPsLUT            },
    { "EMU"  , "EmuDPsLUT" ,  1, ""                                                                              , mainDumpRFChanTREmulPsLUT            },
    { "EMU"  , "EmuDPsLUT" ,  2, "[<adr nbit(0, 8-16)>]                   Dump  RFChanTREmul PsLUT"              , mainDumpRFChanTREmulPsLUT            },
    { "EMU"  , "EmuRHwLUT" ,  1, "<adr>                                   Read  RFChanTREmul HwLUT"              , mainReadRFChanTREmulHwLUT            },
    { "EMU"  , "EmuWHwLUT" ,  2, "<adr val>                               Write RFChanTREmul HwLUT"              , mainWriteRFChanTREmulHwLUT           },
    { "EMU"  , "EmuIHwLUT" ,  3, "<adr nwd val>                           Init  RFChanTREmul HwLUT"              , mainInitRFChanTREmulHwLUT            },
    { "EMU"  , "EmuLHwLUT" ,  0, "                                        Load  RFChanTREmul HwLUT"              , mainLoadRFChanTREmulHwLUT            },
    { "EMU"  , "EmuDHwLUT" ,  0, ""                                                                              , mainDumpRFChanTREmulHwLUT            },
    { "EMU"  , "EmuDHwLUT" ,  1, ""                                                                              , mainDumpRFChanTREmulHwLUT            },
    { "EMU"  , "EmuDHwLUT" ,  2, "[<adr nbit(0, 8-16)>]                   Dump  RFChanTREmul HwLUT"              , mainDumpRFChanTREmulHwLUT            },
    { "EMU"  , "EmuRslv"   ,  1, "<reg>                                   Read  RFChanTREmul Slv Reg"            , mainReadRFChanTREmulSlvReg           },
    { "EMU"  , "EmuWslv"   ,  2, "<reg val>                               Write RFChanTREmul Slv Reg"            , mainWriteRFChanTREmulSlvReg          },
    { "EMU"  , "EmuDslv"   ,  0, ""                                                                              , mainDumpRFChanTREmulSlvReg           },
    { "EMU"  , "EmuDslv"   ,  1, "[<fmt(0=hex 1=bin)>]                    Dump  RFChanTREmul Slv Reg"            , mainDumpRFChanTREmulSlvReg           },
    { "EMU"  , "EmuDmp"    ,  0, "                                        Dump  RFChanTREmul State"              , mainDumpRFChanTREmul                 },
    { "EMU"  , "EmuRst"    ,  0, "                                        Reset RFChanTREmul"                    , mainResetRFChanTREmul                },
    { "EMU"  , ""          ,  0, "",                                                                                                               NULL },
    { "EMS"  , "EmsRslv"   ,  1, "<reg>                                   Read  RFChanTREmulSw Slv Reg"          , mainReadRFChanTREmulSwSlvReg         },
    { "EMS"  , "EmsWslv"   ,  2, "<reg val>                               Write RFChanTREmulSw Slv Reg"          , mainWriteRFChanTREmulSwSlvReg        },
    { "EMS"  , "EmsDslv"   ,  0, ""                                                                              , mainDumpRFChanTREmulSwSlvReg         },
    { "EMS"  , "EmsDslv"   ,  1, "[<fmt(0=hex 1=bin)>]                    Dump  RFChanTREmulSw Slv Reg"          , mainDumpRFChanTREmulSwSlvReg         },
    { "EMS"  , "EmsDmp"    ,  0, "                                        Dump  RFChanTREmulSw State"            , mainDumpRFChanTREmulSw               },
    { "EMS"  , "EmsRst"    ,  0, "                                        Reset RFChanTREmulSw"                  , mainResetRFChanTREmulSw              },
    { "EMS"  , ""          ,  0, "",                                                                                                               NULL },
    { "EMF"  , "EmfRslv"   ,  1, "<reg>                                   Read  RFChanTREmulFreqCtr Slv Reg"     , mainReadRFChanTREmulFreqCtrSlvReg    },
    { "EMF"  , "EmfWslv"   ,  2, "<reg val>                               Write RFChanTREmulFreqCtr Slv Reg"     , mainWriteRFChanTREmulFreqCtrSlvReg   },
    { "EMF"  , "EmfDslv"   ,  0, ""                                                                              , mainDumpRFChanTREmulFreqCtrSlvReg    },
    { "EMF"  , "EmfDslv"   ,  1, "[<fmt(0=hex 1=bin)>]                    Dump  RFChanTREmulFreqCtr Slv Reg"     , mainDumpRFChanTREmulFreqCtrSlvReg    },
    { "EMF"  , "EmfDmp"    ,  0, "                                        Dump  RFChanTREmulFreqCtr State"       , mainDumpRFChanTREmulFreqCtr          },
    { "EMF"  , "EmfRst"    ,  0, "                                        Reset RFChanTREmulFreqCtr"             , mainResetRFChanTREmulFreqCtr         },
    { "EMF"  , "EmfRead"   ,  0, "                                        Read  RFChanTREmulFreqCtr"             , mainReadRFChanTREmulFreqCtr          },
    { "EMF"  , ""          ,  0, "",                                                                                                               NULL },
    { "ACF"  , "AcfRslv"   ,  1, "<reg>                                   Read  AD9364ClkFreqCtr Slv Reg"        , mainReadAD9364ClkFreqCtrSlvReg       },
    { "ACF"  , "AcfWslv"   ,  2, "<reg val>                               Write AD9364ClkFreqCtr Slv Reg"        , mainWriteAD9364ClkFreqCtrSlvReg      },
    { "ACF"  , "AcfDslv"   ,  0, ""                                                                              , mainDumpAD9364ClkFreqCtrSlvReg       },
    { "ACF"  , "AcfDslv"   ,  1, "[<fmt(0=hex 1=bin)>]                    Dump  AD9364ClkFreqCtr Slv Reg"        , mainDumpAD9364ClkFreqCtrSlvReg       },
    { "ACF"  , "AcfDmp"    ,  0, "                                        Dump  AD9364ClkFreqCtr State"          , mainDumpAD9364ClkFreqCtr             },
    { "ACF"  , "AcfRst"    ,  0, "                                        Reset AD9364ClkFreqCtr"                , mainResetAD9364ClkFreqCtr            },
    { "ACF"  , "AcfRead"   ,  0, "                                        Read  AD9364ClkFreqCtr"                , mainReadAD9364ClkFreqCtr             },
    { "ACF"  , ""          ,  0, "",                                                                                                               NULL },
    { "ADF"  , "AdfRslv"   ,  1, "<reg>                                   Read  AD9364DACFreqCtr Slv Reg"        , mainReadAD9364DACFreqCtrSlvReg       },
    { "ADF"  , "AdfWslv"   ,  2, "<reg val>                               Write AD9364DACFreqCtr Slv Reg"        , mainWriteAD9364DACFreqCtrSlvReg      },
    { "ADF"  , "AdfDslv"   ,  0, ""                                                                              , mainDumpAD9364DACFreqCtrSlvReg       },
    { "ADF"  , "AdfDslv"   ,  1, "[<fmt(0=hex 1=bin)>]                    Dump  AD9364DACFreqCtr Slv Reg"        , mainDumpAD9364DACFreqCtrSlvReg       },
    { "ADF"  , "AdfDmp"    ,  0, "                                        Dump  AD9364DACFreqCtr State"          , mainDumpAD9364DACFreqCtr             },
    { "ADF"  , "AdfRst"    ,  0, "                                        Reset AD9364DACFreqCtr"                , mainResetAD9364DACFreqCtr            },
    { "ADF"  , "AdfRead"   ,  0, "                                        Read  AD9364DACFreqCtr"                , mainReadAD9364DACFreqCtr             },
    { "ADF"  , ""          ,  0, "",                                                                                                               NULL },
    { "AAF"  , "AafRslv"   ,  1, "<reg>                                   Read  AD9364ADCFreqCtr Slv Reg"        , mainReadAD9364ADCFreqCtrSlvReg       },
    { "AAF"  , "AafWslv"   ,  2, "<reg val>                               Write AD9364ADCFreqCtr Slv Reg"        , mainWriteAD9364ADCFreqCtrSlvReg      },
    { "AAF"  , "AafDslv"   ,  0, ""                                                                              , mainDumpAD9364ADCFreqCtrSlvReg       },
    { "AAF"  , "AafDslv"   ,  1, "[<fmt(0=hex 1=bin)>]                    Dump  AD9364ADCFreqCtr Slv Reg"        , mainDumpAD9364ADCFreqCtrSlvReg       },
    { "AAF"  , "AafDmp"    ,  0, "                                        Dump  AD9364ADCFreqCtr State"          , mainDumpAD9364ADCFreqCtr             },
    { "AAF"  , "AafRst"    ,  0, "                                        Reset AD9364ADCFreqCtr"                , mainResetAD9364ADCFreqCtr            },
    { "AAF"  , "AafRead"   ,  0, "                                        Read  AD9364ADCFreqCtr"                , mainReadAD9364ADCFreqCtr             },
    { "AAF"  , ""          ,  0, "",                                                                                                               NULL },
    { "NSD"  , "NsdRslv"   ,  1, "<reg>                                   Read  RxNSDet Slv Reg"                 , mainReadRxNSDetSlvReg                },
    { "NSD"  , "NsdWslv"   ,  2, "<reg val>                               Write RxNSDet Slv Reg"                 , mainWriteRxNSDetSlvReg               },
    { "NSD"  , "NsdDslv"   ,  0, ""                                                                              , mainDumpRxNSDetSlvReg                },
    { "NSD"  , "NsdDslv"   ,  1, "[<fmt(0=hex 1=bin)>]                    Dump  RxNSDet Slv Reg"                 , mainDumpRxNSDetSlvReg                },
    { "NSD"  , "NsdDmp"    ,  0, "                                        Dump  RxNSDet State"                   , mainDumpRxNSDet                      },
    { "NSD"  , "NsdRst"    ,  0, "                                        Reset RxNSDet"                         , mainResetRxNSDet                     },
    { "NSD"  , ""          ,  0, "",                                                                                                               NULL },
    { "DSK"  , "DskRBram"  ,  1, "<adr>                                   Read  RxTddAXISink Bram"               , mainReadRxTddAXISinkBram             },
    { "DSK"  , "DskWBram"  ,  2, "<adr val>                               Write RxTddAXISink Bram"               , mainWriteRxTddAXISinkBram            },
    { "DSK"  , "DskIBram"  ,  3, "<adr nwd val>                           Init  RxTddAXISink Bram"               , mainInitRxTddAXISinkBram             },
    { "DSK"  , "DskDBram"  ,  0, ""                                                                              , mainDumpRxTddAXISinkBram             },
    { "DSK"  , "DskDBram"  ,  1, ""                                                                              , mainDumpRxTddAXISinkBram             },
    { "DSK"  , "DskDBram"  ,  2, "[<adr nbit(0, 8-16)>]                   Dump  RxTddAXISink Bram"               , mainDumpRxTddAXISinkBram             },
    { "DSK"  , "DskRslv"   ,  1, "<reg>                                   Read  RxTddAXISink Slv Reg"            , mainReadRxTddAXISinkSlvReg           },
    { "DSK"  , "DskWslv"   ,  2, "<reg val>                               Write RxTddAXISink Slv Reg"            , mainWriteRxTddAXISinkSlvReg          },
    { "DSK"  , "DskDslv"   ,  0, ""                                                                              , mainDumpRxTddAXISinkSlvReg           },
    { "DSK"  , "DskDslv"   ,  1, "[<fmt(0=hex 1=bin)>]                    Dump  RxTddAXISink Slv Reg"            , mainDumpRxTddAXISinkSlvReg           },
    { "DSK"  , "DskDmp"    ,  0, "                                        Dump  RxTddAXISink State"              , mainDumpRxTddAXISink                 },
    { "DSK"  , "DskRst"    ,  0, "                                        Reset RxTddAXISink"                    , mainResetRxTddAXISink                },
    { "DSK"  , ""          ,  0, "",                                                                                                               NULL },
    { "PRS"  , "PrsRBram"  ,  1, "<adr>                                   Read  RxPrsAXISink Bram"               , mainReadRxPrsAXISinkBram             },
    { "PRS"  , "PrsWBram"  ,  2, "<adr val>                               Write RxPrsAXISink Bram"               , mainWriteRxPrsAXISinkBram            },
    { "PRS"  , "PrsIBram"  ,  3, "<adr nwd val>                           Init  RxPrsAXISink Bram"               , mainInitRxPrsAXISinkBram             },
    { "PRS"  , "PrsDBram"  ,  0, ""                                                                              , mainDumpRxPrsAXISinkBram             },
    { "PRS"  , "PrsDBram"  ,  1, ""                                                                              , mainDumpRxPrsAXISinkBram             },
    { "PRS"  , "PrsDBram"  ,  2, "[<adr nbit(0, 8-16)>]                   Dump  RxPrsAXISink Bram"               , mainDumpRxPrsAXISinkBram             },
    { "PRS"  , "PrsRslv"   ,  1, "<reg>                                   Read  RxPrsAXISink Slv Reg"            , mainReadRxPrsAXISinkSlvReg           },
    { "PRS"  , "PrsWslv"   ,  2, "<reg val>                               Write RxPrsAXISink Slv Reg"            , mainWriteRxPrsAXISinkSlvReg          },
    { "PRS"  , "PrsDslv"   ,  0, ""                                                                              , mainDumpRxPrsAXISinkSlvReg           },
    { "PRS"  , "PrsDslv"   ,  1, "[<fmt(0=hex 1=bin)>]                    Dump  RxPrsAXISink Slv Reg"            , mainDumpRxPrsAXISinkSlvReg           },
    { "PRS"  , "PrsDmp"    ,  0, "                                        Dump  RxPrsAXISink State"              , mainDumpRxPrsAXISink                 },
    { "PRS"  , "PrsRst"    ,  0, "                                        Reset RxPrsAXISink"                    , mainResetRxPrsAXISink                },
    { "PRS"  , ""          ,  0, "",                                                                                                               NULL },
    { "DTA"  , "DtaRBram"  ,  1, "<adr>                                   Read  RxDtaAXISink Bram"               , mainReadRxDtaAXISinkBram             },
    { "DTA"  , "DtaWBram"  ,  2, "<adr val>                               Write RxDtaAXISink Bram"               , mainWriteRxDtaAXISinkBram            },
    { "DTA"  , "DtaIBram"  ,  3, "<adr nwd val>                           Init  RxDtaAXISink Bram"               , mainInitRxDtaAXISinkBram             },
    { "DTA"  , "DtaDBram"  ,  0, ""                                                                              , mainDumpRxDtaAXISinkBram             },
    { "DTA"  , "DtaDBram"  ,  1, ""                                                                              , mainDumpRxDtaAXISinkBram             },
    { "DTA"  , "DtaDBram"  ,  2, "[<adr nbit(0, 8-16)>]                   Dump  RxDtaAXISink Bram"               , mainDumpRxDtaAXISinkBram             },
    { "DTA"  , "DtaRslv"   ,  1, "<reg>                                   Read  RxDtaAXISink Slv Reg"            , mainReadRxDtaAXISinkSlvReg           },
    { "DTA"  , "DtaWslv"   ,  2, "<reg val>                               Write RxDtaAXISink Slv Reg"            , mainWriteRxDtaAXISinkSlvReg          },
    { "DTA"  , "DtaDslv"   ,  0, ""                                                                              , mainDumpRxDtaAXISinkSlvReg           },
    { "DTA"  , "DtaDslv"   ,  1, "[<fmt(0=hex 1=bin)>]                    Dump  RxDtaAXISink Slv Reg"            , mainDumpRxDtaAXISinkSlvReg           },
    { "DTA"  , "DtaDmp"    ,  0, "                                        Dump  RxDtaAXISink State"              , mainDumpRxDtaAXISink                 },
    { "DTA"  , "DtaRst"    ,  0, "                                        Reset RxDtaAXISink"                    , mainResetRxDtaAXISink                },
    { "DTA"  , ""          ,  0, "",                                                                                                               NULL },
    { "A00"  , "A00RBram"  ,  1, "<adr>                                   Read  RxNsdA00AXISink Bram"            , mainReadRxNsdA00AXISinkBram          },
    { "A00"  , "A00WBram"  ,  2, "<adr val>                               Write RxNsdA00AXISink Bram"            , mainWriteRxNsdA00AXISinkBram         },
    { "A00"  , "A00IBram"  ,  3, "<adr nwd val>                           Init  RxNsdA00AXISink Bram"            , mainInitRxNsdA00AXISinkBram          },
    { "A00"  , "A00DBram"  ,  0, ""                                                                              , mainDumpRxNsdA00AXISinkBram          },
    { "A00"  , "A00DBram"  ,  1, ""                                                                              , mainDumpRxNsdA00AXISinkBram          },
    { "A00"  , "A00DBram"  ,  2, "[<adr nbit(0, 8-16)>]                   Dump  RxNsdA00AXISink Bram"            , mainDumpRxNsdA00AXISinkBram          },
    { "A00"  , "A00Rslv"   ,  1, "<reg>                                   Read  RxNsdA00AXISink Slv Reg"         , mainReadRxNsdA00AXISinkSlvReg        },
    { "A00"  , "A00Wslv"   ,  2, "<reg val>                               Write RxNsdA00AXISink Slv Reg"         , mainWriteRxNsdA00AXISinkSlvReg       },
    { "A00"  , "A00Dslv"   ,  0, ""                                                                              , mainDumpRxNsdA00AXISinkSlvReg        },
    { "A00"  , "A00Dslv"   ,  1, "[<fmt(0=hex 1=bin)>]                    Dump  RxNsdA00AXISink Slv Reg"         , mainDumpRxNsdA00AXISinkSlvReg        },
    { "A00"  , "A00Dmp"    ,  0, "                                        Dump  RxNsdA00AXISink State"           , mainDumpRxNsdA00AXISink              },
    { "A00"  , "A00Rst"    ,  0, "                                        Reset RxNsdA00AXISink"                 , mainResetRxNsdA00AXISink             },
    { "A00"  , ""          ,  0, "",                                                                                                               NULL },
    { "A01"  , "A01RBram"  ,  1, "<adr>                                   Read  RxNsdA01AXISink Bram"            , mainReadRxNsdA01AXISinkBram          },
    { "A01"  , "A01WBram"  ,  2, "<adr val>                               Write RxNsdA01AXISink Bram"            , mainWriteRxNsdA01AXISinkBram         },
    { "A01"  , "A01IBram"  ,  3, "<adr nwd val>                           Init  RxNsdA01AXISink Bram"            , mainInitRxNsdA01AXISinkBram          },
    { "A01"  , "A01DBram"  ,  0, ""                                                                              , mainDumpRxNsdA01AXISinkBram          },
    { "A01"  , "A01DBram"  ,  1, ""                                                                              , mainDumpRxNsdA01AXISinkBram          },
    { "A01"  , "A01DBram"  ,  2, "[<adr nbit(0, 8-16)>]                   Dump  RxNsdA01AXISink Bram"            , mainDumpRxNsdA01AXISinkBram          },
    { "A01"  , "A01Rslv"   ,  1, "<reg>                                   Read  RxNsdA01AXISink Slv Reg"         , mainReadRxNsdA01AXISinkSlvReg        },
    { "A01"  , "A01Wslv"   ,  2, "<reg val>                               Write RxNsdA01AXISink Slv Reg"         , mainWriteRxNsdA01AXISinkSlvReg       },
    { "A01"  , "A01Dslv"   ,  0, ""                                                                              , mainDumpRxNsdA01AXISinkSlvReg        },
    { "A01"  , "A01Dslv"   ,  1, "[<fmt(0=hex 1=bin)>]                    Dump  RxNsdA01AXISink Slv Reg"         , mainDumpRxNsdA01AXISinkSlvReg        },
    { "A01"  , "A01Dmp"    ,  0, "                                        Dump  RxNsdA01AXISink State"           , mainDumpRxNsdA01AXISink              },
    { "A01"  , "A01Rst"    ,  0, "                                        Reset RxNsdA01AXISink"                 , mainResetRxNsdA01AXISink             },
    { "A01"  , ""          ,  0, "",                                                                                                               NULL },
    { "A15"  , "A15RBram"  ,  1, "<adr>                                   Read  RxNsdA15AXISink Bram"            , mainReadRxNsdA15AXISinkBram          },
    { "A15"  , "A15WBram"  ,  2, "<adr val>                               Write RxNsdA15AXISink Bram"            , mainWriteRxNsdA15AXISinkBram         },
    { "A15"  , "A15IBram"  ,  3, "<adr nwd val>                           Init  RxNsdA15AXISink Bram"            , mainInitRxNsdA15AXISinkBram          },
    { "A15"  , "A15DBram"  ,  0, ""                                                                              , mainDumpRxNsdA15AXISinkBram          },
    { "A15"  , "A15DBram"  ,  1, ""                                                                              , mainDumpRxNsdA15AXISinkBram          },
    { "A15"  , "A15DBram"  ,  2, "[<adr nbit(0, 8-16)>]                   Dump  RxNsdA15AXISink Bram"            , mainDumpRxNsdA15AXISinkBram          },
    { "A15"  , "A15Rslv"   ,  1, "<reg>                                   Read  RxNsdA15AXISink Slv Reg"         , mainReadRxNsdA15AXISinkSlvReg        },
    { "A15"  , "A15Wslv"   ,  2, "<reg val>                               Write RxNsdA15AXISink Slv Reg"         , mainWriteRxNsdA15AXISinkSlvReg       },
    { "A15"  , "A15Dslv"   ,  0, ""                                                                              , mainDumpRxNsdA15AXISinkSlvReg        },
    { "A15"  , "A15Dslv"   ,  1, "[<fmt(0=hex 1=bin)>]                    Dump  RxNsdA15AXISink Slv Reg"         , mainDumpRxNsdA15AXISinkSlvReg        },
    { "A15"  , "A15Dmp"    ,  0, "                                        Dump  RxNsdA15AXISink State"           , mainDumpRxNsdA15AXISink              },
    { "A15"  , "A15Rst"    ,  0, "                                        Reset RxNsdA15AXISink"                 , mainResetRxNsdA15AXISink             },
    { "A15"  , ""          ,  0, "",                                                                                                               NULL },
    { "CSK"  , "CskRBram"  ,  1, "<adr>                                   Read  RxNsdCordicAXISink Bram"         , mainReadRxNsdCordicAXISinkBram       },
    { "CSK"  , "CskWBram"  ,  2, "<adr val>                               Write RxNsdCordicAXISink Bram"         , mainWriteRxNsdCordicAXISinkBram      },
    { "CSK"  , "CskIBram"  ,  3, "<adr nwd val>                           Init  RxNsdCordicAXISink Bram"         , mainInitRxNsdCordicAXISinkBram       },
    { "CSK"  , "CskDBram"  ,  0, ""                                                                              , mainDumpRxNsdCordicAXISinkBram       },
    { "CSK"  , "CskDBram"  ,  1, ""                                                                              , mainDumpRxNsdCordicAXISinkBram       },
    { "CSK"  , "CskDBram"  ,  2, "[<adr nbit(0, 8-16)>]                   Dump  RxNsdCordicAXISink Bram"         , mainDumpRxNsdCordicAXISinkBram       },
    { "CSK"  , "CskRslv"   ,  1, "<reg>                                   Read  RxNsdCordicAXISink Slv Reg"      , mainReadRxNsdCordicAXISinkSlvReg     },
    { "CSK"  , "CskWslv"   ,  2, "<reg val>                               Write RxNsdCordicAXISink Slv Reg"      , mainWriteRxNsdCordicAXISinkSlvReg    },
    { "CSK"  , "CskDslv"   ,  0, ""                                                                              , mainDumpRxNsdCordicAXISinkSlvReg     },
    { "CSK"  , "CskDslv"   ,  1, "[<fmt(0=hex 1=bin)>]                    Dump  RxNsdCordicAXISink Slv Reg"      , mainDumpRxNsdCordicAXISinkSlvReg     },
    { "CSK"  , "CskDmp"    ,  0, "                                        Dump  RxNsdCordicAXISink State"        , mainDumpRxNsdCordicAXISink           },
    { "CSK"  , "CskRst"    ,  0, "                                        Reset RxNsdCordicAXISink"              , mainResetRxNsdCordicAXISink          },
    { "CSK"  , ""          ,  0, "",                                                                                                               NULL },
    { "FSK"  , "FskRBram"  ,  1, "<adr>                                   Read  RxNsdFFTAXISink Bram"            , mainReadRxNsdFFTAXISinkBram          },
    { "FSK"  , "FskWBram"  ,  2, "<adr val>                               Write RxNsdFFTAXISink Bram"            , mainWriteRxNsdFFTAXISinkBram         },
    { "FSK"  , "FskIBram"  ,  3, "<adr nwd val>                           Init  RxNsdFFTAXISink Bram"            , mainInitRxNsdFFTAXISinkBram          },
    { "FSK"  , "FskDBram"  ,  0, ""                                                                              , mainDumpRxNsdFFTAXISinkBram          },
    { "FSK"  , "FskDBram"  ,  1, ""                                                                              , mainDumpRxNsdFFTAXISinkBram          },
    { "FSK"  , "FskDBram"  ,  2, "[<adr nbit(0, 8-16)>]                   Dump  RxNsdFFTAXISink Bram"            , mainDumpRxNsdFFTAXISinkBram          },
    { "FSK"  , "FskRslv"   ,  1, "<reg>                                   Read  RxNsdFFTAXISink Slv Reg"         , mainReadRxNsdFFTAXISinkSlvReg        },
    { "FSK"  , "FskWslv"   ,  2, "<reg val>                               Write RxNsdFFTAXISink Slv Reg"         , mainWriteRxNsdFFTAXISinkSlvReg       },
    { "FSK"  , "FskDslv"   ,  0, ""                                                                              , mainDumpRxNsdFFTAXISinkSlvReg        },
    { "FSK"  , "FskDslv"   ,  1, "[<fmt(0=hex 1=bin)>]                    Dump  RxNsdFFTAXISink Slv Reg"         , mainDumpRxNsdFFTAXISinkSlvReg        },
    { "FSK"  , "FskDmp"    ,  0, "                                        Dump  RxNsdFFTAXISink State"           , mainDumpRxNsdFFTAXISink              },
    { "FSK"  , "FskRst"    ,  0, "                                        Reset RxNsdFFTAXISink"                 , mainResetRxNsdFFTAXISink             },
    { "FSK"  , ""          ,  0, "",                                                                                                               NULL },
    { "B00"  , "B00RBram"  ,  1, "<adr>                                   Read  RxNsdB00AXISink Bram"            , mainReadRxNsdB00AXISinkBram          },
    { "B00"  , "B00WBram"  ,  2, "<adr val>                               Write RxNsdB00AXISink Bram"            , mainWriteRxNsdB00AXISinkBram         },
    { "B00"  , "B00IBram"  ,  3, "<adr nwd val>                           Init  RxNsdB00AXISink Bram"            , mainInitRxNsdB00AXISinkBram          },
    { "B00"  , "B00DBram"  ,  0, ""                                                                              , mainDumpRxNsdB00AXISinkBram          },
    { "B00"  , "B00DBram"  ,  1, ""                                                                              , mainDumpRxNsdB00AXISinkBram          },
    { "B00"  , "B00DBram"  ,  2, "[<adr nbit(0, 8-16)>]                   Dump  RxNsdB00AXISink Bram"            , mainDumpRxNsdB00AXISinkBram          },
    { "B00"  , "B00Rslv"   ,  1, "<reg>                                   Read  RxNsdB00AXISink Slv Reg"         , mainReadRxNsdB00AXISinkSlvReg        },
    { "B00"  , "B00Wslv"   ,  2, "<reg val>                               Write RxNsdB00AXISink Slv Reg"         , mainWriteRxNsdB00AXISinkSlvReg       },
    { "B00"  , "B00Dslv"   ,  0, ""                                                                              , mainDumpRxNsdB00AXISinkSlvReg        },
    { "B00"  , "B00Dslv"   ,  1, "[<fmt(0=hex 1=bin)>]                    Dump  RxNsdB00AXISink Slv Reg"         , mainDumpRxNsdB00AXISinkSlvReg        },
    { "B00"  , "B00Dmp"    ,  0, "                                        Dump  RxNsdB00AXISink State"           , mainDumpRxNsdB00AXISink              },
    { "B00"  , "B00Rst"    ,  0, "                                        Reset RxNsdB00AXISink"                 , mainResetRxNsdB00AXISink             },
    { "B00"  , ""          ,  0, "",                                                                                                               NULL },
    { "B01"  , "B01RBram"  ,  1, "<adr>                                   Read  RxNsdB01AXISink Bram"            , mainReadRxNsdB01AXISinkBram          },
    { "B01"  , "B01WBram"  ,  2, "<adr val>                               Write RxNsdB01AXISink Bram"            , mainWriteRxNsdB01AXISinkBram         },
    { "B01"  , "B01IBram"  ,  3, "<adr nwd val>                           Init  RxNsdB01AXISink Bram"            , mainInitRxNsdB01AXISinkBram          },
    { "B01"  , "B01DBram"  ,  0, ""                                                                              , mainDumpRxNsdB01AXISinkBram          },
    { "B01"  , "B01DBram"  ,  1, ""                                                                              , mainDumpRxNsdB01AXISinkBram          },
    { "B01"  , "B01DBram"  ,  2, "[<adr nbit(0, 8-16)>]                   Dump  RxNsdB01AXISink Bram"            , mainDumpRxNsdB01AXISinkBram          },
    { "B01"  , "B01Rslv"   ,  1, "<reg>                                   Read  RxNsdB01AXISink Slv Reg"         , mainReadRxNsdB01AXISinkSlvReg        },
    { "B01"  , "B01Wslv"   ,  2, "<reg val>                               Write RxNsdB01AXISink Slv Reg"         , mainWriteRxNsdB01AXISinkSlvReg       },
    { "B01"  , "B01Dslv"   ,  0, ""                                                                              , mainDumpRxNsdB01AXISinkSlvReg        },
    { "B01"  , "B01Dslv"   ,  1, "[<fmt(0=hex 1=bin)>]                    Dump  RxNsdB01AXISink Slv Reg"         , mainDumpRxNsdB01AXISinkSlvReg        },
    { "B01"  , "B01Dmp"    ,  0, "                                        Dump  RxNsdB01AXISink State"           , mainDumpRxNsdB01AXISink              },
    { "B01"  , "B01Rst"    ,  0, "                                        Reset RxNsdB01AXISink"                 , mainResetRxNsdB01AXISink             },
    { "B01"  , ""          ,  0, "",                                                                                                               NULL },
    { "B15"  , "B15RBram"  ,  1, "<adr>                                   Read  RxNsdB15AXISink Bram"            , mainReadRxNsdB15AXISinkBram          },
    { "B15"  , "B15WBram"  ,  2, "<adr val>                               Write RxNsdB15AXISink Bram"            , mainWriteRxNsdB15AXISinkBram         },
    { "B15"  , "B15IBram"  ,  3, "<adr nwd val>                           Init  RxNsdB15AXISink Bram"            , mainInitRxNsdB15AXISinkBram          },
    { "B15"  , "B15DBram"  ,  0, ""                                                                              , mainDumpRxNsdB15AXISinkBram          },
    { "B15"  , "B15DBram"  ,  1, ""                                                                              , mainDumpRxNsdB15AXISinkBram          },
    { "B15"  , "B15DBram"  ,  2, "[<adr nbit(0, 8-16)>]                   Dump  RxNsdB15AXISink Bram"            , mainDumpRxNsdB15AXISinkBram          },
    { "B15"  , "B15Rslv"   ,  1, "<reg>                                   Read  RxNsdB15AXISink Slv Reg"         , mainReadRxNsdB15AXISinkSlvReg        },
    { "B15"  , "B15Wslv"   ,  2, "<reg val>                               Write RxNsdB15AXISink Slv Reg"         , mainWriteRxNsdB15AXISinkSlvReg       },
    { "B15"  , "B15Dslv"   ,  0, ""                                                                              , mainDumpRxNsdB15AXISinkSlvReg        },
    { "B15"  , "B15Dslv"   ,  1, "[<fmt(0=hex 1=bin)>]                    Dump  RxNsdB15AXISink Slv Reg"         , mainDumpRxNsdB15AXISinkSlvReg        },
    { "B15"  , "B15Dmp"    ,  0, "                                        Dump  RxNsdB15AXISink State"           , mainDumpRxNsdB15AXISink              },
    { "B15"  , "B15Rst"    ,  0, "                                        Reset RxNsdB15AXISink"                 , mainResetRxNsdB15AXISink             },
    { "B15"  , ""          ,  0, "",                                                                                                               NULL },
    { "CIR"  , "CIRRBramFFT",  1, "<adr>                                   Read  RxPrsCIR BramFFT"               , mainReadRxPrsCIRBramFFT              },
    { "CIR"  , "CIRWBramFFT",  2, "<adr val>                               Write RxPrsCIR BramFFT"               , mainWriteRxPrsCIRBramFFT             },
    { "CIR"  , "CIRIBramFFT",  3, "<adr nwd val>                           Init  RxPrsCIR BramFFT"               , mainInitRxPrsCIRBramFFT              },
    { "CIR"  , "CIRDBramFFT",  0, ""                                                                             , mainDumpRxPrsCIRBramFFT              },
    { "CIR"  , "CIRDBramFFT",  1, ""                                                                             , mainDumpRxPrsCIRBramFFT              },
    { "CIR"  , "CIRDBramFFT",  2, "[<adr nbit(0, 8-16)>]                   Dump  RxPrsCIR BramFFT"               , mainDumpRxPrsCIRBramFFT              },
    { "CIR"  , "CIRRBramABS",  1, "<adr>                                   Read  RxPrsCIR BramABS"               , mainReadRxPrsCIRBramABS              },
    { "CIR"  , "CIRWBramABS",  2, "<adr val>                               Write RxPrsCIR BramABS"               , mainWriteRxPrsCIRBramABS             },
    { "CIR"  , "CIRIBramABS",  3, "<adr nwd val>                           Init  RxPrsCIR BramABS"               , mainInitRxPrsCIRBramABS              },
    { "CIR"  , "CIRDBramABS",  0, ""                                                                             , mainDumpRxPrsCIRBramABS              },
    { "CIR"  , "CIRDBramABS",  1, ""                                                                             , mainDumpRxPrsCIRBramABS              },
    { "CIR"  , "CIRDBramABS",  2, "[<adr nbit(0, 8-16)>]                   Dump  RxPrsCIR BramABS"               , mainDumpRxPrsCIRBramABS              },
    { "CIR"  , "CIRRslv"    ,  1, "<reg>                                   Read  RxPrsCIR Slv Reg"               , mainReadRxPrsCIRSlvReg               },
    { "CIR"  , "CIRWslv"    ,  2, "<reg val>                               Write RxPrsCIR Slv Reg"               , mainWriteRxPrsCIRSlvReg              },
    { "CIR"  , "CIRDslv"    ,  0, ""                                                                             , mainDumpRxPrsCIRSlvReg               },
    { "CIR"  , "CIRDslv"    ,  1, "[<fmt(0=hex 1=bin)>]                    Dump  RxPrsCIR Slv Reg"               , mainDumpRxPrsCIRSlvReg               },
    { "CIR"  , "CIRDmp"     ,  0, "                                        Dump  RxPrsCIR State"                 , mainDumpRxPrsCIR                     },
    { "CIR"  , "CIRRst"     ,  0, "                                        Reset RxPrsCIR"                       , mainResetRxPrsCIR                    },
    { "CIR"  , "CIRDelta"   ,  0, "                                        Calc CIR Deltas"                      , mainRxPrsCIRCalc                     },
    { "CIR"  , "CIRRecalc"  ,  0, "                                        Calc CIR Deltas (further pass)"       , mainRxPrsCIRRecalc                   },
    { "CIR"  , "CIRTable"   ,  0, "                                        Dump CIR Table"                       , mainRxPrsCIRDumpCIRTable             },
    { "CIR"  , "CIRPrsRt"   ,  0, "                                        Dump CIR Rx PRS Time Domain"          , mainRxPrsCIRDumpRxPRSTd              },
    { "CIR"  , "CIRPrsRf"   ,  0, "                                        Dump CIR Rx PRS Freq Domain"          , mainRxPrsCIRDumpRxPRSFd              },
    { "CIR"  , "CIRPrsTf"   ,  0, "                                        Dump CIR Tx PRS Freq Domain"          , mainRxPrsCIRDumpTxPRSFd              },
    { "CIR"  , "CIRPrsRtCT" ,  0, "                                        Dump CIR Rx PRS TD Corr Delta T"      , mainRxPrsCIRDumpRxPRSTdCorDeltaT     },
    { "CIR"  , "CIRPrsRtCf" ,  0, "                                        Dump CIR Rx PRS TD Corr delta f'"     , mainRxPrsCIRDumpRxPRSTdCordeltafp    },
    { "CIR"  , "CIRPrsRfCf" ,  0, "                                        Dump CIR Rx PRS FD Dorr delta f'"     , mainRxPrsCIRDumpRxPRSFdCordeltafp    },
    { "CIR"  , "CIRPrsRfCt" ,  0, "                                        Dump CIR Rx PRS FD Corr delta t'"     , mainRxPrsCIRDumpRxPRSFdCordeltatp    },
    { "CIR"  , "CIRPrsRtCFT",  0, "                                        Dump CIR Rx PRS TD Corr Freq+Time"    , mainRxPrsCIRDumpRxPRSTdCorFT         },
    { "CIR"  , ""           ,  0, "",                                                                                                              NULL },
    { "BPA"  , "BPARslv"    ,  1, "<reg>                                   Read  RxPrsBPA Slv Reg"               , mainReadRxPrsBPASlvReg               },
    { "BPA"  , "BPAWslv"    ,  2, "<reg val>                               Write RxPrsBPA Slv Reg"               , mainWriteRxPrsBPASlvReg              },
    { "BPA"  , "BPADslv"    ,  0, ""                                                                             , mainDumpRxPrsBPASlvReg               },
    { "BPA"  , "BPADslv"    ,  1, "[<fmt(0=hex 1=bin)>]                    Dump  RxPrsBPA Slv Reg"               , mainDumpRxPrsBPASlvReg               },
    { "BPA"  , "BPADmp"     ,  0, "                                        Dump  RxPrsBPA State"                 , mainDumpRxPrsBPA                     },
    { "BPA"  , "BPARst"     ,  0, "                                        Reset RxPrsBPA"                       , mainResetRxPrsBPA                    },
    { "BPA"  , "BPADelta"   ,  0, "                                        RxPrsBPA Calc Deltas"                 , mainRxPrsBPACalcDelta                },
    { "BPA"  , "BPADelta2"  ,  0, "                                        RxPrsBPA Calc Deltas, 2 pass CIR"     , mainRxPrsBPACalcDelta2               },
    { "BPA"  , "BPAdCSFd"   ,  0, "                                        RxPrsBPA Dump Corr Symb Freq Dom"     , mainRxPrsBPADumpCorSymbFd            },
    { "BPA"  , "BPAdCSTd"   ,  0, "                                        RxPrsBPA Dump Corr Symb Time Dom"     , mainRxPrsBPADumpCorSymbTd            },
    { "BPA"  , "BPAdatFd"   ,  0, "                                        RxData BPA Dump Corr Symb Freq Dom"   , mainRxDataBPADumpCorSymbFd           },
    { "BPA"  , "BPAdatTd"   ,  0, "                                        RxDataBPA Dump Corr Symb Time Dom"    , mainRxDataBPADumpCorSymbTd           },
    { "BPA"  , "PrsBPA"     ,  0, "                                        Enter BPA Shell"                      , mainPrsBPAShell                      },
    { "BPA"  , ""           ,  0, "",                                                                                                              NULL },
    { "BPA"  , "DSRSFd"     ,  0, "                                        dump Data Symbol raw FD"              , mainRxDataSymbolDumpRawFd                      },
    { "BPA"  , "DSCFFd"     ,  0, "                                        dump Data Symbol freq corr FD"        , mainRxDataSymbolDumpFCFd                      },
    { "BPA"  , "DSCFTSFd"   ,  0, "                                        dump Data Symbol freq and time FD"    , mainRxDataSymbolDumpFTCFd                      },
    { "BPA"  , ""           ,  0, "",                                                                                                              NULL },
//===========================================
//////// End of script-generated code ///////
//===========================================
    { "EMS"  , "EmsSelChn" ,  1, "<0 = RF 1=Emul>                         Select Emul Switch Chan"               , mainSelChRFChanEmulSw                },
    { "EMS"  , ""          ,  0, "",                                                                                                               NULL },
    { "RF"   , "ad9364"    ,  0, "                                        AD9364 Shell"                          , mainAD9364Shell                      },
    { "RF"   , "ad9364i"   ,  1, "0=26.00, 1=30.72, 2=61.44MHz            AD9364 Init"                           , mainInitAD9364                       },
    { "PSK"  , ""          ,  0, "",                                                                                                               NULL },
    { "NSD"  , "NsdTest"   ,  1, "<0 off 1 ctr 2 1shot>                   RxNSDet Test (set mode, start)"        , mainTestRxNSDet                      },
    { "NSD"  , "NsdTMode"  ,  1, "<0 off 1 ctr 2 1shot>                   Set RxNSDet Test Mode"                 , mainSetTmodeRxNSDet                  },
    { "NSD"  , "NsdStart"  ,  0, "                                        Start RxNSDet"                         , mainStartRxNSDet                     },
    { "NSD"  , "NsdStop"   ,  0, "                                        Stop RxNSDet"                          , mainStopTmodeRxNSDet                 },
    { "NSD"  , "NsdDstate" ,  0, "                                        Set RxNSDet Test Mode"                 , mainDstateRxNSDet                    },
    { "NSD"  , ""          ,  0, "",                                                                                                               NULL },
    { "all"  , "SysParam"  ,  1, "?                                       Show System Params"                    , mainSetSysParams                     },
    { "all"  , "SysParam"  ,  1, "default                                 Set  System Params to Default Values"  , mainSetSysParams                     },
    { "all"  , "SysParam"  , 10, "<Sys Param List>                        Set  System Params"                    , mainSetSysParams                     },
    { "all"  , "RunTest"   ,  3, "<selEmulChan rmean srcmode>             Run Test"                              , mainRunTest                          },
    { "all"  , "TBPA"      ,  0, "                                        Run BPA Test"                          , mainRunTestBPA                       },
    { "all"  , "StopTest"  ,  0, "                                        Stop Test"                             , mainStopTest                         },
    { "all"  , "ChkTsiFrq" ,  0, "                                        Check TSI Tx Frequency"                , mainChkTSITxFreq                     },
    { "all"  , "ChkNsdFrq" ,  0, "                                        Check NSD Rx Frequency"                , mainChkNSDRxFreq                     },
    { "all"  , "InitSys"   ,  0, "                                        Initialise and reset HW"               , mainSysInit                          },
    { "all"  , "ResetSys"  ,  0, "                                        Global HW Reset"                       , mainSysReset                         },
    { "\0" }
};



int main()
{
	Xil_ICacheEnable();
	Xil_DCacheEnable();
	SysTimer   = SysTimerAPI;
	if( ( SysTimerId = SysTimer->InitTimer( SysTimer, SYS_TIMER_BASEADDR, ID_NEW_DEVICE, "System Timer", SYS_TIMER_DEVICE_ID ) ) == ID_DEVICE_ERROR )
	{
		cprintf( "\nError Initialising System Timer: %s./n", SysTimer->DevManErrMsg( SysTimer ) );
	}
	else
	{
		SysTimer->StartTimer( SysTimer, SysTimerId );
		dispSysTimerRunning();
	}
	( void ) clkWizInit( &clkWizDac, RFCHAN_TREMUL_CLKWIZ_BASEADDR, 56, 0x00, " Emulator Clock Generator" );


	FLUSH;

    mainShowBanner( NULL );
	mainSysInit( NULL );
//	( void ) pause("\nPaused - hit any key");

	( void ) ad9364_Init_EMM( DATA_SEL_DMA, FS_SEL_61_44_MHz, 1 );
//	( void ) pause("\nPaused - hit any key");

	//( void ) shellConsole( &mainShell, mainCmdDef, "main shell >", "v 0; loop 10; echo test\\_%l\\n; runtest 0 6 1; wait 0.1; endl; v 1; cirdelta; bpadelta;", 1000 );
	( void ) shellConsole( &mainShell, mainCmdDef, "main shell >", "v 0; loop 3; echo test\\_%l\\n; runtest 0 6 0; wait 0.1; cirdelta;bpadelta;endl; ", 1000 );
	cprintf( "---Exiting main---\n\r" );
	Xil_DCacheDisable();
	Xil_ICacheDisable();
	return 0;
}
/* End of main() */

/************************************************************
 System Testing
 ************************************************************/
void mainStopTest( shell * psh )
{
    StopPLSystem( pEMMSysCtrl );
}

void mainRunTest( shell *psh ) {
    mainStopTest( psh );
    mainSysReset( psh );
    u32 selEmulChan = psh->uarg[ 0 ];
    double rmean    = psh->darg[ 1 ];
    u32 srcmode     = psh->uarg[ 2 ];
	SysRunTest( pEMMSysCtrl, selEmulChan, rmean, srcmode );
//    mainRxPrsCIRCalc( psh );
}

void mainRunTestBPA( shell *psh ) {
    u32 selEmulChan = 2;
    double rmean    = 6.0;
    u32 srcmode     = AXISOURCE_ONESHOT;
    SysRunTest( pEMMSysCtrl, selEmulChan, rmean, srcmode );
    WaitSec( 1.00 );
    mainStopTest( psh );
    mainRxPrsBPACalcDelta( psh );
}


void mainSetSysParams ( shell *psh ) {
    if ( !strcmp ( psh->sarg [ 0 ], "?" ) ) {
        SysSetSysParams ( true, false, pEMMSysCtrl, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 );
    } else if ( !strcmp ( psh->sarg [ 0 ], "default" ) ) {
        SysSetSysParams ( false, true, pEMMSysCtrl, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 );
    } else {
        u32 nCPre         = psh->uarg [ 0 ];
        u32 nFFT          = psh->uarg [ 1 ];
        u32 nSymbData     = psh->uarg [ 2 ];
        u32 nFrameBlock   = psh->uarg [ 3 ];
        u32 nSymbPRS      = psh->uarg [ 4 ];
        u32 nSymbAGC      = psh->uarg [ 5 ];
        u32 nSymbNSD      = psh->uarg [ 6 ];
        u32 TxRxFreq      = psh->uarg [ 7 ];
        u32 BWPercent     = psh->uarg [ 8 ];
        ModType ModType   = psh->uarg [ 9 ];
        SoftBits SoftBits = FOUR_SOFT_BITS;
        SysSetSysParams ( false, false, pEMMSysCtrl,
                nCPre,
                nFFT,
                nSymbData,
                nFrameBlock,
                nSymbPRS,
                nSymbAGC,
                nSymbNSD,
                TxRxFreq,
                BWPercent,
                ModType,
                SoftBits
                );
    }
}


/************************************************************
 User Interface
 ************************************************************/
void mainShowBanner( shell *psh )
{
	cprintf( "\n" );
	cprintf( "*************************************************************\n" );
	cprintf( "User Console for Telemetry System Debug\n" );
	cprintf( "G. Mortimer 2017\n" );
	cprintf( "( c ) Relevant Technologies Ltd\n" );
	cprintf( "( c ) EM Motorsport Ltd\n" );
	cprintf( "*************************************************************\n" );
	cprintf( "\n" );
}

void mainShowHelp( shell *psh )
{
	shellShowHelp( psh );
}


void mainToggleEscSeqDet( shell *psh )
{
	if ( escSeqDet )
	{
		escSeqDet = false;
		cprintf( "ESC sequences inhibited\n" );
	}
	else
	{
		escSeqDet = true;
		cprintf( "ESC sequences enabled\n" );
	}
}

void mainTestKeyboard( shell *psh )
{
	testKB();
}

void mainMemDump( shell *psh )
{
	static u32 addr  = 0x100000;
	static u32 fmt   = DUMP_FMT_HEX;
	static u32 arg1  = 4;
	static u32 arg2  = 4096;
	mainArgs( &addr, &fmt, &arg1, &arg2, psh );

	utilScreenDumpData( DUMP_SEL_PS7_RAM, fmt, DUMP_TYP_MEMORY, addr, arg1, arg2, 0, "PS7 RAM:" );
}


void mainArgs( u32 *arg0, u32 *arg1, u32 *arg2, u32 *arg3, shell *psh )
{
	switch ( psh->narg )
	{
	case 4:
		*arg3 = psh->uarg[ 3 ];
	case 3:
		*arg2 = psh->uarg[ 2 ];
	case 2:
		*arg1 = psh->uarg[ 1 ];
	case 1:
		*arg0 = psh->uarg[ 0 ];
	case 0:
	default:
		break;
	}
}


/************************************************************
 Access to PrsBPA
 ************************************************************/
void mainPrsBPAShell( shell *psh )
{
    shellConsole( &PrsBPAShell, PrsBPACmdDef, "PrsBPA >", "h", 1000 );
}

/************************************************************
 Access to AD9364
 ************************************************************/
void mainAD9364Shell( shell *psh )
{
	shellConsole( &AD9364Shell, AD9364CmdDef, "AD9364 >", "h", 1000 );
}

void mainInitAD9364( shell *psh )
{
	enum dds_data_select  dac_data_sel    = DATA_SEL_DMA;
	enum samp_freq_select sample_freq_sel = *( enum samp_freq_select* )&psh->uarg[ 0 ];
	uint32_t              FIR_enable      = 1;
	( void ) ad9364_Init_EMM( dac_data_sel, sample_freq_sel, FIR_enable );
}


/************************************************************
 Debug Testing
 ************************************************************/
void mainTestArg( shell *psh )
{
	cprintf( "hex    : %08x\nuint   : %u\nint    : %i\nstring : %s\ndouble : %g\n", toUint( psh->uarg[ 0 ] ), toUint( psh->uarg[ 0 ] ), toInt( psh->iarg[ 0 ] ), psh->sarg[ 0 ], psh->darg[ 0 ] );
}

/************************************************************
 Test Buffer Initialisation
 ************************************************************/

/************************************************************
 HW System Reset - Dedicated Module Pending
 ************************************************************/
void mainSysReset( shell *psh )
{
	DVC_RESET( pEMMSysCtrl->pTxSourceAXIBramC );
	DVC_RESET( pEMMSysCtrl->pTxAXISource );
	DVC_RESET( pEMMSysCtrl->pTxFFTC );
	DVC_RESET( pEMMSysCtrl->pRxNsdFFTC );
	DVC_RESET( pEMMSysCtrl->pTxTSInj );
	DVC_RESET( pEMMSysCtrl->pRFChanTREmul );
	DVC_RESET( pEMMSysCtrl->pRFChanTREmulSw );
	DVC_RESET( pEMMSysCtrl->pRFChanTREmulFreqCtr );
	DVC_RESET( pEMMSysCtrl->pAD9364ClkFreqCtr );
	DVC_RESET( pEMMSysCtrl->pAD9364DACFreqCtr );
	DVC_RESET( pEMMSysCtrl->pAD9364ADCFreqCtr );
	DVC_RESET( pEMMSysCtrl->pRxNSDet );
	DVC_RESET( pEMMSysCtrl->pRxTddAXISink );
    DVC_RESET( pEMMSysCtrl->pRxPrsAXISink );
    DVC_RESET( pEMMSysCtrl->pRxDtaAXISink );
	DVC_RESET( pEMMSysCtrl->pRxNsdA00AXISink );
    DVC_RESET( pEMMSysCtrl->pRxNsdA01AXISink );
    DVC_RESET( pEMMSysCtrl->pRxNsdA15AXISink );
    DVC_RESET( pEMMSysCtrl->pRxNsdFFTAXISink );
    DVC_RESET( pEMMSysCtrl->pRxNsdCordicAXISink );
    DVC_RESET( pEMMSysCtrl->pRxNsdB00AXISink );
    DVC_RESET( pEMMSysCtrl->pRxNsdB01AXISink );
    DVC_RESET( pEMMSysCtrl->pRxNsdB15AXISink );
#ifdef verbose
	cprintf( "All Devices Reset\n" );
#endif
}

/************************************************************
 HW System Init
 ************************************************************/
void mainSysInit( shell * psh )
{
	InitPLSystem( pEMMSysCtrl );
	cprintf( "\n" );
}

/************************************************************
 NSD Testing
 ************************************************************/
void mainChkTSITxFreq( shell * psh )
{
    ChkTSITxFreq();
}

void mainChkNSDRxFreq( shell * psh )
{
    ChkNSDRxFreq();
}

void mainTestRxNSDet( shell *psh )
{
	u32 tstmode = *( u32 *)psh->uarg[0];
	DVC_PROC( pEMMSysCtrl->pRxNSDet, SetTstMode, min( tstmode, NSDET_TM_ONE_SHOT ) );
	DVC_PROC( pEMMSysCtrl->pRxNSDet, Start );
}

void mainSetTmodeRxNSDet( shell *psh )
{
	u32 tstmode = *( u32 *)psh->uarg[0];
	DVC_PROC( pEMMSysCtrl->pRxNSDet, SetTstMode, min( tstmode, NSDET_TM_ONE_SHOT ) );
}

void mainStartRxNSDet( shell *psh )
{
	DVC_PROC( pEMMSysCtrl->pRxNSDet, Start );
}

void mainStopTmodeRxNSDet( shell *psh )
{
	DVC_PROC( pEMMSysCtrl->pRxNSDet, Stop );
}

void mainDstateRxNSDet( shell *psh )
{
	DVC_PROC( pEMMSysCtrl->pRxNSDet, Dstate );
}

void mainSelChRFChanEmulSw( shell *psh )
{
	DVC_PROC( pEMMSysCtrl->pRFChanTREmulSw, SelChan, psh->uarg[ 0 ] );
}

/************************************************************
 Timer
 ************************************************************/
void mainResetTimer( shell *psh )
{
	SysTimer->ResetTimer( SysTimer, SysTimerId );
	cprintf( "Timer reset\n" );
}

void mainStartTimer( shell *psh )
{
	SysTimer->StartTimer( SysTimer, SysTimerId );
	cprintf( "Timer started\n" );
}

void mainStopTimer( shell *psh )
{
	SysTimer->StopTimer( SysTimer, SysTimerId );
	cprintf( "Timer stopped\n" );
}

void mainClearTimer( shell *psh )
{
	SysTimer->ClearTimer( SysTimer, SysTimerId );
	cprintf( "Timer cleared\n" );
}

void mainReadTimer( shell *psh )
{
	double f = SysTimer->ReadTimerLo( SysTimer, SysTimerId );
	cprintf( "Time = %.9fs\n", f/PL_FABRIC_CLOCK_FREQ );
}

void mainReadClearTimer( shell *psh )
{
	double f = SysTimer->ReadClearTimer( SysTimer, SysTimerId );
	cprintf( "Time = %.9fs\n", f/PL_FABRIC_CLOCK_FREQ );
}

void mainReadTimerLong( shell *psh )
{
	cprintf( "Time ( long ) = %.9fs\n", SysTimer->ReadTimerLongSec( SysTimer, SysTimerId ) );
}

void mainDumpTimerSlvReg( shell *psh )
{
	SysTimer->DumpTimerSlvRegs( SysTimer, SysTimerId, psh->parg->uarg[ 0 ] );
}

/************************************************************
 AXI BRAM Controller - test only
 ************************************************************/
//static void InitAxiBramCtrlrWrBuf( void )
//{
//	PRNG32 rng = { 0x600D1DEA, 0, prng32s }; PRNG32 *prng = &rng;
//	( void ) prng->prng( prng, true );
//	for( int i = 0; i < AXI_BRAM_CTRLR_BRAM_SIZE; i++ ) axiBramCtrlrWrBuf[ i ] = prng->prng( prng, false );
//}
//
//void mainTestXilAxiBramCtrlr( shell *psh )
//{
//	u32 *prd = axiBramCtrlrWrBuf;
//	u32 *pwr = axiBramCtrlrRdBuf;
//	u32 loops = psh->uarg[ 0 ];
//    InitAxiBramCtrlrWrBuf();
//
//	TIME (
//	for ( int j = 0; j < loops; j++ )
//	{
//		pwr = axiBramCtrlrRdBuf;
//		for ( int i = 0; i < AXI_BRAM_CTRLR_BRAM_SIZE; i++ ) *( u32 * ) ( AXI_BRAM_CTRLR_BASE_ADDR + 4 * i ) = *pwr++;
//	}, loops * AXI_BRAM_CTRLR_BRAM_SIZE, 4 );
//
//	TIME (
//	for ( int j = 0; j < loops; j++ ) {
//		prd = axiBramCtrlrWrBuf;
//		for ( int i = 0; i < AXI_BRAM_CTRLR_BRAM_SIZE; i++ ) {
//			*prd = *( u32 * ) ( AXI_BRAM_CTRLR_BASE_ADDR + 4 * i );
//			if ( *prd++ != *pwr++ ) {
//				cprintf( "\nError at address %08x\n", ( unsigned int )i );
//				break;
//			}
//		}
//	}, loops * AXI_BRAM_CTRLR_BRAM_SIZE, 4 );
//}



/************************************************************
 Device Standard Functions
 ************************************************************/

//===========================================
//// Script-generated code - do not edit ////
//===========================================
//=============
// TxAXISource
//=============
void mainReadTxAXISourceBram( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 data  = pEMMSysCtrl->pTxAXISource->ReadBram( pEMMSysCtrl->pTxAXISource, addr );
    cprintf( " data = %u 0x%08x\n", toUint( data ), toUint( data ) );
}

void mainWriteTxAXISourceBram( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 data  = psh->uarg[ 1 ];
    pEMMSysCtrl->pTxAXISource->WriteBram( pEMMSysCtrl->pTxAXISource, addr, data );
}

void mainInitTxAXISourceBram( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 nwds  = psh->uarg[ 1 ];
    u32 data  = psh->uarg[ 2 ];
    TIME_TFR( pEMMSysCtrl->pTxAXISource->InitBram( pEMMSysCtrl->pTxAXISource, addr, nwds, data ), nwds, 4 );
}

void mainLoadTxAXISourceBram( shell *psh )
{
    TIME_TFR( pEMMSysCtrl->pTxAXISource->LoadBram( pEMMSysCtrl->pTxAXISource ), DVC_PAR( pEMMSysCtrl->pTxAXISource, BramSize ), 4 );
}

void mainDumpTxAXISourceBram( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 nbits = psh->uarg[ 1 ];
    pEMMSysCtrl->pTxAXISource->DumpBram( pEMMSysCtrl->pTxAXISource, addr, nbits );
}

void mainReadTxAXISourceSlvReg( shell *psh )
{
    u32 reg = psh->uarg[ 0 ];
    u32 val = DVC_RD_SLV_REG( pEMMSysCtrl->pTxAXISource, reg );
    cprintf( " data = %u 0x%08x\n", toUint( val ), toUint( val ) );
}

void mainWriteTxAXISourceSlvReg( shell *psh )
{
    u32 reg = psh->uarg[ 0 ];
    u32 val = psh->uarg[ 1 ];
    DVC_WR_SLV_REG( pEMMSysCtrl->pTxAXISource, reg, val );
    cprintf( "\n" );
}

void mainDumpTxAXISourceSlvReg( shell *psh )
{
    u32 bBin = psh->uarg[ 0 ];
    DVC_DUMP_SLV_REG( pEMMSysCtrl->pTxAXISource, bBin );
}

void mainDumpTxAXISource( shell *psh )
{
    DVC_DUMP( pEMMSysCtrl->pTxAXISource );
}

void mainResetTxAXISource( shell *psh )
{
    DVC_RESET( pEMMSysCtrl->pTxAXISource );
    cprintf( "Device TxAXISource Reset\n" );
}

//===============
// TxSourceBramC
//===============
void mainReadTxSourceAXIBramCBram( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 data  = pEMMSysCtrl->pTxSourceAXIBramC->ReadBram( pEMMSysCtrl->pTxSourceAXIBramC, addr );
    cprintf( " data = %u 0x%08x\n", toUint( data ), toUint( data ) );
}

void mainWriteTxSourceAXIBramCBram( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 data  = psh->uarg[ 1 ];
    pEMMSysCtrl->pTxSourceAXIBramC->WriteBram( pEMMSysCtrl->pTxSourceAXIBramC, addr, data );
}

void mainInitTxSourceAXIBramCBram( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 nwds  = psh->uarg[ 1 ];
    u32 data  = psh->uarg[ 2 ];
    TIME_TFR( pEMMSysCtrl->pTxSourceAXIBramC->InitBram( pEMMSysCtrl->pTxSourceAXIBramC, addr, nwds, data ), nwds, 4 );
}

void mainLoadTxSourceAXIBramCBram( shell *psh )
{
    TIME_TFR( pEMMSysCtrl->pTxSourceAXIBramC->LoadBram( pEMMSysCtrl->pTxSourceAXIBramC ), DVC_PAR( pEMMSysCtrl->pTxSourceAXIBramC, BramSize ), 4 );
}

void mainReadTxSourceAXIBramCSlvReg( shell *psh )
{
    u32 reg = psh->uarg[ 0 ];
    u32 val = DVC_RD_SLV_REG( pEMMSysCtrl->pTxSourceAXIBramC, reg );
    cprintf( " data = %u 0x%08x\n", toUint( val ), toUint( val ) );
}

void mainWriteTxSourceAXIBramCSlvReg( shell *psh )
{
    u32 reg = psh->uarg[ 0 ];
    u32 val = psh->uarg[ 1 ];
    DVC_WR_SLV_REG( pEMMSysCtrl->pTxSourceAXIBramC, reg, val );
    cprintf( "\n" );
}

void mainDumpTxSourceAXIBramCSlvReg( shell *psh )
{
    u32 bBin = psh->uarg[ 0 ];
    DVC_DUMP_SLV_REG( pEMMSysCtrl->pTxSourceAXIBramC, bBin );
}

void mainDumpTxSourceAXIBramC( shell *psh )
{
    DVC_DUMP( pEMMSysCtrl->pTxSourceAXIBramC );
}

void mainResetTxSourceAXIBramC( shell *psh )
{
    DVC_RESET( pEMMSysCtrl->pTxSourceAXIBramC );
    cprintf( "Device TxSourceAXIBramC Reset\n" );
}

//========
// TxFFTC
//========
void mainReadTxFFTCSlvReg( shell *psh )
{
    u32 reg = psh->uarg[ 0 ];
    u32 val = DVC_RD_SLV_REG( pEMMSysCtrl->pTxFFTC, reg );
    cprintf( " data = %u 0x%08x\n", toUint( val ), toUint( val ) );
}

void mainWriteTxFFTCSlvReg( shell *psh )
{
    u32 reg = psh->uarg[ 0 ];
    u32 val = psh->uarg[ 1 ];
    DVC_WR_SLV_REG( pEMMSysCtrl->pTxFFTC, reg, val );
    cprintf( "\n" );
}

void mainDumpTxFFTCSlvReg( shell *psh )
{
    u32 bBin = psh->uarg[ 0 ];
    DVC_DUMP_SLV_REG( pEMMSysCtrl->pTxFFTC, bBin );
}

void mainDumpTxFFTC( shell *psh )
{
    DVC_DUMP( pEMMSysCtrl->pTxFFTC );
}

void mainResetTxFFTC( shell *psh )
{
    DVC_RESET( pEMMSysCtrl->pTxFFTC );
    cprintf( "Device TxFFTC Reset\n" );
}

//===========
// RxNsdFFTC
//===========
void mainReadRxNsdFFTCSlvReg( shell *psh )
{
    u32 reg = psh->uarg[ 0 ];
    u32 val = DVC_RD_SLV_REG( pEMMSysCtrl->pRxNsdFFTC, reg );
    cprintf( " data = %u 0x%08x\n", toUint( val ), toUint( val ) );
}

void mainWriteRxNsdFFTCSlvReg( shell *psh )
{
    u32 reg = psh->uarg[ 0 ];
    u32 val = psh->uarg[ 1 ];
    DVC_WR_SLV_REG( pEMMSysCtrl->pRxNsdFFTC, reg, val );
    cprintf( "\n" );
}

void mainDumpRxNsdFFTCSlvReg( shell *psh )
{
    u32 bBin = psh->uarg[ 0 ];
    DVC_DUMP_SLV_REG( pEMMSysCtrl->pRxNsdFFTC, bBin );
}

void mainDumpRxNsdFFTC( shell *psh )
{
    DVC_DUMP( pEMMSysCtrl->pRxNsdFFTC );
}

void mainResetRxNsdFFTC( shell *psh )
{
    DVC_RESET( pEMMSysCtrl->pRxNsdFFTC );
    cprintf( "Device RxNsdFFTC Reset\n" );
}

//=========
// TxTSInj
//=========
void mainReadTxTSInjTsBuf( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 data  = pEMMSysCtrl->pTxTSInj->ReadTsBuf( pEMMSysCtrl->pTxTSInj, addr );
    cprintf( " data = %u 0x%08x\n", toUint( data ), toUint( data ) );
}

void mainWriteTxTSInjTsBuf( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 data  = psh->uarg[ 1 ];
    pEMMSysCtrl->pTxTSInj->WriteTsBuf( pEMMSysCtrl->pTxTSInj, addr, data );
}

void mainInitTxTSInjTsBuf( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 nwds  = psh->uarg[ 1 ];
    u32 data  = psh->uarg[ 2 ];
    TIME_TFR( pEMMSysCtrl->pTxTSInj->InitTsBuf( pEMMSysCtrl->pTxTSInj, addr, nwds, data ), nwds, 4 );
}

void mainLoadTxTSInjTsBuf( shell *psh )
{
    TIME_TFR( pEMMSysCtrl->pTxTSInj->LoadTsBuf( pEMMSysCtrl->pTxTSInj ), DVC_PAR( pEMMSysCtrl->pTxTSInj, TsBufSize ), 4 );
}

void mainDumpTxTSInjTsBuf( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 nbits = psh->uarg[ 1 ];
    pEMMSysCtrl->pTxTSInj->DumpTsBuf( pEMMSysCtrl->pTxTSInj, addr, nbits );
}

void mainReadTxTSInjTxBuf( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 data  = pEMMSysCtrl->pTxTSInj->ReadTxBuf( pEMMSysCtrl->pTxTSInj, addr );
    cprintf( " data = %u 0x%08x\n", toUint( data ), toUint( data ) );
}

void mainDumpTxTSInjTxBuf( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 nbits = psh->uarg[ 1 ];
    pEMMSysCtrl->pTxTSInj->DumpTxBuf( pEMMSysCtrl->pTxTSInj, addr, nbits );
}

void mainWriteTxTSInjRxBuf( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 data  = psh->uarg[ 1 ];
    pEMMSysCtrl->pTxTSInj->WriteRxBuf( pEMMSysCtrl->pTxTSInj, addr, data );
}

void mainInitTxTSInjRxBuf( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 nwds  = psh->uarg[ 1 ];
    u32 data  = psh->uarg[ 2 ];
    TIME_TFR( pEMMSysCtrl->pTxTSInj->InitRxBuf( pEMMSysCtrl->pTxTSInj, addr, nwds, data ), nwds, 4 );
}

void mainReadTxTSInjPsBuf( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 data  = pEMMSysCtrl->pTxTSInj->ReadPsBuf( pEMMSysCtrl->pTxTSInj, addr );
    cprintf( " data = %u 0x%08x\n", toUint( data ), toUint( data ) );
}

void mainDumpTxTSInjPsBuf( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 nbits = psh->uarg[ 1 ];
    pEMMSysCtrl->pTxTSInj->DumpPsBuf( pEMMSysCtrl->pTxTSInj, addr, nbits );
}

void mainReadTxTSInjSlvReg( shell *psh )
{
    u32 reg = psh->uarg[ 0 ];
    u32 val = DVC_RD_SLV_REG( pEMMSysCtrl->pTxTSInj, reg );
    cprintf( " data = %u 0x%08x\n", toUint( val ), toUint( val ) );
}

void mainWriteTxTSInjSlvReg( shell *psh )
{
    u32 reg = psh->uarg[ 0 ];
    u32 val = psh->uarg[ 1 ];
    DVC_WR_SLV_REG( pEMMSysCtrl->pTxTSInj, reg, val );
    cprintf( "\n" );
}

void mainDumpTxTSInjSlvReg( shell *psh )
{
    u32 bBin = psh->uarg[ 0 ];
    DVC_DUMP_SLV_REG( pEMMSysCtrl->pTxTSInj, bBin );
}

void mainDumpTxTSInj( shell *psh )
{
    DVC_DUMP( pEMMSysCtrl->pTxTSInj );
}

void mainResetTxTSInj( shell *psh )
{
    DVC_RESET( pEMMSysCtrl->pTxTSInj );
    cprintf( "Device TxTSInj Reset\n" );
}

//==============
// RFChanTREmul
//==============
void mainReadRFChanTREmulPsLUT( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 data  = pEMMSysCtrl->pRFChanTREmul->ReadPsLUT( pEMMSysCtrl->pRFChanTREmul, addr );
    cprintf( " data = %u 0x%08x\n", toUint( data ), toUint( data ) );
}

void mainWriteRFChanTREmulPsLUT( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 data  = psh->uarg[ 1 ];
    pEMMSysCtrl->pRFChanTREmul->WritePsLUT( pEMMSysCtrl->pRFChanTREmul, addr, data );
}

void mainInitRFChanTREmulPsLUT( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 nwds  = psh->uarg[ 1 ];
    u32 data  = psh->uarg[ 2 ];
    TIME_TFR( pEMMSysCtrl->pRFChanTREmul->InitPsLUT( pEMMSysCtrl->pRFChanTREmul, addr, nwds, data ), nwds, 4 );
}

void mainLoadRFChanTREmulPsLUT( shell *psh )
{
    TIME_TFR( pEMMSysCtrl->pRFChanTREmul->LoadPsLUT( pEMMSysCtrl->pRFChanTREmul ), DVC_PAR( pEMMSysCtrl->pRFChanTREmul, PsLUTSize ), 4 );
}

void mainDumpRFChanTREmulPsLUT( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 nbits = psh->uarg[ 1 ];
    pEMMSysCtrl->pRFChanTREmul->DumpPsLUT( pEMMSysCtrl->pRFChanTREmul, addr, nbits );
}

void mainReadRFChanTREmulHwLUT( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 data  = pEMMSysCtrl->pRFChanTREmul->ReadHwLUT( pEMMSysCtrl->pRFChanTREmul, addr );
    cprintf( " data = %u 0x%08x\n", toUint( data ), toUint( data ) );
}

void mainWriteRFChanTREmulHwLUT( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 data  = psh->uarg[ 1 ];
    pEMMSysCtrl->pRFChanTREmul->WriteHwLUT( pEMMSysCtrl->pRFChanTREmul, addr, data );
}

void mainInitRFChanTREmulHwLUT( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 nwds  = psh->uarg[ 1 ];
    u32 data  = psh->uarg[ 2 ];
    TIME_TFR( pEMMSysCtrl->pRFChanTREmul->InitHwLUT( pEMMSysCtrl->pRFChanTREmul, addr, nwds, data ), nwds, 4 );
}

void mainLoadRFChanTREmulHwLUT( shell *psh )
{
    TIME_TFR( pEMMSysCtrl->pRFChanTREmul->LoadHwLUT( pEMMSysCtrl->pRFChanTREmul ), DVC_PAR( pEMMSysCtrl->pRFChanTREmul, HwLUTSize ), 4 );
}

void mainDumpRFChanTREmulHwLUT( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 nbits = psh->uarg[ 1 ];
    pEMMSysCtrl->pRFChanTREmul->DumpHwLUT( pEMMSysCtrl->pRFChanTREmul, addr, nbits );
}

void mainReadRFChanTREmulSlvReg( shell *psh )
{
    u32 reg = psh->uarg[ 0 ];
    u32 val = DVC_RD_SLV_REG( pEMMSysCtrl->pRFChanTREmul, reg );
    cprintf( " data = %u 0x%08x\n", toUint( val ), toUint( val ) );
}

void mainWriteRFChanTREmulSlvReg( shell *psh )
{
    u32 reg = psh->uarg[ 0 ];
    u32 val = psh->uarg[ 1 ];
    DVC_WR_SLV_REG( pEMMSysCtrl->pRFChanTREmul, reg, val );
    cprintf( "\n" );
}

void mainDumpRFChanTREmulSlvReg( shell *psh )
{
    u32 bBin = psh->uarg[ 0 ];
    DVC_DUMP_SLV_REG( pEMMSysCtrl->pRFChanTREmul, bBin );
}

void mainDumpRFChanTREmul( shell *psh )
{
    DVC_DUMP( pEMMSysCtrl->pRFChanTREmul );
}

void mainResetRFChanTREmul( shell *psh )
{
    DVC_RESET( pEMMSysCtrl->pRFChanTREmul );
    cprintf( "Device RFChanTREmul Reset\n" );
}

//================
// RFChanTREmulSw
//================
void mainReadRFChanTREmulSwSlvReg( shell *psh )
{
    u32 reg = psh->uarg[ 0 ];
    u32 val = DVC_RD_SLV_REG( pEMMSysCtrl->pRFChanTREmulSw, reg );
    cprintf( " data = %u 0x%08x\n", toUint( val ), toUint( val ) );
}

void mainWriteRFChanTREmulSwSlvReg( shell *psh )
{
    u32 reg = psh->uarg[ 0 ];
    u32 val = psh->uarg[ 1 ];
    DVC_WR_SLV_REG( pEMMSysCtrl->pRFChanTREmulSw, reg, val );
    cprintf( "\n" );
}

void mainDumpRFChanTREmulSwSlvReg( shell *psh )
{
    u32 bBin = psh->uarg[ 0 ];
    DVC_DUMP_SLV_REG( pEMMSysCtrl->pRFChanTREmulSw, bBin );
}

void mainDumpRFChanTREmulSw( shell *psh )
{
    DVC_DUMP( pEMMSysCtrl->pRFChanTREmulSw );
}

void mainResetRFChanTREmulSw( shell *psh )
{
    DVC_RESET( pEMMSysCtrl->pRFChanTREmulSw );
    cprintf( "Device RFChanTREmulSw Reset\n" );
}

//=====================
// RFChanTREmulFreqCtr
//=====================
void mainReadRFChanTREmulFreqCtrSlvReg( shell *psh )
{
    u32 reg = psh->uarg[ 0 ];
    u32 val = DVC_RD_SLV_REG( pEMMSysCtrl->pRFChanTREmulFreqCtr, reg );
    cprintf( " data = %u 0x%08x\n", toUint( val ), toUint( val ) );
}

void mainWriteRFChanTREmulFreqCtrSlvReg( shell *psh )
{
    u32 reg = psh->uarg[ 0 ];
    u32 val = psh->uarg[ 1 ];
    DVC_WR_SLV_REG( pEMMSysCtrl->pRFChanTREmulFreqCtr, reg, val );
    cprintf( "\n" );
}

void mainDumpRFChanTREmulFreqCtrSlvReg( shell *psh )
{
    u32 bBin = psh->uarg[ 0 ];
    DVC_DUMP_SLV_REG( pEMMSysCtrl->pRFChanTREmulFreqCtr, bBin );
}

void mainDumpRFChanTREmulFreqCtr( shell *psh )
{
    DVC_DUMP( pEMMSysCtrl->pRFChanTREmulFreqCtr );
}

void mainResetRFChanTREmulFreqCtr( shell *psh )
{
    DVC_RESET( pEMMSysCtrl->pRFChanTREmulFreqCtr );
    cprintf( "Device RFChanTREmulFreqCtr Reset\n" );
}

void mainReadRFChanTREmulFreqCtr( shell *psh )
{
    double g = DVC_FUNC( pEMMSysCtrl->pRFChanTREmulFreqCtr, ReadFreqCtr );
    cprintf( "Frequency: %.6f Mhz\n", g );
}

//==================
// AD9364ClkFreqCtr
//==================
void mainReadAD9364ClkFreqCtrSlvReg( shell *psh )
{
    u32 reg = psh->uarg[ 0 ];
    u32 val = DVC_RD_SLV_REG( pEMMSysCtrl->pAD9364ClkFreqCtr, reg );
    cprintf( " data = %u 0x%08x\n", toUint( val ), toUint( val ) );
}

void mainWriteAD9364ClkFreqCtrSlvReg( shell *psh )
{
    u32 reg = psh->uarg[ 0 ];
    u32 val = psh->uarg[ 1 ];
    DVC_WR_SLV_REG( pEMMSysCtrl->pAD9364ClkFreqCtr, reg, val );
    cprintf( "\n" );
}

void mainDumpAD9364ClkFreqCtrSlvReg( shell *psh )
{
    u32 bBin = psh->uarg[ 0 ];
    DVC_DUMP_SLV_REG( pEMMSysCtrl->pAD9364ClkFreqCtr, bBin );
}

void mainDumpAD9364ClkFreqCtr( shell *psh )
{
    DVC_DUMP( pEMMSysCtrl->pAD9364ClkFreqCtr );
}

void mainResetAD9364ClkFreqCtr( shell *psh )
{
    DVC_RESET( pEMMSysCtrl->pAD9364ClkFreqCtr );
    cprintf( "Device AD9364ClkFreqCtr Reset\n" );
}

void mainReadAD9364ClkFreqCtr( shell *psh )
{
    double g = DVC_FUNC( pEMMSysCtrl->pAD9364ClkFreqCtr, ReadFreqCtr );
    cprintf( "Frequency: %.6f Mhz\n", g );
}

//==================
// AD9364DACFreqCtr
//==================
void mainReadAD9364DACFreqCtrSlvReg( shell *psh )
{
    u32 reg = psh->uarg[ 0 ];
    u32 val = DVC_RD_SLV_REG( pEMMSysCtrl->pAD9364DACFreqCtr, reg );
    cprintf( " data = %u 0x%08x\n", toUint( val ), toUint( val ) );
}

void mainWriteAD9364DACFreqCtrSlvReg( shell *psh )
{
    u32 reg = psh->uarg[ 0 ];
    u32 val = psh->uarg[ 1 ];
    DVC_WR_SLV_REG( pEMMSysCtrl->pAD9364DACFreqCtr, reg, val );
    cprintf( "\n" );
}

void mainDumpAD9364DACFreqCtrSlvReg( shell *psh )
{
    u32 bBin = psh->uarg[ 0 ];
    DVC_DUMP_SLV_REG( pEMMSysCtrl->pAD9364DACFreqCtr, bBin );
}

void mainDumpAD9364DACFreqCtr( shell *psh )
{
    DVC_DUMP( pEMMSysCtrl->pAD9364DACFreqCtr );
}

void mainResetAD9364DACFreqCtr( shell *psh )
{
    DVC_RESET( pEMMSysCtrl->pAD9364DACFreqCtr );
    cprintf( "Device AD9364DACFreqCtr Reset\n" );
}

void mainReadAD9364DACFreqCtr( shell *psh )
{
    double g = DVC_FUNC( pEMMSysCtrl->pAD9364DACFreqCtr, ReadFreqCtr );
    cprintf( "Frequency: %.6f Mhz\n", g );
}

//==================
// AD9364ADCFreqCtr
//==================
void mainReadAD9364ADCFreqCtrSlvReg( shell *psh )
{
    u32 reg = psh->uarg[ 0 ];
    u32 val = DVC_RD_SLV_REG( pEMMSysCtrl->pAD9364ADCFreqCtr, reg );
    cprintf( " data = %u 0x%08x\n", toUint( val ), toUint( val ) );
}

void mainWriteAD9364ADCFreqCtrSlvReg( shell *psh )
{
    u32 reg = psh->uarg[ 0 ];
    u32 val = psh->uarg[ 1 ];
    DVC_WR_SLV_REG( pEMMSysCtrl->pAD9364ADCFreqCtr, reg, val );
    cprintf( "\n" );
}

void mainDumpAD9364ADCFreqCtrSlvReg( shell *psh )
{
    u32 bBin = psh->uarg[ 0 ];
    DVC_DUMP_SLV_REG( pEMMSysCtrl->pAD9364ADCFreqCtr, bBin );
}

void mainDumpAD9364ADCFreqCtr( shell *psh )
{
    DVC_DUMP( pEMMSysCtrl->pAD9364ADCFreqCtr );
}

void mainResetAD9364ADCFreqCtr( shell *psh )
{
    DVC_RESET( pEMMSysCtrl->pAD9364ADCFreqCtr );
    cprintf( "Device AD9364ADCFreqCtr Reset\n" );
}

void mainReadAD9364ADCFreqCtr( shell *psh )
{
    double g = DVC_FUNC( pEMMSysCtrl->pAD9364ADCFreqCtr, ReadFreqCtr );
    cprintf( "Frequency: %.6f Mhz\n", g );
}

//=============
// RxNSDet
//=============
void mainReadRxNSDetSlvReg( shell *psh )
{
    u32 reg = psh->uarg[ 0 ];
    u32 val = DVC_RD_SLV_REG( pEMMSysCtrl->pRxNSDet, reg );
    cprintf( " data = %u 0x%08x\n", toUint( val ), toUint( val ) );
}

void mainWriteRxNSDetSlvReg( shell *psh )
{
    u32 reg = psh->uarg[ 0 ];
    u32 val = psh->uarg[ 1 ];
    DVC_WR_SLV_REG( pEMMSysCtrl->pRxNSDet, reg, val );
    cprintf( "\n" );
}

void mainDumpRxNSDetSlvReg( shell *psh )
{
    u32 bBin = psh->uarg[ 0 ];
    DVC_DUMP_SLV_REG( pEMMSysCtrl->pRxNSDet, bBin );
}

void mainDumpRxNSDet( shell *psh )
{
    DVC_DUMP( pEMMSysCtrl->pRxNSDet );
}

void mainResetRxNSDet( shell *psh )
{
    DVC_RESET( pEMMSysCtrl->pRxNSDet );
    cprintf( "Device RxNSDet Reset\n" );
}

//==============
// RxTddAxiSink
//==============
void mainReadRxTddAXISinkBram( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 data  = pEMMSysCtrl->pRxTddAXISink->ReadBram( pEMMSysCtrl->pRxTddAXISink, addr );
    cprintf( " data = %u 0x%08x\n", toUint( data ), toUint( data ) );
}

void mainWriteRxTddAXISinkBram( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 data  = psh->uarg[ 1 ];
    pEMMSysCtrl->pRxTddAXISink->WriteBram( pEMMSysCtrl->pRxTddAXISink, addr, data );
}

void mainInitRxTddAXISinkBram( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 nwds  = psh->uarg[ 1 ];
    u32 data  = psh->uarg[ 2 ];
    TIME_TFR( pEMMSysCtrl->pRxTddAXISink->InitBram( pEMMSysCtrl->pRxTddAXISink, addr, nwds, data ), nwds, 4 );
}

void mainDumpRxTddAXISinkBram( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 nbits = psh->uarg[ 1 ];
    pEMMSysCtrl->pRxTddAXISink->DumpBram( pEMMSysCtrl->pRxTddAXISink, addr, nbits );
}

void mainReadRxTddAXISinkSlvReg( shell *psh )
{
    u32 reg = psh->uarg[ 0 ];
    u32 val = DVC_RD_SLV_REG( pEMMSysCtrl->pRxTddAXISink, reg );
    cprintf( " data = %u 0x%08x\n", toUint( val ), toUint( val ) );
}

void mainWriteRxTddAXISinkSlvReg( shell *psh )
{
    u32 reg = psh->uarg[ 0 ];
    u32 val = psh->uarg[ 1 ];
    DVC_WR_SLV_REG( pEMMSysCtrl->pRxTddAXISink, reg, val );
    cprintf( "\n" );
}

void mainDumpRxTddAXISinkSlvReg( shell *psh )
{
    u32 bBin = psh->uarg[ 0 ];
    DVC_DUMP_SLV_REG( pEMMSysCtrl->pRxTddAXISink, bBin );
}

void mainDumpRxTddAXISink( shell *psh )
{
    DVC_DUMP( pEMMSysCtrl->pRxTddAXISink );
}

void mainResetRxTddAXISink( shell *psh )
{
    DVC_RESET( pEMMSysCtrl->pRxTddAXISink );
    cprintf( "Device RxTddAXISink Reset\n" );
}

//==============
// RxPrsAxiSink
//==============
void mainReadRxPrsAXISinkBram( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 data  = pEMMSysCtrl->pRxPrsAXISink->ReadBram( pEMMSysCtrl->pRxPrsAXISink, addr );
    cprintf( " data = %u 0x%08x\n", toUint( data ), toUint( data ) );
}

void mainWriteRxPrsAXISinkBram( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 data  = psh->uarg[ 1 ];
    pEMMSysCtrl->pRxPrsAXISink->WriteBram( pEMMSysCtrl->pRxPrsAXISink, addr, data );
}

void mainInitRxPrsAXISinkBram( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 nwds  = psh->uarg[ 1 ];
    u32 data  = psh->uarg[ 2 ];
    TIME_TFR( pEMMSysCtrl->pRxPrsAXISink->InitBram( pEMMSysCtrl->pRxPrsAXISink, addr, nwds, data ), nwds, 4 );
}

void mainDumpRxPrsAXISinkBram( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 nbits = psh->uarg[ 1 ];
    pEMMSysCtrl->pRxPrsAXISink->DumpBram( pEMMSysCtrl->pRxPrsAXISink, addr, nbits );
}

void mainReadRxPrsAXISinkSlvReg( shell *psh )
{
    u32 reg = psh->uarg[ 0 ];
    u32 val = DVC_RD_SLV_REG( pEMMSysCtrl->pRxPrsAXISink, reg );
    cprintf( " data = %u 0x%08x\n", toUint( val ), toUint( val ) );
}

void mainWriteRxPrsAXISinkSlvReg( shell *psh )
{
    u32 reg = psh->uarg[ 0 ];
    u32 val = psh->uarg[ 1 ];
    DVC_WR_SLV_REG( pEMMSysCtrl->pRxPrsAXISink, reg, val );
    cprintf( "\n" );
}

void mainDumpRxPrsAXISinkSlvReg( shell *psh )
{
    u32 bBin = psh->uarg[ 0 ];
    DVC_DUMP_SLV_REG( pEMMSysCtrl->pRxPrsAXISink, bBin );
}

void mainDumpRxPrsAXISink( shell *psh )
{
    DVC_DUMP( pEMMSysCtrl->pRxPrsAXISink );
}

void mainResetRxPrsAXISink( shell *psh )
{
    DVC_RESET( pEMMSysCtrl->pRxPrsAXISink );
    cprintf( "Device RxPrsAXISink Reset\n" );
}

//==============
// RxDtaAxiSink
//==============
void mainReadRxDtaAXISinkBram( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 data  = pEMMSysCtrl->pRxDtaAXISink->ReadBram( pEMMSysCtrl->pRxDtaAXISink, addr );
    cprintf( " data = %u 0x%08x\n", toUint( data ), toUint( data ) );
}

void mainWriteRxDtaAXISinkBram( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 data  = psh->uarg[ 1 ];
    pEMMSysCtrl->pRxDtaAXISink->WriteBram( pEMMSysCtrl->pRxDtaAXISink, addr, data );
}

void mainInitRxDtaAXISinkBram( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 nwds  = psh->uarg[ 1 ];
    u32 data  = psh->uarg[ 2 ];
    TIME_TFR( pEMMSysCtrl->pRxDtaAXISink->InitBram( pEMMSysCtrl->pRxDtaAXISink, addr, nwds, data ), nwds, 4 );
}

void mainDumpRxDtaAXISinkBram( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 nbits = psh->uarg[ 1 ];
    pEMMSysCtrl->pRxDtaAXISink->DumpBram( pEMMSysCtrl->pRxDtaAXISink, addr, nbits );
}

void mainReadRxDtaAXISinkSlvReg( shell *psh )
{
    u32 reg = psh->uarg[ 0 ];
    u32 val = DVC_RD_SLV_REG( pEMMSysCtrl->pRxDtaAXISink, reg );
    cprintf( " data = %u 0x%08x\n", toUint( val ), toUint( val ) );
}

void mainWriteRxDtaAXISinkSlvReg( shell *psh )
{
    u32 reg = psh->uarg[ 0 ];
    u32 val = psh->uarg[ 1 ];
    DVC_WR_SLV_REG( pEMMSysCtrl->pRxDtaAXISink, reg, val );
    cprintf( "\n" );
}

void mainDumpRxDtaAXISinkSlvReg( shell *psh )
{
    u32 bBin = psh->uarg[ 0 ];
    DVC_DUMP_SLV_REG( pEMMSysCtrl->pRxDtaAXISink, bBin );
}

void mainDumpRxDtaAXISink( shell *psh )
{
    DVC_DUMP( pEMMSysCtrl->pRxDtaAXISink );
}

void mainResetRxDtaAXISink( shell *psh )
{
    DVC_RESET( pEMMSysCtrl->pRxDtaAXISink );
    cprintf( "Device RxDtaAXISink Reset\n" );
}

//=================
// RxNsdA00AxiSink
//=================
void mainReadRxNsdA00AXISinkBram( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 data  = pEMMSysCtrl->pRxNsdA00AXISink->ReadBram( pEMMSysCtrl->pRxNsdA00AXISink, addr );
    cprintf( " data = %u 0x%08x\n", toUint( data ), toUint( data ) );
}

void mainWriteRxNsdA00AXISinkBram( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 data  = psh->uarg[ 1 ];
    pEMMSysCtrl->pRxNsdA00AXISink->WriteBram( pEMMSysCtrl->pRxNsdA00AXISink, addr, data );
}

void mainInitRxNsdA00AXISinkBram( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 nwds  = psh->uarg[ 1 ];
    u32 data  = psh->uarg[ 2 ];
    TIME_TFR( pEMMSysCtrl->pRxNsdA00AXISink->InitBram( pEMMSysCtrl->pRxNsdA00AXISink, addr, nwds, data ), nwds, 4 );
}

void mainDumpRxNsdA00AXISinkBram( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 nbits = psh->uarg[ 1 ];
    pEMMSysCtrl->pRxNsdA00AXISink->DumpBram( pEMMSysCtrl->pRxNsdA00AXISink, addr, nbits );
}

void mainReadRxNsdA00AXISinkSlvReg( shell *psh )
{
    u32 reg = psh->uarg[ 0 ];
    u32 val = DVC_RD_SLV_REG( pEMMSysCtrl->pRxNsdA00AXISink, reg );
    cprintf( " data = %u 0x%08x\n", toUint( val ), toUint( val ) );
}

void mainWriteRxNsdA00AXISinkSlvReg( shell *psh )
{
    u32 reg = psh->uarg[ 0 ];
    u32 val = psh->uarg[ 1 ];
    DVC_WR_SLV_REG( pEMMSysCtrl->pRxNsdA00AXISink, reg, val );
    cprintf( "\n" );
}

void mainDumpRxNsdA00AXISinkSlvReg( shell *psh )
{
    u32 bBin = psh->uarg[ 0 ];
    DVC_DUMP_SLV_REG( pEMMSysCtrl->pRxNsdA00AXISink, bBin );
}

void mainDumpRxNsdA00AXISink( shell *psh )
{
    DVC_DUMP( pEMMSysCtrl->pRxNsdA00AXISink );
}

void mainResetRxNsdA00AXISink( shell *psh )
{
    DVC_RESET( pEMMSysCtrl->pRxNsdA00AXISink );
    cprintf( "Device RxNsdA00AXISink Reset\n" );
}

//=================
// RxNsdA01AxiSink
//=================
void mainReadRxNsdA01AXISinkBram( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 data  = pEMMSysCtrl->pRxNsdA01AXISink->ReadBram( pEMMSysCtrl->pRxNsdA01AXISink, addr );
    cprintf( " data = %u 0x%08x\n", toUint( data ), toUint( data ) );
}

void mainWriteRxNsdA01AXISinkBram( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 data  = psh->uarg[ 1 ];
    pEMMSysCtrl->pRxNsdA01AXISink->WriteBram( pEMMSysCtrl->pRxNsdA01AXISink, addr, data );
}

void mainInitRxNsdA01AXISinkBram( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 nwds  = psh->uarg[ 1 ];
    u32 data  = psh->uarg[ 2 ];
    TIME_TFR( pEMMSysCtrl->pRxNsdA01AXISink->InitBram( pEMMSysCtrl->pRxNsdA01AXISink, addr, nwds, data ), nwds, 4 );
}

void mainDumpRxNsdA01AXISinkBram( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 nbits = psh->uarg[ 1 ];
    pEMMSysCtrl->pRxNsdA01AXISink->DumpBram( pEMMSysCtrl->pRxNsdA01AXISink, addr, nbits );
}

void mainReadRxNsdA01AXISinkSlvReg( shell *psh )
{
    u32 reg = psh->uarg[ 0 ];
    u32 val = DVC_RD_SLV_REG( pEMMSysCtrl->pRxNsdA01AXISink, reg );
    cprintf( " data = %u 0x%08x\n", toUint( val ), toUint( val ) );
}

void mainWriteRxNsdA01AXISinkSlvReg( shell *psh )
{
    u32 reg = psh->uarg[ 0 ];
    u32 val = psh->uarg[ 1 ];
    DVC_WR_SLV_REG( pEMMSysCtrl->pRxNsdA01AXISink, reg, val );
    cprintf( "\n" );
}

void mainDumpRxNsdA01AXISinkSlvReg( shell *psh )
{
    u32 bBin = psh->uarg[ 0 ];
    DVC_DUMP_SLV_REG( pEMMSysCtrl->pRxNsdA01AXISink, bBin );
}

void mainDumpRxNsdA01AXISink( shell *psh )
{
    DVC_DUMP( pEMMSysCtrl->pRxNsdA01AXISink );
}

void mainResetRxNsdA01AXISink( shell *psh )
{
    DVC_RESET( pEMMSysCtrl->pRxNsdA01AXISink );
    cprintf( "Device RxNsdA01AXISink Reset\n" );
}

//=================
// RxNsdA15AxiSink
//=================
void mainReadRxNsdA15AXISinkBram( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 data  = pEMMSysCtrl->pRxNsdA15AXISink->ReadBram( pEMMSysCtrl->pRxNsdA15AXISink, addr );
    cprintf( " data = %u 0x%08x\n", toUint( data ), toUint( data ) );
}

void mainWriteRxNsdA15AXISinkBram( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 data  = psh->uarg[ 1 ];
    pEMMSysCtrl->pRxNsdA15AXISink->WriteBram( pEMMSysCtrl->pRxNsdA15AXISink, addr, data );
}

void mainInitRxNsdA15AXISinkBram( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 nwds  = psh->uarg[ 1 ];
    u32 data  = psh->uarg[ 2 ];
    TIME_TFR( pEMMSysCtrl->pRxNsdA15AXISink->InitBram( pEMMSysCtrl->pRxNsdA15AXISink, addr, nwds, data ), nwds, 4 );
}

void mainDumpRxNsdA15AXISinkBram( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 nbits = psh->uarg[ 1 ];
    pEMMSysCtrl->pRxNsdA15AXISink->DumpBram( pEMMSysCtrl->pRxNsdA15AXISink, addr, nbits );
}

void mainReadRxNsdA15AXISinkSlvReg( shell *psh )
{
    u32 reg = psh->uarg[ 0 ];
    u32 val = DVC_RD_SLV_REG( pEMMSysCtrl->pRxNsdA15AXISink, reg );
    cprintf( " data = %u 0x%08x\n", toUint( val ), toUint( val ) );
}

void mainWriteRxNsdA15AXISinkSlvReg( shell *psh )
{
    u32 reg = psh->uarg[ 0 ];
    u32 val = psh->uarg[ 1 ];
    DVC_WR_SLV_REG( pEMMSysCtrl->pRxNsdA15AXISink, reg, val );
    cprintf( "\n" );
}

void mainDumpRxNsdA15AXISinkSlvReg( shell *psh )
{
    u32 bBin = psh->uarg[ 0 ];
    DVC_DUMP_SLV_REG( pEMMSysCtrl->pRxNsdA15AXISink, bBin );
}

void mainDumpRxNsdA15AXISink( shell *psh )
{
    DVC_DUMP( pEMMSysCtrl->pRxNsdA15AXISink );
}

void mainResetRxNsdA15AXISink( shell *psh )
{
    DVC_RESET( pEMMSysCtrl->pRxNsdA15AXISink );
    cprintf( "Device RxNsdA15AXISink Reset\n" );
}

//====================
// RxNsdCordicAXISink
//====================
void mainReadRxNsdCordicAXISinkBram( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 data  = pEMMSysCtrl->pRxNsdCordicAXISink->ReadBram( pEMMSysCtrl->pRxNsdCordicAXISink, addr );
    cprintf( " data = %u 0x%08x\n", toUint( data ), toUint( data ) );
}

void mainWriteRxNsdCordicAXISinkBram( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 data  = psh->uarg[ 1 ];
    pEMMSysCtrl->pRxNsdCordicAXISink->WriteBram( pEMMSysCtrl->pRxNsdCordicAXISink, addr, data );
}

void mainInitRxNsdCordicAXISinkBram( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 nwds  = psh->uarg[ 1 ];
    u32 data  = psh->uarg[ 2 ];
    TIME_TFR( pEMMSysCtrl->pRxNsdCordicAXISink->InitBram( pEMMSysCtrl->pRxNsdCordicAXISink, addr, nwds, data ), nwds, 4 );
}

void mainDumpRxNsdCordicAXISinkBram( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 nbits = psh->uarg[ 1 ];
    pEMMSysCtrl->pRxNsdCordicAXISink->DumpBram( pEMMSysCtrl->pRxNsdCordicAXISink, addr, nbits );
}

void mainReadRxNsdCordicAXISinkSlvReg( shell *psh )
{
    u32 reg = psh->uarg[ 0 ];
    u32 val = DVC_RD_SLV_REG( pEMMSysCtrl->pRxNsdCordicAXISink, reg );
    cprintf( " data = %u 0x%08x\n", toUint( val ), toUint( val ) );
}

void mainWriteRxNsdCordicAXISinkSlvReg( shell *psh )
{
    u32 reg = psh->uarg[ 0 ];
    u32 val = psh->uarg[ 1 ];
    DVC_WR_SLV_REG( pEMMSysCtrl->pRxNsdCordicAXISink, reg, val );
    cprintf( "\n" );
}

void mainDumpRxNsdCordicAXISinkSlvReg( shell *psh )
{
    u32 bBin = psh->uarg[ 0 ];
    DVC_DUMP_SLV_REG( pEMMSysCtrl->pRxNsdCordicAXISink, bBin );
}

void mainDumpRxNsdCordicAXISink( shell *psh )
{
    DVC_DUMP( pEMMSysCtrl->pRxNsdCordicAXISink );
}

void mainResetRxNsdCordicAXISink( shell *psh )
{
    DVC_RESET( pEMMSysCtrl->pRxNsdCordicAXISink );
    cprintf( "Device RxNsdCordicAXISink Reset\n" );
}

//=================
// RxNsdFFTAXISink
//=================
void mainReadRxNsdFFTAXISinkBram( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 data  = pEMMSysCtrl->pRxNsdFFTAXISink->ReadBram( pEMMSysCtrl->pRxNsdFFTAXISink, addr );
    cprintf( " data = %u 0x%08x\n", toUint( data ), toUint( data ) );
}

void mainWriteRxNsdFFTAXISinkBram( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 data  = psh->uarg[ 1 ];
    pEMMSysCtrl->pRxNsdFFTAXISink->WriteBram( pEMMSysCtrl->pRxNsdFFTAXISink, addr, data );
}

void mainInitRxNsdFFTAXISinkBram( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 nwds  = psh->uarg[ 1 ];
    u32 data  = psh->uarg[ 2 ];
    TIME_TFR( pEMMSysCtrl->pRxNsdFFTAXISink->InitBram( pEMMSysCtrl->pRxNsdFFTAXISink, addr, nwds, data ), nwds, 4 );
}

void mainDumpRxNsdFFTAXISinkBram( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 nbits = psh->uarg[ 1 ];
    pEMMSysCtrl->pRxNsdFFTAXISink->DumpBram( pEMMSysCtrl->pRxNsdFFTAXISink, addr, nbits );
}

void mainReadRxNsdFFTAXISinkSlvReg( shell *psh )
{
    u32 reg = psh->uarg[ 0 ];
    u32 val = DVC_RD_SLV_REG( pEMMSysCtrl->pRxNsdFFTAXISink, reg );
    cprintf( " data = %u 0x%08x\n", toUint( val ), toUint( val ) );
}

void mainWriteRxNsdFFTAXISinkSlvReg( shell *psh )
{
    u32 reg = psh->uarg[ 0 ];
    u32 val = psh->uarg[ 1 ];
    DVC_WR_SLV_REG( pEMMSysCtrl->pRxNsdFFTAXISink, reg, val );
    cprintf( "\n" );
}

void mainDumpRxNsdFFTAXISinkSlvReg( shell *psh )
{
    u32 bBin = psh->uarg[ 0 ];
    DVC_DUMP_SLV_REG( pEMMSysCtrl->pRxNsdFFTAXISink, bBin );
}

void mainDumpRxNsdFFTAXISink( shell *psh )
{
    DVC_DUMP( pEMMSysCtrl->pRxNsdFFTAXISink );
}

void mainResetRxNsdFFTAXISink( shell *psh )
{
    DVC_RESET( pEMMSysCtrl->pRxNsdFFTAXISink );
    cprintf( "Device RxNsdFFTAXISink Reset\n" );
}

//====================
// RxNsdBcast00AXISink
//====================
void mainReadRxNsdB00AXISinkBram( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 data  = pEMMSysCtrl->pRxNsdB00AXISink->ReadBram( pEMMSysCtrl->pRxNsdB00AXISink, addr );
    cprintf( " data = %u 0x%08x\n", toUint( data ), toUint( data ) );
}

void mainWriteRxNsdB00AXISinkBram( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 data  = psh->uarg[ 1 ];
    pEMMSysCtrl->pRxNsdB00AXISink->WriteBram( pEMMSysCtrl->pRxNsdB00AXISink, addr, data );
}

void mainInitRxNsdB00AXISinkBram( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 nwds  = psh->uarg[ 1 ];
    u32 data  = psh->uarg[ 2 ];
    TIME_TFR( pEMMSysCtrl->pRxNsdB00AXISink->InitBram( pEMMSysCtrl->pRxNsdB00AXISink, addr, nwds, data ), nwds, 4 );
}

void mainDumpRxNsdB00AXISinkBram( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 nbits = psh->uarg[ 1 ];
    pEMMSysCtrl->pRxNsdB00AXISink->DumpBram( pEMMSysCtrl->pRxNsdB00AXISink, addr, nbits );
}

void mainReadRxNsdB00AXISinkSlvReg( shell *psh )
{
    u32 reg = psh->uarg[ 0 ];
    u32 val = DVC_RD_SLV_REG( pEMMSysCtrl->pRxNsdB00AXISink, reg );
    cprintf( " data = %u 0x%08x\n", toUint( val ), toUint( val ) );
}

void mainWriteRxNsdB00AXISinkSlvReg( shell *psh )
{
    u32 reg = psh->uarg[ 0 ];
    u32 val = psh->uarg[ 1 ];
    DVC_WR_SLV_REG( pEMMSysCtrl->pRxNsdB00AXISink, reg, val );
    cprintf( "\n" );
}

void mainDumpRxNsdB00AXISinkSlvReg( shell *psh )
{
    u32 bBin = psh->uarg[ 0 ];
    DVC_DUMP_SLV_REG( pEMMSysCtrl->pRxNsdB00AXISink, bBin );
}

void mainDumpRxNsdB00AXISink( shell *psh )
{
    DVC_DUMP( pEMMSysCtrl->pRxNsdB00AXISink );
}

void mainResetRxNsdB00AXISink( shell *psh )
{
    DVC_RESET( pEMMSysCtrl->pRxNsdB00AXISink );
    cprintf( "Device RxNsdB00AXISink Reset\n" );
}

//====================
// RxNsdBcast01AXISink
//====================
void mainReadRxNsdB01AXISinkBram( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 data  = pEMMSysCtrl->pRxNsdB01AXISink->ReadBram( pEMMSysCtrl->pRxNsdB01AXISink, addr );
    cprintf( " data = %u 0x%08x\n", toUint( data ), toUint( data ) );
}

void mainWriteRxNsdB01AXISinkBram( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 data  = psh->uarg[ 1 ];
    pEMMSysCtrl->pRxNsdB01AXISink->WriteBram( pEMMSysCtrl->pRxNsdB01AXISink, addr, data );
}

void mainInitRxNsdB01AXISinkBram( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 nwds  = psh->uarg[ 1 ];
    u32 data  = psh->uarg[ 2 ];
    TIME_TFR( pEMMSysCtrl->pRxNsdB01AXISink->InitBram( pEMMSysCtrl->pRxNsdB01AXISink, addr, nwds, data ), nwds, 4 );
}

void mainDumpRxNsdB01AXISinkBram( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 nbits = psh->uarg[ 1 ];
    pEMMSysCtrl->pRxNsdB01AXISink->DumpBram( pEMMSysCtrl->pRxNsdB01AXISink, addr, nbits );
}

void mainReadRxNsdB01AXISinkSlvReg( shell *psh )
{
    u32 reg = psh->uarg[ 0 ];
    u32 val = DVC_RD_SLV_REG( pEMMSysCtrl->pRxNsdB01AXISink, reg );
    cprintf( " data = %u 0x%08x\n", toUint( val ), toUint( val ) );
}

void mainWriteRxNsdB01AXISinkSlvReg( shell *psh )
{
    u32 reg = psh->uarg[ 0 ];
    u32 val = psh->uarg[ 1 ];
    DVC_WR_SLV_REG( pEMMSysCtrl->pRxNsdB01AXISink, reg, val );
    cprintf( "\n" );
}

void mainDumpRxNsdB01AXISinkSlvReg( shell *psh )
{
    u32 bBin = psh->uarg[ 0 ];
    DVC_DUMP_SLV_REG( pEMMSysCtrl->pRxNsdB01AXISink, bBin );
}

void mainDumpRxNsdB01AXISink( shell *psh )
{
    DVC_DUMP( pEMMSysCtrl->pRxNsdB01AXISink );
}

void mainResetRxNsdB01AXISink( shell *psh )
{
    DVC_RESET( pEMMSysCtrl->pRxNsdB01AXISink );
    cprintf( "Device RxNsdB01AXISink Reset\n" );
}

//====================
// RxNsdBcast15AXISink
//====================
void mainReadRxNsdB15AXISinkBram( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 data  = pEMMSysCtrl->pRxNsdB15AXISink->ReadBram( pEMMSysCtrl->pRxNsdB15AXISink, addr );
    cprintf( " data = %u 0x%08x\n", toUint( data ), toUint( data ) );
}

void mainWriteRxNsdB15AXISinkBram( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 data  = psh->uarg[ 1 ];
    pEMMSysCtrl->pRxNsdB15AXISink->WriteBram( pEMMSysCtrl->pRxNsdB15AXISink, addr, data );
}

void mainInitRxNsdB15AXISinkBram( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 nwds  = psh->uarg[ 1 ];
    u32 data  = psh->uarg[ 2 ];
    TIME_TFR( pEMMSysCtrl->pRxNsdB15AXISink->InitBram( pEMMSysCtrl->pRxNsdB15AXISink, addr, nwds, data ), nwds, 4 );
}

void mainDumpRxNsdB15AXISinkBram( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 nbits = psh->uarg[ 1 ];
    pEMMSysCtrl->pRxNsdB15AXISink->DumpBram( pEMMSysCtrl->pRxNsdB15AXISink, addr, nbits );
}

void mainReadRxNsdB15AXISinkSlvReg( shell *psh )
{
    u32 reg = psh->uarg[ 0 ];
    u32 val = DVC_RD_SLV_REG( pEMMSysCtrl->pRxNsdB15AXISink, reg );
    cprintf( " data = %u 0x%08x\n", toUint( val ), toUint( val ) );
}

void mainWriteRxNsdB15AXISinkSlvReg( shell *psh )
{
    u32 reg = psh->uarg[ 0 ];
    u32 val = psh->uarg[ 1 ];
    DVC_WR_SLV_REG( pEMMSysCtrl->pRxNsdB15AXISink, reg, val );
    cprintf( "\n" );
}

void mainDumpRxNsdB15AXISinkSlvReg( shell *psh )
{
    u32 bBin = psh->uarg[ 0 ];
    DVC_DUMP_SLV_REG( pEMMSysCtrl->pRxNsdB15AXISink, bBin );
}

void mainDumpRxNsdB15AXISink( shell *psh )
{
    DVC_DUMP( pEMMSysCtrl->pRxNsdB15AXISink );
}

void mainResetRxNsdB15AXISink( shell *psh )
{
    DVC_RESET( pEMMSysCtrl->pRxNsdB15AXISink );
    cprintf( "Device RxNsdB15AXISink Reset\n" );
}

//=======
// PrsCIR
//=======
void mainReadRxPrsCIRBramFFT( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 data  = pEMMSysCtrl->pRxPrsCIR->ReadBramFFT( pEMMSysCtrl->pRxPrsCIR, addr );
    cprintf( " data = %u 0x%08x\n", toUint( data ), toUint( data ) );
}

void mainWriteRxPrsCIRBramFFT( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 data  = psh->uarg[ 1 ];
    pEMMSysCtrl->pRxPrsCIR->WriteBramFFT( pEMMSysCtrl->pRxPrsCIR, addr, data );
}

void mainInitRxPrsCIRBramFFT( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 nwds  = psh->uarg[ 1 ];
    u32 data  = psh->uarg[ 2 ];
    TIME_TFR( pEMMSysCtrl->pRxPrsCIR->InitBramFFT( pEMMSysCtrl->pRxPrsCIR, addr, nwds, data ), nwds, 4 );
}

void mainDumpRxPrsCIRBramFFT( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 nbits = psh->uarg[ 1 ];
    pEMMSysCtrl->pRxPrsCIR->DumpBramFFT( pEMMSysCtrl->pRxPrsCIR, addr, nbits );
}

void mainReadRxPrsCIRBramABS( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 data  = pEMMSysCtrl->pRxPrsCIR->ReadBramABS( pEMMSysCtrl->pRxPrsCIR, addr );
    cprintf( " data = %u 0x%08x\n", toUint( data ), toUint( data ) );
}

void mainWriteRxPrsCIRBramABS( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 data  = psh->uarg[ 1 ];
    pEMMSysCtrl->pRxPrsCIR->WriteBramABS( pEMMSysCtrl->pRxPrsCIR, addr, data );
}

void mainInitRxPrsCIRBramABS( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 nwds  = psh->uarg[ 1 ];
    u32 data  = psh->uarg[ 2 ];
    TIME_TFR( pEMMSysCtrl->pRxPrsCIR->InitBramABS( pEMMSysCtrl->pRxPrsCIR, addr, nwds, data ), nwds, 4 );
}

void mainDumpRxPrsCIRBramABS( shell *psh )
{
    u32 addr  = psh->uarg[ 0 ];
    u32 nbits = psh->uarg[ 1 ];
    pEMMSysCtrl->pRxPrsCIR->DumpBramABS( pEMMSysCtrl->pRxPrsCIR, addr, nbits );
}

void mainReadRxPrsCIRSlvReg( shell *psh )
{
    u32 reg = psh->uarg[ 0 ];
    u32 val = DVC_RD_SLV_REG( pEMMSysCtrl->pRxPrsCIR, reg );
    cprintf( " data = %u 0x%08x\n", toUint( val ), toUint( val ) );
}

void mainWriteRxPrsCIRSlvReg( shell *psh )
{
    u32 reg = psh->uarg[ 0 ];
    u32 val = psh->uarg[ 1 ];
    DVC_WR_SLV_REG( pEMMSysCtrl->pRxPrsCIR, reg, val );
    cprintf( "\n" );
}

void mainDumpRxPrsCIRSlvReg( shell *psh )
{
    u32 bBin = psh->uarg[ 0 ];
    DVC_DUMP_SLV_REG( pEMMSysCtrl->pRxPrsCIR, bBin );
}

void mainDumpRxPrsCIR( shell *psh )
{
    DVC_DUMP( pEMMSysCtrl->pRxPrsCIR );
}

void mainResetRxPrsCIR( shell *psh )
{
    DVC_RESET( pEMMSysCtrl->pRxPrsCIR );
    cprintf( "Device RxPrsCIR Reset\n" );
}

void mainRxPrsCIRCalc( shell *psh )
{
    DVC_FUNC( pEMMSysCtrl->pRxPrsCIR, Calc );
    DVC_FUNC( pEMMSysCtrl->pRxPrsCIR, Dump );
}

void mainRxPrsCIRRecalc( shell *psh )
{
    DVC_FUNC( pEMMSysCtrl->pRxPrsCIR, Recalc );
    DVC_FUNC( pEMMSysCtrl->pRxPrsCIR, Dump );
}

void mainRxPrsCIRDumpCIRTable( shell *psh )
{
    DVC_FUNC( pEMMSysCtrl->pRxPrsCIR, DumpCIRTable );
}

void mainRxPrsCIRDumpTxPRSFd( shell *psh )
{
    DVC_FUNC( pEMMSysCtrl->pRxPrsCIR, DumpTxPRSFd );
}

void mainRxPrsCIRDumpRxPRSTd( shell *psh )
{
    DVC_FUNC( pEMMSysCtrl->pRxPrsCIR, DumpRxPRSTd );
}

void mainRxPrsCIRDumpRxPRSTdCorFT( shell *psh )
{
    DVC_FUNC( pEMMSysCtrl->pRxPrsCIR, DumpRxPRSTdCorFT );
}

void mainRxPrsCIRDumpRxPRSFd( shell *psh )
{
    DVC_FUNC( pEMMSysCtrl->pRxPrsCIR, DumpRxPRSFd );
}

void mainRxPrsCIRDumpRxPRSTdCorDeltaT( shell *psh )
{
    DVC_FUNC( pEMMSysCtrl->pRxPrsCIR, DumpRxPRSTdCorDeltaT );
}

void mainRxPrsCIRDumpRxPRSTdCordeltafp( shell *psh )
{
    DVC_FUNC( pEMMSysCtrl->pRxPrsCIR, DumpRxPRSTdCordeltafp );
}

void mainRxPrsCIRDumpRxPRSFdCordeltafp( shell *psh )
{
    DVC_FUNC( pEMMSysCtrl->pRxPrsCIR, DumpRxPRSFdCordeltafp );
}

void mainRxPrsCIRDumpRxPRSFdCordeltatp( shell *psh )
{
    DVC_FUNC( pEMMSysCtrl->pRxPrsCIR, DumpRxPRSFdCordeltatp );
}


//=======
// PrsBPA
//=======
void mainReadRxPrsBPASlvReg( shell *psh )
{
    u32 reg = psh->uarg[ 0 ];
    u32 val = DVC_RD_SLV_REG( pEMMSysCtrl->pRxPrsBPA, reg );
    cprintf( " data = %u 0x%08x\n", toUint( val ), toUint( val ) );
}

void mainWriteRxPrsBPASlvReg( shell *psh )
{
    u32 reg = psh->uarg[ 0 ];
    u32 val = psh->uarg[ 1 ];
    DVC_WR_SLV_REG( pEMMSysCtrl->pRxPrsBPA, reg, val );
    cprintf( "\n" );
}

void mainDumpRxPrsBPASlvReg( shell *psh )
{
    u32 bBin = psh->uarg[ 0 ];
    DVC_DUMP_SLV_REG( pEMMSysCtrl->pRxPrsBPA, bBin );
}

void mainDumpRxPrsBPA( shell *psh )
{
    DVC_DUMP( pEMMSysCtrl->pRxPrsBPA );
}

void mainResetRxPrsBPA( shell *psh )
{
    DVC_RESET( pEMMSysCtrl->pRxPrsBPA );
    cprintf( "Device RxPrsBPA Reset\n" );
}

void mainRxPrsBPACalcDelta( shell *psh )
{
    DVC_FUNC( pEMMSysCtrl->pRxPrsCIR, Calc );
    DVC_FUNC( pEMMSysCtrl->pRxPrsBPA, CalcDeltas );
}

void mainRxPrsBPACalcDelta2( shell *psh )
{
    DVC_FUNC( pEMMSysCtrl->pRxPrsCIR, Calc );
    DVC_FUNC( pEMMSysCtrl->pRxPrsCIR, Recalc );
    DVC_FUNC( pEMMSysCtrl->pRxPrsBPA, CalcDeltas );
}

void mainRxPrsBPADumpCorSymbFd( shell *psh )
{
    DVC_FUNC( pEMMSysCtrl->pRxPrsBPA, DumpCorrSymbFd );
}
void mainRxDataBPADumpCorSymbFd( shell *psh )
{
    DVC_FUNC( pEMMSysCtrl->pRxPrsBPA, DumpCorrDataSymbFd );
}
void mainRxDataBPADumpCorSymbTd( shell *psh )
{
    DVC_FUNC( pEMMSysCtrl->pRxPrsBPA, DumpCorrDataSymbTd );
}

void mainRxPrsBPADumpCorSymbTd( shell *psh )
{
    DVC_FUNC( pEMMSysCtrl->pRxPrsBPA, DumpCorrSymbTd );
}

void mainRxDataSymbolDumpRawFd( shell *psh )
{
    DVC_FUNC( pEMMSysCtrl->pRxPrsBPA, RxDataSymbolDumpRawFd );
}

void mainRxDataSymbolDumpFCFd( shell *psh )
{
    DVC_FUNC( pEMMSysCtrl->pRxPrsBPA, RxDataSymbolDumpFCFd );
}

void mainRxDataSymbolDumpFTCFd( shell *psh )
{
    DVC_FUNC( pEMMSysCtrl->pRxPrsBPA, RxDataSymbolDumpFTCFd );
}



//===========================================
//////// End of script-generated code ///////
//===========================================
/************************************************************
 DAC Simulator Clock Generator
 ************************************************************/

void mainSetDacFreq( shell *psh )
{
	u32 freq = psh->uarg[ 0 ];
	if ( clkWizSetFreq( &clkWizDac, freq ) )
		cprintf( "\n%s Frequency set to %u MHz\n", clkWizDac.name, ( unsigned int ) freq );
}

