/*
 * SysData.c
 *
 *  Created on: 11 Aug 2018
 *      Author: JULIAN MORTIMER
 */
#include "include.h"

void SysDataInit( SysDvce *pd,
		u32 RxDataCtrReg,
		u32 RxSymbCtrReg,
		u32 RxFrameCtrReg,
		u32 RxBlockCtrReg,
		u32 TxDataCtrReg,
		u32 TxSymbCtrReg,
		u32 TxFrameCtrReg,
		u32 TxBlockCtrReg,
		u32 RxTVCtrReg,
		u32 TxTVCtrReg
		);

static void Update( SysData *ps );
static void Dump(   SysData *ps );
static void Clear(  SysData *ps );

void Clear( SysData *ps )
{
	ps->RxDataCtr  = 0;
	ps->RxSymbCtr  = 0;
	ps->RxFrameCtr = 0;
	ps->RxBlockCtr = 0;
	ps->TxDataCtr  = 0;
	ps->TxSymbCtr  = 0;
	ps->TxFrameCtr = 0;
	ps->TxBlockCtr = 0;
    ps->RxTVCtr    = 0;
    ps->TxTVCtr    = 0;
}

void Dump( SysData *ps )
{
    ps->Update( ps );
    _Bool nl = false;
	if( ps->RxDataCtrReg  != 0 ) { nl = true; printf( "\nRx data ctr      : %u",   toUint( ps->RxDataCtr  ) );  fflush( stdout ); }
	if( ps->RxSymbCtrReg  != 0 ) { nl = true; printf( "\nRx symbol ctr    : %u",   toUint( ps->RxSymbCtr  ) );  fflush( stdout ); }
	if( ps->RxFrameCtrReg != 0 ) { nl = true; printf( "\nRx frame ctr     : %u",   toUint( ps->RxFrameCtr ) );  fflush( stdout ); }
	if( ps->RxBlockCtrReg != 0 ) { nl = true; printf( "\nRx block ctr     : %u",   toUint( ps->RxBlockCtr ) );  fflush( stdout ); }
	if( ps->TxDataCtrReg  != 0 ) { nl = true; printf( "\nTx data ctr      : %u",   toUint( ps->TxDataCtr  ) );  fflush( stdout ); }
	if( ps->TxSymbCtrReg  != 0 ) { nl = true; printf( "\nTx symbol ctr    : %u",   toUint( ps->TxSymbCtr  ) );  fflush( stdout ); }
	if( ps->TxFrameCtrReg != 0 ) { nl = true; printf( "\nTx frame ctr     : %u",   toUint( ps->TxFrameCtr ) );  fflush( stdout ); }
	if( ps->TxBlockCtrReg != 0 ) { nl = true; printf( "\nTx block ctr     : %u",   toUint( ps->TxBlockCtr ) );  fflush( stdout ); }
    if( ps->RxTVCtrReg    != 0 ) { nl = true; printf( "\nRx tvalid us ctr : %u",   toUint( ps->RxTVCtr    ) );  fflush( stdout ); }
    if( ps->TxTVCtrReg    != 0 ) { nl = true; printf( "\nTx tvalid us ctr : %u",   toUint( ps->TxTVCtr    ) );  fflush( stdout ); }
	if( nl ) { printf( "\n" ); fflush( stdout ); }
}

void Update( SysData *ps )
{
	SysDvce *pd = ps->pSysDvce;
	if( ps->RxDataCtrReg  ) ps->RxDataCtr  = pd->ReadSlvReg( pd, ps->RxDataCtrReg  );
	if( ps->RxSymbCtrReg  ) ps->RxSymbCtr  = pd->ReadSlvReg( pd, ps->RxSymbCtrReg  );
	if( ps->RxFrameCtrReg ) ps->RxFrameCtr = pd->ReadSlvReg( pd, ps->RxFrameCtrReg );
	if( ps->RxBlockCtrReg ) ps->RxBlockCtr = pd->ReadSlvReg( pd, ps->RxBlockCtrReg );
	if( ps->TxDataCtrReg  ) ps->TxDataCtr  = pd->ReadSlvReg( pd, ps->TxDataCtrReg  );
	if( ps->TxSymbCtrReg  ) ps->TxSymbCtr  = pd->ReadSlvReg( pd, ps->TxSymbCtrReg  );
	if( ps->TxFrameCtrReg ) ps->TxFrameCtr = pd->ReadSlvReg( pd, ps->TxFrameCtrReg );
	if( ps->TxBlockCtrReg ) ps->TxBlockCtr = pd->ReadSlvReg( pd, ps->TxBlockCtrReg );
    if( ps->RxTVCtrReg    ) ps->RxTVCtr    = pd->ReadSlvReg( pd, ps->RxTVCtrReg    );
    if( ps->TxTVCtrReg    ) ps->TxTVCtr    = pd->ReadSlvReg( pd, ps->TxTVCtrReg    );
}


void SysDataInit( SysDvce *pd,
		u32 RxDataCtrReg,
		u32 RxSymbCtrReg,
		u32 RxFrameCtrReg,
		u32 RxBlockCtrReg,
		u32 TxDataCtrReg,
		u32 TxSymbCtrReg,
		u32 TxFrameCtrReg,
		u32 TxBlockCtrReg,
		u32 RxTVCtrReg,
		u32 TxTVCtrReg
		)
{
	SysData *ps       =  pd->pSysData    ;
	ps->pSysDvce      =  pd              ;
	ps->RxDataCtrReg  =  RxDataCtrReg    ;
	ps->RxSymbCtrReg  =  RxSymbCtrReg    ;
	ps->RxFrameCtrReg =  RxFrameCtrReg   ;
	ps->RxBlockCtrReg =  RxBlockCtrReg   ;
	ps->TxDataCtrReg  =  TxDataCtrReg    ;
	ps->TxSymbCtrReg  =  TxSymbCtrReg    ;
	ps->TxFrameCtrReg =  TxFrameCtrReg   ;
	ps->TxBlockCtrReg =  TxBlockCtrReg   ;
    ps->RxTVCtrReg    =  RxTVCtrReg      ;
    ps->TxTVCtrReg    =  TxTVCtrReg      ;
	ps->Clear         =  Clear           ;
	ps->Update        =  Update          ;
	ps->Dump          =  Dump            ;
	ps->Clear( ps )                      ;
}


