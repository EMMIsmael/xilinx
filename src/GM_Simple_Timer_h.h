/*
 * GM_Simple_Timer_h.h
 *
 *  Created on: 17 May 2017
 *      Author: JULIAN MORTIMER
 */

#ifndef GM_SIMPLE_TIMER_H_H_
#define GM_SIMPLE_TIMER_H_H_


#define SYS_TIMER_DEVICE_ID     ( 0x16E02506 )
#define SYS_TIMER_N_SLV_REGS    ( 6 )


typedef struct sTimerStatus TimerStatus;
struct sTimerStatus {
    u32   TimerIsRunning;
    u32   ErrorFlags;
    float ClockFrequency;
};



typedef struct sTimerAPI TimerAPI;
struct sTimerAPI {
	void          *pDeviceManager;
    DeviceErrCode (* DevManErr )( TimerAPI *pa );
    const char*   (* DevManErrMsg )( TimerAPI *pa );
    u32           (* DevManInitialised )( TimerAPI *pa );
    DeviceErrCode (* InitTimerAPI )( TimerAPI *pa );
    u32           (* InitTimer )( TimerAPI *pa, u32 baseAddr, u32 id, const char name[], u32 deviceId);
    const char*   (* Name )( TimerAPI *pa, u32 Id );
    void          (* ResetTimer )( TimerAPI *pa, u32 id );
    void          (* StartTimer )( TimerAPI *pa, u32 id );
    void          (* StopTimer )( TimerAPI *pa, u32 id );
    void          (* ClearTimer )( TimerAPI *pa, u32 id );
    u32           (* ReadTimerLo )( TimerAPI *pa, u32 id );
    u32           (* ReadTimerHi )( TimerAPI *pa, u32 id );
    u32           (* ReadClearTimer )( TimerAPI *pa, u32 id );
    double        (* ReadTimerLong )( TimerAPI *pa, u32 id );
    double        (* ReadTimerLongSec )( TimerAPI *pa, u32 id );
    u32           (* TimerIsRunning )( TimerAPI *pa, u32 id );
    void          (* ReadTimerSlvRegs )( TimerAPI *pa, u32 id, u32 *pData );
    void          (* WriteTimerSlvReg )( TimerAPI *pa, u32 id, u32 reg, u32 data );
    u32           (* ReadTimerSlvReg )( TimerAPI *pa, u32 id, u32 reg );
    u32           (* NSlvReg )( TimerAPI *pa, u32 id  );
    u32           (* DeviceId )( TimerAPI *pa, u32 id );
    void          (* DumpTimerSlvRegs )( TimerAPI *pa, u32 id, u32 fmtsel );
};






#endif /* GM_SIMPLE_TIMER_H_H_ */
