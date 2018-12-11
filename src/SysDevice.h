/*
 * SysDevice.h
 *
 *  Created on: 19 Sep 2018
 *      Author: JULIAN MORTIMER
 */

#ifndef SYSDEVICE_H_
#define SYSDEVICE_H_
//===============================
// System Device Common Structure
//===============================
struct sSysDvce {
    char                  Name[ SYS_MAX_OBJ_NAME_LEN + 1 ];
    const SysCtrl         *pSysCtrl;
    const SysDesc         *pSysDesc;
    const SysPtrs         *pSysPtrs;
    SysData               *pSysData;
    u32                   BaseAddr;
    u32                   HwInfoReg;
    u32                   CtrlReg;
    u32                   InfoReg;
    u32                   StatReg;
    u32                   ResetBit;
    u32                   IntAckBit;
    u32                   RdSyncBit;
    u32                   RxDataCtrReg;
    u32                   RxSymbCtrReg;
    u32                   RxFrameCtrReg;
    u32                   RxBlockCtrReg;
    u32                   TxDataCtrReg;
    u32                   TxSymbCtrReg;
    u32                   TxFrameCtrReg;
    u32                   TxBlockCtrReg;
    u32                   RxTVCtrReg;
    u32                   TxTVCtrReg;
    u32                   SimSlvRegs[ SYS_N_SLV_REG_MAX ];
    void (* Reset         ) ( SysDvce *pd );
    u32  (* ReadStatReg   ) ( SysDvce *pd );
    u32  (* ReadSlvReg    ) ( SysDvce *pd, u32 reg );
    u64  (* ReadSlvReg64  ) ( SysDvce *pd, u32 reg );
    void (* WriteSlvReg   ) ( SysDvce *pd, u32 reg, u32 data );
    void (* WriteSlvReg64 ) ( SysDvce *pd, u32 reg, u64 data );
    u32  (* NSlvReg       ) ( SysDvce *pd );
    u32  (* DeviceId      ) ( SysDvce *pd );
    void (* IntAck        ) ( SysDvce *pd );
    void (* Dump          ) ( SysDvce *pd );
    void (* DumpSlvReg    ) ( SysDvce *pd, _Bool bBin );
    u32  (* SimHwReg      ) ( u32 baseAddr, u32 reg, u32 data, _Bool bWrite );
};

typedef struct sSysDvceInfo SysDvceInfo;
struct sSysDvceInfo {
    u32 BaseAddr;
    SysDvce *pSysDvce;
};

#endif /* SYSDEVICE_H_ */
