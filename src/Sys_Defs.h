/*
 * Sys_Defs.h
 *
 *  Created on: 6 Jul 2017
 *      Author: JULIAN MORTIMER
 */

#ifndef SYS_DEFS_H_
#define SYS_DEFS_H_
#define DEVICE_INIT( pSysCtrl, pSysDesc, pSysPtrs, Dvce, Dvcn, Dscn, DVCN ) do { Dvce ## Init( pSysCtrl->p ## Dvcn, Dvcn ## Name, DVCN ## _BASEADDR, pSysCtrl, &Dvcn ## SysDvce, pSysDesc, pSysPtrs, &Dvcn ## SysData, &Dscn ## Desc, &Dvcn ## Data ); } while( 0 )

#define SYS_CTL( Dptr, X ) Dptr->pSysDvce->pSysCtrl->X
#define SYS_DVC( Dptr, X ) Dptr->pSysDvce->X
#define SYS_PAR( Dptr, X ) Dptr->pSysDvce->pSysDesc->X
#define SYS_DAT( Dptr, X ) Dptr->pSysDvce->pSysData->X
#define SYS_PTR( Dptr, X ) Dptr->pSysDvce->pSysPtrs->X
#define DVC_PAR( Dptr, X ) Dptr->pDvcDesc->X
#define DVC_DAT( Dptr, X ) Dptr->pDvcData->X

#define SYSCTL( X ) SYS_CTL( pd, X )
#define SYSDVC( X ) SYS_DVC( pd, X )
#define SYSPAR( X ) SYS_PAR( pd, X )
#define SYSDAT( X ) SYS_DAT( pd, X )
#define SYSPTR( X ) SYS_PTR( pd, X )
#define DVCPAR( X ) DVC_PAR( pd, X )
#define DVCDAT( X ) DVC_DAT( pd, X )




#define DVC_RESET( Dptr )                    do{ SYS_DVC( Dptr, Reset         ( Dptr->pSysDvce           ) ); } while( 0 )
#define DVC_RD_STAT_REG( Dptr )              (   SYS_DVC( Dptr, ReadStatReg   ( Dptr->pSysDvce           ) ) )
#define DVC_RD_SLV_REG( Dptr, n )            (   SYS_DVC( Dptr, ReadSlvReg    ( Dptr->pSysDvce, n        ) ) )
#define DVC_RD_SLV_REG_64( Dptr, nlo, nhi )  (   SYS_DVC( Dptr, ReadSlvReg64  ( Dptr->pSysDvce, nlo, nhi ) ) )
#define DVC_WR_SLV_REG( Dptr, n, u )         do{ SYS_DVC( Dptr, WriteSlvReg   ( Dptr->pSysDvce, n, u     ) ); } while( 0 )
#define DVC_N_SLV_REG( Dptr )                (   SYS_DVC( Dptr, NSlvReg       ( Dptr->pSysDvce           ) ) )
#define DVC_DEVICE_ID( Dptr )                (   SYS_DVC( Dptr, DeviceId      ( Dptr->pSysDvce           ) ) )
#define DVC_INT_ACK( Dptr )                  do{ SYS_DVC( Dptr, DeviceId      ( Dptr->pSysDvce          ) ); } while( 0 )
#define DVC_DUMP( Dptr )                     do{ Dptr->Dump( Dptr )                                        ; } while( 0 )
#define DVC_DUMP_SLV_REG( Dptr, bBin )       do{ SYS_DVC( Dptr, DumpSlvReg    ( Dptr->pSysDvce, bBin    ) ); } while( 0 )
#define DVC_FUNC( Dptr, Func, ... )          ( Dptr->Func( Dptr, ##__VA_ARGS__ ))
#define DVC_PROC( Dptr, Func, ... )          do{ Dptr->Func( Dptr, ##__VA_ARGS__ )                         ; } while( 0 )

#define SYSRESET( Dptr )                     do { Dptr->pSysDvce->Reset( Dptr->pSysDvce ); } while( 0 )
#define SYSDVCEDUMP( Dptr )                  do { Dptr->pSysDvce->Dump(  Dptr->pSysDvce ); } while( 0 )
#define SYSUPDDATA( Dptr )                   do { Dptr->pSysDvce->pSysData->Update( Dptr->pSysDevc->pSysData ); } while( 0 )
#define SYSCLRDATA( Dptr )                   do { Dptr->pSysDvce->pSysData->Clear( Dptr->pSysDevc->pSysData );  } while( 0 )
#define SYSDATADUMP( Dptr )                  do { Dptr->pSysDvce->pSysData->Dump( Dptr->pSysDvce->pSysData );   } while( 0 )
#define SYSRDSLVREG( Dptr, n )               ( Dptr->pSysDvce->ReadSlvReg( Dptr->pSysDvce, n ) )
#define SYSRDSLVREG64( Dptr, nlo )           ( Dptr->pSysDvce->ReadSlvReg64( Dptr->pSysDvce, nlo ) )
#define SYSWRSLVREG( Dptr, n, u )            do { Dptr->pSysDvce->WriteSlvReg( Dptr->pSysDvce, n, u );          } while( 0 )
#define SYSWRSLVREG64( Dptr, nlo, u )        do { Dptr->pSysDvce->WriteSlvReg64( Dptr->pSysDvce, nlo, u );      } while( 0 )
#define SYSRDSTATREG( Dptr )                 ( Dptr->pSysDvce->ReadStatReg( Dptr->pSysDvce ) )

// System Control Access
#define SYSCTRLPARAM( Param )                ( pSysCtrl->Param )
#define SYSCTRLFUNC( Func, ... )             ( pSysCtrl->Func( pSysCtrl, ##__VA_ARGS__ ))
#define SYSCTRLPROC( Func, ... )             do { pSysCtrl->Func( pSysCtrl, ##__VA_ARGS__ ); } while( 0 )


#define SYS_MAX_OBJ_NAME_LEN     ( 100 )
#define SYS_MAX_DEVICES          ( 64 )
#define SYS_MAX_ERR_MSG_LEN      ( 50 )

#define SYS_MAX_CARRIERS         ( 2048 )
#define SYS_MIN_CARRIERS         ( 256 )
#define SYS_N_CARRIER_CFG        ( 4 )
#define SYS_MAX_BITS_PER_CARRIER ( 6 )
#define SYS_N_CPRE_MIN           ( 0 )
#define SYS_N_DSYMB_MIN          ( 1 )
#define SYS_N_DSYMB_MAX          ( 1024 )
#define SYS_N_TSYMB_MIN          ( 0 )
#define SYS_N_TSYMB_MAX          ( 4 )
#define SYS_N_CHAN_PRSBPA        ( 9 )
#define SYS_N_CHAN_PRSCIR        ( 21 )
#define SYS_N_SLV_REG_MAX        ( 128 )
#define SYS_N_SOFT_BITS          ( FOUR_SOFT_BITS )

#define PL_DEVICE_N_SLV_REG_MASK ( 0xFF )
#define PL_DEVICE_ID_LSB         ( 8 )

#define SYS_RF_ADC_N_BITS        ( 12 )
#define SYS_RF_DAC_N_BITS        ( 12 )
#define SYS_RF_ADC_MAX           (   ( 1 << ( SYS_RF_ADC_N_BITS - 1 ) ) - 1 )
#define SYS_RF_ADC_MIN           ( - ( 1 << ( SYS_RF_ADC_N_BITS - 1 ) ) )
#define SYS_RF_DAC_MAX           (   ( 1 << ( SYS_RF_ADC_N_BITS - 1 ) ) - 1 )
#define SYS_RF_DAC_MIN           ( - ( 1 << ( SYS_RF_ADC_N_BITS - 1 ) ) )

#define SYS_N_FFT_CTRLS           ( 6 )

#define TX_SOURCE_AXIBRAMC_BASEADDR          ( 0x43C00000 )
#define TX_AXISOURCE_BASEADDR                ( 0x43C10000 )
#define TX_FFTC_BASEADDR                     ( 0x43C60000 )
#define RX_NSD_FFTC_BASEADDR                 ( 0x43CD0000 )
#define TX_TSINJ_BASEADDR                    ( 0x43C70000 )
#define RFCHAN_TREMUL_BASEADDR               ( 0x43C80000 )
#define RFCHAN_TREMULSW_BASEADDR             ( 0x43CB0000 )
#define RFCHAN_TREMUL_FREQCTR_BASEADDR       ( 0x43CA0000 )
#define AD9364_CLK_FREQCTR_BASEADDR          ( 0x43C30000 )
#define AD9364_DAC_FREQCTR_BASEADDR          ( 0x43C40000 )
#define AD9364_ADC_FREQCTR_BASEADDR          ( 0x43D10000 )
#define RX_NSDET_BASEADDR                    ( 0x43CC0000 )
#define RX_TDD_AXISINK_BASEADDR              ( 0x43CE0000 )
#define RX_PRS_AXISINK_BASEADDR              ( 0x43D90000 )
#define RX_DTA_AXISINK_BASEADDR              ( 0x43CF0000 )
#define RX_NSD_B00_AXISINK_BASEADDR          ( 0x43D40000 )
#define RX_NSD_B01_AXISINK_BASEADDR          ( 0x43D70000 )
#define RX_NSD_B15_AXISINK_BASEADDR          ( 0x43D80000 )
#define RX_NSD_FFT_AXISINK_BASEADDR          ( 0x43D30000 )
#define RX_NSD_CORDIC_AXISINK_BASEADDR       ( 0x43D20000 )
#define RX_NSD_A00_AXISINK_BASEADDR          ( 0x43D00000 )
#define RX_NSD_A01_AXISINK_BASEADDR          ( 0x43D50000 )
#define RX_NSD_A15_AXISINK_BASEADDR          ( 0x43D60000 )

// Simulated devices
#define RX_PRS_BPA_BASEADDR                  ( 0x00000001 )
#define RX_PRS_CIR_BASEADDR                  ( 0x00000002 )
#define RX_PRS_FFTC_BASEADDR                 ( 0x00000003 )
#define RX_CIR_FFTC_BASEADDR                 ( 0x00000004 )
#define RX_BPA_FFTC_BASEADDR                 ( 0x00000005 )
#define RX_FFTC_BASEADDR                     ( 0x00000006 )

#define SYS_TIMER_BASEADDR                   ( 0x43C20000 )
#define RFCHAN_TREMUL_CLKWIZ_BASEADDR        ( 0x43C90000 )


#define AXI_BRAM_CTRLR_BASE_ADDR             ( 0x40000000 )
#define AXI_BRAM_CTRLR_BRAM_SIZE             ( 8192 )

#define DUMMY_DVC_BASE_ADDR                  ( ( INTPTR ) Dummy_Device )
#define DATA_DEMOD_LUT_SIZE                  ( 4096 )

#define DATA_MOD_SCALE                       ( 0.125 )

#define N_LINES_DUMP                         ( 32 )

typedef struct sSysCtrl SysCtrl;



typedef enum __dec {
    DEVICE_NAME_TOO_LONG = 0,
    DEVICE_BAD_ID,
    DEVICE_NOT_FOUND,
    DEVICE_NO_MORE_ALLOWED,
    DEVICE_OUT_OF_MEMORY,
    DEVICE_BAD_DEVICE_ID,
    DEVICE_SUCCESS
} DeviceErrCode;

#endif /* SYS_DEFS_H_ */
