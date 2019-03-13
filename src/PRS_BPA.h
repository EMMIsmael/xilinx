/*
 * PRS_BPA.h
 *
 *  Created on: 12 Sep 2018
 *      Author: JULIAN MORTIMER
 */

#ifndef PRS_BPA_H_
#define PRS_BPA_H_

#define BPAPAR( Param ) DEVICE_LOC_PAR( pd, BPA, Param )
#define BPADAT( Param ) DEVICE_LOC_DAT( pd, BPA, Param )
#define BPA_N_SLV_REGS                 ( 64 )

#define BPA_CTRL_REG                   (  0 )
#define BPA_N_FFT_REG                  (  2 )
#define BPA_N_SYMB_FRAME_REG           (  3 )
#define BPA_N_FRAME_BLOCK_REG          (  4 )

#define BPA_INFO_REG                   (  1 )
#define BPA_STAT_REG                   ( 32 )

#define BPA_TX_DATA_CTR_REG            ( 48 )
#define BPA_TX_SYMB_CTR_REG            ( 49 )
#define BPA_TX_FRAME_CTR_REG           ( 50 )
#define BPA_TX_BLOCK_CTR_REG           ( 51 )

// Control register bits
#define BPA_RESET                      ( 0x00000001 )
#define BPA_INT_ACK                    ( 0x00000002 )

#define BPA_HW_INFO_REG                ( 0xBEBEAA20 )

#define DELTA_PHI_MAX                  ( M_PI )





typedef struct sprsphi PrsPhi;
typedef enum eprsphi {
    PHI_ZERO,
    PHI_PI_BY_2,
    PHI_PI,
    PHI_3_PI_BY_2,
    PHI_UNDEF
} PrsPhase;


struct sprsphi {
    _Bool init;
    u32  phase [ SYS_MAX_CARRIERS ];
    PrsPhase ( *Phase )      ( PrsBPA * pd, u32 n    );
    void     ( *Initialise ) ( PrsBPA * pd           );
};




typedef struct sLinReg LinReg;
struct sLinReg {
    double Alpha;
    double Beta;
};




#endif /* PRS_BPA_H_ */
