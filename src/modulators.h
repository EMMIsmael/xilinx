/*
 * modulators.h
 *
 *  Created on: 29 Nov 2017
 *      Author: JULIAN MORTIMER
 */
#ifndef MODULATORS_H_
#define MODULATORS_H_

#define DMOD_NAME_LEN  ( 10 )

typedef struct sModem       Modem;
typedef struct sDmodTFDesc  DmodTFDesc;
typedef struct sDmodLUT     DmodLUT;

typedef enum eModType {
    MOD_QPSK  = 0,
    MOD_4QAM  = 1,
    MOD_16QAM = 2,
    MOD_64QAM = 3
} ModType;

typedef enum e_softbits {
    FOUR_SOFT_BITS  = 4
} SoftBits;

typedef enum eDmodMode   {
    DMOD_MODE_DATA = 0,
    DMOD_MODE_PHASE
} DmodMode ;

struct sModem
{
	char             Name[ DMOD_NAME_LEN + 1 ];
    ModType          ModType;
    DmodMode         DmodMode;
	u32              SymbBits;
	u32              SoftBits;
	u32              AddrBits;
	u32              DataBits;
	u32              IdxMask;
	u32              *pConstellation;
	DmodTFDesc       *pTFDesc;
	u32              *pLUTf;
	u32              *pLUTg;
	u32              *pLUTh;
	HwModulator      *pHwModulator;
    HwDemodulator    *pHwDemodulator;
	Cplx32 ( * Modulate       ) ( Modem *pd, u32    data, double sf );
	u32    ( * Demodulate     ) ( Modem *pd, Cplx32 data, double sf );
    void   ( * CalcDmodTF     ) ( Modem *pd, u32 bit );
    void   ( * CalcDmodLUTBit ) ( Modem *pd, u32 bit );
    void   ( * MakeDmodLUT    ) ( Modem *pd );
    void   ( * InitModem ) (
                 Modem          *pd,
                 ModType        ModType,
                 u32            *pLUTf,
                 u32            *pLUTg,
                 u32            *pLUTh,
                 HwModulator    *pHwModulator,
                 HwDemodulator  *pHwDemodulator,
                 void           (* CalcDmodTF     ) ( Modem *pd, u32 bit ),
                 void           (* CalcDmodLUTBit ) ( Modem *pd, u32 bit ),
                 void           (* MakeDmodLUT    ) ( Modem *pd )
                 );
};




struct sDmodTFDesc {
    u32 InitVal;
    const u32 TFLength;
    const u32 TFTable[ 4 ];
};


#endif /* MODULATORS_H_ */
