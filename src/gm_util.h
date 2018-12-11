/*
 * gm_util.h
 *
 *  Created on: 15 Mar 2016
 *      Author: RTL3
 */

#ifndef GM_UTIL_H_
#define GM_UTIL_H_

#define CHAR_RX_BUF_SIZE ( 10 )


#define MAX_FFT_SIZE ( 32768 )




// PRNG
typedef struct sprng32 PRNG32;
struct sprng32 {
	u32 seed;
	u32 state;
	u32 (* prng )( PRNG32 *p, u32 load );
};

typedef struct sprng64 PRNG64;
struct sprng64 {
    u64 seed;
    u64 state;
    u64 (* prng )( PRNG64 *p, u32 load );
};


// Dump to Console
#define RD_DVC_DATA( ptr, srce, addr ) ( ptr->Read ## srce( ptr, addr ))
typedef enum
{
    DUMP_SEL_PS7_RAM,
	DUMP_SEL_SOURCE,
	DUMP_SEL_SINK,
    DUMP_SEL_TI_BUF,
    DUMP_SEL_TI_MAP,
    DUMP_SEL_CE_BUF,
    DUMP_SEL_CE_MAP,
    DUMP_SEL_FI_TX_BUF,
    DUMP_SEL_FI_INTLVR,
    DUMP_SEL_FI_MOD,
    DUMP_SEL_FFT_CTRLR,
    DUMP_SEL_FREQ_INTLVR,
    DUMP_SEL_TSINJ_PSBUF,
    DUMP_SEL_TSINJ_TXBUF,
    DUMP_SEL_TSINJ_TSBUF,
    DUMP_SEL_DDM_UTIL_DMOD_LUT,
    DUMP_SEL_DDM_UTIL_RAW_LUT,
    DUMP_SEL_TRE_PS_LUT,
    DUMP_SEL_TRE_HW_LUT,
    DUMP_SEL_PRS_BRAM,
    DUMP_SEL_CIR_BRAM_FFT,
    DUMP_SEL_CIR_BRAM_ABS,
    DUMP_SEL_BPA_BRAM,
} dmpSelector;

typedef enum
{
	DUMP_FMT_HEX          = 0,
	DUMP_FMT_DEC          = 1,
	DUMP_FMT_IQ_SIGNED    = 2,
	DUMP_FMT_SOFT_BITS    = 3,
	DUMP_FMT_SOFT_BITS_HW = 4,
	DUMP_FMT_HEX_DEC      = 5,
	DUMP_FMT_BIN32        = 6,
	DUMP_FMT_DEC_SIGNED   = 7,
	DUMP_FMT_GRAPH        = 8,
	DUMP_FMT_HEX_DEMOD    = 9,
	DUMP_FMT_CPLX         = 10,
} dmpFormat;

typedef enum
{
	DUMP_TYP_MEMORY       = 0,
  	DUMP_TYP_SLV_REG      = 1,
  	DUMP_TYP_SLV_REG_LOC  = 2,
} dmpType;


typedef u32 dmpAddr;
typedef u32 dmpArg;
typedef char dmpTitle[];

typedef struct ctrlch CtrlChar;
struct ctrlch {
    char eseq[ESEQ_RX_BUF_SIZE];
    char code;
    char disp[20];
};




#endif /* GM_UTIL_H_ */

