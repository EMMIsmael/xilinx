#include "include.h"

Modem      *GetModulator( ModType m );
const char *GetModTypeStr( ModType m );

static Cplx32 Modulate   ( Modem *pm, u32    data, double sf );
static u32    Demodulate ( Modem *pm, Cplx32 data, double sf );
void  InitModem (
        Modem *pd,
        ModType ModType,
        u32   *pLUT,
        HwModulator   *pHwModulator,
        HwDemodulator *pHwDemodulator,
        void  (* CalcDmodTF     ) ( Modem *pd, u32 bit ),
        void  (* CalcDmodLUTBit ) ( Modem *pd, u32 bit ),
        void  (* MakeDmodLUT    ) ( Modem *pd )
);


static u32 constQPSK[4] = {
		0x7FFF0000,
		0x00008001,
		0x00007FFF,
		0x80010000
};

static u32 const4QAM[ 4 ] = {
		0x5A82A57E,
		0xA57EA57E,
		0x5A825A82,
		0xA57E5A82
};


static u32 const16QAM[ 16 ] = {
		0x5A81A57F,
		0x1E2BA57F,
		0xA57FA57F,
		0xE1D5A57F,
		0x5A81E1D5,
		0x1E2BE1D5,
		0xA57FE1D5,
		0xE1D5E1D5,
		0x5A815A81,
		0x1E2B5A81,
		0xA57F5A81,
		0xE1D55A81,
		0x5A811E2B,
		0x1E2B1E2B,
		0xA57F1E2B,
		0xE1D51E2B
};


static u32 const64QAM[ 64 ] = {
		0x5A82A57E,
		0x40A6A57E,
		0x0CEEA57E,
		0x26CAA57E,
		0xA57EA57E,
		0xBF5AA57E,
		0xF312A57E,
		0xD936A57E,
		0x5A82BF5A,
		0x40A6BF5A,
		0x0CEEBF5A,
		0x26CABF5A,
		0xA57EBF5A,
		0xBF5ABF5A,
		0xF312BF5A,
		0xD936BF5A,
		0x5A82F312,
		0x40A6F312,
		0x0CEEF312,
		0x26CAF312,
		0xA57EF312,
		0xBF5AF312,
		0xF312F312,
		0xD936F312,
		0x5A82D936,
		0x40A6D936,
		0x0CEED936,
		0x26CAD936,
		0xA57ED936,
		0xBF5AD936,
		0xF312D936,
		0xD936D936,
		0x5A825A82,
		0x40A65A82,
		0x0CEE5A82,
		0x26CA5A82,
		0xA57E5A82,
		0xBF5A5A82,
		0xF3125A82,
		0xD9365A82,
		0x5A8240A6,
		0x40A640A6,
		0x0CEE40A6,
		0x26CA40A6,
		0xA57E40A6,
		0xBF5A40A6,
		0xF31240A6,
		0xD93640A6,
		0x5A820CEE,
		0x40A60CEE,
		0x0CEE0CEE,
		0x26CA0CEE,
		0xA57E0CEE,
		0xBF5A0CEE,
		0xF3120CEE,
		0xD9360CEE,
		0x5A8226CA,
		0x40A626CA,
		0x0CEE26CA,
		0x26CA26CA,
		0xA57E26CA,
		0xBF5A26CA,
		0xF31226CA,
		0xD93626CA
};




static DmodTFDesc tfDesc4QAM[ 1 ] = {
        {  1, 2896, {  599,     0,     0,     0  }},
};

static DmodTFDesc tfDesc16QAM[ 2 ] = {
        {  1,  964, {  599,  2531,     0,     0  }},
        {  1,  964, { 1565,     0,     0,     0  }}
};

static DmodTFDesc tfDesc64QAM[ 3 ] = {
        {   1, 414, {  599,  1427,  2255,  3082  }},
        {   1, 414, { 1014,  2668,     0,     0  }},
        {   1, 414, { 1840,     0,     0,     0  }}
};


static Modem modQPSK   = { "QPSK  ", MOD_QPSK , DMOD_MODE_DATA, 2, 4, 12, 12, 0x03, constQPSK , tfDesc4QAM,  NULL, NULL, NULL, Modulate, Demodulate, NULL, NULL, NULL, InitModem };
static Modem mod4QAM   = { "4-QAM ", MOD_4QAM , DMOD_MODE_DATA, 2, 4, 12, 12, 0x03, const4QAM , tfDesc4QAM,  NULL, NULL, NULL, Modulate, Demodulate, NULL, NULL, NULL, InitModem };
static Modem mod16QAM  = { "16-QAM", MOD_16QAM, DMOD_MODE_DATA, 4, 4, 12, 12, 0x0f, const16QAM, tfDesc16QAM, NULL, NULL, NULL, Modulate, Demodulate, NULL, NULL, NULL, InitModem };
static Modem mod64QAM  = { "64-QAM", MOD_16QAM, DMOD_MODE_DATA, 6, 4, 12, 12, 0x3f, const64QAM, tfDesc64QAM, NULL, NULL, NULL, Modulate, Demodulate, NULL, NULL, NULL, InitModem };

Modem UtilModem        = { ""      , MOD_QPSK ,              0, 0, 0,  0,  0,    0,       NULL,        NULL, NULL, NULL, NULL,     NULL,       NULL, NULL, NULL, NULL, InitModem };
Modem TxModem          = { ""      , MOD_QPSK ,              0, 0, 0,  0,  0,    0,       NULL,        NULL, NULL, NULL, NULL,     NULL,       NULL, NULL, NULL, NULL, InitModem };
Modem RxModem          = { ""      , MOD_QPSK ,              0, 0, 0,  0,  0,    0,       NULL,        NULL, NULL, NULL, NULL,     NULL,       NULL, NULL, NULL, NULL, InitModem };

static Modem *modems[4] =
{
		&modQPSK,
		&mod4QAM,
		&mod16QAM,
		&mod64QAM
};


// Modulate, dividing by sf, which will be rounded up to a power of two
static Cplx32 Modulate( Modem *pm, u32 data, double sf )
{
	u32 idx = data & pm->IdxMask;
	return( ScaleCplx32( U32toCplx32( pm->pConstellation[ idx ] ), sf ));
}

// Demodulate, must be exact. Redefine according to scheme used
static u32 Demodulate( Modem *p, Cplx32 data, double sf )
{
	for ( int i = 0; i < ( 1 << p->SymbBits ); i++ )
	{
		if ( Cplx32Eq( data, ScaleCplx32( U32toCplx32( p->pConstellation[ i ] ), sf ))) return i;
	}
	return ( 0 );
}

Modem *GetModem( ModType m )
{
	return ( modems[ m ] );
}

const char *GetModemName( ModType m )
{
    return ( modems[ m ]->Name );
}

void  InitModem (
        Modem   *pd,
        ModType ModType,
        u32     *pLUT,
        HwModulator   *pHwModulator,
        HwDemodulator *pHwDemodulator,
        void    (* CalcDmodTF     ) ( Modem *pd, u32 bit ),
        void    (* CalcDmodLUTBit ) ( Modem *pd, u32 bit ),
        void    (* MakeDmodLUT    ) ( Modem *pd )
        )
{
    Modem *pMod = GetModem( ModType );
    *pd = *pMod;
    pd->pLUT            = pLUT;
    pd->pHwModulator    = pHwModulator;
    pd->pHwDemodulator  = pHwDemodulator;
    pd->CalcDmodTF      = CalcDmodTF;
    pd->CalcDmodLUTBit  = CalcDmodLUTBit;
    pd->MakeDmodLUT     = MakeDmodLUT;
}
