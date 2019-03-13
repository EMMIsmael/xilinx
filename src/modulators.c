#include "include.h"

Modem      *GetModulator( ModType m );
const char *GetModTypeStr( ModType m );

static Cplx32 Modulate   ( Modem *pm, u32    data, double sf );
static u32    Demodulate ( Modem *pm, Cplx32 data, double sf );
void  InitModem (
        Modem *pd,
        ModType ModType,
        u32   *pLUTf,
        u32   *pLUTg,
        u32   *pLUTh,
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

static u32 const4QAM[ 4 ] = { //[imag,real] ToDo: Table Only Valid for 1024-FFT
		0x09000900,
		0xF7000900,
		0x0900F700,
		0xF700F700
};


static u32 const16QAM[ 16 ] = {//[imag,real] ToDo: Table Only Valid for 1024-FFT
        0x03000300,
        0x09000300,
        0xFD000300,
        0xF7000300,
        0x03000900,
        0x09000900,
        0xFD000900,
        0xF7000900,
        0x0300FD00,
        0x0900FD00,
        0xFD00FD00,
        0xF700FD00,
        0x0300F700,
        0x0900F700,
        0xFD00F700,
        0xF700F700,
};


static u32 const64QAM[ 64 ] = {//[imag,real] ToDo: Table Only Valid for 1024-FFT
        0x01400140,
        0x03C00140,
        0x08C00140,
        0x06400140,
        0xFEC00140,
        0xFC400140,
        0xF7400140,
        0xF9C00140,
        0x014003C0,
        0x03C003C0,
        0x08C003C0,
        0x064003C0,
        0xFEC003C0,
        0xFC4003C0,
        0xF74003C0,
        0xF9C003C0,
        0x014008C0,
        0x03C008C0,
        0x08C008C0,
        0x064008C0,
        0xFEC008C0,
        0xFC4008C0,
        0xF74008C0,
        0xF9C008C0,
        0x01400640,
        0x03C00640,
        0x08C00640,
        0x06400640,
        0xFEC00640,
        0xFC400640,
        0xF7400640,
        0xF9C00640,
        0x0140FEC0,
        0x03C0FEC0,
        0x08C0FEC0,
        0x0640FEC0,
        0xFEC0FEC0,
        0xFC40FEC0,
        0xF740FEC0,
        0xF9C0FEC0,
        0x0140FC40,
        0x03C0FC40,
        0x08C0FC40,
        0x0640FC40,
        0xFEC0FC40,
        0xFC40FC40,
        0xF740FC40,
        0xF9C0FC40,
        0x0140F740,
        0x03C0F740,
        0x08C0F740,
        0x0640F740,
        0xFEC0F740,
        0xFC40F740,
        0xF740F740,
        0xF9C0F740,
        0x0140F9C0,
        0x03C0F9C0,
        0x08C0F9C0,
        0x0640F9C0,
        0xFEC0F9C0,
        0xFC40F9C0,
        0xF740F9C0,
        0xF9C0F9C0
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


static Modem modQPSK   = { "QPSK  ", MOD_QPSK , DMOD_MODE_DATA, 2, 4, 12, 12, 0x03, constQPSK , tfDesc4QAM,  NULL, NULL, NULL, NULL, NULL, Modulate, Demodulate, NULL, NULL, NULL, InitModem };
static Modem mod4QAM   = { "4-QAM ", MOD_4QAM , DMOD_MODE_DATA, 2, 4, 12, 12, 0x03, const4QAM , tfDesc4QAM,  NULL, NULL, NULL, NULL, NULL, Modulate, Demodulate, NULL, NULL, NULL, InitModem };
static Modem mod16QAM  = { "16-QAM", MOD_16QAM, DMOD_MODE_DATA, 4, 4, 12, 12, 0x0f, const16QAM, tfDesc16QAM, NULL, NULL, NULL, NULL, NULL, Modulate, Demodulate, NULL, NULL, NULL, InitModem };
static Modem mod64QAM  = { "64-QAM", MOD_64QAM, DMOD_MODE_DATA, 6, 4, 12, 12, 0x3f, const64QAM, tfDesc64QAM, NULL, NULL, NULL, NULL, NULL, Modulate, Demodulate, NULL, NULL, NULL, InitModem };

Modem UtilModem        = { ""      , MOD_QPSK ,              0, 0, 0,  0,  0,    0,       NULL,        NULL, NULL, NULL, NULL, NULL, NULL,     NULL,       NULL, NULL, NULL, NULL, InitModem };
Modem TxModem          = { ""      , MOD_QPSK ,              0, 0, 0,  0,  0,    0,       NULL,        NULL, NULL, NULL, NULL, NULL, NULL,     NULL,       NULL, NULL, NULL, NULL, InitModem };
Modem RxModem          = { ""      , MOD_QPSK ,              0, 0, 0,  0,  0,    0,       NULL,        NULL, NULL, NULL, NULL, NULL, NULL,     NULL,       NULL, NULL, NULL, NULL, InitModem };

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
	//u32 idx = data & pm->IdxMask;
	//return( ScaleCplx32( U32toCplx32( pm->pConstellation[ idx ] ), sf ));
	u32 idx = pm->pConstellation [(data & pm->IdxMask)];
	return( ScaleCplx32( U32toCplx32( idx ), sf ));
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
        u32     *pLUTf,
        u32     *pLUTg,
        u32     *pLUTh,
        HwModulator   *pHwModulator,
        HwDemodulator *pHwDemodulator,
        void    (* CalcDmodTF     ) ( Modem *pd, u32 bit ),
        void    (* CalcDmodLUTBit ) ( Modem *pd, u32 bit ),
        void    (* MakeDmodLUT    ) ( Modem *pd )
        )
{
    Modem *pMod = GetModem( ModType );
    *pd = *pMod;
    pd->pLUTf           = pLUTf;
    pd->pLUTg           = pLUTg;
    pd->pLUTh           = pLUTh;
    pd->pHwModulator    = pHwModulator;
    pd->pHwDemodulator  = pHwDemodulator;
    pd->CalcDmodTF      = CalcDmodTF;
    pd->CalcDmodLUTBit  = CalcDmodLUTBit;
    pd->MakeDmodLUT     = MakeDmodLUT;
}
