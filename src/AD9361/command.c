/**************************************************************************//**
 *   @file   command.c
 *   @brief  Implementation of AD9361 Command Driver.
 *   @author DBogdan ( dragos.bogdan@analog.com )
 *******************************************************************************
 * Copyright 2013( c ) Analog Devices, Inc.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *  - Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *  - Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *  - Neither the name of Analog Devices, Inc. nor the names of its
 *    contributors may be used to endorse or promote products derived
 *    from this software without specific prior written permission.
 *  - The use of this software may or may not infringe the patent rights
 *    of one or more patent holders.  This license does not release you
 *    from the requirement that you obtain separate licenses from these
 *    patent holders to use this software.
 *  - Use of the software either in source or binary form, must be run
 *    on or directly connected to an Analog Devices Inc. component.
 *
 * THIS SOFTWARE IS PROVIDED BY ANALOG DEVICES "AS IS" AND ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, NON-INFRINGEMENT,
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL ANALOG DEVICES BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES ( INCLUDING, BUT NOT
 * LIMITED TO, INTELLECTUAL PROPERTY RIGHTS, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION ) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT ( INCLUDING NEGLIGENCE OR OTHERWISE ) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*******************************************************************************/

/******************************************************************************/
/***************************** Include Files **********************************/
/******************************************************************************/
#include "command.h"
#include "console.h"
#include "ad9361_api.h"
#include "dac_core.h"
#include "..\include.h"

shell AD9364Shell;

shellCommandDef AD9364CmdDef[ ] = {
    { "all"  , ""                         ,  0, "",                                                                                            NULL },
    { "all"  , "q"                        ,  0, "              Exit this shell"                               , shellConsoleExit                    },
    { "all"  , ""                         ,  0, "",                                                                                            NULL },
    { "all"  , "register"                 ,  1, "<reg>         Get the specified register value"              , get_register                        },
    { "all"  , "tx_lo_freq"               ,  0, "              Get current TX LO frequency [ MHz ]"           , get_tx_lo_freq                      },
    { "all"  , "tx_lo_freq"               ,  1, "<freq>        Set the TX LO frequency [ MHz ]"               , set_tx_lo_freq                      },
    { "all"  , "tx_samp_freq"             ,  0, "              Get current TX sampling frequency [ Hz ]"      , get_tx_samp_freq                    },
    { "all"  , "tx_samp_freq"             ,  1, "<freq>        Set the TX sampling frequency [ Hz ]"          , set_tx_samp_freq                    },
    { "all"  , "tx_rf_bandwidth"          ,  0, "              Get current TX RF bandwidth [ Hz ]"            , get_tx_rf_bandwidth                 },
    { "all"  , "tx_rf_bandwidth"          ,  1, "<bandwidth>   Set the TX RF bandwidth [ Hz ]"                , set_tx_rf_bandwidth                 },
    { "AD1"  , "tx1_attenuation"          ,  0, "              Get current TX1 attenuation [ mdB ]"           , get_tx1_attenuation                 },
    { "AD1"  , "tx1_attenuation"          ,  1, "<atten>       Set the TX1 attenuation [ mdB ]"               , set_tx1_attenuation                 },
    { "AD2"  , "tx2_attenuation"          ,  0, "              Get current TX2 attenuation [ mdB ]"           , get_tx2_attenuation                 },
    { "AD2"  , "tx2_attenuation"          ,  1, "<atten>       Set the TX2 attenuation [ mdB ]"               , set_tx2_attenuation                 },
    { "all"  , "tx_fir_en"                ,  0, "              Get current TX FIR state"                      , get_tx_fir_en                       },
    { "all"  , "tx_fir_en"                ,  1, "<1=en 0=inh>  Set the TX FIR state"                          , set_tx_fir_en                       },
    { "all"  , "rx_lo_freq"               ,  0, "              Get current RX LO frequency [ MHz ]"           , get_rx_lo_freq                      },
    { "all"  , "rx_lo_freq"               ,  1, "<freq>        Set the RX LO frequency [ MHz ]"               , set_rx_lo_freq                      },
    { "all"  , "rx_samp_freq"             ,  0, "              Get current RX sampling frequency [ Hz ]"      , get_rx_samp_freq                    },
    { "all"  , "rx_samp_freq"             ,  1, "<freq>        Set the RX sampling frequency [ Hz ]"          , set_rx_samp_freq                    },
    { "all"  , "rx_rf_bandwidth"          ,  0, "              Get current RX RF bandwidth [ Hz ]"            , get_rx_rf_bandwidth                 },
    { "all"  , "rx_rf_bandwidth"          ,  1, "<bandwidth>   Set the RX RF bandwidth [ Hz ]"                , set_rx_rf_bandwidth                 },
    { "AD1"  , "rx1_gc_mode"              ,  0, "              Get current RX1 GC mode"                       , get_rx1_gc_mode                     },
    { "AD1"  , "rx1_gc_mode"              ,  1, "<mode>        Set the RX1 GC mode"                           , set_rx1_gc_mode                     },
    { "AD2"  , "rx2_gc_mode"              ,  0, "              Get current RX2 GC mode"                       , get_rx2_gc_mode                     },
    { "AD2"  , "rx2_gc_mode"              ,  1, "<mode>        Set the RX2 GC mode"                           , set_rx2_gc_mode                     },
    { "AD1"  , "rx1_rf_gain"              ,  0, "              Get current RX1 RF gain"                       , get_rx1_rf_gain                     },
    { "AD1"  , "rx1_rf_gain"              ,  1, "<gain>        Set the RX1 RF gain"                           , set_rx1_rf_gain                     },
    { "AD2"  , "rx2_rf_gain"              ,  0, "              Get current RX2 RF gain"                       , get_rx2_rf_gain                     },
    { "AD2"  , "rx2_rf_gain"              ,  1, "<gain>        Set the RX2 RF gain"                           , set_rx2_rf_gain                     },
    { "all"  , "rx_fir_en"                ,  0, "              Get current RX FIR state"                      , get_rx_fir_en                       },
    { "all"  , "rx_fir_en"                ,  1, "<1=en 0=inh>  Set the RX FIR state"                          , set_rx_fir_en                       },
    { "AD1"  , "dds_tx1_tone1_freq"       ,  0, "              Get current DDS TX1 Tone 1 frequency [ Hz ]"   , get_dds_tx1_tone1_freq              },
    { "AD1"  , "dds_tx1_tone1_freq"       ,  1, "<freq>        Set the DDS TX1 Tone 1 frequency [ Hz ]"       , set_dds_tx1_tone1_freq              },
    { "AD1"  , "dds_tx1_tone2_freq"       ,  0, "              Get current DDS TX1 Tone 2 frequency [ Hz ]"   , get_dds_tx1_tone2_freq              },
    { "AD1"  , "dds_tx1_tone2_freq"       ,  1, "<freq>        Set the DDS TX1 Tone 2 frequency [ Hz ]"       , set_dds_tx1_tone2_freq              },
    { "AD1"  , "dds_tx1_tone1_phase"      ,  0, "              Get current DDS TX1 Tone 1 phase [ degrees ]"  , get_dds_tx1_tone1_phase             },
    { "AD1"  , "dds_tx1_tone1_phase"      ,  1, "<phase>       Set the DDS TX1 Tone 1 phase [ degrees ]"      , set_dds_tx1_tone1_phase             },
    { "AD1"  , "dds_tx1_tone2_phase"      ,  0, "              Get current DDS TX1 Tone 2 phase [ degrees ]"  , get_dds_tx1_tone2_phase             },
    { "AD1"  , "dds_tx1_tone2_phase"      ,  1, "<phase>       Set the DDS TX1 Tone 2 phase [ degrees ]"      , set_dds_tx1_tone2_phase             },
    { "AD1"  , "dds_tx1_tone1_scale"      ,  0, "              Get current DDS TX1 Tone 1 scale"              , get_dds_tx1_tone1_scale             },
    { "AD1"  , "dds_tx1_tone1_scale"      ,  1, "<scale>       Set the DDS TX1 Tone 1 scale"                  , set_dds_tx1_tone1_scale             },
    { "AD1"  , "dds_tx1_tone2_scale"      ,  0, "              Get current DDS TX1 Tone 2 scale"              , get_dds_tx1_tone2_scale             },
    { "AD1"  , "dds_tx1_tone2_scale"      ,  1, "<scale>       Set the DDS TX1 Tone 2 scale"                  , set_dds_tx1_tone2_scale             },
    { "AD2"  , "dds_tx2_tone1_freq"       ,  0, "              Get current DDS TX2 Tone 1 frequency [ Hz ]"   , get_dds_tx2_tone1_freq              },
    { "AD2"  , "dds_tx2_tone1_freq"       ,  1, "<freq>        Set the DDS TX2 Tone 1 frequency [ Hz ]"       , set_dds_tx2_tone1_freq              },
    { "AD2"  , "dds_tx2_tone2_freq"       ,  0, "              Get current DDS TX2 Tone 2 frequency [ Hz ]"   , get_dds_tx2_tone2_freq              },
    { "AD2"  , "dds_tx2_tone2_freq"       ,  1, "<freq>        Set the DDS TX2 Tone 2 frequency [ Hz ]"       , set_dds_tx2_tone2_freq              },
    { "AD2"  , "dds_tx2_tone1_phase"      ,  0, "              Get current DDS TX2 Tone 1 phase [ degrees ]"  , get_dds_tx2_tone1_phase             },
    { "AD2"  , "dds_tx2_tone1_phase"      ,  1, "<phase>       Set the DDS TX2 Tone 1 phase [ degrees ]"      , set_dds_tx2_tone1_phase             },
    { "AD2"  , "dds_tx2_tone2_phase"      ,  0, "              Get current DDS TX2 Tone 2 phase [ degrees ]"  , get_dds_tx2_tone2_phase             },
    { "AD2"  , "dds_tx2_tone2_phase"      ,  1, "<phase>       Set the DDS TX2 Tone 2 phase [ degrees ]"      , set_dds_tx2_tone2_phase             },
    { "AD2"  , "dds_tx2_tone1_scale"      ,  0, "              Get current DDS TX2 Tone 1 scale"              , get_dds_tx2_tone1_scale             },
    { "AD2"  , "dds_tx2_tone1_scale"      ,  1, "<scale>       Set the DDS TX2 Tone 1 scale"                  , set_dds_tx2_tone1_scale             },
    { "AD2"  , "dds_tx2_tone2_scale"      ,  0, "              Get current DDS TX2 Tone 2 scale"              , dds_tx2_tone2_scale                 },
    { "AD2"  , "dds_tx2_tone2_scale"      ,  1, "<scale>       Set the DDS TX2 Tone 2 scale"                  , set_dds_tx2_tone2_scale             },
    { "all"  , "ovflow"                   ,  1, "<chan (0|1)>  Overflow Monitor"                              , overflow_monitor                    },
    { "all"  , "bbrfoffs"                 ,  0, "              Run Offset Calibration"                        , calib_dc_offset                     },
    { "all"  , "rx_timing"                ,  0, "              Rx channel timing analysis"                    , rx_timing_analysis                  },
    { "all"  , "timing"                   ,  0, "              Tx and Rx channel timing analysis"             , txrx_timing_analysis                },
    { "\0" }
};


/******************************************************************************/
/************************ Variables Definitions *******************************/
/******************************************************************************/
extern struct dds_state dds_st;
extern struct ad9361_rf_phy *ad9361_phy;


/**************************************************************************//***
 * @brief Overflow Monitor
 *
 * @return None.
*******************************************************************************/
void overflow_monitor( shell *psh )
{
    uint16_t reg_addr = ( uint16_t ) MIN( psh->uarg[ 0 ], 1 ) + REG_CH_1_OVERFLOW;
    struct spi_device spi;
    spi.id_no = 0;

    uint8_t  reg_val   = 0;
    uint8_t  reg_acc   = 0;
    uint32_t ovf [ 7 ] = { 0, 0, 0, 0, 0, 0, 0 };
    for ( int i = 0; i < 10000; i++ ) {
       reg_val  = ad9361_spi_read( &spi, reg_addr );
       reg_acc |= reg_val;
       for ( int j = 0; j < 6; j++ ) {
           ovf [ j ] += ( reg_val & ( 1 << j ) ? 1 : 0 );
       }
    }
    cprintf( "\n" );
    cprintf( "Ch%u Overflow Reg ( 10k reads ) =\n", toUint( reg_addr - REG_CH_1_OVERFLOW + 1 ) );
    cprintf( "| BBPLL Lock | TxInt3 | TxHb3  | TxHb2  | TxQEC  | TxHB1  | TxFIR  | RxFIR  |\n");
    cprintf( "|            | Filt   | Filt   | Filt   | Filt   | Filt   | Filt   | Filt   |\n");
    cprintf( "| %-10s ", u32bit( reg_val, 7 ) ? "Locked" : "Not locked" );
    for ( int i = 6; i >= 1; i-- ) {
        cprintf( "| %-6u ", toUint( ovf [ i ] ) );
    }
    cprintf( "| %-6u |\n", toUint( ovf [ 0 ] ) );
    cprintf( "\n" );
}
/**************************************************************************//***
 * @brief Rx digital port timing analysis
 *
 * @return None.
 *******************************************************************************/
void rx_timing_analysis( shell *psh )
{
    char buf[ 2048 ];
    int32_t len = 0, i = 0;
    len = ad9361_dig_interface_timing_analysis( ad9361_phy, buf, sizeof( buf ) );
    if ( len > 0 ) {
        for ( i = 0; i < len; i++ )
            if ( buf[ i ] == '\0' ) printf( "\\0" ); else printf( "%c", buf[ i ] );
    }
    printf( "\n" );
}

/**************************************************************************//***
 * @brief Tx and Rx digital port timing analysis
 *
 * @return None.
 *******************************************************************************/
void txrx_timing_analysis( shell *psh )
{
    char buf[ 4096 ];
    int32_t len = 0, i = 0;
    len = ad9361_dig_tune_analyse( ad9361_phy, buf, sizeof( buf ) );
    if ( len > 0 )
    {
        for ( i = 0; i < len; i++ )
            if ( buf[ i ] == '\0' ) printf( "\\0" ); else printf( "%c", buf[ i ] );
    }
    printf( "\n" );
}

/**************************************************************************//***
 * @brief Autotune digital port timing
 *
 * @return None.
 *******************************************************************************/
void autotuneDigitalInterface( shell *psh )
{

}

/**************************************************************************//***
 * @brief Displays all available commands.
 *
 * @return None.
 *******************************************************************************/
void get_register( shell *psh ) // "register?" command
{
    uint16_t reg_addr;
    uint8_t reg_val;
    struct spi_device spi;

    spi.id_no = 0;
    reg_addr = ( uint16_t ) psh->uarg[ 0 ];
    reg_val = ad9361_spi_read( &spi, reg_addr );
    //      cprintf( "register[ 0x%x ]=0x%x\n", reg_addr, reg_val );
    cprintf( "register[ 0x%04x ] = 0x%02x\n", toUint( reg_addr ), toUint( reg_val ) );
}

/**************************************************************************//***
 * @brief Gets current TX LO frequency [ MHz ].
 *
 * @return None.
 *******************************************************************************/
void get_tx_lo_freq( shell *psh ) // "tx_lo_freq?" command
{
    uint64_t lo_freq_hz;

    ad9361_get_tx_lo_freq( ad9361_phy, &lo_freq_hz );
    lo_freq_hz /= 1000000;
    cprintf( "tx_lo_freq = %lu\n", toULint( lo_freq_hz ) );
}

/**************************************************************************//***
 * @brief Sets the TX LO frequency [ MHz ].
 *
 * @return None.
 *******************************************************************************/
void set_tx_lo_freq( shell *psh ) // "tx_lo_freq=" command
{
    uint64_t lo_freq_hz;

    lo_freq_hz = ( uint64_t ) psh->darg[ 0 ];
    lo_freq_hz *= 1000000;
    ad9361_set_tx_lo_freq( ad9361_phy, lo_freq_hz );
    lo_freq_hz /= 1000000;
    cprintf( "tx_lo_freq = %lu\n", toULint( lo_freq_hz ) );
}

/**************************************************************************//***
 * @brief Gets current sampling frequency [ Hz ].
 *
 * @return None.
 *******************************************************************************/
void get_tx_samp_freq( shell *psh ) // "tx_samp_freq?" command
{
    uint32_t sampling_freq_hz;

    ad9361_get_tx_sampling_freq( ad9361_phy, &sampling_freq_hz );
    cprintf( "tx_samp_freq = %u\n", toUint( sampling_freq_hz ) );
}

/**************************************************************************//***
 * @brief Sets the sampling frequency [ Hz ].
 *
 * @return None.
 *******************************************************************************/
void setTxSampFreq( uint32_t fs ) // "tx_samp_freq=" command
{
    ad9361_set_tx_sampling_freq( ad9361_phy, fs );
    ad9361_get_tx_sampling_freq( ad9361_phy, &fs );
    dds_update( ad9361_phy );
    cprintf( "tx_samp_freq = %u\n", toUint( fs ) );
}

void set_tx_samp_freq( shell *psh ) // "tx_samp_freq=" command
{
    setTxSampFreq( psh->uarg[ 0 ] );
}

/**************************************************************************//***
 * @brief Gets current TX RF bandwidth [ Hz ].
 *
 * @return None.
 *******************************************************************************/
void get_tx_rf_bandwidth( shell *psh ) // "tx_rf_bandwidth?" command
{
    uint32_t bandwidth_hz;

    ad9361_get_tx_rf_bandwidth( ad9361_phy, &bandwidth_hz );
    cprintf( "tx_rf_bandwidth = %u\n", toUint( bandwidth_hz ) );
}

/**************************************************************************//***
 * @brief Sets the TX RF bandwidth [ Hz ].
 *
 * @return None.
 *******************************************************************************/
void setTxRfBandwidth( uint32_t bw )
{
    ad9361_set_tx_rf_bandwidth( ad9361_phy, bw );
    ad9361_get_tx_rf_bandwidth( ad9361_phy, &bw );
    cprintf( "tx_rf_bandwidth = %u\n", toUint( bw ) );
}

void set_tx_rf_bandwidth( shell *psh ) // "tx_rf_bandwidth=" command
{
    setTxRfBandwidth( psh->uarg[ 0 ] );
}


/**************************************************************************//***
 * @brief Gets current TX1 attenuation [ mdB ].
 *
 * @return None.
 *******************************************************************************/
void get_tx1_attenuation( shell *psh ) // "tx1_attenuation?" command
{
    uint32_t attenuation_mdb;

    ad9361_get_tx_attenuation( ad9361_phy, 0, &attenuation_mdb );
    cprintf( "tx1_attenuation = %u.%u dB\n", toUint( attenuation_mdb / 1000 ), toUint( attenuation_mdb % 1000 )  );
}

/**************************************************************************//***
 * @brief Sets the TX1 attenuation [ mdB ].
 *
 * @return None.
 *******************************************************************************/
void set_tx1_attenuation( shell *psh ) // "tx1_attenuation=" command
{
    uint32_t attenuation_mdb;

    attenuation_mdb = psh->uarg[ 0 ];
    ad9361_set_tx_attenuation( ad9361_phy, 0, attenuation_mdb );
    get_tx1_attenuation( NULL );
}

/**************************************************************************//***
 * @brief Gets current TX2 attenuation [ mdB ].
 *
 * @return None.
 *******************************************************************************/
void get_tx2_attenuation( shell *psh ) // "tx1_attenuation?" command
{
    uint32_t attenuation_mdb;

    ad9361_get_tx_attenuation( ad9361_phy, 1, &attenuation_mdb );
    cprintf( "tx2_attenuation = %u.%u dB\n", toUint( attenuation_mdb / 1000 ), toUint( attenuation_mdb % 1000 )  );
}

/**************************************************************************//***
 * @brief Sets the TX2 attenuation [ mdB ].
 *
 * @return None.
 *******************************************************************************/
void set_tx2_attenuation( shell *psh ) // "tx1_attenuation=" command
{
    uint32_t attenuation_mdb;

    attenuation_mdb = psh->uarg[ 0 ];
    ad9361_set_tx_attenuation( ad9361_phy, 1, attenuation_mdb );
    get_tx2_attenuation( NULL );
}

/**************************************************************************//***
 * @brief Gets current TX FIR state.
 *
 * @return None.
 *******************************************************************************/
void get_tx_fir_en( shell *psh ) // "tx_fir_en?" command
{
    uint8_t en_dis;

    ad9361_get_tx_fir_en_dis( ad9361_phy, &en_dis );
    cprintf( "tx_fir_en = %u\n", toUint( en_dis ) );
}

/**************************************************************************//***
 * @brief Sets the TX FIR state.
 *
 * @return None.
 *******************************************************************************/
void set_tx_fir_en( shell *psh ) // "tx_fir_en=" command
{
    uint8_t en_dis;

    en_dis = ( uint8_t ) psh->uarg[ 0 ];
    ad9361_set_tx_fir_en_dis( ad9361_phy, en_dis );
    cprintf( "tx_fir_en = %d\n", toUint( en_dis ) );
}

/**************************************************************************//***
 * @brief Gets current RX LO frequency [ MHz ].
 *
 * @return None.
 *******************************************************************************/
void get_rx_lo_freq( shell *psh ) // "rx_lo_freq?" command
{
    uint64_t lo_freq_hz;

    ad9361_get_rx_lo_freq( ad9361_phy, &lo_freq_hz );
    lo_freq_hz /= 1000000;
    cprintf( "rx_lo_freq = %lu\n", toULint( lo_freq_hz ) );
}

/**************************************************************************//***
 * @brief Sets the RX LO frequency [ MHz ].
 *
 * @return None.
 *******************************************************************************/
void set_rx_lo_freq( shell *psh ) // "rx_lo_freq=" command
{
    uint64_t lo_freq_hz;

    lo_freq_hz = ( uint64_t ) psh->darg[ 0 ];
    lo_freq_hz *= 1000000;
    ad9361_set_rx_lo_freq( ad9361_phy, lo_freq_hz );
    lo_freq_hz /= 1000000;
    cprintf( "rx_lo_freq = %lu\n", toULint( lo_freq_hz ) );
}

/**************************************************************************//***
 * @brief Gets current RX sampling frequency [ Hz ].
 *
 * @return None.
 *******************************************************************************/

u32 AD9364RxSampleFreq( void )
{
    uint32_t sampling_freq_hz = 0;
    ad9361_get_rx_sampling_freq( ad9361_phy, &sampling_freq_hz );
    return ( sampling_freq_hz );
}

void get_rx_samp_freq( shell *psh ) // "rx_samp_freq?" command
{
    uint32_t sampling_freq_hz;

    ad9361_get_rx_sampling_freq( ad9361_phy, &sampling_freq_hz );
    cprintf( "rx_samp_freq=%u\n", toUint( sampling_freq_hz ) );
}

/**************************************************************************//***
 * @brief Sets the RX sampling frequency [ Hz ].
 *
 * @return None.
 *******************************************************************************/
void setRxSampFreq( uint32_t fs )
{
    ad9361_set_rx_sampling_freq( ad9361_phy, fs );
    ad9361_get_rx_sampling_freq( ad9361_phy, &fs );
    dds_update( ad9361_phy );
    cprintf( "rx_samp_freq = %u\n", toUint( fs ) );
}

void set_rx_samp_freq( shell *psh ) // "rx_samp_freq=" command
{
    setRxSampFreq( psh->uarg[ 0 ] );
}

/**************************************************************************//***
 * @brief Gets current RX RF bandwidth [ Hz ].
 *
 * @return None.
 *******************************************************************************/
void get_rx_rf_bandwidth( shell *psh ) // "rx_rf_bandwidth?" command
{
    uint32_t bandwidth_hz;

    ad9361_get_rx_rf_bandwidth( ad9361_phy, &bandwidth_hz );
    cprintf( "rx_rf_bandwidth = %u\n", toUint( bandwidth_hz ) );
}

/**************************************************************************//***
 * @brief Sets the RX RF bandwidth [ Hz ].
 *
 * @return None.
 *******************************************************************************/
void setRxRfBandwidth( uint32_t bw )
{
    ad9361_set_rx_rf_bandwidth( ad9361_phy, bw );
    ad9361_get_rx_rf_bandwidth( ad9361_phy, &bw );
    cprintf( "rx_rf_bandwidth = %u\n", toUint( bw ) );
}

void set_rx_rf_bandwidth( shell *psh ) // "rx_rf_bandwidth=" command
{
    setRxRfBandwidth( psh->uarg[ 0 ] );
}

/**************************************************************************//***
 * @brief Gets current RX1 GC mode.
 *
 * @return None.
 *******************************************************************************/
void get_rx1_gc_mode( shell *psh ) // "rx1_gc_mode?" command
{
    uint8_t gc_mode;

    ad9361_get_rx_gain_control_mode( ad9361_phy, 0, &gc_mode );
    cprintf( "rx1_gc_mode = %u\n", toUint( gc_mode ) );
}

/**************************************************************************//***
 * @brief Sets the RX1 GC mode.
 *
 * @return None.
 *******************************************************************************/
void set_rx1_gc_mode( shell *psh ) // "rx1_gc_mode=" command
{
    uint8_t gc_mode;

    gc_mode = ( uint8_t ) psh->uarg[ 0 ];
    ad9361_set_rx_gain_control_mode( ad9361_phy, 0, gc_mode );
    ad9361_get_rx_gain_control_mode( ad9361_phy, 0, &gc_mode );
    cprintf( "rx1_gc_mode = %u\n", toUint( gc_mode ) );
}

/**************************************************************************//***
 * @brief Gets current RX2 GC mode.
 *
 * @return None.
 *******************************************************************************/
void get_rx2_gc_mode( shell *psh ) // "rx2_gc_mode?" command
{
    uint8_t gc_mode;

    ad9361_get_rx_gain_control_mode( ad9361_phy, 1, &gc_mode );
    cprintf( "rx2_gc_mode = %u\n", toUint( gc_mode ) );
}

/**************************************************************************//***
 * @brief Sets the RX2 GC mode.
 *
 * @return None.
 *******************************************************************************/
void set_rx2_gc_mode( shell *psh ) // "rx2_gc_mode=" command
{
    uint8_t gc_mode;

    gc_mode = ( uint8_t ) psh->uarg[ 0 ];
    ad9361_set_rx_gain_control_mode( ad9361_phy, 1, gc_mode );
    ad9361_get_rx_gain_control_mode( ad9361_phy, 1, &gc_mode );
    cprintf( "rx2_gc_mode = %u\n", toUint( gc_mode ) );
}

/**************************************************************************//***
 * @brief Gets current RX1 RF gain.
 *
 * @return None.
 *******************************************************************************/
void get_rx1_rf_gain( shell *psh ) // "rx1_rf_gain?" command
{
    int32_t gain_db;

    ad9361_get_rx_rf_gain ( ad9361_phy, 0, &gain_db );
    cprintf( "rx1_rf_gain = %i\n", toInt( gain_db ) );
}

/**************************************************************************//***
 * @brief Sets the RX1 RF gain.
 *
 * @return None.
 *******************************************************************************/
void set_rx1_rf_gain( shell *psh ) // "rx1_rf_gain=" command
{
    int32_t gain_db;

    gain_db = ( int32_t ) psh->darg[ 0 ];
    ad9361_set_rx_rf_gain ( ad9361_phy, 0, gain_db );
    ad9361_get_rx_rf_gain ( ad9361_phy, 0, &gain_db );
    cprintf( "rx1_rf_gain = %i\n", toInt( gain_db ) );
}

/**************************************************************************//***
 * @brief Gets current RX2 RF gain.
 *
 * @return None.
 *******************************************************************************/
void get_rx2_rf_gain( shell *psh ) // "rx2_rf_gain?" command
{
    int32_t gain_db;

    ad9361_get_rx_rf_gain ( ad9361_phy, 1, &gain_db );
    cprintf( "rx2_rf_gain = %i\n", toInt( gain_db ) );
}

/**************************************************************************//***
 * @brief Sets the RX2 RF gain.
 *
 * @return None.
 *******************************************************************************/
void set_rx2_rf_gain( shell *psh ) // "rx2_rf_gain=" command
{
    int32_t gain_db;

    gain_db = ( int32_t ) psh->darg[ 0 ];
    ad9361_set_rx_rf_gain ( ad9361_phy, 1, gain_db );
    ad9361_get_rx_rf_gain ( ad9361_phy, 1, &gain_db );
    cprintf( "rx2_rf_gain = %i\n", toInt( gain_db ));
}

/**************************************************************************//***
 * @brief Gets current RX FIR state.
 *
 * @return None.
 *******************************************************************************/
void get_rx_fir_en( shell *psh ) // "rx_fir_en?" command
{
    uint8_t en_dis;

    ad9361_get_rx_fir_en_dis( ad9361_phy, &en_dis );
    cprintf( "rx_fir_en = %u\n", toUint( en_dis ) );
}

/**************************************************************************//***
 * @brief Sets the RX FIR state.
 *
 * @return None.
 *******************************************************************************/
void set_rx_fir_en( shell *psh ) // "rx_fir_en=" command
{
    uint8_t en_dis;

    en_dis = ( uint8_t ) psh->uarg[ 0 ];
    ad9361_set_rx_fir_en_dis( ad9361_phy, en_dis );
    cprintf( "rx_fir_en = %u\n", toUint( en_dis ) );
}

/**************************************************************************//***
 * @brief Gets current DDS TX1 Tone 1 frequency [ Hz ].
 *
 * @return None.
 *******************************************************************************/
void get_dds_tx1_tone1_freq( shell *psh )   // dds_tx1_tone1_freq?
{
    uint32_t freq = dds_st.cached_freq[ DDS_CHAN_TX1_I_F1 ];

    cprintf( "dds_tx1_tone1_freq = %u\n", toUint( freq ) );
}

/**************************************************************************//***
 * @brief Sets the DDS TX1 Tone 1 frequency [ Hz ].
 *
 * @return None.
 *******************************************************************************/
void set_dds_tx1_tone1_freq( shell *psh )   // dds_tx1_tone1_freq=
{
    uint32_t freq = psh->uarg[ 0 ];

    dds_set_frequency( ad9361_phy, DDS_CHAN_TX1_I_F1, freq );
    dds_set_frequency( ad9361_phy, DDS_CHAN_TX1_Q_F1, freq );
    cprintf( "dds_tx1_tone1_freq = %u\n", toUint( freq ) );
}

/**************************************************************************//***
 * @brief Gets current DDS TX1 Tone 2 frequency [ Hz ].
 *
 * @return None.
 *******************************************************************************/
void get_dds_tx1_tone2_freq( shell *psh )   // dds_tx1_tone2_freq?
{
    uint32_t freq = dds_st.cached_freq[ DDS_CHAN_TX1_I_F2 ];

    cprintf( "dds_tx1_tone2_freq = %u\n", toUint( freq ) );
}

/**************************************************************************//***
 * @brief Sets the DDS TX1 Tone 2 frequency [ Hz ].
 *
 * @return None.
 *******************************************************************************/
void set_dds_tx1_tone2_freq( shell *psh )   // dds_tx1_tone2_freq=
{
    uint32_t freq = psh->uarg[ 0 ];

    dds_set_frequency( ad9361_phy, DDS_CHAN_TX1_I_F2, freq );
    dds_set_frequency( ad9361_phy, DDS_CHAN_TX1_Q_F2, freq );
    cprintf( "dds_tx1_tone2_freq = %u\n", toUint( freq ) );
}

/**************************************************************************//***
 * @brief Gets current DDS TX1 Tone 1 phase [ degrees ].
 *
 * @return None.
 *******************************************************************************/
void get_dds_tx1_tone1_phase( shell *psh )  // dds_tx1_tone1_phase?
{
    uint32_t phase = dds_st.cached_phase[ DDS_CHAN_TX1_I_F1 ];

    phase /= 1000;
    cprintf( "dds_tx1_tone1_phase = %u\n", toUint( phase ) );
}

/**************************************************************************//***
 * @brief Sets the DDS TX1 Tone 1 phase [ degrees ].
 *
 * @return None.
 *******************************************************************************/
void set_dds_tx1_tone1_phase( shell *psh )  // dds_tx1_tone1_phase=
{
    int32_t phase = ( int32_t ) psh->darg[ 0 ];

    dds_set_phase( ad9361_phy, DDS_CHAN_TX1_I_F1, ( uint32_t )( phase * 1000 ) );
    if ( ( phase - 90 ) < 0 )
        phase += 360;
    dds_set_phase( ad9361_phy, DDS_CHAN_TX1_Q_F1, ( uint32_t )( ( phase - 90 ) * 1000 ) );
    phase = dds_st.cached_phase[ DDS_CHAN_TX1_I_F1 ] / 1000;
    cprintf( "dds_tx1_tone1_phase = %i\n", toInt( phase ) );
}

/**************************************************************************//***
 * @brief Gets current DDS TX1 Tone 2 phase [ degrees ].
 *
 * @return None.
 *******************************************************************************/
void get_dds_tx1_tone2_phase( shell *psh )  // dds_tx1_tone2_phase?
{
    uint32_t phase = dds_st.cached_phase[ DDS_CHAN_TX1_I_F2 ];

    phase /= 1000;
    cprintf( "dds_tx1_tone2_phase = %u\n", toUint( phase ) );
}

/**************************************************************************//***
 * @brief Sets the DDS TX1 Tone 2 phase [ degrees ].
 *
 * @return None.
 *******************************************************************************/
void set_dds_tx1_tone2_phase( shell *psh )  // dds_tx1_tone2_phase=
{
    int32_t phase = ( int32_t ) psh->darg[ 0 ];

    dds_set_phase( ad9361_phy, DDS_CHAN_TX1_I_F2, ( uint32_t )( phase * 1000 ) );
    if ( ( phase - 90 ) < 0 )
        phase += 360;
    dds_set_phase( ad9361_phy, DDS_CHAN_TX1_Q_F2, ( uint32_t )( ( phase - 90 ) * 1000 ) );
    phase = dds_st.cached_phase[ DDS_CHAN_TX1_I_F2 ] / 1000;
    cprintf( "dds_tx1_tone2_phase = %i\n", toInt( phase ) );
}

/**************************************************************************//***
 * @brief Gets current DDS TX1 Tone 1 scale.
 *
 * @return None.
 *******************************************************************************/
void get_dds_tx1_tone1_scale( shell *psh )  // dds_tx1_tone1_scale?
{
    int32_t scale = dds_st.cached_scale[ DDS_CHAN_TX1_I_F1 ];

    cprintf( "dds_tx1_tone1_scale = %i\n", toInt( scale ) );
}

/**************************************************************************//***
 * @brief Sets the DDS TX1 Tone 1 scale.
 *
 * @return None.
 *******************************************************************************/
void set_dds_tx1_tone1_scale( shell *psh )  // dds_tx1_tone1_scale=
{
    int32_t scale = ( int32_t ) psh->darg[ 0 ];

    dds_set_scale( ad9361_phy, DDS_CHAN_TX1_I_F1, scale );
    dds_set_scale( ad9361_phy, DDS_CHAN_TX1_Q_F1, scale );
    scale = dds_st.cached_scale[ DDS_CHAN_TX1_I_F1 ];
    cprintf( "dds_tx1_tone1_scale=%i\n", toInt( scale ) );
}

/**************************************************************************//***
 * @brief Gets current DDS TX1 Tone 2 scale.
 *
 * @return None.
 *******************************************************************************/
void get_dds_tx1_tone2_scale( shell *psh )  // dds_tx1_tone2_scale?
{
    int32_t scale = dds_st.cached_scale[ DDS_CHAN_TX1_I_F2 ];

    cprintf( "dds_tx1_tone2_scale=%i\n", toInt( scale ) );
}

/**************************************************************************//***
 * @brief Sets the DDS TX1 Tone 2 scale.
 *
 * @return None.
 *******************************************************************************/
void set_dds_tx1_tone2_scale( shell *psh )  // dds_tx1_tone2_scale=
{
    int32_t scale = ( int32_t ) psh->darg[ 0 ];

    dds_set_scale( ad9361_phy, DDS_CHAN_TX1_I_F2, scale );
    dds_set_scale( ad9361_phy, DDS_CHAN_TX1_Q_F2, scale );
    scale = dds_st.cached_scale[ DDS_CHAN_TX1_I_F2 ];
    cprintf( "dds_tx1_tone2_scale = %i\n", toInt( scale ) );
}

/**************************************************************************//***
 * @brief Gets current DDS TX2 Tone 1 frequency [ Hz ].
 *
 * @return None.
 *******************************************************************************/
void get_dds_tx2_tone1_freq( shell *psh )   // dds_tx2_tone1_freq?
{
    uint32_t freq = dds_st.cached_freq[ DDS_CHAN_TX2_I_F1 ];

    cprintf( "dds_tx2_tone1_freq = %i\n", toUint( freq ) );
}

/**************************************************************************//***
 * @brief Sets the DDS TX2 Tone 1 frequency [ Hz ].
 *
 * @return None.
 *******************************************************************************/
void set_dds_tx2_tone1_freq( shell *psh )   // dds_tx2_tone1_freq=
{
    uint32_t freq = psh->uarg[ 0 ];

    dds_set_frequency( ad9361_phy, DDS_CHAN_TX2_I_F1, freq );
    dds_set_frequency( ad9361_phy, DDS_CHAN_TX2_Q_F1, freq );
    cprintf( "dds_tx2_tone1_freq = %u\n", toUint( freq ) );
}

/**************************************************************************//***
 * @brief Gets current DDS TX2 Tone 2 frequency [ Hz ].
 *
 * @return None.
 *******************************************************************************/
void get_dds_tx2_tone2_freq( shell *psh )   // dds_tx2_tone2_freq?
{
    uint32_t freq = dds_st.cached_freq[ DDS_CHAN_TX2_I_F2 ];

    cprintf( "dds_tx2_tone2_freq = %u\n", toUint( freq ) );
}

/**************************************************************************//***
 * @brief Sets the DDS TX2 Tone 2 frequency [ Hz ].
 *
 * @return None.
 *******************************************************************************/
void set_dds_tx2_tone2_freq( shell *psh )   // dds_tx2_tone2_freq=
{
    uint32_t freq = psh->uarg[ 0 ];

    dds_set_frequency( ad9361_phy, DDS_CHAN_TX2_I_F2, freq );
    dds_set_frequency( ad9361_phy, DDS_CHAN_TX2_Q_F2, freq );
    cprintf( "dds_tx2_tone2_freq = %u\n", toUint( freq ) );
}

/**************************************************************************//***
 * @brief Gets current DDS TX2 Tone 1 phase [ degrees ].
 *
 * @return None.
 *******************************************************************************/
void get_dds_tx2_tone1_phase( shell *psh )  // dds_tx2_tone1_phase?
{
    uint32_t phase = dds_st.cached_phase[ DDS_CHAN_TX2_I_F1 ];

    phase /= 1000;
    cprintf( "dds_tx2_tone1_phase = %u\n", toUint( phase ) );
}

/**************************************************************************//***
 * @brief Sets the DDS TX2 Tone 1 phase [ degrees ].
 *
 * @return None.
 *******************************************************************************/
void set_dds_tx2_tone1_phase( shell *psh )  // dds_tx2_tone1_phase=
{
    int32_t phase = psh->uarg[ 0 ];

    dds_set_phase( ad9361_phy, DDS_CHAN_TX2_I_F1, ( uint32_t )( phase * 1000 ) );
    if ( ( phase - 90 ) < 0 )
        phase += 360;
    dds_set_phase( ad9361_phy, DDS_CHAN_TX2_Q_F1, ( uint32_t )( ( phase - 90 ) * 1000 ) );
    phase = dds_st.cached_phase[ DDS_CHAN_TX2_I_F1 ] / 1000;
    cprintf( "dds_tx2_tone1_phase = %i\n", toInt( phase ) );
}

/**************************************************************************//***
 * @brief Gets current DDS TX2 Tone 2 phase [ degrees ].
 *
 * @return None.
 *******************************************************************************/
void get_dds_tx2_tone2_phase( shell *psh )  // dds_tx2_tone2_phase?
{
    uint32_t phase = dds_st.cached_phase[ DDS_CHAN_TX2_I_F2 ];

    phase /= 1000;
    cprintf( "dds_tx2_f2_phase = %u\n", toUint( phase ) );
}

/**************************************************************************//***
 * @brief Sets the DDS TX2 Tone 2 phase [ degrees ].
 *
 * @return None.
 *******************************************************************************/
void set_dds_tx2_tone2_phase( shell *psh )  // dds_tx2_tone2_phase=
{
    int32_t phase = ( int32_t ) psh->darg[ 0 ];

    dds_set_phase( ad9361_phy, DDS_CHAN_TX2_I_F2, ( uint32_t )( phase * 1000 ) );
    if ( ( phase - 90 ) < 0 )
        phase += 360;
    dds_set_phase( ad9361_phy, DDS_CHAN_TX2_Q_F2, ( uint32_t )( ( phase - 90 ) * 1000 ) );
    phase = dds_st.cached_phase[ DDS_CHAN_TX2_I_F2 ] / 1000;
    cprintf( "dds_tx2_tone2_phase = %i\n", toInt( phase ) );
}

/**************************************************************************//***
 * @brief Gets current DDS TX2 Tone 1 scale.
 *
 * @return None.
 *******************************************************************************/
void get_dds_tx2_tone1_scale( shell *psh )  // dds_tx2_tone1_scale?
{
    int32_t scale = dds_st.cached_scale[ DDS_CHAN_TX2_I_F1 ];

    cprintf( "dds_tx2_tone1_scale = %i\n", toInt( scale ) );
}

/**************************************************************************//***
 * @brief Sets the DDS TX2 Tone 1 scale.
 *
 * @return None.
 *******************************************************************************/
void set_dds_tx2_tone1_scale( shell *psh )  // dds_tx2_tone1_scale=
{
    int32_t scale = ( int32_t ) psh->darg[ 0 ];

    dds_set_scale( ad9361_phy, DDS_CHAN_TX2_I_F1, scale );
    dds_set_scale( ad9361_phy, DDS_CHAN_TX2_Q_F1, scale );
    scale = dds_st.cached_scale[ DDS_CHAN_TX2_I_F1 ];
    cprintf( "dds_tx2_tone1_scale = %i\n", toInt( scale ) );
}

/**************************************************************************//***
 * @brief Gets current DDS TX2 Tone 2 scale.
 *
 * @return None.
 *******************************************************************************/
void dds_tx2_tone2_scale( shell *psh )  // dds_tx2_tone2_scale?
{
    int32_t scale = dds_st.cached_scale[ DDS_CHAN_TX2_I_F2 ];

    cprintf( "dds_tx2_tone2_scale = %i\n", toInt( scale ) );
}

/**************************************************************************//***
 * @brief Sets the DDS TX2 Tone 2 scale.
 *
 * @return None.
 *******************************************************************************/
void set_dds_tx2_tone2_scale( shell *psh )  // dds_tx2_tone2_scale=
{
    int32_t scale = ( int32_t ) psh->darg[ 0 ];

    dds_set_scale( ad9361_phy, DDS_CHAN_TX2_I_F2, scale );
    dds_set_scale( ad9361_phy, DDS_CHAN_TX2_Q_F2, scale );
    scale = dds_st.cached_scale[ DDS_CHAN_TX2_I_F2 ];
    cprintf( "dds_tx2_tone2_scale = %i\n", toInt( scale ) );
}



void calib_dc_offset( shell *psh )
{
    if ( !ad9361_calibrate( ad9361_phy, BBDC_CAL ) )
    {
        cprintf( "baseband DC offset calibrated\n");
    }
    else
    {
        cprintf( "baseband DC offset calibration failed\n");

    }
    if ( !ad9361_calibrate( ad9361_phy, RFDC_CAL ) )
    {
        cprintf( "RF DC offset calibrated\n");
    }
    else
    {
        cprintf( "RF DC offset calibration failed\n");
    }
}
