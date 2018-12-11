/*
 * GM_AXI_BRAM_Controller_h.h
 *
 *  Created on: 20 May 2017
 *      Author: JULIAN MORTIMER
 */

#ifndef GM_AXI_BRAM_CONTROLLER_0_H_H_
#define GM_AXI_BRAM_CONTROLLER_0_H_H_

#define AXIBRAMC_N_SLAVE_REGS          ( 8 )



#define AXIBRAMC_CTRL_REG              ( 0 )
#define AXIBRAMC_WR_ADDR_REG           ( 2 )
#define AXIBRAMC_WR_DATA_REG           ( 3 )
#define AXIBRAMC_RD_ADDR_REG           ( 4 )

#define AXIBRAMC_INFO_REG              ( 1 )
#define AXIBRAMC_BRAM_SIZE_REG         ( 5 )
#define AXIBRAMC_STAT_REG              ( 6 )
#define AXIBRAMC_RD_DATA_REG           ( 7 )

#define AXIBRAMC_RESET        ( 0x00000001 )
#define AXIBRAMC_WR_DATA      ( 0x00000002 )
#define AXIBRAMC_END_WR       ( 0x00000004 )
#define AXIBRAMC_RD_DATA      ( 0x00000008 )
#define AXIBRAMC_END_RD       ( 0x00000010 )

#define AXIBRAMC_HW_INFO_REG  ( 0x1BEB8C08 )

typedef struct sAxiBRAMC AxiBRAMC;

typedef enum e_axibramcerr     { ABC_NAME_TOO_LONG,
                                 ABC_SUCCESS
                               } axiBRAMCErr;

#endif /* GM_AXI_BRAM_CONTROLLER_H_H_ */
