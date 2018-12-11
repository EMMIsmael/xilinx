/*
 * Clk_Wiz.h
 *
 *  Created on: 13 Jun 2017
 *      Author: JULIAN MORTIMER
 */

#ifndef CLK_WIZ_H_
#define CLK_WIZ_H_



#define CLK_WIZ_SRR_REG_OFS   ( 0x000 )
#define CLK_WIZ_STAT_REG_OFS  ( 0x004 )
#define CLK_WIZ_CFG0_REG_OFS  ( 0x200 )
#define CLK_WIZ_CFG2_REG_OFS  ( 0x208 )
#define CLK_WIZ_CFG23_REG_OFS ( 0x25c )

#define CLK_WIZ_SOFT_RESET    ( 0x00a )
#define CLK_WIZ_LOCK          ( 0x001 )
#define CLK_WIZ_SET_LD_SEN    ( 0x007 )
#define CLK_WIZ_CLR_LD_SEN    ( 0x002 )




typedef struct sClkWiz    tClkWiz;
typedef struct sClkWizCfg tClkWizCfg;

struct sClkWiz {
	char name[SYS_MAX_OBJ_NAME_LEN+1];
	u32 id;
    u32 baseAddr;
    const tClkWizCfg *pcfg;
};

struct sClkWizCfg {
	u32 frqMHz;
    u32 divclkDivide;
    u32 clkfboutMult;
    u32 clkfboutFrac;
    u32 clkout0Divide;
    u32 clkout0Frac;
};




#endif /* CLK_WIZ_H_ */
