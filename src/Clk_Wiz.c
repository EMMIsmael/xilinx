/*
 * Clk_Wiz.c
 *
 *  Created on: 13 Jun 2017
 *      Author: JULIAN MORTIMER
 */
#include "include.h"

void clkWizReset(tClkWiz *ps);
void clkWizWaitLock(tClkWiz *ps);
void clkWizConfig(tClkWiz *ps);
int clkWizSetFreq(tClkWiz *ps, u32 frqMhz);
int clkWizInit(tClkWiz *ps, u32 baseAddr, u32 frqMHz, u32 id, const char name[]);

tClkWiz clkWizDac;

const tClkWizCfg clkWizCfgTab[] =
{
  { 16, 1, 10,   0, 62, 500 },
  { 20, 1,  8, 500, 42, 500 },
  { 24, 5, 50, 250, 41, 875 },
  { 28, 5, 49, 875, 35, 625 },
  { 40, 1, 10,   0, 25,   0 },
  { 56, 5, 50, 750, 18, 125 },
  { 0 }
};


void clkWizReset(tClkWiz *ps)
{
	hwWriteReg(ps->baseAddr, CLK_WIZ_SRR_REG_OFS, CLK_WIZ_SOFT_RESET);
    clkWizWaitLock(ps);
}

void clkWizWaitLock(tClkWiz *ps)
{
    u32 t = 100000;
	while (((hwReadReg(ps->baseAddr, CLK_WIZ_STAT_REG_OFS) & CLK_WIZ_LOCK) == 0) && (t > 0))
		t--;
	if (t == 0) printf("\nTimed Out Waiting For Clock Wizard %u (""%s"") to Lock\n", (unsigned int)ps->id, ps->name);
}

void clkWizConfig(tClkWiz *ps)
{
	const tClkWizCfg *pc = ps->pcfg;
    u32 clkConfig0 = (pc->divclkDivide  & 0x0ff) + ((pc->clkfboutMult & 0x0ff) << 8) + ((pc->clkfboutFrac & 0x3ff) << 16);
    u32 clkConfig2 = (pc->clkout0Divide & 0x0ff) + ((pc->clkout0Frac  & 0x3ff) << 8);
    if (pc->clkfboutFrac != 0) clkConfig0 += (1 << 26);
    if (pc->clkout0Frac != 0)  clkConfig2 += (1 << 18);
    clkWizReset(ps);
    hwWriteReg(ps->baseAddr, CLK_WIZ_CFG0_REG_OFS,  clkConfig0);
    hwWriteReg(ps->baseAddr, CLK_WIZ_CFG2_REG_OFS,  clkConfig2);
    hwWriteReg(ps->baseAddr, CLK_WIZ_CFG23_REG_OFS, CLK_WIZ_SET_LD_SEN);
    hwWriteReg(ps->baseAddr, CLK_WIZ_CFG23_REG_OFS, CLK_WIZ_CLR_LD_SEN);
    clkWizWaitLock(ps);
}

int clkWizSetFreq(tClkWiz *ps, u32 frqMHz)
{
    const tClkWizCfg *pc = clkWizCfgTab;
    while (pc->frqMHz != 0)
    {
    	if (pc->frqMHz == frqMHz)
    	{
    		ps->pcfg = pc;
    		break;
    	}
    	pc++;
    }
    if (pc->frqMHz == 0)
    {
    	printf("\nFrequency %u For Clock Wizard %u (""%s"") not found\n", (unsigned int)frqMHz, (unsigned int)ps->id, ps->name);
    	return FALSE;
    }
    clkWizConfig(ps);
    return TRUE;
}

int clkWizInit(tClkWiz *ps, u32 baseAddr, u32 frqMHz, u32 id, const char name[])
{
    ps->id = id;
    ps->baseAddr = baseAddr;

    if (strlen(name) <= SYS_MAX_OBJ_NAME_LEN)
    {
    	strcpy(ps->name, name);
    }
    else
    {
    	printf("\nName %s too long for Clock Wizard object\n", name);
    	return FALSE;
    }

    return (clkWizSetFreq(ps, frqMHz));
}

