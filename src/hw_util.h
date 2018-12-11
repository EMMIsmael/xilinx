/*
 * hw_util.h
 *
 *  Created on: 16 May 2017
 *      Author: JULIAN MORTIMER
 */

#ifndef HW_UTIL_H_
#define HW_UTIL_H_

#define hwWriteReg(BaseAddress, RegOffset, Data) \
  	Xil_Out32((BaseAddress) + (RegOffset), (u32)(Data))

#define hwReadReg(BaseAddress, RegOffset) \
    Xil_In32((BaseAddress) + (RegOffset))

#define ID_NEW_DEVICE   ( 0x5AA51CBA )
#define ID_DEVICE_ERROR ( 0xFAECE5 )

#define HW_READ_BRAM_TIMEOUT ( 0.001 )
#endif /* HW_UTIL_H_ */



