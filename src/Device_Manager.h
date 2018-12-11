/*
 * Device_Manager.h
 *
 *  Created on: 16 Apr 2018
 *      Author: JULIAN MORTIMER
 */

#ifndef DEVICE_MANAGER_H_
#define DEVICE_MANAGER_H_



typedef struct sDeviceManager DeviceManager;
typedef struct sDeviceAPI     DeviceAPI;
typedef struct sDevice        Device;
typedef struct sDeviceInfo    DeviceInfo;
extern  const  DeviceManager DEVICE_MANAGER_TEMPLATE;



//Always Overloaded
struct sDeviceAPI {
	DeviceManager *pDeviceManager;
};

//Always Overloaded
struct sDevice {
	char     name[SYS_MAX_OBJ_NAME_LEN + 1];
	u32      baseAddr;
	u32      id;
};

struct sDeviceManager {
  DeviceAPI         *pAPI;
  Device**          DeviceArray;
  DeviceErrCode     errCode;
  u32               isInitialised;
  DeviceErrCode     (* InitDeviceManager )( DeviceManager *pm  );
  u32               (* IsInitialised )( DeviceManager *pm );
  u32               (* GetDeviceNewId )( Device **array );
  Device*           (* AllocateDevice )( DeviceManager *pm );
  Device*           (* AddDevice )( DeviceManager *pm );
  void              (* SetDeviceErr )( DeviceManager *p, DeviceErrCode e );
  DeviceErrCode     (* GetDeviceErr )( DeviceManager *pm );
  Device*           (* GetDevicePtr )( DeviceManager *p, u32 id );
  const char *      (* GetDeviceErrMsg )( DeviceManager *pm  );
};

#endif /* DEVICE_MANAGER_H_ */
