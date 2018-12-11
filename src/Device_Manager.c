/*
 * Device_Manager.c

 *
 *  Created on: 16 Apr 2018
 *      Author: JULIAN MORTIMER
 */
#include "include.h"
#include "device_manager.h"


static const char deviceErrMsgs[][SYS_MAX_ERR_MSG_LEN] = {
		"Name too long",
		"Bad identifier",
		"Device not found",
		"Too many devices of this type",
		"Out of memory",
		"Device ID does not match",
		"No errors found"
		"\0"
};

static DeviceErrCode InitDeviceManager( DeviceManager *pm  );
static u32 IsInitialised( DeviceManager *pm );
static u32 GetDeviceNewId( Device **array );
static Device* AllocateDevice( DeviceManager *pm );
static Device* AddDevice( DeviceManager *pm );
static void SetDeviceErr( DeviceManager *pm, DeviceErrCode e );
static DeviceErrCode GetDeviceErr( DeviceManager *pm );
static Device* GetDevicePtr( DeviceManager *pm, u32 id );
static const char *GetDeviceErrMsg( DeviceManager *pm  );


// This copied and overloaded
const DeviceManager DEVICE_MANAGER_TEMPLATE = {
    NULL,
    NULL,
    DEVICE_SUCCESS,
    false,
    InitDeviceManager,
    IsInitialised,
    GetDeviceNewId,
    AllocateDevice,
    AddDevice,
    SetDeviceErr,
    GetDeviceErr,
    GetDevicePtr,
    GetDeviceErrMsg,
};


// This never used
static const Device DEVICE_TEMPLATE = {
	"\0",
	0,
	0
};




// Accessible only from this module


// little do do here. Clear the pointer array.
static DeviceErrCode InitDeviceManager( DeviceManager *pm  )
{
	for (u32 i = 0; i < SYS_MAX_DEVICES; i++ )
    {
    	pm->DeviceArray[i] = NULL;
    }
    pm->isInitialised = true;
	return( DEVICE_SUCCESS );
}

static u32 IsInitialised( DeviceManager *pm  )
{
	return pm->isInitialised;
}

// add a device
static u32 GetDeviceNewId( Device **array )
{
	u32 i = 0;
	while(i < SYS_MAX_DEVICES )
	{
		if ( !array[i] )
		{
			return( i );
		}
	}
	return ( -1 );
}

//This is always overloaded.
static Device* AllocateDevice( DeviceManager *pm )
{
	Device *pd  = ( Device* ) malloc(sizeof( Device ));
	if( !pd ) return ( Device * )NULL;
    *pd = DEVICE_TEMPLATE;
	return ( Device *) pd;
}


// This is always the same
static Device* AddDevice( DeviceManager *pm )
{
    u32 id = pm->GetDeviceNewId( pm->DeviceArray );
    if (id < 0 )
    {
    	pm->SetDeviceErr( pm, DEVICE_NO_MORE_ALLOWED);
    	return( NULL );
    }
    Device *pd = pm->AllocateDevice( pm );
    if (!pd)
    {
	    pm->SetDeviceErr( pm, DEVICE_OUT_OF_MEMORY );
	    return( NULL );
    }
  	pd->id = id;
    pm->DeviceArray[ id ] = pd;
    return( pd );
}


static void SetDeviceErr( DeviceManager *pm, DeviceErrCode e )
{
	pm->errCode = e;
}

static DeviceErrCode GetDeviceErr( DeviceManager *pm )
{
	DeviceErrCode e = pm->errCode;
	pm->errCode = DEVICE_SUCCESS;
	return( e );
}

// Obtain a pointer to a device - Should be recast by consumer
static Device* GetDevicePtr( DeviceManager *pm, u32 id )
{
    if (( id < 0 ) || ( id >= SYS_MAX_DEVICES - 1 )) {
    	pm->SetDeviceErr( pm, DEVICE_BAD_ID ); return ( NULL );
    }
    if ( pm->DeviceArray[ id ] == NULL) {
    	pm->SetDeviceErr( pm, DEVICE_NOT_FOUND ); return ( NULL );
    }
    return( pm->DeviceArray[ id ]);
}


const char *GetDeviceErrMsg( DeviceManager *pm  )
{
	return( deviceErrMsgs[ GetDeviceErr( pm ) ]);
}

