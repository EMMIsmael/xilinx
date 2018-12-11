/***************************************************************************//**
 *   @file   console.c
 *   @brief  Implementation of Console Driver.
 *   @author DBogdan (dragos.bogdan@analog.com)
********************************************************************************
 * Copyright 2013(c) Analog Devices, Inc.
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
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, INTELLECTUAL PROPERTY RIGHTS, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*******************************************************************************/

/******************************************************************************/
/***************************** Include Files **********************************/
/******************************************************************************/
#include "stdarg.h"
#include "stdlib.h"
#include "stdio.h"
#include "console.h"


/***************************************************************************//**
 * @brief Reads one character from UART.
 *
 * @return received_char   - Read character.
 *
 * Note: Blocking function - Waits until get a valid data.
*******************************************************************************/
void uart_read_char(char * data)
{
	*data = getchar();
}


/***************************************************************************//**
 * @brief Converts an integer number to a string of ASCII characters string.
 *
 * @param number - Integer number.
 * @param base   - Numerical base used to represent the integer number as
 *                 string.
 *
 * @return Pointer to the string of ASCII characters.
*******************************************************************************/
char *int_to_str(long number, char base)
{
	unsigned long pos_number = 0;
	char		  neg_sign   = 0;
	const char	  digits[17] = "0123456789ABCDEF";
	static char   buffer[17] = "                ";
	char*		  buffer_ptr = &buffer[16];

	if((number < 0) && (base == 10))
	{
		neg_sign = 1;
		pos_number = -1 * number;
	}
	else
	{
		pos_number = (unsigned long)number;
	}
	do
	{
		*buffer_ptr-- = digits[pos_number % base];
		pos_number /= base;
	}
	while(pos_number != 0);

	if(neg_sign)
	{
		*buffer_ptr-- = '-';
	}

    buffer_ptr++;

	return buffer_ptr;
}

