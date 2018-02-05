/* VECTORBLOX MXP SOFTWARE DEVELOPMENT KIT
 *
 * Copyright (C) 2012-2017 VectorBlox Computing Inc., Vancouver, British Columbia, Canada.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *
 *     * Neither the name of VectorBlox Computing Inc. nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * This agreement shall be governed in all respects by the laws of the Province
 * of British Columbia and by the laws of Canada.
 *
 * This file is part of the VectorBlox MXP Software Development Kit.
 *
 */

#include "scalar_vec_copy.h"

// Scalar reverse byte
void scalar_vec_copy_ubyte(uint8_t *output, uint8_t *input, const int N)
{
	int i;
    for (i = 0; i < N; i++) {
        output[i]=input[i];
	}
}

// Scalar reverse half
void scalar_vec_copy_uhalf(uint16_t *output, uint16_t *input, const int N)
{
	int i;
    for (i = 0; i < N; i++) {
        output[i]=input[i];
	}
}

// Scalar reverse word
void scalar_vec_copy_uword(uint32_t *output, uint32_t *input, const int N)
{
	int i;
    for (i = 0; i < N; i++) {
        output[i]=input[i];
	}
}

void scalar_vec_copy_byte(int8_t *output, int8_t *input, const int N)
{
    scalar_vec_copy_ubyte((uint8_t *)output, (uint8_t *)input, N);
}

void scalar_vec_copy_half(int16_t *output, int16_t *input, const int N)
{
    scalar_vec_copy_uhalf((uint16_t *)output, (uint16_t *)input, N);
}

void scalar_vec_copy_word(int32_t *output, int32_t *input, const int N)
{
    scalar_vec_copy_uword((uint32_t *)output, (uint32_t *)input, N);
}

