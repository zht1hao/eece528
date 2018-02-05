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

#include "scalar_vec_fir.h"

#define ACCUM_BITS 40

// Scalar FIR filter byte
void scalar_vec_fir_byte(int8_t *output, int8_t *input, const int8_t *const coeffs, const int32_t sample_size, const int32_t num_taps)
{
	int32_t i, j;

	for (j = 0; j <= (sample_size - num_taps); j++) {
		int64_t accum = 0;
		for (i = 0; i < num_taps; i++) {
			int8_t mul_result = input[i] * coeffs[i];
			accum += mul_result;
		}
		accum = accum << (64-ACCUM_BITS);
		accum = accum >> (64-ACCUM_BITS);
		int8_t truncated_accum = accum;
		if(((int64_t)truncated_accum) == accum){
			*output = accum;
		} else {
			*output = 0x7F ^ (accum >> (ACCUM_BITS-1));
		}
		output++;
		input++;
	}
}

// Scalar FIR filter half
void scalar_vec_fir_half(int16_t *output, int16_t *input, const int16_t *const coeffs, const int32_t sample_size, const int32_t num_taps)
{
	int32_t i, j;

	for (j = 0; j <= (sample_size - num_taps); j++) {
		int64_t accum = 0;
		for (i = 0; i < num_taps; i++) {
			int16_t mul_result = input[i] * coeffs[i];
			accum += mul_result;
		}
		accum = accum << (64-ACCUM_BITS);
		accum = accum >> (64-ACCUM_BITS);
		int16_t truncated_accum = accum;
		if(((int64_t)truncated_accum) == accum){
			*output = accum;
		} else {
			*output = 0x7FFF ^ (accum >> (ACCUM_BITS-1));
		}
		output++;
		input++;
	}
}

// Scalar FIR filter word
void scalar_vec_fir_word(int32_t *output, int32_t *input, const int32_t *const coeffs, const int32_t sample_size, const int32_t num_taps)
{
	int32_t i, j;

	for (j = 0; j <= (sample_size - num_taps); j++) {
		int64_t accum = 0;
		for (i = 0; i < num_taps; i++) {
			int32_t mul_result = input[i] * coeffs[i];
			accum += mul_result;
		}
		accum = accum << (64-ACCUM_BITS);
		accum = accum >> (64-ACCUM_BITS);
		int32_t truncated_accum = accum;
		if(((int64_t)truncated_accum) == accum){
			*output = accum;
		} else {
			*output = 0x7FFFFFFF ^ (accum >> (ACCUM_BITS-1));
		}
		output++;
		input++;
	}
}
// Scalar FIR filter unsigned byte
void scalar_vec_fir_ubyte(uint8_t *output, uint8_t *input, const uint8_t *const coeffs, const int32_t sample_size, const int32_t num_taps)
{
	int32_t i, j;

	for (j = 0; j <= (sample_size - num_taps); j++) {
		uint64_t accum = 0;
		for (i = 0; i < num_taps; i++) {
			uint8_t mul_result = input[i] * coeffs[i];
			accum += mul_result;
		}
		accum = accum << (64-ACCUM_BITS);
		accum = accum >> (64-ACCUM_BITS);
		uint8_t truncated_accum = accum;
		if(((uint64_t)truncated_accum) == accum){
			*output = accum;
		} else {
			*output = 0xFF;
		}
		output++;
		input++;
	}
}

// Scalar FIR filter unsigned half
void scalar_vec_fir_uhalf(uint16_t *output, uint16_t *input, const uint16_t *const coeffs, const int32_t sample_size, const int32_t num_taps)
{
	int32_t i, j;

	for (j = 0; j <= (sample_size - num_taps); j++) {
		uint64_t accum = 0;
		for (i = 0; i < num_taps; i++) {
			uint16_t mul_result = input[i] * coeffs[i];
			accum += mul_result;
		}
		accum = accum << (64-ACCUM_BITS);
		accum = accum >> (64-ACCUM_BITS);
		uint16_t truncated_accum = accum;
		if(((uint64_t)truncated_accum) == accum){
			*output = accum;
		} else {
			*output = 0xFFFF;
		}
		output++;
		input++;
	}
}

// Scalar FIR filter unsigned word
void scalar_vec_fir_uword(uint32_t *output, uint32_t *input, const uint32_t *const coeffs, const int32_t sample_size, const int32_t num_taps)
{
	int32_t i, j;

	for (j = 0; j <= (sample_size - num_taps); j++) {
		uint64_t accum = 0;
		for (i = 0; i < num_taps; i++) {
			uint32_t mul_result = input[i] * coeffs[i];
			accum += mul_result;
		}
		accum = accum << (64-ACCUM_BITS);
		accum = accum >> (64-ACCUM_BITS);
		uint32_t truncated_accum = accum;
		if(((uint64_t)truncated_accum) == accum){
			*output = accum;
		} else {
			*output = 0xFFFFFFFF;
		}
		output++;
		input++;
	}
}
