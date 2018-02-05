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

/**@file*/
#include "scalar_mtx_fir.h"

/** 2D FIR filter scalar byte.
 *	Scan through the input multiplying
 * 	and accumulating with the filter
 *
 * @param[in] input.
 * @param[out] output.
 * @param[in] coeffs.
 * @param[in] num_row.
 * @param[in] num_column.
 * @param[in] ntaps_row.
 * @param[in] ntaps_column.
 *
 */
void scalar_mtx_2Dfir_byte(
		int8_t *output, int8_t *input, const int32_t *coeffs,
		const int num_row, const int num_column, const int ntaps_row, const int ntaps_column )
{
	int l,k,j,i;
	int32_t temp1, temp2;
	for( l = 0 ; l < num_row-ntaps_row; l++) {
		for( k = 0; k < num_column-ntaps_column; k++) {
			temp1 = 0;
			for( j = 0; j < ntaps_row; j++) {
				temp2 = 0;
				for( i = 0; i < ntaps_column; i++) {
					temp2 += input[(l+j)*num_column+(k+i)] * coeffs[j*ntaps_column+i];
				}
				temp1 += temp2;
			}
			output[l*num_column+k] = temp1 >> 8;
		}
	}
}

/** 2D FIR filter scalar half.
 *	Scan through the input multiplying
 * 	and accumulating with the filter
 *
 * @param[in] input.
 * @param[out] output.
 * @param[in] coeffs.
 * @param[in] num_row.
 * @param[in] num_column.
 * @param[in] ntaps_row.
 * @param[in] ntaps_column.
 *
 */
void scalar_mtx_2Dfir_half(
		int16_t *output, int16_t *input, const int32_t *coeffs,
		const int num_row, const int num_column, const int ntaps_row, const int ntaps_column )
{
	int l,k,j,i;
	int32_t temp1, temp2;
	for( l = 0 ; l < num_row-ntaps_row; l++) {
		for( k = 0; k < num_column-ntaps_column; k++) {
			temp1 = 0;
			for( j = 0; j < ntaps_row; j++) {
				temp2 = 0;
				for( i = 0; i < ntaps_column; i++) {
					temp2 += input[(l+j)*num_column+(k+i)] * coeffs[j*ntaps_column+i];
				}
				temp1 += temp2;
			}
			output[l*num_column+k] = temp1 >> 8;
		}
	}
}

/** 2D FIR filter scalar word.
 *	Scan through the input multiplying
 * 	and accumulating with the filter
 *
 * @param[in] input.
 * @param[out] output.
 * @param[in] coeffs.
 * @param[in] num_row.
 * @param[in] num_column.
 * @param[in] ntaps_row.
 * @param[in] ntaps_column.
 *
 */
void scalar_mtx_2Dfir_word(
		int32_t *output, int32_t *input, const int32_t *coeffs,
		const int num_row, const int num_column, const int ntaps_row, const int ntaps_column )
{
	int l,k,j,i;
	int32_t temp1, temp2;
	for( l = 0 ; l < num_row-ntaps_row; l++) {
		for( k = 0; k < num_column-ntaps_column; k++) {
			temp1 = 0;
			for( j = 0; j < ntaps_row; j++) {
				temp2 = 0;
				for( i = 0; i < ntaps_column; i++) {
					temp2 += input[(l+j)*num_column+(k+i)] * coeffs[j*ntaps_column+i];
				}
				temp1 += temp2;
			}
			output[l*num_column+k] = temp1 >> 8;
		}
	}
}

/** 2D FIR filter scalar byte.
 *	Scan through the input multiplying
 * 	and accumulating with the filter
 *
 * @param[in] input.
 * @param[out] output.
 * @param[in] coeffs.
 * @param[in] num_row.
 * @param[in] num_column.
 * @param[in] ntaps_row.
 * @param[in] ntaps_column.
 *
 */
void scalar_mtx_2Dfir_ubyte(
		uint8_t *output, uint8_t *input, const uint32_t *coeffs,
		const int num_row, const int num_column, const int ntaps_row, const int ntaps_column )
{
	int l,k,j,i;
	uint32_t temp1, temp2;
	for( l = 0 ; l < num_row-ntaps_row; l++) {
		for( k = 0; k < num_column-ntaps_column; k++) {
			temp1 = 0;
			for( j = 0; j < ntaps_row; j++) {
				temp2 = 0;
				for( i = 0; i < ntaps_column; i++) {
					temp2 += input[(l+j)*num_column+(k+i)] * coeffs[j*ntaps_column+i];
				}
				temp1 += temp2;
			}
			output[l*num_column+k] = temp1 >> 8;
		}
	}
}

/** 2D FIR filter scalar half.
 *	Scan through the input multiplying
 * 	and accumulating with the filter
 *
 * @param[in] input.
 * @param[out] output.
 * @param[in] coeffs.
 * @param[in] num_row.
 * @param[in] num_column.
 * @param[in] ntaps_row.
 * @param[in] ntaps_column.
 *
 */
void scalar_mtx_2Dfir_uhalf(
		uint16_t *output, uint16_t *input, const uint32_t *coeffs,
		const int num_row, const int num_column, const int ntaps_row, const int ntaps_column )
{
	int l,k,j,i;
	uint32_t temp1, temp2;
	for( l = 0 ; l < num_row-ntaps_row; l++) {
		for( k = 0; k < num_column-ntaps_column; k++) {
			temp1 = 0;
			for( j = 0; j < ntaps_row; j++) {
				temp2 = 0;
				for( i = 0; i < ntaps_column; i++) {
					temp2 += input[(l+j)*num_column+(k+i)] * coeffs[j*ntaps_column+i];
				}
				temp1 += temp2;
			}
			output[l*num_column+k] = temp1 >> 8;
		}
	}
}

/** 2D FIR filter scalar word.
 *	Scan through the input multiplying
 * 	and accumulating with the filter
 *
 * @param[in] input.
 * @param[out] output.
 * @param[in] coeffs.
 * @param[in] num_row.
 * @param[in] num_column.
 * @param[in] ntaps_row.
 * @param[in] ntaps_column.
 *
 */
void scalar_mtx_2Dfir_uword(
		uint32_t *output, uint32_t *input, const uint32_t *coeffs,
		const int num_row, const int num_column, const int ntaps_row, const int ntaps_column )
{
	int l,k,j,i;
	uint32_t temp1, temp2;
	for( l = 0 ; l < num_row-ntaps_row; l++) {
		for( k = 0; k < num_column-ntaps_column; k++) {
			temp1 = 0;
			for( j = 0; j < ntaps_row; j++) {
				temp2 = 0;
				for( i = 0; i < ntaps_column; i++) {
					temp2 += input[(l+j)*num_column+(k+i)] * coeffs[j*ntaps_column+i];
				}
				temp1 += temp2;
			}
			output[l*num_column+k] = temp1 >> 8;
		}
	}
}

/**@}*/
