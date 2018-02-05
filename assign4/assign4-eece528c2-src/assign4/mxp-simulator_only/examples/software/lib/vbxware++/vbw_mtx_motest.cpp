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

#include "Vector.hpp"
#include "vbw_mtx_motest.hpp"

namespace vbw{

	int mtx_motest_byte(uint32_t *res, uint8_t *x, uint8_t *y,
	                        const int32_t search_height, const int32_t search_width,
	                        const int32_t block_height, const int32_t block_width, const int32_t image_width)
	{
		using namespace VBX;

		Vector<vbx_ubyte_t,2> ref_block(block_width,block_height,block_width);
		Vector<vbx_ubyte_t,2> search_block(block_width+search_width,block_height+search_height,block_width+search_width);
		Vector<vbx_ubyte_t,2>  sub_block(block_width,block_height+search_height,block_width);
		Vector<vbx_uword_t,2> result(search_width,search_height,search_width);



		ref_block.dma_read(x,image_width);

		search_block.dma_read(y,image_width);

		for(int col=0;col<search_width;col++){
			//copy the sub window into a packed vector
			sub_block=search_block[col upto block_width+col , 0 upto block_height+search_height];

			result[col upto col+1,0 upto search_height] =
				accumulate(absdiff(sub_block.to3D(block_width*block_height,1,0,search_height,block_width),
				                   ref_block.to3D(block_width*block_height,1,0,search_height,0)));
		}
		result.dma_write(res,search_width);
		return 0;
	}

	int mtx_motest_byte(output_type *result, input_type *x, input_type *y, vbw_motest_t *m)
	{
		return mtx_motest_byte( (uint32_t*)result, (uint8_t*)x, (uint8_t*)y, m->search_height, m->search_width, m->block_height,
		                        m->block_width, m->image_width);

	}

}
