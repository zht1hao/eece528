/* VECTORBLOX MXP SOFTWARE DEVELOPMENT KIT
 *
 * Copyright (C) 2012-2016 VectorBlox Computing Inc., Vancouver, British Columbia, Canada.
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

#include "lbp_detect.h"

#define VLBPLUT VCUSTOM0
#define VLBPPAT VCUSTOM1

#define VLBPLUTWR VCUSTOM0
#define VLBPLUTZERO VCUSTOM1
#define VLBPLUTINC VCUSTOM2
#define VLBPCELL1 VCUSTOM3
#define VLBPCELL2 VCUSTOM4
#define VLBPCELL4 VCUSTOM5
#define DOUBLE_LUT 1

unsigned short** ScalarLBPRestrictedSums(unsigned short* img, const unsigned width, const unsigned height, const unsigned log)
{
	int i, j, l, n, m, offset, sum;
    unsigned short **sums = (unsigned short**)malloc((log+1)*sizeof(unsigned short*));

    for (l=0; l < log+1; l++) {
        sums[l] = (unsigned short*)vbx_shared_malloc(height*width*sizeof(unsigned short));
    }

    for(l = 0; l < log+1; l++) {
        offset = (1 << l);
        for (j = 0; j < height - (offset-1); j++) {
            for (i = 0; i < width - (offset-1); i++) {
                sum = 0;
                for (n = 0; n < offset; n++) {
                    for (m = 0; m < offset; m++) {
                        sum += img[(j+n)*width+(i+m)];
                    }
                }
                sums[l][j*width + i] = sum;
            }
        }
    }
    return sums;
}


unsigned char** ScalarLBPRestrictedPatterns(unsigned short **sums, const unsigned width, const unsigned height, const unsigned log)
{
	int i, j, l, n, cell, center, neighbor;
    unsigned char pattern;

    int coeff_x[8] = {0, 1, 2, 2, 2, 1, 0, 0};
    int coeff_y[8] = {0, 0, 0, 1, 2, 2, 2, 1};

    unsigned char **patterns = (unsigned char**)malloc((log+1)*sizeof(unsigned char*));

    for (l=0; l<log+1; l++) {
        patterns[l] = (unsigned char*)vbx_shared_malloc(height*width*sizeof(unsigned char));
    }

    for(l = 0; l < log+1; l++) {
        cell = 1 << l;

        for (j = 0; j < height - cell*2; j++) {
            for (i = 0; i < width - cell*2; i++) {
                pattern = 0;
                center = sums[l][(j+cell)*width + (i+cell)];

                for (n = 0; n < 8; n++) {
                    neighbor = sums[l][(j+coeff_y[n]*cell)*width + (i+coeff_x[n]*cell)];

                    if (neighbor >= center) {
                        pattern += 1 << (7-n);
                    }
                }
                patterns[l][j*width+i] = pattern;
            }
        }
    }
    return patterns;
}


unsigned char** LBPRestrictedSums(unsigned short* img, const unsigned width, const unsigned height, const unsigned log)
{
	vbx_mxp_t *this_mxp = VBX_GET_THIS_MXP();
	int i, j, offset, rows, max_rows;


    unsigned short *tmp = (unsigned short*)vbx_shared_malloc(height*width*sizeof(unsigned short));
    unsigned char **sums = (unsigned char**)malloc((log+1)*sizeof(unsigned char*));
    for (i=0; i<log+1; i++) {
        sums[i] = (unsigned char*)vbx_shared_malloc(height*width*sizeof(unsigned char));
    }

    max_rows = this_mxp->scratchpad_size/(width*(sizeof(vbx_uhalf_t)+sizeof(vbx_ubyte_t)));
    if (max_rows > height) {
        max_rows = height;
    }
    vbx_uhalf_t *v_x = (vbx_uhalf_t*)vbx_sp_malloc(max_rows*width*sizeof(vbx_uhalf_t));
    vbx_ubyte_t *v_b = (vbx_ubyte_t*)vbx_sp_malloc(max_rows*width*sizeof(vbx_ubyte_t));

    for(i = 0; i < log; i++) {
        offset = 1<<i;
        rows = max_rows - offset;

        for (j = 0; j < height; j = j+rows) {
            if(j+rows >= height) {
                rows = height - j;
            }
            if (i == 0) {
                vbx_dma_to_vector(v_x, img + j*width, (rows+offset)*width*sizeof(unsigned short));
                vbx_set_vl(width*rows);
                vbx(VVHBU, VMOV, v_b, v_x, 0);
                vbx_dma_to_host(sums[0]+ j*width, v_b, rows*width*sizeof(vbx_ubyte_t));
            } else {
                vbx_dma_to_vector(v_x, tmp + j*width, (rows+offset)*width*sizeof(unsigned short));
            }
            vbx_set_vl(width);
            /* requires loading in *offset* extra rows */
            vbx_set_2D(rows+offset, width*sizeof(vbx_half_t), width*sizeof(vbx_half_t), width*sizeof(vbx_half_t));
            vbx_2D(VVHU, VADD, v_x, v_x, v_x + offset);
            vbx_set_2D(rows, width*sizeof(vbx_half_t), width*sizeof(vbx_half_t), width*sizeof(vbx_half_t));
            vbx_2D(VVHU, VADD, v_x, v_x, v_x + offset*width);
            if (i < log - 1) {
                vbx_dma_to_host(tmp + j*width, v_x, rows*width*sizeof(vbx_uhalf_t));
            }
            vbx_set_2D(rows, width*sizeof(vbx_byte_t), 0, width*sizeof(vbx_half_t));
            vbx_2D(SVHBU, VSHR, v_b, 2*offset, v_x);
            vbx_dma_to_host(sums[i+1]+ j*width, v_b, rows*width*sizeof(vbx_ubyte_t));
        }
    }
    vbx_sync();
    vbx_sp_free();
    vbx_shared_free(tmp);

    return sums;
}

unsigned char** LBPRestrictedPatterns(unsigned char **sums, const unsigned width, const unsigned height, const unsigned log)
{
	vbx_mxp_t *this_mxp = VBX_GET_THIS_MXP();
	int i, j, n, rows, max_rows, cell;

    int coeff_x[8] = {0, 1, 2, 2, 2, 1, 0, 0};
    int coeff_y[8] = {0, 0, 0, 1, 2, 2, 2, 1};

    unsigned char **patterns = (unsigned char**)malloc((log+1)*sizeof(unsigned char*));

    for (i=0; i<log+1; i++) {
        patterns[i] = (unsigned char*)vbx_shared_malloc(height*width*sizeof(unsigned char));
    }

    max_rows = this_mxp->scratchpad_size/(4*width*sizeof(vbx_ubyte_t));
    if (max_rows > height) {
        max_rows = height;
    }

    vbx_ubyte_t *v_row =   (vbx_ubyte_t*)vbx_sp_malloc(max_rows*width*sizeof(vbx_ubyte_t));
    vbx_ubyte_t *v_cmp =   (vbx_ubyte_t*)vbx_sp_malloc(max_rows*width*sizeof(vbx_ubyte_t));
    vbx_ubyte_t *v_shift = (vbx_ubyte_t*)vbx_sp_malloc(max_rows*width*sizeof(vbx_ubyte_t));
    vbx_ubyte_t *v_lbp =   (vbx_ubyte_t*)vbx_sp_malloc(max_rows*width*sizeof(vbx_ubyte_t));

    for(i = 0; i < log+1; i++) {
        cell = 1<<i;
        rows = max_rows - 2*(1<<i);
        for (j = 0; j < height; j = j + rows) {
            if(j+rows >= height) {
                rows = height - j;
            }
            vbx_dma_to_vector(v_row, sums[i]+j*width, (rows+2*cell)*width*sizeof(vbx_ubyte_t));
            vbx_set_vl(width);
            vbx_set_2D(rows, width*sizeof(vbx_byte_t), width*sizeof(vbx_byte_t), width*sizeof(vbx_byte_t));

            vbx_2D(SVBU, VMOV, v_lbp, 0, 0);

            for (n = 0; n < 8; n++) {
                vbx_2D(SVBU, VMOV, v_shift, 0, 0);
                vbx_2D(VVBU, VSUB, v_cmp, v_row + coeff_x[n]*cell + coeff_y[n]*cell*width, v_row + width*cell + cell);
                vbx_2D(SVBU, VCMV_GEZ, v_shift, 1<<(7-n), v_cmp);
                vbx_2D(VVBU, VADD, v_lbp, v_lbp, v_shift);
            }
            vbx_dma_to_host(patterns[i]+j*width, v_lbp, rows*width*sizeof(vbx_ubyte_t));
        }
    }
    vbx_sync();
    vbx_sp_free();

    return patterns;
}

unsigned short** LBPRestrictedSums2(unsigned short* img, const unsigned width, const unsigned height, const unsigned log)
{
	int i, j, k, offset;

    unsigned short **sums = (unsigned short**)malloc((log+1)*sizeof(unsigned short*));
    for (i=0; i<log+1; i++) {
        sums[i] = (unsigned short*)vbx_shared_malloc(height*width*sizeof(unsigned short));
    }
    vbx_uhalf_t **v_x = (vbx_uhalf_t**)vbx_shared_malloc(3*sizeof(vbx_uhalf_t*));
    for (i=0; i<log+1; i++) {
        v_x[i] = (vbx_uhalf_t*)vbx_sp_malloc(width*sizeof(vbx_uhalf_t));
    }
    vbx_uhalf_t *v_in = (vbx_uhalf_t*)vbx_sp_malloc(width*sizeof(vbx_uhalf_t));
    vbx_uhalf_t *v_out = (vbx_uhalf_t*)vbx_sp_malloc(width*sizeof(vbx_uhalf_t));
    vbx_uhalf_t *v_tmp;

    for(i = 0; i < log; i++) {
        offset = 1<<i;
        vbx_set_vl(width - (offset-1));

        for(k=0; k < offset+1; k++) {
            if(i == 0) {
                vbx_dma_to_vector(v_x[k], img + k*width, width*sizeof(unsigned short));
            } else {
                vbx_dma_to_vector(v_x[k], sums[i] + k*width, width*sizeof(unsigned short));
            }
        }

        for(k=0; k < offset+1; k++) {
            vbx(VVHU, VADD, v_x[k], v_x[k], v_x[k] + offset);
        }

        for (j = 0; j < height-(offset-1); j++) {
            if(i == 0) {
                vbx_dma_to_vector(v_in, img + (j+offset+1)*width, width*sizeof(unsigned short));
            } else {
                vbx_dma_to_vector(v_in, sums[i] + (j+offset+1)*width, width*sizeof(unsigned short));
            }
            vbx(VVHU, VADD, v_x[0], v_x[0], v_x[offset]);
            vbx(VVHU, VADD, v_in, v_in, v_in + offset);
            vbx_dma_to_host(sums[i+1] + j*width, v_x[0], width*sizeof(vbx_uhalf_t));

            v_tmp = v_out;
            v_out = v_x[0];
            for(k=0; k < offset; k++) {
                v_x[k] = v_x[k+1];
            }
            v_x[offset] = v_in;
            v_in = v_out;
        }
    }
    vbx_sync();
    vbx_sp_free();
    vbx_shared_free(v_x);

    return sums;
}

unsigned char** LBPRestrictedPatterns2(unsigned short* img, unsigned short **sums, const unsigned width, const unsigned height, const unsigned log)
{
	int l, j, n, cell;

    int coeff_x[8] = {0, 1, 2, 2, 2, 1, 0, 0};
    int coeff_y[8] = {0, 0, 0, 1, 2, 2, 2, 1};

    unsigned char **patterns = (unsigned char**)malloc((log+1)*sizeof(unsigned char*));

    for (l=0; l<log+1; l++) {
        patterns[l] = (unsigned char*)vbx_shared_malloc(height*width*sizeof(unsigned char));
    }

    vbx_ubyte_t *v_cmp =   (vbx_ubyte_t*)vbx_sp_malloc(width*sizeof(vbx_ubyte_t));
    vbx_ubyte_t *v_shift = (vbx_ubyte_t*)vbx_sp_malloc(width*sizeof(vbx_ubyte_t));
    vbx_ubyte_t *v_lbp =   (vbx_ubyte_t*)vbx_sp_malloc(width*sizeof(vbx_ubyte_t));
    vbx_uhalf_t *v_row =   (vbx_uhalf_t*)vbx_sp_malloc(12*width*sizeof(vbx_uhalf_t));

    for(l = 0; l < log+1; l++) {
        cell = 1 << l;
        for (j = 0; j < height - (3*cell-1); j++) {
            if(l == 0) {
                vbx_dma_to_vector(v_row, img+j*width, (3*cell)*width*sizeof(vbx_uhalf_t));
            } else {
                vbx_dma_to_vector(v_row, sums[l]+j*width, (3*cell)*width*sizeof(vbx_uhalf_t));
            }
            vbx_set_vl(width - (3*cell-1));

            vbx(SVBU, VMOV, v_lbp, 0, 0);

            for (n = 0; n < 8; n++) {
                vbx(SVBU, VMOV, v_shift, 0, 0);
                vbx(VVHBU, VSUB, v_cmp, v_row + coeff_x[n]*cell + coeff_y[n]*cell*width, v_row + width*cell + cell);
                vbx(SVBU, VCMV_GEZ, v_shift, 1<<(7-n), v_cmp);
                vbx(VVBU, VADD, v_lbp, v_lbp, v_shift);
            }
            vbx_dma_to_host(patterns[l]+j*width, v_lbp, width*sizeof(vbx_ubyte_t));
        }
    }
    vbx_sync();
    vbx_sp_free();

    return patterns;
}

unsigned char** LBPRestricted(unsigned short *img, const unsigned width, const unsigned height, const unsigned log)
{
    int coeff_x[8] = {0, 1, 2, 2, 2, 1, 0, 0};
    int coeff_y[8] = {0, 0, 0, 1, 2, 2, 2, 1};

    int i, j, l, n, x, y;

    unsigned char **patterns = (unsigned char**)malloc((log+1)*sizeof(unsigned char*));
    for (l=0; l<log+1; l++) {
        patterns[l] = (unsigned char*)vbx_shared_malloc(height*width*sizeof(unsigned char));
    }


    int num_in, num_dbl, num_row;
    num_in = 4;

    vbx_uhalf_t ***v_dbl, ***v_row, **v_in, *v_tmp;
    vbx_ubyte_t *v_lbp, *v_cmp, *v_shift;

    v_in = (vbx_uhalf_t**)malloc(num_in*sizeof(vbx_uhalf_t*));
    v_dbl = (vbx_uhalf_t***)malloc(log*sizeof(vbx_uhalf_t**));
    v_row = (vbx_uhalf_t***)malloc(log*sizeof(vbx_uhalf_t**));

    v_in[0] = (vbx_uhalf_t*)vbx_sp_malloc(num_in*width*sizeof(vbx_uhalf_t));
    for(j=1; j< num_in; j++){
        v_in[j] = v_in[0] + j*width;
    }

    for(l=0; l<log; l++){
        num_dbl = 3 + l; //2x2 needs 0 1, 2   4x4 needs 0 1 2, 3
        v_dbl[l] = (vbx_uhalf_t**)malloc(num_dbl*sizeof(vbx_uhalf_t*));
        v_dbl[l][0] = (vbx_uhalf_t*)vbx_sp_malloc(num_dbl*width*sizeof(vbx_uhalf_t));
        for(j=1; j<num_dbl; j++){
            v_dbl[l][j] = v_dbl[l][0] + j*width;
        }
        num_row = 6 + 4*l; //2x2 needs 6   4x4 needs 10
        v_row[l] = (vbx_uhalf_t**)malloc(num_row*sizeof(vbx_uhalf_t*));
        v_row[l][0] = (vbx_uhalf_t*)vbx_sp_malloc(num_row*width*sizeof(vbx_uhalf_t));
        for(j=1; j<num_row; j++){
            v_row[l][j] = v_row[l][0] + j*width;
        }
    }

    v_lbp = (vbx_ubyte_t*)vbx_sp_malloc((log+1)*width*sizeof(vbx_ubyte_t));
    v_cmp = (vbx_ubyte_t*)vbx_sp_malloc((log+1)*width*sizeof(vbx_ubyte_t));
    v_shift = (vbx_ubyte_t*)vbx_sp_malloc((log+1)*width*sizeof(vbx_ubyte_t));

    /* setup initial rows */
    vbx_dma_to_vector(v_in[0], img, 3*width*sizeof(vbx_uhalf_t));

    vbx_set_vl(width);
    vbx_set_2D(2, width*sizeof(vbx_ubyte_t), width*sizeof(vbx_ubyte_t), width*sizeof(vbx_ubyte_t));

    for(j = 0; j < 16; j++) {
        if (j >= 0){
            if (j < 2){
                vbx(VVHU, VADD, v_dbl[0][j], v_in[j], v_in[j] + 1);
            } else {
                vbx(VVHU, VADD, v_dbl[0][2], v_in[2], v_in[2] + 1);
            }
        }

        if (j >= 2) {
            vbx_dma_to_vector(v_in[3], img+(j+1)*width, width*sizeof(vbx_uhalf_t));

            if(j < 7){
                vbx(SVBU, VMOV, v_lbp, 0, 0);
                for (n = 0; n < 8; n++) {
                    y = coeff_y[n];
                    x = coeff_x[n];
                    vbx(SVBU, VMOV, v_shift, 0, 0);
                    vbx(VVHBU, VSUB, v_cmp, v_in[y] + x, v_in[1] + 1);
                    vbx(SVBU, VCMV_GEZ, v_shift, 1<<(7-n), v_cmp);
                    vbx(VVBU, VADD, v_lbp, v_lbp, v_shift);
                }
            }

            if (j < 7) {
                vbx(VVHU, VADD, v_row[0][j-2], v_dbl[0][0], v_dbl[0][0+1]);
            } else {
                vbx(VVHU, VADD, v_row[0][5], v_dbl[0][0], v_dbl[0][0+1]);
            }
        }

        if (j >= 4) {
            if (j < 7) {
                vbx(VVHU, VADD, v_dbl[1][j-4], v_row[0][j-4], v_row[0][j-4] + 2);
            } else {
                vbx(VVHU, VADD, v_dbl[1][3], v_row[0][3], v_row[0][3] + 2);
            }
        }

        if (j >= 7) {
            vbx(VVHU, VADD, v_row[1][j-7], v_dbl[1][0], v_dbl[1][0+2]);

            vbx_2D(SVBU, VMOV, v_lbp, 0, 0);
            for (n = 0; n < 8; n++) {
                y = coeff_y[n];
                x = coeff_x[n];
                vbx_2D(SVBU, VMOV, v_shift, 0, 0);
                vbx(VVHBU, VSUB, v_cmp,       v_in[y] + x, v_in[1] + 1);
                vbx(VVHBU, VSUB, v_cmp+width, v_row[0][y*2] + x*2, v_row[0][2] + 2);
                vbx_2D(SVBU, VCMV_GEZ, v_shift, 1<<(7-n), v_cmp);
                vbx_2D(VVBU, VADD, v_lbp, v_lbp, v_shift);
            }
        }
        /* pointer swaps and dma out results*/
        if (j >= 2){
            vbx_dma_to_host(patterns[0]+(j-2)*width, v_lbp, width*sizeof(vbx_ubyte_t));

            v_tmp = v_in[0];
            v_in[0] = v_in[1];
            v_in[1] = v_in[2];
            v_in[2] = v_in[3];
            v_in[3] = v_tmp;

            v_tmp = v_dbl[0][0];
            v_dbl[0][0] = v_dbl[0][1];
            v_dbl[0][1] = v_dbl[0][2];
            v_dbl[0][2] = v_tmp;
        }
        if (j >= 7){
            vbx_dma_to_host(patterns[1]+(j-7)*width, v_lbp+width, width*sizeof(vbx_ubyte_t));

            v_tmp = v_row[0][0];
            for(i = 0; i < 5; i++){
                v_row[0][i] = v_row[0][i+1];
            }
            v_row[0][5] = v_tmp;

            v_tmp = v_dbl[1][0];
            v_dbl[1][0] = v_dbl[1][1];
            v_dbl[1][1] = v_dbl[1][2];
            v_dbl[1][2] = v_dbl[1][3];
            v_dbl[1][3] = v_tmp;
        }
    }

    vbx_set_2D(3, width*sizeof(vbx_ubyte_t), width*sizeof(vbx_ubyte_t), width*sizeof(vbx_ubyte_t));

    for (j = 16; j < height + 5; j++) {
        if(j < height - 1) {
        vbx_dma_to_vector(v_in[3], img+(j+1)*width, width*sizeof(vbx_uhalf_t));
        }

        vbx_2D(SVBU, VMOV, v_lbp, 0, 0); 
        for (n = 0; n < 8; n++) {
            y = coeff_y[n];
            x = coeff_x[n];
            vbx_2D(SVBU, VMOV, v_shift, 0, 0); 
            vbx(VVHBU, VSUB, v_cmp+2*width,       v_in[y] + x, v_in[1] + 1); 
            vbx(VVHBU, VSUB, v_cmp+width, v_row[0][y*2] + x*2, v_row[0][2] + 2); 
            vbx(VVHBU, VSUB, v_cmp, v_row[1][y*4] + x*4, v_row[1][4] + 4); 
            vbx_2D(SVBU, VCMV_GEZ, v_shift, 1<<(7-n), v_cmp); 
            vbx_2D(VVBU, VADD, v_lbp, v_lbp, v_shift); 
        }

        vbx(VVHU, VADD, v_dbl[0][2], v_in[2], v_in[2] + 1);
        vbx(VVHU, VADD, v_row[0][5], v_dbl[0][0], v_dbl[0][0+1]);
        vbx(VVHU, VADD, v_dbl[1][3], v_row[0][3], v_row[0][3] + 2);
        vbx(VVHU, VADD, v_row[1][9], v_dbl[1][0], v_dbl[1][0+2]);

        if(j < height){
            vbx_dma_to_host(patterns[0]+(j-2)*width, v_lbp+2*width, width*sizeof(vbx_ubyte_t));
        }
        if(j < height + 2){
            vbx_dma_to_host(patterns[1]+(j-7)*width, v_lbp+width, width*sizeof(vbx_ubyte_t));
        }
        vbx_dma_to_host(patterns[2]+(j-16)*width, v_lbp, width*sizeof(vbx_ubyte_t));

        /* swap pointers */
        v_tmp = v_in[0];
        v_in[0] = v_in[1];
        v_in[1] = v_in[2];
        v_in[2] = v_in[3];
        v_in[3] = v_tmp;

        v_tmp = v_dbl[0][0];
        v_dbl[0][0] = v_dbl[0][1];
        v_dbl[0][1] = v_dbl[0][2];
        v_dbl[0][2] = v_tmp;

        v_tmp = v_row[0][0];
        for(i = 0; i < 5; i++){
            v_row[0][i] = v_row[0][i+1];
        }
        v_row[0][5] = v_tmp;

        v_tmp = v_row[1][0];
        for(i = 0; i < 9; i++){
            v_row[1][i] = v_row[1][i+1];
        }
        v_row[1][9] = v_tmp;

        v_tmp = v_dbl[1][0];
        v_dbl[1][0] = v_dbl[1][1];
        v_dbl[1][1] = v_dbl[1][2];
        v_dbl[1][2] = v_dbl[1][3];
        v_dbl[1][3] = v_tmp;
    }

    vbx_sync();
    vbx_sp_free();
    for(l=0; l<log; l++){
        free(v_dbl[l]);
        free(v_row[l]);
    }
    free(v_in);
    free(v_dbl);
    free(v_row);

    return patterns;
}

#if 0
unsigned char** LBPRestricted8H(unsigned char *img, const unsigned width, const unsigned height, const unsigned log)
{
    int coeff_x[8] = {0, 1, 2, 2, 2, 1, 0, 0};
    int coeff_y[8] = {0, 0, 0, 1, 2, 2, 2, 1};

    int i, j, l, n, x, y;

    unsigned char **patterns = (unsigned char**)malloc((log+1)*sizeof(unsigned char*));
    for (l=0; l<log+1; l++) {
        patterns[l] = (unsigned char*)vbx_shared_malloc(height*width*sizeof(unsigned char));
    }


    int num_in, num_dbl, num_row;
    num_in = 4;

    vbx_uhalf_t ***v_dbl, ***v_row, **v_in, *v_tmp;
    vbx_ubyte_t *v_lbp, *v_cmp, *v_shift;

    v_in = (vbx_uhalf_t**)malloc(num_in*sizeof(vbx_uhalf_t*));
    v_dbl = (vbx_uhalf_t***)malloc(log*sizeof(vbx_uhalf_t**));
    v_row = (vbx_uhalf_t***)malloc(log*sizeof(vbx_uhalf_t**));

    v_in[0] = (vbx_uhalf_t*)vbx_sp_malloc(num_in*width*sizeof(vbx_uhalf_t));
    for(j=1; j< num_in; j++){
        v_in[j] = v_in[0] + j*width;
    }

    for(l=0; l<log; l++){
        num_dbl = 3 + l; //2x2 needs 0 1, 2   4x4 needs 0 1 2, 3
        v_dbl[l] = (vbx_uhalf_t**)malloc(num_dbl*sizeof(vbx_uhalf_t*));
        v_dbl[l][0] = (vbx_uhalf_t*)vbx_sp_malloc(num_dbl*width*sizeof(vbx_uhalf_t));
        for(j=1; j<num_dbl; j++){
            v_dbl[l][j] = v_dbl[l][0] + j*width;
        }
        num_row = 6 + 4*l; //2x2 needs 6   4x4 needs 10
        v_row[l] = (vbx_uhalf_t**)malloc(num_row*sizeof(vbx_uhalf_t*));
        v_row[l][0] = (vbx_uhalf_t*)vbx_sp_malloc(num_row*width*sizeof(vbx_uhalf_t));
        for(j=1; j<num_row; j++){
            v_row[l][j] = v_row[l][0] + j*width;
        }
    }

    v_lbp = (vbx_ubyte_t*)vbx_sp_malloc((log+1)*width*sizeof(vbx_ubyte_t));
    v_cmp = (vbx_ubyte_t*)vbx_sp_malloc((log+1)*width*sizeof(vbx_ubyte_t));
    v_shift = (vbx_ubyte_t*)vbx_sp_malloc((log+1)*width*sizeof(vbx_ubyte_t));

    /* setup initial rows */
    vbx_dma_to_vector(v_in[0], img, 3*width*sizeof(vbx_uhalf_t));

    vbx_set_vl(width);
    vbx_set_2D(2, width*sizeof(vbx_ubyte_t), width*sizeof(vbx_ubyte_t), width*sizeof(vbx_ubyte_t));

    for(j = 0; j < 16; j++) {
        if (j >= 0){
            if (j < 2){
                vbx(VVHU, VADD, v_dbl[0][j], v_in[j], v_in[j] + 1);
            } else {
                vbx(VVHU, VADD, v_dbl[0][2], v_in[2], v_in[2] + 1);
            }
        }

        if (j >= 2) {
            vbx_dma_to_vector(v_in[3], img+(j+1)*width, width*sizeof(vbx_uhalf_t));

            if(j < 7){
                vbx(SVBU, VMOV, v_lbp, 0, 0); 
                for (n = 0; n < 8; n++) {
                    y = coeff_y[n];
                    x = coeff_x[n];
                    vbx(SVBU, VMOV, v_shift, 0, 0);
                    vbx(VVHBU, VSUB, v_cmp, v_in[y] + x, v_in[1] + 1);
                    vbx(SVBU, VCMV_GEZ, v_shift, 1<<(7-n), v_cmp);
                    vbx(VVBU, VADD, v_lbp, v_lbp, v_shift);
                }
            }

            if (j < 7) {
                vbx(VVHU, VADD, v_row[0][j-2], v_dbl[0][0], v_dbl[0][0+1]);
            } else {
                vbx(VVHU, VADD, v_row[0][5], v_dbl[0][0], v_dbl[0][0+1]);
            }
        }

        if (j >= 4) {
            if (j < 7) {
                vbx(VVHU, VADD, v_dbl[1][j-4], v_row[0][j-4], v_row[0][j-4] + 2);
            } else {
                vbx(VVHU, VADD, v_dbl[1][3], v_row[0][3], v_row[0][3] + 2);
            }
        }

        if (j >= 7) {
            vbx(VVHU, VADD, v_row[1][j-7], v_dbl[1][0], v_dbl[1][0+2]);

            vbx_2D(SVBU, VMOV, v_lbp, 0, 0);
            for (n = 0; n < 8; n++) {
                y = coeff_y[n];
                x = coeff_x[n];
                vbx_2D(SVBU, VMOV, v_shift, 0, 0);
                vbx(VVHBU, VSUB, v_cmp,       v_in[y] + x, v_in[1] + 1);
                vbx(VVHBU, VSUB, v_cmp+width, v_row[0][y*2] + x*2, v_row[0][2] + 2);
                vbx_2D(SVBU, VCMV_GEZ, v_shift, 1<<(7-n), v_cmp);
                vbx_2D(VVBU, VADD, v_lbp, v_lbp, v_shift);
            }
        }
        /* pointer swaps and dma out results*/
        if (j >= 2){
            vbx_dma_to_host(patterns[0]+(j-2)*width, v_lbp, width*sizeof(vbx_ubyte_t));

            v_tmp = v_in[0];
            v_in[0] = v_in[1];
            v_in[1] = v_in[2];
            v_in[2] = v_in[3];
            v_in[3] = v_tmp;

            v_tmp = v_dbl[0][0];
            v_dbl[0][0] = v_dbl[0][1];
            v_dbl[0][1] = v_dbl[0][2];
            v_dbl[0][2] = v_tmp;
        }
        if (j >= 7){
            vbx_dma_to_host(patterns[1]+(j-7)*width, v_lbp+width, width*sizeof(vbx_ubyte_t));

            v_tmp = v_row[0][0];
            for(i = 0; i < 5; i++){
                v_row[0][i] = v_row[0][i+1];
            }
            v_row[0][5] = v_tmp;

            v_tmp = v_dbl[1][0];
            v_dbl[1][0] = v_dbl[1][1];
            v_dbl[1][1] = v_dbl[1][2];
            v_dbl[1][2] = v_dbl[1][3];
            v_dbl[1][3] = v_tmp;
        }
    }

    vbx_set_2D(3, width*sizeof(vbx_ubyte_t), width*sizeof(vbx_ubyte_t), width*sizeof(vbx_ubyte_t));

    for (j = 16; j < height + 5; j++) {
        if(j < height - 1) {
        vbx_dma_to_vector(v_in[3], img+(j+1)*width, width*sizeof(vbx_uhalf_t));
        }

        vbx_2D(SVBU, VMOV, v_lbp, 0, 0); 
        for (n = 0; n < 8; n++) {
            y = coeff_y[n];
            x = coeff_x[n];
            vbx_2D(SVBU, VMOV, v_shift, 0, 0); 
            vbx(VVHBU, VSUB, v_cmp+2*width,       v_in[y] + x, v_in[1] + 1); 
            vbx(VVHBU, VSUB, v_cmp+width, v_row[0][y*2] + x*2, v_row[0][2] + 2); 
            vbx(VVHBU, VSUB, v_cmp, v_row[1][y*4] + x*4, v_row[1][4] + 4); 
            vbx_2D(SVBU, VCMV_GEZ, v_shift, 1<<(7-n), v_cmp); 
            vbx_2D(VVBU, VADD, v_lbp, v_lbp, v_shift); 
        }

        vbx(VVHU, VADD, v_dbl[0][2], v_in[2], v_in[2] + 1);
        vbx(VVHU, VADD, v_row[0][5], v_dbl[0][0], v_dbl[0][0+1]);
        vbx(VVHU, VADD, v_dbl[1][3], v_row[0][3], v_row[0][3] + 2);
        vbx(VVHU, VADD, v_row[1][9], v_dbl[1][0], v_dbl[1][0+2]);

        if(j < height){
            vbx_dma_to_host(patterns[0]+(j-2)*width, v_lbp+2*width, width*sizeof(vbx_ubyte_t));
        }
        if(j < height + 2){
            vbx_dma_to_host(patterns[1]+(j-7)*width, v_lbp+width, width*sizeof(vbx_ubyte_t));
        }
        vbx_dma_to_host(patterns[2]+(j-16)*width, v_lbp, width*sizeof(vbx_ubyte_t));

        /* swap pointers */
        v_tmp = v_in[0];
        v_in[0] = v_in[1];
        v_in[1] = v_in[2];
        v_in[2] = v_in[3];
        v_in[3] = v_tmp;

        v_tmp = v_dbl[0][0];
        v_dbl[0][0] = v_dbl[0][1];
        v_dbl[0][1] = v_dbl[0][2];
        v_dbl[0][2] = v_tmp;

        v_tmp = v_row[0][0];
        for(i = 0; i < 5; i++){
            v_row[0][i] = v_row[0][i+1];
        }
        v_row[0][5] = v_tmp;

        v_tmp = v_row[1][0];
        for(i = 0; i < 9; i++){
            v_row[1][i] = v_row[1][i+1];
        }
        v_row[1][9] = v_tmp;

        v_tmp = v_dbl[1][0];
        v_dbl[1][0] = v_dbl[1][1];
        v_dbl[1][1] = v_dbl[1][2];
        v_dbl[1][2] = v_dbl[1][3];
        v_dbl[1][3] = v_tmp;
    }

    vbx_sync();
    vbx_sp_free();
    for(l=0; l<log; l++){
        free(v_dbl[l]);
        free(v_row[l]);
    }
    free(v_in);
    free(v_dbl);
    free(v_row);

    return patterns;
}
#endif

unsigned char** LBPRestricted8(unsigned char *img, const unsigned width, const unsigned height, const unsigned log)
{
    int coeff_x[8] = {0, 1, 2, 2, 2, 1, 0, 0};
    int coeff_y[8] = {0, 0, 0, 1, 2, 2, 2, 1};

    int i, j, l, n, x, y;

    unsigned char **patterns = (unsigned char**)malloc((log+1)*sizeof(unsigned char*));
    for (l=0; l<log+1; l++) {
        patterns[l] = (unsigned char*)vbx_shared_malloc(height*width*sizeof(unsigned char));
    }


    int num_in, num_dbl, num_row;
    num_in = 4;

    vbx_uhalf_t ***v_dbl, ***v_row, *v_tmp;
    vbx_ubyte_t **v_in, *v_lbp, *v_cmp, *v_shift, *v_tmp8;

    v_in = (vbx_ubyte_t**)malloc(num_in*sizeof(vbx_ubyte_t*));
    v_dbl = (vbx_uhalf_t***)malloc(log*sizeof(vbx_uhalf_t**));
    v_row = (vbx_uhalf_t***)malloc(log*sizeof(vbx_uhalf_t**));

    /* v_in[0] = (vbx_ubyte_t*)vbx_sp_malloc(num_in*width*sizeof(vbx_ubyte_t)); */
    for(j=0; j< num_in; j++){
        v_in[j] = (vbx_ubyte_t*)vbx_sp_malloc(width*sizeof(vbx_ubyte_t));
        /* v_in[j] = v_in[0] + j*width; */
    }

    for(l=0; l<log; l++){
        num_dbl = 3 + l; //2x2 needs 0 1, 2   4x4 needs 0 1 2, 3
        v_dbl[l] = (vbx_uhalf_t**)malloc(num_dbl*sizeof(vbx_uhalf_t*));
        v_dbl[l][0] = (vbx_uhalf_t*)vbx_sp_malloc(num_dbl*width*sizeof(vbx_uhalf_t));
        for(j=1; j<num_dbl; j++){
            v_dbl[l][j] = v_dbl[l][0] + j*width;
        }
        num_row = 6 + 4*l; //2x2 needs 6   4x4 needs 10
        v_row[l] = (vbx_uhalf_t**)malloc(num_row*sizeof(vbx_uhalf_t*));
        v_row[l][0] = (vbx_uhalf_t*)vbx_sp_malloc(num_row*width*sizeof(vbx_uhalf_t));
        for(j=1; j<num_row; j++){
            v_row[l][j] = v_row[l][0] + j*width;
        }
    }

    v_lbp = (vbx_ubyte_t*)vbx_sp_malloc((log+1)*width*sizeof(vbx_ubyte_t));
    v_cmp = (vbx_ubyte_t*)vbx_sp_malloc((log+1)*width*sizeof(vbx_ubyte_t));
    v_shift = (vbx_ubyte_t*)vbx_sp_malloc((log+1)*width*sizeof(vbx_ubyte_t));

    /* setup initial rows */
    vbx_dma_to_vector_2D(v_in[0], img, width, 3, v_in[1]-v_in[0], width);

    vbx_set_vl(width);
    vbx_set_2D(2, width*sizeof(vbx_ubyte_t), width*sizeof(vbx_ubyte_t), width*sizeof(vbx_ubyte_t));

    for(j = 0; j < 16; j++) {
        if (j >= 0){
            if (j < 2){
                vbx(VVBHU, VADD, v_dbl[0][j], v_in[j], v_in[j] + 1);
            } else {
                vbx(VVBHU, VADD, v_dbl[0][2], v_in[2], v_in[2] + 1);
            }
        }

        if (j >= 2) {
            vbx_dma_to_vector(v_in[3], img+(j+1)*width, width*sizeof(vbx_ubyte_t));

            if(j < 7){
                vbx(SVBU, VMOV, v_lbp, 0, 0); 
                for (n = 0; n < 8; n++) {
                    y = coeff_y[n];
                    x = coeff_x[n];
                    vbx(SVBU, VMOV, v_shift, 0, 0); 
                    vbx(VVBU, VSUB, v_cmp, v_in[y] + x, v_in[1] + 1); 
                    vbx(SVBU, VCMV_GEZ, v_shift, 1<<(7-n), v_cmp); 
                    vbx(VVBU, VADD, v_lbp, v_lbp, v_shift); 
                }
            }

            if (j < 7) {
                vbx(VVHU, VADD, v_row[0][j-2], v_dbl[0][0], v_dbl[0][0+1]);
            } else {
                vbx(VVHU, VADD, v_row[0][5], v_dbl[0][0], v_dbl[0][0+1]);
            }
        }

        if (j >= 4) {
            if (j < 7) {
                vbx(VVHU, VADD, v_dbl[1][j-4], v_row[0][j-4], v_row[0][j-4] + 2);
            } else {
                vbx(VVHU, VADD, v_dbl[1][3], v_row[0][3], v_row[0][3] + 2);
            }
        }

        if (j >= 7) {
            vbx(VVHU, VADD, v_row[1][j-7], v_dbl[1][0], v_dbl[1][0+2]);

            vbx_2D(SVBU, VMOV, v_lbp, 0, 0); 
            for (n = 0; n < 8; n++) {
                y = coeff_y[n];
                x = coeff_x[n];
                vbx_2D(SVBU, VMOV, v_shift, 0, 0); 
                vbx(VVBU, VSUB, v_cmp,       v_in[y] + x, v_in[1] + 1); 
                vbx(VVHBU, VSUB, v_cmp+width, v_row[0][y*2] + x*2, v_row[0][2] + 2); 
                vbx_2D(SVBU, VCMV_GEZ, v_shift, 1<<(7-n), v_cmp); 
                vbx_2D(VVBU, VADD, v_lbp, v_lbp, v_shift); 
            }
        }
        /* pointer swaps and dma out results*/
        if (j >= 2){
            vbx_dma_to_host(patterns[0]+(j-2)*width, v_lbp, width*sizeof(vbx_ubyte_t));

            v_tmp8 = v_in[0];
            v_in[0] = v_in[1];
            v_in[1] = v_in[2];
            v_in[2] = v_in[3];
            v_in[3] = v_tmp8;

            v_tmp = v_dbl[0][0];
            v_dbl[0][0] = v_dbl[0][1];
            v_dbl[0][1] = v_dbl[0][2];
            v_dbl[0][2] = v_tmp;
        }
        if (j >= 7){
            vbx_dma_to_host(patterns[1]+(j-7)*width, v_lbp+width, width*sizeof(vbx_ubyte_t));

            v_tmp = v_row[0][0];
            for(i = 0; i < 5; i++){
                v_row[0][i] = v_row[0][i+1];
            }
            v_row[0][5] = v_tmp;

            v_tmp = v_dbl[1][0];
            v_dbl[1][0] = v_dbl[1][1];
            v_dbl[1][1] = v_dbl[1][2];
            v_dbl[1][2] = v_dbl[1][3];
            v_dbl[1][3] = v_tmp;
        }
    }

    vbx_set_2D(3, width*sizeof(vbx_ubyte_t), width*sizeof(vbx_ubyte_t), width*sizeof(vbx_ubyte_t));

    for (j = 16; j < height + 5; j++) {
        if(j < height - 1) {
        vbx_dma_to_vector(v_in[3], img+(j+1)*width, width*sizeof(vbx_ubyte_t));
        }

        vbx_2D(SVBU, VMOV, v_lbp, 0, 0); 
        for (n = 0; n < 8; n++) {
            y = coeff_y[n];
            x = coeff_x[n];
            vbx_2D(SVBU, VMOV, v_shift, 0, 0); 
            vbx(VVBU, VSUB, v_cmp+2*width,       v_in[y] + x, v_in[1] + 1); 
            vbx(VVHBU, VSUB, v_cmp+width, v_row[0][y*2] + x*2, v_row[0][2] + 2); 
            vbx(VVHBU, VSUB, v_cmp, v_row[1][y*4] + x*4, v_row[1][4] + 4); 
            vbx_2D(SVBU, VCMV_GEZ, v_shift, 1<<(7-n), v_cmp); 
            vbx_2D(VVBU, VADD, v_lbp, v_lbp, v_shift); 
        }

        vbx(VVBHU, VADD, v_dbl[0][2], v_in[2], v_in[2] + 1);
        vbx(VVHU, VADD, v_row[0][5], v_dbl[0][0], v_dbl[0][0+1]);
        vbx(VVHU, VADD, v_dbl[1][3], v_row[0][3], v_row[0][3] + 2);
        vbx(VVHU, VADD, v_row[1][9], v_dbl[1][0], v_dbl[1][0+2]);

        if(j < height){
            vbx_dma_to_host(patterns[0]+(j-2)*width, v_lbp+2*width, width*sizeof(vbx_ubyte_t));
        }
        if(j < height + 2){
            vbx_dma_to_host(patterns[1]+(j-7)*width, v_lbp+width, width*sizeof(vbx_ubyte_t));
        }
        vbx_dma_to_host(patterns[2]+(j-16)*width, v_lbp, width*sizeof(vbx_ubyte_t));

        /* swap pointers */
        v_tmp8 = v_in[0];
        v_in[0] = v_in[1];
        v_in[1] = v_in[2];
        v_in[2] = v_in[3];
        v_in[3] = v_tmp8;

        v_tmp = v_dbl[0][0];
        v_dbl[0][0] = v_dbl[0][1];
        v_dbl[0][1] = v_dbl[0][2];
        v_dbl[0][2] = v_tmp;

        v_tmp = v_row[0][0];
        for(i = 0; i < 5; i++){
            v_row[0][i] = v_row[0][i+1];
        }
        v_row[0][5] = v_tmp;

        v_tmp = v_row[1][0];
        for(i = 0; i < 9; i++){
            v_row[1][i] = v_row[1][i+1];
        }
        v_row[1][9] = v_tmp;

        v_tmp = v_dbl[1][0];
        v_dbl[1][0] = v_dbl[1][1];
        v_dbl[1][1] = v_dbl[1][2];
        v_dbl[1][2] = v_dbl[1][3];
        v_dbl[1][3] = v_tmp;
    }

    vbx_sync();
    vbx_sp_free();
    for(l=0; l<log; l++){
        free(v_dbl[l]);
        free(v_row[l]);
    }
    free(v_in);
    free(v_dbl);
    free(v_row);

    return patterns;
}

unsigned char** LBPRestrictedCI3(unsigned char *img, const unsigned width, const unsigned height, const unsigned log)
{
    int i, l, y;
	vbx_mxp_t *this_mxp = VBX_GET_THIS_MXP();
	int vci_lanes = this_mxp->vcustom1_lanes;

    unsigned char **patterns = (unsigned char**)malloc((log+1)*sizeof(unsigned char*));
    for (l=0; l<log+1; l++) {
        patterns[l] = (unsigned char*)vbx_shared_malloc(height*width*sizeof(unsigned char));
    }

#if 1
    int max_rows = this_mxp->scratchpad_size / (width*2*sizeof(vbx_ubyte_t));
    vbx_ubyte_t* v_in =  (vbx_ubyte_t*)vbx_sp_malloc(max_rows*width*sizeof(vbx_ubyte_t));
    vbx_ubyte_t* v_out =  (vbx_ubyte_t*)vbx_sp_malloc(max_rows*width*sizeof(vbx_ubyte_t));
    if (v_out == NULL) {
        printf("not enough scratch\n");
    }

    vbx_set_vl(vci_lanes*4);

    for (y = 0; y < height; y = y + (max_rows-12)) {
        int rows = max_rows;
        if (rows > height - y) {
            rows = height - y;
        } 
        vbx_dma_to_vector(v_in, img + y * width, rows*width*sizeof(vbx_ubyte_t));
        vbx_sync();

        vbx_set_2D(rows, width, 0, width);
        for (l=0; l<log+1; l++) {
            for (i = 0; i < width - (vci_lanes*4 - 3*4-1); i += vci_lanes*4 - (3*4-1)) {
                vbx_2D(SVBU, VLBPPAT, v_out+i, l, v_in+i);
            }
            vbx_dma_to_host(patterns[l] + y * width, v_out, rows*width*sizeof(vbx_ubyte_t));
            vbx_sync();
        }
    }
    vbx_sp_free();
#else
    vbx_ubyte_t* v_in =  (vbx_ubyte_t*)vbx_sp_malloc(vci_lanes*4*(height+12)*sizeof(vbx_ubyte_t));
    vbx_ubyte_t* v_out = (vbx_ubyte_t*)vbx_sp_malloc(vci_lanes*4*(height+12)*sizeof(vbx_ubyte_t));
    vbx_ubyte_t* v_temp;
    if (v_out == NULL) {
        printf("out of scratch\n");
    }

    for (i=0; i<width; i += vci_lanes*4-(3*4-1)) {
        vbx_dma_to_vector_2D(v_in, img+i, vci_lanes*4, height, vci_lanes*4, width);
        vbx_sync();
        for (l=0; l<log+1; l++) {

            vbx_set_vl(vci_lanes*4*height);
            vbx(SVBU, VLBPPAT, v_out, l, v_in);

            vbx_dma_to_host_2D(patterns[l]+i, v_out, vci_lanes*4-(3*4-1), height-(2*(1<<l)), width, vci_lanes*4);
            vbx_sync();
        }
    }

    vbx_sp_free();
#endif

    return patterns;
}

unsigned char** LBPRestrictedCI2(unsigned short *img, const unsigned width, const unsigned height, const unsigned log)
{
    int i, j, l, n, x, y;
	vbx_mxp_t *this_mxp = VBX_GET_THIS_MXP();
	int vci_lanes = this_mxp->vcustom1_lanes;
    int cvi_delay = vci_lanes * 3;

    unsigned char **patterns = (unsigned char**)malloc((log+1)*sizeof(unsigned char*));
    for (l=0; l<log+1; l++) {
        patterns[l] = (unsigned char*)vbx_shared_malloc(height*width*sizeof(unsigned char));
    }

    vbx_uhalf_t* v_mem = (vbx_uhalf_t*)vbx_sp_malloc((13+11+9)*width*sizeof(vbx_uhalf_t));
    vbx_uhalf_t** v_1 = (vbx_uhalf_t**)malloc(14*sizeof(vbx_uhalf_t*));
    vbx_uhalf_t** v_2 = (vbx_uhalf_t**)malloc(13*sizeof(vbx_uhalf_t*));
    vbx_uhalf_t** v_4 = (vbx_uhalf_t**)malloc(9*sizeof(vbx_uhalf_t*));

    vbx_ubyte_t* v_1o = (vbx_ubyte_t*)vbx_sp_malloc(width*sizeof(vbx_ubyte_t));
    vbx_ubyte_t* v_2o = (vbx_ubyte_t*)vbx_sp_malloc(width*sizeof(vbx_ubyte_t));
    vbx_ubyte_t* v_4o = (vbx_ubyte_t*)vbx_sp_malloc(width*sizeof(vbx_ubyte_t));

    for (i=0; i<13; i++) {
        v_1[i]= v_mem + (0+i)*width;
    }

    vbx_dma_to_vector(v_1[0], img, 12*width*sizeof(vbx_uhalf_t));

    for (i=0; i<11; i++) {
        v_2[i]= v_mem + (13+i)*width;
    }
    for (i=0; i<9; i++) {
        v_4[i]= v_mem + (13+11+i)*width;
    }

    vbx_uhalf_t* v_tmp;

    vbx_set_vl(width);
    
    for(i=0; i<10; i++) {
        vbx(VVH, VADD, v_2[i], v_1[i], v_1[i+1]);
    }
    for(i=0; i<8; i++) {
        vbx(VVH, VADD, v_4[i], v_2[i], v_2[i+2]);
    }

    vbx_set_vl(width+cvi_delay);
    for(i=0; i<4; i++) {
        vbx(VVHU, VLBPPAT, v_1[i], v_1[i], v_1[i+1]);
        vbx(VVHS, VLBPPAT, v_2[i], v_2[i], v_2[i+2]);
        vbx(VVHS, VLBPPAT, v_4[i], v_4[i], v_4[i+4]);
    }
    for(i=4; i<8; i++) {
        vbx(VVHU, VLBPPAT, v_1[i], v_1[i], v_1[i+1]);
        vbx(VVHS, VLBPPAT, v_2[i], v_2[i], v_2[i+2]);
    }
    for(i=8; i<10; i++) {
        vbx(VVHU, VLBPPAT, v_1[i], v_1[i], v_1[i+1]);
    }
    vbx_set_vl(width);

    for(i=0; i<9; i++) {
        vbx(VVHB, VADD, v_1o, ((vbx_ubyte_t*)v_1[i] + 1), v_1[i+1]);
        vbx_dma_to_host(patterns[0]+i*width, v_1o, width*sizeof(vbx_ubyte_t));
    }
    
    for(i=0; i<6; i++) {
        vbx(VVHB, VADD, v_2o, ((vbx_ubyte_t*)v_2[i] + 1), v_2[i+2]);
        vbx_dma_to_host(patterns[1]+i*width, v_2o, width*sizeof(vbx_ubyte_t));
    }
    
    for (j = 0; j < height - 2 - 9; j++) {
        // sum-custom-interleave
        vbx(VVH, VADD, v_2[10], v_1[10], v_1[11]);
        vbx(VVH, VADD, v_4[8], v_2[8], v_2[10]);

        vbx_set_vl(width+cvi_delay);
        vbx(VVHU, VLBPPAT, v_1[10], v_1[10], v_1[11]);
        vbx(VVHS, VLBPPAT, v_2[8], v_2[8], v_2[10]);
        vbx(VVHS, VLBPPAT, v_4[4], v_4[4], v_4[8]);
        vbx_set_vl(width);

        vbx_dma_to_vector(v_1[12], img+(12+j)*width, width*sizeof(vbx_uhalf_t));

        vbx(VVHB, VADD, v_1o, ((vbx_ubyte_t*)v_1[9] + 1), v_1[10]);
        vbx(VVHB, VADD, v_2o, ((vbx_ubyte_t*)v_2[6] + 1), v_2[8]);
        vbx(VVHB, VADD, v_4o, ((vbx_ubyte_t*)v_4[0] + 1), v_4[4]);

        //swap pointers
        v_tmp = v_1[8];
        for (i=8; i<12; i++){
            v_1[i] = v_1[i+1];
        }
        v_1[12] = v_tmp;

        v_tmp = v_2[6];
        for (i=6; i<10; i++) {
            v_2[i] = v_2[i+1];
        }
        v_2[10] = v_tmp;

        v_tmp = v_4[0];
        for (i=0; i<8; i++) {
            v_4[i] = v_4[i+1];
        }
        v_4[8] = v_tmp;

        // dma out results
        vbx_dma_to_host(patterns[0] + (9+j)*width, v_1o, width*sizeof(vbx_byte_t));
        vbx_dma_to_host(patterns[1] + (6+j)*width, v_2o, width*sizeof(vbx_byte_t));
        vbx_dma_to_host(patterns[2] + (0+j)*width, v_4o, width*sizeof(vbx_byte_t));
    };

    vbx_sync();
    vbx_sp_free();
    free(v_1);
    free(v_2);
    free(v_4);

    return patterns;
}

#if 0
unsigned char** LBPRestrictedCI28(unsigned char *img, const unsigned width, const unsigned height, const unsigned log)
{
    int i, j, l, n, x, y;
	vbx_mxp_t *this_mxp = VBX_GET_THIS_MXP();
	int vci_lanes = this_mxp->vcustom1_lanes;
    int cvi_delay = vci_lanes * 3;
    cvi_delay = 0;

    unsigned char **patterns = (unsigned char**)malloc((log+1)*sizeof(unsigned char*));
    for (l=0; l<log+1; l++) {
        patterns[l] = (unsigned char*)vbx_shared_malloc(height*width*sizeof(unsigned char));
    }

    vbx_ubyte_t* v_mem8 = (vbx_ubyte_t*)vbx_sp_malloc(13*width*sizeof(vbx_ubyte_t));
    vbx_ubyte_t** v_in = (vbx_ubyte_t**)malloc(14*sizeof(vbx_ubyte_t*));

    vbx_uhalf_t* v_mem = (vbx_uhalf_t*)vbx_sp_malloc((13+11+9)*width*sizeof(vbx_uhalf_t));
    vbx_uhalf_t** v_1 = (vbx_uhalf_t**)malloc(14*sizeof(vbx_uhalf_t*));
    vbx_uhalf_t** v_2 = (vbx_uhalf_t**)malloc(13*sizeof(vbx_uhalf_t*));
    vbx_uhalf_t** v_4 = (vbx_uhalf_t**)malloc(9*sizeof(vbx_uhalf_t*));

    vbx_ubyte_t* v_1o = (vbx_ubyte_t*)vbx_sp_malloc(width*sizeof(vbx_ubyte_t));
    vbx_ubyte_t* v_2o = (vbx_ubyte_t*)vbx_sp_malloc(width*sizeof(vbx_ubyte_t));
    vbx_ubyte_t* v_4o = (vbx_ubyte_t*)vbx_sp_malloc(width*sizeof(vbx_ubyte_t));

    for (i=0; i<13; i++) {
        v_in[i]= v_mem8 + (0+i)*width;
    }

    for (i=0; i<13; i++) {
        v_1[i]= v_mem + (0+i)*width;
    }

    vbx_dma_to_vector(v_in[0], img, 12*width*sizeof(vbx_ubyte_t));
    /* vbx_dma_to_vector(v_1[0], img, 12*width*sizeof(vbx_uhalf_t)); */

    for (i=0; i<11; i++) {
        v_2[i]= v_mem + (13+i)*width;
    }
    for (i=0; i<9; i++) {
        v_4[i]= v_mem + (13+11+i)*width;
    }

    vbx_uhalf_t* v_tmp;
    vbx_ubyte_t* v_tmp8;

    vbx_set_vl(width);
#if 0
    vbx_set_2D(10, width*sizeof(vbx_uhalf_t), width*sizeof(vbx_ubyte_t), width*sizeof(vbx_ubyte_t));
    vbx_2D(VVBHU, VADD, v_2[i], v_in[i], v_in[i+1]);
#else
    
    for(i=0; i<10; i++) {
        vbx(VVBHU, VADD, v_2[i], v_in[i], v_in[i+1]);
    }
#endif


#if 0
    vbx_set_2D(8, width*sizeof(vbx_uhalf_t), width*sizeof(vbx_uhalf_t), width*sizeof(vbx_uhalf_t));
    vbx_2D(VVHU, VADD, v_4[i], v_2[i], v_2[i+2]);
#else
    for(i=0; i<8; i++) {
        vbx(VVHU, VADD, v_4[i], v_2[i], v_2[i+2]);
    }
#endif

    vbx_set_vl(width+cvi_delay);
#if 0
    vbx_set_2D(4, width*sizeof(vbx_uhalf_t), width*sizeof(vbx_ubyte_t), width*sizeof(vbx_ubyte_t));
    vbx(VVBHU, VLBPPAT, v_1[i], v_in[i], v_in[i+1]);

    vbx_set_2D(4, width*sizeof(vbx_uhalf_t), width*sizeof(vbx_uhalf_t), width*sizeof(vbx_uhalf_t));
    vbx(VVHS,  VLBPPAT, v_2[i], v_2[i], v_2[i+2]);
    vbx(VVHS,  VLBPPAT, v_4[i], v_4[i], v_4[i+4]);
#else
    for(i=0; i<4; i++) {
        vbx(VVBHU, VLBPPAT, v_1[i], v_in[i], v_in[i+1]);
        vbx(VVHS,  VLBPPAT, v_2[i], v_2[i], v_2[i+2]);
        vbx(VVHS,  VLBPPAT, v_4[i], v_4[i], v_4[i+4]);
    }
#endif
    for(i=4; i<8; i++) {
        vbx(VVBHU, VLBPPAT, v_1[i], v_in[i], v_in[i+1]);
        vbx(VVHS,  VLBPPAT, v_2[i], v_2[i], v_2[i+2]);
    }
    for(i=8; i<10; i++) {
        vbx(VVBHU, VLBPPAT, v_1[i], v_in[i], v_in[i+1]);
    }
    vbx_set_vl(width);

    for(i=0; i<9; i++) {
        vbx(VVHBU, VADD, v_1o, ((vbx_ubyte_t*)v_1[i] + 1), v_1[i+1]);
        vbx_dma_to_host(patterns[0]+i*width, v_1o, width*sizeof(vbx_ubyte_t));
    }
    
    for(i=0; i<6; i++) {
        vbx(VVHBU, VADD, v_2o, ((vbx_ubyte_t*)v_2[i] + 1), v_2[i+2]);
        vbx_dma_to_host(patterns[1]+i*width, v_2o, width*sizeof(vbx_ubyte_t));
    }
    
    for (j = 0; j < height - 2 - 9; j++) {
        // sum-custom-interleave
        vbx(VVBHU, VADD, v_2[10], v_in[10], v_in[11]);
        /* vbx(VVH, VADD, v_2[10], v_1[10], v_1[11]); */
        vbx(VVHU, VADD, v_4[8], v_2[8], v_2[10]);

        vbx_set_vl(width+cvi_delay);
        /* vbx(VVHU, VLBPPAT, v_1[10], v_1[10], v_1[11]); */
        vbx(VVBHU, VLBPPAT, v_1[10], v_in[10], v_in[11]);
        vbx(VVHS, VLBPPAT, v_2[8], v_2[8], v_2[10]);
        vbx(VVHS, VLBPPAT, v_4[4], v_4[4], v_4[8]);
        vbx_set_vl(width);

        /* vbx_dma_to_vector(v_1[12], img+(12+j)*width, width*sizeof(vbx_uhalf_t)); */
        vbx_dma_to_vector(v_in[12], img+(12+j)*width, width*sizeof(vbx_ubyte_t));

        vbx(VVHBU, VADD, v_1o, ((vbx_ubyte_t*)v_1[9] + 1), v_1[10]);
        vbx(VVHBU, VADD, v_2o, ((vbx_ubyte_t*)v_2[6] + 1), v_2[8]);
        vbx(VVHBU, VADD, v_4o, ((vbx_ubyte_t*)v_4[0] + 1), v_4[4]);

        //swap pointers
        v_tmp8 = v_in[8];
        for (i=8; i<12; i++){
            v_in[i] = v_in[i+1];
        }
        v_in[12] = v_tmp8;

        v_tmp = v_1[8];
        for (i=8; i<12; i++){
            v_1[i] = v_1[i+1];
        }
        v_1[12] = v_tmp;

        v_tmp = v_2[6];
        for (i=6; i<10; i++) {
            v_2[i] = v_2[i+1];
        }
        v_2[10] = v_tmp;

        v_tmp = v_4[0];
        for (i=0; i<8; i++) {
            v_4[i] = v_4[i+1];
        }
        v_4[8] = v_tmp;

        // dma out results
        vbx_dma_to_host(patterns[0] + (9+j)*width, v_1o, width*sizeof(vbx_byte_t));
        vbx_dma_to_host(patterns[1] + (6+j)*width, v_2o, width*sizeof(vbx_byte_t));
        vbx_dma_to_host(patterns[2] + (0+j)*width, v_4o, width*sizeof(vbx_byte_t));
    };

    vbx_sync();
    vbx_sp_free();
    free(v_in);
    free(v_1);
    free(v_2);
    free(v_4);

    return patterns;
}
#else
unsigned char** LBPRestrictedCI28(unsigned char **patterns, unsigned char *img, const unsigned width, const unsigned height, const unsigned log)
{
    int j, l, n, x, y;
	vbx_mxp_t *this_mxp = VBX_GET_THIS_MXP();
    int vci_lanes;
    if(VLBPPAT == VCUSTOM0) {
        vci_lanes = this_mxp->vcustom0_lanes;
    } else {
        vci_lanes = this_mxp->vcustom1_lanes;
    }

#if 0
    unsigned char **patterns = (unsigned char**)malloc((log+1)*sizeof(unsigned char*));
    for (l=0; l<log+1; l++) {
        patterns[l] = (unsigned char*)vbx_shared_malloc(height*width*sizeof(unsigned char));
    }
#endif

    vbx_ubyte_t* v_memb = (vbx_ubyte_t*)vbx_sp_malloc((13*width)*sizeof(vbx_ubyte_t));
    vbx_ubyte_t* v_in_0 = v_memb +  0 * width;
    vbx_ubyte_t* v_in_1 = v_memb +  1 * width;
    vbx_ubyte_t* v_in_2 = v_memb +  2 * width;
    vbx_ubyte_t* v_in_3 = v_memb +  3 * width;
    vbx_ubyte_t* v_in_4 = v_memb +  4 * width;
    vbx_ubyte_t* v_in_5 = v_memb +  5 * width;
    vbx_ubyte_t* v_in_6 = v_memb +  6 * width;
    vbx_ubyte_t* v_in_7 = v_memb +  7 * width;
    vbx_ubyte_t* v_in_8 = v_memb +  8 * width;
    vbx_ubyte_t* v_in_9 = v_memb +  9 * width;
    vbx_ubyte_t* v_in_10= v_memb + 10 * width;
    vbx_ubyte_t* v_in_11= v_memb + 11 * width;
    vbx_ubyte_t* v_in_12= v_memb + 12 * width;

    vbx_uhalf_t* v_mem_1 = (vbx_uhalf_t*)vbx_sp_malloc((13*width)*sizeof(vbx_uhalf_t));
    vbx_uhalf_t* v_1_0 = v_mem_1 + 0 * width;
    vbx_uhalf_t* v_1_1 = v_mem_1 + 1 * width;
    vbx_uhalf_t* v_1_2 = v_mem_1 + 2 * width;
    vbx_uhalf_t* v_1_3 = v_mem_1 + 3 * width;
    vbx_uhalf_t* v_1_4 = v_mem_1 + 4 * width;
    vbx_uhalf_t* v_1_5 = v_mem_1 + 5 * width;
    vbx_uhalf_t* v_1_6 = v_mem_1 + 6 * width;
    vbx_uhalf_t* v_1_7 = v_mem_1 + 7 * width;
    vbx_uhalf_t* v_1_8 = v_mem_1 + 8 * width;
    vbx_uhalf_t* v_1_9 = v_mem_1 + 9 * width;
    vbx_uhalf_t* v_1_10= v_mem_1 + 10 * width;
    vbx_uhalf_t* v_1_11= v_mem_1 + 11 * width;
    vbx_uhalf_t* v_1_12= v_mem_1 + 12 * width;

    vbx_uhalf_t* v_mem_2 = (vbx_uhalf_t*)vbx_sp_malloc((11*width)*sizeof(vbx_uhalf_t));
    vbx_uhalf_t* v_2_0 = v_mem_2 + 0 * width;
    vbx_uhalf_t* v_2_1 = v_mem_2 + 1 * width;
    vbx_uhalf_t* v_2_2 = v_mem_2 + 2 * width;
    vbx_uhalf_t* v_2_3 = v_mem_2 + 3 * width;
    vbx_uhalf_t* v_2_4 = v_mem_2 + 4 * width;
    vbx_uhalf_t* v_2_5 = v_mem_2 + 5 * width;
    vbx_uhalf_t* v_2_6 = v_mem_2 + 6 * width;
    vbx_uhalf_t* v_2_7 = v_mem_2 + 7 * width;
    vbx_uhalf_t* v_2_8 = v_mem_2 + 8 * width;
    vbx_uhalf_t* v_2_9 = v_mem_2 + 9 * width;
    vbx_uhalf_t* v_2_10= v_mem_2 + 10 * width;

    vbx_uhalf_t* v_mem_4 = (vbx_uhalf_t*)vbx_sp_malloc((9*width)*sizeof(vbx_uhalf_t));
    vbx_uhalf_t* v_4_0 = v_mem_4 + 0 * width;
    vbx_uhalf_t* v_4_1 = v_mem_4 + 1 * width;
    vbx_uhalf_t* v_4_2 = v_mem_4 + 2 * width;
    vbx_uhalf_t* v_4_3 = v_mem_4 + 3 * width;
    vbx_uhalf_t* v_4_4 = v_mem_4 + 4 * width;
    vbx_uhalf_t* v_4_5 = v_mem_4 + 5 * width;
    vbx_uhalf_t* v_4_6 = v_mem_4 + 6 * width;
    vbx_uhalf_t* v_4_7 = v_mem_4 + 7 * width;
    vbx_uhalf_t* v_4_8 = v_mem_4 + 8 * width;

    vbx_ubyte_t* v_1o = (vbx_ubyte_t*)vbx_sp_malloc(9*width*sizeof(vbx_ubyte_t));
    vbx_ubyte_t* v_2o = (vbx_ubyte_t*)vbx_sp_malloc(6*width*sizeof(vbx_ubyte_t));
    vbx_ubyte_t* v_4o = (vbx_ubyte_t*)vbx_sp_malloc(1*width*sizeof(vbx_ubyte_t));
    if (v_4o == NULL) {
        printf("scratchpad full!\n");
    }

    vbx_dma_to_vector(v_memb, img, 12*width*sizeof(vbx_ubyte_t));

    vbx_uhalf_t* v_tmp;
    vbx_ubyte_t* v_tmp8;

    vbx_set_vl(10*width);
    vbx(VVBHU, VADD, v_mem_2, v_memb, v_memb + 1*width);

    vbx_set_vl(8*width);
    vbx(VVHU, VADD, v_mem_4, v_mem_2, v_mem_2 + 2*width);

    vbx_set_vl(10*width);
    vbx(VVBHU, VLBPPAT, v_mem_1, v_memb, v_memb+1*width);

    vbx_set_vl(8*width);
    vbx(VVHS,  VLBPPAT, v_mem_2, v_mem_2, v_mem_2+2*width);

    vbx_set_vl(4*width);
    vbx(VVHS,  VLBPPAT, v_mem_4, v_mem_4, v_mem_4+4*width);

    vbx_set_vl(9*width);
    vbx(VVHBU, VADD, v_1o, ((vbx_ubyte_t*)v_mem_1 + 1), v_mem_1+1*width);
    vbx_dma_to_host(patterns[0], v_1o, 9*width*sizeof(vbx_ubyte_t));

    vbx_set_vl(6*width);
    vbx(VVHBU, VADD, v_2o, ((vbx_ubyte_t*)v_mem_2 + 1), v_mem_2+2*width);
    vbx_dma_to_host(patterns[1], v_2o, 6*width*sizeof(vbx_ubyte_t));

    vbx_set_vl(width);
    for (j = 0; j < height - 2 - 9; j++) {
        // sum-custom-interleave
        vbx(VVBHU, VADD, v_2_10, v_in_10, v_in_11);
        vbx(VVHU, VADD, v_4_8, v_2_8, v_2_10);

        /* vbx_set_vl(width); */
        vbx(VVBHU, VLBPPAT, v_1_10, v_in_10, v_in_11);
        vbx(VVHS,  VLBPPAT, v_2_8, v_2_8, v_2_10);
        vbx(VVHS,  VLBPPAT, v_4_4, v_4_4, v_4_8);

        if (height > 12 + j) {
            vbx_dma_to_vector(v_in_12, img+(12+j)*width, width*sizeof(vbx_ubyte_t));
        }

        /* vbx_set_vl(width); */
        vbx(VVHBU, VADD, v_1o, ((vbx_ubyte_t*)v_1_9 + 1), v_1_10);
        vbx(VVHBU, VADD, v_2o, ((vbx_ubyte_t*)v_2_6 + 1), v_2_8);
        vbx(VVHBU, VADD, v_4o, ((vbx_ubyte_t*)v_4_0 + 1), v_4_4);

        //swap pointers
        v_tmp8  = v_in_8;
        v_in_8  = v_in_9;
        v_in_9  = v_in_10;
        v_in_10 = v_in_11;
        v_in_11 = v_in_12;
        v_in_12 = v_tmp8;

        v_tmp  = v_1_8;
        v_1_8  = v_1_9;
        v_1_9  = v_1_10;
        v_1_10 = v_1_11;
        v_1_11 = v_1_12;
        v_1_12 = v_tmp;

        v_tmp  = v_2_6;
        v_2_6  = v_2_7;
        v_2_7  = v_2_8;
        v_2_8  = v_2_9;
        v_2_9  = v_2_10;
        v_2_10 = v_tmp;

        v_tmp = v_4_0;
        v_4_0 = v_4_1;
        v_4_1 = v_4_2;
        v_4_2 = v_4_3;
        v_4_3 = v_4_4;
        v_4_4 = v_4_5;
        v_4_5 = v_4_6;
        v_4_6 = v_4_7;
        v_4_7 = v_4_8;
        v_4_8 = v_tmp;

        // dma out results
        vbx_dma_to_host(patterns[0] + (9+j)*width, v_1o, width*sizeof(vbx_byte_t));
        vbx_dma_to_host(patterns[1] + (6+j)*width, v_2o, width*sizeof(vbx_byte_t));
        vbx_dma_to_host(patterns[2] + (0+j)*width, v_4o, width*sizeof(vbx_byte_t));
    };

    vbx_sync();
    vbx_sp_free();

    return patterns;
}
#endif

void LBPRestricted_CI_column_8(unsigned char **patterns, unsigned char *img, const unsigned width, const unsigned height, const unsigned log)
{
    int l, w;
    vbx_mxp_t *this_mxp = VBX_GET_THIS_MXP();
    int vci_lanes;
    if (VLBPCELL1 == VCUSTOM0) {
        vci_lanes = this_mxp->vcustom0_lanes;
    } else if (VLBPCELL1 == VCUSTOM3) {
        vci_lanes = this_mxp->vcustom3_lanes;
    } else {
        printf("custom instructions need to be specified correctly for column");
    }

    int v_instr[3];
    v_instr[0] = VLBPCELL1;
    v_instr[1] = VLBPCELL2;
    v_instr[2] = VLBPCELL4;


    vbx_ubyte_t* v_tmp;
    vbx_ubyte_t* v_a    = (vbx_ubyte_t*)vbx_sp_malloc(height*(4*vci_lanes)*sizeof(vbx_ubyte_t));
    vbx_ubyte_t* v_b    = (vbx_ubyte_t*)vbx_sp_malloc(height*(4*vci_lanes)*sizeof(vbx_ubyte_t));
    vbx_ubyte_t* v_next_in = (vbx_ubyte_t*)vbx_sp_malloc(height*(4*vci_lanes)*sizeof(vbx_ubyte_t));

    vbx_ubyte_t* v_out[3];
    for(l = 0; l < log+1; l++){
        v_out[l] = (vbx_ubyte_t*)vbx_sp_malloc(height*(4*vci_lanes)*sizeof(vbx_ubyte_t));
    }

    vbx_dma_to_vector_2D(v_a,    img + 0*(4*vci_lanes), (4*vci_lanes), height, (4*vci_lanes), width);
    vbx_dma_to_vector_2D(v_b,    img + 1*(4*vci_lanes), (4*vci_lanes), height, (4*vci_lanes), width);
    
    int full_wavefronts = width / (4*vci_lanes);
    int last_wavefront = width % (4*vci_lanes);

    vbx_set_vl((4+height)*(4*vci_lanes));
    for(w=0; w < full_wavefronts; w++){
        vbx_dma_to_vector_2D(v_next_in, img + (w+2)*(4*vci_lanes), (4*vci_lanes), height, (4*vci_lanes), width);

        for(l=0; l < log+1; l++) {
            vbx(VVBU, v_instr[l], v_out[l], v_a, v_b);
        }

        for(l=0; l < log+1; l++) {
            vbx_dma_to_host_2D(patterns[l] + w*(4*vci_lanes), v_out[l], (4*vci_lanes), height, width, (4*vci_lanes));
        }

        v_tmp = v_a;
        v_a = v_b;
        v_b = v_next_in;
        v_next_in = v_tmp;

    }

    if (last_wavefront) {
        for(l=0; l < log+1; l++){
            vbx(VVBU, v_instr[l], v_out[l], v_a, v_b);
            vbx_dma_to_host_2D(patterns[l] + full_wavefronts*(4*vci_lanes), v_out[l], last_wavefront, height, width, (4*vci_lanes));
        }
    }

    vbx_sync();
    vbx_sp_free();
}

void LBPRestricted_CI_column_8_scratch(unsigned char **patterns, unsigned char *img, const unsigned width, const unsigned height, const unsigned log)
{
    int l, w;
    vbx_mxp_t *this_mxp = VBX_GET_THIS_MXP();
    int vci_lanes;
    if (VLBPCELL1 == VCUSTOM0) {
        vci_lanes = this_mxp->vcustom0_lanes;
    } else if (VLBPCELL1 == VCUSTOM3) {
        vci_lanes = this_mxp->vcustom3_lanes;
    } else {
        printf("custom instructions need to be specified correctly for column");
    }

    int v_instr[3];
    v_instr[0] = VLBPCELL1;
    v_instr[1] = VLBPCELL2;
    v_instr[2] = VLBPCELL4;

    /* vbx_ubyte_t* v_tmp; */
    vbx_ubyte_t* v_in = (vbx_ubyte_t*)vbx_sp_malloc(height*width*sizeof(vbx_ubyte_t));
    vbx_ubyte_t* v_out[3];

    for(l = 0; l < log+1; l++){
        v_out[l] = (vbx_ubyte_t*)vbx_sp_malloc(height*width*sizeof(vbx_ubyte_t));
    }

    vbx_dma_to_vector(v_in, img, height*width);
    
    int full_wavefronts = width / (4*vci_lanes);
    int last_wavefront = width % (4*vci_lanes);

    vbx_set_vl(4*vci_lanes);
    vbx_set_2D((4+height), width, width, width);

    for(l=0; l < log+1; l++){
        for(w=0; w < full_wavefronts; w++){
            vbx_2D(VVBU, v_instr[l], v_out[l]+w*(4*vci_lanes), v_in+w*(4*vci_lanes), v_in+(w+1)*(4*vci_lanes));
        }
        if (last_wavefront) {
            vbx_set_vl(last_wavefront);
            vbx_2D(VVBU, v_instr[l], v_out[l]+w*(4*vci_lanes), v_in+w*(4*vci_lanes), v_in+(w+1)*(4*vci_lanes));
            vbx_set_vl(4*vci_lanes);
        }

        vbx_dma_to_host(patterns[l], v_out[l], height*width);
    }

    vbx_sync();
    vbx_sp_free();
}

unsigned char** LBPRestrictedCI(unsigned short *img, const unsigned width, const unsigned height, const unsigned log)
{
    int i, j, l, n, x, y;
	vbx_mxp_t *this_mxp = VBX_GET_THIS_MXP();
	int vci_lanes = this_mxp->vcustom1_lanes;

    unsigned char **patterns = (unsigned char**)malloc((log+1)*sizeof(unsigned char*));
    for (l=0; l<log+1; l++) {
        patterns[l] = (unsigned char*)vbx_shared_malloc(height*width*sizeof(unsigned char));
    }

    vbx_uhalf_t* v_mem = (vbx_uhalf_t*)vbx_sp_malloc((13+11+9)*width*sizeof(vbx_uhalf_t));
    vbx_uhalf_t** v_1 = (vbx_uhalf_t**)malloc(14*sizeof(vbx_uhalf_t*));
    vbx_uhalf_t** v_2 = (vbx_uhalf_t**)malloc(13*sizeof(vbx_uhalf_t*));
    vbx_uhalf_t** v_4 = (vbx_uhalf_t**)malloc(9*sizeof(vbx_uhalf_t*));

    vbx_ubyte_t* v_1o = (vbx_ubyte_t*)vbx_sp_malloc(width*sizeof(vbx_ubyte_t));
    vbx_ubyte_t* v_2o = (vbx_ubyte_t*)vbx_sp_malloc(width*sizeof(vbx_ubyte_t));
    vbx_ubyte_t* v_4o = (vbx_ubyte_t*)vbx_sp_malloc(width*sizeof(vbx_ubyte_t));

    for (i=0; i<13; i++) {
        v_1[i]= v_mem + (0+i)*width;
    }

    vbx_dma_to_vector(v_1[0], img, 12*width*sizeof(vbx_uhalf_t));

    for (i=0; i<11; i++) {
        v_2[i]= v_mem + (13+i)*width;
    }
    for (i=0; i<9; i++) {
        v_4[i]= v_mem + (13+11+i)*width;
    }

    vbx_uhalf_t* v_tmp;

    vbx_set_vl(width);
    
    for(i=0; i<10; i++) {
        vbx(VVH, VADD, v_2[i], v_1[i], v_1[i+1]);
    }
    for(i=0; i<8; i++) {
        vbx(VVH, VADD, v_4[i], v_2[i], v_2[i+2]);
    }

    for(i=0; i<10; i++) {
        vbx(VVH, VCUSTOM1, v_1[i], v_1[i], v_1[i+1]);
    }
    vbx_set_vl(width+(vci_lanes*2*3));
    for(i=0; i<8; i++) {
        vbx(VVH, VCUSTOM2, v_2[i], v_2[i], v_2[i+2]);
    }
    vbx_set_vl(width+(vci_lanes*2*6));
    for(i=0; i<4; i++) {
        vbx(VVH, VCUSTOM3, v_4[i], v_4[i], v_4[i+4]);
    }
    vbx_set_vl(width);

    for(i=0; i<9; i++) {
        vbx(VVHB, VADD, v_1o, ((vbx_ubyte_t*)v_1[i] + 1), v_1[i+1]);
        vbx_dma_to_host(patterns[0]+i*width, v_1o, width*sizeof(vbx_ubyte_t));
    }
    
    for(i=0; i<6; i++) {
        vbx(VVHB, VADD, v_2o, ((vbx_ubyte_t*)v_2[i] + 1), v_2[i+2]);
        vbx_dma_to_host(patterns[1]+i*width, v_2o, width*sizeof(vbx_ubyte_t));
    }
    
    for (j = 0; j < height - 2 - 9; j++) {
        // sum-custom-interleave
        vbx(VVH, VADD, v_2[10], v_1[10], v_1[11]);
        vbx(VVH, VADD, v_4[8], v_2[8], v_2[10]);

        vbx(VVH, VCUSTOM1, v_1[10], v_1[10], v_1[11]);
        vbx_set_vl(width+(vci_lanes*2*3));
        vbx(VVH, VCUSTOM2, v_2[8], v_2[8], v_2[10]);
        vbx_set_vl(width+(vci_lanes*2*6));
        vbx(VVH, VCUSTOM3, v_4[4], v_4[4], v_4[8]);
        vbx_set_vl(width);

        vbx_dma_to_vector(v_1[12], img+(12+j)*width, width*sizeof(vbx_uhalf_t));

        vbx(VVHB, VADD, v_1o, ((vbx_ubyte_t*)v_1[9] + 1), v_1[10]);
        vbx(VVHB, VADD, v_2o, ((vbx_ubyte_t*)v_2[6] + 1), v_2[8]);
        vbx(VVHB, VADD, v_4o, ((vbx_ubyte_t*)v_4[0] + 1), v_4[4]);

        //swap pointers
        v_tmp = v_1[8];
        for (i=8; i<12; i++){
            v_1[i] = v_1[i+1];
        }
        v_1[12] = v_tmp;

        v_tmp = v_2[6];
        for (i=6; i<10; i++) {
            v_2[i] = v_2[i+1];
        }
        v_2[10] = v_tmp;

        v_tmp = v_4[0];
        for (i=0; i<8; i++) {
            v_4[i] = v_4[i+1];
        }
        v_4[8] = v_tmp;

        // dma out results
        vbx_dma_to_host(patterns[0] + (9+j)*width, v_1o, width*sizeof(vbx_byte_t));
        vbx_dma_to_host(patterns[1] + (6+j)*width, v_2o, width*sizeof(vbx_byte_t));
        vbx_dma_to_host(patterns[2] + (0+j)*width, v_4o, width*sizeof(vbx_byte_t));
    };

    vbx_sync();
    vbx_sp_free();
    free(v_1);
    free(v_2);
    free(v_4);

    return patterns;
}

int LBPPassCascade(image_t img, lbp_stage_t *cascade, pair_t p0, int max_stage)
{
	int i;

	for (i = 0; i < max_stage; i++) {
        if (!LBPPassStage(img, cascade[i], p0)) {
            return 0;
        }
    }

    return 1;
}

int LBPPassStage(image_t img, lbp_stage_t stage, pair_t p0)
{
	int i, sum_stage = 0;

	for (i = 0; i < stage.count; i++) {
		if (CheckLUT((unsigned int *)stage.feats[i].lut, SATBinaryPattern(img, &stage.feats[i], p0))) {
			sum_stage = sum_stage + stage.feats[i].fail;
        } else {
			sum_stage = sum_stage + stage.feats[i].pass;
        }
	}

	if (sum_stage >= stage.stageThreshold) {
		return 1;
	}

	return 0;
}

/* compute all rectangular neighbor sections of a given feature */
unsigned char SATBinaryPattern(image_t img, lbp_feat_t *feature, const pair_t p0)
{
    unsigned int n, neighbor, center;
    unsigned char pattern = 0;
    pair_t coeff[8] = {{0, 0}, {1, 0}, {2, 0}, {2, 1}, {2, 2}, {1, 2}, {0, 2}, {0, 1}};

    pair_t pc = {p0.x + feature->pos.size.x,
                 p0.y + feature->pos.size.y};
    center = SATValue(img, feature->pos, pc);

    for (n = 0; n < 8; n++) {
        pair_t p = {p0.x + coeff[n].x * feature->pos.size.x,
                    p0.y + coeff[n].y * feature->pos.size.y};
        neighbor = SATValue(img, feature->pos, p);

        if (neighbor >= center) {
            pattern += 1 << (7-n);
        }
    }

    return pattern;
}

/*
 *  a ---- b
 *  |      |
 *  d ---- c
 */
unsigned int SATValue(image_t img, const cell_t cell, const pair_t p0)
{
    unsigned int a, b, c, d;
    int left = p0.x + cell.src.x - 1;
    int top = p0.y + cell.src.y - 1;

    if (left < 0 || top < 0) {
        a = 0;
    } else {
        a = img.src[top * img.size.x + left];
    }

    if (top < 0) {
        b = 0;
    } else {
        b = img.src[top * img.size.x + left + cell.size.x];
    }

    c = img.src[(top + cell.size.y) * img.size.x + left + cell.size.x];

    if (left < 0) {
        d = 0;
    } else {
        d = img.src[(top + cell.size.y) * img.size.x + left];
    }

    return (a + c) - (b + d);
}

unsigned int SATNormalizedValue(image_t img, const cell_t cell, const pair_t p0)
{
    return SATValue(img, cell, p0) / (cell.size.x * cell.size.y);
}

int CheckLUT(unsigned int *lut, unsigned char value)
{
    unsigned char group = value >> 5;
    unsigned char index = value & 0x1f;

    return (lut[group] >> index) & 1;
}

vptr_word vector_row_lbp_2D(vptr_word v_int, vptr_word v_tmp, int search_width, int image_width, int vector_2D, lbp_stage_t *cascade, short max_stage)
{
	vptr_word v_add   = v_tmp + 0*image_width*vector_2D; // Holds pass or fail values to be added to stage sum
	vptr_word v_stage = v_tmp + 1*image_width*vector_2D; // Holds sum of features in a stages
	vptr_word v_final = v_tmp + 2*image_width*vector_2D; // Holds binary value if passed all stages
	vptr_word v_accum = v_tmp + 3*image_width*vector_2D; // Holds accumulated binary values if stages have been passed, used to exit early
	vptr_uword v_lut   = (vptr_uword)v_tmp + 4*image_width*vector_2D;
	vptr_uword v_lbp   = (vptr_uword)v_tmp + 5*image_width*vector_2D;
	vptr_word v_top   = v_tmp + 6*image_width*vector_2D;
	vptr_word v_bot   = v_tmp + 7*image_width*vector_2D;
	vptr_word v_p0    = v_tmp + 8*image_width*vector_2D;

    /* check for stage threshold */
	vbx_set_vl(search_width);
	vbx_set_2D(vector_2D, image_width*sizeof(vbx_word_t), image_width*sizeof(vbx_word_t), image_width*sizeof(vbx_word_t));

	//Zero components
	vbx_2D(SVW, VMOV, v_final, 1, NULL);

	//Run through stages
	int stage;
	for (stage=0; stage < max_stage; stage++) {

		//Zero out temporary binary stage pass
        vbx_2D(SVW, VMOV, v_stage, 0, NULL);

		int n, f;
		for (f = 0; f < cascade[stage].count; f++) {
            lbp_feat_t feat = cascade[stage].feats[f];

            vbx_2D(SVWU, VMOV, v_lbp,  0, NULL);
            vbx_2D(SVWU, VMOV, v_lut,  0, NULL);
			/* initalize values to be added to default fail value */
			vbx_2D(SVW, VMOV, v_add, (int)feat.fail, NULL);

            int dx = feat.pos.src.x;
            int dy = feat.pos.src.y;
            int dw = feat.pos.size.x;
            int dh = feat.pos.size.y;
            vptr_word v_a = v_int + image_width * (dy + dh - 1)   + (dx + dw - 1);
            vptr_word v_b = v_int + image_width * (dy + dh - 1)   + (dx + 2*dw - 1);
            vptr_word v_c = v_int + image_width * (dy + 2*dh - 1) + (dx + 2*dw - 1);
            vptr_word v_d = v_int + image_width * (dy + 2*dh - 1) + (dx + dw - 1);

            /* Aliases */
            vptr_uword v_shift = (vptr_uword)v_top;
            vptr_word v_sel = v_top;
            vptr_word v_group = v_bot;
            vptr_uword v_idx = (vptr_uword)v_bot;
            vptr_word v_px = (vptr_word)v_lut;

            int coeff_x[8] = {0, 1, 2, 2, 2, 1, 0, 0};
            int coeff_y[8] = {0, 0, 0, 1, 2, 2, 2, 1};
            int x, y;

            /* calc center */
            /* p0 = (a + c) - (b + d) */
            vbx_2D(VVW, VADD, v_top, v_a, v_c);
            vbx_2D(VVW, VADD, v_bot, v_b, v_d);
            vbx_2D(VVW, VSUB, v_p0, v_top, v_bot);

            /* calc each neighbor with center, and shift in to lbp */
            for (n = 0; n < 8; n++) {
                x = (dx + (coeff_x[n])*dw - 1);
                y = (dy + (coeff_y[n])*dh - 1);

                v_a = v_int + image_width * y + x;
                v_b = v_int + image_width * y + (x + dw);
                v_c = v_int + image_width * (y + dh) + (x + dw);
                v_d = v_int + image_width * (y + dh) + x;

                if(x >= 0 && y >= 0) {
                    vbx_2D(VVW, VADD, v_top, v_a, v_c);
                    vbx_2D(VVW, VADD, v_bot, v_b, v_d);
                    vbx_2D(VVW, VSUB, v_px, v_top, v_bot);
                } else if (x >= 0 && y < 0)  {
                    vbx_2D(VVW, VSUB, v_px, v_c, v_d);
                } else if (x < 0 && y >= 0)  {
                    vbx_2D(VVW, VSUB, v_px, v_c, v_b);
                } else {
                    vbx_2D(VVW, VMOV, v_px, v_c, NULL);
                }
                vbx_2D(VVW, VSUB, v_px, v_px, v_p0);

                vbx_2D(SVWU, VMOV, v_shift, 0, NULL);
                vbx_2D(SVW, VCMV_GEZ, (vptr_word)v_shift, 1 << (7-n), v_px);
                vbx_2D(VVWU, VADD, v_lbp, v_lbp, v_shift);
            }

            /* check if pattern is in lut */
            vbx_2D(SVWU, VSHR, (vptr_uword)v_group, 5, v_lbp);
            for (n = 0; n < 8; n++) {
                vbx_2D(SVW, VADD, v_sel, -n, v_group);
                vbx_2D(SVW, VCMV_Z, (vptr_word)v_lut, feat.lut[n], v_sel);
            }

            vbx_2D(SVWU, VAND, v_idx, 0x1f, v_lbp);
            vbx_2D(VVWU, VSHR, v_lut, v_idx, v_lut);
            vbx_2D(SVWU, VAND, v_lut, 1, v_lut);

			/* add either pass or fail sum to running stage total */
			vbx_2D(SVW, VCMV_LEZ, v_add, (int)feat.pass, (vptr_word)v_lut);
			vbx_2D(VVW, VADD, v_stage, v_stage, v_add);
		}

		/* final stage result */
		vbx_2D(SVW, VSUB, v_stage, cascade[stage].stageThreshold, v_stage);
		vbx_2D(SVW, VCMV_GTZ, v_final, 0, v_stage);

		/* exit early if entire group of rows has failed */
		vbx_acc_2D(VVW, VMOV, v_accum, v_final, NULL);
		vbx_sync();
		int accumulated = 0;
		for (n = 0; n < vector_2D; n++) {
			accumulated  = accumulated + v_accum[n*image_width];
		}
#if DEBUG
		if (!accumulated) {
			stage_count[stage] = stage_count[stage]+1;
			break;
		} else if (stage == max_stage-1) {
			stage_count[stage] = stage_count[stage]+1;
		}
#else
		if (!accumulated) break;
#endif
	}

	//Accumulate if the row has any valid pixels in the last pixel of the row
	//(the last 'window' pixels are unused)
	vbx_acc_2D(VVW, VMOV, v_final + image_width - 1, v_final, NULL);
	return v_final;
}

vptr_word vector_row_lbp_masked(vptr_word v_int, vptr_word v_tmp, int search_width, int image_width, int vector_2D, lbp_stage_t *cascade, short max_stage)
{
	vptr_word v_add   = v_tmp + 0*image_width*vector_2D; // Holds pass or fail values to be added to stage sum
	vptr_word v_stage = v_tmp + 1*image_width*vector_2D; // Holds sum of features in a stages
	vptr_word v_final = v_tmp + 2*image_width*vector_2D; // Holds binary value if passed all stages
	/* vptr_word v_accum = v_tmp + 3*image_width*vector_2D; // Holds accumulated binary values if stages have been passed, used to exit early */
	vptr_word v_lut   = v_tmp + 4*image_width*vector_2D;
	vptr_word v_lbp   = v_tmp + 5*image_width*vector_2D;
	vptr_word v_top   = v_tmp + 6*image_width*vector_2D;
	vptr_word v_bot   = v_tmp + 7*image_width*vector_2D;
	vptr_word v_p0    = v_tmp + 8*image_width*vector_2D;

	// clear mask status register in case previous valid data somehow
	int mask_status;
	vbx_get_mask_status(&mask_status);

	//create mask; nothing set in the image_width-win area
	vbx_set_2D(vector_2D, image_width*sizeof(vbx_word_t), image_width*sizeof(vbx_word_t), image_width*sizeof(vbx_word_t));
	vbx_set_vl(image_width-search_width);
	vbx_2D(SVW, VMOV, v_final+search_width, 1, NULL);
	vbx_set_vl(search_width);
	vbx_2D(SVW, VMOV,   v_final, 0, NULL);
	vbx_set_vl(image_width*vector_2D);
	vbx_setup_mask(VCMV_Z, v_final);

	//run through stages
	int stage;
	for(stage=0; stage < max_stage; stage++){


		//Zero out temporary binary stage pass
        vbx_masked(SVW, VMOV, v_stage, 0, NULL);

		int n, f;
		for (f = 0; f < cascade[stage].count; f++) {
            lbp_feat_t feat = cascade[stage].feats[f];

            vbx_masked(SVW, VMOV, v_lbp,  0, NULL);
            vbx_masked(SVW, VMOV, v_lut,  0, NULL);
			/* initalize values to be added to default fail value */
			vbx_masked(SVW, VMOV, v_add, (int)feat.fail, NULL);

            int dx = feat.pos.src.x;
            int dy = feat.pos.src.y;
            int dw = feat.pos.size.x;
            int dh = feat.pos.size.y;
            vptr_word v_a = v_int + image_width * (dy + dh - 1)   + (dx + dw - 1);
            vptr_word v_b = v_int + image_width * (dy + dh - 1)   + (dx + 2*dw - 1);
            vptr_word v_c = v_int + image_width * (dy + 2*dh - 1) + (dx + 2*dw - 1);
            vptr_word v_d = v_int + image_width * (dy + 2*dh - 1) + (dx + dw - 1);

            /* Aliases */
            vptr_word v_shift = v_top;
            vptr_word v_sel = v_top;
            vptr_word v_group = v_bot;
            vptr_word v_idx = v_bot;
            vptr_word v_px = v_lut;

            int coeff_x[8] = {0, 1, 2, 2, 2, 1, 0, 0};
            int coeff_y[8] = {0, 0, 0, 1, 2, 2, 2, 1};
            int x, y;

            /* calc center */
            /* p0 = (a + c) - (b + d) */
            vbx_masked(VVW, VADD, v_top, v_a, v_c);
            vbx_masked(VVW, VADD, v_bot, v_b, v_d);
            vbx_masked(VVW, VSUB, v_p0, v_top, v_bot);

            /* calc each neighbor with center, and shift in to lbp */
            for (n = 0; n < 8; n++) {
                x = (dx + (coeff_x[n])*dw - 1);
                y = (dy + (coeff_y[n])*dh - 1);

                v_a = v_int + image_width * y + x;
                v_b = v_int + image_width * y + (x + dw);
                v_c = v_int + image_width * (y + dh) + (x + dw);
                v_d = v_int + image_width * (y + dh) + x;

                if(x >= 0 && y >= 0) {
                    vbx_masked(VVW, VADD, v_top, v_a, v_c);
                    vbx_masked(VVW, VADD, v_bot, v_b, v_d);
                    vbx_masked(VVW, VSUB, v_px, v_top, v_bot);
                } else if (x >= 0 && y < 0)  {
                    vbx_masked(VVW, VSUB, v_px, v_c, v_d);
                } else if (x < 0 && y >= 0)  {
                    vbx_masked(VVW, VSUB, v_px, v_c, v_b);
                } else {
                    vbx_masked(VVW, VMOV, v_px, v_c, NULL);
                }
                vbx_masked(VVW, VSUB, v_px, v_px, v_p0);

                vbx_masked(SVW, VMOV, v_shift, 0, NULL);
                vbx_masked(SVW, VCMV_GEZ, v_shift, 1 << (7-n), v_px);
                vbx_masked(VVWU, VADD, (vptr_uword)v_lbp, (vptr_uword)v_lbp, (vptr_uword)v_shift);
            }

            vbx_masked(SVW, VSHR, v_group, 5, v_lbp);
            for (n = 0; n < 8; n++) {
                vbx_masked(SVW, VADD, v_sel, -n, v_group);
                vbx_masked(SVW, VCMV_Z, v_lut, feat.lut[n], v_sel);
            }

            vbx_masked(SVW, VAND, v_idx, 0x1f, v_lbp);
            vbx_masked(VVW, VSHR, v_lut, v_idx, v_lut);
            vbx_masked(SVW, VAND, v_lut, 1, v_lut);

			/* add either pass or fail sum to running stage total */
			vbx_masked(SVW, VCMV_LEZ, v_add, (int)feat.pass, v_lut);
			vbx_masked(VVW, VADD, v_stage, v_stage, v_add);
		}

		//final stage result
		vbx_masked(SVW, VSUB, v_stage, cascade[stage].stageThreshold, v_stage);
		//update mask with new existant values
		vbx_setup_mask_masked(VCMV_LEZ, v_stage);

		//exit early if entire group of rows has failed
		vbx_sync();
		vbx_get_mask_status(&mask_status);

#if DEBUG
		if (!mask_status) {
			stage_count[stage] = stage_count[stage]+1;
			break;
		} else if (stage == max_stage-1) {
			stage_count[stage] = stage_count[stage]+1;
		}
#else
		if (!mask_status) break;
#endif
	}
	//set to 1 anything still left
	vbx_masked(SVW, VMOV, v_final, 1, NULL);

	//accumulate if the row has any valid pixels in the last pixel of the row
	//(the last 'window' pixels are unused)
	vbx_set_vl(search_width);
	vbx_set_2D(vector_2D, image_width*sizeof(vbx_word_t), image_width*sizeof(vbx_word_t), image_width*sizeof(vbx_word_t));
	vbx_acc_2D(VVW, VMOV, v_final+image_width-1, v_final, NULL);
	return v_final;
}

vptr_word vector_row_lbp_restricted_masked(vptr_ubyte *v_lbp, vptr_word v_tmp, int offset, int search_width, int image_width, int vector_2D, lbp_stage_t *cascade, short max_stage)
{
#if !USE_BYTES
	vptr_word v_add   = v_tmp + 0*image_width*vector_2D; // Holds pass or fail values to be added to stage sum
	vptr_word v_stage = v_tmp + 1*image_width*vector_2D; // Holds sum of features in a stages
	vptr_word v_final = v_tmp + 2*image_width*vector_2D; // Holds binary value if passed all stages
	vptr_word v_lut   = v_tmp + 3*image_width*vector_2D;
#if !LUT_CI
	vptr_word v_idx   = v_tmp + 4*image_width*vector_2D;
	vptr_word v_sel   = v_tmp + 5*image_width*vector_2D;
    /* Aliases */
    vptr_word v_group = v_idx;
#endif
#else
    vptr_byte v_tmp_b = (vptr_byte)v_tmp;
	vptr_byte v_final = v_tmp_b + 0*image_width*vector_2D; // Holds binary value if passed all stages
	vptr_byte v_stage = v_tmp_b + 1*image_width*vector_2D; // Holds sum of features in a stages
	vptr_byte v_lut   = v_tmp_b + 2*image_width*vector_2D;
#if !LUT_CI
	vptr_byte v_add   = v_tmp_b + 6*image_width*vector_2D; // Holds pass or fail values to be added to stage sum
	vptr_byte v_idx   = v_tmp_b + 7*image_width*vector_2D;
	vptr_byte v_sel   = v_tmp_b + 11*image_width*vector_2D;
    /* Aliases */
    vptr_word v_group = (vptr_word)v_idx;
#endif
#endif
    vptr_ubyte v_pattern;

    lbp_feat_t feat;
    int f, dx, dy, dw, stage, total = 0;

	// clear mask status register in case previous valid data somehow
	int mask_status;
	vbx_get_mask_status(&mask_status);

	//create mask; nothing set in the image_width-win area
#if !USE_BYTES
	vbx_set_2D(vector_2D, image_width*sizeof(vbx_word_t), image_width*sizeof(vbx_word_t), image_width*sizeof(vbx_word_t));
	vbx_set_vl(image_width-search_width);
	vbx_2D(SVW, VMOV, v_final+search_width, 1, NULL);
	vbx_set_vl(search_width);
	vbx_2D(SVW, VMOV,   v_final, 0, NULL);
	vbx_set_vl(image_width*vector_2D);
	vbx_setup_mask(VCMV_Z, v_final);
#else
	vbx_set_2D(vector_2D, image_width*sizeof(vbx_byte_t), image_width*sizeof(vbx_byte_t), image_width*sizeof(vbx_byte_t));
	vbx_set_vl(image_width-search_width);
	vbx_2D(SVB, VMOV, ((vptr_byte)v_final)+search_width, 1, NULL);
	vbx_set_vl(search_width);
	vbx_2D(SVB, VMOV, v_final, 0, NULL);
	vbx_set_vl(image_width*vector_2D);
	vbx_setup_mask(VCMV_Z, (vptr_byte)v_final);
#endif

	//run through stages
	for(stage=0; stage < max_stage; stage++){

		//Zero out temporary binary stage pass
#if !USE_BYTES
#if BAD_NO_MASK
        vbx(SVW, VMOV, v_stage, 0, NULL);
#else
        vbx_masked(SVW, VMOV, v_stage, 0, NULL);
#endif
#else
#if BAD_NO_MASK
        vbx(SVB, VMOV, v_stage, 0, NULL);
#else
        vbx_masked(SVB, VMOV, v_stage, 0, NULL);
#endif
#endif

		for (f = 0; f < cascade[stage].count; f++) {
            feat = cascade[stage].feats[f];
            dx = feat.pos.src.x;
            dy = feat.pos.src.y;
            dw = feat.pos.size.x;

            /* quick hack, won't hold for 2^3 */
            v_pattern = v_lbp[dw>>1]+(dy*image_width+dx+offset);

			/* initalize values to be added to default fail value */
#if !USE_BYTES
#if BAD_NO_MASK
			vbx(SVW, VMOV, v_add, (int)feat.fail, NULL);
#else
			vbx_masked(SVW, VMOV, v_add, (int)feat.fail, NULL);
#endif
#else
#if !LUT_CI
#if BAD_NO_MASK
			vbx(SVB, VMOV, v_add, feat.fail, NULL);
#else
			vbx_masked(SVB, VMOV, v_add, feat.fail, NULL);
#endif
#endif
#endif
#if LUT_CI
#if BAD_NO_MASK
            if((total+f) < 256) {
                vbx(SVBU, VLBPLUT, (vptr_ubyte)v_lut, total+f, v_pattern);
            } else {
                vbx(SVBS, VLBPLUT, (vptr_ubyte)v_lut, (total+f)-256, v_pattern);
            }
#else
            if((total+f) < 256) {
                vbx_masked(SVBU, VLBPLUT, (vptr_ubyte)v_lut, total+f, v_pattern);
            } else {
                vbx_masked(SVBS, VLBPLUT, (vptr_ubyte)v_lut, (total+f)-256, v_pattern);
            }
#endif
#else
#if BAD_NO_MASK
            /* check if pattern is in lut */
            vbx(SVBU, VSHR, (vptr_ubyte)v_group, 5, v_pattern);
            int n;
            for (n = 0; n < 8; n++) {
                vbx(SVB, VADD, (vptr_byte)v_sel, -n, (vptr_byte)v_group);
                vbx(SVBW, VCMV_Z, (vptr_word)v_lut, feat.lut[n], (vptr_byte)v_sel);
            }

            vbx(SVBWU, VAND, (vptr_uword)v_idx, 0x1f, v_pattern);
            vbx(VVWB, VSHR, (vptr_byte)v_lut, (vptr_word)v_idx, (vptr_word)v_lut);
            vbx(SVB, VAND, (vptr_byte)v_lut, 1, (vptr_byte)v_lut);
#else
            /* check if pattern is in lut */
            vbx_masked(SVBU, VSHR, (vptr_ubyte)v_group, 5, v_pattern);
            int n;
            for (n = 0; n < 8; n++) {
                vbx_masked(SVB, VADD, (vptr_byte)v_sel, -n, (vptr_byte)v_group);
                vbx_masked(SVBW, VCMV_Z, (vptr_word)v_lut, feat.lut[n], (vptr_byte)v_sel);
            }

            vbx_masked(SVBWU, VAND, (vptr_uword)v_idx, 0x1f, v_pattern);
            vbx_masked(VVWB, VSHR, (vptr_byte)v_lut, (vptr_word)v_idx, (vptr_word)v_lut);
            vbx_masked(SVB, VAND, (vptr_byte)v_lut, 1, (vptr_byte)v_lut);
#endif
#endif
			/* add either pass or fail sum to running stage total */
#if !USE_BYTES
#if BAD_NO_MASK
			vbx(SVBW, VCMV_LEZ, v_add, (int)feat.pass, (vptr_byte)v_lut);
			vbx(VVW, VADD, v_stage, v_stage, v_add);
#else
			vbx_masked(SVBW, VCMV_LEZ, v_add, (int)feat.pass, (vptr_byte)v_lut);
			vbx_masked(VVW, VADD, v_stage, v_stage, v_add);
#endif
#else
#if !LUT_CI
#if BAD_NO_MASK
			vbx(SVB, VCMV_LEZ, v_add, feat.pass, v_lut);
#else
			vbx_masked(SVB, VCMV_LEZ, v_add, feat.pass, v_lut);
#endif
#if BAD_NO_MASK
			vbx(VVB, VADD, v_stage, v_stage, v_add);
#else
			vbx_masked(VVB, VADD, v_stage, v_stage, v_add);
#endif
#else
#if BAD_NO_MASK
			vbx(VVB, VADD, v_stage, v_stage, v_lut);
#else
			vbx_masked(VVB, VADD, v_stage, v_stage, v_lut);
#endif
#endif
#endif
		}
        total += cascade[stage].count;

		//final stage result
#if !USE_BYTES
#if BAD_NO_MASK
		vbx(SVW, VSUB, v_stage, cascade[stage].stageThreshold, v_stage);
#else
		vbx_masked(SVW, VSUB, v_stage, cascade[stage].stageThreshold, v_stage);
#endif

		//update mask with new existant values
		vbx_setup_mask_masked(VCMV_LEZ, v_stage);
#else
		//update mask with new existant values
		vbx_setup_mask_masked(VCMV_GEZ, ((vptr_byte)v_stage));
#endif

		//exit early if entire group of rows has failed
		vbx_sync();
		vbx_get_mask_status(&mask_status);

#if DEBUG
		if (!mask_status) {
			stage_count[stage] = stage_count[stage]+1;
			break;
		} else if (stage == max_stage-1) {
			stage_count[stage] = stage_count[stage]+1;
		}
#else
		if (!mask_status) break;
#endif
	}
#if !USE_BYTES
	//set to 1 anything still left
	vbx_masked(SVW, VMOV, v_final, 1, NULL);

	//accumulate if the row has any valid pixels in the last pixel of the row
	//(the last 'window' pixels are unused)
	vbx_set_vl(search_width);
	vbx_set_2D(vector_2D, image_width*sizeof(vbx_word_t), image_width*sizeof(vbx_word_t), image_width*sizeof(vbx_word_t));
	vbx_acc_2D(VVW, VMOV, v_final+image_width-1, v_final, NULL);
#else
	vbx_masked(SVB, VMOV, (vptr_byte)v_final, 1, NULL);
	vbx_set_vl(search_width);
	vbx_set_2D(vector_2D, image_width*sizeof(vbx_byte_t), image_width*sizeof(vbx_byte_t), image_width*sizeof(vbx_byte_t));
	vbx_acc_2D(VVB, VMOV, ((vptr_byte)v_final)+image_width-1, (vptr_byte)v_final, NULL);
#endif
	return (vptr_word)v_final;
}

vptr_word vector_row_lbp_restricted_2D(vptr_ubyte *v_lbp, vptr_word v_tmp, int offset, int search_width, int image_width, int vector_2D, lbp_stage_t *cascade, short max_stage)
{
	vptr_word v_add   = v_tmp + 0*image_width*vector_2D; // Holds pass or fail values to be added to stage sum
	vptr_word v_stage = v_tmp + 1*image_width*vector_2D; // Holds sum of features in a stages
	vptr_word v_final = v_tmp + 2*image_width*vector_2D; // Holds binary value if passed all stages
	vptr_word v_accum = v_tmp + 3*image_width*vector_2D; // Holds accumulated binary values if stages have been passed, used to exit early
	vptr_word v_lut   = v_tmp + 4*image_width*vector_2D;
#if !LUT_CI
	vptr_word v_idx   = v_tmp + 5*image_width*vector_2D;
	vptr_word v_sel   = v_tmp + 6*image_width*vector_2D;
    /* Aliases */
    vptr_word v_group = v_idx;
#endif
    vptr_ubyte v_pattern;

    lbp_feat_t feat;
    int n, f, dx, dy, dw, dh, stage, total = 0;

    /* check for stage threshold */
	vbx_set_vl(search_width);
	vbx_set_2D(vector_2D, image_width*sizeof(vbx_word_t), image_width*sizeof(vbx_word_t), image_width*sizeof(vbx_word_t));

	//Zero components
	vbx_2D(SVW, VMOV, v_final, 1, NULL);

	//Run through stages
	for (stage=0; stage < max_stage; stage++) {

		//Zero out temporary binary stage pass
        vbx_2D(SVW, VMOV, v_stage, 0, NULL);

		for (f = 0; f < cascade[stage].count; f++) {
            feat = cascade[stage].feats[f];
            dx = feat.pos.src.x;
            dy = feat.pos.src.y;
            dw = feat.pos.size.x;
            dh = feat.pos.size.y;

            /* quick hack, won't hold for >= 2^3 */
            v_pattern = v_lbp[dw>>1]+(dy*image_width+dx+offset);

			/* initalize values to be added to default fail value */
			vbx_2D(SVW, VMOV, v_add, (int)feat.fail, NULL);
#if LUT_CI
            vbx_2D(SVBU, VLBPLUT, (vptr_ubyte)v_lut, total+f, v_pattern);
#else
            /* check if pattern is in lut */
            vbx_2D(SVBU, VSHR, (vptr_ubyte)v_group, 5, v_pattern);
            for (n = 0; n < 8; n++) {
                vbx_2D(SVB, VADD, (vptr_byte)v_sel, -n, (vptr_byte)v_group);
                vbx_2D(SVBW, VCMV_Z, v_lut, feat.lut[n], (vptr_byte)v_sel);
            }

            vbx_2D(SVBWU, VAND, (vptr_uword)v_idx, 0x1f, v_pattern);
            vbx_2D(VVWB, VSHR, (vptr_byte)v_lut, v_idx, v_lut);
            vbx_2D(SVBU, VAND, (vptr_ubyte)v_lut, 1, (vptr_ubyte)v_lut);
#endif
			/* add either pass or fail sum to running stage total */
			vbx_2D(SVBW, VCMV_LEZ, v_add, (int)feat.pass, (vptr_byte)v_lut);
			vbx_2D(VVW, VADD, v_stage, v_stage, v_add);
		}
        total += cascade[stage].count;

		/* final stage result */
		vbx_2D(SVW, VSUB, v_stage, cascade[stage].stageThreshold, v_stage);
		vbx_2D(SVW, VCMV_GTZ, v_final, 0, v_stage);

		/* exit early if entire group of rows has failed */
		vbx_acc_2D(VVW, VMOV, v_accum, v_final, NULL);
		vbx_sync();
		int accumulated = 0;
		for (n = 0; n < vector_2D; n++) {
			accumulated  = accumulated + v_accum[n*image_width];
		}
#if DEBUG
		if (!accumulated) {
			stage_count[stage] = stage_count[stage]+1;
			break;
		} else if (stage == max_stage-1) {
			stage_count[stage] = stage_count[stage]+1;
		}
#else
		if (!accumulated) break;
#endif
	}

	//Accumulate if the row has any valid pixels in the last pixel of the row
	//(the last 'window' pixels are unused)
	vbx_acc_2D(VVW, VMOV, v_final + image_width - 1, v_final, NULL);
	return v_final;
}


int assign_lbp_lut_ci2(){

    vbx_mxp_t *this_mxp = VBX_GET_THIS_MXP();
    int vci_lanes = this_mxp->vcustom0_lanes;

    unsigned int a, b, c, d, x;
    int i, j, k, l, s, f, total = 0;
    

    /* pflx xxxx 0000 0000 0000 0000 0000 0000 */
    /* addressing scheme for loading the custom instruction */
    /* pass_base is bit 31 */
    /* fail_base is bit 30 */
    /* lut_base is bit 29 */
    /* which of the 32 bytes of 256 lut represented by bits 28-24 */

    unsigned int pass_base = 1 << 23;
    unsigned int fail_base = 1 << 22;
    unsigned int lut_base  = 1 << 21;
    unsigned int lut_start  = 1 << 16;

    int num_features = 0;
    for (s = 0; s < face_lbp_max_stage; s++) {
        num_features += face_lbp[s].count;
    }

    unsigned int *pass        = (unsigned int*)vbx_shared_malloc(num_features*sizeof(unsigned int));
    unsigned int *pass_addr   = (unsigned int*)vbx_shared_malloc(num_features*sizeof(unsigned int));
    unsigned int *fail        = (unsigned int*)vbx_shared_malloc(num_features*sizeof(unsigned int));
    unsigned int *fail_addr   = (unsigned int*)vbx_shared_malloc(num_features*sizeof(unsigned int));
    unsigned int *lut         = (unsigned int*)vbx_shared_malloc(32*num_features*sizeof(unsigned int));
    unsigned int *lut_addr    = (unsigned int*)vbx_shared_malloc(32*num_features*sizeof(unsigned int));
    unsigned int *l_pass      = (unsigned int*)vbx_shared_malloc(vci_lanes*num_features*sizeof(unsigned int));
    unsigned int *l_pass_addr = (unsigned int*)vbx_shared_malloc(vci_lanes*num_features*sizeof(unsigned int));
    unsigned int *l_fail      = (unsigned int*)vbx_shared_malloc(vci_lanes*num_features*sizeof(unsigned int));
    unsigned int *l_fail_addr = (unsigned int*)vbx_shared_malloc(vci_lanes*num_features*sizeof(unsigned int));
    unsigned int *l_lut       = (unsigned int*)vbx_shared_malloc(32*vci_lanes*num_features*sizeof(unsigned int));
    unsigned int *l_lut_addr  = (unsigned int*)vbx_shared_malloc(32*vci_lanes*num_features*sizeof(unsigned int));

    vbx_uword_t* v_data = (vbx_uword_t*)vbx_sp_malloc((num_features*vci_lanes)*sizeof(vbx_uword_t));
    vbx_uword_t* v_addr = (vbx_uword_t*)vbx_sp_malloc((num_features*vci_lanes)*sizeof(vbx_uword_t));
    if (v_addr == NULL) {
        printf("Failed\n");
    }

    for (s = 0; s < face_lbp_max_stage; s++) {
        for (f = 0; f < face_lbp[s].count; f++) {
            lbp_feat_t feat = face_lbp[s].feats[f];

            pass[f+total] = feat.pass;
            fail[f+total] = feat.fail;
            pass_addr[f+total] = pass_base + f+total;
            fail_addr[f+total] = fail_base + f+total;

            for(j = 0; j < 8; j++){
                x = feat.lut[j];

                a = x & 0xff;
                b = (x >> 8) & 0xff;
                c = (x >> 16) & 0xff;
                d = (x >> 24) & 0xff;

                lut[(f+total)*32+(j*4)+0] = a;
                lut[(f+total)*32+(j*4)+1] = b;
                lut[(f+total)*32+(j*4)+2] = c;
                lut[(f+total)*32+(j*4)+3] = d;

                for (k = 0; k < 4; k++) {
                    lut_addr[(f+total)*32+(j*4)+k] = lut_base + ((j*4+k) * lut_start) + f+total;
                }
            }
        }
        total += face_lbp[s].count;
    }

    for(i=0; i<num_features; i++){
        for(l=0; l<vci_lanes; l++){
            l_pass[i*vci_lanes+l] = pass[i];
            l_fail[i*vci_lanes+l] = fail[i];
            l_pass_addr[i*vci_lanes+l] = pass_addr[i];
            l_fail_addr[i*vci_lanes+l] = fail_addr[i];
        }
    }
    for(i=0; i<num_features*32; i++){
        for(l=0; l<vci_lanes; l++){
            l_lut[i*vci_lanes+l] = lut[i];
            l_lut_addr[i*vci_lanes+l] = lut_addr[i];
        }
    }
    vbx_set_vl((num_features*vci_lanes));

    vbx_dma_to_vector(v_data, l_pass, (num_features*vci_lanes)*sizeof(unsigned int));
    vbx_dma_to_vector(v_addr, l_pass_addr, (num_features*vci_lanes)*sizeof(unsigned int));
    vbx(VVWU, VLBPLUTWR, v_data, v_addr, v_data);

    vbx_dma_to_vector(v_data, l_fail, (num_features*vci_lanes)*sizeof(unsigned int));
    vbx_dma_to_vector(v_addr, l_fail_addr, (num_features*vci_lanes)*sizeof(unsigned int));
    vbx(VVWU, VLBPLUTWR, v_data, v_addr, v_data);


    for (s = 0; s < 32; s++) {
        vbx_dma_to_vector(v_data, l_lut + s*(num_features*vci_lanes), (num_features*vci_lanes)*sizeof(unsigned int));
        vbx_dma_to_vector(v_addr, l_lut_addr + s*(num_features*vci_lanes), (num_features*vci_lanes)*sizeof(unsigned int));
        vbx(VVWU, VLBPLUTWR, v_data, v_addr, v_data);
    }
    vbx_sync();

    vbx_sp_free();
    vbx_shared_free(l_pass);
    vbx_shared_free(l_fail);
    vbx_shared_free(l_lut);
    vbx_shared_free(l_pass_addr);
    vbx_shared_free(l_fail_addr);
    vbx_shared_free(l_lut_addr);
    return 0;
}

int assign_lbp_lut_ci(){
    vbx_mxp_t *this_mxp = VBX_GET_THIS_MXP();
    int vci_lanes = this_mxp->vcustom0_lanes;

    unsigned int a, b, c, d, x;
    int i, j, k, l, s, f, total = 0;
    

    /* pflx xxxx 0000 0000 0000 0000 0000 0000 */
    /* addressing scheme for loading the custom instruction */
    /* pass_base is bit 31 */
    /* fail_base is bit 30 */
    /* lut_base is bit 29 */
    /* which of the 32 bytes of 256 lut represented by bits 28-24 */

    unsigned int pass_base = 0x80000000;
    unsigned int fail_base = 0x40000000;
    unsigned int lut_base  = 0x20000000;
    unsigned int lut_start  = 24;

    int num_features = 0;
    for (s = 0; s < face_lbp_max_stage; s++) {
        num_features += face_lbp[s].count;
    }

    unsigned int *pass        = (unsigned int*)vbx_shared_malloc(num_features*sizeof(unsigned int));
    unsigned int *pass_addr   = (unsigned int*)vbx_shared_malloc(num_features*sizeof(unsigned int));
    unsigned int *fail        = (unsigned int*)vbx_shared_malloc(num_features*sizeof(unsigned int));
    unsigned int *fail_addr   = (unsigned int*)vbx_shared_malloc(num_features*sizeof(unsigned int));
    unsigned int *lut         = (unsigned int*)vbx_shared_malloc(32*num_features*sizeof(unsigned int));
    unsigned int *lut_addr    = (unsigned int*)vbx_shared_malloc(32*num_features*sizeof(unsigned int));
    unsigned int *l_pass      = (unsigned int*)vbx_shared_malloc(vci_lanes*num_features*sizeof(unsigned int));
    unsigned int *l_pass_addr = (unsigned int*)vbx_shared_malloc(vci_lanes*num_features*sizeof(unsigned int));
    unsigned int *l_fail      = (unsigned int*)vbx_shared_malloc(vci_lanes*num_features*sizeof(unsigned int));
    unsigned int *l_fail_addr = (unsigned int*)vbx_shared_malloc(vci_lanes*num_features*sizeof(unsigned int));
    unsigned int *l_lut       = (unsigned int*)vbx_shared_malloc(32*vci_lanes*num_features*sizeof(unsigned int));
    unsigned int *l_lut_addr  = (unsigned int*)vbx_shared_malloc(32*vci_lanes*num_features*sizeof(unsigned int));

    vbx_uword_t* v_data = (vbx_uword_t*)vbx_sp_malloc((num_features*vci_lanes)*sizeof(vbx_uword_t));
    vbx_uword_t* v_addr = (vbx_uword_t*)vbx_sp_malloc((num_features*vci_lanes)*sizeof(vbx_uword_t));
    if (v_addr == NULL) {
        printf("Failed\n");
    }

    for (s = 0; s < face_lbp_max_stage; s++) {
        for (f = 0; f < face_lbp[s].count; f++) {
            lbp_feat_t feat = face_lbp[s].feats[f];

            pass[f+total] = feat.pass & 0xff;
            fail[f+total] = feat.fail & 0xff;
            pass_addr[f+total] = pass_base + f+total;
            fail_addr[f+total] = fail_base + f+total;

            for(j = 0; j < 8; j++){
                x = feat.lut[j];

                a = x & 0xff;
                b = (x >> 8) & 0xff;
                c = (x >> 16) & 0xff;
                d = (x >> 24) & 0xff;

                lut[(f+total)*32+(j*4)+0] = a;
                lut[(f+total)*32+(j*4)+1] = b;
                lut[(f+total)*32+(j*4)+2] = c;
                lut[(f+total)*32+(j*4)+3] = d;

                for (k = 0; k < 4; k++) {
                    lut_addr[(f+total)*32+(j*4)+k] = lut_base + ((j*4+k) << lut_start) + f+total;
                }
            }
        }
        total += face_lbp[s].count;
    }

    for(i=0; i<num_features; i++){
        for(l=0; l<vci_lanes; l++){
            l_pass[i*vci_lanes+l] = pass[i];
            l_fail[i*vci_lanes+l] = fail[i];
            l_pass_addr[i*vci_lanes+l] = pass_addr[i];
            l_fail_addr[i*vci_lanes+l] = fail_addr[i];
        }
    }
    for(i=0; i<num_features*32; i++){
        for(l=0; l<vci_lanes; l++){
            l_lut[i*vci_lanes+l] = lut[i];
            l_lut_addr[i*vci_lanes+l] = lut_addr[i];
        }
    }
    vbx_set_vl((num_features*vci_lanes));

    vbx_dma_to_vector(v_data, l_pass, (num_features*vci_lanes)*sizeof(unsigned int));
    vbx_dma_to_vector(v_addr, l_pass_addr, (num_features*vci_lanes)*sizeof(unsigned int));
    vbx(VVWU, VLBPLUT, v_data, v_addr, v_data);

    vbx_dma_to_vector(v_data, l_fail, (num_features*vci_lanes)*sizeof(unsigned int));
    vbx_dma_to_vector(v_addr, l_fail_addr, (num_features*vci_lanes)*sizeof(unsigned int));
    vbx(VVWU, VLBPLUT, v_data, v_addr, v_data);


    for (s = 0; s < 32; s++) {
        vbx_dma_to_vector(v_data, l_lut + s*(num_features*vci_lanes), (num_features*vci_lanes)*sizeof(unsigned int));
        vbx_dma_to_vector(v_addr, l_lut_addr + s*(num_features*vci_lanes), (num_features*vci_lanes)*sizeof(unsigned int));
        vbx(VVWU, VLBPLUT, v_data, v_addr, v_data);
    }
    vbx_sync();

    vbx_sp_free();
    vbx_shared_free(l_pass);
    vbx_shared_free(l_fail);
    vbx_shared_free(l_lut);
    vbx_shared_free(l_pass_addr);
    vbx_shared_free(l_fail_addr);
    vbx_shared_free(l_lut_addr);
    return 0;
}

void vector_row_lbp_restricted_max(unsigned char *output, vptr_ubyte *v_lbp, vptr_byte v_tmp, int x_offset, int search_width, int image_width, int vector_2D, lbp_stage_t *cascade, short max_stage)
{
    lbp_feat_t feat, feat_a, feat_b;
    int f, dx, dy, dw, stage, feat_total = 0;

    /* pass/fail value for a given feature and pattern */
    vptr_byte v_lut   = v_tmp + 0*image_width*vector_2D; 
    /* sum of pass/fail values for features in a given stage */
    vptr_byte v_stage = v_tmp + 1*image_width*vector_2D;
    /* binary value if currently passing stages (1) or has failed (0) */
    vptr_byte v_final = v_tmp + 2*image_width*vector_2D; 
    /* set to the current features cell size and position within search window */
    vptr_byte v_pattern, v_pattern_a, v_pattern_b;

    /* clear mask status register in case previous valid data somehow */
    int mask_status;
    vbx_get_mask_status(&mask_status);

    /* create mask; nothing set in the image_width-win area */
    /* vbx_set_2D remains the same in the final accumation step of this function */
    vbx_set_2D(vector_2D, image_width, 0, 0);
    vbx_set_vl(image_width-search_width);
    vbx_2D(SVB, VMOV, v_final+search_width, 1, NULL);
    vbx_set_vl(search_width);
    vbx_2D(SVB, VMOV, v_final, 0, NULL);
    vbx_set_vl(image_width*vector_2D);
    vbx_setup_mask(VCMV_Z, v_final);

    /* run through stages */
    for(stage=0; stage < max_stage; stage++){
			
        /* zero out temporary binary stage pass */
#if BAD_NO_MASK
        vbx(SVB, VMOV, v_stage, 0, NULL);
#else
        vbx_masked(SVB, VMOV, v_stage, 0, NULL);
#endif

#if DOUBLE_LUT
        for (f = 0; f < cascade[stage].count/2; f++) {
            feat_a = cascade[stage].feats[f*2];
            dx = feat_a.pos.src.x;
            dy = feat_a.pos.src.y;
            dw = feat_a.pos.size.x;

            v_pattern_a = v_lbp[dw >> 1]+(dy * image_width + dx + x_offset);

            feat_b = cascade[stage].feats[f*2+1];
            dx = feat_b.pos.src.x;
            dy = feat_b.pos.src.y;
            dw = feat_b.pos.size.x;

            v_pattern_b = v_lbp[dw >> 1]+(dy * image_width + dx + x_offset);

			/* return pass or fail value, using the lbp pattern as index into the feature LUT */
            /* unsigned is used to access first 256 features, signed for next 256 features */
            if(feat_total == 0) {
#if BAD_NO_MASK
                vbx(VVBU, VLBPLUTZERO, (vptr_ubyte)v_lut, (vptr_ubyte)v_pattern_a, (vptr_ubyte)v_pattern_b);
#else
                vbx_masked(VVBU, VLBPLUTZERO, (vptr_ubyte)v_lut, (vptr_ubyte)v_pattern_a, (vptr_ubyte)v_pattern_b);
#endif

            } else {
#if BAD_NO_MASK
                vbx(VVBU, VLBPLUTINC, (vptr_ubyte)v_lut, (vptr_ubyte)v_pattern_a, (vptr_ubyte)v_pattern_b);
#else
                vbx_masked(VVBU, VLBPLUTINC, (vptr_ubyte)v_lut, (vptr_ubyte)v_pattern_a, (vptr_ubyte)v_pattern_b);
#endif
            }
			/* add the returned pass or fail value to the stage total */
#if BAD_NO_MASK
                vbx(VVB, VADD, v_stage, v_stage, v_lut);
#else
                vbx_masked(VVB, VADD, v_stage, v_stage, v_lut);
#endif
                feat_total++;
            }
#else
            for (f = 0; f < cascade[stage].count; f++) {
            feat = cascade[stage].feats[f];
            dx = feat.pos.src.x;
            dy = feat.pos.src.y;
            dw = feat.pos.size.x;

            /* quick hack to select which lbp pattern based on cell size */
            /* shifting dw won't hold for cell size >= 2^3, should use log2 */
            v_pattern = v_lbp[dw >> 1]+(dy * image_width + dx + x_offset);

			/* return pass or fail value, using the lbp pattern as index into the feature LUT */
            /* unsigned is used to access first 256 features, signed for next 256 features */
            if(feat_total < 256) {
#if BAD_NO_MASK
                vbx(SVBU, VLBPLUT, (vptr_ubyte)v_lut, feat_total % 256, (vptr_ubyte)v_pattern);
#else
                vbx_masked(SVBU, VLBPLUT, (vptr_ubyte)v_lut, feat_total % 256, (vptr_ubyte)v_pattern);
#endif

            } else {
#if BAD_NO_MASK
                vbx(SVBS, VLBPLUT,  (vptr_byte)v_lut, feat_total % 256,  (vptr_byte)v_pattern);
#else
                vbx_masked(SVBS, VLBPLUT,  (vptr_byte)v_lut, feat_total % 256,  (vptr_byte)v_pattern);
#endif
            }
			/* add the returned pass or fail value to the stage total */
#if BAD_NO_MASK
			vbx(VVB, VADD, v_stage, v_stage, v_lut);
#else
			vbx_masked(VVB, VADD, v_stage, v_stage, v_lut);
#endif
            feat_total++;
		}
#endif

		/* final stage result - update mask with new existant values */
		vbx_setup_mask_masked(VCMV_GEZ, v_stage);

		/* exit early if entire group of rows has failed */
		/* vbx_sync(); */
		vbx_get_mask_status(&mask_status);
#if !DEBUG
		if (!mask_status) break;
#else
		if (!mask_status) {
			stage_count[stage] = stage_count[stage]+1;
			break;
		} else if (stage == max_stage-1) {
			stage_count[stage] = stage_count[stage]+1;
		}
#endif
    }
    /* set to 1 anything still left */
    vbx_masked(SVB, VMOV, v_final, 1, NULL);
    /* accumulate if the row has any valid pixels in the last pixel of the row */
    /* the last 'window' width of pixels remain unused */

#if 1
    vbx_set_vl(search_width);
    vbx_set_2D(vector_2D, image_width, image_width, 0);
    vbx_acc_2D(VVB, VMOV, v_final+image_width-1, v_final, NULL);
#else
    int stride = 100;
    int num_strides = search_width / 100;
    vbx_set_vl(100);
    vbx_set_2D(scalar_width/100, 100, 100, 0);
    vbx_set_3D(vector_2D, image_width, image_width, 0);
    vbx_acc_3D(VVB, VMOV, v_check, v_final, NULL);
    vbx_dma_to_host(check, v_check, search_width / 100 * vector_2D);
    if (search_width%stride) {
        vbx_set_vl(100);
        vbx_set_2D(scalar_width/100, 100, 100, 0);
        vbx_set_3D(vector_2D, image_width, image_width, 0);
        vbx_acc_3D(VVB, VMOV, v_check, v_final, NULL);
        vbx_dma_to_host(check, v_check, search_width / 100 * vector_2D);
    }
#endif
    vbx_dma_to_host(output, v_final, image_width*vector_2D);
    vbx_sync();
}
