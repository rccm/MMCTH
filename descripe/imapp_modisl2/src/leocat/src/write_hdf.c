/*$Id: write_hdf.c,v 1.1.1.2 2006/12/05 14:27:49 mpav Exp $*/

#include <math.h>
#include <hdf.h>
#include <mfhdf.h>
#include "common_leocat.h"
#include "output_leocat.h"
#include <string.h>

void write_hdf (hdf_output hdf, void *var)

{

    char *rout = {"write_hdf"};

    int8
    *i8_buf;

    int16
    *i16_buf;

    int32
    sds_id,
    dim_id,
    status,
    n,
    *i32_buf,
    *start,
    *edge,
    npts;

    float32
    range_min=0.0,
    range_max=0.0;

    float32
    *f32_buf;


    float64
    scale,
    offset;
    float32 *sampleBuffer=NULL;

#if 0
 /*########################################################################### */
 #define SAMPLE 5
    char *name;
    int sampler=0;
    int samplerB=0;
    int sampleIdx=0;

    int sizea=0;
    int sizeb=0;

    name=hdf.sds_name;

printf("hdfName %s\n",hdf.sds_name);

    if (strstr(name,"Latitude") != NULL) {

        int aSize =  hdf.dimen[0]*hdf.dimen[1];
printf("latitude\n");
        sizea=hdf.dimen[0];
        sizeb=hdf.dimen[1];
        hdf.dimen[0] = hdf.dimen[0]/SAMPLE;
        hdf.dimen[1] = hdf.dimen[1]/SAMPLE;


        if ((sampleBuffer = (float32 *) malloc(aSize*sizeof(int32))) == NULL)
            error_allo(rout,"sampleBuffer");

        sampleIdx=0;
        for (sampler=0;sampler<hdf.dimen[0]*SAMPLE;sampler+=SAMPLE) {

            for (samplerB=0;samplerB<hdf.dimen[1]*SAMPLE;samplerB+=SAMPLE) {
                sampleBuffer[sampleIdx]=
                    ((float*)var)[((sampler)*sizeb)+samplerB];


                sampleIdx++;
            }

        }
        var = sampleBuffer;

    } else if (strstr(name,"Longitude") != NULL) {

        int aSize =  hdf.dimen[0]*hdf.dimen[1];

        sizea=hdf.dimen[0];
        sizeb=hdf.dimen[1];
        hdf.dimen[0] = hdf.dimen[0]/SAMPLE;
        hdf.dimen[1] = hdf.dimen[1]/SAMPLE;


        if ((sampleBuffer = (float32 *) malloc(aSize*sizeof(int32))) == NULL)
            error_allo(rout,"sampleBuffer");

        sampleIdx=0;
        for (sampler=0;sampler<hdf.dimen[0]*SAMPLE;sampler+=SAMPLE) {

            for (samplerB=0;samplerB<hdf.dimen[1]*SAMPLE;samplerB+=SAMPLE) {
                sampleBuffer[sampleIdx]=
                    ((float*)var)[((sampler)*sizeb)+samplerB];

                sampleIdx++;
            }

        }
        var = sampleBuffer;

    }

 /*########################################################################### */#define SAMPLE 5

#endif

//  printf("\nInside write_hdf(), for %s\n",hdf.sds_name);	// GPC

    sds_id = SDcreate(hdf.sd_id,hdf.sds_name,hdf.scaled_type,hdf.rank,hdf.dimen);
    if (sds_id == FAIL) {
        fprintf(stderr,"%sCannot create HDF SDS, %s - aborting\n",EXE_PROMPT,hdf.sds_name);
        exit(EXIT_FAILURE);
    }

    if ((start = (int32 *) calloc(hdf.rank,sizeof(int32))) == NULL)
        error_allo(rout,"start");
    if ((edge = (int32 *) calloc(hdf.rank,sizeof(int32))) == NULL)
        error_allo(rout,"edge");

    npts = 1;
    for (n=0; n<hdf.rank; n++) {
        dim_id = SDgetdimid(sds_id,n);
        if (dim_id == FAIL) {
            fprintf(stderr,"%sCannot create dimension %d for HDF SDS, %s - aborting\n",EXE_PROMPT,n,hdf.sds_name);
            exit(EXIT_FAILURE);
        }

//        fprintf(stderr,"%sdata0: %d %d\n", EXE_PROMPT, n, dim_id);
        status = SDsetdimname(dim_id,hdf.dim_name[n]);
        if (status == FAIL) {
            fprintf(stderr,"%sCannot set name for dimension %d for HDF SDS, %s - aborting\n",EXE_PROMPT,n,hdf.sds_name);
            exit(EXIT_FAILURE);
        }
        start[n] = 0;
        edge[n] = hdf.dimen[n];
        npts *= hdf.dimen[n];
    }

    status = SDsetattr(sds_id,"algorithm_index",DFNT_INT16,1,(void *) &hdf.algo_num);
    if (status == FAIL) {
        fprintf(stderr,"%sCannot create algorithm_index attribute for SDS, %s - aborting\n",EXE_PROMPT,hdf.sds_name);
        exit(EXIT_FAILURE);
    }

    status = SDsetattr(sds_id,"reference",DFNT_CHAR8,strlen(hdf.reference),(void *) &hdf.reference[0]);
    if (status == FAIL) {
        fprintf(stderr,"%sCannot create reference attribute for SDS, %s - aborting\n",EXE_PROMPT,hdf.sds_name);
        exit(EXIT_FAILURE);
    }

    status = SDsetattr(sds_id,"units",DFNT_CHAR8,strlen(hdf.units),(void *) &hdf.units[0]);
    if (status == FAIL) {
        fprintf(stderr,"%sCannot create units attribute for SDS, %s - aborting\n",EXE_PROMPT,hdf.sds_name);
        exit(EXIT_FAILURE);
    }

    status = SDsetattr(sds_id,"scaling_method",DFNT_INT8,1,(void *) &hdf.scaled_flg);
    if (status == FAIL) {
        fprintf(stderr,"%sCannot create scaling_method attribute for SDS, %s - aborting\n",EXE_PROMPT,hdf.sds_name);
        exit(EXIT_FAILURE);
    }

    switch (hdf.scaled_type) {
    case DFNT_INT8:
        hdf.scaled_min = scaled_min_i8;
        hdf.scaled_max = scaled_max_i8;
        hdf.scaled_missing = scaled_missing_i8;
        break;
    case DFNT_CHAR8:
        hdf.scaled_min = scaled_min_i8;
        hdf.scaled_max = scaled_max_i8;
        hdf.scaled_missing = scaled_missing_i8;
        break;
    case DFNT_INT16:
        hdf.scaled_min = scaled_min_i16;
        hdf.scaled_max = scaled_max_i16;
        hdf.scaled_missing = scaled_missing_i16;
        break;
    case DFNT_INT32:
        hdf.scaled_missing = MISSING_INT;
        break;
    case DFNT_FLOAT32:
        hdf.scaled_missing = MISSING_FLOAT;
        break;
    case DFNT_FLOAT64:
        hdf.scaled_missing = MISSING_FLOAT;
        break;
    default:
        fprintf(stderr,"%sScaled data type, %d is not valid - aborting",EXE_PROMPT,hdf.scaled_type);
        exit(EXIT_FAILURE);
        break;
    }

//      printf("\nhdf.scaled_flg = %d\n",hdf.scaled_flg);	// GPC
//      printf("\nhdf.scaled_type = %d\n",hdf.scaled_type);	// GPC
    if (hdf.scaled_flg >= 1 && hdf.scaled_flg <=3) {
//  	printf("\nhdf.scaled_flg >= 1 && hdf.scaled_flg <=3\n");	// GPC

        switch (hdf.scaled_flg) {
        case 1:
//      	printf("\trange = hdf.range\n");	// GPC
            range_min = hdf.range_min;
            range_max = hdf.range_max;
            break;
        case 2:
//      	printf("\trange = log10(hdf.range)\n");	// GPC
            range_min = log10(hdf.range_min);
            range_max = log10(hdf.range_max);
            break;
        case 3:
//      	printf("\trange = sqrt(hdf.range)\n");	// GPC
            range_min = sqrt(hdf.range_min);
            range_max = sqrt(hdf.range_max);
            break;
        }

        scale = (range_max - range_min)/(hdf.scaled_max - hdf.scaled_min);
        offset = (float64) hdf.scaled_min - (range_min/scale);
        offset = (-1.0)*offset*scale;

//  	printf("\tscale = %f\n\toffset = %f\n",scale,offset);	// GPC

        status = SDsetattr(sds_id,"scale_factor",DFNT_FLOAT64,1,(void *) &scale);
        if (status == FAIL) {
            fprintf(stderr,"%sCannot create scale_factor attribute for SDS, %s - aborting\n",EXE_PROMPT,hdf.sds_name);
            exit(EXIT_FAILURE);
        }

        status = SDsetattr(sds_id,"add_offset",DFNT_FLOAT64,1,(void *) &offset);
        if (status == FAIL) {
            fprintf(stderr,"%sCannot create add_offset attribute for SDS, %s - aborting\n",EXE_PROMPT,hdf.sds_name);
            exit(EXIT_FAILURE);
        }

        status = SDsetattr(sds_id,"_FillValue",hdf.scaled_type,1,(void *) &hdf.scaled_missing);
        if (status == FAIL) {
            fprintf(stderr,"%sCannot create _FillValue attribute for SDS, %s - aborting\n",EXE_PROMPT,hdf.sds_name);
            exit(EXIT_FAILURE);
        }

        if (hdf.scaled_type >= 1) {
            offset = (-1.0)*offset/scale;
        }

        f32_buf = (float32 *) var;

        switch (hdf.scaled_type) {
        case DFNT_INT8:
            if ((i8_buf = (int8 *) malloc(npts*sizeof(int8))) == NULL)
                error_allo(rout,"i8_buf");
            for (n=0; n<npts; n++) {
                if (f32_buf[n] == MISSING_FLOAT) {
                    i8_buf[n] = hdf.scaled_missing;
                } else {
                    f32_buf[n] = max(hdf.range_min,min(hdf.range_max,f32_buf[n]));
                    switch (hdf.scaled_flg) {
                    case 1:
                        i8_buf[n] = nint(f32_buf[n]/scale + offset);
                        break;
                    case 2:
                        i8_buf[n] = (f32_buf[n]);
                        if (f32_buf[n] > 0) i8_buf[n] = nint(log10(f32_buf[n])/scale + offset);
                        break;
                    case 3:
                        i8_buf[n] = (f32_buf[n]);
                        if (f32_buf[n] > 0) i8_buf[n] = nint(sqrt(f32_buf[n])/scale + offset);
                        break;
                    }
                }
            }
            status = SDwritedata(sds_id,start,NULL,edge,(void *) i8_buf);
            if (status == FAIL) {
                fprintf(stderr,"%sError writing DFNT_INT8 HDF SDS, %s - aborting\n",EXE_PROMPT,hdf.sds_name);
                exit(EXIT_FAILURE);
            }
            free(i8_buf);
            break;
        case DFNT_CHAR8:
            if ((i8_buf = (int8 *) malloc(npts*sizeof(int8))) == NULL)
                error_allo(rout,"i8_buf");
            for (n=0; n<npts; n++) {
                if (f32_buf[n] == MISSING_FLOAT) {
                    i8_buf[n] = hdf.scaled_missing;
                } else {
                    f32_buf[n] = max(hdf.range_min,min(hdf.range_max,f32_buf[n]));
                    switch (hdf.scaled_flg) {
                    case 1:
                        i8_buf[n] = nint(f32_buf[n]/scale + offset);
                        break;
                    case 2:
                        i8_buf[n] = (f32_buf[n]);
                        if (f32_buf[n] > 0) i8_buf[n] = nint(log10(f32_buf[n])/scale + offset);
                        break;
                    case 3:
                        i8_buf[n] = (f32_buf[n]);
                        if (f32_buf[n] > 0) i8_buf[n] = nint(sqrt(f32_buf[n])/scale + offset);
                        break;
                    }
                }
            }
            status = SDwritedata(sds_id,start,NULL,edge,(void *) i8_buf);
            if (status == FAIL) {
                fprintf(stderr,"%sError writing DFNT_CHAR8 HDF SDS, %s - aborting\n",EXE_PROMPT,hdf.sds_name);
                exit(EXIT_FAILURE);
            }
            free(i8_buf);
            break;
        case DFNT_INT16:
            if ((i16_buf = (int16 *) malloc(npts*sizeof(int16))) == NULL)
                error_allo(rout,"i16_buf");
            for (n=0; n<npts; n++) {
                if (f32_buf[n] == MISSING_FLOAT) {
                    i16_buf[n] = hdf.scaled_missing;
                } else {
                    f32_buf[n] = max(hdf.range_min,min(hdf.range_max,f32_buf[n]));
                    switch (hdf.scaled_flg) {
                    case 1:
                        i16_buf[n] = nint(f32_buf[n]/scale + offset);
                        break;
                    case 2:
                        i16_buf[n] = (f32_buf[n]);
                        if (f32_buf[n] > 0) i16_buf[n] = nint(log10(f32_buf[n])/scale + offset);
                        break;
                    case 3:
                        i16_buf[n] = (f32_buf[n]);
                        if (f32_buf[n] > 0) i16_buf[n] = nint(sqrt(f32_buf[n])/scale + offset);
                        break;
                    }
                }
            }
            status = SDwritedata(sds_id,start,NULL,edge,(void *) i16_buf);
            if (status == FAIL) {
                fprintf(stderr,"%sError writing DFNT_INT16 HDF SDS, %s - aborting\n",EXE_PROMPT,hdf.sds_name);
                exit(EXIT_FAILURE);
            }
            free(i16_buf);
            break;
        case DFNT_INT32:
            if ((i32_buf = (int32 *) malloc(npts*sizeof(int32))) == NULL)
                error_allo(rout,"i32_buf");
            for (n=0; n<npts; n++) {
                if (f32_buf[n] == MISSING_FLOAT) {
                    i32_buf[n] = hdf.scaled_missing;
                } else {
                    f32_buf[n] = max(hdf.range_min,min(hdf.range_max,f32_buf[n]));
                    switch (hdf.scaled_flg) {
                    case 1:
                        i32_buf[n] = nint(f32_buf[n]/scale + offset);
                        break;
                    case 2:
                        i32_buf[n] = (f32_buf[n]);
                        if (f32_buf[n] > 0) i32_buf[n] = nint(log10(f32_buf[n])/scale + offset);
                        break;
                    case 3:
                        i32_buf[n] = (f32_buf[n]);
                        if (f32_buf[n] > 0) i32_buf[n] = nint(sqrt(f32_buf[n])/scale + offset);
                        break;
                    }
                }
            }
            status = SDwritedata(sds_id,start,NULL,edge,(void *) i32_buf);
            if (status == FAIL) {
                fprintf(stderr,"%sError writing DFNT_INT32 HDF SDS, %s - aborting\n",EXE_PROMPT,hdf.sds_name);
                exit(EXIT_FAILURE);
            }
            free(i32_buf);
            break;
        default:
            fprintf(stderr,"%sScaled data type, %d is not a 8, 16, or 32 bit signed integer - aborting",EXE_PROMPT,hdf.scaled_type);
            exit(EXIT_FAILURE);
            break;
        }

    } else {

        scale = 1.0;
        offset = 0.0;

        status = SDsetattr(sds_id,"scale_factor",DFNT_FLOAT64,1,(void *) &scale);
        if (status == FAIL) {
            fprintf(stderr,"%sCannot create scale_factor attribute for SDS, %s - aborting\n",EXE_PROMPT,hdf.sds_name);
            exit(EXIT_FAILURE);
        }

        status = SDsetattr(sds_id,"add_offset",DFNT_FLOAT64,1,(void *) &offset);
        if (status == FAIL) {
            fprintf(stderr,"%sCannot create add_offset attribute for SDS, %s - aborting\n",EXE_PROMPT,hdf.sds_name);
            exit(EXIT_FAILURE);
        }

        status = SDsetattr(sds_id,"_FillValue",hdf.scaled_type,1,(void *) &hdf.scaled_missing);
        if (status == FAIL) {
            fprintf(stderr,"%sCannot create _FillValue attribute for SDS, %s - aborting\n",EXE_PROMPT,hdf.sds_name);
            exit(EXIT_FAILURE);
        }

        status = SDwritedata(sds_id,start,NULL,edge,var);
        if (status == FAIL) {
            fprintf(stderr,"%sError writinge HDF SDS, %s - aborting\n",EXE_PROMPT,hdf.sds_name);
            exit(EXIT_FAILURE);
        }
    }

    status = SDendaccess(sds_id);
    if (status == FAIL) {
        fprintf(stderr,"%sError ending access to HDF SDS, %s - aborting\n",EXE_PROMPT,hdf.sds_name);
        exit(EXIT_FAILURE);
    }
if ( sampleBuffer != NULL) free(sampleBuffer);


    free(start);
    free(edge);

}

void write_hdf_multidim (hdf_output hdf, int32 *start, int32 *edge, int32 *sds_id, void *var)

{

    char *rout = {"write_hdf_2d"};

    int8
    write_attr = 0;

    int8
    *i8_buf;

    int16
    *i16_buf;

    int32
    dim_id,
    status,
    n,
    *i32_buf,
    npts;

    float32
    range_min=0.0,
    range_max=0.0;

    float32
    *f32_buf;


    float64
    scale,
    offset;

    npts = 1;
    for (n=0; n<hdf.rank; n++) npts *= edge[n];

    if (*sds_id < 0) {

        write_attr = 1;

        *sds_id = SDcreate(hdf.sd_id,hdf.sds_name,hdf.scaled_type,hdf.rank,hdf.dimen);
        if (*sds_id == FAIL) {
            fprintf(stderr,"%sCannot create HDF SDS, %s - aborting\n",EXE_PROMPT,hdf.sds_name);
            exit(EXIT_FAILURE);
        }

        for (n=0; n<hdf.rank; n++) {
            dim_id = SDgetdimid(*sds_id,n);
            if (dim_id == FAIL) {
                fprintf(stderr,"%sCannot create dimension %d for HDF SDS, %s - aborting\n",EXE_PROMPT,n,hdf.sds_name);
                exit(EXIT_FAILURE);
            }
//            fprintf(stderr,"%sdata0: %d %d\n", EXE_PROMPT, n, dim_id);
            status = SDsetdimname(dim_id,hdf.dim_name[n]);
            if (status == FAIL) {
                fprintf(stderr,"%sCannot set name for dimension %d for HDF SDS, %s - aborting\n",EXE_PROMPT,n,hdf.sds_name);
                fprintf(stderr,"%sdata: %d %s\n", EXE_PROMPT, dim_id, hdf.dim_name[n]);
                fprintf(stderr,"%sdata: %d %d %d %d %d %d\n", EXE_PROMPT, start[0],start[1],start[2],edge[0],edge[1],edge[2]);
                fprintf(stderr,"%sdata: %d \n", EXE_PROMPT, *(sds_id));
                exit(EXIT_FAILURE);
            }
        }

        status = SDsetattr(*sds_id,"algorithm_index",DFNT_INT16,1,(void *) &hdf.algo_num);
        if (status == FAIL) {
            fprintf(stderr,"%sCannot create algorithm_index attribute for SDS, %s - aborting\n",EXE_PROMPT,hdf.sds_name);
            exit(EXIT_FAILURE);
        }

        status = SDsetattr(*sds_id,"reference",DFNT_CHAR8,strlen(hdf.reference),(void *) &hdf.reference[0]);
        if (status == FAIL) {
            fprintf(stderr,"%sCannot create reference attribute for SDS, %s - aborting\n",EXE_PROMPT,hdf.sds_name);
            exit(EXIT_FAILURE);
        }

        status = SDsetattr(*sds_id,"units",DFNT_CHAR8,strlen(hdf.units),(void *) &hdf.units[0]);
        if (status == FAIL) {
            fprintf(stderr,"%sCannot create units attribute for SDS, %s - aborting\n",EXE_PROMPT,hdf.sds_name);
            exit(EXIT_FAILURE);
        }

        status = SDsetattr(*sds_id,"scaling_method",DFNT_INT8,1,(void *) &hdf.scaled_flg);
        if (status == FAIL) {
            fprintf(stderr,"%sCannot create scaling_method attribute for SDS, %s - aborting\n",EXE_PROMPT,hdf.sds_name);
            exit(EXIT_FAILURE);
        }
    }

    switch (hdf.scaled_type) {
    case DFNT_INT8:
        hdf.scaled_min = scaled_min_i8;
        hdf.scaled_max = scaled_max_i8;
        hdf.scaled_missing = scaled_missing_i8;
        break;
    case DFNT_CHAR8:
        hdf.scaled_min = scaled_min_i8;
        hdf.scaled_max = scaled_max_i8;
        hdf.scaled_missing = scaled_missing_i8;
        break;
    case DFNT_INT16:
        hdf.scaled_min = scaled_min_i16;
        hdf.scaled_max = scaled_max_i16;
        hdf.scaled_missing = scaled_missing_i16;
        break;
    case DFNT_INT32:
        hdf.scaled_missing = MISSING_INT;
        break;
    case DFNT_FLOAT32:
        hdf.scaled_missing = MISSING_FLOAT;
        break;
    case DFNT_FLOAT64:
        hdf.scaled_missing = MISSING_FLOAT;
        break;
    default:
        fprintf(stderr,"%sScaled data type, %d is not valid - aborting",EXE_PROMPT,hdf.scaled_type);
        exit(EXIT_FAILURE);
        break;
    }

    if (hdf.scaled_flg >= 1 && hdf.scaled_flg <=3) {

        switch (hdf.scaled_flg) {
        case 1:
            range_min = hdf.range_min;
            range_max = hdf.range_max;
            break;
        case 2:
            range_min = log10(hdf.range_min);
            range_max = log10(hdf.range_max);
            break;
        case 3:
            range_min = sqrt(hdf.range_min);
            range_max = sqrt(hdf.range_max);
            break;
        }

        scale = (range_max - range_min)/(hdf.scaled_max - hdf.scaled_min);
        offset = (float64) hdf.scaled_min - (range_min/scale);
        offset = (-1.0)*offset*scale;

        if (write_attr) {
            status = SDsetattr(*sds_id,"scale_factor",DFNT_FLOAT64,1,(void *) &scale);
            if (status == FAIL) {
                fprintf(stderr,"%sCannot create scale_factor attribute for SDS, %s - aborting\n",EXE_PROMPT,hdf.sds_name);
                exit(EXIT_FAILURE);
            }

            status = SDsetattr(*sds_id,"add_offset",DFNT_FLOAT64,1,(void *) &offset);
            if (status == FAIL) {
                fprintf(stderr,"%sCannot create add_offset attribute for SDS, %s - aborting\n",EXE_PROMPT,hdf.sds_name);
                exit(EXIT_FAILURE);
            }

            status = SDsetattr(*sds_id,"_FillValue",hdf.scaled_type,1,(void *) &hdf.scaled_missing);
            if (status == FAIL) {
                fprintf(stderr,"%sCannot create _FillValue attribute for SDS, %s - aborting\n",EXE_PROMPT,hdf.sds_name);
                exit(EXIT_FAILURE);
            }
        }

        if (hdf.scaled_type >= 1) {
            offset = (-1.0)*offset/scale;
        }

        f32_buf = (float32 *) var;

        switch (hdf.scaled_type) {
        case DFNT_INT8:
            if ((i8_buf = (int8 *) malloc(npts*sizeof(int8))) == NULL)
                error_allo(rout,"i8_buf");
            for (n=0; n<npts; n++) {
                if (f32_buf[n] == MISSING_FLOAT) {
                    i8_buf[n] = hdf.scaled_missing;
                } else {
                    f32_buf[n] = max(hdf.range_min,min(hdf.range_max,f32_buf[n]));
                    switch (hdf.scaled_flg) {
                    case 1:
                        i8_buf[n] = nint(f32_buf[n]/scale + offset);
                        break;
                    case 2:
                        i8_buf[n] = (f32_buf[n]);
                        if (f32_buf[n] > 0) i8_buf[n] = nint(log10(f32_buf[n])/scale + offset);
                        break;
                    case 3:
                        i8_buf[n] = (f32_buf[n]);
                        if (f32_buf[n] > 0) i8_buf[n] = nint(sqrt(f32_buf[n])/scale + offset);
                        break;
                    }
                }
            }
            status = SDwritedata(*sds_id,start,NULL,edge,(void *) i8_buf);
            if (status == FAIL) {
                fprintf(stderr,"%sError writingf HDF SDS, %s - aborting\n",EXE_PROMPT,hdf.sds_name);
                exit(EXIT_FAILURE);
            }
            free(i8_buf);
            break;
        case DFNT_CHAR8:
            if ((i8_buf = (int8 *) malloc(npts*sizeof(int8))) == NULL)
                error_allo(rout,"i8_buf");
            for (n=0; n<npts; n++) {
                if (f32_buf[n] == MISSING_FLOAT) {
                    i8_buf[n] = hdf.scaled_missing;
                } else {
                    f32_buf[n] = max(hdf.range_min,min(hdf.range_max,f32_buf[n]));
                    switch (hdf.scaled_flg) {
                    case 1:
                        i8_buf[n] = nint(f32_buf[n]/scale + offset);
                        break;
                    case 2:
                        i8_buf[n] = (f32_buf[n]);
                        if (f32_buf[n] > 0) i8_buf[n] = nint(log10(f32_buf[n])/scale + offset);
                        break;
                    case 3:
                        i8_buf[n] = (f32_buf[n]);
                        if (f32_buf[n] > 0) i8_buf[n] = nint(sqrt(f32_buf[n])/scale + offset);
                        break;
                    }
                }
            }
            status = SDwritedata(*sds_id,start,NULL,edge,(void *) i8_buf);
            if (status == FAIL) {
                fprintf(stderr,"%sError writingg HDF SDS, %s - aborting\n",EXE_PROMPT,hdf.sds_name);
                exit(EXIT_FAILURE);
            }
            free(i8_buf);
            break;
        case DFNT_INT16:
            if ((i16_buf = (int16 *) malloc(npts*sizeof(int16))) == NULL)
                error_allo(rout,"i16_buf");
            for (n=0; n<npts; n++) {
                if (f32_buf[n] == MISSING_FLOAT) {
                    i16_buf[n] = hdf.scaled_missing;
                } else {
                    f32_buf[n] = max(hdf.range_min,min(hdf.range_max,f32_buf[n]));
                    switch (hdf.scaled_flg) {
                    case 1:
                        i16_buf[n] = nint(f32_buf[n]/scale + offset);
                        break;
                    case 2:
                        i16_buf[n] = (f32_buf[n]);
                        if (f32_buf[n] > 0) i16_buf[n] = nint(log10(f32_buf[n])/scale + offset);
                        break;
                    case 3:
                        i16_buf[n] = (f32_buf[n]);
                        if (f32_buf[n] > 0) i16_buf[n] = nint(sqrt(f32_buf[n])/scale + offset);
                        break;
                    }
                }
            }
            status = SDwritedata(*sds_id,start,NULL,edge,(void *) i16_buf);
            if (status == FAIL) {
                fprintf(stderr,"%sError writingh HDF SDS, %s - aborting\n",EXE_PROMPT,hdf.sds_name);
                exit(EXIT_FAILURE);
            }
            free(i16_buf);
            break;
        case DFNT_INT32:
            if ((i32_buf = (int32 *) malloc(npts*sizeof(int32))) == NULL)
                error_allo(rout,"i32_buf");
            for (n=0; n<npts; n++) {
                if (f32_buf[n] == MISSING_FLOAT) {
                    i32_buf[n] = hdf.scaled_missing;
                } else {
                    f32_buf[n] = max(hdf.range_min,min(hdf.range_max,f32_buf[n]));
                    switch (hdf.scaled_flg) {
                    case 1:
                        i32_buf[n] = nint(f32_buf[n]/scale + offset);
                        break;
                    case 2:
                        i32_buf[n] = (f32_buf[n]);
                        if (f32_buf[n] > 0) i32_buf[n] = nint(log10(f32_buf[n])/scale + offset);
                        break;
                    case 3:
                        i32_buf[n] = (f32_buf[n]);
                        if (f32_buf[n] > 0) i32_buf[n] = nint(sqrt(f32_buf[n])/scale + offset);
                        break;
                    }
                }
            }
            status = SDwritedata(*sds_id,start,NULL,edge,(void *) i32_buf);
            if (status == FAIL) {
                fprintf(stderr,"%sError writingi HDF SDS, %s - aborting\n",EXE_PROMPT,hdf.sds_name);
                exit(EXIT_FAILURE);
            }
            free(i32_buf);
            break;
        default:
            fprintf(stderr,"%sScaled data type, %d is not a 8, 16, or 32 bit signed integer - aborting",EXE_PROMPT,hdf.scaled_type);
            exit(EXIT_FAILURE);
            break;
        }

    } else {

        if (write_attr) {
            scale = 1.0;
            offset = 0.0;

            status = SDsetattr(*sds_id,"scale_factor",DFNT_FLOAT64,1,(void *) &scale);
            if (status == FAIL) {
                fprintf(stderr,"%sCannot create scale_factor attribute for SDS, %s - aborting\n",EXE_PROMPT,hdf.sds_name);
                exit(EXIT_FAILURE);
            }

            status = SDsetattr(*sds_id,"add_offset",DFNT_FLOAT64,1,(void *) &offset);
            if (status == FAIL) {
                fprintf(stderr,"%sCannot create add_offset attribute for SDS, %s - aborting\n",EXE_PROMPT,hdf.sds_name);
                exit(EXIT_FAILURE);
            }

            status = SDsetattr(*sds_id,"_FillValue",hdf.scaled_type,1,(void *) &hdf.scaled_missing);
            if (status == FAIL) {
                fprintf(stderr,"%sCannot create _FillValue attribute for SDS, %s - aborting\n",EXE_PROMPT,hdf.sds_name);
                exit(EXIT_FAILURE);
            }
        }

        status = SDwritedata(*sds_id,start,NULL,edge,var);
        if (status == FAIL) {
            fprintf(stderr,"%sError writingj HDF SDS, %s - aborting\n",EXE_PROMPT,hdf.sds_name);
            exit(EXIT_FAILURE);
        }
    }

}
