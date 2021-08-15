#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <errno.h>
#include <float.h>
#include <inttypes.h>
#include <errno.h>
#include <stdlib.h>
//#define DEBUG_MODEL_PRINT 1
#include "error.h"
#include "model.h"
#include "f5c.h"

uint32_t set_model(model_t* model, uint32_t model_id) {

    uint32_t kmer_size=0;
    uint32_t num_kmer=0;
    float *inbuilt_model=NULL;

    if(model_id==MODEL_ID_DNA_NUCLEOTIDE){
        kmer_size=6;
        num_kmer=4096;
        inbuilt_model=r9_4_450bps_nucleotide_6mer_template_model_builtin_data;
        assert(num_kmer == (uint32_t)(1 << 2*kmer_size)); //num_kmer should be 4^kmer_size
    }

    else if(model_id==MODEL_ID_RNA_NUCLEOTIDE){
        kmer_size=5;
        num_kmer=1024;
        inbuilt_model=r9_4_70bps_u_to_t_rna_5mer_template_model_builtin_data;
        assert(num_kmer == (uint32_t)(1 << 2*kmer_size)); //num_kmer should be 4^kmer_size
    }
    else{
        assert(0);
    }

    uint32_t i = 0;
    for (i = 0; i < num_kmer; i++) {
        model[i].level_mean = inbuilt_model[i * 4 + 0];
        model[i].level_stdv = inbuilt_model[i * 4 + 1];
    #ifdef LOAD_SD_MEANSSTDV
        model[i].sd_mean = inbuilt_model[i * 4 + 2];
        model[i].sd_stdv = inbuilt_model[i * 4 + 3];
    #endif
    #ifdef CACHED_LOG
        model[i].level_log_stdv=log(model[i].level_stdv);
    #endif
    }

#ifdef DEBUG_MODEL_PRINT
    i = 0;
    fprintf(stderr, "level_mean\tlevel_stdv\tsd_mean\tsd_stdv\n");
    for (i = 0; i < num_kmer; i++) {
        fprintf(stderr, "%f\t%f\t%f\t%f\n", model[i].level_mean,
                model[i].level_stdv, 0.0, 0.0);
    }
#endif

    return kmer_size;
}
