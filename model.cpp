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

//this function can be made more efficient by setting the address to the global variable
void set_model(model_t* model) {
    uint32_t i = 0;
    for (i = 0; i < NUM_KMER; i++) {
        model[i].level_mean =
            r9_4_450bps_nucleotide_6mer_template_model_builtin_data[i * 4 + 0];
        model[i].level_stdv =
            r9_4_450bps_nucleotide_6mer_template_model_builtin_data[i * 4 + 1];
    #ifdef LOAD_SD_MEANSSTDV
        model[i].sd_mean =
            r9_4_450bps_nucleotide_6mer_template_model_builtin_data[i * 4 + 2];
        model[i].sd_stdv =
            r9_4_450bps_nucleotide_6mer_template_model_builtin_data[i * 4 + 3];
    #endif
    #ifdef CACHED_LOG
        model[i].level_log_stdv=log(model[i].level_stdv);
    #endif
    }
#ifdef DEBUG_MODEL_PRINT
    i = 0;
    fprintf(stderr, "level_mean\tlevel_stdv\tsd_mean\tsd_stdv\n");
    for (i = 0; i < NUM_KMER; i++) {
        fprintf(stderr, "%f\t%f\t%f\t%f\n", model[i].level_mean,
                model[i].level_stdv, model[i].sd_mean, model[i].sd_stdv);
    }
#endif
}
