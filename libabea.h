#include "f5c.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

typedef struct {
    int8_t align_success;
    int32_t size_of_arrays;
    int32_t *base_index;
    int64_t *raw_start_index;
    int64_t *raw_end_index;
} abea_out_t;


void  run_abea_on_read(abea_out_t *output,int32_t read_len, char *read, int64_t n_samples, float *samples, float digitisation, float offset, float range, float sample_rate, int8_t debug);