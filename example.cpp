#include "f5c.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "example.h"
#include "libabea.h"

int main(int argc, char* argv[]) {
  
    abea_out_t *output = (abea_out_t *)malloc(sizeof(abea_out_t));
    output->base_index = (int32_t *)malloc(sizeof(int32_t)*READ_LEN);
    output->raw_start_index = (int64_t *)malloc(sizeof(int64_t)*READ_LEN);
    output->raw_end_index = (int64_t *)malloc(sizeof(int64_t)*READ_LEN);
    
    run_abea_on_read(output, READ_LEN, READ, N_SAMPLES, SAMPLES, DIGITISATION, OFFSET, RANGE, SAMPLE_RATE, 0);

    fprintf(stdout,"base_index\tstart_raw_index(inclusive)\tend_raw_index(non_inclusive)\n");
    for(int i=0; i<output->size_of_arrays; i++){
        fprintf(stdout,"%d\t%ld\t%ld\n",output->base_index[i],output->raw_start_index[i],output->raw_end_index[i]);
    }

    free(output->base_index);
    free(output->raw_start_index);
    free(output->raw_end_index);
    free(output);

    return 0;
}
