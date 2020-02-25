import numpy as np
cimport numpy as np
from libc.stdlib cimport *
from libc.stdint cimport *

np.import_array()
from abea cimport *

def abea_python(output,READ_LEN, READ, N_SAMPLES, SAMPLES, DIGITISATION, OFFSET, RANGE, SAMPLE_RATE):
    run_abea_on_read(<abea_out_t *> output, <int32_t> READ_LEN, <char *> READ, <int64_t> N_SAMPLES, <float *> SAMPLES, <float> DIGITISATION, <float> OFFSET, <float> RANGE, <float> SAMPLE_RATE, <int8_t> 0)



