
# distutils: language = c++
# cython: language_level=3

import numpy as np
cimport numpy as np
from libc.stdlib cimport *
from libc.stdint cimport *
from libc.stdlib cimport malloc, free
from libabea cimport *

np.import_array()

def abea_python(read, samples, digitisation, offset, range, sample_rate):

    # this won't work for batching, but code in motif2sig.py can be ported
    # to support batching
    samples_floats = [float(i) for i in samples]
    cdef np.ndarray[np.float32_t,ndim=1] samples_array
    samples_array = np.ascontiguousarray(samples_floats, dtype=np.float32)

    cdef read_array = str.encode(read)
    READ_LEN = len(read)

    cdef abea_out_t *output = <abea_out_t *> malloc(sizeof(abea_out_t))
    output.base_index = <int32_t *> malloc(sizeof(int32_t)*READ_LEN)
    output.raw_start_index = <int64_t *> malloc(sizeof(int64_t)*READ_LEN)
    output.raw_end_index = <int64_t *> malloc(sizeof(int64_t)*READ_LEN)

    # code above this to create the batch of data as run_abea_on_batch with expect
    # alternatively we can have 2 function calls, one for single, one for batch

    run_abea_on_read(<abea_out_t *> output, <int32_t> READ_LEN, <char *> read_array, <int64_t> samples_array.shape[0], <float *> samples_array.data, <float> digitisation, <float> offset, <float> range, <float> sample_rate, <int8_t> 0)

    ret = {}
    if output.align_success:
        for i in xrange(output.size_of_arrays):
            ret[i] = [output.base_index[i], output.raw_start_index[i], output.raw_end_index[i]]
    else:
        ret['FAIL'] = ''

    free(output.base_index)
    free(output.raw_start_index)
    free(output.raw_end_index)
    free(output)

    return ret
