cdef extern from "libabea.h":

    
    ctypedef signed int int32_t
    ctypedef signed long int64_t
    ctypedef signed char int8_t

    ctypedef struct abea_out_t:
        int8_t align_success;
        int32_t size_of_arrays;
        int32_t *base_index;
        int64_t *raw_start_index;
        int64_t *raw_end_index;


    void  run_abea_on_read(abea_out_t *output,int32_t read_len, char *read, int64_t n_samples, float *samples, float digitisation, float offset, float range, float sample_rate, int8_t debug)

