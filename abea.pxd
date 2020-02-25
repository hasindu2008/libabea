cdef extern from "abea.h":

    void  run_abea_on_read(abea_out_t *output,int32_t read_len, char *read, int64_t n_samples, float *samples, float digitisation, float offset, float range, float sample_rate, int8_t debug);

