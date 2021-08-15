cdef extern from "libabea.h":

    ctypedef signed int int32_t
    ctypedef signed long int64_t
    ctypedef signed char int8_t


# struct that returns an alignment of a base-called read to the signal
# the signal interval for the base index given by base_index[i] are sample indices from raw_start_index[i] to raw_end_index[i],
# where raw_start_index[i] is included in the interval and raw_end_index[i] is not
# int8_t align_success          if the alignment was successful or not
# int32_t size_of_arrays        the number of elements in each of following arrays
# int32_t *base_index           index of the base (0-indexing) in the read
# int64_t *raw_start_index      index of the sample of the raw signal that indicates the start of the base (0-indexing).
# int64_t *raw_end_index        index of the sample of the raw signal that indicates the end of the base (0-indexing and non-inclusive).

    ctypedef struct abea_out_t:
        int8_t align_success
        int32_t size_of_arrays
        int32_t *base_index
        int64_t *raw_start_index
        int64_t *raw_end_index

# runs ABEA on a given read and returns the alignment in the space pointed by output
# NOTE: The output should be malloced and freed later by the user : see abea.pyx
    void  run_abea_on_read(abea_out_t *output,int32_t read_len, char *read, int64_t n_samples, float *samples, float digitisation, float offset, float range, float sample_rate, int8_t debug, int8_t rna)
