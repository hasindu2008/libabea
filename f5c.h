/* @f5c
**
** f5c interface 
** @author: Hasindu Gamaarachchi (hasindu@unsw.edu.au)
** @@
******************************************************************************/

#ifndef F5C_H
#define F5C_H

#include <stdint.h>
#include <unistd.h>
#include <stdbool.h>

/* hard coded numbers*/
#define KMER_SIZE 6 //hard coded for now; todo : change to dynamic?
#define NUM_KMER 4096   //num k-mers for 6-mers DNA
#define NUM_KMER_METH 15625 //number k-mers for 6-mers with methylated C
//#define HAVE_CUDA 1 //if compiled for CUDA or not
#define ALN_BANDWIDTH 100 // the band size in adaptive_banded_dynamic_alignment

/*flags for a read status (related to db_t->read_stat_flag)*/
#define FAILED_CALIBRATION 0x001 //if the calibration failed
#define FAILED_ALIGNMENT 0x002 //if the alignment failed
#define FAILED_QUALITY_CHK  0x004 //if the quality check failed


/* other hard coded options */


//#define ALIGN_2D_ARRAY 1 //for CPU whether to use a 1D array or a 2D array
// note : (2D arrays are very slow due to mallocs when the number of threads is high)

#define CACHED_LOG 1 //if the log values of scalings and the model k-mers are cached
//#define LOAD_SD_MEANSSTDV 1 //if the sd_mean and the sd_stdv is to be loaded (they are unused anyway)

#define MIN_CALIBRATION_VAR 2.5
#define MAX_EVENT_TO_BP_RATIO 20

//user options
typedef struct {

    int32_t batch_size;         //max reads loaded at once

} opt_t;

/* data structures */


typedef struct {
    float* rawptr;   // raw signal (float is not the best datatype type though)
    int64_t nsample; // number of samples

    //	Information for scaling raw data from ADC values to pA (are these duplicates?)
    float digitisation;
    float offset;
    float range;
    float sample_rate;

    // computed scaling paramtersd
    float scale;
    float shift;
    float drift;
    float var;
    float scale_sd;
    float var_sd;

    // derived parameters that are cached for efficiency. do we need these?
    float log_var;
    float scaled_var;
    float log_scaled_var;

} fast5_t;

// a single event : adapted from taken from scrappie
typedef struct {
    uint64_t start; //start index
    float length; //todo : cant be made int? end = start+length : end is not inclusive
    float mean;
    float stdv;
    //int32_t pos;   //todo : always -1 can be removed
    //int32_t state; //todo : always -1 can be removed
} event_t;

// event table : adapted from scrappie
typedef struct {
    size_t n;     //todo : int32_t not enough?
    size_t start; //todo : always 0?
    size_t end;   //todo : always equal to n?
    event_t* event;
} event_table;

//k-mer model
typedef struct {
    float level_mean;
    float level_stdv;

#ifdef CACHED_LOG
    //calculated for efficiency
    float level_log_stdv;
#endif

#ifdef LOAD_SD_MEANSSTDV
    //float sd_mean;
    //float sd_stdv;
    //float weight;
#endif
} model_t;

//scaling parameters for the signal : taken from nanopolish
typedef struct {
    // direct parameters that must be set
    float scale;
    float shift;
    //float drift; = 0 always?
    float var; // set later when calibrating
    //float scale_sd;
    //float var_sd;

    // derived parameters that are cached for efficiency
#ifdef CACHED_LOG
    float log_var;
#endif
    //float scaled_var;
    //float log_scaled_var;
} scalings_t;

//from nanopolish
typedef struct {
        int event_idx;
        int kmer_idx;
} EventKmerPair;

//from nanopolish
typedef struct {
    int ref_pos;
    int read_pos;
} AlignedPair;

//from nanopolish
typedef struct {
    int32_t start;
    int32_t stop; // inclusive
} index_pair_t;

//from nanopolish
typedef struct {
    // ref data
    //char* ref_name;
    char ref_kmer[KMER_SIZE + 1];
    int32_t ref_position;

    // event data
    int32_t read_idx;
    //int32_t strand_idx;
    int32_t event_idx;
    bool rc;

    // hmm data
    char model_kmer[KMER_SIZE + 1];
    char hmm_state;
} event_alignment_t;

// a data batch (dynamic data based on the reads)
typedef struct {

    int32_t capacity_bam_rec; // will these overflow?
    int32_t n_bam_rec;

    //read sequence //todo : optimise by grabbing it from bam seq. is it possible due to clipping?
    char** read;
    int32_t* read_len;

    // fast5 file //should flatten this to reduce mallocs
    fast5_t** f5;

    //event table
    event_table* et;

    //scalings
    scalings_t* scalings;

    //aligned pairs
    AlignedPair** event_align_pairs;
    int32_t* n_event_align_pairs;

    //event alignments
    event_alignment_t** event_alignment;
    int32_t* n_event_alignment;
    double* events_per_base; //todo : do we need double?

    index_pair_t** base_to_event_map;

    int32_t* read_stat_flag;


} db_t;


//core data structure (mostly static data throughout the program lifetime)
typedef struct {

    // models
    model_t* model; //dna model

    // options
    opt_t opt;


} core_t;


/* Function prototype for major functions */

db_t* init_db(core_t* core);
core_t* init_core(opt_t opt);
void align_single(core_t* core, db_t* db, int32_t i);
void free_core(core_t* core,opt_t opt);
void free_db_tmp(db_t* db);
void free_db(db_t* db);
void init_opt(opt_t* opt);
void process_single(core_t* core, db_t* db,int32_t i);

event_table getevents(size_t nsample, float* rawptr);
void set_model(model_t* model);
scalings_t estimate_scalings_using_mom(char* sequence, int32_t sequence_len,
                                       model_t* pore_model, event_table et);
int32_t align(AlignedPair* out_2, char* sequence, int32_t sequence_len,
              event_table events, model_t* models, scalings_t scaling,
              float sample_rate);
int32_t postalign(event_alignment_t* alignment, index_pair_t* base_to_event_map, double* events_per_base,
                  char* sequence, int32_t n_kmers, AlignedPair* event_alignment,
                  int32_t n_events);
bool recalibrate_model(model_t* pore_model, event_table et,
                       scalings_t* scallings,
                       const event_alignment_t* alignment_output,
                       int32_t num_alignments, bool scale_var);


#endif
