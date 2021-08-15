/* @f5c
**
** f5c interface
** @author: Hasindu Gamaarachchi (hasindu@unsw.edu.au)
** @@
******************************************************************************/
#ifdef  NDEBUG
    #undef NDEBUG
#endif
#include <assert.h>
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "f5c.h"
#include "error.h"

/*
todo :
Error counter for consecutive failures in the skip unreadable mode
not all the memory allocations are needed for eventalign mode
*/



core_t* init_core(opt_t opt) {
    core_t* core = (core_t*)malloc(sizeof(core_t));
    MALLOC_CHK(core);

    //model
    core->model = (model_t*)malloc(sizeof(model_t) * NUM_KMER); //4096 is 4^6 which is hardcoded now
    MALLOC_CHK(core->model);

    //load the model
    if(opt.rna){
        INFO("%s","builtin RNA nucleotide model loaded");
        uint32_t kmer_size=set_model(core->model, MODEL_ID_RNA_NUCLEOTIDE);
        assert(kmer_size==5);
    }
    else{
        uint32_t kmer_size=set_model(core->model, MODEL_ID_DNA_NUCLEOTIDE);
        assert(kmer_size==6);
    }

    core->opt = opt;

    return core;
}

void free_core(core_t* core,opt_t opt) {
    free(core->model);
    free(core);
}

db_t* init_db(core_t* core) {
    db_t* db = (db_t*)(malloc(sizeof(db_t)));
    MALLOC_CHK(db);

    db->capacity_bam_rec = core->opt.batch_size;
    db->n_bam_rec = 0;

    db->read = (char**)(malloc(sizeof(char*) * db->capacity_bam_rec));
    MALLOC_CHK(db->read);
    db->read_len = (int32_t*)(malloc(sizeof(int32_t) * db->capacity_bam_rec));
    MALLOC_CHK(db->read_len);

    db->f5 = (fast5_t**)malloc(sizeof(fast5_t*) * db->capacity_bam_rec);
    MALLOC_CHK(db->f5);

    db->et = (event_table*)malloc(sizeof(event_table) * db->capacity_bam_rec);
    MALLOC_CHK(db->et);

    db->scalings =
        (scalings_t*)malloc(sizeof(scalings_t) * db->capacity_bam_rec);
    MALLOC_CHK(db->scalings);

    db->event_align_pairs =
        (AlignedPair**)malloc(sizeof(AlignedPair*) * db->capacity_bam_rec);
    MALLOC_CHK(db->event_align_pairs);
    db->n_event_align_pairs =
        (int32_t*)malloc(sizeof(int32_t) * db->capacity_bam_rec);
    MALLOC_CHK(db->n_event_align_pairs);

    db->event_alignment = (event_alignment_t**)malloc(
        sizeof(event_alignment_t*) * db->capacity_bam_rec);
    MALLOC_CHK(db->event_alignment);
    db->n_event_alignment =
        (int32_t*)malloc(sizeof(int32_t) * db->capacity_bam_rec);
    MALLOC_CHK(db->n_event_alignment);

    db->events_per_base =
        (double*)malloc(sizeof(double) * db->capacity_bam_rec);
    MALLOC_CHK(db->events_per_base);

    db->base_to_event_map =
        (index_pair_t**)malloc(sizeof(index_pair_t*) * db->capacity_bam_rec);
    MALLOC_CHK(db->base_to_event_map);

    db->read_stat_flag = (int32_t *)malloc(sizeof(int32_t) * db->capacity_bam_rec);
    MALLOC_CHK(db->read_stat_flag);


    return db;
}


void free_db_tmp(db_t* db) {
    int32_t i = 0;
    for (i = 0; i < db->n_bam_rec; ++i) {

        free(db->read[i]);
        free(db->f5[i]->rawptr);
        free(db->f5[i]);
        free(db->et[i].event);
        free(db->event_align_pairs[i]);
        free(db->base_to_event_map[i]);

    }
}

void free_db(db_t* db) {
    //int32_t i = 0;

    free(db->read);
    free(db->read_len);
    free(db->et);
    free(db->f5);
    free(db->scalings);
    free(db->event_align_pairs);
    free(db->n_event_align_pairs);
    free(db->event_alignment);
    free(db->n_event_alignment);
    free(db->events_per_base);
    free(db->base_to_event_map);
    free(db->read_stat_flag);


    free(db);
}

void event_single(core_t* core,db_t* db, int32_t i) {

    float* rawptr = db->f5[i]->rawptr;
    float range = db->f5[i]->range;
    float digitisation = db->f5[i]->digitisation;
    float offset = db->f5[i]->offset;
    int32_t nsample = db->f5[i]->nsample;

    // convert to pA
    float raw_unit = range / digitisation;
    for (int32_t j = 0; j < nsample; j++) {
        rawptr[j] = (rawptr[j] + offset) * raw_unit;
    }

    int8_t rna=core->opt.rna;
    db->et[i] = getevents(db->f5[i]->nsample, rawptr, rna);

    // if(db->et[i].n/(float)db->read_len[i] > 20){
    //     fprintf(stderr,"%s\tevents_per_base\t%f\tread_len\t%d\n",bam_get_qname(db->bam_rec[i]), db->et[i].n/(float)db->read_len[i],db->read_len[i]);
    // }

    //get the scalings
    db->scalings[i] = estimate_scalings_using_mom(
        db->read[i], db->read_len[i], core->model, db->et[i]);

    //If sequencing RNA, reverse the events to be 3'->5'
    if (rna){
        event_t *events = db->et[i].event;
        size_t n_events = db->et[i].n;
        for (size_t i = 0; i < n_events/2; ++i) {
            event_t tmp_event = events[i];
            events[i]=events[n_events-1-i];
            events[n_events-1-i]=tmp_event;
        }
    }

}




void align_single(core_t* core, db_t* db, int32_t i) {
    db->n_event_align_pairs[i] = align(
            db->event_align_pairs[i], db->read[i], db->read_len[i], db->et[i],
            core->model, db->scalings[i], db->f5[i]->sample_rate);
        //fprintf(stderr,"readlen %d,n_events %d\n",db->read_len[i],n_event_align_pairs);
}




void process_single(core_t* core, db_t* db,int32_t i) {

    event_single(core,db,i);

    db->event_align_pairs[i] = (AlignedPair*)malloc(
        sizeof(AlignedPair) * db->et[i].n * 2); //todo : find a good heuristic to save memory //todo : save memory by freeing here itself
    MALLOC_CHK(db->event_align_pairs[i]);

    align_single(core, db,i);

    db->event_alignment[i] = NULL;
    db->n_event_alignment[i] = 0;
    db->events_per_base[i] = 0; //todo : is double needed? not just float?

    int32_t n_kmers = db->read_len[i] - KMER_SIZE + 1;
    db->base_to_event_map[i]=(index_pair_t*)(malloc(sizeof(index_pair_t) * n_kmers));
    MALLOC_CHK(db->base_to_event_map[i]);

    if (db->n_event_align_pairs[i] > 0) {
        // prepare data structures for the final calibration

        db->event_alignment[i] = (event_alignment_t*)malloc(
            sizeof(event_alignment_t) * db->n_event_align_pairs[i]);
        MALLOC_CHK(db->event_alignment[i]);

        // for (int j = 0; j < n_event_align_pairs; ++j) {
        //     fprintf(stderr, "%d-%d\n",event_align_pairs[j].ref_pos,event_align_pairs[j].read_pos);
        // }


        //todo : verify if this n is needed is needed
        db->n_event_alignment[i] = postalign(
            db->event_alignment[i],db->base_to_event_map[i], &db->events_per_base[i], db->read[i],
            n_kmers, db->event_align_pairs[i], db->n_event_align_pairs[i]);

        //fprintf(stderr,"n_event_alignment %d\n",n_events);

        // run recalibration to get the best set of scaling parameters and the residual
        // between the (scaled) event levels and the model.

        // internally this function will set shift/scale/etc of the pore model
        bool calibrated = recalibrate_model(
            core->model, db->et[i], &db->scalings[i],
            db->event_alignment[i], db->n_event_alignment[i], 1);

        // QC calibration
        if (!calibrated || db->scalings[i].var > MIN_CALIBRATION_VAR) {
            //     events[strand_idx].clear();
            free(db->event_alignment[i]);
            //free(db->event_align_pairs[i]);
            db->read_stat_flag[i] |= FAILED_CALIBRATION;
            return;
        }

        free(db->event_alignment[i]);

    } else {
        // Could not align, fail this read
        // this->events[strand_idx].clear();
        // this->events_per_base[strand_idx] = 0.0f;
        //free(db->event_align_pairs[i]);
        db->read_stat_flag[i] |= FAILED_ALIGNMENT;
        return;
    }

    // Filter poor quality reads that have too many "stays"

    if (db->events_per_base[i] > 5.0) {
        //     events[0].clear();
        //     events[1].clear();
        //free(db->event_align_pairs[i]);
        db->read_stat_flag[i] |= FAILED_QUALITY_CHK;
        return;
    }


}

void init_opt(opt_t* opt) {
    memset(opt, 0, sizeof(opt_t));
    opt->batch_size = 1;
}
