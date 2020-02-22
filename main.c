#include "f5c.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "example.h"

void print(core_t *core, db_t *db){
 
        int32_t i = 0;
        for (i = 0; i < db->n_bam_rec; i++) {
            printf(">%s\tLN:%d\tEVENTSTART:%d\tEVENTEND:%d\t{raw_start_index,length,mean.stdv}\n",
                   "readid", (int)db->et[i].n,
                   (int)db->et[i].start, (int)db->et[i].end);
            uint32_t j = 0;
            for (j = 0; j < db->et[i].n; j++) {
                printf("{%d,%d,%f,%f}\t", (int)db->et[i].event[j].start,
                       (int)db->et[i].event[j].length, db->et[i].event[j].mean,
                       db->et[i].event[j].stdv);
            }
            printf("\n");
        }
    
        for (i = 0; i < db->n_bam_rec; i++) {
            if((db->read_stat_flag[i]) & FAILED_ALIGNMENT){
                continue;
            }
            printf(">%s\tN_ALGN_PAIR:%d\t{kmer_index,event_index}\n",
                   "readid",
                   (int)db->n_event_align_pairs[i]);
            AlignedPair* event_align_pairs = db->event_align_pairs[i];
            int32_t j = 0;
            for (j = 0; j < db->n_event_align_pairs[i]; j++) {
                printf("{%d,%d}\t", event_align_pairs[j].ref_pos,
                       event_align_pairs[j].read_pos);
            }
            printf("\n");
        }
    
        printf("read\tshift\tscale\tvar\n");

        for (i = 0; i < db->n_bam_rec; i++) {
            if((db->read_stat_flag[i])&(FAILED_ALIGNMENT|FAILED_CALIBRATION)){
                continue;
            }
            printf("%s\t%.2lf\t%.2lf\t%.2lf\n", "readid",
                   db->scalings[i].shift, db->scalings[i].scale,
                   db->scalings[i].var);
        }

        for (i = 0; i < db->n_bam_rec; i++) {
            if((db->read_stat_flag[i])&(FAILED_ALIGNMENT)){
                continue;
            }
            int32_t n_kmers = db->read_len[i] - KMER_SIZE + 1;
            printf("base_index\tstart_raw_index(inclusive)\tend_raw_index(non inclusive)\n");
            for(int j=0; j<n_kmers; j++){
                index_pair_t elem = db->base_to_event_map[i][j];
                int32_t event_start = elem.start; //eventindex start (inclusive)
                int32_t event_end = elem.stop; //eventindex end (inclusive)

                int32_t raw_start = -1 ;
                int32_t raw_end = -1;

                if(event_start!=-1){
                    raw_start = (int32_t)db->et[i].event[event_start].start; //raw signal inde start
                }
                

                if(event_end!=-1){
                    raw_end = (int32_t)db->et[i].event[event_end].start + (int32_t)db->et[i].event[event_end].length ; //raw signal end (not inclusive)
                }
                printf("%d\t%d\t%d\n",j+2,raw_start,raw_end);
            }
            
        }

    
}

int main(int argc, char* argv[]) {
  
    /************************** init - done once **********************************/
    opt_t opt;
    init_opt(&opt); //initialise options to defaults

    //initialise the core data structure
    core_t* core = init_core(opt);

    //initialise a databatch
    db_t* db = init_db(core);

    //set the batch size to 1 (for now until expansion)
    db->n_bam_rec=1;
    int i=0; //index of the read in batch - always 0 for now until expansion

	assert(READ_LEN==strlen(READ));
	assert(N_SAMPLES*sizeof(float)==sizeof(SAMPLES));


/***************** Process a single read  *********************************/
// this part can be looped for multiple reads manually, but keep i as 0 always, 

    db->read_stat_flag[i] = 0; //reset the flag

    //set the fastq
    db->read_len[i] = READ_LEN;
    db->read[i] = (char*)malloc(READ_LEN + 1); 
    strcpy(db->read[i],READ);
	
    //set the fast5
    db->f5[i] = (fast5_t*)calloc(1, sizeof(fast5_t));
    db->f5[i]->nsample= N_SAMPLES;
    db->f5[i]->rawptr = (float*)calloc(N_SAMPLES, sizeof(float));
    memcpy(db->f5[i]->rawptr,SAMPLES,N_SAMPLES*sizeof(float));
    db->f5[i]->digitisation = DIGITISATION;
    db->f5[i]->offset = OFFSET;
    db->f5[i]->range = RANGE;
    db->f5[i]->sample_rate = SAMPLE_RATE;

    //do abea
    process_single(core,db,0);

    //print the alignment
    print(core, db);

    //free temporary; free  the db->read[i], db->f5[i] as well including other internally malloced stuff
    free_db_tmp(db);


/********************************final clean up - done once*********************************/

    //free the databatch
    free_db(db);

    //free the core data structure
    free_core(core,opt);

    return 0;
}
