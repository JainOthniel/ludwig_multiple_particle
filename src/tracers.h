#ifndef TRACERS_H
#define TRACERS_H

#include <mpi.h>
#include "coords.h"
#include "hydro.h"
//temporary
#define N_TRACERS 5

typedef enum {XY00=0, XY10, XY11, XY01, GRID_NEIGHBOUR_COUNT} XYcorner;

typedef struct tracer_data trac_d;

struct tracer_data{
    double tr_u[3];
    double tr_actual_pos[3];
};


typedef struct tracer trac;

struct tracer{
    int rank_current_ts; /* current rank of the position of the tracer*/
    int rank_previous_ts;
    int tracer_id;

    double intial_pos[3];
    // trac_d * tr_d;
    double actual_pos[3];
    double local_pos[3];
    int local_grid_arr[GRID_NEIGHBOUR_COUNT][3];
    int rel_local_coords[3];

    double tracer_u[3];
    double tracer_u_arr[GRID_NEIGHBOUR_COUNT][3];

    MPI_File tr_file;
    int file_op; /* 0 if closed and 1 if open*/
    
    trac *next;// not intersted in linked list
};



typedef struct tracers_info trs_info;

struct tracers_info{

    int Ntracers; // no of tracers total in domain
    int ntracers_local; // no of tracers local to a rank

    int tracers_io_freq;
    int Ntracers_io_no;
    int tracer_seed;

    trac *tr_array;
     // array of ids of tracers
};

// sub functions
__host__  int tracer_pos_rank(cs_t * cs, double domain_pos[3]);
__host__ void tracers_parse_input(rt_t *rt, trs_info *tinfo);
__host__ void tracers_random_pos(trs_info *tinfo, pe_t *pe, cs_t *cs, double (*initial_domain_pos)[3]);


//initialisation of struct
__host__ int tracers_create(cs_t *cs,  pe_t *pe, rt_t *rt, trs_info **trsinfo);



#endif