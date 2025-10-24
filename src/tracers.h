#ifndef TRACERS_H
#define TRACERS_H

#include <mpi.h>
#include "coords.h"
#include "hydro.h"
#include "ran.h"
#include "runtime.h"
#include "util_fopen.h"

#define MPI_TAG_TRAC_COUNT 55050
#define MPI_TAG_TRAC_COMM_BUF 77050


typedef enum {XY00=0, XY10, XY11, XY01, GRID_NEIGHBOUR_COUNT} XYcorner;

typedef struct tracer{
    int tracer_id;
    double intial_pos[3];
    double actual_pos[3];
    double local_pos[3];
    double tracer_u[3];
    int file_op;// if closed and 1 if open */
    int sel_for_writing; // 1 if selected for writing */

    int rel_local_coords[3];
    int local_grid_arr[GRID_NEIGHBOUR_COUNT][3];
    double tracer_u_arr[GRID_NEIGHBOUR_COUNT][3];

    
    FILE * tr_file;         
    
} trac;


typedef struct  tracers_info{

    int Ntracers; // no of tracers total in domain
    int ntracers_local; // no of tracers local to a rank
    int ntracer_capacity;

    int tracers_io_freq;
    int Ntracers_io_no;
    int tracer_seed;

    int tr_nbr[3][3][3];

    trac *tr_array;
}trs_info;

//mpi data-type
extern MPI_Datatype MPI_TRAC_TYPE;

__host__ int create_mpi_trac_datatype(MPI_Datatype * MPI_TRAC_TYPE);
__host__ int tracers_destruct(trs_info **tinfo, MPI_Datatype * MPI_TRAC_TYPE);


// sub functions
__host__  int tracer_pos_rank(cs_t * cs, double domain_pos[3]);
__host__ int tracers_parse_input(rt_t *rt, trs_info *tinfo);
__host__ int tracers_random_pos(trs_info *tinfo, pe_t *pe, cs_t *cs, double (*initial_domain_pos)[3]);
__host__  int get_velocity_at_grid(cs_t *cs, hydro_t *hydro, int local_grid_pos[3], double tracer_u_grid[3]);
__host__ int set_tr_nbr(cs_t *cs, trs_info * tinfo);
__host__  int tracer_periodic_local_pos_update(cs_t * cs, trac * tr);


//initialisation of struct
__host__ int tracers_create(cs_t *cs,  pe_t *pe, rt_t *rt, trs_info **trsinfo);

//main functions
__host__ int tracers_main( cs_t *cs, hydro_t *hydro, trs_info *tinfo);
__host__  int tracer_position_update(cs_t *cs, hydro_t *hydro, trs_info *tinfo);
__host__  int tracer_pos_euler_integ(cs_t * cs,  trac * tr);
__host__ int tracer_vel_update(cs_t *cs, hydro_t *hydro, trs_info *tinfo);
__host__ int tracer_pos_grids( trac * tr, cs_t *cs);
__host__  int tracer_grid_velocities(cs_t *cs, hydro_t *hydro, trac *tr);
__host__  int bilinear_interp_velocity(trac *tr);

//tracer particle exchange
__host__ int tracer_particle_exchange(cs_t *cs, trs_info * tinfo);
__host__ int tracer_send_recv_count(cs_t *cs, trs_info *tinfo, 
                                                    int countSend[3][3][3], int countRecv[3][3][3]);
                                    
__host__ int tracer_allocate_buffer(cs_t *cs, trac *buffSend[3][3][3], 
            trac *buffRecv[3][3][3], int countSend[3][3][3], int countRecv[3][3][3]);  
__host__ int tracer_update_trac_array_capacity(cs_t * cs, trs_info * tinfo, 
                    int countTotalSend, int countTotalRecv);
__host__ int tracer_pack_buffer(cs_t *cs, trs_info *tinfo, trac *buffSend[3][3][3]);   

__host__ int tracer_send_recv_particles(cs_t * cs, trs_info * tinfo, trac * buffSend[3][3][3], trac * buffRecv[3][3][3],
                int countSend[3][3][3], int countRecv[3][3][3], MPI_Datatype MPI_TRAC_TYPE);
__host__ int tracer_unpack_recv_buffer(cs_t *cs, trs_info * tinfo, trac *buffRecv[3][3][3], int countRecv[3][3][3]);   

//tracer write functions
__host__ int tracer_open_file(cs_t * cs, trac *tr);
__host__ int tracer_close_file(trac *tr);
__host__ int tracer_write_file(cs_t * cs,  trs_info * tinfo, int step);
__host__ int tracer_init_file(cs_t * cs, trs_info *tinfo);
__host__ int tracer_write_num_distribute(cs_t *cs, trs_info *tinfo);
__host__ int tracer_id_select_writing(cs_t *cs, trs_info *tinfo, int tracers_to_each_rank);

#endif