/* Defined by saranya*/

#ifndef TRACER_H
#define TRACER_H
#include <mpi.h>
// #include </usr/lib/x86_64-linux-gnu/openmpi/include/mpi.h>

#define N_TRACERS 1 // total no of tracers in simaulation
typedef enum {XY00=0, XY10, XY11, XY01, GRID_NEIGHBOUR_COUNT} XYcorner;

typedef struct tracers_info trs_in;

struct tracers_info{

    int Ntracers; // no of tracers total in domain
    int ntracers_local; // no of tracers local to a rank

    trac *head;
     // array of ids of tracers
};

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
    // double actual_pos[3];
    trac_d * tr_d;
    double local_pos[3];
    int local_grid_arr[GRID_NEIGHBOUR_COUNT][3];
    int rel_local_coords[3];

    // double tracer_u[3];
    double tracer_u_arr[GRID_NEIGHBOUR_COUNT][3];

    MPI_File tr_file;
    int file_op; /* 0 if closed and 1 if open*/
    
};


/* tracer create*/
__host__ int tracer_create(cs_t *cs, trac **ptr);
__host__   void tracer_init(trac *tr, cs_t * cs, double initial_domain_pos[3]);
__host__   int tracer_pos_rank(cs_t * cs, double domain_pos[3]);

/*tracer updation*/
__host__   int tracer_position_update(cs_t *cs, hydro_t *hydro, trac *tr);

__host__   void tracer_pos_grids( trac * tr);
__host__   void tracer_grid_velocities(cs_t *cs, hydro_t *hydro, trac *tr);
__host__   void bilinear_interp_velocity(trac *tr);
__host__   void tracer_pos_euler_integ(cs_t * cs,  trac * tr);
__host__  void tracer_update_rank(cs_t *cs, trac *tr);    

/* send recv data for tracer across ranks*/
__host__   void tracer_receive_buffer(cs_t *cs, trac *tr, MPI_Comm comm);
__host__   void tracer_send_buffer(trac *tr, MPI_Comm comm);
__host__ void tracer_send_recv(trac * tr, cs_t * cs, int working_rank);
                                                                                                             
/* write data*/ //will do later
__host__ void tracer_create_file(cs_t *cs, trac *tr, trac_d *tr_d);
__host__ void tracer_init_file(cs_t * cs, trac *tr);
__host__ void tracer_open_file(cs_t * cs, trac *tr);
__host__ void tracer_write_file(cs_t * cs, trac *tr, int step);
__host__ void tracer_close_file(trac *tr);


__host__   void tracer_periodic_local_pos_update(cs_t * cs, double local_pos[3]);
__host__   void get_velocity_at_grid(cs_t *cs, hydro_t *hydro,int local_grid_pos[3], double tracer_u_grid[3]);
                                                                                                               
#endif



