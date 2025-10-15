/* Defined by saranaya for tracer particles*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>  
#include <mpi.h>
#include "coords.h"
#include "hydro.h"
#include "tracer.h"
#include <string.h>


//////////////////////////////////////////////////////////////////

// tracer struct creater

//////////////////////////////////////////////////////////////////
__host__ int tracer_create(cs_t *cs, trac **ptr){

    trac *tr = NULL;
    trac_d *tr_d = NULL;
    assert(cs);
    assert(ptr);
    
    tr = (trac *) calloc(1, sizeof(trac));
    tr_d = (trac_d *) calloc(1, sizeof(trac_d));
    assert(tr);
    assert(tr_d);
    if (tr == NULL) pe_fatal(cs->pe, "calloc(trac) failed \n");
    if (tr_d == NULL) pe_fatal(cs->pe, "calloc(trac) failed \n");

    tr->tr_d = tr_d;

    tr->rank_current_ts = -1;
    tr->rank_previous_ts = -1;
    tr->tracer_id = -1;
    for (int dim=X; dim < NHDIM; dim++){
        tr->tr_d->tr_actual_pos[dim] = 1.0;
        tr->tr_d->tr_u[dim] = 0.0;

        tr->intial_pos[dim] = 1.0;
        tr->local_pos[dim] = 1.0;
        tr->rel_local_coords[dim] = 0;
    }

    for (XYcorner corner = XY00; corner < GRID_NEIGHBOUR_COUNT; corner++){
        for (int dim=X; dim <NHDIM; dim++){
            tr->tracer_u_arr[corner][dim] = 0.0;
            tr->local_grid_arr[corner][dim] = 1;
        }
    }
    tr->file_op = 0;
    *ptr = tr;
    return 0;
}

//////////////////////////////////////////////////////////////////

// tracer initialisation

//////////////////////////////////////////////////////////////////

/*****************************************************************************
 * 
 * tracer_init
 *
 * 
 *****************************************************************************/
__host__  void tracer_init(trac * tr, cs_t * cs, double initial_domain_pos[3]) {

    tr->rank_current_ts = tracer_pos_rank(cs,initial_domain_pos);
    tr->rank_previous_ts = tr->rank_current_ts;
    tr->tracer_id = 1; // must include somelogic to divise the id for multiple tracers.
    int rank = cs_cart_rank(cs);
    if(rank == tr->rank_current_ts){
        int noffset[3];
        cs_nlocal_offset(cs, noffset);
        int nlocal[3];
        cs_nlocal(cs, nlocal);
        for(int dim=X; dim<3; dim++){
            tr->intial_pos[dim] = initial_domain_pos[dim];
            tr->tr_d->tr_actual_pos[dim] = initial_domain_pos[dim];
            
            tr->local_pos[dim] = tr->tr_d->tr_actual_pos[dim] - noffset[dim];
            tr->rel_local_coords[dim] = (int)(tr->local_pos[dim]/nlocal[dim]);
        }
    }

}


//////////////////////////////////////////////////////////////////

// tracer_main function

//////////////////////////////////////////////////////////////////

/*****************************************************************************
 * 
 * tracer_position_update
 *
 * 
 *****************************************************************************/

__host__  int tracer_position_update(cs_t *cs, hydro_t *hydro, trac *tr) {

    int rank = cs_cart_rank(cs);
    // hydro_u_halo(hydro); // collect halo for all ranks
    if(rank == tr->rank_current_ts){
        tracer_pos_euler_integ(cs, tr);
        tracer_pos_grids(tr);
        tracer_grid_velocities(cs, hydro, tr);
        bilinear_interp_velocity(tr);
        tracer_update_rank(cs, tr);
    }

    tracer_send_recv(tr, cs, rank);
   
    return 0;
}

//////////////////////////////////////////////////////////////////

// tracer_main communication

//////////////////////////////////////////////////////////////////

/*****************************************************************************
 * 
 * tracer_rank_update
 *
 * 
 *****************************************************************************/

__host__ void tracer_update_rank(cs_t *cs, trac *tr){


    if(tr->rel_local_coords[X] !=0 || tr->rel_local_coords[Y] != 0 || tr->rel_local_coords[Z] != 0){
        int new_coords[3];
        MPI_Comm trac_comm;
        cs_cart_comm(cs, &trac_comm);

        for(int dim =X; dim<NHDIM; dim++){
            new_coords[dim] = (cs->param->mpi_cartcoords[dim] + tr->rel_local_coords[dim]) % cs->param->mpi_cartsz[dim];
        }
        // rank update
        tr->rank_previous_ts = tr->rank_current_ts;
        MPI_Cart_rank(trac_comm,new_coords,&(tr->rank_current_ts));   // update the current rank  
    }

}

/*****************************************************************************
 * 
 * tracer_send_recv for ranks
 *
 * rank argument here is the working  rank in trac_position update
 *****************************************************************************/
__host__ void tracer_send_recv(trac * tr, cs_t * cs, int working_rank){

    MPI_Bcast(&tr->rank_current_ts, 1, MPI_INT, tr->rank_previous_ts, cs->commcart);

    if (tr->rank_current_ts != tr->rank_previous_ts){
        if (working_rank == tr->rank_previous_ts ){
                
            tracer_send_buffer(tr, cs->commcart);
        }

        if(working_rank == tr->rank_current_ts){
            tracer_receive_buffer(cs, tr, cs->commcart);
            // convert the local which is out of bound to the actual local for the current rank
            tracer_periodic_local_pos_update(cs, tr->local_pos);
        }

        if(working_rank != tr->rank_current_ts){
            for (int i=0; i <3; i++){
            tr->tr_d->tr_actual_pos[i] = 0.0;
            tr->local_pos[i] = 0.0;
            tr->tr_d->tr_u[i] = 0.0;
            }
        } 
        tr->rank_previous_ts = tr->rank_current_ts; // set the previous to current as well for all the ranks   
    }
    

}
/*****************************************************************************
 * 
 * tracer_send_buffer
 *
 * 
 *****************************************************************************/
__host__  void tracer_send_buffer(trac *tr, MPI_Comm comm){

    double send_buffer[9]; // for 3 arrays of size 3 each

    for(int i=0; i<3;i++){
        send_buffer[i] = tr->tr_d->tr_actual_pos[i];
        send_buffer[i+3] = tr->local_pos[i];
        send_buffer[i+6] = tr->tr_d->tr_u[i];
    }

    // at this point the rank update has happened , data should be send from previous rank to cuurent
    MPI_Send(send_buffer, 9, MPI_DOUBLE, tr->rank_current_ts, 13, comm);
}

/*****************************************************************************
 * 
 * tracer_send_buffer
 *
 * 
 *****************************************************************************/
__host__  void tracer_receive_buffer(cs_t *cs, trac *tr, MPI_Comm comm){

    double rec_buffer[9];// for 3 arrays of size 3 each
    
    MPI_Recv(rec_buffer, 9, MPI_DOUBLE, tr->rank_previous_ts,13, comm, MPI_STATUS_IGNORE);

    for (int i=0; i <3; i++){
        tr->tr_d->tr_actual_pos[i] = rec_buffer[i];
        tr->local_pos[i] = rec_buffer[i+3];
        tr->tr_d->tr_u[i] = rec_buffer[i+6];
    }
    
}

//////////////////////////////////////////////////////////////////

// step 1 perform euler integartion

//////////////////////////////////////////////////////////////////

/*****************************************************************************
 * added by jain
 * tracer_pos_euler_integ
 *
 * perform euler integration
 *
 *****************************************************************************/
__host__  void tracer_pos_euler_integ(cs_t * cs,  trac * tr){

  
    for (int dim = X; dim < NHDIM; dim++){
        tr->local_pos[dim] += tr->tr_d->tr_u[dim]*1.0;
        tr->tr_d->tr_actual_pos[dim] += tr->tr_d->tr_u[dim]*1.0;
        // printf("\ntr_u%d = %e, pos %d = %e\n", dim, tr->tracer_u[dim], dim, tr->actual_pos[dim]);

        tr->rel_local_coords[dim] = (int)floor((tr->local_pos[dim] / cs->param->nlocal[dim])); // zero mean same rankas previous step
    }
  
  /* periodic check*/
}


/*****************************************************************************
 * added by jain
 * tracer_periodic_pos_update
 *
 * makes sure local positions always stays inside the rank
 *
 *****************************************************************************/
__host__  void tracer_periodic_local_pos_update(cs_t * cs, double local_pos[3]){

  int nlocal[3];
  cs_nlocal(cs, nlocal);
   local_pos[X] = fmod( local_pos[X], nlocal[X]);
   local_pos[X] = ( local_pos[X] <= 0.0) ?  local_pos[X] + nlocal[X] :  local_pos[X]; 

   local_pos[Y] = fmod( local_pos[Y], nlocal[Y]);
   local_pos[Y] = ( local_pos[Y] <= 0.0) ?  local_pos[Y] + nlocal[Y] :  local_pos[Y]; 
  
   local_pos[Z] = fmod( local_pos[Z], nlocal[Z]);
   local_pos[Z] = ( local_pos[Z] <= 0.0) ?  local_pos[Z] + nlocal[Z]:  local_pos[Z]; 

}

//////////////////////////////////////////////////////////////////

// step 2 get grid point velocity

//////////////////////////////////////////////////////////////////

/*****************************************************************************
 * added by jain
 * tracer_pos_grids
 *
 * gets the immediate grid points around the global position of tracer.
 *         3* * * * * * *2 
 *         *             *          3(0,1)= floor(x),floor(y) + 1   2(1,1)= floor(x) + 1,floor(y) + 1          
 *         *     ##      * ## = x,y
 *         *             *             
 *         0* * * * * * *1          0(0,0)= floor(x),floor(y)       1(1,0)= floor(x) + 1,floor(y)
 *****************************************************************************/
__host__ void tracer_pos_grids( trac * tr) {

    
    // find immediate grid neighbours
    for (XYcorner corner=XY00; corner<GRID_NEIGHBOUR_COUNT; corner++){
        tr->local_grid_arr[corner][X] = ((corner == XY00 || corner == XY01)) ? 
                                     (int)(tr->local_pos[X]) : (int)(tr->local_pos[X]) + 1;
        tr->local_grid_arr[corner][Y] = ((corner == XY00 || corner == XY10)) ? 
                                     (int)(tr->local_pos[Y]) : (int)(tr->local_pos[Y]) + 1;
        tr->local_grid_arr[corner][Z]  = (int)tr->local_pos[Z];
        
    }
    
}

/*****************************************************************************
 * added by jain
 * get_velocity_at_grid
 *
 * fetches the lattice velocity corresponding to a grid posint
 *****************************************************************************/
__host__  void get_velocity_at_grid(cs_t *cs, hydro_t *hydro, int local_grid_pos[3], double tracer_u_grid[3]) {

    int index = cs_index(cs, local_grid_pos[X], local_grid_pos[Y], local_grid_pos[Z]);
    hydro_u(hydro, index, tracer_u_grid);

    // Ensure all processes see uu (broadcast from the owning rank if needed)
    //MPI_Bcast(uu, 3, MPI_DOUBLE, 0, comm);
}

/*****************************************************************************
 * added by jain
 * tracer_grid_velocities
 *
 * makes sure global positions always stays inside the simulation box
 *
 *****************************************************************************/

__host__  void tracer_grid_velocities(cs_t *cs, hydro_t *hydro, trac *tr){
    
    for (XYcorner corner = XY00; corner < GRID_NEIGHBOUR_COUNT; corner++){
        get_velocity_at_grid(cs, hydro, tr->local_grid_arr[corner], tr->tracer_u_arr[corner]);
        // printf("\ncorner = %d, %e, %e, %e\n",corner, tr->tracer_u_arr[corner][X], tr->tracer_u_arr[corner][Y],tr->tracer_u_arr[corner][Z] );
    }

         
}

//////////////////////////////////////////////////////////////////

// step 3 do bilinear interpolation and get tracer velocity

//////////////////////////////////////////////////////////////////

/*****************************************************************************
 * 
 * bilinear_interp_velocity
 *
 * get tracer velocity - interpolated from the immediate grid point velocities
 *
 *****************************************************************************/

__host__  void bilinear_interp_velocity(trac *tr) {

    double tx, ty;
    double v_y0[3], v_y1[3];

    tx = (tr->local_pos[X] -  tr->local_grid_arr[XY00][X]) / 
                                            (tr->local_grid_arr[XY10][X] - tr->local_grid_arr[XY00][X]);
    ty = (tr->local_pos[Y] -  tr->local_grid_arr[XY00][Y]) / 
                                            (tr->local_grid_arr[XY01][Y] - tr->local_grid_arr[XY00][Y]);
    // printf("\ntx = %e, ty = %e\n", tx, ty);

    /* Interpolate along x */
    v_y0[Z] = tr->tracer_u_arr[XY00][Z];
    v_y1[Z] = tr->tracer_u_arr[XY00][Z];
    for (int i = 0; i < 2; i++) {
        v_y0[i] = (1 - tx) * tr->tracer_u_arr[XY00][i] + tx * tr->tracer_u_arr[XY10][i];
        v_y1[i] = (1 - tx) * tr->tracer_u_arr[XY01][i] + tx * tr->tracer_u_arr[XY11][i];
        // printf("\n%d : vy0 = %e, v_y1 = %e\n",i, v_y0[i], v_y1[i]);
    }

    /* Interpolate along y */
    for (int i = 0; i < 2; i++) {
        tr->tr_d->tr_u[i] = (1 - ty) * v_y0[i] + ty * v_y1[i];
    }
    // printf("\ntr_ux = %e, tr_uy = %e\n", tr->tracer_u[X], tr->tracer_u[Y]);
    tr->tr_d->tr_u[Z] = tr->tracer_u_arr[XY00][Z];
    // Ensure all processes see uu (broadcast from the owning rank if needed)
    //MPI_Bcast(f_xy, 3, MPI_DOUBLE, 0, comm);
}


//////////////////////////////////////////////////////////////////

// other sub functions

//////////////////////////////////////////////////////////////////


/*****************************************************************************
 * added by jain
 * tracer_pos_rank
 *
 *returns the rank where global position of tracer
 *
 *****************************************************************************/

__host__  int tracer_pos_rank(cs_t * cs, double domain_pos[3]){

  assert(cs);
  int nlocal[3], coords[3];
  int rank;
  
  cs_nlocal(cs, nlocal);
  MPI_Comm trac_comm;
  
  cs_cart_comm(cs, &trac_comm);

  for(int dim=X; dim<3; dim++){
    coords[dim] = (int)(domain_pos[dim]/nlocal[dim]);
  }

  MPI_Cart_rank(trac_comm, coords, &rank);
  return rank;
}

//////////////////////////////////////////////////////////////////

// Writing functions

//////////////////////////////////////////////////////////////////

/*****************************************************************************
 * added by jain
 * tracer_open_file
 *
 *opens tracer dat file
 *
 *****************************************************************************/

__host__ void tracer_open_file(cs_t * cs, trac *tr){

    char filename[64];
    snprintf(filename, sizeof(filename), "tracer_%05d.dat", tr->tracer_id);
    // tr->tr_file = fopen(filename, "a");
    // if (!tr->tr_file) {
    //     fprintf(stderr, "Error: Could not open file %s\n", filename);
    //     MPI_Abort(cs->commcart, 1);
    // }
    MPI_File_open(cs->commcart, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &tr->tr_file);
    tr->file_op = 1; 

}

/*****************************************************************************
 * added by jain
 * tracer_close_file
 *
 *closes tracer dat file
 *
 *****************************************************************************/

__host__ void tracer_close_file(trac *tr){

    MPI_File_close(&tr->tr_file);
    tr->file_op = 0;
}

/*****************************************************************************
 * added by jain
 * tracer_init_file
 *
 *initiate tracer dat file
 *
 *****************************************************************************/

__host__ void tracer_init_file(cs_t * cs, trac *tr){

    char head[64];
    int rank = cs_cart_rank(cs);
    tracer_open_file(cs, tr);
    
    if (rank == tr->rank_current_ts){        
      snprintf(head, sizeof(head), "time x y z vx vy vz\n"); 
      MPI_Offset offset;
      MPI_File_get_size(tr->tr_file, &offset);
      MPI_File_write_at(tr->tr_file, offset, head, strlen(head), MPI_CHAR, MPI_STATUS_IGNORE);
    } 
    
}

/*****************************************************************************
 * added by jain
 * tracer_init_file
 *
 *initiate tracer dat file
 *
 *****************************************************************************/

__host__ void tracer_write_file(cs_t * cs,  trac * tr, int step){

    // char filename[64];
    // sprintf(filename, sizeof(filename), "tracer_%04d.dat", tr->tracer_id);
    char line[256];
    int rank = cs_cart_rank(cs);
    // tracer_open_file(cs, tr);
    if (rank == tr->rank_current_ts){
        snprintf(line, sizeof(line), "%d %.8e %.8e %.8e %.8e %.8e %.8e\n", step,
        tr->tr_d->tr_actual_pos[X], tr->tr_d->tr_actual_pos[Y], tr->tr_d->tr_actual_pos[Z],
        tr->tr_d->tr_u[X], tr->tr_d->tr_u[Y], tr->tr_d->tr_u[Z]);
        MPI_Offset offset;
        MPI_File_get_size(tr->tr_file, &offset);
        MPI_File_write_at(tr->tr_file, offset, line,strlen(line), MPI_CHAR, MPI_STATUS_IGNORE);
    }
    // MPI_Barrier(cs->commcart);
    // tracer_close_file(tr);

}

// __host__ void create_mpi_array(cs_t *cs, pe_t *pe, trac *tr){

//     // MPI_Int
//     offsetof()
// }