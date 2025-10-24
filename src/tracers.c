#include <stdio.h>
#include <stdlib.h>
#include<stddef.h>
#include <math.h> 
#include "tracers.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

// tracers main function

//////////////////////////////////////////////////////////////////

/*****************************************************************************
 * 
 * tracers_main 
 *
 * 
 *****************************************************************************/
__host__ int tracers_main( cs_t *cs, hydro_t *hydro, trs_info *tinfo){

  // collect halo
  // update the position and get the new velocity
  tracer_position_update(cs, hydro, tinfo);

  //handle communication
  tracer_particle_exchange(cs, tinfo);

  //find current velocity
  tracer_vel_update(cs, hydro, tinfo);
  
  return 0;
}

/*****************************************************************************
 * 
 * tracers_position_update
 *
 * 
 *****************************************************************************/

 __host__  int tracer_position_update(cs_t *cs, hydro_t *hydro, trs_info *tinfo){

  int nlocal_tracers = tinfo->ntracers_local;
  trac *tr_arry = tinfo->tr_array; 
  int nlocal[3];
  cs_nlocal(cs, nlocal);
    
  for(int nt =0; nt < nlocal_tracers; nt++ ){
    // tracer_pos_euler_integ(cs, &tr_arry[ncount]);
    for(int dim = 0; dim<NHDIM; dim++){
      tr_arry[nt].actual_pos[dim] += tr_arry[nt].tracer_u[dim];
      tr_arry[nt].local_pos[dim] += tr_arry[nt].tracer_u[dim];
      tr_arry[nt].rel_local_coords[dim] = (tr_arry[nt].local_pos[dim] >= nlocal[dim]) - (tr_arry[nt].local_pos[dim] < 0);
      
      //make the relative coordinates to 0 in those directions where there is no domain split
      if(nlocal[dim] == 1){tr_arry[nt].rel_local_coords[dim] = 0;}                                   
    }
    // tracer periodic update , accesses the neighbouring ranks nlocal to update the local pos
    tracer_periodic_local_pos_update(cs, &tr_arry[nt]);

  }

  return 0;
 }

 /*****************************************************************************
 * 
 * tracers_vel_update
 *
 * 
 *****************************************************************************/
__host__ int tracer_vel_update(cs_t *cs, hydro_t *hydro, trs_info *tinfo){

  int nlocal_tracers = tinfo->ntracers_local, ncount;
  trac *tr_arry = tinfo->tr_array; 
    
  for(ncount =0; ncount < nlocal_tracers; ncount++ ){
    tracer_pos_grids(&tr_arry[ncount], cs);
    tracer_grid_velocities(cs, hydro, &tr_arry[ncount]);
    bilinear_interp_velocity(&tr_arry[ncount]);
  }

  return 0;

}

/******************************************************************************************************
 * added by jain
 * tracer_pos_grids
 *
 * gets the immediate grid points around the global position of tracer.
 *         3* * * * * * *2 
 *         *             *          3(0,1)= floor(x),floor(y) + 1   2(1,1)= floor(x) + 1,floor(y) + 1          
 *         *     ##      * ## = x,y
 *         *             *             
 *         0* * * * * * *1          0(0,0)= floor(x),floor(y)       1(1,0)= floor(x) + 1,floor(y)
 *********************************************************************************************************/

__host__ int tracer_pos_grids( trac * tr, cs_t *cs) {
  // find immediate grid neighbours
  int nlocal[3];
  cs_nlocal(cs, nlocal);
  for (XYcorner corner=XY00; corner<GRID_NEIGHBOUR_COUNT; corner++){

    tr->local_grid_arr[corner][X] = ((corner == XY00 || corner == XY01)) ? 
                                  (int)floor(tr->local_pos[X]) : (int)floor(tr->local_pos[X]) + 1;

    tr->local_grid_arr[corner][Y] = ((corner == XY00 || corner == XY10)) ? 
                                  (int)floor(tr->local_pos[Y]) : (int)floor(tr->local_pos[Y]) + 1;
    
    tr->local_grid_arr[corner][Z]  = (int)floor(tr->local_pos[Z]);      
  }    
  return 0;
}

/*****************************************************************************
 * added by jain
 * tracer_grid_velocities
 *
 * makes sure global positions always stays inside the simulation box
 *
 *****************************************************************************/

__host__  int tracer_grid_velocities(cs_t *cs, hydro_t *hydro, trac *tr){
    
  for (XYcorner corner = XY00; corner < GRID_NEIGHBOUR_COUNT; corner++){
    get_velocity_at_grid(cs, hydro, tr->local_grid_arr[corner], tr->tracer_u_arr[corner]);
  } 
  return 0;        
}

/*****************************************************************************
 * 
 * bilinear_interp_velocity
 *
 * get tracer velocity - interpolated from the immediate grid point velocities
 *
 *****************************************************************************/

__host__  int bilinear_interp_velocity(trac *tr) {

  double tx, ty;
  double v_y0[3], v_y1[3];

  tx = (tr->local_pos[X] -  tr->local_grid_arr[XY00][X]) / 
                                          (tr->local_grid_arr[XY10][X] - tr->local_grid_arr[XY00][X]);
  ty = (tr->local_pos[Y] -  tr->local_grid_arr[XY00][Y]) / 
                                          (tr->local_grid_arr[XY01][Y] - tr->local_grid_arr[XY00][Y]);

  /* Interpolate along x */
  v_y0[Z] = tr->tracer_u_arr[XY00][Z];
  v_y1[Z] = tr->tracer_u_arr[XY00][Z];
  for (int i = 0; i < 2; i++) {
    v_y0[i] = (1 - tx) * tr->tracer_u_arr[XY00][i] + tx * tr->tracer_u_arr[XY10][i];
    v_y1[i] = (1 - tx) * tr->tracer_u_arr[XY01][i] + tx * tr->tracer_u_arr[XY11][i];
  }

  /* Interpolate along y */
  for (int i = 0; i < 2; i++) {
    tr->tracer_u[i] = (1 - ty) * v_y0[i] + ty * v_y1[i];
  }
  tr->tracer_u[Z] = tr->tracer_u_arr[XY00][Z];

  return 0;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

// tracer main communication

//////////////////////////////////////////////////////////////////


/*****************************************************************************
 * 
 * tracer_particle_exchange
 *
 * 
 *****************************************************************************/

__host__ int tracer_particle_exchange(cs_t *cs, trs_info * tinfo){
  
  // int nlocal_tracers = tinfo->ntracers_local;
  int countSend[3][3][3] = {{{0}}}, countRecv[3][3][3] = {{{0}}};

  //send the count of particles leaving and recevi the count entering the rank
  //(1,1,1) will have total incoming and outgoing particles for countSend and count Recv respectively.

  tracer_send_recv_count(cs, tinfo, countSend, countRecv);

  //prepare the buffer for sending and receving
  //an array of 3*3*3 shape pointer pointing to an array of shape 3
  trac *buffSend[3][3][3], *buffRecv[3][3][3];
  memset(buffSend, 0, sizeof(buffSend));
  memset(buffRecv, 0, sizeof(buffRecv));

  //allocate the space for send and receive buffer
  tracer_allocate_buffer(cs, buffSend, buffRecv, countSend, countRecv);
  
  //extend the trac array if needed
  tracer_update_trac_array_capacity(cs, tinfo, countSend[1][1][1], countRecv[1][1][1]);
  
  //prepare the send buffer and reorder the tr_array
  tracer_pack_buffer(cs, tinfo, buffSend);

  //send and receive the buffer
  tracer_send_recv_particles(cs, tinfo, buffSend, buffRecv, countSend, countRecv, MPI_TRAC_TYPE);

  // unwrap the incomin buffer
  tracer_unpack_recv_buffer(cs, tinfo, buffRecv, countRecv);

  //free the buffer

  for(int dx = -1; dx<=1; dx++){
    for(int dy = -1; dy<=1; dy++){
      for(int dz = -1; dz<=1; dz++){

        free(buffSend[dx+1][dy+1][dz+1]);
        free(buffRecv[dx+1][dy+1][dz+1]);
      }
    }
  }

  return 0;

}

/*****************************************************************************
 * 
 * tracer_send_recv_count
 * handles particle exchange in all 26 directions 
 * 
 *****************************************************************************/

__host__ int tracer_send_recv_count(cs_t *cs, trs_info *tinfo, 
                         int countSend[3][3][3], int countRecv[3][3][3]){
            
  int nlocal_tracers = tinfo->ntracers_local;            
  // count the particles leaving through each faces.
  for(int ncount = 0; ncount < nlocal_tracers; ncount++){
    int *dxyz = tinfo->tr_array[ncount].rel_local_coords;
    //counts for everyone expcept 0,0,0 which means it didnt go past the local domain
    if (dxyz[X] || dxyz[Y] || dxyz[Z]){
      countSend[dxyz[X] + 1][dxyz[Y] + 1][dxyz[Z] + 1]++;
      countSend[1][1][1]++;// record the total sending count
    }
    
  }

  int dxo, dyo, dzo, tag;
  for(int dx = -1; dx <= 1; dx++){
    for(int dy = -1; dy <= 1; dy++){
      for(int dz = -1; dz <= 1; dz++){

        if (dx == 0 && dy == 0 && dz == 0) continue;// skip sending for (1,1,1) which is the own rank
       
        dxo = -dx; dyo = -dy; dzo = -dz;
        tag = MPI_TAG_TRAC_COUNT + ((dz + 1) + (dy + 1) *3 + (dx + 1)*3*3);

        int send_rank = tinfo->tr_nbr[dx+1][dy+1][dz+1];
        int recv_rank = tinfo->tr_nbr[dxo+1][dyo+1][dzo+1];

        if(send_rank == MPI_PROC_NULL && recv_rank == MPI_PROC_NULL) continue;
        if(recv_rank == MPI_PROC_NULL) countRecv[dxo+1][dyo+1][dzo+1] = 0;


        MPI_Sendrecv(&countSend[dx+1][dy+1][dz+1], 1, MPI_INT, send_rank, tag,
                    &countRecv[dxo+1][dyo+1][dzo+1], 1, MPI_INT, recv_rank, tag,
                     cs->commcart, MPI_STATUS_IGNORE);
      }
    }
  }

// compute total incoming particles after accumulation
countRecv[1][1][1] = 0;
for(int dx=0; dx<3; dx++){
  for(int dy=0; dy<3; dy++){
    for(int dz=0; dz<3; dz++){
      if(dx!=1 || dy!=1 || dz!=1){
        countRecv[1][1][1] += countRecv[dx][dy][dz];
      }
    }
  }
}


  return 0;
}

/*****************************************************************************
 * 
 * tracer_allocate_buffer
 * allocate memory for send and recieve buffer based on the count of incoming and oitgoing particles
 * 
 *****************************************************************************/

__host__ int tracer_allocate_buffer(cs_t *cs, trac *buffSend[3][3][3], 
            trac *buffRecv[3][3][3], int countSend[3][3][3], int countRecv[3][3][3]){

  // allocate memory for buffer
  //send buffer
  for(int dx = -1; dx <= 1; dx++){
    for(int dy = -1; dy <= 1; dy++){
      for(int dz = -1; dz <= 1; dz++){

        int nsend = countSend[dx+1][dy+1][dz+1]; 
        if (nsend <= 0 || (dx == 0 && dy == 0 && dz == 0)) continue;//skips(1,1,1)
        buffSend[dx+1][dy+1][dz+1] = calloc(nsend, sizeof(trac));

        if(!buffSend[dx+1][dy+1][dz+1]) pe_fatal(cs->pe, "Send buffer allocation failed for buffSend[%d][%d][%d]"
          , dx+1, dy+1, dz+1);

      }
    }
  }

  //Recv buffer
  for(int dx = -1; dx <= 1; dx++){
    for(int dy = -1; dy <= 1; dy++){
      for(int dz = -1; dz <= 1; dz++){

        int nrecv = countRecv[dx+1][dy+1][dz+1]; //skips(1,1,1)
        if (nrecv <= 0 || (dx == 0 && dy == 0 && dz == 0)) continue;
        buffRecv[dx+1][dy+1][dz+1] = calloc(nrecv, sizeof(trac));

        if(!buffRecv[dx+1][dy+1][dz+1]) pe_fatal(cs->pe, "Receive buffer allocation failed for buffRecv[%d][%d][%d]"
          , dx+1, dy+1, dz+1);
        
      }
    }
  }
return 0;
}

/*****************************************************************************
 * 
 * tracer_update_buffer_capacity
 *
 * 
 *****************************************************************************/
__host__ int tracer_update_trac_array_capacity(cs_t * cs, trs_info * tinfo, 
                    int countTotalSend, int countTotalRecv){

  int free_space = tinfo->ntracer_capacity - (tinfo->ntracers_local - countTotalSend);
  
  if(free_space < countTotalRecv){

    int req_space = countTotalRecv - free_space;
    int new_capacity = tinfo->ntracer_capacity + req_space;
    trac *new_arry = realloc(tinfo->tr_array, new_capacity * sizeof(trac));

    if(!new_arry){
      pe_fatal(cs->pe, "Failed intialization of %d capacity to tracer array", new_capacity);
    }

    tinfo->tr_array = new_arry;
    tinfo->ntracer_capacity = new_capacity;

  }

  return 0;  
}

/*****************************************************************************
 * 
 * tracer_pack_buffer
 *
 * 
 *****************************************************************************/
__host__ int tracer_pack_buffer(cs_t *cs, trs_info *tinfo, trac *buffSend[3][3][3]){

  int sendBuff_idx[3][3][3] = {0};
  int nlocal_tracers = tinfo->ntracers_local;
  trac *temp_arry = malloc(nlocal_tracers*sizeof(trac));
  int new_local_count = 0;

  for(int count = 0; count< nlocal_tracers; count++){
    int *dxyz = tinfo->tr_array[count].rel_local_coords;

    if(dxyz[X]|| dxyz[Y] || dxyz[Z]){
      int i = dxyz[X] + 1, j = dxyz[Y] + 1, k = dxyz[Z] + 1;
      int idx = sendBuff_idx[i][j][k];

      buffSend[i][j][k][idx] = tinfo->tr_array[count];
      sendBuff_idx[i][j][k]++;

    }else{

      temp_arry[new_local_count] = tinfo->tr_array[count];
      new_local_count++;
    }
  }

  memcpy(tinfo->tr_array, temp_arry, new_local_count*sizeof(trac));
  tinfo->ntracers_local = new_local_count;//local count after sending
  free(temp_arry);

  return 0;

}


/*****************************************************************************
 * 
 * tracer_send_recv_particles
 *
 * 
 *****************************************************************************/
__host__ int tracer_send_recv_particles(cs_t * cs, trs_info * tinfo, trac * buffSend[3][3][3], trac * buffRecv[3][3][3],
                int countSend[3][3][3], int countRecv[3][3][3], MPI_Datatype MPI_TRAC_TYPE){

  //send and receive the buffer
  MPI_Request reqs[26 * 2];
  int req_idx =0;
  for(int dx = -1; dx<=1; dx++){
    for(int dy = -1; dy<=1; dy++){
      for(int dz = -1; dz<=1; dz++){

        if (dx == 0 && dy ==0 && dz ==0)continue;
        int i = dx +1; int j = dy + 1; int k = dz +1;
        int tag = MPI_TAG_TRAC_COMM_BUF + (dz+1) + (dy+1)*3 + (dx+1)*3*3;// unique tag for all 26 directions
        int send_rank = tinfo->tr_nbr[i][j][k];
        int recv_rank = tinfo->tr_nbr[-dx+1][-dy+1][-dz+1];

        if (send_rank != MPI_PROC_NULL){
            MPI_Isend(buffSend[i][j][k],countSend[i][j][k], MPI_TRAC_TYPE, send_rank, tag, cs->commcart, &reqs[req_idx++]);
        }

        if(recv_rank != MPI_PROC_NULL){
          //2 - i = 2 -(dx +1) = -dx + 1
          MPI_Irecv(buffRecv[2-i][2-j][2-k], countRecv[2-i][2-j][2-k], MPI_TRAC_TYPE, recv_rank, tag, cs->commcart, &reqs[req_idx++]);
        }

      }
    }
  }

  MPI_Waitall(req_idx, reqs, MPI_STATUSES_IGNORE);

  return 0;
}

/*****************************************************************************
 * 
 * tracer_copy_recv_buffer
 *
 * 
 *****************************************************************************/
__host__ int tracer_unpack_recv_buffer(cs_t *cs, trs_info * tinfo, trac *buffRecv[3][3][3], int countRecv[3][3][3]){

  int recv_size =0;
  trac *temp_arr = malloc(countRecv[1][1][1]*sizeof(trac));
   for(int dx = -1; dx<=1; dx++){
    for(int dy = -1; dy<=1; dy++){
      for(int dz = -1; dz<=1; dz++){

        if(dx ==0 && dy ==0 && dz ==0) continue;
        int i = dx+1; int j = dy +1; int k = dz +1;
        memcpy(&temp_arr[recv_size], buffRecv[i][j][k],countRecv[i][j][k] * sizeof(trac));
        
        recv_size+=countRecv[i][j][k]; 
      }
    }
  }

  if (recv_size != countRecv[1][1][1]){
    pe_fatal(cs->pe, "no of recvied particle = %d , mismatched with the count of total particles expected = %d",
    recv_size, countRecv[1][1][1]);
  }

  memcpy(&tinfo->tr_array[tinfo->ntracers_local], temp_arr, recv_size*sizeof(trac));
  tinfo->ntracers_local += recv_size;
  free(temp_arr);

  return 0;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

// tracer structs creater and initialize

//////////////////////////////////////////////////////////////////

/*****************************************************************************
 * 
 * tracers_create
 *
 * 
 *****************************************************************************/
__host__ int tracers_create(cs_t *cs,  pe_t *pe, rt_t *rt, trs_info **trsinfo){

  // make trs_info struct
  trs_info *tinfo;    
  tinfo = (trs_info *) calloc(1,sizeof(trs_info));
  tracers_parse_input(rt, tinfo);

  // no of tracers trajectories requested shouldn't exceed total no of tracers
  assert(tinfo->Ntracers >= tinfo->Ntracers_io_no); 
  if (tinfo->Ntracers <= 0){
    pe_info(cs->pe, "no tracers specified - Tracer number %d", tinfo->Ntracers);
    return 0;
  }

  //make an array of initial positions for tracers
  double (*initial_domain_pos)[3];
  initial_domain_pos = calloc(tinfo->Ntracers, sizeof(*initial_domain_pos));

   // generate random positions  for tracers
  tracers_random_pos(tinfo, pe, cs, initial_domain_pos);

  int working_rank = cs_cart_rank(cs);
  int tracer_rank;
  int nlocal_tracer_count = 0;

  // find the total no of tracers local to a rank
  for(int nt=0; nt < tinfo->Ntracers; nt++){                
    tracer_rank = tracer_pos_rank(cs, initial_domain_pos[nt]);
    if (working_rank != tracer_rank) continue;
    nlocal_tracer_count++;//count the no of tracers local to a rank
  }
  
  tinfo->ntracers_local = nlocal_tracer_count;
  tinfo->ntracer_capacity = tinfo->ntracers_local + tinfo->Ntracers/5;

  //find the neighbouring ranks
  set_tr_nbr(cs, tinfo);

//////////////////////////////////////////////////////////////////////////////
//initialize and allocate memory for tracer struct
//////////////////////////////////////////////////////////////////////////////

  //define an mpi struct data type for tracer struct(for communication)
  create_mpi_trac_datatype(&MPI_TRAC_TYPE);
  
  // allocate memory for the tracer struct
  trac *tr;    
  tr = (trac *) calloc(tinfo->ntracer_capacity, sizeof(trac));

  // initialize/ build the tracers for each rank
  int noffset[3];
  cs_nlocal_offset(cs, noffset);
  tinfo->tr_array = tr;

  nlocal_tracer_count =0; //again set to zero
  for(int i=0; i<tinfo->Ntracers; i++){      
    tracer_rank = tracer_pos_rank(cs, initial_domain_pos[i]);
    if (tracer_rank != working_rank) continue;

    for(int dim=X; dim<NHDIM; dim++){
      tr[nlocal_tracer_count].intial_pos[dim] = initial_domain_pos[i][dim];
      tr[nlocal_tracer_count].actual_pos[dim] = initial_domain_pos[i][dim];
      tr[nlocal_tracer_count].local_pos[dim] = (initial_domain_pos[i][dim] - noffset[dim]);
    }
    tr[nlocal_tracer_count].tracer_id = i+1;
    nlocal_tracer_count++;
  }
  //pass the created struct to the global scope
  *trsinfo = tinfo;
  free(initial_domain_pos);

  return 0;
}

/*****************************************************************************
 * 
 * tracers_parse_input
 *
 * 
 *****************************************************************************/
__host__ int tracers_parse_input(rt_t *rt, trs_info *tinfo){

  rt_int_parameter(rt, "Ntracers", &tinfo->Ntracers); //total no of tracers requested
  rt_int_parameter(rt, "tracer_io_freq", &tinfo->tracers_io_freq); // freq at which data should be written for tracerrrr
  rt_int_parameter(rt, "tracer_io_no", &tinfo->Ntracers_io_no); //total no of tracers trajec requested
  rt_int_parameter(rt, "tracer_seed", &tinfo->tracer_seed);

  return 0;
}

/*****************************************************************************
 * 
 * tracers_random_pos
 * generate random positions for tracers inside the domain
 * 
 *****************************************************************************/
__host__ int tracers_random_pos(trs_info *tinfo, pe_t *pe, cs_t *cs, double (*initial_domain_pos)[3]){

  int dh = 0;// initiate tracers at least 5 grid points from the boundary
  double ntotal[3];
  int rand_seed = tinfo->tracer_seed;// include time module to get time(NULL) to generate different random seed
  double l[3], lmin, lmax;
  ran_init_seed(pe, rand_seed);

  cs_lmin(cs, l);
  cs_ltot(cs, ntotal);
  for (int nt =0; nt < tinfo->Ntracers; nt++){
    for(int dim = 0; dim<2; dim++){
      lmin = l[dim] + dh;
      lmax = ntotal[dim] - dh;
      assert(lmax >= lmin);
      initial_domain_pos[nt][dim] = lmin + (lmax - lmin)*ran_serial_uniform();
    }
    initial_domain_pos[nt][Z] = 1.0;
  }
  return 0;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
  int noffset[3];
  
  cs_nlocal(cs, nlocal);
  MPI_Comm trac_comm;
  
  cs_cart_comm(cs, &trac_comm);

  cs_nlocal_offset(cs, noffset);
  for(int dim=X; dim<3; dim++){
    if ((domain_pos[dim] - noffset[dim] <= 0.0) || (domain_pos[dim] > (noffset[dim] + nlocal[dim])) ) return -1;
    coords[dim] = (int)floor(domain_pos[dim]/nlocal[dim]);
  }

  MPI_Cart_rank(trac_comm, coords, &rank);
  return rank;
}

/*****************************************************************************
 * added by jain
 * get_velocity_at_grid
 *
 * fetches the lattice velocity corresponding to a grid posint
 *****************************************************************************/
__host__  int get_velocity_at_grid(cs_t *cs, hydro_t *hydro, int local_grid_pos[3], double tracer_u_grid[3]) {

  int index = cs_index(cs, local_grid_pos[X], local_grid_pos[Y], local_grid_pos[Z]);
  hydro_u(hydro, index, tracer_u_grid);

  return 0;

}

/*****************************************************************************
 * added by jain
 * tracer_periodic_pos_update
 * 
 ************************* must do the periodic update before particle exchange **** 
 * 
 * makes sure local positions always stays inside the rank
 * if the particle was in a different rank in the previous timestep
 * 
 *****************************************************************************/
__host__  int tracer_periodic_local_pos_update(cs_t * cs, trac * tr){
  double nlocal[3];

  for (int dim = X; dim < NHDIM; dim++){

    if (tr->rel_local_coords[dim] == 0) continue;
    int coords = cs->param->mpi_cartcoords[dim];
    int size = cs->param->mpi_cartsz[dim];
    
    // makes sure it wraps around
    int idm = (coords + tr->rel_local_coords[dim] + size) % size;
    nlocal[dim] = cs->listnlocal[dim][idm];
    if(nlocal[dim] == 0) continue;

    tr->local_pos[dim] = fmod(tr->local_pos[dim], nlocal[dim]) + (nlocal[dim] * (tr->local_pos[dim] < 0.0));
    
  }
   return 0;
}

/*****************************************************************************
 * added by jain
 * set_tr_nbr
 *
 * finds all the neighbours of a rank
 *****************************************************************************/

  __host__ int set_tr_nbr(cs_t *cs, trs_info * tinfo){

    int own_rank = cs_cart_rank(cs);
    int coords[3], size[3], periodic[3];
    int nbr_coords[3];

    cs_cart_coords(cs, coords);
    cs_cartsz(cs,size);
    cs_periodic(cs, periodic);

    for(int i = -1; i<=1; i++){
      for(int j = -1; j<=1; j++){
        for(int k = -1; k<=1; k++){

          // own rank goes to 1,1,1
          if (i==0 && j==0 && k==0) {tinfo->tr_nbr[i+1][j+1][k+1] = own_rank; continue;}

          // find the neibhouring ranks coordinates, wrap around if periodic else keep the values
          nbr_coords[X] = (periodic[X] == 1) ? (coords[X] + i + size[X]) % size[X] : (coords[X] + i);
          nbr_coords[Y] = (periodic[Y] == 1) ? (coords[Y] + j + size[Y]) % size[Y] : (coords[Y] + i);
          nbr_coords[Z] = (periodic[Z] == 1) ? (coords[Z] + k + size[Z]) % size[Z] : (coords[Z] + i);
          
          // checks whether any direction breaks periodicity,and puts rank to null
          if (nbr_coords[X] < 0 || nbr_coords[X] > size[X] - 1){tinfo->tr_nbr[i+1][j+1][k+1] = MPI_PROC_NULL; continue;}
          if (nbr_coords[Y] < 0 || nbr_coords[Y] > size[Y] - 1){tinfo->tr_nbr[i+1][j+1][k+1] = MPI_PROC_NULL; continue;}
          if (nbr_coords[Z] < 0 || nbr_coords[Z] > size[Z] - 1){tinfo->tr_nbr[i+1][j+1][k+1] = MPI_PROC_NULL; continue;}
                                                          
          MPI_Cart_rank(cs->commcart, nbr_coords, &tinfo->tr_nbr[i+1][j+1][k+1]);                         
      }
    }
  }
  return 0;
}

/*****************************************************************************
 * added by jain
 * create_mpi_trac_datatype
 *
 * initializes the mpi_data type for the tracer struct stricv\ctky for sending
 *****************************************************************************/
__host__ int create_mpi_trac_datatype(MPI_Datatype * MPI_TRAC_TYPE){

  const int num_items = 7;
  int blocks[7] = {1, 3, 3, 3, 3, 1,1};
  
  MPI_Datatype types[7] = {MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT, MPI_INT};
  MPI_Aint offsets[7];

  offsets[0] = offsetof(trac, tracer_id);
  offsets[1] = offsetof(trac, intial_pos);
  offsets[2] = offsetof(trac, actual_pos);
  offsets[3] = offsetof(trac, local_pos);
  offsets[4] = offsetof(trac, tracer_u);
  offsets[5] = offsetof(trac, file_op);
  offsets[6] = offsetof(trac, sel_for_writing);

  MPI_Type_create_struct(num_items, blocks, offsets, types, MPI_TRAC_TYPE);
  MPI_Type_commit(MPI_TRAC_TYPE);

  return 0;

}

/*****************************************************************************
 * added by jain
 *tracer_destruct
 *
 * initializes the mpi_data type for the tracer struct stricv\ctky for sending
 * also removes tr_Array and tinfo srtuct
 *****************************************************************************/

__host__ int tracers_destruct(trs_info **tinfo, MPI_Datatype * MPI_TRAC_TYPE){

  if(MPI_TRAC_TYPE && *MPI_TRAC_TYPE != MPI_DATATYPE_NULL){
     MPI_Type_free(MPI_TRAC_TYPE);
        *MPI_TRAC_TYPE = MPI_DATATYPE_NULL;
  }

  if (tinfo && *tinfo) {
        if ((*tinfo)->tr_array) {
            free((*tinfo)->tr_array);
            (*tinfo)->tr_array = NULL;
        }
        free(*tinfo);
        *tinfo = NULL;
  }
  return 0;
}

///////////////////////////////////////////////////////////////////////////////////

// Writing functions

///////////////////////////////////////////////////////////////////////////////////

/*****************************************************************************
 * added by jain
 * tracer_open_file
 *
 *opens tracer dat file
 *
 *****************************************************************************/

__host__ int tracer_open_file(cs_t * cs, trac *tr){

  char filename[64];
  snprintf(filename, sizeof(filename), "tracer_%07d.dat", tr->tracer_id);
  tr->tr_file = fopen(filename,"w");
  fprintf(tr->tr_file, "timestep x y z ux uy uz\n");
  tr->file_op = 1; 
  fclose(tr->tr_file);
  tr->file_op = 0;

  return 0;
}

/*****************************************************************************
 * added by jain
 * tracer_close_file
 *
 *closes tracer dat file
 *
 *****************************************************************************/

__host__ int tracer_close_file(trac *tr){

  if (tr->file_op == 1 && (tr->rel_local_coords[X] != 0 || tr->rel_local_coords[Y] != 0 
                                          || tr->rel_local_coords[Z] != 0)){
  // MPI_File_close(&tr->tr_file);
  fclose(tr->tr_file);
  tr->file_op = 0;
  }
  return 0;
}

/*****************************************************************************
 * added by jain
 * tracer_init_file
 *
 *initiate tracer dat file
 *
 *****************************************************************************/

__host__ int  tracer_init_file(cs_t * cs, trs_info *tinfo){

  // char head[64];
  // equally distribute tracers to ranks for wrtiting trajectories.
  int tracers_to_each_rank = tracer_write_num_distribute(cs, tinfo);
  
  if(tracers_to_each_rank != 0){
    
    tracer_id_select_writing(cs, tinfo, tracers_to_each_rank);
    for(int nlt =0; nlt < tinfo->ntracers_local; nlt++){

      if (tinfo->tr_array[nlt].sel_for_writing != 1){ continue;}
      tracer_open_file(cs, &tinfo->tr_array[nlt]);

    }
  }

  return 0; 
}

/*****************************************************************************
 * added by jain
 * 
 * tracer_write_file
 *
 * write tracer dat file
 *
 *****************************************************************************/

__host__ int tracer_write_file(cs_t * cs,  trs_info * tinfo, int step){

  

  for(int nlt =0; nlt< tinfo->ntracers_local; nlt++){

    if (tinfo->tr_array[nlt].sel_for_writing != 1) continue;
    char filename[256];
    snprintf(filename, sizeof(filename), "tracer_%07d.dat", tinfo->tr_array[nlt].tracer_id);
    tinfo->tr_array[nlt].tr_file = fopen(filename,"a");

    // fseek(tinfo->tr_array[nlt].tr_file, 0, SEEK_END);
    fprintf(tinfo->tr_array[nlt].tr_file,"%d %.8e %.8e %.8e %.8e %.8e %.8e\n", step,
    tinfo->tr_array[nlt].actual_pos[X], tinfo->tr_array[nlt].actual_pos[Y], tinfo->tr_array[nlt].actual_pos[Z],
    tinfo->tr_array[nlt].tracer_u[X], tinfo->tr_array[nlt].tracer_u[Y], tinfo->tr_array[nlt].tracer_u[Z]);
    
    /*if(step % 2000 == 0){
      fflush(tinfo->tr_array[nlt].tr_file);
    }*/    
    fclose(tinfo->tr_array[nlt].tr_file);

  }
  return 0;

}

/*****************************************************************************
 * added by jain
 * tracer_write_num_distribute
 *
 * calculate the required no of tracers for writing from the requested number
 * such that the requested number is equally distributed among the ranks
 *
 *****************************************************************************/
__host__ int tracer_write_num_distribute(cs_t *cs, trs_info *tinfo){

  int working_rank = cs_cart_rank(cs);
  int total_size = cs->param->mpi_cartsz[X] * cs->param->mpi_cartsz[Y] * cs->param->mpi_cartsz[Z];

  //divide the tracer trajectories requested  equally to ranks. if the requested number cannot be eqally distributed
  // then each extra will be divided to ranks starting from 0.
  int tracers_to_each_rank = (tinfo->Ntracers_io_no / total_size) + (tinfo->Ntracers_io_no % total_size > working_rank);
  
  return tracers_to_each_rank;
}

/*****************************************************************************
 * added by jain
 * tracer_id_select_writing
 *
 * selecting the tracers for writing at equal intervals from the local tracer array
 * based on the set number of trajectories for the rank.
 *****************************************************************************/

 __host__ int tracer_id_select_writing(cs_t *cs, trs_info *tinfo, int tracers_to_each_rank){

  //if required no is not present in one rank then write for all locally available tracers
  int count_space = (tinfo->ntracers_local >= tracers_to_each_rank) ? tinfo->ntracers_local / tracers_to_each_rank : 1;
  
  int count_inclu = 0;
  for(int nlt=0; nlt < tinfo->ntracers_local; nlt++){

    // select tracers for writing, at equal intervals from the local tracers array
    int idx = (nlt % count_space == 0) ? nlt : -1; 
    if (idx >= 0){
      tinfo->tr_array[idx].sel_for_writing = 1; // set the tracer for writing
      count_inclu++;
    } 

    // break after selecting required number
    if (count_inclu == tracers_to_each_rank ) {
      break;
    }
  }
  return 0;  
 }