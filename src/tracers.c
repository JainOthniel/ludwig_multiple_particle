#include <stdio.h>
#include <stdlib.h>
#include <math.h> 
#include <string.h> 
// #include <mpi.h>
#include "tracers.h"
#include "ran.h"
#include "runtime.h"



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
    //  printf("passed by rank %d", rank);
    tracers_parse_input(rt, tinfo);

    //make an array of initial positions for tracers
    double (*initial_domain_pos)[3];
    initial_domain_pos = calloc(tinfo->Ntracers, sizeof(*initial_domain_pos));
  
    // generate random positions
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
    int extra_mem = nlocal_tracer_count;
    
    // allocate memory for the tracer struct
    trac *tr;    
    tr = (trac *) calloc(tinfo->ntracers_local + extra_mem,sizeof(trac));

    // initialize/ build the tracers for each rank
    int noffset[3];
    cs_nlocal_offset(cs, noffset);
    tinfo->tr_array = tr;

    nlocal_tracer_count =0;
    for(int i=0; i<tinfo->Ntracers; i++){      
      tracer_rank = tracer_pos_rank(cs, initial_domain_pos[i]);
      if (tracer_rank != working_rank) continue;

      for(int dim=X; dim<NHDIM; dim++){
        tr[nlocal_tracer_count].intial_pos[dim] = initial_domain_pos[i][dim];
        tr[nlocal_tracer_count].actual_pos[dim] = initial_domain_pos[i][dim];
        tr[nlocal_tracer_count].tracer_u[dim] = 0.0;
        tr[nlocal_tracer_count].rel_local_coords[dim] = 0;
        tr[nlocal_tracer_count].local_pos[dim] = (int)(initial_domain_pos[i][dim] - noffset[dim]);
      }
      tr[nlocal_tracer_count].tracer_id = i+1;
      nlocal_tracer_count++;
    }
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
__host__ void tracers_parse_input(rt_t *rt, trs_info *tinfo){

  rt_int_parameter(rt, "Ntracers", &tinfo->Ntracers); //total no of tracers requested
  rt_int_parameter(rt, "tracer_io_freq", &tinfo->tracers_io_freq); // freq at which data should be written for tracerrrr
  rt_int_parameter(rt, "tracer_io_no", &tinfo->Ntracers_io_no); //total no of tracers trajec requested
  rt_int_parameter(rt, "tracer_seed", &tinfo->tracer_seed);

}

/*****************************************************************************
 * 
 * tracers_random_pos
 * generate random positions for tracers inside the domain
 * 
 *****************************************************************************/
__host__ void tracers_random_pos(trs_info *tinfo, pe_t *pe, cs_t *cs, double (*initial_domain_pos)[3]){

  int dh = 0.5;// initiate tracers at least 5 grid points from the boundary
  int ntotal[3];
  int rand_seed = tinfo->tracer_seed;// include time module to get time(NULL) to generate different random seed
  double l[3], lmin, lmax;
  ran_init_seed(pe, rand_seed);

  cs_lmin(cs, l);
  cs_ntotal(cs, ntotal);
  for (int nt =0; nt < tinfo->Ntracers; nt++){
    for(int dim = 0; dim<2; dim++){
      lmin = l[dim] + dh;
      lmax = l[dim] +  ntotal[dim] - dh;
      assert(lmax >= lmin);
      initial_domain_pos[nt][dim] = lmin + (lmax - lmin)*ran_serial_uniform();
    }
    initial_domain_pos[nt][Z] = 1.0;
  }

}