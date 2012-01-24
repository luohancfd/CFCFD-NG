/*Copyright (C) Joseph Tang, 2007
                                                                                                        
    This file is part of OctVCE.
                                                                                                        
    OctVCE is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
                                                                                                        
    OctVCE is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
                                                                                                        
    You should have received a copy of the GNU General Public License
    along with OctVCE.  If not, see <http://www.gnu.org/licenses/>.
*/

/**\file Source file for reading/writing data to files*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "ov_kernel.h"
#include "ov_io.h"
#include "ov_lists.h"
#include "ov_adapt.h"
#include "ov_recon.h"
#include "ov_thermo.h"

#define PRINT_TEMP 1 
#define PRINT_VEL 1

#define NOTINCLUDE 0
#define PRINT_SCHL 1
#define RECON_VEL 1
#define TRY_PLOT_SHOCK 0 /*Use the adaptation parameters as a guide for shock tracking*/
#define SHOCK_PLOT_PARAM 0.085
#define SERIAL_WRITE 0

extern short int twoD;
extern double twoD_thickness;
extern Cart_cell Root;
extern short int only_get_boundary_cells_flow;
extern int num_timesteps_dump_data;
extern double time_interval_write;
extern short int num_times_refine_root;
extern short int min_refinement_level;
extern short int max_refinement_level;
extern short int min_intersect_refine_level;
extern short int IC_level;
extern short int allow_multiple_refinements;
extern short int want_to_adapt;
extern short int adapt_type;
extern double adapt_param;
extern double refine_param;
extern double coarsen_param;
extern double adapt_noise;
extern int num_timesteps_adapt;
extern short int use_geom_engine_always;
extern short int level_before_can_interrogate;
extern short int area_subcell_num;
extern short int volume_subcell_num;
extern short int wall_rep;
extern double root_shift[3];
extern double root_scale;
extern char many_limit;
extern char always_many_limit;
extern Gas_dat gas_data;
extern short int flux_type;
extern double max_schl;
extern History_loc Histories;
extern char base_hist_filename[128];
extern Jwlb_coeff JWLB_coeffs;
extern Cfl_dat CFL_data;
extern Deton_pts Deton; 
extern double global_dt;
extern int number_of_threads;
extern short int adapt_in_parallel;
extern short int set_global_dt;
extern int num_timesteps_run;
extern double finish_time;
extern char base_soln_filename[128];
extern short int output_grid;
extern short int output_grid_last;

#if NOTINCLUDE
extern double IC_bbox[2][3];
extern State_vector IC_inside_state;
extern State_vector IC_outside_state;
#else
extern State_vector IC_ambient;
extern State_vector **IC_states;
extern Body *IC_regions;
extern double ***IC_bboxes;
extern int num_ICs;
#endif

extern List_leaf Leaves;
extern List_leaf Leaves_tail;
extern List_leaf *Leaf_heads;
extern List_leaf *Leaf_tails;
extern int * num_cells_per_thread;
extern List_vtx_glob Vtxlist;
extern List_vtx_glob Vtxlist_tail;
extern IC_heat ICheat;
extern BC_data BC_est;
extern BC_data BC_wst;
extern BC_data BC_nth;
extern BC_data BC_sth;
extern BC_data BC_upr;
extern BC_data BC_lwr;

#if !NOTINCLUDE
extern int stepnow;
#endif

/*------------------------------------------------------------------*/

/**\brief Read parameters for global co-ordinate system, adaptation and geometry engine interrogation for 
    octree VCE cartesian cells*/

short int read_ov_param(char file[])
{
  FILE *fin;
  FILE *fout;
  char temp[128];
  char buffer[512];
  short int i;
  int j;
  double root_ctr[3];
  double rscale;
  short int initial_root_refine;
  short int minrefl;
  short int maxrefl;
  short int minintl;
  short int ICl;
#if NOTINCLUDE /*These variables inessential for now*/
  char multipleref;
  char usegeomalways;
#endif
  short int level_interrogate;
  char only_flow_boundary;
  short int num_sb_area;
  short int num_sb_vol;
  double dt;
  char adaptinpar;
  int total_steps;
  double param;
  int deltdump, deltadapt;

  param = adapt_param; /*This must be initialized*/

  fin = fopen(file, "r");

  if(fin == NULL)
    {
      printf("read_ov_param() could not open OctVCE paramter file \"%s\".  Exiting ...\n", file);
      return(ERROR);
    }
  
  for(i = 0; i < 7; i++)
    fgets(buffer, sizeof(buffer), fin); /*Ignore first 7 lines*/

  fgets(buffer, sizeof(buffer), fin);
  sscanf(buffer, "Centroid location: %lf %lf %lf", &root_ctr[0], &root_ctr[1], &root_ctr[2]);

  fgets(buffer, sizeof(buffer), fin);
  sscanf(buffer, "Edge length scale: %lf", &rscale); 
  if(root_scale < 0)
    {
      printf("read_ov_param(): Root edge length scale < 0.  Exiting ...\n");
      return(ERROR);
    }

  for(i = 0; i < 5; i++)
    fgets(buffer, sizeof(buffer), fin); /*Ignore next 5 lines*/

  fgets(buffer, sizeof(buffer), fin);
  strncpy(temp, "\0", 128);
  sscanf(buffer, "Type of flux solver: %s", temp);

  if(strcmp(temp, "EFM") == 0)
    flux_type = EFM;
  else if(strcmp(temp, "AUSM") == 0)
    flux_type = AUSM;
  else if(strcmp(temp, "AUSMDV") == 0)
    flux_type = AUSMDV;
  else if(strcmp(temp, "ADAPTIVE") == 0)
    flux_type = ADAPTIVE;
  else {
    printf("read_ov_param(): Unknown flux calculator - enter \"EFM\" or \"AUSM\" or \"AUSMDV\" or \"ADAPTIVE\".  Exiting ...\n");
    return(ERROR);
  }

  fgets(buffer, sizeof(buffer), fin);
  sscanf(buffer, "Use multiple limiters? %c\n", &always_many_limit);

  if((always_many_limit != 'y') && (always_many_limit != 'n')) {
    printf("read_ov_param(): For use multiple limiters, please enter 'y' or 'n'.  Exiting ...\n");
    return(ERROR);
  }

  if(always_many_limit == 'n')
    many_limit = 'n';
  else many_limit = 'y';

  fgets(buffer, sizeof(buffer), fin);
  sscanf(buffer, "Global timestep (set to 0 if don't want it fixed): %lf", &dt);

  if(dt < 0)
    {
      printf("read_ov_param(): Global timestep unknown.  Exiting ...\n");
      return(ERROR);
    }
  
  fgets(buffer, sizeof(buffer), fin);
  sscanf(buffer, "Max CFL (must be between 0 and 1): %lf", &(CFL_data.CFL_max));
  
  if((dt == 0) && (CFL_data.CFL_max < 0.0))
    {
      printf("read_ov_param(): CFL doesn't seem right.  Exiting ...\n");
      return(ERROR);
    }

  fgets(buffer, sizeof(buffer), fin);
  
  fgets(buffer, sizeof(buffer), fin);
  sscanf(buffer, "  By timestep: %d", &total_steps);

  if(total_steps < 0) {
    printf("read_ov_param(): You can't have finish timesteps < 0!  Exiting ...\n");
    return(ERROR);
  }
  else {
    fgets(buffer, sizeof(buffer), fin);
    
    if(total_steps == 0) {
      sscanf(buffer, "  By flow time: %lf", &finish_time);
      
      if(finish_time <= 0) {
	printf("read_ov_param(): You can't have a finish time <= 0!  Exiting ...\n");
	return(ERROR);
      }
    }
  }
  
  for(i = 0; i < 5; i++)
    fgets(buffer, sizeof(buffer), fin); /*Ignore next 5 lines*/

  fgets(buffer, sizeof(buffer), fin);
  sscanf(buffer, "Adapt in parallel? %c", &adaptinpar);
 
  if((adaptinpar != 'n') && (adaptinpar != 'y'))
    {
      printf("read_ov_param(): Unknown value for adapt in parallel.  Exiting ...\n");
      return(ERROR);
    } 
  
  for(i = 0; i < 5; i++) 
    fgets(buffer, sizeof(buffer), fin); /*Ignore next 5 lines*/

  fgets(buffer, sizeof(buffer), fin);
  sscanf(buffer, "Adapt every how many time steps (enter '0' for no adaptation): %d", &deltadapt);

  if(deltadapt < 0)
    {
      printf("read_ov_param(): adapt over how many time steps not right.  Exiting ...\n");
      return(ERROR);
    }

  fgets(buffer, sizeof(buffer), fin);
  sscanf(buffer, "Type of error indicator: %hd", &adapt_type);

  if((adapt_type != 1) && (adapt_type != 2) && (adapt_type != 3))
    {
      printf("read_ov_param(): Error indicator type must be 1 (using velocity differences) or 2 (density derivatives) or 3 (both) only.  Exiting ...\n");
      return(ERROR);
    }

  fgets(buffer, sizeof(buffer), fin);

  if((adapt_type == 1) || (adapt_type == 3))
    {
      fgets(buffer, sizeof(buffer), fin);
      sscanf(buffer, "  If error indicator type is 1 or 3 - compression threshold: %lf", &param);

      if(param < 0.0)
	{
	  printf("read_ov_param(): Error threshold not right.  Exiting ...\n");
	  return(ERROR);
	}

      if(adapt_type == 3) {
	for(i = 0 ; i < 3; i++)
	  fgets(buffer, sizeof(buffer), fin);
      }
      else {
	for(i = 0 ; i < 6; i++)
	  fgets(buffer, sizeof(buffer), fin);
      }
    }

  if((adapt_type == 2) || (adapt_type == 3))
    {
      if(adapt_type == 2) {
	for(i = 0 ; i < 4; i++)
	  fgets(buffer, sizeof(buffer), fin);
      }

      fgets(buffer, sizeof(buffer), fin);
      sscanf(buffer, "    Refinement threshold: %lf", &refine_param);

      if(refine_param <= 0) 
	{
	  printf("read_ov_param(): Refinement threshold <= 0.  Must be > 0.  Exiting ...\n");
	  return(ERROR);
	}

      fgets(buffer, sizeof(buffer), fin);
      sscanf(buffer, "    Coarsening threshold: %lf", &coarsen_param);

      if(coarsen_param <= 0) 
	{
	  printf("read_ov_param(): Coarsening threshold <= 0.  Must be > 0.  Exiting ...\n");
	  return(ERROR);
	}

      fgets(buffer, sizeof(buffer), fin);
      sscanf(buffer, "    Noise filter value: %lf", &adapt_noise);

      if(adapt_noise <= 0) 
	{
	  printf("read_ov_param(): Noise filter value <= 0.  Must be > 0.  Exiting ...\n");
	  return(ERROR);
	}
    }
    
  for(i = 0; i < 5; i++)
    fgets(buffer, sizeof(buffer), fin); /*Ignore next 5 lines*/

  fgets(buffer, sizeof(buffer), fin);
  sscanf(buffer, "No. times to refine root initially: %hd", &initial_root_refine);
  if(initial_root_refine < 0)
    {
      printf("read_ov_param(): No. times to refine root initially < 0.  Exiting ...\n");
      return(ERROR);
    }

  fgets(buffer, sizeof(buffer), fin);
  sscanf(buffer, "Max refinement level: %hd", &maxrefl);
  if(maxrefl < 0)
    {
      printf("read_ov_param(): Max refinement level < 0.  Exiting ...\n");
      return(ERROR);
    }

  if(twoD != 0) {
    if(deltadapt == 0) {
      printf("For 2D modes, your gap thickness should thus be %e.  Please ensure this is correct\n", 
	     rscale*pow(0.5, (double) initial_root_refine+1.0));
      twoD_thickness = rscale*pow(0.5, (double) initial_root_refine+1.0);
    }
    else {
      printf("For 2D modes, your gap thickness should thus be %e.  Please ensure this is correct\n", 
	     rscale*pow(0.5, (double) maxrefl+1.0));
      twoD_thickness = rscale*pow(0.5, (double) maxrefl+1.0);
    }
  }

  fgets(buffer, sizeof(buffer), fin);
  sscanf(buffer, "Min refinement level: %hd", &minrefl);
  if((minrefl < 0) || (minrefl > maxrefl))
    {
      printf("read_ov_param(): Min refinement level in error.  Exiting ...\n");
      return(ERROR);
    }

  fgets(buffer, sizeof(buffer), fin);
  sscanf(buffer, "Min intersect refinement level: %hd", &minintl);
  if((minintl < 0) || (minintl > maxrefl))
    {
      printf("read_ov_param(): Min intersect refinement level in error.  Exiting ...\n");
      return(ERROR);
    }
  if(minintl < minrefl)
    {
      printf("read_ov_param(): Warning.  Min intersect refinement level < min refinement level.  Adjusting so they equal ...\n");
      minintl = minrefl;
    }

  fgets(buffer, sizeof(buffer), fin);
  sscanf(buffer, "Level of cell IC intersection: %hd", &ICl);
  if(deltadapt == 0) {
    ICl = initial_root_refine;
    printf("read_ov_param(): Warning.  As you're not doing adaptive, level of IC cell intersected = no. times to refine root initially\n");
  }
  else {
    if((ICl < 0) || (ICl > maxrefl))
      {
	printf("read_ov_param(): Level of cell IC intersection in error; it's \"%hd\".  Exiting ...\n", ICl);
	return(ERROR);
      }
    if(ICl < minrefl)
      {
	printf("read_ov_param(): Warning.  Level of cell IC intersection < min refinement level.  Adjusting so they equal ...\n");
	ICl = minrefl;
      }
  }

#if NOTINCLUDE
  /*This variable not necessary for now*/

  fgets(buffer, sizeof(buffer), fin);
  sscanf(buffer, "Allow multiple refinements? %c", &multipleref);
  if((multipleref != 'y') && (multipleref != 'n'))
    {
      printf("read_ov_param(): For allow multiple refine - please enter y or n.  Exiting ...\n");
      return(ERROR);
    }
#endif
  
  for(i = 0; i < 5; i++)
    fgets(buffer, sizeof(buffer), fin); /*Ignore next 5 lines*/

#if NOTINCLUDE
  /*This variable not necessary for now*/

  fgets(buffer, sizeof(buffer), fin);
  sscanf(buffer, "Interrogate geometry engine always? %c", &usegeomalways);
  if((usegeomalways != 'y') && (usegeomalways != 'n'))
    {
      printf("read_ov_param(): For interrogating engine always - please enter y or n.  Exiting ...\n");
      return(ERROR);
    }
#endif

  fgets(buffer, sizeof(buffer), fin);
  sscanf(buffer, "Interrogate geometry engine only after what level: %hd", &level_interrogate);
  if(level_interrogate < 0)
    {
      printf("read_ov_param(): For interrogating engine after certain level - level < 0.  Exiting ...\n");
      return(ERROR);
    }
  else if(level_interrogate > minrefl) 
    {
      printf("read_ov_param(): Warning.  Minimum refinement level < minimum interrogation level.  Adjusting so min refine level = min interrogate level ...\n");
      minrefl = level_interrogate;
    }
    
#if 0 /*Is this really necessary???*/
  if(twoD != 0) {
    if((deltadapt > 0) && (level_interrogate < minrefl)) {
      printf("read_ov_param():  Warning.  For 2D adaptive modes, minimum refinement level must equal minimum interrogation level.  Adjusting so this is the case ...\n");
      level_interrogate = minrefl;
    }
  }
#endif

  fgets(buffer, sizeof(buffer), fin);
  sscanf(buffer, "No. area subcells along edge (must be even): %hd", &num_sb_area);
  if((num_sb_area < 0) || (num_sb_area%2 != 0))
    {
      printf("read_ov_param(): No. area subcells along edge in error.  Exiting ...\n");
      return(ERROR);
    }
 
  fgets(buffer, sizeof(buffer), fin);
  sscanf(buffer, "No. volume subcells along edge (must be even): %hd", &num_sb_vol);
  if((num_sb_vol < 0) || (num_sb_vol%2 != 0))
    {
      printf("read_ov_param(): No. volume subcells along edge in error.  Exiting ...\n");
      return(ERROR);
    }

  fgets(buffer, sizeof(buffer), fin);
  sscanf(buffer, "Wall representation ('0' for staircased, '1' for smoothed): %hd", &wall_rep);
  if((wall_rep != 0) && (wall_rep != 1))
    {
      printf("read_ov_param(): Wall representation must either be 0 (staircased) or 1 (smoothed).  Exiting ...\n");
      return(ERROR);
    }

  if((twoD == 2) && (wall_rep == 0)) {
    printf("read_ov_param(): Warning - in 2D axisymmetric mode, only smooth walls (representation 1) are allowed.  Changing to 1 ...\n");
    wall_rep = 1;
  }  

  for(i = 0; i < 5; i++)
    fgets(buffer, sizeof(buffer), fin); /*Ignore next 5 lines*/

  fgets(buffer, sizeof(buffer), fin);
  sscanf(buffer, "Write flow soln of only intersected cells? %c", &only_flow_boundary);
  if((only_flow_boundary != 'y') && (only_flow_boundary != 'n'))
    {
      printf("read_ov_param(): Write flow soln of only intersected cells unknown.  Exiting ...\n");
      return(ERROR);
    }

  if((twoD != 0) && (only_flow_boundary == 'y')) {
    printf("read_ov_param(): Warning.  Sorry, but cannot write only boundary cells in 2D mode.  Writing all cell info ...\n");
    return(ERROR);
  }

  fgets(buffer, sizeof(buffer), fin);

  fgets(buffer, sizeof(buffer), fin);
  sscanf(buffer, "  By timestep: %d", &deltdump);

  if(deltdump < 0)
    { 
      printf("read_ov_param(): Dump data every how many time steps not right.  Exiting ...\n");
      return(ERROR);
    }
  else {
    fgets(buffer, sizeof(buffer), fin);

    if(deltdump == 0) {
      sscanf(buffer, "  By flow time: %lf", &time_interval_write);

      if(time_interval_write <= 0) {
	printf("read_ov_param(): You can't have a solution time write interval <= 0!  Exiting ...\n");
	return(ERROR);
      }
    }
  }

  strncpy(base_soln_filename, "\0", 128);
  fgets(buffer, sizeof(buffer), fin);
  if(sscanf(buffer, "Base solution file name: %s", base_soln_filename) <= 0)
    {
      printf("read_ov_param(): Please enter a base solution filename or the simulation won't proceed.  Exiting ...\n");
      return(ERROR);
    }  

  for(i = 0; i < 5; i++)
    fgets(buffer, sizeof(buffer), fin); /*Ignore next 5 lines*/

  fgets(buffer, sizeof(buffer), fin);
  sscanf(buffer, "No. history locations: %d", &(Histories.num_locs));
  if(Histories.num_locs < 0)
    {
      printf("read_ov_param(): No. history locations < 0.  Enter 0 or more history locations.  Exiting ...\n");
      return(ERROR);
    }

  if(Histories.num_locs > 0)
    {
      fgets(buffer, sizeof(buffer), fin);
      sscanf(buffer, "Dump history locations every how many time steps: %d", &(Histories.num_timesteps_history));
      if(Histories.num_timesteps_history <= 0)
	{
	  printf("read_ov_param(): No. timesteps history files are dumped <= 0.  Exiting ...\n");
	  return(ERROR);
	}

      strncpy(base_hist_filename, "\0", 128);
      fgets(buffer, sizeof(buffer), fin);
      if(sscanf(buffer, "Base history file name: %s", base_hist_filename) <= 0) 
	{
	  printf("read_ov_param(): Please enter a base history filename or the simulation won't proceed.  Exiting ...\n");
	  return(ERROR);
	}

      for(i = 0; i < 3; i++)
	fgets(buffer, sizeof(buffer), fin); /*Ignore next 3 lines*/

      Histories.locs = malloc(sizeof(double *)*(Histories.num_locs));
      if(Histories.locs == NULL) {
	printf("read_ov_param(): Error allocating history files data structure.  Exiting ...\n");
	return(ERROR);
      }
      
      for(j = 0; j < Histories.num_locs; j++) {
	Histories.locs[j] = malloc(3*sizeof(double));
	if(Histories.locs[j] == NULL) {
	  printf("read_ov_param(): Error allocating history files data structure.  Exiting ...\n");
	  return(ERROR);
	}
			
	fgets(buffer, sizeof(buffer), fin);
	if(sscanf(buffer, "%lf %lf %lf", &(Histories.locs[j][0]), &(Histories.locs[j][1]), &(Histories.locs[j][2])) != 3) {
	  printf("read_ov_param(): History location %d seems to be in error.  Exiting ...\n", j);
	  return(ERROR);
	}

	strncpy(temp, "\0", 128);
	sprintf(temp, "%s%d.dat", base_hist_filename, j);

	fout = fopen(temp, "w");
	fprintf(fout, "# Pressure history file %hd for location (%e %e %e) - 1st col time, 2nd col overpressure\n", j, Histories.locs[j][0], Histories.locs[j][1], Histories.locs[j][2]);
	fclose(fout);
      }
    }

  fclose(fin);

  /*All parameters input correctly*/

  for(i = 0; i < 3; i++)
    root_shift[i] = root_ctr[i];
  root_scale = rscale;
  num_times_refine_root = initial_root_refine;
  min_refinement_level = minrefl;
  max_refinement_level = maxrefl;
  min_intersect_refine_level = minintl;
  IC_level = ICl;

#if NOTINCLUDE
  if(multipleref == 'y')
    allow_multiple_refinements = TRUE;
  else allow_multiple_refinements = FALSE;

  if(usegeomalways == 'y')
    use_geom_engine_always = TRUE;
  else use_geom_engine_always = FALSE;
#endif

  level_before_can_interrogate = level_interrogate;

  area_subcell_num = num_sb_area;
  volume_subcell_num = num_sb_vol;
  global_dt = dt;  
  num_timesteps_run = total_steps;

  if(dt > 0)
    set_global_dt = TRUE;
  
  if(number_of_threads == 1) 
    adapt_in_parallel = FALSE;
  else
    {
      if(adaptinpar == 'y')
	adapt_in_parallel = TRUE;
      else adapt_in_parallel = FALSE;
    } 
  
  if(deltadapt == 0)
    want_to_adapt = FALSE;
  else want_to_adapt = TRUE;

  adapt_param = param;
  num_timesteps_adapt = deltadapt;
  
  if(only_flow_boundary == 'y')
    only_get_boundary_cells_flow = TRUE;
  else only_get_boundary_cells_flow = FALSE;

  num_timesteps_dump_data = deltdump;

  return(NOERROR);
}

/*------------------------------------------------------------------*/

/**\brief Read in gas properties*/

short int read_ov_gas(char file[])
{
  FILE *fin;
  char buffer[512], temp[32];
  short int i;

  fin = fopen(file, "r");
  if(fin == NULL) {
    printf("Read_ov_gas(): Couldn't open gas properties file \"%s\".  Exiting ...\n", file);
    return(ERROR);
  }

  for(i = 0; i < 3; i++)
    fgets(buffer, sizeof(buffer), fin); /*Ignore first 3 lines*/

  fgets(buffer, sizeof(buffer), fin);
  sscanf(buffer, "No. species: %d", &(gas_data.num_species));
  if((gas_data.num_species <= 0) || (gas_data.num_species > 2)) {
    printf("Read_ov_gas(): Can only support 1 or 2 species at the moment.  Exiting ...\n");
    return(ERROR);
  }

  fgets(buffer, sizeof(buffer), fin);
  sscanf(buffer, "Ambient gas Cv: %lf", &(gas_data.Cv_amb));
  if(gas_data.Cv_amb <= 0) {
    printf("Read_ov_gas(): Cv of ambient gas <= 0.  Exiting ...\n");
    return(ERROR);
  }

  fgets(buffer, sizeof(buffer), fin);
  sscanf(buffer, "Ambient gas gamma: %lf", &(gas_data.gam_amb));
  if(gas_data.gam_amb <= 0) {
    printf("Read_ov_gas(): Gamma of ambient gas <= 0.  Exiting ...\n");
    return(ERROR);
  }

  gas_data.Cv_prod = gas_data.Cv_amb;
  gas_data.gam_prod = gas_data.gam_amb; /*Initialize these values*/

  if(gas_data.num_species == 1) {
    gas_data.eos_prod = PERFECT_GAS; /*By default*/
  }
  else {
    for(i = 0; i < 3; i++)
      fgets(buffer, sizeof(buffer), fin); /*Ignore next 3 lines*/

    fgets(buffer, sizeof(buffer), fin);
    sscanf(buffer, "  EOS type for products: %s", temp);

    if(strncmp(temp, "JWLB", 4) == 0)
      gas_data.eos_prod = JWLB;
    else if(strncmp(temp, "Ideal", 5) == 0)
      gas_data.eos_prod = PERFECT_GAS;
    else {
      printf("Read_ov_gas(): Don't know what equation of state \"%s\" is (only 'Ideal' or 'JWLB' allowed).  Exiting ...\n", temp);
      return(ERROR);
    }

    fgets(buffer, sizeof(buffer), fin);
    sscanf(buffer, "  Products gas Cv: %lf", &(gas_data.Cv_prod));
    if(gas_data.Cv_prod <= 0) {
      printf("Read_ov_gas(): Cv of explosion products gas <= 0.  Exiting ...\n");
      return(ERROR);
    }

    if(gas_data.eos_prod == PERFECT_GAS) {
      for(i = 0; i < 3; i++)
	fgets(buffer, sizeof(buffer), fin);

      fgets(buffer, sizeof(buffer), fin);
      sscanf(buffer, "    Products gas gamma: %lf", &(gas_data.gam_prod));
      if(gas_data.gam_prod <= 0) {
	printf("Read_ov_gas(): Gamma of explosion products gas <= 0.  Exiting ...\n");
	return(ERROR);
      }
    }
    else if(gas_data.eos_prod == JWLB) {
      for(i = 0; i < 5; i++)
	fgets(buffer, sizeof(buffer), fin); /*Ignore next 5 lines*/

      fgets(buffer, sizeof(buffer), fin);
      sscanf(buffer, "    If EOS is JWLB - want to treat low products density as air? %c", &(gas_data.jwlb_cutoff));

      if((gas_data.jwlb_cutoff != 'y') && (gas_data.jwlb_cutoff != 'n')) {
	printf("Read_ov_gas(): Put 'y' or 'n' only if you want to treat low products density as air.  Exiting ...\n");
	return(ERROR);
      }
      
      for(i = 0; i < 3; i++)
	fgets(buffer, sizeof(buffer), fin); /*Ignore next 3 lines*/

      fgets(buffer, sizeof(buffer), fin);
      if(sscanf(buffer, "    A1: %lf", &(JWLB_coeffs.A[0])) <= 0) {
	printf("Read_ov_gas(): Can't read in JWLB A1.  Exiting ...\n");
	return(ERROR);
      }

      fgets(buffer, sizeof(buffer), fin);
      if(sscanf(buffer, "    A2: %lf", &(JWLB_coeffs.A[1])) <= 0) {
	printf("Read_ov_gas(): Can't read in JWLB A2.  Exiting ...\n");
	return(ERROR);
      }

      fgets(buffer, sizeof(buffer), fin);
      if(sscanf(buffer, "    A3: %lf", &(JWLB_coeffs.A[2])) <= 0) {
	printf("Read_ov_gas(): Can't read in JWLB A3.  Exiting ...\n");
	return(ERROR);
      }
            
      fgets(buffer, sizeof(buffer), fin);
      if(sscanf(buffer, "    A4: %lf", &(JWLB_coeffs.A[3])) <= 0) {
	printf("Read_ov_gas(): Can't read in JWLB A4.  Exiting ...\n");
	return(ERROR);
      }

      fgets(buffer, sizeof(buffer), fin);
      if(sscanf(buffer, "    A5: %lf", &(JWLB_coeffs.A[4])) <= 0) {
	printf("Read_ov_gas(): Can't read in JWLB A5.  Exiting ...\n");
	return(ERROR);
      }

      fgets(buffer, sizeof(buffer), fin);
      if(sscanf(buffer, "    R1: %lf", &(JWLB_coeffs.R[0])) <= 0) {
	printf("Read_ov_gas(): Can't read in JWLB R1.  Exiting ...\n");
	return(ERROR);
      }

      fgets(buffer, sizeof(buffer), fin);
      if(sscanf(buffer, "    R2: %lf", &(JWLB_coeffs.R[1])) <= 0) {
	printf("Read_ov_gas(): Can't read in JWLB R2.  Exiting ...\n");
	return(ERROR);
      }

      fgets(buffer, sizeof(buffer), fin);
      if(sscanf(buffer, "    R3: %lf", &(JWLB_coeffs.R[2])) <= 0) {
	printf("Read_ov_gas(): Can't read in JWLB R3.  Exiting ...\n");
	return(ERROR);
      }

      fgets(buffer, sizeof(buffer), fin);
      if(sscanf(buffer, "    R4: %lf", &(JWLB_coeffs.R[3])) <= 0) {
	printf("Read_ov_gas(): Can't read in JWLB R4.  Exiting ...\n");
	return(ERROR);
      }

      fgets(buffer, sizeof(buffer), fin);
      if(sscanf(buffer, "    R5: %lf", &(JWLB_coeffs.R[4])) <= 0) {
	printf("Read_ov_gas(): Can't read in JWLB R5.  Exiting ...\n");
	return(ERROR);
      }

      fgets(buffer, sizeof(buffer), fin);
      if(sscanf(buffer, "    AL1: %lf", &(JWLB_coeffs.AL[0])) <= 0) {
	printf("Read_ov_gas(): Can't read in JWLB AL1.  Exiting ...\n");
	return(ERROR);
      }

      fgets(buffer, sizeof(buffer), fin);
      if(sscanf(buffer, "    AL2: %lf", &(JWLB_coeffs.AL[1])) <= 0) {
	printf("Read_ov_gas(): Can't read in JWLB AL2.  Exiting ...\n");
	return(ERROR);
      }

      fgets(buffer, sizeof(buffer), fin);
      if(sscanf(buffer, "    AL3: %lf", &(JWLB_coeffs.AL[2])) <= 0) {
	printf("Read_ov_gas(): Can't read in JWLB AL3.  Exiting ...\n");
	return(ERROR);
      }

      fgets(buffer, sizeof(buffer), fin);
      if(sscanf(buffer, "    AL4: %lf", &(JWLB_coeffs.AL[3])) <= 0) {
	printf("Read_ov_gas(): Can't read in JWLB AL4.  Exiting ...\n");
	return(ERROR);
      }

      fgets(buffer, sizeof(buffer), fin);
      if(sscanf(buffer, "    AL5: %lf", &(JWLB_coeffs.AL[4])) <= 0) {
	printf("Read_ov_gas(): Can't read in JWLB AL5.  Exiting ...\n");
	return(ERROR);
      }

      fgets(buffer, sizeof(buffer), fin);
      if(sscanf(buffer, "    BL1: %lf", &(JWLB_coeffs.BL[0])) <= 0) {
	printf("Read_ov_gas(): Can't read in JWLB BL1.  Exiting ...\n");
	return(ERROR);
      }

      fgets(buffer, sizeof(buffer), fin);
      if(sscanf(buffer, "    BL2: %lf", &(JWLB_coeffs.BL[1])) <= 0) {
	printf("Read_ov_gas(): Can't read in JWLB BL2.  Exiting ...\n");
	return(ERROR);
      }

      fgets(buffer, sizeof(buffer), fin);
      if(sscanf(buffer, "    BL3: %lf", &(JWLB_coeffs.BL[2])) <= 0) {
	printf("Read_ov_gas(): Can't read in JWLB BL3.  Exiting ...\n");
	return(ERROR);
      }

      fgets(buffer, sizeof(buffer), fin);
      if(sscanf(buffer, "    BL4: %lf", &(JWLB_coeffs.BL[3])) <= 0) {
	printf("Read_ov_gas(): Can't read in JWLB BL4.  Exiting ...\n");
	return(ERROR);
      }

      fgets(buffer, sizeof(buffer), fin);
      if(sscanf(buffer, "    BL5: %lf", &(JWLB_coeffs.BL[4])) <= 0) {
	printf("Read_ov_gas(): Can't read in JWLB BL5.  Exiting ...\n");
	return(ERROR);
      }

      fgets(buffer, sizeof(buffer), fin);
      if(sscanf(buffer, "    RL1: %lf", &(JWLB_coeffs.RL[0])) <= 0) {
	printf("Read_ov_gas(): Can't read in JWLB RL1.  Exiting ...\n");
	return(ERROR);
      }

      fgets(buffer, sizeof(buffer), fin);
      if(sscanf(buffer, "    RL2: %lf", &(JWLB_coeffs.RL[1])) <= 0) {
	printf("Read_ov_gas(): Can't read in JWLB RL2.  Exiting ...\n");
	return(ERROR);
      }

      fgets(buffer, sizeof(buffer), fin);
      if(sscanf(buffer, "    RL3: %lf", &(JWLB_coeffs.RL[2])) <= 0) {
	printf("Read_ov_gas(): Can't read in JWLB RL3.  Exiting ...\n");
	return(ERROR);
      }

      fgets(buffer, sizeof(buffer), fin);
      if(sscanf(buffer, "    RL4: %lf", &(JWLB_coeffs.RL[3])) <= 0) {
	printf("Read_ov_gas(): Can't read in JWLB RL4.  Exiting ...\n");
	return(ERROR);
      }

      fgets(buffer, sizeof(buffer), fin);
      if(sscanf(buffer, "    RL5: %lf", &(JWLB_coeffs.RL[4])) <= 0) {
	printf("Read_ov_gas(): Can't read in JWLB RL5.  Exiting ...\n");
	return(ERROR);
      }

      fgets(buffer, sizeof(buffer), fin);
      if(sscanf(buffer, "    C: %lf", &(JWLB_coeffs.C)) <= 0) {
	printf("Read_ov_gas(): Can't read in JWLB C.  Exiting ...\n");
	return(ERROR);
      }

      fgets(buffer, sizeof(buffer), fin);
      if(sscanf(buffer, "    OMEGA: %lf", &(JWLB_coeffs.omega)) <= 0) {
	printf("Read_ov_gas(): Can't read in JWLB omega.  Exiting ...\n");
	return(ERROR);
      }

      fgets(buffer, sizeof(buffer), fin);
      if(sscanf(buffer, "    v0: %lf", &(JWLB_coeffs.v0)) <= 0) {
	printf("Read_ov_gas(): Can't read in JWLB v0.  Exiting ...\n");
	return(ERROR);
      }
      else if((JWLB_coeffs.v0) <= 0) {
	printf("Read_ov_gas(): v0 <= 0.  Exiting ...\n");
	return(ERROR);
      }

      if(JWLB_coeffs.RL[0] > 0)
	precompute_JRL();
      get_smallest_JWLB_terms();
    }
  }

  return(NOERROR);
}

/*------------------------------------------------------------------*/

/**\brief Read in special initial condition requirements i.e. where specific regions have high
   temperature, pressure etc.*/

short int read_IC(char file[])
{
  FILE *fin;
  char temp1[128], temp2[128], temp3[128];
  int str_length = 128;
  int i, j;
  int c = 'A';
  int count_gas;
  int count_vel;
  int region_seen;  
  int ambient_now, specific_now, dealt_ambient = FALSE;
  State_vector *S;    
  double rho, p, q, T, U, V, W, R, v, term1, term3;

  double expRv[5] = {0};
  double expRLv[5] = {0};
  double expRL[5] = {0};

  IC_ambient.Rhof = 0;
  IC_ambient.Pres = 0;
  IC_ambient.RhoE = 0;
  U = 0; V = 0; W = 0;
  fin = fopen(file, "r");
  if(fin == NULL) {
    printf("Read_IC(): Couldn't open initial condition file \"%s\".  Exiting ...\n", file);
    return(ERROR);
  }

  while(c != EOF) { /*Keep going through the file until EOF returned*/
    c = getc(fin); 
    
    if(isalpha(c)) {
      i = 0;
      strncpy(temp1, "\0", str_length);
      
      while(isalpha(c)) {
	temp1[i] = (char) c;
	i++;
	c = getc(fin);

	if(c == EOF) {
	  printf("read_IC(): Inappropriate end of file.  Exiting ...\n"); 
	  return(ERROR);
	}

	if(i > str_length) {
	  printf("read_IC(): String overflow.  Exiting ...\n");
	  return(ERROR);
	}
      }

      if((strncmp(temp1, "Ambient", strlen(temp1)) == 0) || (strncmp(temp1, "Products", strlen(temp1)) == 0)) {

	/*We're now supposedly entering a section which describes some sort of IC*/

	while(isspace(c))
	  c = getc(fin);

	if(isalpha(c)) {
	  i = 0;
	  strncpy(temp2, "\0", str_length);
	  while(isalpha(c)) {
	    temp2[i] = (char) c;
	    i++;
	    c = getc(fin);

	    if(c == EOF) {
	      printf("read_IC(): Inappropriate end of file.  Exiting ...\n"); 
	      return(ERROR);
	    }

	    if(i > str_length) {
	      printf("read_IC(): String overflow.  Exiting ...\n");
	      return(ERROR);
	    }   
	  }
	}

	if(strncmp(temp2, "conditions", strlen(temp2)) == 0) { /*Now read in the next 5 values*/
	  
	  ambient_now = FALSE; 
	  specific_now = FALSE;
	  S = NULL;
	  count_gas = 0;
	  count_vel = 0;
	  rho = 0;
	  p = 0;
	  T = 0;

	  if(strncmp(temp1, "Ambient", strlen(temp1)) == 0) {
	    ambient_now = TRUE;
	    region_seen = TRUE;
	  }
	  else {
	    specific_now = TRUE;
	    region_seen = FALSE;
	  }
	  
	  if(specific_now == TRUE)
	    S = malloc(sizeof(struct state_vector));
	  
	  while((count_gas < 2) || (count_vel < 3) || (region_seen == FALSE)) {
	    
	    while((c != EOF) && !isalpha(c)) 
	      c = getc(fin); 	    
	    
	    if(c == EOF) {
	      printf("read_IC(): Inappropriate end of file.  Exiting ...\n"); 
	      return(ERROR);
	    }
	  
	    i = 0;
	    strncpy(temp2, "\0", str_length);
	    while(isalpha(c)) {
	      temp2[i] = (char) c;
	      i++;
	      c = getc(fin);
	      if(c == EOF) {
		printf("read_IC(): Inappropriate end of file.  Exiting ...\n");
		return(ERROR);
	      }

	      if(i > str_length) {
		printf("read_IC(): String overflow.  Exiting ...\n");
		return(ERROR);
	      }
	    } 
	    
	    if((strncmp(temp2, "Density", 7) == 0) && (count_gas < 2)) { 
	      while(c != EOF) {
		if(isdigit(c) || (c == '.'))
		  break;
		c = getc(fin); 
	      }
	      
	      i = 0;
	      strncpy(temp2, "\0", str_length);	      
	      while(!isspace(c)) {
		temp2[i] = (char) c;
		i++;
		if(i > str_length) {
		  printf("read_IC(): String overflow.  Exiting ...\n");
		  return(ERROR);
		}
		c = getc(fin);
	      } 
	      
	      rho = strtod(temp2, NULL); 
	      if(rho == 0) {
		printf("read_IC(): Non-physical density value of %s was seen.  Exiting ...\n", temp2);
		return(ERROR);
	      }
	      
	      count_gas++;
	    }
	    else if((strncmp(temp2, "Pressure", 8) == 0) && (count_gas < 2)) {
	      while(c != EOF) {
		if(isdigit(c) || (c == '.'))
		  break;
		c = getc(fin); 
	      }

	      i = 0;
	      strncpy(temp2, "\0", str_length);	      
	      while(!isspace(c)) {
		temp2[i] = (char) c;
		i++;
		c = getc(fin);
		if(i > str_length) {
		  printf("read_IC(): String overflow.  Exiting ...\n");
		  return(ERROR);
		}
	      } 

	      p = strtod(temp2, NULL); 
	      if(p == 0) {
		printf("read_IC(): Non-physical pressure value of %s was seen.  Exiting ...\n", temp2);
		return(ERROR);
	      }

	      count_gas++;
	    }
	    else if((strncmp(temp2, "Temperature", 11) == 0) && (count_gas < 2)) {
	      while(c != EOF) {
		if(isdigit(c) || (c == '.'))
		  break;
		c = getc(fin); 
	      }

	      i = 0;
	      strncpy(temp2, "\0", str_length);	      
	      while(!isspace(c)) {
		temp2[i] = (char) c;
		i++;
		c = getc(fin);
		if(i > str_length) {
		  printf("read_IC(): String overflow.  Exiting ...\n");
		  return(ERROR);
		}
	      }

	      T = strtod(temp2, NULL); 
	      if(T == 0) {
		printf("read_IC(): Non-physical temperature value of %s was seen.  Exiting ...\n", temp2);
		return(ERROR);
	      }

	      count_gas++;
	    }
	    else if(strncmp(temp2, "U", 1) == 0) { 
	      while(c != EOF) {
		if(isdigit(c) || (c == '.') || (c == '-'))
		  break;
		c = getc(fin); 
	      }
	      
	      i = 0;
	      strncpy(temp2, "\0", str_length);	      
	      while(!isspace(c)) {
		temp2[i] = (char) c;
		i++;
		c = getc(fin);
		if(i > str_length) {
		  printf("read_IC(): String overflow.  Exiting ...\n");
		  return(ERROR);
		}
	      }

	      U = strtod(temp2, NULL); 	

	      count_vel++;
	    }
	    else if(strncmp(temp2, "V", 1) == 0) {
	      while(c != EOF) {
		if(isdigit(c) || (c == '.') || (c == '-'))
		  break;
		c = getc(fin); 
	      }

	      i = 0;
	      strncpy(temp2, "\0", str_length);	      
	      while(!isspace(c)) {
		temp2[i] = (char) c;
		i++;
		c = getc(fin);
		if(i > str_length) {
		  printf("read_IC(): String overflow.  Exiting ...\n");
		  return(ERROR);
		}
	      }

	      V = strtod(temp2, NULL);	      

	      count_vel++;
	    }
	    else if(strncmp(temp2, "W", 1) == 0) {
	      while(c != EOF) {
		if(isdigit(c) || (c == '.') || (c == '-'))
		  break;
		c = getc(fin); 
	      }
	      
	      i = 0;
	      strncpy(temp2, "\0", str_length);	      
	      while(!isspace(c)) {
		temp2[i] = (char) c;
		i++;
		c = getc(fin);
		if(i > str_length) {
		  printf("read_IC(): String overflow.  Exiting ...\n");
		  return(ERROR);
		}
	      }

	      W = strtod(temp2, NULL);	      

	      count_vel++;
	    }
	    else if(strncmp(temp2, "File", 4) == 0) { 
	      if(strncmp(temp1, "Products", strlen(temp1)) != 0) {
		printf("read_IC(): Found \"File\" keyword in descripton of conditions not in explosion products gas.  Exiting ...\n");	
		return(ERROR);
	      }
	      
	      if(region_seen == TRUE) {
		printf("read_IC(): More than 1 region file listed.  Exiting ...\n");
		return(ERROR);
	      }

	      while(!isalpha(c) && (c != '.') && !isdigit(c) && (c != '_'))
		c = getc(fin);

	      i = 0;
	      strncpy(temp3, "\0", str_length);
	      while(isalpha(c) || (c == '.') || isdigit(c) || (c == '_')) {
		temp3[i] = (char) c;
		i++;
		c = getc(fin);
		if(i > str_length) {
		  printf("read_IC(): String overflow.  Exiting ...\n");
		  return(ERROR);
		}
	      }

	      region_seen = TRUE;
	      num_ICs++;
	      
	      if(num_ICs == 1) {
		IC_regions = malloc(sizeof(Body)*num_ICs);
		if(IC_regions == NULL) {
		  printf("read_IC(): Couldn't allocate array of specific regions.  Exiting ...\n");
		  return(ERROR);
		} 

		IC_bboxes = malloc(sizeof(double **)*num_ICs);
		if(IC_bboxes== NULL) {
		  printf("read_IC(): Couldn't allocate array of bounding boxes.  Exiting ...\n");
		  return(ERROR);
		}

		IC_bboxes[0] = malloc(sizeof(double *)*2);
		if(IC_bboxes[0] == NULL) {
		  printf("read_IC(): Couldn't allocate array of bounding boxes.  Exiting ...\n");
		  return(ERROR);
		}

		IC_bboxes[0][0] = malloc(sizeof(double)*3);
		if(IC_bboxes[0][0] == NULL) {
		  printf("read_IC(): Couldn't allocate array of bounding boxes.  Exiting ...\n");
		  return(ERROR);
		}

		IC_bboxes[0][1] = malloc(sizeof(double)*3);
		if(IC_bboxes[0][1] == NULL) {
		  printf("read_IC(): Couldn't allocate array of bounding boxes.  Exiting ...\n");
		  return(ERROR);
		}
	      }
	      else {
		IC_regions = realloc(IC_regions ,sizeof(Body)*num_ICs);
		if(IC_regions == NULL) {
		  printf("read_IC(): Couldn't allocate array of specific regions.  Exiting ...\n");
		  return(ERROR);
		}

		IC_bboxes = realloc(IC_bboxes, sizeof(double **)*num_ICs);
		if(IC_bboxes == NULL) {
		  printf("read_IC(): Couldn't allocate array of bounding boxes.  Exiting ...\n");
		  return(ERROR);
		}

		IC_bboxes[num_ICs-1] = malloc(sizeof(double *)*2);
		if(IC_bboxes[0] == NULL) {
		  printf("read_IC(): Couldn't allocate array of bounding boxes.  Exiting ...\n");
		  return(ERROR);
		}

		IC_bboxes[num_ICs-1][0] = malloc(sizeof(double)*3);
		if(IC_bboxes[num_ICs-1][0] == NULL) {
		  printf("read_IC(): Couldn't allocate array of bounding boxes.  Exiting ...\n");
		  return(ERROR);
		}

		IC_bboxes[num_ICs-1][1] = malloc(sizeof(double)*3);
		if(IC_bboxes[num_ICs-1][1] == NULL) {
		  printf("read_IC(): Couldn't allocate array of bounding boxes.  Exiting ...\n");
		  return(ERROR);
		}		
	      }
	      
	      /*Process the body here given by file in temp3*/

	      if(process_body(temp3, TRUE, &(IC_regions[num_ICs-1]), IC_bboxes[num_ICs-1]) == ERROR)
		return(ERROR); 
	    }
	  }
 
	  /*So we've finished reading the 5 values*/

	  if(ambient_now == TRUE) {

	    R = (gas_data.Cv_amb)*((gas_data.gam_amb) - 1.0);

	    if(rho == 0) 
	      rho = p/(R*T);	    
	    else if(p == 0) 
	      p = rho*R*T;

	    IC_ambient.Rhof = rho;
	    IC_ambient.Pres = p;
	    IC_ambient.T = p/(rho*R);
	    IC_ambient.RhoU = rho*U;
	    IC_ambient.RhoV = rho*V;
	    IC_ambient.RhoW = rho*W;
	    IC_ambient.RhoE = p/((gas_data.gam_amb)-1.0) + 0.5*(SQR(IC_ambient.RhoU) + SQR(IC_ambient.RhoV) + SQR(IC_ambient.RhoW))/rho;
	    IC_ambient.rhof_prod = 0; /*Explosion products gas not initially mixed with ambient gas*/
#if AFTERBURN
	    IC_ambient.rhof_prod_ub = 0;
#endif
	    	    	    
	    if(IC_ambient.Rhof <= 0) {
	      printf("Read_IC(): Non-physical density in ambient conditions.  Exiting ...\n");
	      return(ERROR);
	    }
	    else if(IC_ambient.Pres <= 0) {
	      printf("Read_IC(): Non-physical pressure in ambient conditions.  Exiting ...\n");
	      return(ERROR);
	    }
	    else if(IC_ambient.RhoE <= 0) {
	      printf("Read_IC(): Non-physical energy in ambient conditions.  Exiting ...\n");
	      return(ERROR);
	    }

	    dealt_ambient = TRUE;
	  }
	  else {
	    if(region_seen == FALSE) {
	      printf("Read_IC(): No input region specified.  Exiting ...\n");
	      return(ERROR);
	    }
	    
	    /*Explosion products properties*/
	    	    
	    if(gas_data.eos_prod == PERFECT_GAS) {
	      R = (gas_data.Cv_prod)*((gas_data.gam_prod) - 1.0);
	      if(rho == 0) 
		rho = p/(R*T);	    
	      else if(p == 0) 
		p = rho*R*T;
	      T = p/(rho*R); /*Compute it anyway*/
	    }
	    else if(gas_data.eos_prod == JWLB) {
	      if(rho == 0) {
		printf("Read_IC(): JWLB equation of state hard to invert to find rho = rho(p, T).  Use rho,P or rho,T.  Exiting ...\n");
		return(ERROR);
	      }
	      else { /*We need to find either temperature or pressure*/
		v = (1.0/rho)/(JWLB_coeffs.v0); /*Relative specific volume - which should be initial CJ value*/
		
		/*Precompute some quantities*/
		for(j = 0; j < 5; j++) {
		  if(JWLB_coeffs.R[j] > 0)
		    expRv[j] = exp(-(JWLB_coeffs.R[j])*v);

		  if(JWLB_coeffs.RL[j] > 0) 
		    expRLv[j] = exp(-(JWLB_coeffs.RL[j])*v);
		}
		q = JWLB_q(v, expRLv, expRL);

		R = JWLB_coeffs.omega;
		for(j = 0; j < 5; j++)
		  if(JWLB_coeffs.RL[j] > 0)
		    R += ((JWLB_coeffs.AL[j])*v + (JWLB_coeffs.BL[j]))*expRLv[j]; /*Get lambda value*/

		term1 = 0;
		for(j = 0; j < 5; j++)
		  if(JWLB_coeffs.R[j] > 0)
		    term1 += (JWLB_coeffs.A[j])*expRv[j];

		term3 = (JWLB_coeffs.C)*(1.0 - R*q/(JWLB_coeffs.omega))*pow(v, -(JWLB_coeffs.omega)-1.0);

		if(p == 0)
		  p = term1 + R*(gas_data.Cv_prod)*T/((JWLB_coeffs.v0)*v) + term3;
		else T = (p - term1 - term3)*((JWLB_coeffs.v0)*v)/(R*(gas_data.Cv_prod));
	      }
	    }
	    
	    S -> Rhof = rho;
	    S -> Pres = p;
	    S -> T = T;
	    S -> RhoU = rho*U;
	    S -> RhoV = rho*V;
	    S -> RhoW = rho*W;
	    S -> rhof_prod = rho;
#if AFTERBURN
	    S -> rhof_prod_ub = rho;
#endif

	    if(gas_data.eos_prod == PERFECT_GAS) 
	      S -> RhoE = p/((gas_data.gam_prod)-1.0) + 0.5*(SQR(S -> RhoU) + SQR(S -> RhoV) + SQR(S -> RhoW))/rho;
	    else if(gas_data.eos_prod == JWLB) {
	      v = (1.0/rho)/(JWLB_coeffs.v0);

	      R = JWLB_coeffs.omega;
	      for(j = 0; j < 5; j++)
		if(JWLB_coeffs.RL[j] > 0)
		  R += ((JWLB_coeffs.AL[j])*v + (JWLB_coeffs.BL[j]))*exp(-(JWLB_coeffs.RL[j])*v);

	      term1 = 0;
	      for(j = 0; j < 5; j++)
		if(JWLB_coeffs.R[j] > 0)
		  term1 += (JWLB_coeffs.A[j])*(1.0 - R/(v*JWLB_coeffs.R[j]))*exp(-(JWLB_coeffs.R[j])*v);

	      term3 = (JWLB_coeffs.C)*(1.0 - R/JWLB_coeffs.omega)*pow(v, -(JWLB_coeffs.omega)-1.0);

	      S -> RhoE =  (p - term1 - term3)/R + 0.5*(SQR(S -> RhoU) + SQR(S -> RhoV) + SQR(S -> RhoW))/rho;
	    }
	    	    
	    if(S -> Rhof <= 0) {
	      printf("Read_IC(): Non-physical density in explosion products gas conditions.  Exiting ...\n");
	      return(ERROR);
	    }
	    else if(S -> Pres <= 0) {
	      printf("Read_IC(): Non-physical pressure in explosion products gas conditions.  Exiting ...\n");
	      return(ERROR);
	    }
	    else if(S -> RhoE <= 0) {
	      printf("Read_IC(): Non-physical energy in explosion products gas conditions.  Exiting ...\n");
	      return(ERROR);
	    }
	    
	    if(num_ICs == 1)
	      IC_states = malloc(sizeof(State_vector *)*num_ICs);
	    else IC_states = realloc(IC_states, sizeof(State_vector *)*num_ICs);
 
	    IC_states[num_ICs-1] = S;
	  }
	}
      }
      else if(strncmp(temp1, "Heat", strlen(temp1)) == 0) { /*Read in unsteady heat addition to explosion products volume*/
	while(isspace(c))
	  c = getc(fin);

	if(isalpha(c)) {
	  i = 0;
	  strncpy(temp2, "\0", str_length);
	  while(isalpha(c)) {
	    temp2[i] = (char) c;
	    i++;
	    c = getc(fin);

	    if(c == EOF) {
	      printf("read_IC(): Inappropriate end of file.  Exiting ...\n"); 
	      return(ERROR);
	    }

	    if(i > str_length) {
	      printf("read_IC(): String overflow.  Exiting ...\n");
	      return(ERROR);
	    }   
	  }
	}

	if(strncmp(temp2, "file", strlen(temp2)) == 0) { /*Read in the file containing transient heat addition values*/
	  
	  if(gas_data.num_species > 1) {

	    while(!isalpha(c) && (c != '.') && !isdigit(c))
	      c = getc(fin);

	    /*Get the file name*/
	    i = 0;
	    strncpy(temp3, "\0", str_length);
	    while(isalpha(c) || (c == '.') || isdigit(c) || (c == '_')) {
	      temp3[i] = (char) c;
	      i++;
	      c = getc(fin);
	      if(i > str_length) {
		printf("read_IC(): String overflow.  Exiting ...\n");
		return(ERROR);
	      }
	    }

	    if(get_unsteady_IC(temp3, &ICheat) == ERROR)
	      return(ERROR);
	  }
	  else printf("read_IC(): Warning - no. species = 1, thus no heat addition to explosion prodcuts (must be different species)\n");
	}
      }      
    }
  }

  fclose(fin);

  if(dealt_ambient == FALSE) {
    printf("Read_IC(): Didn't see proper declaration of ambient IC.  Exiting ...\n");
    return(ERROR);
  }

  return(NOERROR);
}

/*------------------------------------------------------------------*/

/**\brief Read in boundary conditions for border faces of root octree cell*/

short int read_BC(char file[])
{
  FILE *fin;
  char temp1[256], temp2[256], temp3[256];
  int c = 'A';
  int i;
  int str_length = 256;
  short int count_gas, count_vel;
  short int dirs_seen = 0;
  double val;
  BC_data *BC_dat;
  BC_dat = 0;
  fin = fopen(file, "r");
  if(fin == NULL) {
    printf("read_BC(): Couldn't open boundary condition file \"%s\".  Exiting ...\n", file);
    return(ERROR);
  }

  while((dirs_seen < 6) && (c != EOF)) { /*Expect to cycle through 6 directions*/
    
    while(c != EOF) {
      c = getc(fin);
      if(isalpha(c))
	break;
    }

    i = 0;
    strncpy(temp1, "\0", str_length);
    while(isalpha(c)) {
      temp1[i] = (char) c; 
      i++;
      c = getc(fin);
      if(i > str_length) {
	printf("Read_BC(): String overflow.  Exiting ...\n"); 
	return(ERROR);
      }
      if(c == EOF) {
	printf("Read_BC(): Inappropriate end of file.  Exiting ...\n");
	return(ERROR);
      }
    }
    
    if((strncmp(temp1, "East", 4) == 0) || (strncmp(temp1, "West", 4) == 0) || 
       (strncmp(temp1, "North", 5) == 0) || (strncmp(temp1, "South", 5) == 0) || 
       (strncmp(temp1, "Upper", 5) == 0) || (strncmp(temp1, "Lower", 5) == 0)) {
      
      /*We are in 1 block*/

      while(!isalpha(c) && (c != EOF)) /*We want to get the next word*/
	c = getc(fin);

      if(c == EOF) {
	printf("Read_BC(): Inappropriate end of file.  Exiting ...\n");
	return(ERROR);
      }

      i = 0;
      strncpy(temp2, "\0", str_length);
      while(isalpha(c)) {
	temp2[i] = (char) c;
	i++;
	c = getc(fin);
	if(i > str_length) {        
	  printf("Read_BC(): String overflow.  Exiting ...\n");
	  return(ERROR);
	}
	if(c == EOF) {
	  printf("Read_BC(): Inappropriate end of file.  Exiting ...\n");
	  return(ERROR);
	}
      }

      if(strncmp(temp2, "boundary", strlen(temp2)) != 0) {
	printf("Read_BC(): Found qualifier \"%s\" (1st word in line), but next word, which should be \"boundary\", is \"%s\".  Exiting ...\n", 
	       temp1, temp2); 
	return(ERROR);
      }
      else { 
	if(strncmp(temp1, "East", strlen(temp1)) == 0)
	  BC_dat = &BC_est;
	else if(strncmp(temp1, "West", strlen(temp1)) == 0)
	  BC_dat = &BC_wst;
	else if(strncmp(temp1, "North", strlen(temp1)) == 0)
	  BC_dat = &BC_nth;
	else if(strncmp(temp1, "South", strlen(temp1)) == 0)
	  BC_dat = &BC_sth;
	else if(strncmp(temp1, "Upper", strlen(temp1)) == 0)
	  BC_dat = &BC_upr;
	else if(strncmp(temp1, "Lower", strlen(temp1)) == 0)
	  BC_dat = &BC_lwr;

	BC_dat -> type = UNKNOWN;
	BC_dat -> rhof_stat[0] = UNKNOWN;
	BC_dat -> pres_stat[0] = UNKNOWN;
	BC_dat -> T_stat[0] = UNKNOWN;
	BC_dat -> U_stat[0] = UNKNOWN;
	BC_dat -> V_stat[0] = UNKNOWN;
	BC_dat -> W_stat[0] = UNKNOWN;
	BC_dat -> num_pts = 0;
	BC_dat -> step = 0;
	BC_dat -> stat = NULL;

	/*Must first read in the BC type*/

	while(!isalpha(c) && (c != EOF))
	  c = getc(fin);

	if(c == EOF) {
	  printf("Read_BC(): Inappropriate end of file.  Exiting ...\n");
	  return(ERROR);
	}

	i = 0;
	strncpy(temp2, "\0", str_length);
	while(isalpha(c)) {
	  temp2[i] = (char) c;
	  i++;
	  c = getc(fin);
	  if(i > str_length) {
	    printf("Read_BC(): String overflow.  Exiting ...\n");
	    return(ERROR);
	  }
	  if(c == EOF) {
	    printf("Read_BC(): Inappropriate end of file.  Exiting ...\n");
	    return(ERROR);
	  }
	}

	while(!isalpha(c) && (c != EOF))
	  c = getc(fin);

	if(c == EOF) {
	  printf("Read_BC(): Inappropriate end of file.  Exiting ...\n");
	  return(ERROR);
	}

	i = 0;
	strncpy(temp3, "\0", str_length);
	while(isalpha(c)) {
	  temp3[i] = (char) c;
	  i++;
	  c = getc(fin);
	  if(i > str_length) {
	    printf("Read_BC(): String overflow.  Exiting ...\n");
	    return(ERROR);
	  }
	  if(c == EOF) {
	    printf("Read_BC(): Inappropriate end of file.  Exiting ...\n");
	    return(ERROR);
	  }
	}
	
	if((strncmp(temp2, "BC", strlen(temp2)) != 0) || (strncmp(temp3, "type", 4) != 0)) {
	  printf("Read_BC(): Next qualifiers after boundary declaration should be \"BC type\", but found \"%s %s\".  Exiting ...\n",
		 temp2,temp3);
	  return(ERROR);
	}
	else {
	  while(!isalpha(c) && (c != EOF))
	    c = getc(fin);

	  if(c == EOF) {
	    printf("Read_BC(): Inappropriate end of file.  Exiting ...\n");
	    return(ERROR);
	  }

	  i = 0;
	  strncpy(temp2, "\0", str_length);
	  while(!isspace(c) && (c != EOF)) {
	    temp2[i] = (char) c;
	    i++;
	    c = getc(fin);
	    if(i > str_length) {
	      printf("Read_BC(): String overflow.  Exiting ...\n");
	      return(ERROR);
	    }
	    if(c == EOF) {
	      printf("Read_BC(): Inappropriate end of file.  Exiting ...\n");
	      return(ERROR);
	    }
	  }

	  if((strncmp(temp2, "Non-reflecting", strlen(temp2)) == 0) || (strncmp(temp2, "Specific", strlen(temp2)) == 0) ||
	     (strncmp(temp2, "Transient", strlen(temp2)) == 0)) {
	    
	    if(strncmp(temp2, "Non-reflecting", strlen(temp2)) == 0) 
	      BC_dat -> type = NONREFLECTING; /*Will try and implement the non-reflecting boundary conditions*/ 
	    else if(strncmp(temp2, "Transient", strlen(temp2)) == 0) { /*Get table of values*/
	      BC_dat -> type = TRANSIENT;

	      while(!isalpha(c) && (c != EOF))
		c = getc(fin);
	      if(c == EOF) {
		printf("Read_BC(): Inappropriate end of file.  Exiting ...\n");
		return(ERROR);
	      }

	      i = 0;
	      strncpy(temp3, "\0", str_length);
	      while(isalpha(c)) {
		temp3[i] = (char) c;
		i++;
		c = getc(fin);
		if(i > str_length) {
		  printf("Read_BC(): String overflow.  Exiting ...\n");
		  return(ERROR);
		}
		if(c == EOF) {
		  printf("Read_BC(): Inappropriate end of file.  Exiting ...\n");
		  return(ERROR);
		}
	      } 

	      if(strncmp(temp3, "File", 4) == 0) {
		while(c != EOF) {
		  if((c == '.') || isalpha(c) || (c == '_') || isdigit(c))
		    break;
		  c = getc(fin);
		}
		  
		if(c == EOF) {
		  printf("Read_BC(): Inappropriate end of file.  Exiting ...\n");
		  return(ERROR);
		}

		/*Try and get the filename*/
		i = 0;
		strncpy(temp3, "\0", str_length);
		while(!isspace(c) && (c != EOF)) {
		  temp3[i] = (char) c; 
		  i++;
		  c = getc(fin);
		  if(i > str_length) {
		    printf("Read_BC(): String overflow.  Exiting ...\n");
		    return(ERROR);
		  }
		  if(c == EOF) {
		    printf("Read_BC(): Inappropriate end of file.  Exiting ...\n");
		    return(ERROR);
		  }
		}

		/*Now process this file*/ 
		if(get_unsteady_BC(temp3, BC_dat) == ERROR)
		  return(ERROR);		
	      }
	      else {
		printf("Read_BC(): If BC type is unsteady, qualifier \"File\" should be seen, but insteady found \"%s\".  Exiting ...\n",
		       temp3);
		return(ERROR);
	      }
	    }
	    else {
	      count_gas = 0;
	      count_vel = 0;

	      BC_dat -> type = SPECIFIC;

	      while((count_gas < 2) || (count_vel < 3)) { /*Need 2 state variables and 3 velocities*/ 
		while(!isalpha(c) && (c != EOF))
		  c = getc(fin);
		if(c == EOF) {
		  printf("Read_BC(): Inappropriate end of file.  Exiting ...\n"); 
		  return(ERROR);
		}
	  
		i = 0;
		strncpy(temp2, "\0", str_length);
		while(isalpha(c)) {
		  temp2[i] = (char) c;
		  i++;
		  c = getc(fin);
		  if(i > str_length) {
		    printf("Read_BC(): String overflow.  Exiting ...\n");
		    return(ERROR);
		  }
		  if(c == EOF) {
		    printf("Read_BC(): Inappropriate end of file.  Exiting ...\n");
		    return(ERROR);
		  }
		}
		
		if((strncmp(temp2,"Density",7) == 0) || (strncmp(temp2,"Pressure",8) == 0) || (strncmp(temp2,"Temperature",11) == 0) ||
		   (strncmp(temp2, "U", 1) == 0) || (strncmp(temp2, "V", 1) == 0) || (strncmp(temp2, "W", 1) == 0)) { 
		  /*Now try and read in the values*/
		  
		  while(c != EOF) {
		    if(isdigit(c) || (c == '.') || isalpha(c) || (c == '-'))
		      break;
		    c = getc(fin);
		  } 

		  /*Now read the value in*/
		  i = 0;
		  strncpy(temp3, "\0", str_length);
		  while((!isspace(c)) && (c != EOF)) {
		    temp3[i] = (char) c; 
		    i++;
		    c = getc(fin);
		    if(i > str_length) {
		      printf("Read_BC(): String overflow.  Exiting ...\n");
		      return(ERROR);
		    }
		    if(c == EOF) {
		      printf("Read_BC(): Inappropriate end of file.  Exiting ...\n");
		      return(ERROR);
		    }
		  }

		  if(strncmp(temp3, "Extrapolated", strlen(temp3)) == 0) { /*Want this variable to be extrapolated*/
		    if((strncmp(temp2, "Density", 7) == 0) && (count_gas < 2)) {
		      count_gas++; 
		      BC_dat -> rhof_stat[0] = EXTRAPOLATED;		      
		    }
		    else if((strncmp(temp2, "Pressure", 8) == 0) && (count_gas < 2)) {
		      count_gas++; 
		      BC_dat -> pres_stat[0] = EXTRAPOLATED;
		    }
		    else if((strncmp(temp2, "Temperature", 11) == 0) && (count_gas < 2)) {
		      count_gas++; 
		      BC_dat -> T_stat[0] = EXTRAPOLATED;
		    }
		    else if((strncmp(temp2, "U", 1) == 0) && (count_vel < 3)) {
		      count_vel++;
		      BC_dat -> U_stat[0] = EXTRAPOLATED;
		    }
		    else if((strncmp(temp2, "V", 1) == 0) && (count_vel < 3)) {
		      count_vel++;
		      BC_dat -> V_stat[0] = EXTRAPOLATED;
		    }
		    else if((strncmp(temp2, "W", 1) == 0) && (count_vel < 3)) {
		      count_vel++;
		      BC_dat -> W_stat[0] = EXTRAPOLATED;
		    }
		  }
		  else {
		    val = strtod(temp3, NULL);

		    if(((strncmp(temp2, "Density", 7) == 0) || (strncmp(temp2, "Pressure", 8) == 0) || 
			(strncmp(temp2, "Temperature", 11) == 0)) && (val <= 0)) {
		      printf("Read_BC(): Found non-physical %s value, which is \"%s\".  Exiting ...\n", temp2, temp3);
		      return(ERROR);
		    }
		    
		    if((strncmp(temp2, "Density", 7) == 0) && (count_gas < 2)) {
		      count_gas++; 
		      BC_dat -> rhof_stat[0] = FIXED; 
		      BC_dat -> rhof_stat[1] = val;
		    }
		    else if((strncmp(temp2, "Pressure", 8) == 0) && (count_gas < 2)) {
		      count_gas++; 
		      BC_dat -> pres_stat[0] = FIXED; 
		      BC_dat -> pres_stat[1] = val;
		    }
		    else if((strncmp(temp2, "Temperature", 11) == 0) && (count_gas < 2)) {
		      count_gas++; 
		      BC_dat -> T_stat[0] = FIXED; 
		      BC_dat -> T_stat[1] = val;
		    }
		    else if((strncmp(temp2, "U", 1) == 0) && (count_vel < 3)) {
		      count_vel++;
		      BC_dat -> U_stat[0] = FIXED;
		      BC_dat -> U_stat[1] = val;
		    }
		    else if((strncmp(temp2, "V", 1) == 0) && (count_vel < 3)) {
		      count_vel++;
		      BC_dat -> V_stat[0] = FIXED;
		      BC_dat -> V_stat[1] = val;
		    }
		    else if((strncmp(temp2, "W", 1) == 0) && (count_vel < 3)) {
		      count_vel++;
		      BC_dat -> W_stat[0] = FIXED;
		      BC_dat -> W_stat[1] = val;
		    }
		  }
		}
		else {
		  printf("Read_BC(): Unknown flow parameter \"%s\".  Exiting ...\n", temp2);
		  return(ERROR);
		}
	      }	 

	      if((gas_data.eos_prod == JWLB) && 
#if FLOAT_EQ_STABLE
		 EQ(BC_dat -> rhof_stat[0], UNKNOWN)
#else
		 (BC_dat -> rhof_stat[0] == UNKNOWN)
#endif
		 ) {
		printf("Read_BC(): JWLB equation of state difficult to find rho = rho(p, T).  Exiting ...\n");
		return(ERROR);
	      }

	      if((count_gas < 2) || (count_vel < 3)) {
		printf("Read_BC(): 2 gas parameters and 3 velocity parameters not seen.  Exiting ...\n"); 
		return(ERROR);
	      }
	    }

	    dirs_seen++;
	  }
	  else {
	    printf("Read_BC(): Unknown BC type.  Exiting ...\n");
	    return(ERROR);
	  }
	}
      }
    }
  }
  
  if(dirs_seen != 6) {
    printf("Read_BC(): Not all 6 borders appropriately specified.  Exiting ...\n"); 
    return(ERROR);
  }

  return(NOERROR); 
}


/*------------------------------------------------------------------*/

/**\brief Get table of transient values for unsteady heat addition to IC volume*/

short int get_unsteady_IC(char file_name[], IC_heat *ICh)
{
  int i;
  char buffer[512];
  FILE *fin = fopen(file_name, "r");
  if(fin == NULL) {
    printf("get_unsteady_IC(): Couldn't open file %s.  Exiting ...\n", file_name);
    return(ERROR);
  }

  /*First line should be no. of points*/
  
  fgets(buffer, sizeof(buffer), fin);
  if(sscanf(buffer, "No. points: %d", &(ICh -> num_pts)) != 1) {
    printf("get_unsteady_IC(): Didn't read an appropriate no. time points in file %s.  Exiting ...\n", file_name);
    return(ERROR);
  }
  
  ICh -> heat = malloc(sizeof(double *)*(ICh -> num_pts));
  if(ICh -> heat == NULL) {
    printf("get_unsteady_IC(): Couldn't allocate time matrix.  Exiting ...\n");
    return(ERROR);
  }

  for(i = 0; i < (ICh -> num_pts); i++) {
    ICh -> heat[i] = malloc(sizeof(double)*4); /*3 columns - time, fuel consumption rate (kg/s), energy release rate (J/(kg.s)),
						 ignition temperature (K)*/
    if(ICh -> heat[i] == NULL) {
      printf("get_unsteady_IC(): Couldn't allocate time matrix.  Exiting ...\n");
      return(ERROR);
    }
  }

  /*Now proceed to read in transient values*/

  for(i = 0; i < (ICh -> num_pts); i++) {
    fgets(buffer, sizeof(buffer), fin);
    
    if(sscanf(buffer, "%lf %lf %lf %lf", &(ICh -> heat[i][0]), &(ICh -> heat[i][1]), &(ICh -> heat[i][2]), &(ICh -> heat[i][3])) != 2) {
      printf("get_unsteady_IC(): Couldn't read in 4 transient values in line %d.  Exiting ...\n", i+2);
      return(ERROR);
    }

    if(i == 0) {
      if(ICh -> heat[0][0] != 0) {
	printf("get_unsteady_IC(): Time history must begin at time 0 (but is %e).  Exiting ...\n",ICh -> heat[0][0]);
	return(ERROR);
      }
    }
    else if(ICh -> heat[i][0] <= ICh -> heat[i-1][0]) {
      printf("get_unsteady_IC(): Found stationary/decreasing time in line %d.  Exiting ...\n", i+2);
      return(ERROR);
    }
  }

  return(NOERROR);
}

/*------------------------------------------------------------------*/

/**\brief Get table of transient values for unsteady BCs*/

short int get_unsteady_BC(char file_name[], BC_data *BC_dat)
{
  int i, j;
  int str_length = 72;
  char buffer[512];
  char temp1[30], temp2[30];
  double t, q1, q2;
  FILE *fin = fopen(file_name, "r");
  if(fin == NULL) {
    printf("get_unsteady_BC(): Couldn't open file %s.  Exiting ...\n", file_name);
    return(ERROR);
  }

  /*First line should be no. of points*/
  
  fgets(buffer, sizeof(buffer), fin);
  if(sscanf(buffer, "%d", &(BC_dat -> num_pts)) != 1) {
    printf("get_unsteady_BC(): Didn't read an appropriate no. time points in file %s.  Exiting ...\n", file_name);
    return(ERROR);
  }
  
  BC_dat -> stat = malloc(sizeof(char **)*(BC_dat -> num_pts));
  if(BC_dat -> stat == NULL) {
    printf("get_unsteady_BC(): Couldn't allocate time matrix.  Exiting ...\n");
    return(ERROR);
  }
  
  for(i = 0; i < (BC_dat -> num_pts); i++) {
    BC_dat -> stat[i] = malloc(sizeof(char *)*6);
    if(BC_dat -> stat[i] == NULL) {
      printf("get_unsteady_BC(): Couldn't allocate time matrix.  Exiting ...\n");
      return(ERROR);
    }

    for(j = 0; j < 6; j++) {
      BC_dat -> stat[i][j] = malloc(sizeof(char)*str_length); /*Values stored in string format as some might be extrapolated*/
      if(BC_dat -> stat[i][j] == NULL) {
	printf("get_unsteady_BC(): Couldn't allocate time matrix.  Exiting ...\n");
	return(ERROR);
      }

      strncpy(BC_dat -> stat[i][j], "\0", str_length);
    }
  }
  
  /*2nd line says which gas variables are specified*/

  fgets(buffer, sizeof(buffer), fin);
  if(sscanf(buffer, "%s %s", temp1, temp2) != 2) {
    printf("get_unsteady_BC(): Didn't read which gas variables are specified in 1st two entries.  Exiting ...\n");
    return(ERROR);
  }

  BC_dat -> temp_spec = 0; 
  BC_dat -> pres_spec = 0;
  BC_dat -> rhof_spec = 0;
  
  if(strcmp(temp1, "Temperature") == 0)
    BC_dat -> temp_spec = 1;
  else if(strcmp(temp1, "Pressure") == 0)
    BC_dat -> pres_spec = 1;
  else if(strcmp(temp1, "Density") == 0)
    BC_dat -> rhof_spec = 1;
  else {
    printf("get_unsteady_BC(): Unknown gas variable \"%s\"; expected \"Temperature\", \"Pressure\" or \"Density\".  Exiting ...\n", temp1);
    return(ERROR);
  }
  
  if(strcmp(temp2, "Temperature") == 0) {
    if(BC_dat -> temp_spec == 1) {
      printf("get_unsteady_BC(): Temperature already specified.  Exiting ...\n");
      return(ERROR);
    }
    BC_dat -> temp_spec = 2;
  }
  else if(strcmp(temp2, "Pressure") == 0) {
    if(BC_dat -> pres_spec == 1) {
      printf("get_unsteady_BC(): Pressure already specified.  Exiting ...\n");
      return(ERROR);
    }
    BC_dat -> pres_spec = 2;
  }
  else if(strcmp(temp2, "Density") == 0) {
    if(BC_dat -> rhof_spec == 1) {
      printf("get_unsteady_BC(): Density already specified.  Exiting ...\n");
      return(ERROR);
    }
    BC_dat -> rhof_spec = 2;
  }
  else {
    printf("get_unsteady_BC(): Unknown gas variable \"%s\"; expected \"Temperature\", \"Pressure\" or \"Density\".  Exiting ...\n", temp2);
    return(ERROR);
  }  

  if((gas_data.eos_prod == JWLB) && (BC_dat -> rhof_spec == 0)) {
    printf("get_unsteady_BC(): JWLB equation of state difficult to find rho = rho(p, T).  Please fix rho.  Exiting ...\n");
    return(ERROR);
  }

  /*
    Now we should read in the transient values - for each row (1 line in the file), 
    col 0 = time
    col 1-2 = gas variables in order they were written in previous line
    col 3 = u velocity
    col 4 = v velocity
    col 5 = w velocity 
  */
  
  for(i = 0; i < (BC_dat -> num_pts); i++) {
    fgets(buffer, sizeof(buffer), fin);
    
    if(sscanf(buffer, "%s %s %s %s %s %s", (BC_dat->stat[i][0]), (BC_dat->stat[i][1]), (BC_dat->stat[i][2]),
	      (BC_dat->stat[i][3]), (BC_dat->stat[i][4]), (BC_dat->stat[i][5])) != 6) {
      printf("get_unsteady_BC(): Couldn't read in 6 transient values in line %d.  Exiting ...\n", i+2);
      return(ERROR);
    }

    t = strtod(BC_dat -> stat[i][0], NULL);
    q1 = strtod(BC_dat -> stat[i][1], NULL);
    q2 = strtod(BC_dat -> stat[i][2], NULL);

    if(i == 0) {
      if(t != 0) {
	printf("get_unsteady_BC(): Time history must begin at time 0.  Exiting ...\n");
	return(ERROR);
      }
    }
    else if(t <= strtod(BC_dat -> stat[i-1][0], NULL)) {
      printf("get_unsteady_BC(): Found stationary/decreasing time in line %d.  Exiting ...\n", i+2);
      return(ERROR);
    }

    if((q1 <= 0) || (q2 <= 0)) {
      printf("get_unsteady_BC(): Non-physical value in gas parameters seen in line %d.  Exiting ...\n", i+2);
      return(ERROR);
    }
  }

  fclose(fin);
  return(NOERROR);
}

/*------------------------------------------------------------------*/

/**\brief Function to write in flow solution of all (or just intersected) cells in serial/parallel*/

short int write_soln(int run_timesteps_so_far, int num_cells_flow, double current_time)
{ 
  int num_cells;
  char tmp[128];
  int thread_num;
  short int error_found = FALSE;
  List_leaf Head, Tail;
  Head = Leaves;
  Tail = Leaves_tail;

#if SERIAL_WRITE
  List_leaf L;
#endif
  thread_num = 0;
#if 0
  L = Leaves;
  Cart_cell C;
  while(L != NULL) {
    C = L -> cell_loc;
    printf("C @ %g %g, cn %s, p %e, r %e, ru %e, rv %e, re %e\n",C->centroid[0],C->centroid[1],C->cell_num,C->Flow_data.State_vec.Pres,
	   C->Flow_data.State_vec.Rhof,C->Flow_data.State_vec.RhoU,C->Flow_data.State_vec.RhoV,C->Flow_data.State_vec.RhoE);
    L=L->next;
  }
#endif

#if SERIAL_WRITE
  nt=number_of_threads;
  number_of_threads=1;
#endif

#if PRINT_SCHL
  get_max_schl(); /*Get maximum density gradient for later schlieren normalization*/
#endif
  
  sprintf(tmp, "%s.%d", base_soln_filename, run_timesteps_so_far);
	  
  if(number_of_threads > 1)
    if(link_solution_files(tmp, current_time) == ERROR)
      return(ERROR);

   if(only_get_boundary_cells_flow == TRUE) /*Only want flow soln next to bodies*/
    {
      if(number_of_threads > 1)
	{
	  #pragma omp parallel private(Head, Tail, num_cells, thread_num)
	  {
	    #ifdef _OPENMP
	    thread_num = omp_get_thread_num();
	    #endif

	    Head = Leaf_heads[thread_num];
	    Tail = Leaf_tails[thread_num];

	    num_cells = get_num_boundary_cells(Head, Tail);

#if TRY_PLOT_SHOCK
	    compute_adapt_param(Head, Tail);
	    #pragma omp barrier
#endif
	    
	    if(error_found == FALSE)
	      if(write_oct_flow_soln(Head, Tail, TRUE, num_cells, current_time, tmp) == ERROR)
		{
		  printf("write_soln(): Error writing flow solution\n");
		  error_found = FALSE;
                  #pragma omp flush(error_found)
		}
	  }
	}
      else {
	num_cells = get_num_boundary_cells(Head, Tail);

#if TRY_PLOT_SHOCK
	    compute_adapt_param(Head, Tail);
	    #pragma omp barrier
#endif

	if(write_oct_flow_soln(Head, Tail, TRUE, num_cells, current_time, tmp) == ERROR)
	  {
	    printf("write_soln(): Error writing flow solution\n");
	    error_found = FALSE;            
	  }
      }
    }
  else 
    {   
      /*I do things this way because the portland compiler behaves funny if I compile using OpenMP 
	but only use 1 CPU*/

      if(number_of_threads > 1)
	{
	  #pragma omp parallel private(Head, Tail, num_cells, thread_num)
	  {
	    #ifdef _OPENMP
	    thread_num = omp_get_thread_num();
	    #endif

	    Head = Leaf_heads[thread_num];
	    Tail = Leaf_tails[thread_num];

	    num_cells = num_cells_per_thread[thread_num];

#if TRY_PLOT_SHOCK
	    compute_adapt_param(Head, Tail);
	    #pragma omp barrier
#endif
	    if(error_found == FALSE)
	      if(write_oct_flow_soln(Head, Tail, FALSE, num_cells, current_time, tmp) == ERROR)
		{
		  printf("write_soln(): Error writing flow solution\n");
		  error_found = FALSE;
                  #pragma omp flush(error_found)
		}     
	  }
	}
      else {
	num_cells = num_cells_flow;

#if SERIAL_WRITE
	    L=Head;
	    while(L != NULL) {
	      L->thread_num = 0;
	      L=L->next;
	    }
#endif

#if TRY_PLOT_SHOCK
	    compute_adapt_param(Head, Tail);
	    #pragma omp barrier
#endif

	if(write_oct_flow_soln(Head, Tail, FALSE, num_cells, current_time, tmp) == ERROR)
	  {
	    printf("write_soln(): Error writing flow solution\n");
	    error_found = FALSE;            
	  }
      }
    }
   


#if SERIAL_WRITE
  number_of_threads=nt;
#endif

#if SOLNSPLIT
  if(output_grid == TRUE) {
    if((run_timesteps_so_far > 0) && ((output_grid_last == FALSE) || 
       ((output_grid_last == TRUE) && (((num_timesteps_run > 0) && (run_timesteps_so_far >= num_timesteps_run)) || 
				       ((finish_time > 0) && (current_time >= finish_time)))))) {
      write_mesh_and_soln(Root, run_timesteps_so_far, current_time);
    }
  }
#endif

   if(error_found == TRUE)
     return(ERROR);
   else return(NOERROR);
}

/*------------------------------------------------------------------*/

/**\brief Link solution files (in ascii VTK XML format) when they are outputted in parallel*/

short int link_solution_files(char base_name[], double current_time)
{
  FILE *fout;
  int i;
  char oname[64];

  strncpy(oname, base_name, 64);
  strcat(oname, ".pvtu");

  fout = fopen(oname, "w");
  if(fout == NULL) 
    {
      printf("link_solution_files(): could not open file %s, exiting ...\n", oname);
      return(ERROR);
    }

  fprintf(fout, "<?xml version=\"1.0\"?>\n");
  fprintf(fout, "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
  fprintf(fout, "<note> Flow soln, time %e s </note>\n", current_time);
  fprintf(fout, "<PUnstructuredGrid GhostLevel=\"0\">\n");
  fprintf(fout, "<PPoints>\n");
  fprintf(fout, "<PDataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\"/>\n");
  fprintf(fout, "</PPoints>\n");
  fprintf(fout, "<PPointData Scalars=\"Scalars (SI units)\" Vectors=\"Vectors\">\n");
  fprintf(fout, "<PDataArray type=\"Float64\" Name=\"Pressure-Pa\" format=\"ascii\"/>\n");
  fprintf(fout, "<PDataArray type=\"Float64\" Name=\"Density-kg/m^3\" format=\"ascii\"/>\n");

#if PRINT_SCHL
  fprintf(fout, "<PDataArray type=\"Float64\" Name=\"Schlieren\" format=\"ascii\"/>\n");
#endif

#if PRINT_TEMP
  fprintf(fout, "<DataArray type=\"Float64\" Name=\"Temperature-K\" format=\"ascii\"/>\n");
#endif

  if(gas_data.num_species > 1)
    fprintf(fout, "<DataArray type=\"Float64\" Name=\"Products m_frac\" format=\"ascii\"/>\n");

#if PRINT_VEL
  fprintf(fout, "<PDataArray type=\"Float64\" Name=\"Velocity.x-m/s\" format=\"ascii\"/>\n");
  fprintf(fout, "<PDataArray type=\"Float64\" Name=\"Velocity.y-m/s\" format=\"ascii\"/>\n");
  fprintf(fout, "<PDataArray type=\"Float64\" Name=\"Velocity.z-m/s\" format=\"ascii\"/>\n");
  fprintf(fout, "<DataArray type=\"Float64\" Name=\"Velocity-m/s\" NumberOfComponents=\"3\" format=\"ascii\"/>\n");
  fprintf(fout, "<PDataArray type=\"Float64\" Name=\"Speed-m/s\" format=\"ascii\"/>\n");
  fprintf(fout, "<PDataArray type=\"Float64\" Name=\"M\" format=\"ascii\"/>\n");
#endif
  fprintf(fout, "</PPointData>\n");

  fprintf(fout, "<PDataArray type=\"Float64\" Name=\"Shocks\" format=\"ascii\"/>\n");

  for(i = 1; i <= number_of_threads; i++)
    fprintf(fout, "<Piece Source=\"%s.part%dof%d.vtu\"/>\n", base_name, i, number_of_threads); /*Collate individual pieces together*/

  fprintf(fout, "</PUnstructuredGrid>\n");
  fprintf(fout, "</VTKFile>\n");

  fclose(fout);
  return(NOERROR);
}

/*------------------------------------------------------------------*/

/**\brief Write flow solution to VTK file using cell-centred data using vertex (VTK voxel) values of each cell
   (reconstructed to cell corners for each cell)*/

short int write_oct_flow_soln(List_leaf L, List_leaf Tail, short int only_get_boundary_cells, int num_cells, double current_time, 
			      char base_name[])
{
  FILE *fout;
  List_leaf node;
  List_vtx_glob Vl;
  Cart_cell C;
  int j, thread_num;
  short int i, num_vtx_cells;
  double avg, q, prod_m_frac;
  char oname[64];
  int num_points = 0;
  char tmp[48];

#if NOTINCLUDE
  if(stepnow == 200) {
    printf("begin write_oct_flow_soln\n");
  }
#endif
  
  thread_num = L -> thread_num; /*Only deal with cells with this thread number*/

  strncpy(oname, base_name, 64);
  
  if(number_of_threads > 1)
    sprintf(tmp, ".part%dof%d.vtu",thread_num+1, number_of_threads);
  else sprintf(tmp, ".vtu");

  strcat(oname, tmp);

  fout = fopen(oname, "w");
  if(fout == NULL) 
    {
      printf("write_oct_flow_soln(): could not open file %s, exiting ...\n", oname);
      return(ERROR);
    }

  fprintf(fout, "<?xml version=\"1.0\"?>\n");
  fprintf(fout, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
  if(number_of_threads > 1)
    fprintf(fout, "<note> Partial flow soln, time %e s </note>\n", current_time);
  else fprintf(fout, "<note> Flow soln, time %e s </note>\n", current_time);
  fprintf(fout,"<UnstructuredGrid>\n");

  /*Count and number all verticies first*/
   
  Vl = Vtxlist;

  if(only_get_boundary_cells == FALSE)
    {
      while(Vl != NULL)
	{
	  /*Due to refinement a vertex that was previously plotted may be shared by all solid cells - disregard verticies
	    shared by only solid cells*/
	 
	  for(i = 0; i < MAX_NUM_CHILDREN; i++)
	    { 
	      if((Vl -> vtx_loc -> Leaf_cells[i] != NULL) && (Vl -> vtx_loc -> Leaf_cells[i] -> cell_type != SOLID) &&
#if DETON
		 (Vl -> vtx_loc -> Leaf_cells[i] -> un_det != TRUE) && 
#endif
		 (Vl -> vtx_loc -> Leaf_cells[i] -> Leaf_list_loc -> thread_num == thread_num)) 
		{ /*Not shared solely by solid cells and at least one cell has the same thread number*/

		  Vl -> vtx_loc -> vtx_nums[thread_num] = num_points;
		  num_points++;		
		  break;
		}
	    } 
	  
	  if(i == MAX_NUM_CHILDREN)
	    Vl -> vtx_loc -> vtx_nums[thread_num] = -1; /*Special flag for verticies which need not be outputted*/
	  	  	  
	  Vl = Vl -> next;
	} 
    }
  else /*See if vertex shared by any intersected or (non-pure) fluid cell*/
    {
      while(Vl != NULL)
	{
	  /*Due to refinement a vertex that was previously plotted may no longer be shared by any intersected cell*/

	  for(i = 0; i < MAX_NUM_CHILDREN; i++)
	    {
	      if((Vl -> vtx_loc -> Leaf_cells[i] != NULL) &&
#if DETON
		 (Vl -> vtx_loc -> Leaf_cells[i] -> un_det != TRUE) &&
#endif
		 ((Vl -> vtx_loc -> Leaf_cells[i] -> cell_type == INTERSECTED) || 
		  ((Vl -> vtx_loc -> Leaf_cells[i] -> Leaf_list_loc != NULL) && 
		   (Vl -> vtx_loc -> Leaf_cells[i] -> Leaf_list_loc -> pure_fluid == FALSE))) && 
		 (Vl -> vtx_loc -> Leaf_cells[i] -> Leaf_list_loc -> thread_num == thread_num))
		{
		  Vl -> vtx_loc -> vtx_nums[thread_num] = num_points;
		  num_points++;
		  break;
		}
	    }

	  if(i == MAX_NUM_CHILDREN)
	    Vl -> vtx_loc -> vtx_nums[thread_num] = -1;	    
	  	  	  
	  Vl = Vl -> next;
	}
    }
  
  fprintf(fout, "<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", num_points, num_cells);  


  /*Write the points*/
  fprintf(fout, "<Points>\n");
  fprintf(fout,"<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n");

  Vl = Vtxlist;
  while(Vl != NULL)
    {
      if(Vl -> vtx_loc -> vtx_nums[thread_num] != -1)
	fprintf(fout, "%20.12e %20.12e %20.12e\n", Vl -> vtx_loc -> loc[0], Vl -> vtx_loc -> loc[1], Vl -> vtx_loc -> loc[2]);
	      
      Vl = Vl -> next;
    }

  fprintf(fout, "</DataArray>\n");
  fprintf(fout, "</Points>\n");


  /*Write the cells*/
  fprintf(fout, "<Cells>\n");
  fprintf(fout, "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");
        
  node = L;
  if(only_get_boundary_cells == TRUE)
    {
      while(node != NULL)
	{
	  C = node -> cell_loc;

#if DETON
	  if(C -> un_det != TRUE)
#endif
	    {
	      if((C -> cell_type == INTERSECTED) || (node -> pure_fluid == FALSE))
		{
		  fprintf(fout, "%d ", ABS(C -> verticies[0] -> vtx_nums[thread_num]));
		  fprintf(fout, "%d ", ABS(C -> verticies[2] -> vtx_nums[thread_num]));
		  fprintf(fout, "%d ", ABS(C -> verticies[1] -> vtx_nums[thread_num]));
		  fprintf(fout, "%d ", ABS(C -> verticies[3] -> vtx_nums[thread_num]));
		  fprintf(fout, "%d ", ABS(C -> verticies[4] -> vtx_nums[thread_num]));
		  fprintf(fout, "%d ", ABS(C -> verticies[6] -> vtx_nums[thread_num]));
		  fprintf(fout, "%d ", ABS(C -> verticies[5] -> vtx_nums[thread_num]));
		  fprintf(fout, "%d ", ABS(C -> verticies[7] -> vtx_nums[thread_num]));
	  
		  fprintf(fout, "\n");
		}
	    }
      
	  if(node == Tail)
	    break;
	  else node = node -> next;
	}
    }
  else
    {
      while(node != NULL)
	{
	  C = node -> cell_loc;

#if DETON
	  if(C -> un_det != TRUE)
#endif
	    {
	      fprintf(fout, "%d ", ABS(C -> verticies[0] -> vtx_nums[thread_num]));
	      fprintf(fout, "%d ", ABS(C -> verticies[2] -> vtx_nums[thread_num]));
	      fprintf(fout, "%d ", ABS(C -> verticies[1] -> vtx_nums[thread_num]));
	      fprintf(fout, "%d ", ABS(C -> verticies[3] -> vtx_nums[thread_num]));
	      fprintf(fout, "%d ", ABS(C -> verticies[4] -> vtx_nums[thread_num]));
	      fprintf(fout, "%d ", ABS(C -> verticies[6] -> vtx_nums[thread_num]));
	      fprintf(fout, "%d ", ABS(C -> verticies[5] -> vtx_nums[thread_num]));
	      fprintf(fout, "%d ", ABS(C -> verticies[7] -> vtx_nums[thread_num]));

	      fprintf(fout, "\n");
	    }
	  	      
	  if(node == Tail)
	    break;
	  else node = node -> next;
	}
    }
  fprintf(fout, "</DataArray>\n");

  fprintf(fout, "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");
  for(j = 1; j <= num_cells; j++)
    fprintf(fout, "%d\n", j*8); /*Offsets start at a positive number, and always 8 verticies per cell*/
  fprintf(fout, "</DataArray>\n");

  fprintf(fout, "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
  for(j = 0; j < num_cells; j++)
    fprintf(fout, "%d\n", 11); /* 11 == VTK_VOXEL */
  fprintf(fout, "</DataArray>\n");

  fprintf(fout, "</Cells>\n");

  /*Now write some scalar/vector point data*/

  fprintf(fout, "<PointData Scalars=\"Scalars (SI units)\" Vectors=\"Vectors\">\n");
  
  fprintf(fout, "<DataArray type=\"Float64\" Name=\"Pressure-Pa\" format=\"ascii\">\n");
  Vl = Vtxlist;
  while(Vl != NULL)
    {
      if(Vl -> vtx_loc -> vtx_nums[thread_num] != -1)
	{	  
	  avg = 0; num_vtx_cells = 0;
	  
	  /*Get the average nodal value for all computational cells*/
	  
	  for(i = 0; i < MAX_NUM_CHILDREN; i++)
	    {
	      if((Vl -> vtx_loc -> Leaf_cells[i] != NULL) && (Vl -> vtx_loc -> Leaf_cells[i] -> cell_type != SOLID)
#if DETON
		 && (Vl -> vtx_loc -> Leaf_cells[i] -> un_det != TRUE)
#endif
		 )
		{
		  reconstruct_flow(Vl -> vtx_loc -> Leaf_cells[i], IRRELEVANT, IRRELEVANT, 7-i, 1, 'P', &q, NULL, 'n');
		  avg += q;

		  num_vtx_cells++;
		}
	    } 
	  
	  fprintf(fout, "%g\n", avg/num_vtx_cells);	  
	}	     
      
      Vl = Vl -> next;
    }
  fprintf(fout, "</DataArray>\n");
  
  fprintf(fout, "<DataArray type=\"Float64\" Name=\"Density-kg/m^3\" format=\"ascii\">\n");
  Vl = Vtxlist;    
  while(Vl != NULL)
    {
      if(Vl -> vtx_loc -> vtx_nums[thread_num] != -1)
	{	  
	  avg = 0; num_vtx_cells = 0;
	  
	  for(i = 0; i < MAX_NUM_CHILDREN; i++)
	    {
	      if((Vl -> vtx_loc -> Leaf_cells[i] != NULL) && (Vl -> vtx_loc -> Leaf_cells[i] -> cell_type != SOLID)
#if DETON
		 && (Vl -> vtx_loc -> Leaf_cells[i] -> un_det != TRUE)
#endif
		 )
		{
		  reconstruct_flow(Vl -> vtx_loc -> Leaf_cells[i], IRRELEVANT, IRRELEVANT, 7-i, 1, 'R', &q, NULL, 'n');
		  avg += q;
		  
		  num_vtx_cells++;
		}
	    } 
	  
	  fprintf(fout, "%g\n", avg/num_vtx_cells);
	}	     
      
      Vl = Vl -> next;
    }    
  fprintf(fout, "</DataArray>\n");

#if PRINT_SCHL
  fprintf(fout, "<DataArray type=\"Float64\" Name=\"Schlieren\" format=\"ascii\">\n");
  Vl = Vtxlist;    
  while(Vl != NULL)
    {
      if(Vl -> vtx_loc -> vtx_nums[thread_num] != -1)
	{	  
	  avg = 0; num_vtx_cells = 0;
	  
	  for(i = 0; i < MAX_NUM_CHILDREN; i++)
	    {
	      if((Vl -> vtx_loc -> Leaf_cells[i] != NULL) && (Vl -> vtx_loc -> Leaf_cells[i] -> cell_type != SOLID)
#if DETON
		 && (Vl -> vtx_loc -> Leaf_cells[i] -> un_det != TRUE)
#endif
		 )
		{
		  /*Use cell-centered schlieren for now*/
		  if(max_schl > 0)
		    q = (SQRT(SQR(Vl->vtx_loc->Leaf_cells[i]->Flow_data.Grads.grad_rhof[0]) + 
			      SQR(Vl->vtx_loc->Leaf_cells[i]->Flow_data.Grads.grad_rhof[1]) + 
			      SQR(Vl->vtx_loc->Leaf_cells[i]->Flow_data.Grads.grad_rhof[2])))/max_schl;
		  else q = 0;
		  
		  avg += q;
		  
		  num_vtx_cells++;
		}
	    } 
	  
	  fprintf(fout, "%g\n", avg/num_vtx_cells);
	}	     
      
      Vl = Vl -> next;
    }    
  fprintf(fout, "</DataArray>\n");
#endif

#if PRINT_TEMP 
  fprintf(fout, "<DataArray type=\"Float64\" Name=\"Temperature-K\" format=\"ascii\">\n");
  Vl = Vtxlist;
  while(Vl != NULL)
    {
      if(Vl -> vtx_loc -> vtx_nums[thread_num] != -1)
	{
	  avg = 0; num_vtx_cells = 0;
	 
	  for(i = 0; i < MAX_NUM_CHILDREN; i++)
	    {
	      if((Vl -> vtx_loc -> Leaf_cells[i] != NULL) && (Vl -> vtx_loc -> Leaf_cells[i] -> cell_type != SOLID)
#if DETON
		 && (Vl -> vtx_loc -> Leaf_cells[i] -> un_det != TRUE)
#endif
		 )
		{
		  q = Vl -> vtx_loc -> Leaf_cells[i] -> Flow_data.State_vec.T;
		  avg += q;
		  num_vtx_cells++;
		}
	    }
 	  
	  fprintf(fout, "%g\n", avg/num_vtx_cells);
	}	     
      
      Vl = Vl -> next;
    }
  fprintf(fout, "</DataArray>\n");
#endif

  if(gas_data.num_species > 1)
    {
      fprintf(fout, "<DataArray type=\"Float64\" Name=\"Products m_frac\" format=\"ascii\">\n");
      Vl = Vtxlist;
      while(Vl != NULL)
	{
	  if(Vl -> vtx_loc -> vtx_nums[thread_num] != -1)
	    {	  
	      avg = 0; num_vtx_cells = 0;
	  
	      for(i = 0; i < MAX_NUM_CHILDREN; i++)
		{
		  if((Vl -> vtx_loc -> Leaf_cells[i] != NULL) && (Vl -> vtx_loc -> Leaf_cells[i] -> cell_type != SOLID)
#if DETON
		     && (Vl -> vtx_loc -> Leaf_cells[i] -> un_det != TRUE)
#endif
		     )
		    {
		      /*reconstruct_flow(Vl -> vtx_loc -> Leaf_cells[i], IRRELEVANT, IRRELEVANT, 7-i, 1, 'R', &q, NULL);*/
		      prod_m_frac = (Vl -> vtx_loc -> Leaf_cells[i] -> Flow_data.State_vec.rhof_prod)/
			(Vl -> vtx_loc -> Leaf_cells[i] -> Flow_data.State_vec.Rhof);

		      avg += prod_m_frac; /*Just use cell-centred mass fraction to avoid more reconstruction*/

		      num_vtx_cells++;
		    }
		} 
	  
	      fprintf(fout, "%g\n", avg/num_vtx_cells);
	    }	     
      
	  Vl = Vl -> next;
	}
      fprintf(fout, "</DataArray>\n");
    }

#if PRINT_VEL
  fprintf(fout, "<DataArray type=\"Float64\" Name=\"Velocity.x-m/s\" format=\"ascii\">\n");
  Vl = Vtxlist;
  while(Vl != NULL)
    {
      if(Vl -> vtx_loc -> vtx_nums[thread_num] != -1)
	{
	  avg = 0; num_vtx_cells = 0;
	  	
	  for(i = 0; i < MAX_NUM_CHILDREN; i++)
	    {
	      if((Vl -> vtx_loc -> Leaf_cells[i] != NULL) && (Vl -> vtx_loc -> Leaf_cells[i] -> cell_type != SOLID)
#if DETON
		 && (Vl -> vtx_loc -> Leaf_cells[i] -> un_det != TRUE)
#endif
		 )
		{
#if RECON_VEL
		  reconstruct_flow(Vl -> vtx_loc -> Leaf_cells[i], IRRELEVANT, IRRELEVANT, 7-i, 1, 'U', &q, NULL, 'n');
#else
	q = (Vl->vtx_loc->Leaf_cells[i]->Flow_data.State_vec.RhoU)/(Vl->vtx_loc->Leaf_cells[i]->Flow_data.State_vec.Rhof);
#endif
		  avg += q;

		  num_vtx_cells++;
		}
	    }
  
	  fprintf(fout, "%g\n", avg/num_vtx_cells);
	}	     
      
      Vl = Vl -> next;
    }
  fprintf(fout, "</DataArray>\n");

  fprintf(fout, "<DataArray type=\"Float64\" Name=\"Velocity.y-m/s\" format=\"ascii\">\n");
  Vl = Vtxlist;
  while(Vl != NULL)
    {
      if(Vl -> vtx_loc -> vtx_nums[thread_num] != -1)
	{
	  avg = 0; num_vtx_cells = 0;
	  	
	  for(i = 0; i < MAX_NUM_CHILDREN; i++)
	    {
	      if((Vl -> vtx_loc -> Leaf_cells[i] != NULL) && (Vl -> vtx_loc -> Leaf_cells[i] -> cell_type != SOLID)
#if DETON
		 && (Vl -> vtx_loc -> Leaf_cells[i] -> un_det != TRUE)
#endif
		 )
		{
#if RECON_VEL
		  reconstruct_flow(Vl -> vtx_loc -> Leaf_cells[i], IRRELEVANT, IRRELEVANT, 7-i, 1, 'V', &q, NULL, 'n');
#else
q = (Vl->vtx_loc->Leaf_cells[i]->Flow_data.State_vec.RhoV)/(Vl->vtx_loc->Leaf_cells[i]->Flow_data.State_vec.Rhof);
#endif
		  avg += q;

		  num_vtx_cells++;
		}
	    }
  
	  fprintf(fout, "%g\n", avg/num_vtx_cells);
	}	     
      
      Vl = Vl -> next;
    }
  fprintf(fout, "</DataArray>\n");

  fprintf(fout, "<DataArray type=\"Float64\" Name=\"Velocity.z-m/s\" format=\"ascii\">\n");
  Vl = Vtxlist;
  while(Vl != NULL)
    {
      if(Vl -> vtx_loc -> vtx_nums[thread_num] != -1)
	{
	  avg = 0; num_vtx_cells = 0;
	  	
	  for(i = 0; i < MAX_NUM_CHILDREN; i++)
	    {
	      if((Vl -> vtx_loc -> Leaf_cells[i] != NULL) && (Vl -> vtx_loc -> Leaf_cells[i] -> cell_type != SOLID)
#if DETON
		 && (Vl -> vtx_loc -> Leaf_cells[i] -> un_det != TRUE)
#endif
		 )
		{
#if RECON_VEL
		  reconstruct_flow(Vl -> vtx_loc -> Leaf_cells[i], IRRELEVANT, IRRELEVANT, 7-i, 1, 'W', &q, NULL, 'n');
#else
		  q = (Vl->vtx_loc->Leaf_cells[i]->Flow_data.State_vec.RhoW)/(Vl->vtx_loc->Leaf_cells[i]->Flow_data.State_vec.Rhof);
#endif
		  avg += q;

		  num_vtx_cells++;
		}
	    }
  
	  fprintf(fout, "%g\n", avg/num_vtx_cells);
	}	     
      
      Vl = Vl -> next;
    }
  fprintf(fout, "</DataArray>\n");

  fprintf(fout, "<DataArray type=\"Float64\" Name=\"Velocity-m/s\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  Vl = Vtxlist;
  while(Vl != NULL)
    {
      if(Vl -> vtx_loc -> vtx_nums[thread_num] != -1)
	{
	  avg = 0; num_vtx_cells = 0;
 	  for(i = 0; i < MAX_NUM_CHILDREN; i++)
	    {
	      if((Vl -> vtx_loc -> Leaf_cells[i] != NULL) && (Vl -> vtx_loc -> Leaf_cells[i] -> cell_type != SOLID)
#if DETON
		 && (Vl -> vtx_loc -> Leaf_cells[i] -> un_det != TRUE)
#endif
		 )
		{
#if RECON_VEL
		  reconstruct_flow(Vl -> vtx_loc -> Leaf_cells[i], IRRELEVANT, IRRELEVANT, 7-i, 1, 'U', &q, NULL, 'n');
#else
q = (Vl->vtx_loc->Leaf_cells[i]->Flow_data.State_vec.RhoU)/(Vl->vtx_loc->Leaf_cells[i]->Flow_data.State_vec.Rhof);
#endif
		  avg += q;

		  num_vtx_cells++;
		}
	    }
	  fprintf(fout, "%g ", avg/num_vtx_cells);

	  avg = 0; 
 	  for(i = 0; i < MAX_NUM_CHILDREN; i++)
	    {
	      if((Vl -> vtx_loc -> Leaf_cells[i] != NULL) && (Vl -> vtx_loc -> Leaf_cells[i] -> cell_type != SOLID)
#if DETON
		 && (Vl -> vtx_loc -> Leaf_cells[i] -> un_det != TRUE)
#endif
		 )
		{
#if RECON_VEL
		  reconstruct_flow(Vl -> vtx_loc -> Leaf_cells[i], IRRELEVANT, IRRELEVANT, 7-i, 1, 'V', &q, NULL, 'n');
#else
q = (Vl->vtx_loc->Leaf_cells[i]->Flow_data.State_vec.RhoV)/(Vl->vtx_loc->Leaf_cells[i]->Flow_data.State_vec.Rhof);
#endif
		  avg += q;
		}		  
	    }
	  fprintf(fout, "%g ", avg/num_vtx_cells);

	  avg = 0; 
 	  for(i = 0; i < MAX_NUM_CHILDREN; i++)
	    {
	      if((Vl -> vtx_loc -> Leaf_cells[i] != NULL) && (Vl -> vtx_loc -> Leaf_cells[i] -> cell_type != SOLID)
#if DETON
		 && (Vl -> vtx_loc -> Leaf_cells[i] -> un_det != TRUE)
#endif
		 )
		{
#if RECON_VEL
		  reconstruct_flow(Vl -> vtx_loc -> Leaf_cells[i], IRRELEVANT, IRRELEVANT, 7-i, 1, 'W', &q, NULL, 'n');
#else
q = (Vl->vtx_loc->Leaf_cells[i]->Flow_data.State_vec.RhoW)/(Vl->vtx_loc->Leaf_cells[i]->Flow_data.State_vec.Rhof);
#endif
		  avg += q;
		}		
	    }
	  fprintf(fout, "%g\n", avg/num_vtx_cells);
	}	     
      
      Vl = Vl -> next;
    }
  fprintf(fout, "</DataArray>\n");

  fprintf(fout, "<DataArray type=\"Float64\" Name=\"Speed-m/s\" format=\"ascii\">\n");
  Vl = Vtxlist;
  while(Vl != NULL)
    {
      if(Vl -> vtx_loc -> vtx_nums[thread_num] != -1)
	{
	  avg = 0; num_vtx_cells = 0;
	  	
	  for(i = 0; i < MAX_NUM_CHILDREN; i++)
	    {
	      if((Vl -> vtx_loc -> Leaf_cells[i] != NULL) && (Vl -> vtx_loc -> Leaf_cells[i] -> cell_type != SOLID)
#if DETON
		 && (Vl -> vtx_loc -> Leaf_cells[i] -> un_det != TRUE)
#endif
		 )
		{
#if RECON_VEL
		  q = SQRT((SQR(Vl->vtx_loc->Leaf_cells[i]->Flow_data.State_vec.RhoU) + 
			    SQR(Vl->vtx_loc->Leaf_cells[i]->Flow_data.State_vec.RhoV) +
			    SQR(Vl->vtx_loc->Leaf_cells[i]->Flow_data.State_vec.RhoW))/
			   SQR(Vl->vtx_loc->Leaf_cells[i]->Flow_data.State_vec.Rhof));
#else
		  q = SQRT((SQR(Vl->vtx_loc->Leaf_cells[i]->Flow_data.State_vec.RhoU) + 
			    SQR(Vl->vtx_loc->Leaf_cells[i]->Flow_data.State_vec.RhoV) +
			    SQR(Vl->vtx_loc->Leaf_cells[i]->Flow_data.State_vec.RhoW))/
			   SQR(Vl->vtx_loc->Leaf_cells[i]->Flow_data.State_vec.Rhof));
#endif
		  avg += q;

		  num_vtx_cells++;
		}
	    }
  
	  fprintf(fout, "%g\n", avg/num_vtx_cells);
	}	     
      
      Vl = Vl -> next;
    }
  fprintf(fout, "</DataArray>\n");

fprintf(fout, "<DataArray type=\"Float64\" Name=\"M\" format=\"ascii\">\n");
  Vl = Vtxlist;
  while(Vl != NULL)
    {
      if(Vl -> vtx_loc -> vtx_nums[thread_num] != -1)
	{
	  avg = 0; num_vtx_cells = 0;
	  	
	  for(i = 0; i < MAX_NUM_CHILDREN; i++)
	    {
	      if((Vl -> vtx_loc -> Leaf_cells[i] != NULL) && (Vl -> vtx_loc -> Leaf_cells[i] -> cell_type != SOLID)
#if DETON
		 && (Vl -> vtx_loc -> Leaf_cells[i] -> un_det != TRUE)
#endif
		 )
		{
#if RECON_VEL
		  q = SQRT((SQR(Vl->vtx_loc->Leaf_cells[i]->Flow_data.State_vec.RhoU) + 
			    SQR(Vl->vtx_loc->Leaf_cells[i]->Flow_data.State_vec.RhoV) +
			    SQR(Vl->vtx_loc->Leaf_cells[i]->Flow_data.State_vec.RhoW))/
			   SQR(Vl->vtx_loc->Leaf_cells[i]->Flow_data.State_vec.Rhof))/
		    SQRT((gas_data.gam_amb)*(Vl->vtx_loc->Leaf_cells[i]->Flow_data.State_vec.Pres)/
			 (Vl->vtx_loc->Leaf_cells[i]->Flow_data.State_vec.Rhof));
#else
		  q = SQRT((SQR(Vl->vtx_loc->Leaf_cells[i]->Flow_data.State_vec.RhoU) + 
			    SQR(Vl->vtx_loc->Leaf_cells[i]->Flow_data.State_vec.RhoV) +
			    SQR(Vl->vtx_loc->Leaf_cells[i]->Flow_data.State_vec.RhoW))/
			   SQR(Vl->vtx_loc->Leaf_cells[i]->Flow_data.State_vec.Rhof))/
		    SQRT((gas_data.gam_amb)*(Vl->vtx_loc->Leaf_cells[i]->Flow_data.State_vec.Pres)/
			 (Vl->vtx_loc->Leaf_cells[i]->Flow_data.State_vec.Rhof));
#endif
		  avg += q;

		  num_vtx_cells++;
		}
	    }
  
	  fprintf(fout, "%g\n", avg/num_vtx_cells);
	}	     
      
      Vl = Vl -> next;
    }
  fprintf(fout, "</DataArray>\n");

#endif

#if TRY_PLOT_SHOCK
  fprintf(fout, "<DataArray type=\"Float64\" Name=\"Shocks\" format=\"ascii\">\n");
  Vl = Vtxlist;
  while(Vl != NULL)
    {
      if(Vl -> vtx_loc -> vtx_nums[thread_num] != -1)
	{	
	  avg = 0; num_vtx_cells = 0;
	
	  for(i = 0; i < MAX_NUM_CHILDREN; i++)
	    {
	      if((Vl -> vtx_loc -> Leaf_cells[i] != NULL) && (Vl -> vtx_loc -> Leaf_cells[i] -> cell_type != SOLID)
#if DETON
		 && (Vl -> vtx_loc -> Leaf_cells[i] -> un_det != TRUE)
#endif
		 )
		{
		  avg += Vl -> vtx_loc -> Leaf_cells[i] -> adapt_flag;
		  num_vtx_cells++;

#if 0
		  if(q == REFINE) /*This vertex belongs to a cell that needs refinement; I suppose we say this is a 'shock' vertex*/
		    break;
#endif
		}
	    }
  
#if 0
	  if(i < MAX_NUM_CHILDREN)
	    fprintf(fout, "1.0\n"); /*'1' is for shocks*/
	  else fprintf(fout, "0.0\n");
#else
	  fprintf(fout, "%g\n", avg/num_vtx_cells);
#endif
	}	     
      
      Vl = Vl -> next;
    }
  fprintf(fout, "</DataArray>\n");
#endif

  fprintf(fout, "</PointData>\n");
  fprintf(fout, "</Piece>\n");   
  fprintf(fout,"</UnstructuredGrid>\n");
  fprintf(fout, "</VTKFile>\n");
  fclose(fout);

  return(NOERROR);
}

/*------------------------------------------------------------------*/

/**\brief Output history files of flow variables (for the moment, just pressure) at a location - by
   traversing the root cell*/

short int write_history_loc(Cart_cell C, double loc[], double time, char filename[])
{
  FILE *fout;
  short int i;
  double p;
        
  if(C -> children[0] == NULL) {
    if((C -> cell_type != SOLID)
#if DETON
       && (C -> un_det != TRUE)
#endif
       ) { /*Use this cell's values and reconstruct to the location*/
      fout = fopen(filename, "a"); 
      if(fout == NULL) {
	printf("write_history_loc(): Can't open file \"%s\".  Exiting ...\n", filename);
	return(ERROR);
      } 

      reconstruct_flow(C, IRRELEVANT, IRRELEVANT, IRRELEVANT, IRRELEVANT, 'P', &p, loc, 'n');
      
      fprintf(fout, "%e %e\n", time, p-IC_ambient.Pres);
      fclose(fout);
    }      
  }
  else { /*Recursively traverse children*/
    for(i = 0; i < MAX_NUM_CHILDREN; i++) {
      if((C->children[i]->cell_type != SOLID) && 
#if DETON
	 (C->children[i]->un_det != TRUE) &&
#endif

#if 0
	 (C->children[i]->verticies[0]->loc[0] <= loc[0]) && (C->children[i]->verticies[0]->loc[1] <= loc[1]) && 
	 (C->children[i]->verticies[0]->loc[2] <= loc[2]) && 
	 (C->children[i]->verticies[7]->loc[0] >= loc[0]) &&
	 (C->children[i]->verticies[7]->loc[1] >= loc[1]) && (C->children[i]->verticies[7]->loc[2] >= loc[2])
#else
	 ((C->children[i]->verticies[0]->loc[0]-1e-10) <= loc[0]) && ((C->children[i]->verticies[0]->loc[1]-1e-10) <= loc[1]) && 
	 ((C->children[i]->verticies[0]->loc[2]-1e-10) <= loc[2]) && 
	 ((C->children[i]->verticies[7]->loc[0]+1e-10) >= loc[0]) &&
	 ((C->children[i]->verticies[7]->loc[1]+1e-10) >= loc[1]) && ((C->children[i]->verticies[7]->loc[2]+1e-10) >= loc[2])

	 /*Give a slight offset so that points may not align exactly with grid lines*/
#endif
	 ) { 
	/*Overlap with bounding box of cell*/
  
	if(write_history_loc(C->children[i], loc, time, filename) == ERROR)
	  return(ERROR);
	break; /*Just traverse 1 branch - it should be sufficient*/
      }
    }

    if(i == MAX_NUM_CHILDREN) {

#if 0
      printf("C level %hd at %g %g %g\n",C->cell_level,C->centroid[0],C->centroid[1],C->centroid[2]);
      printf("loc %g %g %g\n",loc[0],loc[1],loc[2]);
      for(i = 0; i < MAX_NUM_CHILDREN; i++) {
	if(C->children[i]->cell_type != SOLID) {
	  printf("child type %hd, child %hd, child level %hd, vtx0 %6.5e %6.5e %6.5e, vtx7 %6.5e %6.5e %6.5e\n",
		 C->children[i]->cell_type,i,C->children[i]->cell_level,C->children[i]->verticies[0]->loc[0],
		 C->children[i]->verticies[0]->loc[1],C->children[i]->verticies[0]->loc[2],C->children[i]->verticies[7]->loc[0],
		 C->children[i]->verticies[7]->loc[1],C->children[i]->verticies[7]->loc[2]);

	  printf("lwr logic - %d %d %d\n",(C->children[i]->verticies[0]->loc[0] <= loc[0]),
		 (C->children[i]->verticies[0]->loc[1] <= loc[1]),(C->children[i]->verticies[0]->loc[2] <= loc[2]));
	  printf("upr logic - %d %d %d\n",(C->children[i]->verticies[7]->loc[0] >= loc[0]),
		 (C->children[i]->verticies[7]->loc[1] >= loc[1]),
		 (C->children[i]->verticies[7]->loc[2] >= loc[2]));
	  if(i == 3) {
	    printf("vtx7 y co-ord is %20.19e, loc[1] is %20.19e\n",C->children[i]->verticies[7]->loc[1],loc[1]);
	    if(C->children[i]->verticies[7]->loc[1] > loc[1])
	      printf("In fact the child vtx is greater by %12.11e\n",C->children[i]->verticies[7]->loc[1]-loc[1]);
	    else if(C->children[i]->verticies[7]->loc[1] < loc[1])
	      printf("In fact the child vtx is smaller by %12.11e\n",loc[1]-C->children[i]->verticies[7]->loc[1]);
	  }
	}
      }
      abort();
#endif

      fout = fopen(filename, "a"); 
      if(fout == NULL) {
	printf("write_history_loc(): Can't open file \"%s\".  Exiting ...\n", filename);
	return(ERROR);
      } 
      fprintf(fout, "# Sorry, but no cell exists for this history location at time %e\n", time);
      fclose(fout);
    }
  }

  return(NOERROR);
}

/*------------------------------------------------------------------*/
/**\brief Output history locations at a given timestep*/

short int output_history_loc(double time)
{
  int i;
  char filename[128];
  
  for(i = 0; i < Histories.num_locs; i++) {
    strncpy(filename, "\0", 128);
    sprintf(filename, "%s%d.dat", base_hist_filename, i);

    if(write_history_loc(Root, Histories.locs[i], time, filename) == ERROR)
      return(ERROR);  
  }   

 return(NOERROR); 
}

/*------------------------------------------------------------------*/
/**\brief Finally make history files m files readable by Octave or Matlab (for easy manipulation/editing)*/

void close_history_files(void)
{
  char filename[128];
  int i;
  
  for(i = 0; i < Histories.num_locs; i++) {
    strncpy(filename, "\0", 128);
    sprintf(filename, "%s%d.dat", base_hist_filename, i);
    
#if 0
    fout = fopen(filename, "a");

    while(c != EOF)
      c = getc(fout);

    fprintf(fout, "];\n");
    fclose(fout);
#endif
  }
}

/*------------------------------------------------------------------*/

/**\brief Helper function for mesh and flow soln writing routines - if we want only boundary cells, returns the no. cells that
   need to be plotted and also checks if a FLUID cell touches a wall*/

int get_num_boundary_cells(List_leaf Head, List_leaf Tail)
{
  short int i, j;
  int num_boundary_cells = 0;
  List_leaf L = Head;
  Cart_cell C;
  double basis_quad;

  while(L != NULL)
    {
      C = L -> cell_loc;
      
#if DETON
      if(C -> un_det != TRUE)
#endif
	{
	  if(C -> cell_type == INTERSECTED)
	    num_boundary_cells++;
	  else
	    {
	      basis_quad = CALC_CELL_FLUID_AREA(C)/MAX_NUM_NEIGHBOURS;
	      L -> pure_fluid = TRUE;
	  
	      for(i = 0; i < NUM_FACES; i++)
		{
		  for(j = 0; j < MAX_NUM_NEIGHBOURS; j++)
		    {
		      if(C -> flux_area[i][j] < basis_quad) /*Obstructed area here*/
			break;
		    }
	      
		  if(j < MAX_NUM_NEIGHBOURS) /*Obstructed area found*/
		    {
		      L -> pure_fluid = FALSE;
		      break;
		    }
		}
	  
	      if(L -> pure_fluid == FALSE)
		num_boundary_cells++;
	    }
	}

      if(L == Tail)
	break;
      else L = L -> next;
    }

  return(num_boundary_cells);
}

/*------------------------------------------------------------------*/

/**\brief Read in detonation info*/

short int read_deton(char file_name[])
{
  int i, j;
  char buffer[512];
  FILE *fin = fopen(file_name, "r");
  if(fin == NULL) {
    printf("read_deton(): Couldn't open file %s.  Exiting ...\n", file_name);
    return(ERROR);
  }

  for(i = 0; i < 3; i++) /*Ignore first 3 lines*/
    fgets(buffer, sizeof(buffer), fin);

  fgets(buffer, sizeof(buffer), fin);
  if(sscanf(buffer, "No. detonation points: %d", &(Deton.num_deton_pts)) == 0) {
    printf("read_deton(): Didn't read an appropriate no. detonation points in file %s.  Exiting ...\n", file_name);
    return(ERROR);
  }

  if(Deton.num_deton_pts < 0) {
    printf("read_deton(): No. detonation points < 0.  Exiting ...\n");
    return(ERROR);
  }

  if(Deton.num_deton_pts == 0) {
    printf("read_deton(): Warning - no. detonation points is 0, so instantaneous detonation model assumed\n");
    return(NOERROR);
  }

#if DETON
  if(Deton.num_deton_pts > 0) {
    IC_level = max_refinement_level;
    Deton.deton_active = TRUE;
  }
#endif

  Deton.dpts = malloc(sizeof(double *)*(Deton.num_deton_pts));
  if(Deton.dpts == NULL) {
    printf("read_deton(): Couldn't allocate data structure.  Exiting ...\n");
    return(ERROR);
  }

  for(i = 0; i < Deton.num_deton_pts; i++) {
    Deton.dpts[i] = malloc(sizeof(double)*3);
    if(Deton.dpts[i] == NULL) {
      printf("read_deton(): Couldn't allocate array.  Exiting ...\n");
      return(ERROR);
    }
  }

  Deton.dcjvels = malloc(sizeof(double)*(Deton.num_deton_pts));
  if(Deton.dcjvels == NULL) {
    printf("read_deton(): Couldn't allocate data structure.  Exiting ...\n");
    return(ERROR);
  }

  for(i = 0; i < 4; i++) /*Ignore next 4 lines*/
    fgets(buffer, sizeof(buffer), fin);

  j = 9; /*9 lines so far*/

  for(i = 0; i < Deton.num_deton_pts; i++) {
    fgets(buffer, sizeof(buffer), fin);
    if(sscanf(buffer, "Detonation point: %lf %lf %lf", &(Deton.dpts[i][0]), &(Deton.dpts[i][1]), &(Deton.dpts[i][2])) != 3) {
      printf("read_deton(): Couldn't read in detonation point on line %d.  Exiting ...\n", j);
      return(ERROR);
    }

    j++;

    fgets(buffer, sizeof(buffer), fin);
    if(sscanf(buffer, "Detonation velocity: %lf", &(Deton.dcjvels[i])) != 1) {
      printf("read_deton(): Couldn't read in VoD on line %d.  Exiting ...\n", j);
      return(ERROR);
    }

    if(Deton.dcjvels[i] <= 0) {
      printf("read_deton(): detonation velocity no. %d <= 0.  Exiting ...\n",i);
      return(ERROR);
    }

    fgets(buffer, sizeof(buffer), fin); /*Space in between entries*/
    j++;
  }  

  return(NOERROR);
}

/*------------------------------------------------------------------*/
/**\brief Find maximum density gradient for schlieren normalization.  Can't be bothered with parallel*/

void get_max_schl(void)
{
  List_leaf L;
  Cart_cell C;
  double mag_dp;
  L = Leaves;
  
  max_schl = 0;

  while(L != NULL) {
    C = L -> cell_loc;
    
    mag_dp = SQRT(SQR(C->Flow_data.Grads.grad_rhof[0]) + SQR(C->Flow_data.Grads.grad_rhof[1]) + SQR(C->Flow_data.Grads.grad_rhof[2]));
    /*Magnitude of density gradient*/

    if(mag_dp > max_schl)
      max_schl = mag_dp;

    L = L -> next;
  }  
}

/*------------------------------------------------------------------*/
/**\brief Output the whole solution and mesh data for later solution continuation*/

void write_mesh_and_soln(Cart_cell C, int run_timesteps_so_far, double current_time)
{
  FILE *fout;
  short int i, j;
  char file_name[256], tmp[128];
  
  if(C -> children[0] != NULL) {
    /*First give info on cell level of parent and its cell no. */

    sprintf(tmp, ".%d.ov_out", run_timesteps_so_far);
    strncpy(file_name, base_soln_filename, 128);
    strcat(file_name, tmp);

    fout = fopen(file_name, "a");

    fprintf(fout, "%hd ", C -> cell_level);
    for(i = 0; i < ((C -> cell_level)+1); i++)
      fprintf(fout, "%c", C -> cell_num[i]);
    fprintf(fout, " ");

    /*Need to also know if any children are leaf children*/
    for(i = 0; i < MAX_NUM_CHILDREN; i++) {
      if(C -> children[i] -> children[0] != NULL)
	fprintf(fout, "n ");
      else fprintf(fout, "y ");
    }
    fprintf(fout, "\n");

    /*Now the next 8 lines printout properties of the children - 
    cell_type, additional_refine_times, is_small_cell, flux_areas, face_r_centroid, r_centroid, wall_area, wall_norm, cell_volume,
    shocked_cell, Rhof, RhoU, RhoV, RhoW, RhoE, Pres, T, rhof_prod, rhof_prod_ub (if applicable)
    But not all the above info is necessary for non-leaf cells (only the 1st 3)
    */
    
    for(i = 0; i < MAX_NUM_CHILDREN; i++) {
      if(C -> children[i] -> children[0] != NULL)
	fprintf(fout, "%hd %hd %hd\n", C -> children[i] -> cell_type, C -> children[i] -> additional_refine_times,
		C -> children[i] -> is_small_cell);
      else {
	fprintf(fout, "%hd %hd %hd ", C -> children[i] -> cell_type, C -> children[i] -> additional_refine_times,
		C -> children[i] -> is_small_cell);
	for(j = 0; j < 6; j++)
	  fprintf(fout, "%33.32e %33.32e %33.32e %33.32e ", C->children[i]->flux_area[j][0],C->children[i]->flux_area[j][1],
		  C->children[i]->flux_area[j][2],C->children[i]->flux_area[j][3]);
	fprintf(fout, "%33.32e %33.32e %33.32e %33.32e ", C->children[i]->face_r_centroid[0][0], C->children[i]->face_r_centroid[0][1],
		C->children[i]->face_r_centroid[1][0],C->children[i]->face_r_centroid[1][1]);
	fprintf(fout, "%33.32e %33.32e %33.32e %33.32e %33.32e %33.32e %c ", C->children[i]->r_centroid, C->children[i]->wall_area,
		C->children[i]->wall_norm[0], C->children[i]->wall_norm[1], C->children[i]->wall_norm[2], C->children[i]->cell_volume,
		C->children[i]->shocked_cell);
	fprintf(fout, "%33.32e %33.32e %33.32e %33.32e %33.32e %33.32e %33.32e %33.32e ",
		C->children[i]->Flow_data.State_vec.Rhof,C->children[i]->Flow_data.State_vec.RhoU,
		C->children[i]->Flow_data.State_vec.RhoV,C->children[i]->Flow_data.State_vec.RhoW,
		C->children[i]->Flow_data.State_vec.RhoE,C->children[i]->Flow_data.State_vec.Pres,C->children[i]->Flow_data.State_vec.T,
		C->children[i]->Flow_data.State_vec.rhof_prod);
#if AFTERBURN
	fprintf(fout, "%33.32e ", C->children[i]->Flow_data.State_vec.rhof_prod_ub);
#endif
	fprintf(fout, "\n");
      }
    }
    

    fclose(fout);

    for(i = 0; i < MAX_NUM_CHILDREN; i++)
      write_mesh_and_soln(C -> children[i], run_timesteps_so_far, current_time);
  }
}

/*------------------------------------------------------------------*/
/**\brief Read in a mesh solution to generate it.  Useful for continuing solutions after where they left off*/

short int read_mesh(char file[])
{
  FILE *fin;
  
  Cart_cell C;
  char buffer[2048], cell_num[256], shocked_cell, is_leaf[MAX_NUM_CHILDREN];
  short int cell_type, cell_level, additional_refine_times, is_small_cell;
  double flux_area[NUM_FACES][MAX_NUM_NEIGHBOURS], face_r_centroid[2][2], r_centroid, wall_area, wall_norm[3], cell_volume;
  State_vector sv; 

  int i, j, line;
  short int root_seen = FALSE;
  line = 0;
  
  fin = fopen(file, "r");
  if(fin == NULL) {
    printf("read_mesh() couldn't open file '%s'.  Exiting ...\n", file);
    return(ERROR);
  }
  
  while(1) {
    strncpy(cell_num, "\0", 256); /*Clear over it*/
    fgets(buffer, sizeof(buffer), fin);  /*1st 2 values read should be parent level, and parent cell_num*/
    if(sscanf(buffer, "%hd %s %c %c %c %c %c %c %c %c", &cell_level, cell_num, &is_leaf[0], &is_leaf[1], &is_leaf[2],
	      &is_leaf[3], &is_leaf[4], &is_leaf[5], &is_leaf[6], &is_leaf[7]) != 10) {
      printf("read_mesh() expected to see cell_level, cell_num and is_leaf data on line %d.  Exiting ...\n",line+1);
      break;
    }
    
    line++;
    if((root_seen == TRUE) && (cell_level == 0))
      break;
    if(cell_level == 0)
      root_seen = TRUE;

    /*Refine the parent cell (must find it first)*/
    C = NULL;
    refine_specific(Root, cell_level, cell_num, &C);
    if(C == NULL) {
      printf("read_mesh() had trouble constructing the grid.  Exiting ...\n");
      return(ERROR);
    }

    for(i = 0; i < MAX_NUM_CHILDREN; i++) {
      fgets(buffer, sizeof(buffer), fin);  
      line++;

      /*Next 8 lines l each should contain for child l -
	cell_type, additional_refine_times, is_small_cell, flux_areas, face_r_centroid, r_centroid, wall_area, wall_norm, cell_volume,
	shocked_cell, Rhof, RhoU, RhoV, RhoW, RhoE, Pres, T, rhof_prod (no rhof_prod_ub for now)
	But not all the above info is necessary for non-leaf cells (only the 1st 3)
       */
      /***rhof_prod_ub needs fixing - just involving 1 extra data point to read in but I can't be bothered**/

	if(is_leaf[i] == 'n') {
	  if(sscanf(buffer, "%hd %hd %hd", &cell_type, &additional_refine_times, &is_small_cell) != 3) {
	    printf("read_mesh() encountered error with file '%s' on line %d\n", file, line);
	    return(ERROR);
	  }
	}
	else {
	  if(sscanf(buffer, "%hd %hd %hd %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %c %lf %lf %lf %lf %lf %lf %lf %lf", &cell_type, &additional_refine_times, &is_small_cell,
		    &flux_area[0][0], &flux_area[0][1], &flux_area[0][2], &flux_area[0][3],
		    &flux_area[1][0], &flux_area[1][1], &flux_area[1][2], &flux_area[1][3],
		    &flux_area[2][0], &flux_area[2][1], &flux_area[2][2], &flux_area[2][3],
		    &flux_area[3][0], &flux_area[3][1], &flux_area[3][2], &flux_area[3][3],
		    &flux_area[4][0], &flux_area[4][1], &flux_area[4][2], &flux_area[4][3],
		    &flux_area[5][0], &flux_area[5][1], &flux_area[5][2], &flux_area[5][3],
		    &face_r_centroid[0][0], &face_r_centroid[0][1], &face_r_centroid[1][0], &face_r_centroid[1][1],
		    &r_centroid, &wall_area, &wall_norm[0], &wall_norm[1], &wall_norm[2], &cell_volume, &shocked_cell,
		    &(sv.Rhof), &(sv.RhoU), &(sv.RhoV), &(sv.RhoW), &(sv.RhoE), &(sv.Pres), &(sv.T), &(sv.rhof_prod)) != 46) { 
	    /*There're 46 data points to read in (unfortunately)*/
	    printf("read_mesh() encountered error with file '%s' on line %d\n", file, line);
	    return(ERROR);
	  }
   
	  for(j = 0; j < NUM_FACES; j++) {
	    C -> children[i] -> flux_area[j][0] = flux_area[j][0];
	    C -> children[i] -> flux_area[j][1] = flux_area[j][1];
	    C -> children[i] -> flux_area[j][2] = flux_area[j][2];
	    C -> children[i] -> flux_area[j][3] = flux_area[j][3];
	  }

	  C -> children[i] -> face_r_centroid[0][0] = face_r_centroid[0][0]; 
	  C -> children[i] -> face_r_centroid[0][1] = face_r_centroid[0][1];
	  C -> children[i] -> face_r_centroid[1][0] = face_r_centroid[1][0];
	  C -> children[i] -> face_r_centroid[1][1] = face_r_centroid[1][1];
	  C -> children[i] -> r_centroid = r_centroid;
	
	  C -> children[i] -> wall_area = wall_area;
	  C -> children[i] -> wall_norm[0] = wall_norm[0]; C -> children[i] -> wall_norm[1] = wall_norm[1];
	  C -> children[i] -> wall_norm[2] = wall_norm[2];
	  C -> children[i] -> cell_volume = cell_volume;

	  C -> children[i] -> shocked_cell = shocked_cell;
	
	  C -> children[i] -> Flow_data.State_vec.Rhof = sv.Rhof;
	  C -> children[i] -> Flow_data.State_vec.RhoU = sv.RhoU;
	  C -> children[i] -> Flow_data.State_vec.RhoV = sv.RhoV;
	  C -> children[i] -> Flow_data.State_vec.RhoW = sv.RhoW;
	  C -> children[i] -> Flow_data.State_vec.RhoE = sv.RhoE;
	  C -> children[i] -> Flow_data.State_vec.Pres = sv.Pres;
	  C -> children[i] -> Flow_data.State_vec.T = sv.T;
	  C -> children[i] -> Flow_data.State_vec.rhof_prod = sv.rhof_prod;

	  C->children[i]->Flow_data.Grads.grad_rhof[0] = 0.0; C->children[i]->Flow_data.Grads.grad_rhof[1] = 0.0; 
	  C->children[i]->Flow_data.Grads.grad_rhof[2] = 0.0;
	  C->children[i]->Flow_data.Grads.grad_E[0] = 0.0; C->children[i]->Flow_data.Grads.grad_E[1] = 0.0; 
	  C->children[i]->Flow_data.Grads.grad_E[2] = 0.0;
	  C->children[i]->Flow_data.Grads.grad_u[0] = 0.0; C->children[i]->Flow_data.Grads.grad_u[1] = 0.0; 
	  C->children[i]->Flow_data.Grads.grad_u[2] = 0.0;
	  C->children[i]->Flow_data.Grads.grad_v[0] = 0.0; C->children[i]->Flow_data.Grads.grad_v[1] = 0.0; 
	  C->children[i]->Flow_data.Grads.grad_v[2] = 0.0;
	  C->children[i]->Flow_data.Grads.grad_w[0] = 0.0; C->children[i]->Flow_data.Grads.grad_w[1] = 0.0; 
	  C->children[i]->Flow_data.Grads.grad_w[2] = 0.0;  
	}

      /*Now initialize these values*/
      C -> children[i] -> cell_type = cell_type;
      C -> children[i] -> additional_refine_times = additional_refine_times;
      C -> children[i] -> is_small_cell = is_small_cell;

      /*We need to do this step finally (normally done in refine())*/
      if((C -> children[i] -> Leaf_list_loc != NULL) && (cell_type == SOLID)) 
	delete_from_leaf_list(C -> children[i], &Leaves, &Leaves_tail);
    }
  }

  return(NOERROR);
}


/*------------------------------------------------------------------*/
