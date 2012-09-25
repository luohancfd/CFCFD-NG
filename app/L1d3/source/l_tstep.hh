// l_tstep.hh

#ifndef L_TSTEP_HH
#define L_TSTEP_HH

int L_encode_conserved_for_one_cell(struct L_cell* c);
int L_encode_conserved(struct slug_data* A);
int L_decode_conserved(struct slug_data* A);
int L_set_chemistry_timestep(struct slug_data* A, double dt);
int L_set_thermal_timestep(struct slug_data* A, double dt);
int L_chemical_increment(struct slug_data* A, double dt);
int L_source_vector(struct slug_data* A);
int L_axial_heat_flux(struct slug_data* A, double k);
int L_adjust_end_cells(struct slug_data* A);
int L_apply_rivp(struct slug_data* A);
int L_time_derivatives(struct slug_data* A, int time_level);
int L_record_slug_state(struct slug_data* A);
int L_restore_slug_state(struct slug_data* A);
int L_predictor_step(struct slug_data* A);
int L_corrector_step(struct slug_data* A);
int L_check_cells(struct slug_data* A, int js);
int L_check_cfl(struct slug_data* A);
int P_time_derivatives(struct piston_data* B, int time_level, double sim_time);
int P_record_piston_state(struct piston_data* B);
int P_restore_piston_state(struct piston_data* B);
int P_predictor_step(struct piston_data* B);
int P_corrector_step(struct piston_data* B);

#endif
