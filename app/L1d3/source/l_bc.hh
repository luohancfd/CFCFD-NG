// l_bc.hh

#ifndef L_BC_HH
#define L_BC_HH

int L_bc_left_velocity(struct slug_data *A, double v);
int L_bc_left_reflect(struct slug_data *A);
int L_bc_left_free(struct slug_data *A);
int L_bc_right_velocity(struct slug_data *A, double v);
int L_bc_right_reflect(struct slug_data *A);
int L_bc_right_free(struct slug_data *A);
int L_exchange_bc_data(struct slug_data *A, struct slug_data *B);
int L_apply_bc(struct slug_data *A);

#endif
