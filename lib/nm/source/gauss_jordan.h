void gaussj(struct cell_center *cell, int n, int m);
void inverse_product(struct cell_center *cell);
void matrix_inverse(struct cell_center *cell, int n);
void gaussj_pivoting(double **a, int n, double **b, int m);
void nrerror(char error_text[]);
int *ivector(long nl, long nh);
void free_ivector(int *v);

