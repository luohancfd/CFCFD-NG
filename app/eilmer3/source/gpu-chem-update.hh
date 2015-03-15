// Author: Kyle Damm and Rowan Gollan
// Date: 22-Jan-2015

#ifndef GPU_CHEM_UPDATE_HH
#define GPU_CHEM_UPDATE_HH

#include <vector>
#include <CL/opencl.h>
#include "cell.hh"

struct gpu_work_arrays {
    double *kf, *kb, *MM;
    double *Y, *Yp, *Ypo, *Yc;
    double *h;
    double *alpha, *q, *p, *qp, *pp;
    double *q_bar, *p_bar, *alpha_bar;
    double *cond, *sigma, *debugging;
};

struct opencl_container {
    cl_int status;
    cl_uint numPlatforms;
    cl_platform_id* platforms;
    cl_uint numDevices;
    cl_device_id* devices;
    cl_context context;
    cl_command_queue cmdQueue;

    cl_mem bufkf, bufkb, bufMM;
    cl_mem bufY, bufYp, bufYpo, bufYc;
    cl_mem bufh;
    cl_mem bufalpha, bufq, bufp, bufqp, bufpp;
    cl_mem bufq_bar, bufp_bar, bufalpha_bar;
    cl_mem bufcond, bufsigma, bufdebugging;

    cl_program program;
    cl_kernel kernel;

    size_t nreac_datasize;
    size_t datasize;
    size_t ncell;
    size_t debugnum;
};

struct gpu_chem {
    int nspec, nreac, numcells;
    gpu_work_arrays w_arrays;
    opencl_container cl; 
};
     
gpu_chem* init_gpu_module(int ncells);
void clear_gpu_module(gpu_chem *gchem);
int update_chemistry(gpu_chem &gchem, double dt_global, std::vector<FV_Cell*>& cells);

#endif
