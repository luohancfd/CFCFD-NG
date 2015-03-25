// Author: Kyle Damm and Rowan Gollan
// Date: 22-Jan-2015

#include "gpu-chem-update.hh"


#define MAX_SOURCE_SIZE (0x100000)


using namespace std;

#include "kernel.hh"
#include "../../../lib/gas/kinetics/chemical-kinetic-ODE-update.hh"
#include "../../../lib/gas/kinetics/chemical-kinetic-system.hh"



gpu_chem* init_gpu_module(int ncells) 
{

  // Load the kernel source code into the array programSource                                                                                                                 
    FILE *fp;
    char *programSource;
    size_t source_size;

    fp = fopen("/home/qpsf01/kdamm/cfcfd3/app/eilmer3/source/alpha_qss_kernel.cl", "r");
    if (fp ==  NULL) {
        fprintf(stderr, "Failed to load kernel.\n");
        exit(1);
    }
    programSource = (char*)malloc(MAX_SOURCE_SIZE);
    source_size = fread( programSource, 1, MAX_SOURCE_SIZE, fp);
    fclose( fp );



  gpu_chem* gchem = new gpu_chem();
    opencl_container& oc = gchem->cl;
    Gas_model *gmodel = get_gas_model_ptr();
    Chemical_kinetic_ODE_update *ckupdate = dynamic_cast<Chemical_kinetic_ODE_update*>(get_reaction_update_ptr());
    Chemical_kinetic_system *cks = ckupdate->get_cks_pointer();
    
    gchem->nspec = gmodel->get_number_of_species();
    gchem->nreac = cks->get_n_reactions();


    int numcells = ncells;
    int numdebug = 0;
    const int nreac_elements = numcells*gchem->nreac;
    const int nspec_elements = numcells*gchem->nspec;
    oc.nreac_datasize=sizeof(double)*nreac_elements;
    oc.datasize=sizeof(double)*nspec_elements;
    oc.ncell = sizeof(double)*numcells;
    oc.debugnum = sizeof(double)*numdebug;

    std::cout << "number of cells: " << ncells << "\n";
    std::cout << "number of species: " << gchem->nspec << "\n";
    std::cout << "number of reactions: " << gchem->nreac << "\n";


    
    gpu_work_arrays& w = gchem->w_arrays;

    w.kf = (double*)malloc(oc.nreac_datasize);
    w.kb = (double*)malloc(oc.nreac_datasize);
    w.MM = (double*)malloc(oc.nreac_datasize);
    w.Y = (double*)malloc(oc.datasize);
    w.Yc = (double*)malloc(oc.datasize);
    w.Yp = (double*)malloc(oc.datasize);
    w.Ypo = (double*)malloc(oc.datasize);
    w.h = (double*)malloc(oc.ncell);
    w.alpha = (double*)malloc(oc.datasize);
    w.q = (double*)malloc(oc.datasize);
    w.p = (double*)malloc(oc.datasize);
    w.qp = (double*)malloc(oc.datasize);
    w.pp = (double*)malloc(oc.datasize);
    w.q_bar = (double*)malloc(oc.datasize);
    w.p_bar = (double*)malloc(oc.datasize);
    w.debugging = (double*)malloc(oc.debugnum);

    
    //Find available devices----------------------------------------------------------------------------
    //retrieve number of platforms
    oc.status=clGetPlatformIDs(0,NULL, &(oc.numPlatforms));
 
    //allocate enough space for each platform
    oc.platforms=(cl_platform_id*)malloc(oc.numPlatforms*sizeof(cl_platform_id));
 
    //fill in platforms
    oc.status = clGetPlatformIDs(oc.numPlatforms, oc.platforms, NULL);
 
    //retrieve number of devices, note for GPU change CPU to GPU
    oc.status=clGetDeviceIDs(oc.platforms[0], CL_DEVICE_TYPE_GPU, 0, NULL, &(oc.numDevices));
 
    //allocate enough space for each device
    oc.devices=(cl_device_id*)malloc(oc.numDevices*sizeof(cl_device_id));
 
    //fill in the devices, for for GPU change CPU to GPU
    oc.status=clGetDeviceIDs(oc.platforms[0], CL_DEVICE_TYPE_GPU, oc.numDevices, oc.devices, NULL);

    //D. Create Context and Command Queue for each device for use---------------------------------------------
 
    oc.context = clCreateContext(NULL, oc.numDevices, oc.devices, NULL, NULL, &(oc.status));
 
    oc.cmdQueue=clCreateCommandQueue(oc.context, oc.devices[0], 0, &(oc.status));
 
    //E. Create buffers which will contain data on device-----------------------------------------------------

    oc.bufkf=clCreateBuffer(oc.context, CL_MEM_READ_ONLY, oc.nreac_datasize, NULL, &(oc.status));
    oc.bufkb=clCreateBuffer(oc.context, CL_MEM_READ_ONLY, oc.nreac_datasize, NULL, &(oc.status));
    oc.bufMM=clCreateBuffer(oc.context, CL_MEM_READ_WRITE, oc.nreac_datasize, NULL, &(oc.status));
    oc.bufY=clCreateBuffer(oc.context, CL_MEM_READ_WRITE, oc.datasize, NULL, &(oc.status));
    oc.bufYp=clCreateBuffer(oc.context, CL_MEM_READ_WRITE, oc.datasize, NULL, &(oc.status));
    oc.bufYpo=clCreateBuffer(oc.context, CL_MEM_READ_WRITE, oc.datasize, NULL, &(oc.status));
    oc.bufYc=clCreateBuffer(oc.context, CL_MEM_READ_WRITE, oc.datasize, NULL, &(oc.status));
    oc.bufh=clCreateBuffer(oc.context, CL_MEM_READ_WRITE, oc.ncell, NULL, &(oc.status));
    oc.bufalpha=clCreateBuffer(oc.context, CL_MEM_READ_WRITE, oc.datasize, NULL, &(oc.status));
    oc.bufq=clCreateBuffer(oc.context, CL_MEM_READ_WRITE, oc.datasize, NULL, &(oc.status));
    oc.bufp=clCreateBuffer(oc.context, CL_MEM_READ_WRITE, oc.datasize, NULL, &(oc.status));
    oc.bufqp=clCreateBuffer(oc.context, CL_MEM_READ_WRITE, oc.datasize, NULL, &(oc.status));
    oc.bufpp=clCreateBuffer(oc.context, CL_MEM_READ_WRITE, oc.datasize, NULL, &(oc.status));
    oc.bufq_bar=clCreateBuffer(oc.context, CL_MEM_READ_WRITE, oc.datasize, NULL, &(oc.status));
    oc.bufp_bar=clCreateBuffer(oc.context, CL_MEM_READ_WRITE, oc.datasize, NULL, &(oc.status));
    oc.bufalpha_bar=clCreateBuffer(oc.context, CL_MEM_READ_WRITE, oc.datasize, NULL, &(oc.status));
    oc.bufdebugging=clCreateBuffer(oc.context, CL_MEM_READ_WRITE, oc.debugnum, NULL, &(oc.status));

    //create a program with source code
    oc.program=clCreateProgramWithSource(oc.context, 1, (const char**)&programSource,  (const size_t *)&source_size, &(oc.status));
     
    //build (compile) the program for the device
    oc.status = clBuildProgram(oc.program, oc.numDevices, oc.devices, NULL, NULL, NULL);
     
    //create the kernel
    oc.kernel = clCreateKernel(oc.program, "new_qss", &(oc.status));
    
    return gchem;
}

void clear_gpu_module(gpu_chem* gchem)
{
    gpu_work_arrays& w = gchem->w_arrays;
    free(w.kf);
    free(w.kb);
    free(w.MM);
    free(w.Y);
    free(w.Yp);
    free(w.Ypo);
    free(w.Yc);
    free(w.h);
    free(w.alpha);
    free(w.q);
    free(w.p);
    free(w.qp);
    free(w.pp);
    free(w.q_bar);
    free(w.p_bar);
    free(w.alpha_bar);
    free(w.debugging);

    opencl_container& oc = gchem->cl;
    clReleaseCommandQueue(oc.cmdQueue);
    clReleaseMemObject(oc.bufkf);
    clReleaseMemObject(oc.bufkb);
    clReleaseMemObject(oc.bufMM);
    clReleaseMemObject(oc.bufY);
    clReleaseMemObject(oc.bufYp);
    clReleaseMemObject(oc.bufYpo);
    clReleaseMemObject(oc.bufYc);
    clReleaseMemObject(oc.bufh);
    clReleaseMemObject(oc.bufalpha);
    clReleaseMemObject(oc.bufq);
    clReleaseMemObject(oc.bufp);
    clReleaseMemObject(oc.bufqp);
    clReleaseMemObject(oc.bufpp);
    clReleaseMemObject(oc.bufq_bar);
    clReleaseMemObject(oc.bufp_bar);
    clReleaseMemObject(oc.bufalpha_bar);
    clReleaseMemObject(oc.bufdebugging);
    clReleaseContext(oc.context);
    clReleaseKernel(oc.kernel);
    clReleaseProgram(oc.program);
    
    delete gchem;
}

int update_chemistry(gpu_chem &gchem, double dt_global, vector<FV_Cell*>& cells)
{
 
    // HARD-CODED ARRAY SIZES
    // FIXME: when generalising for different GPUs
    // end FIXME
    global_data &G = *get_global_data_ptr();
    Gas_model *gmodel = get_gas_model_ptr();
    Chemical_kinetic_ODE_update *ckupdate = dynamic_cast<Chemical_kinetic_ODE_update*>(get_reaction_update_ptr());
    Chemical_kinetic_system *cks = ckupdate->get_cks_pointer();
    vector<double> conc;
    conc.resize(gchem.nspec);
    gpu_work_arrays& w = gchem.w_arrays;
    int cellnum = cells.size();

    for ( size_t i = 0; i < cells.size(); ++i ) {
	FV_Cell *cp = cells[i];
	convert_massf2conc(cp->fs->gas->rho, cp->fs->gas->massf, gmodel->M(), conc);
	// In this loop, we transfer the arrangement of cell-stored mass fraction data
	// to the large vector storage of data for the GPU. In actuality, we convert
	// to concentration rather than mass fraction for use in the ODE solver.
	for ( int isp = 0; isp < gchem.nspec; ++isp ) {
	  //  w.Y[i*gchem.nspec+isp] = conc[isp] ;
	  w.Y[i+isp*cellnum] = conc[isp]; 
	}
	
	// In this loop, we populate the large kf and kb vectors for the GPU
	// using the lib/gas routines to compute reaction rate coefficients.
	for ( int ir = 0; ir < gchem.nreac; ++ir ) {
	    Reaction* reac = cks->get_reaction(ir);
	    reac->compute_k_f(*(cp->fs->gas));
	    reac->compute_k_b(*(cp->fs->gas));
	    // w.kf[i*gchem.nreac+ir] = reac->k_f();
	    // w.kb[i*gchem.nreac+ir] = reac->k_b();
	    w.kf[i+ir*cellnum] = reac->k_f();
	    w.kb[i+ir*cellnum] = reac->k_b();

	}

	// Populate the stepsize (h) array
	if ( cp->dt_chem <= 0.0 ) {
	    cks->set_gas_data_ptr(*(cp->fs->gas));
	    w.h[i] = cks->stepsize_select(conc);
	}
	else {
	    w.h[i] = cp->dt_chem;
	}

    }

    //    std::cout << "Species Concentration Check (before kernel): " << w.Y[0] << "\t";
    //    std::cout << "kf (before kernel): " << w.kf[31] << "\t";
    //    std::cout << "kb (before kernel): " << w.kb[31] << "\t";

    //    std::cout << "Y[0] =   " << w.h[0] << "\n";

    //F. Write reaction rates to device buffers----------------------------------------------------------------
    opencl_container& oc = gchem.cl;
    oc.status=clEnqueueWriteBuffer(oc.cmdQueue, oc.bufkf, CL_FALSE, 0, oc.nreac_datasize,(void *)w.kf, 0, NULL, NULL);
    oc.status=clEnqueueWriteBuffer(oc.cmdQueue, oc.bufkb, CL_FALSE, 0, oc.nreac_datasize,(void *)w.kb, 0, NULL, NULL);
    oc.status=clEnqueueWriteBuffer(oc.cmdQueue, oc.bufY, CL_FALSE, 0, oc.datasize,(void *)w.Y, 0, NULL, NULL);
    oc.status=clEnqueueWriteBuffer(oc.cmdQueue, oc.bufh, CL_FALSE, 0, oc.ncell,(void*)w.h, 0, NULL, NULL);
    /*
    std::cout << "Y[0]b =   " << w.Y[0] << "\n";
    std::cout << "Y[1]b =   " << w.Y[4096] << "\n";
    std::cout << "Y[2]b =   " << w.Y[8192] << "\n";
    std::cout << "Y[3]b =   " << w.Y[12288] << "\n";
    std::cout << "Y[4]b =   " << w.Y[16384] << "\n";
    std::cout << "Y[5]b =   " << w.Y[20480] << "\n";
    std::cout << "Y[6]b =   " << w.Y[24576] << "\n";
    */
    //set kernel arguments
    oc.status = clSetKernelArg(oc.kernel, 0, sizeof(cl_mem),(void *)&(oc.bufkf));
    oc.status = clSetKernelArg(oc.kernel, 1, sizeof(cl_mem),(void *)&(oc.bufkb));
    oc.status = clSetKernelArg(oc.kernel, 2, sizeof(cl_mem),(void *)&(oc.bufMM));
    oc.status = clSetKernelArg(oc.kernel, 3, sizeof(cl_mem),(void *)&(oc.bufY));
    oc.status = clSetKernelArg(oc.kernel, 4, sizeof(cl_mem),(void *)&(oc.bufYc));
    oc.status = clSetKernelArg(oc.kernel, 5, sizeof(cl_mem),(void *)&(oc.bufYp));
    oc.status = clSetKernelArg(oc.kernel, 6, sizeof(cl_mem),(void *)&(oc.bufh));
    oc.status = clSetKernelArg(oc.kernel, 7, sizeof(cl_mem),(void *)&(oc.bufalpha));
    oc.status = clSetKernelArg(oc.kernel, 8, sizeof(cl_mem),(void *)&(oc.bufq));
    oc.status = clSetKernelArg(oc.kernel, 9, sizeof(cl_mem),(void *)&(oc.bufp));
    oc.status = clSetKernelArg(oc.kernel, 10, sizeof(cl_mem),(void *)&(oc.bufqp));
    oc.status = clSetKernelArg(oc.kernel, 11, sizeof(cl_mem),(void *)&(oc.bufpp));
    oc.status = clSetKernelArg(oc.kernel, 12, sizeof(cl_mem),(void *)&(oc.bufq_bar));
    oc.status = clSetKernelArg(oc.kernel, 13, sizeof(cl_mem),(void *)&(oc.bufp_bar));
    oc.status = clSetKernelArg(oc.kernel, 14, sizeof(cl_mem),(void *)&(oc.bufalpha_bar));
    oc.status = clSetKernelArg(oc.kernel, 15, sizeof(cl_mem),(void *)&(oc.bufYpo));
    oc.status = clSetKernelArg(oc.kernel, 16, sizeof(cl_int),(void *)&(gchem.nspec));
    oc.status = clSetKernelArg(oc.kernel, 17, sizeof(cl_int),(void *)&(gchem.nreac));
    oc.status = clSetKernelArg(oc.kernel, 18, sizeof(cl_int),(void *)&(cellnum));
    oc.status = clSetKernelArg(oc.kernel, 19, sizeof(cl_double),(void *)&(dt_global));
    oc.status = clSetKernelArg(oc.kernel, 20, sizeof(cl_mem),(void *)&(oc.bufdebugging));


       
    //define an index space (global work size) of work items for execution
    size_t globalWorkSize[1];
    size_t localWorkSize[1]; 

    globalWorkSize[0]=cells.size();
    //localWorkSize[0]= 64;

    //execute the kernel
    oc.status = clEnqueueNDRangeKernel(oc.cmdQueue, oc.kernel, 1, NULL, globalWorkSize, NULL, 0 , NULL, NULL);
     
    //read output off device back to CPU
    clEnqueueReadBuffer(oc.cmdQueue, oc.bufY, CL_TRUE, 0, oc.datasize,(void *)w.Y, 0, NULL, NULL);
    //    clEnqueueReadBuffer(oc.cmdQueue, oc.bufdebugging, CL_TRUE, 0, oc.debugnum,(void *)w.debugging, 0, NULL, NULL);
    clEnqueueReadBuffer(oc.cmdQueue, oc.bufh, CL_TRUE, 0, oc.ncell, (void *)w.h, 0, NULL, NULL);
    /*
    std::cout << "Y[0]a =   " << w.Y[0] << "\n";
    std::cout << "Y[1]a =   " << w.Y[4096] << "\n";
    std::cout << "Y[2]a =   " << w.Y[8192] << "\n";
    std::cout << "Y[3]a =   " << w.Y[12288] << "\n";
    std::cout << "Y[4]a =   " << w.Y[16384] << "\n";
    std::cout << "Y[5]a =   " << w.Y[20480] << "\n";
    std::cout << "Y[6]a =   " << w.Y[24576] << "\n";
    */
    //    std::cout << "Y[0] =   " << w.Y[0] << "\n";
  
    //    std::cout << "Species Concentration Check (after kernel): " << w.Y[0] << "\n";
    //    std::cout << "debug value: " << w.debugging[0] << "\n";

    // Now reverse step 1. Get concentrations from large Y vector and populate cells with mass fractions.
    for ( size_t i = 0; i < cells.size(); ++i ) {
	FV_Cell *cp = cells[i];
	for ( int isp = 0; isp < gchem.nspec; ++isp ) {
	    conc[isp] = w.Y[i+isp*cellnum];
	}
	convert_conc2massf(cp->fs->gas->rho, conc, gmodel->M(), cp->fs->gas->massf);
	// Store a copy of the chemical step size.
	cp->dt_chem = w.h[i];
	// Enforce thermodynamic constraint of fixed mass, fixed energy.
	gmodel->eval_thermo_state_rhoe(*(cp->fs->gas));
	// If we are doing a viscous sim, we'll need to ensure
	// viscous properties are up-to-date
	if ( G.viscous ) gmodel->eval_transport_coefficients(*(cp->fs->gas));
	if ( G.diffusion ) gmodel->eval_diffusion_coefficients(*(cp->fs->gas));
	// ...but we have to manually update the conservation quantities
	// for the gas-dynamics time integration.
	// Species densities: mass of species isp per unit volume.
	for ( int isp = 0; isp < gchem.nspec; ++isp )
	    cp->U[0]->massf[isp] = cp->fs->gas->rho * cp->fs->gas->massf[isp];
    }

      
    return SUCCESS;
}
