#include <Python.h>
#include "numpy/arrayobject.h"
#include "sginfo.h"

PyObject *exception = NULL;

PyArrayObject* SeitzToPython(T_RTMx *seitz, int nValues){
  npy_intp twod[] = {nValues, 4, 4};
  float lastrow[] = {0.0, 0.0, 0.0, 1.0};
  float *p;  
  PyArrayObject* SeitzMx;

  SeitzMx = (PyArrayObject *)PyArray_SimpleNew(3, twod, NPY_FLOAT);

  for(twod[0]=0; twod[0] < nValues; twod[0]++){
    for(twod[1]=0; twod[1] < 3; twod[1]++){
      for(twod[2]=0; twod[2] < 3; twod[2]++){
	p = PyArray_GetPtr(SeitzMx, twod); 
	*p = (float)(seitz->s.R[(twod[1] * 3) + twod[2]]);
      }
    }
    
    twod[2] = 3;
    for(twod[1]=0; twod[1] < 3; twod[1]++){
      p = PyArray_GetPtr(SeitzMx, twod); 
      *p = (float)(seitz->s.T[twod[1]]) / 12.0;
    }

    twod[1] = 3;
    for(twod[2]=0; twod[2]<4; twod[2]++){
      p = PyArray_GetPtr(SeitzMx, twod); 
      *p = lastrow[twod[2]];
    }
    
    seitz++;
  }

  return SeitzMx;
}

PyArrayObject* TrVectorToPython(const T_LatticeInfo *lattice, int nValues){
  npy_intp twod[] = {nValues, 3};
  float *p;  
  PyArrayObject* TrVector;

  TrVector = (PyArrayObject *)PyArray_SimpleNew(2, twod, PyArray_FLOAT);
  for(twod[0] = 0; twod[0] < nValues; twod[0]++){
    for(twod[1] = 0; twod[1] < 3; twod[1]++){
      p = PyArray_GetPtr(TrVector, twod); 
      *p = (float)(lattice->TrVector[(twod[0] * 3) + twod[1]]) / 12.0;
    }
  }

  return TrVector;
}

int parse_sg(T_SgInfo* sg, const char* name){
  const T_TabSgName  *tsgn;
  const char *Name = NULL;
  
  tsgn = NULL;
  tsgn = FindTabSgNameEntry(name, 'A');
  if(tsgn != NULL)
    Name = tsgn->HallSymbol;
  
  sg->MaxList = 192; /* absolute maximum number of symops */
  sg->ListSeitzMx
    = malloc(sg->MaxList * sizeof (*sg->ListSeitzMx));
  
  if (sg->ListSeitzMx == NULL) {
    SetSgError("Not enough core");
    return 0;
  }

  sg->ListRotMxInfo
    = malloc(sg->MaxList * sizeof (*sg->ListRotMxInfo));
  
  if (sg->ListRotMxInfo == NULL) {
    SetSgError("Not enough core");
    return 0;
  }
  
  /* Initialize the SgInfo structure
   */
  
  InitSgInfo(sg);
  sg->TabSgName = tsgn; /* in case we know the table entry */
  
  /* Translate the Hall symbol and generate the whole group
   */
  
  if(Name != NULL){
    ParseHallSymbol(Name, sg);
    if (SgError != NULL) return -2;
  
    /* Do some book-keeping and derive crystal system, point group,
       and - if not already set - find the entry in the internal
       table of space group symbols
    */
  
    CompleteSgInfo(sg);
    return -1;
  } else {
    return -2;
  }
}

static PyObject *
sginfo_getsg(PyObject *self, PyObject *args)
{
  const char *command;
  int iList;
  int rtn;
  T_SgInfo  SgInfo;

  if (!PyArg_ParseTuple(args, "s", &command))
    return NULL;
  rtn = parse_sg(&SgInfo, command);
  if(rtn == 0){
    PyErr_SetString(PyExc_RuntimeError, "Unable to initialize SGINFO.");
    return NULL;
  } else if(rtn == -2) {
    PyErr_SetString(exception, "Unable to parse spacegroup.");
    return NULL;
  } 

  iList = PG_Index(SgInfo.PointGroup);
  
  PyObject *tSgInfo;
  PyObject *SgName;
  PyArrayObject *SeitzMx;
  PyArrayObject *TrVector;
  
  SeitzMx = SeitzToPython(SgInfo.ListSeitzMx, SgInfo.nList);
  TrVector = TrVectorToPython(SgInfo.LatticeInfo, SgInfo.LatticeInfo->nTrVector);
  
  SgName = Py_BuildValue("{s:i, s:s, s:s, s:s}",
			 "Number", SgInfo.TabSgName->SgNumber, 
			 "Label", SgInfo.TabSgName->SgLabels,
			 "Hall", SgInfo.TabSgName->HallSymbol,
			 "Setting", SgInfo.TabSgName->Extension);

  tSgInfo = Py_BuildValue("{s:i,s:s,s:s,s:s,s:s,s:i,s:i,s:i,s:c}", 
			  "Centric", SgInfo.Centric, 
			  "PointGroup", PG_Names[iList],
			  "LaueGroup", PG_Names[PG_Index(LG_Code_of_PG_Index[iList])],
			  "System", XS_Name[SgInfo.XtalSystem],
			  "Extra", EI_Name[SgInfo.ExtraInfo],
			  "InversionOffOrigin", SgInfo.InversionOffOrigin,
			  "OrderL", SgInfo.OrderL,
			  "OrderP", SgInfo.OrderP,
			  "LatticeCode", SgInfo.LatticeInfo->Code);

  /* ListSgInfo(&SgInfo, 1, 0, stderr); */
  return Py_BuildValue("(O,O,O,O)", SgName, tSgInfo, TrVector, SeitzMx);
}

static PyMethodDef SginfoMethods[] = {
    {"getsg",  sginfo_getsg, METH_VARARGS, "Execute a shell command."},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

PyMODINIT_FUNC
initsginfo(void)
{
    PyObject *m;

    m = Py_InitModule("sginfo", SginfoMethods);
    if (m == NULL)
        return;

    import_array();

    exception = PyErr_NewException("sginfo.error", NULL, NULL);
    Py_INCREF(exception);
    PyModule_AddObject(m, "error", exception);
}
