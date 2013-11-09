//
// Python/C interface file for basic subroutines
//

#include <Python.h>
#include "numpy/arrayobject.h"

#include <cfloat>
#include <cstring>
#include <iostream>
#include "trm/subs.h"

// Voigt profile routine

static PyObject* 
subs_voigt(PyObject *self, PyObject *args)
{
    double a;
    PyObject *vo = NULL;
    double eps = 1.e-8;
    if(!PyArg_ParseTuple(args, "dO|d:voigt", &a, &vo, &eps))
	return NULL;

    if(!vo){
	PyErr_SetString(PyExc_ValueError, "subs.voigt: v not defined");
	return NULL;
    } 

    double v;
    if(PyFloat_Check(vo)){
	v = PyFloat_AS_DOUBLE(vo);
	double f = Subs::voigt(a,v,eps);
	return Py_BuildValue("f", f);

    }else if(PyArray_Check(vo)){
	int nd = PyArray_NDIM(vo);
	if(nd != 1){
	    PyErr_SetString(PyExc_ValueError, "subs.voigt: v must be a 1D array or float");
	    return NULL;
	}

	npy_intp nelem = PyArray_Size(vo);
	PyObject *arr = PyArray_FROM_OTF(vo, NPY_DOUBLE, NPY_IN_ARRAY);
	if(arr == NULL) return NULL;
	double *vptr = (double *)PyArray_DATA(arr);

	npy_intp dim[1] = {nelem};
	PyArrayObject *oarr = (PyArrayObject*) PyArray_SimpleNew(1, dim, PyArray_DOUBLE);
	if(oarr == NULL){
	    Py_DECREF(arr);
	    return NULL;
	}
	double *optr = (double*) oarr->data;

	for(npy_intp i=0; i<nelem; i++)
	    optr[i] = Subs::voigt(a,vptr[i],eps);

	Py_DECREF(arr);
	return PyArray_Return(oarr);

    }else{
	PyErr_SetString(PyExc_ValueError, "subs.voigt: v must be a 1D array or float");
	return NULL;
    }
};

// Incomplete gamma function

static PyObject* 
subs_gammq(PyObject *self, PyObject *args)
{
    double a;
    PyObject *xo = NULL;
    if(!PyArg_ParseTuple(args, "dO:gammq", &a, &xo))
	return NULL;

    if(!xo){
	PyErr_SetString(PyExc_ValueError, "subs.gammq: x not defined");
	return NULL;
    } 

    double x;
    if(PyFloat_Check(xo)){
	x = PyFloat_AS_DOUBLE(xo);
	double f = Subs::gammq(a,x);
	return Py_BuildValue("f", f);

    }else if(PyArray_Check(xo)){
	int nd = PyArray_NDIM(xo);
	if(nd != 1){
	    PyErr_SetString(PyExc_ValueError, "subs.gammq: x must be a 1D array or float");
	    return NULL;
	}

	npy_intp nelem = PyArray_Size(xo);
	PyObject *arr = PyArray_FROM_OTF(xo, NPY_DOUBLE, NPY_IN_ARRAY);
	if(arr == NULL) return NULL;
	double *xptr = (double *)PyArray_DATA(arr);

	npy_intp dim[1] = {nelem};
	PyArrayObject *oarr = (PyArrayObject*) PyArray_SimpleNew(1, dim, PyArray_DOUBLE);
	if(oarr == NULL){
	    Py_DECREF(arr);
	    return NULL;
	}
	double *optr = (double*) oarr->data;

	for(npy_intp i=0; i<nelem; i++)
	    optr[i] = Subs::gammq(a,xptr[i]);

	Py_DECREF(arr);
	return PyArray_Return(oarr);

    }else{
	PyErr_SetString(PyExc_ValueError, "subs.gammq: x must be a 1D array or float");
	return NULL;
    }
};

//----------------------------------------------------------------------------------------
// The methods

static PyMethodDef SubsMethods[] = {

    {"voigt", subs_voigt, METH_VARARGS, 
     "fv = voigt(a,v,eps=1.e-8) -- returns Voigt function. v can be an array or a float"},

    {"gammq", subs_gammq, METH_VARARGS, 
     "g = gammq(a,x) -- returns incomplete gamma function. x can be an array or a float"},

    {NULL, NULL, 0, NULL} /* Sentinel */
};

PyMODINIT_FUNC
init_subs(void)
{
    (void) Py_InitModule("_subs", SubsMethods);
    import_array();
}
