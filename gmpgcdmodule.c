#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <gmp.h>

static PyObject *gmpgcdError;



int _asdf(mpz_t pr[], int length) {

/*	mpz_t a;
	mpz_init(a);
	mpz_set_ui(a,1);
*/
	unsigned long int b;
	b=2;
	mpz_t temp;
	long out;
	long mult;
	int index;

	for (index = 0; index < length; index++)
		out = mpz_get_si(pr[index]);
		printf("pr[%d] = %ld\n", index, out);
		mpz_mul_si(temp,pr[index],b);
		mult = mpz_get_si(temp);
		printf("pr[%d] * 2 = %ld\n", index, mult);

	return 0;
}






static PyObject* gmpgcd(PyObject* self, PyObject *args){

/*	unsigned long int b;
	unsigned long int a;
	unsigned long int out;

	mpz_t aa;
	mpz_t bb;
	mpz_t d;
	mpz_init(aa);
	mpz_init(bb);
	mpz_init(d);
*/
	PyObject *long_list;
	int pr_length;
	mpz_t *pr;
	int index;

	if (!PyArg_ParseTuple(args, "O", &long_list))
		return NULL;
	pr_length = PyObject_Length(long_list);
	if (pr_length < 0)
		return NULL;
	pr = (mpz_t *) malloc(sizeof(mpz_t *) * pr_length);
	if (pr == NULL)
		return NULL;
	for (index = 0; index < pr_length; index++) {
		PyObject *item;
		item = PyList_GetItem(long_list, index);
		mpz_init(pr[index]);
		if (!PyLong_Check(item))
			mpz_set_ui(pr[index], 0);
		mpz_set_ui(pr[index], PyLong_AsLong(item));
	}
	return Py_BuildValue("i", _asdf(pr, pr_length));



/*	mpz_set_ui(aa,a);
	mpz_set_ui(bb,b);



	mpz_add(d, aa, bb);

	out = mpz_get_ui(d);

	return PyLong_FromLong(out);



	mpz_clear(aa);
	mpz_clear(d);
*/
}

static PyMethodDef gmpgcdMethods[] = {
	{"gmpgcd", gmpgcd, METH_VARARGS, "compute gcd of two integers using GMP lib"},
  {NULL,NULL,0,NULL} /*Sentinel*/
};

static struct PyModuleDef gmpgcdmodule = {
	PyModuleDef_HEAD_INIT,
	"gmpgcd",   /* name of module */
	NULL, /* module documentation, may be NULL */
	-1,       /* size of per-interpreter state of the module,
				 or -1 if the module keeps state in global variables. */
	gmpgcdMethods
};

PyMODINIT_FUNC PyInit_gmpgcd(void){
	PyObject *m;
	m = PyModule_Create(&gmpgcdmodule);
	if(m == NULL)
		return NULL;

	gmpgcdError = PyErr_NewException("gcd.error", NULL, NULL);
	Py_XINCREF(gmpgcdError);

	if (PyModule_AddObject(m, "error", gmpgcdError) < 0) {
		Py_XDECREF(gmpgcdError);
		Py_CLEAR(gmpgcdError);
		Py_DECREF(m);
		return NULL;
	}
	return m;
}

int
main(int argc, char *argv[])
{
	wchar_t *program = Py_DecodeLocale(argv[0], NULL);
	if (program == NULL) {
		fprintf(stderr, "Fatal error: cannot decode argv[0]\n");
		exit(1);
	}

	/* Add a built-in module, before Py_Initialize */
	PyImport_AppendInittab("gmpgcd", PyInit_gmpgcd);

	/* Pass argv[0] to the Python interpreter */
	Py_SetProgramName(program);

	/* Initialize the Python interpreter.  Required. */
	Py_Initialize();

	/* Optionally import the module; alternatively,
	   import can be deferred until the embedded script
	   imports it. */
	PyImport_ImportModule("gmpgcd");


	PyMem_RawFree(program);
	return 0;
}






