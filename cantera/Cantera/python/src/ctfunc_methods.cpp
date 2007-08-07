
static PyObject*
py_func_new(PyObject *self, PyObject *args)
{
    int type, n;
    PyObject* c;
    if (!PyArg_ParseTuple(args, "iiO:func_new", &type, &n, &c))
        return NULL;
    PyArrayObject* coeffs = (PyArrayObject*)c;
    double* xd = (double*)coeffs->data;
    int lenc = coeffs->dimensions[0];
    int nn = func_new(type, n, lenc, xd);
    if (nn < 0) return reportError(nn);
    return Py_BuildValue("i",nn);
}

static PyObject*
py_func_newcombo(PyObject *self, PyObject *args)
{
    int type, n, m;
    if (!PyArg_ParseTuple(args, "iii:func_newcombo", &type, &n, &m))
        return NULL;
    int nn = func_new(type, n, m, 0);
    if (nn < 0) return reportError(nn);
    return Py_BuildValue("i",nn);
}

static PyObject*
py_func_del(PyObject *self, PyObject *args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:func_del", &n))
        return NULL;
    int iok = func_del(n);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",0);
}

static PyObject*
py_func_value(PyObject *self, PyObject *args)
{
    int n;
    double t;
    if (!PyArg_ParseTuple(args, "id:func_value", &n, &t))
        return NULL;
    double r = func_value(n, t);
    return Py_BuildValue("d",r);
}

