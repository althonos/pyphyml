from .utilities cimport option as t_option

cdef extern from "cl.h" nogil:
    bint Read_Command_Line(t_option* input, int argc, char** argv)