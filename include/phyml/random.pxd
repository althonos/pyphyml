from libc.stdint cimport uint32_t


cdef extern from "random.h" nogil:

    cdef struct __MTRand:
        pass
    ctypedef __MTRand MTRand

    void MTRand_Seed(MTRand* rand, uint32_t seed)