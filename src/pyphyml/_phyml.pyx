from cpython.unicode cimport PyUnicode_FromString

from libc.stdlib cimport malloc, calloc, free, srand
from libc.string cimport strdup, strndup, strcmp
from libc.time cimport time, time_t
from libc.float cimport DECIMAL_DIG
from libc.stdio cimport (
    printf as PhyML_Printf,
    fprintf as PhyML_Fprintf,
    fflush,
    rewind,
)

cimport phyml.utilities
cimport phyml.make
cimport phyml.models
cimport phyml.cl
cimport phyml.io
cimport phyml.init
cimport phyml.free
cimport phyml.lk
cimport phyml.spr
cimport phyml.optimiz
cimport phyml.pars
cimport phyml.ancestral

from phyml.utilities cimport (
    phydbl,
    align as t_align,
    option as t_option,
    t_mod,
    t_opt,
    t_tree,
    t_topo,
    t_type,
)

import os

# --- Alignment ----------------------------------------------------------------

cdef class Alignment:
    """A multiple sequence alignment.
    """

    cdef t_align** _data
    cdef int       _n_otu
    cdef int       _datatype

    def __cinit__(self):
        self._data = NULL
        self._n_otu = -1
        self._datatype = t_type.UNDEFINED

    def __init__(
        self,
        object names not None,
        object sequences not None,
        *,
        bint protein = False,
    ):
        cdef int i
        cdef int n = len(sequences)

        cdef str name
        cdef str seq

        #
        self._datatype = <int> t_type.AA if protein else <int> t_type.NT

        if len(names) != n:
            raise ValueError("dimension mismatch")

        self._allocate(n)
        self._n_otu = n

        for i, (name, seq) in enumerate(zip(names, sequences)):
            # copy name
            self._data[i].name = strndup(name.encode('utf-8'), len(name)) # FIXME: allow bytes
            if self._data[i].name is NULL:
                raise MemoryError
            # copy seq
            self._data[i].len = len(seq)
            self._data[i].state = strndup(seq.encode('utf-8'), len(seq)) # FIXME: allow bytes
            if self._data[i] is NULL:
                raise MemoryError
            # check length consistency
            if self._data[i].len != self._data[0].len:
                raise ValueError("dimension mismatch")

        if n > 0:
            self._post_process()

    def __dealloc__(self):
        phyml.free.Free_Seq(self._data, self._n_otu)

    cdef void _allocate(self, int n) except *:

        self._data = <t_align**> calloc(n, sizeof(t_align*))
        if self._data is NULL:
            raise MemoryError

        for i in range(n):
            self._data[i] = <t_align*> calloc(1, sizeof(t_align))
            if self._data[i] is NULL:
                raise MemoryError

    cdef void _post_process(self) except *:
        # NOTE: Adapted from `Post_Process_Data` in `io.c` to avoid requiring
        #       a `t_option*` argument. This function ensures the sequence data
        #       is using a consistent indeterminate symbol and that sequences
        #       are sorted by sequence name.

        cdef int      i
        cdef int      j
        cdef bint     swap
        cdef t_align* data_buff

        for j in range(self._n_otu):
            for i in range(self._data[0].len):
                if self._data[j].state[i] == ord('*') or self._data[j].state[i] == ord('?') or self._data[j].state[i] == ord('-'):
                    self._data[j].state[i] = ord('X')
                elif self._datatype == t_type.NT and self._data[j].state[i] == ord('N'):
                    self._data[j].state[i] = ord('X')
                elif self._data[j].state[i] == ord('U'):
                    self._data[j].state[i] = ord('T')

        # FIXME: maybe not needed?
        for i in range(self._n_otu):
            self._data[i].len = self._data[0].len

        # Sequences are to be ordered alphabetically
        data_buff = NULL
        swap      = True
        # FIXME: use qsort instead here?
        while swap:
            swap = False
            for i in range(self._n_otu):
                for j in range(i+1, self._n_otu):
                    if strcmp(self._data[i].name, self._data[j].name) < 0:
                        swap        = True
                        data_buff   = self._data[i]
                        self._data[i] = self._data[j]
                        self._data[j] = data_buff


# --- Tree ---------------------------------------------------------------------

cdef class Tree:
    """A phylogenetic tree.
    """
    cdef t_tree* _tree

    def __cinit__(self):
        self._tree = NULL

    def __init__(self):
        raise TypeError("cannot instantiate Tree")

    def __dealloc__(self):
        if self._tree is not NULL:
            phyml.free.Free_Tree(self._tree)

    cpdef str dumps(self):
        assert self._tree is not NULL

        cdef char* data
        cdef str   s

        data = phyml.io.Write_Tree(self._tree)
        try:
            s = PyUnicode_FromString(data)
        finally:
            free(data)

        return s



# --- Main using alignment -----------------------------------------------------

cdef t_option* Get_Input_Default():
    cdef t_option* io     = NULL
    cdef t_mod*    mod    = NULL
    cdef t_opt*    s_opt  = NULL
    cdef int       rv     = 1

    io    = phyml.make.Make_Input()
    mod   = phyml.make.Make_Model_Basic()
    s_opt = phyml.make.Make_Optimiz()

    phyml.utilities.Set_Defaults_Input(io)
    phyml.utilities.Set_Defaults_Model(mod)
    phyml.utilities.Set_Defaults_Optimiz(s_opt)

    io.mod    = mod
    mod.io    = io
    mod.s_opt = s_opt

    # switch (argc)
    # {
    # case 1:
    # {
    #     Launch_Interface(io);
    #     break;
    # }
    # default:
    # {
    # rv = phyml.cl.Read_Command_Line(io, argc, argv);
    # }
    # }
    #endif
    # if not rv:
    #     phyml.free.Free_Optimiz(s_opt)
    #     phyml.free.Free_Model(mod)
    #     phyml.free.Free_Input(io)
    #     return NULL

    return io

def main_ali(
    Alignment ali not None,
    *,
    bint debug = False,
):

    # cdef calign *cdata;
    cdef t_option* io                       = NULL
    cdef t_tree*   tree                     = NULL
    cdef t_tree*   dum                      = NULL
    cdef int       num_data_set
    cdef int       num_tree
    cdef int       num_rand_tree
    # cdef t_mod  *mod;
    cdef time_t    time_beg
    cdef time_t    time_end
    cdef phydbl    best_lnL                 = phyml.utilities.UNLIKELY
    cdef int       r_seed
    cdef char*     s
    cdef char*     most_likely_tree         = NULL
    cdef bint      orig_random_input_tree

    #ifdef QUIET
    # setvbuf(stdout, NULL, _IOFBF, 2048);
    #endif

    io = Get_Input_Default()
    if io is NULL:
        raise ValueError("Failed parsing input")

    #
    io.quiet = not debug

    # Seed global RNG
    # FIXME: see if the global RNG can be replaced with a local type
    # r_seed = time(NULL) if io.r_seed < 0 else io.r_seed
    r_seed = 42
    srand(r_seed)
    io.r_seed = r_seed

    #
    if io.in_tree == 2:
        # Test_Multiple_Data_Set_Format(io)
        raise NotImplementedError("Test_Multiple_Data_Set_Format")
    else:
        io.n_trees = 1

    #
    if io.n_trees == 0 and io.in_tree == 2:
        raise ValueError("invalid format of input tree")

    #
    if io.n_data_sets > 1 and io.n_trees > 1:
        io.n_data_sets = min(io.n_trees, io.n_data_sets)
        io.n_trees     = min(io.n_trees, io.n_data_sets)

    #
    for num_data_set in range(0, io.n_data_sets):

        best_lnL = phyml.utilities.UNLIKELY

        # phyml.io.Get_Seq(io)
        io.n_otu = ali._n_otu
        io.data = ali._data
        io.datatype = ali._datatype


        phyml.make.Make_Model_Complete(io.mod)
        phyml.utilities.Set_Model_Name(io.mod)
        phyml.io.Print_Settings(io)

        mod                    = io.mod
        orig_random_input_tree = io.mod.s_opt.random_input_tree

        if io.data is not NULL:

            if io.n_data_sets > 1:
                PhyML_Printf("\n. Data set [#%d]\n", num_data_set + 1)

            cdata = phyml.utilities.Compact_Data(io.data, io)

            # phyml.free.Free_Seq(io.data, cdata.n_otu)

            for num_tree in range(0 if io.n_trees == 1 else num_data_set, io.n_trees):
                if not io.mod.s_opt.random_input_tree:
                    io.mod.s_opt.n_rand_starts = 1

                if orig_random_input_tree and io.n_trees > 1:
                    PhyML_Printf("\n== Cannot combine random starting trees with multiple input trees.")
                    raise RuntimeError("EXIT")

                for num_rand_tree in range(io.mod.s_opt.n_rand_starts):
                # {
                    if (io.mod.s_opt.random_input_tree) and (io.mod.s_opt.topo_search != t_topo.NNI_MOVE):
                        if not io.quiet:
                            PhyML_Printf("\n\n. [Random start %3d/%3d]", num_rand_tree + 1, io.mod.s_opt.n_rand_starts)

                    phyml.init.Init_Model(cdata, mod, io)
                    phyml.models.Set_Model_Parameters(mod)

                    # #ifdef M4
                    # if (io->mod->use_m4mod) M4_Init_Model(mod->m4mod, cdata, mod);
                    # #endif

                    # Make the initial tree
                    if io.in_tree == 0 or io.in_tree == 1:
                        tree = phyml.utilities.Dist_And_BioNJ(cdata, mod, io)
                        if io.print_mat_and_exit:
                            phyml.io.Print_Mat(phyml.lk.ML_Dist(cdata, mod))
                            raise RuntimeError("exit", -1)
                    elif io.in_tree == 2:
                        tree = phyml.io.Read_User_Tree(cdata, mod, io)
                        assert tree is not NULL # FIXME?

                    if io.mod.s_opt.opt_topo:
                        phyml.utilities.Remove_Duplicates(cdata, io, tree)

                    if io.fp_in_constraint_tree is not NULL:

                        PhyML_Printf("\n. Reading constraint tree file...")

                        io.cstr_tree = phyml.io.Read_Tree_File_Phylip(io.fp_in_constraint_tree)
                        if io.cstr_tree.n_root is not NULL:
                            PhyML_Printf("\n== The constraint tree file must be unrooted")
                            raise RuntimeError("exit")

                        s = phyml.utilities.Add_Taxa_To_Constraint_Tree(io.fp_in_constraint_tree, cdata)
                        fflush(NULL)
                        phyml.free.Free_Tree(tree)
                        tree = phyml.io.Read_Tree(&s)
                        io.in_tree = 2
                        phyml.free.Free(s)
                        s = NULL

                        phyml.utilities.Check_Constraint_Tree_Taxa_Names(io.cstr_tree, cdata)
                        phyml.utilities.Alloc_Bip(io.cstr_tree)
                        phyml.utilities.Get_Bip(io.cstr_tree.a_nodes[0], io.cstr_tree.a_nodes[0].v[0], io.cstr_tree)
                        if not tree.has_branch_lengths:
                            phyml.utilities.Add_BioNJ_Branch_Lengths(tree, cdata, mod, NULL)

                    if tree is NULL:
                        continue

                    time(&time_beg)
                    time(&(tree.t_beg))

                    tree.mod          = mod
                    tree.io           = io
                    tree.data         = cdata
                    tree.n_root       = NULL
                    tree.e_root       = NULL
                    tree.n_tot_bl_opt = 0

                    phyml.utilities.Set_Both_Sides(True, tree)

                    # check memory requirement on first run
                    if num_data_set == 0 and num_tree == 0 and num_rand_tree == 0:
                        phyml.utilities.Check_Memory_Amount(tree)

                    if io.cstr_tree is not NULL and not phyml.utilities.Check_Topo_Constraints(tree, io.cstr_tree):
                        PhyML_Printf("\n\n== The initial tree does not satisfy the topological constraint.")
                        PhyML_Printf("\n== Please use the user input tree option with an adequate tree topology.")
                        raise RuntimeError("EXIT")

                    phyml.utilities.Connect_CSeqs_To_Nodes(tree.data, tree.io, tree)
                    phyml.make.Make_Tree_For_Pars(tree)
                    phyml.make.Make_Tree_For_Lk(tree)
                    phyml.make.Make_Spr(tree)
                    phyml.utilities.Br_Len_Not_Involving_Invar(tree)
                    phyml.utilities.Unscale_Br_Len_Multiplier_Tree(tree)

                    # #ifdef BEAGLE
                    # if (mod->bootstrap == YES)
                    # {
                    #     PhyML_Printf("\n== PhyML-BEAGLE does not support bootstrap "
                    #                 "analysis yet... ");
                    #     Exit("\n");
                    # }
                    # if (mod->ras->invar == YES)
                    # {
                    #     PhyML_Printf("\n== PhyML-BEAGLE does not support invariant site "
                    #                 "models yet... ");
                    #     Exit("\n");
                    # }
                    # #endif

                    if tree.io.print_json_trace:
                        phyml.io.JSON_Tree_Io(tree, tree.io.fp_out_json_trace)

                    phyml.utilities.Set_Update_Eigen(True, tree.mod)
                    phyml.lk.Lk(NULL, tree)
                    phyml.utilities.Set_Update_Eigen(False, tree.mod)

                    PhyML_Printf("\n. Init log-likelihood: %f", tree.c_lnL)

                    if (tree.mod.s_opt.opt_topo):
                        phyml.spr.Global_Spr_Search(tree)
                        if (tree.n_root):
                            phyml.utilities.Add_Root(tree.a_edges[0], tree)
                    else:
                        # #ifdef BEAGLE
                        # tree->b_inst = create_beagle_instance(tree, io->quiet, io);
                        # #endif
                        if tree.mod.s_opt.opt_subst_param or tree.mod.s_opt.opt_bl_one_by_one:
                            phyml.optimiz.Round_Optimize(tree, phyml.utilities.ROUND_MAX)

                    # if(tree->mod->gamma_mgf_bl)
                    # Best_Root_Position_IL_Model(tree);

                    phyml.utilities.Set_Both_Sides(True, tree)
                    phyml.lk.Lk(NULL, tree)
                    phyml.pars.Pars(NULL, tree)
                    phyml.utilities.Get_Tree_Size(tree)
                    PhyML_Printf("\n\n. Log likelihood of the current tree: %.*f.", DECIMAL_DIG, tree.c_lnL)

                    if tree.io.ancestral:
                        phyml.ancestral.Ancestral_Sequences(tree, True)

                    phyml.utilities.Check_Br_Lens(tree)
                    phyml.utilities.Br_Len_Involving_Invar(tree)
                    phyml.utilities.Rescale_Br_Len_Multiplier_Tree(tree)

                    if tree.n_root is NULL:
                        phyml.utilities.Get_Best_Root_Position(tree)

                    # Print the tree estimated using the current random (or BioNJ) starting tree
                    # if(io->mod->s_opt->n_rand_starts > 1)
                    if orig_random_input_tree:
                        phyml.io.Print_Tree(io.fp_out_trees, tree)
                        fflush(NULL)

                    # Record the most likely tree in a string of characters
                    # FIXME: avoid serialization here?
                    if tree.c_lnL > best_lnL:
                        best_lnL = tree.c_lnL
                        if most_likely_tree is not NULL:
                            phyml.free.Free(most_likely_tree)
                        most_likely_tree = phyml.io.Write_Tree(tree)

                        time(&time_end)

                        # FIXME: we don't have FP out but maybe we could
                        #        capture the stats here rather than outputting
                        #        them
                        # phyml.io.Print_Fp_Out(
                        #     io.fp_out_stats,
                        #     time_beg,
                        #     time_end,
                        #     tree,
                        #     io,
                        #     num_data_set + 1,
                        #     num_rand_tree if orig_random_input_tree else num_tree,
                        #     num_rand_tree == io.mod.s_opt.n_rand_starts - 1,
                        #     io.precision
                        # )

                        if tree.io.print_site_lnl:
                            phyml.io.Print_Site_Lk(tree, io.fp_out_lk)

                    # Start from BioNJ tree
                    if num_rand_tree == io.mod.s_opt.n_rand_starts - 1 and tree.mod.s_opt.random_input_tree:
                        # Do one more iteration in the loop, but don't randomize the tree
                        tree.mod.s_opt.n_rand_starts += 1
                        tree.mod.s_opt.random_input_tree = False

                    # #ifdef BEAGLE
                    # finalize_beagle_instance(tree);
                    # #endif
                    phyml.free.Free_Best_Spr(tree)
                    phyml.free.Free_Spr_List_One_Edge(tree)
                    phyml.free.Free_Spr_List_All_Edge(tree)
                    phyml.free.Free_Tree_Pars(tree)
                    phyml.free.Free_Tree_Lk(tree)
                    phyml.free.Free_Tree(tree)
                # Tree done

                if io.n_data_sets == 1:
                    if io.fp_out_tree is not NULL:
                        rewind(io.fp_out_tree)

                if most_likely_tree is not NULL:
                    if io.fp_out_tree is not NULL:
                        PhyML_Fprintf(io.fp_out_tree, "%s\n", most_likely_tree)

                # Launch bootstrap analysis
                if io.do_boot or io.do_tbe:
                    if not io.quiet:
                        PhyML_Printf("\n\n. Launch bootstrap analysis on the most likely tree...")
                    most_likely_tree = phyml.utilities.Bootstrap_From_String(most_likely_tree, cdata, mod, io)
                    PhyML_Printf("\n\n. Completed the bootstrap analysis succesfully.")
                    fflush(NULL)
                elif io.ratio_test != 0.0:
                    # Launch aLRT
                    most_likely_tree = phyml.utilities.aLRT_From_String(most_likely_tree, cdata, mod, io)

                # Print the most likely tree in the output file
                if not io.quiet:
                    PhyML_Printf("\n\n. Printing the most likely tree in file '%s'.", phyml.utilities.Basename(io.out_tree_file));
                if io.n_data_sets == 1:
                    if io.fp_out_tree is not NULL:
                        rewind(io.fp_out_tree)

                dum      = phyml.io.Read_Tree(&most_likely_tree)
                dum.data = cdata
                dum.mod  = mod
                dum.io   = io
                phyml.utilities.Connect_CSeqs_To_Nodes(cdata, io, dum)
                phyml.utilities.Insert_Duplicates(dum)
                phyml.free.Free(most_likely_tree)
                most_likely_tree = phyml.io.Write_Tree(dum)
                phyml.free.Free_Tree(dum)

                if io.fp_out_tree is not NULL:
                    PhyML_Fprintf(io.fp_out_tree, "%s\n", most_likely_tree)

                if io.n_trees > 1 and io.n_data_sets > 1:
                    break

            phyml.free.Free_Calign(cdata)

        else:
            PhyML_Printf("\n== No data was found.\n")
            raise RuntimeError("Exit")

        phyml.free.Free_Model_Complete(mod)

    # if most_likely_tree is not NULL:
        # phyml.free.Free(most_likely_tree)

    if mod.s_opt.n_rand_starts > 1:
        PhyML_Printf("\n. Best log likelihood: %f\n", best_lnL)

    phyml.free.Free_Optimiz(mod.s_opt)
    phyml.free.Free_Model_Basic(mod)

    # if (io->fp_in_constraint_tree) fclose(io->fp_in_constraint_tree);
    # if (io->fp_in_align) fclose(io->fp_in_align);
    # if (io->fp_in_tree) fclose(io->fp_in_tree);
    # if (io->fp_out_lk) fclose(io->fp_out_lk);
    # if (io->fp_out_tree) fclose(io->fp_out_tree);
    # if (io->fp_out_trees) fclose(io->fp_out_trees);
    # if (io->fp_out_stats) fclose(io->fp_out_stats);
    # if (io->fp_out_trace) fclose(io->fp_out_trace);
    # if (io->fp_out_json_trace) fclose(io->fp_out_json_trace);

    # if (io->fp_in_constraint_tree != NULL) Free_Tree(io->cstr_tree);
    phyml.free.Free_Input(io)

    time(&time_end)
    phyml.io.Print_Time_Info(time_beg, time_end)

    if most_likely_tree is NULL:
        return None

    # FIXME: read, ocnnect cseqs, insert duplciates, free
    # return most_likely_tree.decode('utf-8')

    cdef Tree t = Tree.__new__(Tree)
    t._tree = phyml.io.Read_Tree(&most_likely_tree)
    phyml.free.Free(most_likely_tree)
    return t
