from cpython.unicode cimport PyUnicode_FromString

from libc.stdlib cimport malloc, calloc, free, srand
from libc.string cimport strdup, strcmp
from libc.time cimport time, time_t
from libc.float cimport DECIMAL_DIG
from libc.math cimport NAN
from libc.stdio cimport (
    printf as PhyML_Printf,
    fprintf as PhyML_Fprintf,
    fflush,
    rewind,
)

cimport phyml.utilities
cimport phyml.make
cimport phyml.models
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
    calign as t_calign,
    option as t_option,
    t_mod,
    t_opt,
    t_tree,
    t_topo,
    t_type,
    t_whichmodel,
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
            self._data[i].name = strdup(name.encode('utf-8')) # FIXME: allow bytes
            if self._data[i].name is NULL:
                raise MemoryError
            # copy seq
            self._data[i].len = len(seq)
            self._data[i].state = strdup(seq.encode('utf-8')) # FIXME: allow bytes
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

# --- CompressedAlignment ------------------------------------------------------

cdef class CompressedAlignment:
    cdef t_calign* _data

    def __dealloc__(self):
        if self._data is not NULL:
            phyml.free.Free_Calign(self._data)

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


# --- Model --------------------------------------------------------------------

cdef dict _NT_MODELS = {
    "JC69": <int> t_whichmodel.JC69,
    "K80": <int> t_whichmodel.K80,
    "F81": <int> t_whichmodel.F81,
    "HKY85": <int> t_whichmodel.HKY85,
    "F84": <int> t_whichmodel.F84,
    "TN93": <int> t_whichmodel.TN93,
    "GTR": <int> t_whichmodel.GTR,
}

cdef dict _AA_MODELS = {
    "WAG": <int> t_whichmodel.WAG,
}

cdef class ModelPrototype:
    cdef t_mod* _mod
    cdef int    _type

    def __cinit__(self):
        self._mod = NULL
        self._type = <int> t_type.UNDEFINED

    def __init__(self):
        raise TypeError("cannot instantiate a Model")

    def __dealloc__(self):
        if self._mod is not NULL:
            phyml.free.Free_Model_Complete(self._mod)
            phyml.free.Free_Model_Basic(self._mod)

    @classmethod
    def from_name(
        cls,
        str name,
        *,
        int rate_categories = 4,
        object kappa = None,
        object lambda_ = None,
        object alpha = 1.0,
    ):
        """Create a new substitution model prototype.
        """
        cdef ModelPrototype model

        # validate parameters
        if rate_categories < 1:
            raise ValueError("number of rate categories must be a positive integer")
        if name not in _NT_MODELS and name not in _AA_MODELS:
            raise ValueError(f"unknown model: {name!r}")
        if alpha is not None and alpha < 1e-10:
            raise ValueError(f"alpha must be >=1E-10 (got {alpha!r})")
        if kappa is not None and kappa < 0:
            raise ValueError(f"kappa must be positive (got {kappa!r})")

        # allocate data
        model = ModelPrototype.__new__(ModelPrototype)
        model._mod = phyml.make.Make_Model_Basic()
        if model._mod is NULL:
            raise MemoryError

        # initialize PhyML defaults
        phyml.utilities.Set_Defaults_Model(model._mod)

        # set model kind and name
        if name in _NT_MODELS:
            model._type = <int> t_type.NT
            model._mod.ns = 4
            model._mod.whichmodel = <int> _NT_MODELS[name]
        elif name in _AA_MODELS:
            raise NotImplementedError
        else:
            raise ValueError(name)
        phyml.utilities.Set_Model_Name(model._mod)

        # set alpha value
        if alpha is None:
            model._mod.ras.alpha.v = 1.0
            model._mod.ras.alpha.optimize = True
        else:
            model._mod.ras.alpha.v = alpha
            model._mod.ras.alpha.optimize = False

        # set lambda and kappa values
        if model._mod.whichmodel != t_whichmodel.JC69 and model._mod.whichmodel != t_whichmodel.F81 and model._mod.whichmodel != t_whichmodel.GTR:
            if kappa is None:
                model._mod.kappa.v = 4.0
                model._mod.kappa.optimize = True
                if model._mod.whichmodel == t_whichmodel.TN93:
                    model._mod.lambda_.optimize = True
            else:
                model._mod.kappa.optimize = False
                model._mod.lambda_.optimize = False
                if model._mod.whichmodel == t_whichmodel.TN93:
                    raise NotImplementedError
                else:
                    model._mod.kappa.v = kappa
        if model._mod.whichmodel != t_whichmodel.K80 and model._mod.whichmodel != t_whichmodel.HKY85 and model._mod.whichmodel != t_whichmodel.F84 and model._mod.whichmodel != t_whichmodel.TN93:
            model._mod.kappa.optimize = False

        # set number of rate categories
        model._mod.ras.n_catg = rate_categories
        if rate_categories == 1:
            model._mod.ras.alpha.optimize = False

        # finish initialization
        phyml.make.Make_Model_Complete(model._mod)
        return model


# --- Main using alignment -----------------------------------------------------


cdef class Result:
    cdef readonly Tree                tree
    cdef readonly Alignment           alignment
    cdef readonly CompressedAlignment compressed
    cdef readonly phydbl              log_likelihood

    def __init__(
        self,
        Tree tree not None,
        Alignment alignment not None,
        CompressedAlignment compressed not None,
        phydbl log_likelihood,
    ):
        self.tree = tree
        self.alignment = alignment
        self.compressed = compressed
        self.log_likelihood = log_likelihood


cdef class TreeBuilder:
    cdef readonly ModelPrototype model
    cdef readonly int            seed

    def __init__(
        self,
        ModelPrototype model = None,
        *,
        int seed = 0,
    ):
        """Create a new `TreeBuilder` with the given parameters.

        Keyword Arguments:
            seed (`int`): The seed to initialize the random number generator
                with. If a negative number is given, the seed will be
                initialized using the system clock.
            model (`ModelPrototype` or `None`): The substitution model
                parameters to use. Pass `None` to use the default *HKY85*
                (for nucleotide) or *LG* (for protein) models with the
                transition/transversion ratio estimated by maximum likelihood.

        """
        self.seed = seed
        self.model = model or ModelPrototype.from_name("HKY85")

    # ---

    cdef t_option* _create_default_options(self) except NULL:
        cdef t_option* io     = NULL
        cdef t_mod*    mod    = NULL
        cdef t_opt*    s_opt  = NULL
        cdef int       rv     = 1

        io = phyml.make.Make_Input()
        if io is NULL:
            raise MemoryError

        mod = phyml.utilities.Copy_Model(self.model._mod)
        if mod is NULL:
            raise MemoryError

        s_opt = phyml.make.Make_Optimiz()
        if s_opt is NULL:
            raise MemoryError

        phyml.utilities.Set_Defaults_Input(io)
        phyml.utilities.Set_Defaults_Optimiz(s_opt)

        io.mod    = mod
        mod.io    = io
        mod.s_opt = s_opt

        return io

    # FIXME: see if the global RNG can be replaced with a local type
    # r_seed = time(NULL) if io.r_seed < 0 else io.r_seed

    cdef void _seed_rng(self, t_option* io) noexcept nogil:
        cdef int r_seed = self.seed
        if r_seed < 0:
            r_seed = time(NULL)
        io.r_seed = r_seed
        srand(r_seed)

    cdef void _initialize_options(self, t_option* io) except *:
        # initialize RNG (--r_seed flag)
        self._seed_rng(io)

    # ---

    cdef t_tree* _bootstrap(
        self,
        t_calign* cdata,
        t_mod* mod,
        t_option* io,
        char* most_likely_tree,
    ):
        cdef t_tree* dum

        # FIXME: Move bootstrap to another post-process function?
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
            PhyML_Printf("\n\n. Printing the most likely tree in file '%s'.", phyml.utilities.Basename(io.out_tree_file))
        if io.n_data_sets == 1:
            if io.fp_out_tree is not NULL:
                rewind(io.fp_out_tree)

        # Recover most likely tree and add removed duplicate sequences
        # FIXME: avoid serialization / deserialization?
        dum      = phyml.io.Read_Tree(&most_likely_tree)
        dum.data = cdata
        dum.mod  = mod
        dum.io   = io
        phyml.utilities.Connect_CSeqs_To_Nodes(cdata, io, dum)
        phyml.utilities.Insert_Duplicates(dum)
        phyml.free.Free(most_likely_tree)

        return dum

    cpdef Result build(
        self,
        Alignment alignment,
        Tree start_tree = None,
        bint debug = False,
    ):
        """Build a tree from the given alignment.

        Arguments:
            alignment (`Alignment`): The multiple sequence alignment to
                estimate a phylogeny for.
            tree (`Tree` or `None`): A starting tree to use, or `None`
                to compute an initial tree with the BioNJ algorithm.

        """
        # TODO: allow multiple start tree

        cdef t_calign* cdata                    = NULL
        cdef t_option* io                       = NULL
        cdef t_tree*   tree                     = NULL
        cdef t_tree*   dum                      = NULL
        cdef t_mod*    mod                      = NULL
        cdef int       num_tree
        cdef int       num_rand_tree
        cdef time_t    time_beg
        cdef time_t    time_end
        cdef phydbl    best_lnL                 = phyml.utilities.UNLIKELY
        cdef char*     most_likely_tree         = NULL
        cdef bint      orig_random_input_tree

        cdef Tree                out_tree
        cdef CompressedAlignment out_compressed

        # Validate alignmentphyml.utilities.Connect_CSeqs_To_Nodes
        if alignment is None:
            raise TypeError("expected Alignment, found None")
        assert alignment._data is not NULL

        #ifdef QUIET
        # setvbuf(stdout, NULL, _IOFBF, 2048);
        #endif

        # Create a new `t_option` struct to handle PhyML configuration
        io = self._create_default_options()
        if io is NULL:
            raise ValueError("Failed parsing input")

        # Enable verbose output in debug mode
        io.quiet = not debug
        # Initialize options and seed global RNG
        self._initialize_options(io)

        # Use user-defined start tree
        if start_tree is not None:
            io.in_tree = 2  # in_tree=2 for user-defined tree, =1 for parsimony, =0 by default

        # if io.in_tree == 2:
        #     # Test_Multiple_Data_Set_Format(io)
        #     raise NotImplementedError("Test_Multiple_Data_Set_Format")
        # else:
        io.n_trees = 1

        # if io.n_trees == 0 and io.in_tree == 2:
        #     raise ValueError("invalid format of input tree")

        #
        # if io.n_data_sets > 1 and io.n_trees > 1:
        #     io.n_data_sets = min(io.n_trees, io.n_data_sets)
        #     io.n_trees     = min(io.n_trees, io.n_data_sets)

        # phyml.io.Get_Seq(io)
        io.n_otu = alignment._n_otu
        io.data = alignment._data
        io.datatype = alignment._datatype

        # phyml.make.Make_Model_Complete(io.mod)
        # phyml.utilities.Set_Model_Name(io.mod)
        if not io.quiet:
            phyml.io.Print_Settings(io)

        mod                    = io.mod
        orig_random_input_tree = io.mod.s_opt.random_input_tree

        cdata = phyml.utilities.Compact_Data(io.data, io)
        # phyml.free.Free_Seq(io.data, cdata.n_otu)

        for num_tree in range(io.n_trees):
            if not io.mod.s_opt.random_input_tree:
                io.mod.s_opt.n_rand_starts = 1

            if orig_random_input_tree and io.n_trees > 1:
                PhyML_Printf("\n== Cannot combine random starting trees with multiple input trees.")
                raise RuntimeError("EXIT")

            for num_rand_tree in range(io.mod.s_opt.n_rand_starts):

                if io.mod.s_opt.random_input_tree and io.mod.s_opt.topo_search != t_topo.NNI_MOVE:
                    if not io.quiet:
                        PhyML_Printf("\n\n. [Random start %3d/%3d]", num_rand_tree + 1, io.mod.s_opt.n_rand_starts)

                phyml.init.Init_Model(cdata, mod, io)
                phyml.models.Set_Model_Parameters(mod)

                # Make the initial tree
                # in_tree == 0 or in_tree == 1 means BioNJ
                # in_tree == 2 --> user-provided tree
                if io.in_tree == 0 or io.in_tree == 1:
                    tree = phyml.utilities.Dist_And_BioNJ(cdata, mod, io)
                    # if io.print_mat_and_exit:
                    #     phyml.io.Print_Mat(phyml.lk.ML_Dist(cdata, mod))
                    #     raise RuntimeError("exit", -1)
                elif io.in_tree == 2:
                    # NOTE: User-provided tree, make a copy to be safe
                    #       (but do not read from file like original code)
                    #tree = phyml.io.Read_User_Tree(cdata, mod, io)
                    tree = phyml.utilities.Duplicate_Tree(start_tree._tree)
                assert tree is not NULL # FIXME?

                if io.mod.s_opt.opt_topo:
                    phyml.utilities.Remove_Duplicates(cdata, io, tree)

                # TODO: Support constraint tree
                # if io.fp_in_constraint_tree is not NULL:
                #     PhyML_Printf("\n. Reading constraint tree file...")
                #
                #     io.cstr_tree = phyml.io.Read_Tree_File_Phylip(io.fp_in_constraint_tree)
                #     if io.cstr_tree.n_root is not NULL:
                #         PhyML_Printf("\n== The constraint tree file must be unrooted")
                #         raise RuntimeError("exit")
                #
                #     s = phyml.utilities.Add_Taxa_To_Constraint_Tree(io.fp_in_constraint_tree, cdata)
                #     fflush(NULL)
                #     phyml.free.Free_Tree(tree)
                #     tree = phyml.io.Read_Tree(&s)
                #     io.in_tree = 2
                #     phyml.free.Free(s)
                #     s = NULL
                #
                #     phyml.utilities.Check_Constraint_Tree_Taxa_Names(io.cstr_tree, cdata)
                #     phyml.utilities.Alloc_Bip(io.cstr_tree)
                #     phyml.utilities.Get_Bip(io.cstr_tree.a_nodes[0], io.cstr_tree.a_nodes[0].v[0], io.cstr_tree)
                #     if not tree.has_branch_lengths:
                #         phyml.utilities.Add_BioNJ_Branch_Lengths(tree, cdata, mod, NULL)

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
                if num_tree == 0 and num_rand_tree == 0:
                    phyml.utilities.Check_Memory_Amount(tree)

                # TODO: Support constraint tree
                # if io.cstr_tree is not NULL and not phyml.utilities.Check_Topo_Constraints(tree, io.cstr_tree):
                #     PhyML_Printf("\n\n== The initial tree does not satisfy the topological constraint.")
                #     PhyML_Printf("\n== Please use the user input tree option with an adequate tree topology.")
                #     raise RuntimeError("EXIT")

                phyml.utilities.Connect_CSeqs_To_Nodes(tree.data, tree.io, tree)
                phyml.make.Make_Tree_For_Pars(tree)
                phyml.make.Make_Tree_For_Lk(tree)
                phyml.make.Make_Spr(tree)
                phyml.utilities.Br_Len_Not_Involving_Invar(tree)
                phyml.utilities.Unscale_Br_Len_Multiplier_Tree(tree)

                # NOTE: Unused (I/O)
                # if tree.io.print_json_trace:
                #     phyml.io.JSON_Tree_Io(tree, tree.io.fp_out_json_trace)

                phyml.utilities.Set_Update_Eigen(True, tree.mod)
                phyml.lk.Lk(NULL, tree)
                phyml.utilities.Set_Update_Eigen(False, tree.mod)
                if not io.quiet:
                    PhyML_Printf("\n. Init log-likelihood: %f", tree.c_lnL)

                if tree.mod.s_opt.opt_topo:
                    phyml.spr.Global_Spr_Search(tree)
                    if tree.n_root:
                        phyml.utilities.Add_Root(tree.a_edges[0], tree)
                else:
                    if tree.mod.s_opt.opt_subst_param or tree.mod.s_opt.opt_bl_one_by_one:
                        phyml.optimiz.Round_Optimize(tree, phyml.utilities.ROUND_MAX)

                phyml.utilities.Set_Both_Sides(True, tree)
                phyml.lk.Lk(NULL, tree)
                phyml.pars.Pars(NULL, tree)
                phyml.utilities.Get_Tree_Size(tree)
                if not io.quiet:
                    PhyML_Printf("\n\n. Log likelihood of the current tree: %.*f.", DECIMAL_DIG, tree.c_lnL)

                if tree.io.ancestral:
                    phyml.ancestral.Ancestral_Sequences(tree, True)

                phyml.utilities.Check_Br_Lens(tree)
                phyml.utilities.Br_Len_Involving_Invar(tree)
                phyml.utilities.Rescale_Br_Len_Multiplier_Tree(tree)

                if tree.n_root is NULL:
                    phyml.utilities.Get_Best_Root_Position(tree)

                # Print the tree estimated using the current random (or BioNJ) starting tree
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

                    # FIXME: we don't have FP out but we could
                    #        capture the stats here rather than outputting
                    #        them to a file
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

                # Deallocate memory for the current tree
                phyml.free.Free_Best_Spr(tree)
                phyml.free.Free_Spr_List_One_Edge(tree)
                phyml.free.Free_Spr_List_All_Edge(tree)
                phyml.free.Free_Tree_Pars(tree)
                phyml.free.Free_Tree_Lk(tree)
                phyml.free.Free_Tree(tree)

            # Run bootstrap
            dum = self._bootstrap(cdata, mod, io, most_likely_tree)
            most_likely_tree = phyml.io.Write_Tree(dum)
            phyml.free.Free_Tree(dum)

        # Recover the most likely tree
        if most_likely_tree is not NULL:
            out_tree = Tree.__new__(Tree)
            out_tree._tree      = phyml.io.Read_Tree(&most_likely_tree)
            out_tree._tree.data = cdata
            out_tree._tree.mod  = mod
            out_tree._tree.io   = io
            phyml.utilities.Connect_CSeqs_To_Nodes(cdata, io, out_tree._tree)
            phyml.utilities.Insert_Duplicates(out_tree._tree)
            phyml.free.Free(most_likely_tree)
        else:
            out_tree = None

        # Keep the compressed alignment
        out_compressed = CompressedAlignment.__new__(CompressedAlignment)
        out_compressed._data = cdata

        # Free model
        phyml.free.Free_Model_Complete(mod)

        if mod.s_opt.n_rand_starts > 1:
            if not io.quiet:
                PhyML_Printf("\n. Best log likelihood: %f\n", best_lnL)

        # Free remaining data
        phyml.free.Free_Optimiz(mod.s_opt)
        phyml.free.Free_Model_Basic(mod)
        phyml.free.Free_Input(io)

        # Record and display finish time
        time(&time_end)
        if not io.quiet:
            phyml.io.Print_Time_Info(time_beg, time_end)

        # Return result structure
        return Result(
            tree=out_tree,
            compressed=out_compressed,
            alignment=alignment,
            log_likelihood=best_lnL,
        )
