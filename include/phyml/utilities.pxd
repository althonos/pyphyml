from libc.stdio cimport FILE
from libc.time cimport time_t

cdef extern from "utilities.h" nogil:

    cdef enum t_topo:
        NNI_MOVE = 0
        SPR_MOVE = 1
        BEST_OF_NNI_AND_SPR = 2

    # double-precision
    ctypedef double phydbl
    cdef phydbl LOG(phydbl)
    cdef phydbl POW(phydbl)
    cdef phydbl EXP(phydbl)
    cdef phydbl FABS(phydbl)
    cdef phydbl SQRT(phydbl)
    cdef phydbl CEIL(phydbl)
    cdef phydbl FLOOR(phydbl)
    cdef phydbl RINT(phydbl)
    cdef phydbl ROUND(phydbl)
    cdef phydbl TRUNC(phydbl)
    cdef phydbl COS(phydbl)
    cdef phydbl SIN(phydbl)
    cdef phydbl TAN(phydbl)
    
    const phydbl SMALL
    const phydbl BIG
    const phydbl LOGBIG
    const phydbl LOGSMALL
    const phydbl UNLIKELY

    const int    N_MAX_INSERT
    const int    N_MAX_OTU
    const int    ROUND_MAX
    const int    T_MAX_ALPHABET

    cdef struct __Generic_LL:
        void*         v
        __Generic_LL* next
        __Generic_LL* prev
        __Generic_LL* tail
        __Generic_LL* head
    ctypedef __Generic_LL t_ll

    cdef struct __Scalar_Int:
        int           v
        bint          optimize
        __Scalar_Int* next
        __Scalar_Int* prev
    ctypedef __Scalar_Int scalar_int

    cdef struct __Scalar_Dbl:
        phydbl        v
        bint          onoff
        bint          optimize
        bint          print
        __Scalar_Dbl* next
        __Scalar_Dbl* prev
    ctypedef __Scalar_Dbl scalar_dbl

    cdef struct __Vect_Int:
        int*        v
        int         len
        bint        optimize
        __Vect_Int* next
        __Vect_Int* prev
    ctypedef __Vect_Int vect_int

    cdef struct __Vect_Dbl:
        phydbl*     v
        int         len
        bint        optimize
        __Vect_Dbl* next
        __Vect_Dbl* prev
    ctypedef __Vect_Dbl vect_dbl

    cdef struct __String:
        char*     s
        int       len
        __String* next
        __String* prev
    ctypedef __String t_string

    cdef struct __Node:
        __Node**  v         # table of pointers to neighbor nodes. Dimension = 3 
        __Node*** bip_node; # three lists of pointer to tip nodes. One list for each direction
        # struct __Edge **b;         /*! table of pointers to neighbor branches */
        # struct __Node  *anc;   /*! direct ancestor t_node (for rooted tree only) */
        # struct __Edge  *b_anc; /*! edge between this node and its direct ancestor (for
        #                             rooted tree only) */
        # struct __Node  *ext_node;
        # struct __Node  *match_node;
        # struct __Align *c_seq;     /*! corresponding compressed sequence */
        # struct __Align *c_seq_anc; /*! corresponding compressed ancestral sequence */

        # struct __Node *next; /*! tree->a_nodes[i]->next <=> tree->next->a_nodes[i] */
        # struct __Node *prev; /*! See above */
        # struct __Node *next_mixt; /*! Next mixture tree*/
        # struct __Node *prev_mixt; /*! Parent mixture tree */
        # struct __Calibration *
        #     *cal; /*! List of calibration constraints attached to this node */
        # struct __Lindisk_Node
        #     *ldsk; /*! Used in PhyREX. Lineage/Disk this node corresponds to */
        # struct __Node *rk_next; /*! Next node in the list of ranked nodes (from oldest
        #                             to youngest) */
        # struct __Node *rk_prev; /*! Previous node in the list of ranked nodes (from
        #                             oldest to youngest) */
        # struct __Label *label;

        # int *bip_size; /*! Size of each of the three lists from bip_node */
        # int  num;      /*! t_node number */
        # int  tax;      /*! tax = 1 -> external node, else -> internal t_node */

        # char *name;     /*! taxon name (if exists) */
        # char *ori_name; /*! taxon name (if exists) */
        # int   n_cal;    /*! Number of calibration constraints */

        # phydbl *score; /*! score used in BioNJ to determine the best pair of nodes to
        #                     agglomerate */
        # phydbl dist_to_root; /*! distance to the root t_node */

        # short int common;
        # phydbl    y_rank;
        # phydbl    y_rank_ori;
        # phydbl    y_rank_min;
        # phydbl    y_rank_max;

        # int *s_ingrp;  /*! does the subtree beneath belong to the ingroup */
        # int *s_outgrp; /*! does the subtree beneath belong to the outgroup */

        # int id_rank; /*! order taxa alphabetically and use id_rank to store the ranks
        #                 */
        int rank
        int rank_max

        # struct __Geo_Coord *coord;
        # struct __Geo_Veloc *veloc;
    ctypedef __Node t_node

    cdef struct __Edge:
        __Node *left, # t_node on the left side of the t_edge
        __Node *rght; # t_node on the right side of the t_edge
        # short int      l_r, r_l, l_v1, l_v2, r_v1, r_v2;
        # /*! these are directions (i.e., 0, 1 or 2): */
        # /*! l_r (left to right) -> left[b_fcus->l_r] = right */
        # /*! r_l (right to left) -> right[b_fcus->r_l] = left */
        # /*! l_v1 (left t_node to first t_node != from right) -> left[b_fcus->l_v1] =
        # * left_1 */
        # /*! l_v2 (left t_node to secnd t_node != from right) -> left[b_fcus->l_v2] =
        # * left_2 */
        # /*! r_v1 (right t_node to first t_node != from left) -> right[b_fcus->r_v1] =
        # * right_1 */
        # /*! r_v2 (right t_node to secnd t_node != from left) -> right[b_fcus->r_v2] =
        # * right_2 */

        # struct __NNI   *nni;
        # struct __Edge  *next;
        # struct __Edge  *prev;
        # struct __Edge  *next_mixt;
        # struct __Edge  *prev_mixt;
        # struct __Label *label;

        # int         num;       /*! branch number */
        # scalar_dbl *l;         /*! branch length */
        # scalar_dbl *best_l;    /*! best branch length found so far */
        # scalar_dbl *l_old;     /*! old branch length */
        # scalar_dbl *l_var;     /*! variance of edge length */
        # scalar_dbl *l_var_old; /*! variance of edge length (previous value) */

        # int bip_score; /*! score of the bipartition generated by the corresponding
        # edge bip_score = 1 iif the branch is found in both trees to be compared,
        # bip_score = 0 otherwise. */
        # phydbl tdist_score; /*! average transfer distance over bootstrap trees */

        # int num_st_left; /*! number of the subtree on the left side */
        # int num_st_rght; /*! number of the subtree on the right side */

        # phydbl *Pij_rr;  /*! matrix of change probabilities and its first and secnd
        #                     derivates (rate*state*state) */
        # phydbl *tPij_rr; /*! transpose matrix of change probabilities and its first
        #                     and secnd derivates (rate*state*state) */
        # #ifdef BEAGLE
        # int Pij_rr_idx;
        # #endif

        # short int *div_post_pred_left; /*! posterior prediction of nucleotide/aa
        #                                     diversity (left-hand subtree) */
        # short int *div_post_pred_rght; /*! posterior prediction of nucleotide/aa
        #                                     diversity (rght-hand subtree) */
        # short int does_exist;

        # phydbl *p_lk_left,
        #     *p_lk_rght; /*! likelihoods of the subtree on the left and right side (for
        #                     each site and each relative rate category) */
        # phydbl *p_lk_tip_r, *p_lk_tip_l;

        # #ifdef BEAGLE
        # int p_lk_left_idx, p_lk_rght_idx;
        # int p_lk_tip_idx;
        # #endif

        # int *patt_id_left;
        # int *patt_id_rght;
        # int *p_lk_loc_left;
        # int *p_lk_loc_rght;

        # int *pars_l, *pars_r; /*! parsimony of the subtree on the left and right sides
        #                         (for each site) */
        # int *ui_l, *ui_r; /*! union - intersection vectors used in Fitch's parsimony
        #                     algorithm */
        # int *p_pars_l, *p_pars_r; /*! conditional parsimony vectors */

        # /*! Below are the likelihood scaling factors (used in functions
        #     `Get_All_Partial_Lk_Scale' in lk.c. */
        # /*
        #     For every site, every subtree and every rate class, PhyML maintains
        #     a`sum_scale_pow' value where sum_scale_pow = sum_scale_pow_v1 +
        #     sum_scale_pow_v2 + curr_scale_pow' sum_scale_pow_v1 and sum_scale_pow_v2 are
        #     sum_scale_pow of the left and right subtrees. curr_scale_pow is an integer
        #     greater than one. The smaller the partials, the larger curr_scale_pow.

        #     Now the partials for this subtree are scaled by *multiplying* each of
        #     them by 2^curr_scale_pow. The reason for doing the scaling this way is
        #     that multiplications by 2^x (x an integer) can be done in an 'exact'
        #     manner (i.e., there is no loss of numerical precision)

        #     At the root edge, the log-likelihood is then
        #     logL = logL' - (sum_scale_pow_left + sum_scale_pow_right)log(2),
        #     where L' is the scaled likelihood.
        # */

        # int *sum_scale_left_cat;
        # int *sum_scale_rght_cat;
        # int *sum_scale_left;
        # int *sum_scale_rght;

        # phydbl bootval; /*! bootstrap value (if exists) */

        # short int is_alive; /*! is_alive = 1 if this t_edge is used in a tree */

        # phydbl dist_btw_edges;
        # int    topo_dist_btw_edges;

        # int has_zero_br_len;

        phydbl ratio_test;     # approximate likelihood ratio test 
        phydbl alrt_statistic; # aLRT statistic 
        phydbl support_val

        int n_jumps; # number of jumps of substitution rates 

        # int *n_diff_states_l; # Number of different states found in the subtree on the left of this edge 
        # int *n_diff_states_r; # Number of different states found in the subtree on the right of this edge

        # phydbl bin_cod_num;

        # short int update_partial_lk_left;
        # short int update_partial_lk_rght;
    ctypedef __Edge t_edge

    cdef struct __Tree:
        __Node*   n_root   # root t_node
        __Edge*   e_root   # t_edge on which lies the root
        __Node**  a_nodes; # array of nodes that defines the tree topology 
        __Edge**  a_edges; # array of edges 
        __Model*  mod      # substitution model
        __Calign* data     # sequences 
        __Tree*   next     # set to NULL by default. Used for mixture models 
        __Tree*   prev     # set to NULL by default. Used for mixture models 
        # __Tree   *next_mixt; /*! set to NULL by default. Used for mixture models */
        # __Tree   *prev_mixt; /*! set to NULL by default. Used for mixture models */
        # __Tree   *mixt_tree; /*! set to NULL by default. Used for mixture models */
        # __Tree   **aux_tree;   /*! set to NULL by default. Used as a latent variable in molecular dating */
        __Option* io       # input/output 
        __Matrix* mat      # pairwise distance matrix 
        # __Node  **curr_path; /*! list of nodes that form a path in the tree */
        # __SPR   **spr_list_one_edge;
        # __SPR   **spr_list_all_edge;
        # __SPR    *best_spr;

        # __Tdraw  *ps_tree; /*! structure for drawing trees in postscript format */
        # __T_Rate       *rates; /*! structure for handling rates of evolution */
        # __T_Time       *times; /*! structure for handling node ages */
        # __Tmcmc        *mcmc;
        # __Phylogeo     *geo;
        # __Migrep_Model *mmod;
        # __Disk_Event   *young_disk;  /*! Youngest disk (i.e., disk which age is
        #                                         the closest to present). Used in PhyREX */
        # __Disk_Event *old_samp_disk; /*! Oldest sampled disk. Used in PhyREX */
        # __XML_node   *xml_root;
        # __Generic_LL *edge_list;
        # struct __Generic_LL *node_list;
        # struct __Independent_Contrasts
        #     *ctrst; /*! Pointer to data structure used for independent contrasts */
        # struct __Continuous_Model *contmod;

        # #if (defined(__AVX__) || defined(__AVX2__))
        # __m256d *_tPij1, *_tPij2, *_pmat1plk1, *_pmat2plk2, *_plk0, *_l_ev, *_r_ev,
        #     *_prod_left, *_prod_rght;
        # #elif (defined(__SSE__) || defined(__SSE2__) || defined(__SSE3__) ||           \
        #     defined(__ARM_NEON))
        # __m128d *_tPij1, *_tPij2, *_pmat1plk1, *_pmat2plk2, *_plk0, *_l_ev, *_r_ev,
        #     *_prod_left, *_prod_rght;
        # #endif

        # phydbl *p_lk_left_pi, *l_ev;
        # phydbl *big_lk_array;
        # int     big_lk_array_pos;

        # short int eval_alnL; /*! Evaluate likelihood for genetic data */
        # short int eval_rlnL; /*! Evaluate likelihood for rates along the tree */
        # short int eval_glnL; /*! Evaluate likelihood under phylogeo model */
        # short int eval_tlnL; /*! Evaluate likelihood under tree-generating process */
        # short int scaling_method;
        # short int fully_nni_opt;
        # short int numerical_warning;

        # short int use_eigen_lr;
        # int       is_mixt_tree;
        # int       tree_num;  /*! tree number. Used for mixture models */
        # int ps_page_number;  /*! when multiple trees are printed, this variable give
        #                         the current page number */
        # int depth_curr_path; /*! depth of the t_node path defined by curr_path */
        bint has_bip;          # has_bip=1, then the structure to compare tree topologies is allocated, has_bip=0 otherwise 
        int n_otu;    # number of taxa 
        int curr_site # current site of the alignment to be processed
        int curr_catg # current class of the discrete gamma rate distribution
        int n_swap    # number of NNIs performed 
        bint has_branch_lengths; # =1 iff input tree displays branch lengths 
        # short int both_sides;   /*! both_sides=1 -> a pre-order and a post-order tree
        #     traversals are required to compute the likelihood
        #     of every subtree in the phylogeny*/
        # int num_curr_branch_available; /*!gives the number of the next cell in a_edges
        #                                     that is free to receive a pointer to a
        #                                     branch */
        # short int *t_dir;
        # int        n_moves;
        # int        verbose;

        # int dp;        /*! Data partition */
        # int s_mod_num; /*! Substitution model number */
        # int lock_topo; /*! = 1 any subsequent topological modification will be
        #                     banished */
        # int       print_labels;
        # int       write_br_lens;
        # int      *mutmap; /*! Mutational map */
        # int       json_num;
        # short int update_eigen_lr;
        # int       tip_root; /*! Index of tip node used as the root */
        # phydbl   *dot_prod;

        # phydbl *expl;

        phydbl  init_lnL
        phydbl  best_lnL  # highest value of the loglikelihood found so far
        int     best_pars # highest value of the parsimony found so far
        phydbl  c_lnL;     # loglikelihood *
        phydbl  p_lnL;     # loglikelihood (previous value)
        # phydbl  old_lnL;   /*! old loglikelihood */
        # phydbl  sum_min_sum_scale; /*! common factor of scaling factors */
        # phydbl *c_lnL_sorted;      /*! used to compute c_lnL by adding sorted terms to
        #                                 minimize CPU errors */
        # phydbl *cur_site_lk;       /*! vector of loglikelihoods at individual sites */
        # phydbl *old_site_lk;       /*! vector of likelihoods at individual sites */
        # phydbl  annealing_temp;    /*! annealing temperature in simulated annealing
        #                                 optimization algo */
        # phydbl c_dlnL;  /*! First derivative of the log-likelihood with respect to the
        #                     length of a branch */
        # phydbl c_d2lnL; /*! Second derivative of the log-likelihood with respect to
        #                     the length of a branch */

        # phydbl *unscaled_site_lk_cat; /*! partially scaled site likelihood at
        #                                 individual sites */

        # phydbl *site_lk_cat; /*! loglikelihood at a single site and for each class of
        #                         rate*/
        # phydbl unconstraint_lk; /*! unconstrained (or multinomial) log-likelihood  */
        # phydbl composite_lk;    /*! composite log-likelihood  */
        # int   *fact_sum_scale;
        # phydbl **log_lks_aLRT; /*! used to compute several branch supports */
        # phydbl   n_root_pos;   /*! position of the root on its t_edge */
        # phydbl   size;         /*! tree size */
        # int     *site_pars;
        # int      c_pars;
        # int     *step_mat;

        # int size_spr_list_one_edge;
        # int size_spr_list_all_edge;
        # int perform_spr_right_away;

        time_t t_beg
        time_t t_current

        # int bl_from_node_stamps; /*! == 1 -> Branch lengths are determined by t_node
        #                             times */
        # phydbl sum_y_dist_sq;
        # phydbl sum_y_dist;
        # phydbl tip_order_score;
        # int    write_tax_names;
        # int    update_alias_subpatt;

        # phydbl geo_mig_sd; /*! standard deviation of the migration step random
        #                         variable */
        # phydbl geo_lnL;    /*! log likelihood of the phylo-geography model */

        # int bl_ndigits;

        # phydbl *short_l;   /*! Vector of short branch length values */
        # int     n_short_l; /*! Length of short_l */
        # phydbl  norm_scale;

        # short int br_len_recorded;

        # short int apply_lk_scaling; /*! Applying scaling of likelihoods. YES/NO */

        # phydbl *K; /*! a vector of the norm.constants for the node times prior. */

        # short int ignore_root;
        # short int ignore_mixt_info;
        # #ifdef BEAGLE
        # int b_inst; /*! The BEAGLE instance id associated with this tree. */
        # #endif

        # // Extra partial lk structure for bookkeeping
        # short int *div_post_pred_extra_0;
        # int       *sum_scale_cat_extra_0;
        # int       *sum_scale_extra_0;
        # phydbl    *p_lk_extra_0;
        # phydbl    *p_lk_tip_extra_0;
        # int       *patt_id_extra_0;

        # short int *div_post_pred_extra_1;
        # int       *sum_scale_cat_extra_1;
        # int       *sum_scale_extra_1;
        # phydbl    *p_lk_extra_1;
        # phydbl    *p_lk_tip_extra_1;
        # int       *patt_id_extra_1;

        # int n_edges_traversed;
        int n_tot_bl_opt

        bint opt_topo
    ctypedef __Tree t_tree

    cdef struct __Align:
        char*      name      # sequence name
        int        len       # sequence length
        char*      state     # sequence itself 
        short int* d_state   # sequence itself (digits)
        short int* is_ambigu # is_ambigu[site] = 1 if state[site] is an ambiguous character. 0 otherwise
        short int  is_duplicate
        int        num
    ctypedef __Align align

    cdef struct __Calign:
        __Align** c_seq     # compressed sequences
        __Align** c_seq_rm  # removed sequences
        __Option* io        # input/output
        phydbl*       obs_state_frq; # observed state frequencies
        short int*    invar;         # < 0 -> polymorphism observed
        phydbl*       wght;          # # of each site in c_align
        short int*    ambigu; # ambigu[i]=1 is one or more of the sequences at site i display an ambiguous character 
        
        phydbl obs_pinvar
        int    n_otu;        # number of taxa
        int    n_rm;         # number of taxa removed
        int    clean_len;    # uncrunched sequences lenghts without gaps
        int    n_pattern;    # crunched sequences lengths
        int    init_len;     # length of the uncompressed sequences
        int*   sitepatt;     # this array maps the position of the patterns in the compressed alignment to the positions in the uncompressed one
        int         format;  # 0 (default): PHYLIP. 1: NEXUS. 
        scalar_dbl *io_wght; # weight of each *site* (not pattern) given as input 

        int       n_masked;   # Number of masked positions (or columns, depending on mask_type)
        short int mask_type;  # MASK_TYPE_POSITION or MASK_TYPE_COLUMN
        int*      masked_pos; # Vector of masked positions/columns
    ctypedef __Calign calign

    cdef struct __Matrix:
                                 # mostly used in BIONJ
        phydbl **P, **Q, **dist; # observed proportions of transition, transverion and  distances between pairs of sequences
        t_tree *tree # tree
        int *on_off  # on_off[i]=1 if column/line i corresponds to a t_node that has not been agglomerated yet
        int             n_otu # number of taxa
        char**          name  # sequence names
        int             r     # number of nodes that have not been agglomerated yet 
        __Node**        tip_node  # array of pointer to the leaves of the tree 
        int             curr_int  # used in the NJ/BIONJ algorithms
        int             method    # if method=1->NJ method is used, BIONJ otherwise
    ctypedef __Matrix matrix

    cdef struct __RateMatrix:
        int n_diff_rr; # number of different relative substitution rates in the custom model
        vect_dbl *rr;  # relative rate parameters of the GTR or custom model (rescaled)
        vect_dbl *rr_val; # log of relative rate parameters of the GTR or custom model (unscaled) 
        vect_int *rr_num
        vect_int *n_rr_per_cat # number of rate parameters in each category
        vect_dbl *qmat
        vect_dbl *qmat_buff
        bint optimize
        __RateMatrix *next
        __RateMatrix *prev
    ctypedef __RateMatrix t_rmat

    cdef struct __Option:
        __Model*  mod  # pointer to a substitution model
        __Tree*   tree;      # pointer to the current tree 
        __Align** data # pointer to the uncompressed sequences
        __Tree*   cstr_tree; # pointer to a constraint tree (can be a multifurcating one) 
        __Calign* cdata; # pointer to the compressed sequences
        # struct __Super_Tree *st;    /*! pointer to supertree */
        # struct __Tnexcom   **nex_com_list;
        # struct __List_Tree  *treelist; /*! list of trees. */
        # struct __Option     *next;
        # struct __Option     *prev;
        # struct __Tmcmc      *mcmc;
        # struct __T_Rate     *rates;
        # struct __T_Time     *times;

        int interleaved # interleaved or sequential sequence file format ?
        int in_tree     # =1 iff a user input tree is used as input

        char *in_align_file # alignment file name
        FILE *fp_in_align   # pointer to the alignment file

        char *in_tree_file; # input tree file name 
        FILE *fp_in_tree;   # pointer to the input tree file 

        char *in_constraint_tree_file; # input constraint tree file name 
        FILE *fp_in_constraint_tree;   # pointer to the input constraint tree file 

        char *out_tree_file # name of the tree file 
        FILE *fp_out_tree

        char *weight_file;   # name of the file containing site weights
        FILE *fp_weight_file

        char *out_trees_file; # name of the tree file 
        FILE *fp_out_trees;   # pointer to the tree file containing all the trees estimated using random starting trees 

        char *out_boot_tree_file; # name of the tree file 
        FILE *fp_out_boot_tree;   # pointer to the bootstrap tree file

        char *out_boot_stats_file; # name of the tree file 
        FILE *fp_out_boot_stats;   # pointer to the statistics file 

        char *out_stats_file # name of the statistics file
        FILE *fp_out_stats

        char *out_trace_file # name of the file in which the trace is written
        FILE *fp_out_trace

        char *out_json_trace_file; #! name of the file in which json trace is written 
        FILE *fp_out_json_trace

        char *out_lk_file; # name of the file in which the likelihood of the model is written 
        FILE *fp_out_lk

        char *out_summary_file; # name of the file in which summary statistics are written 
        FILE *fp_out_summary

        char *out_ps_file # name of the file in which tree(s) is(are) written 
        FILE *fp_out_ps

        char *out_ancestral_seq_file   # name of the file containing the ancestral sequences 
        char *out_ancestral_tree_file; # name of the file containing the tree with internal node labelled according to refs in ancestral_seq_file 

        FILE *fp_out_ancestral_seq   # pointer to the file containing the ancestral sequences 
        FILE *fp_out_ancestral_tree  # pointer to the file containing the tree with labels on internal nodes  

        char *in_xml_file
        FILE *fp_in_xml   # pointer to the file containing XML formatted data */

        char *in_coord_file # name of input file containing coordinates 
        FILE *fp_in_coord   # pointer to the file containing coordinates 

        char *out_file; # name of the output file 

        # char *clade_list_file;

        # int datatype;         /*! 0->DNA, 1->AA */
        # int print_boot_trees; /*! =1 if the bootstrapped trees are printed in output
        #                         */
        # int   out_stats_file_open_mode; /*! opening file mode for statistics file */
        # int   out_tree_file_open_mode;  /*! opening file mode for tree file */
        int   n_data_sets              # number of data sets to be analysed
        int   n_trees                  # number of trees
        int   init_len;                # sequence length 
        int   n_otu;                   # number of taxa 
        int   n_data_set_asked;        # number of bootstrap replicates 
        char* nt_or_cd;                # nucleotide or codon data ? (not used)
        bint  multigene;               # if=1 -> analyse several partitions. 
        int   config_multigene
        # int   n_part; /*! number of data partitions */
        # int   curr_gt;
        int   ratio_test # from 1 to 4 for specific branch supports, 0 of not 
        # int   ready_to_go;
        # int   data_file_format; /*! Data format: Phylip or Nexus */
        # int   tree_file_format; /*! Tree format: Phylip or Nexus */
        # int   state_len;

        # int curr_interface;
        int r_seed  # random seed
        bint collapse_boot; # 0 -> branch length on bootstrap trees are not collapsed if too small
        bint random_boot_seq_order; # !0 -> sequence order in bootstrapped data set is random
        bint print_trace
        bint print_json_trace
        bint print_site_lnl
        # int m4_model;
        bint rm_ambigu; # 0 is the default. 1: columns with ambiguous characters are discarded prior further analysis
        # int   colalias;
        # int   append_run_ID;
        # char *run_id_string;
        bint   quiet; # 0 is the default. 1: no interactive question (for batch mode)

        # int    lk_approx; /* EXACT or NORMAL */
        char **  alphabet
        int      codpos
        int      mutmap
        bint     use_xml

        char **long_tax_names
        char **short_tax_names
        int    size_tax_names

        phydbl *z_scores
        phydbl *lat
        phydbl *lon

        int boot_prog_every
        int mem_question
        int do_alias_subpatt

        # #ifdef BEAGLE
        # int beagle_resource;
        # #endif

        bint ancestral
        bint has_io_weights
        int tbe_bootstrap # Replace standard bootstrap with tbe bootstrap (only when b>0) 

        int leave_duplicates   # Leave duplicated sequences 
        int precision          # Decimal output precision for values in stats file 
        int n_boot_replicates

        bint print_mat_and_exit 
        bint print_node_num     # print node numbers if print_node_num=1
        bint print_support_val

        bint do_tbe
        bint do_boot
        bint do_alrt

        # short int edge_len_unit;

        # short int mcmc_output_times;
        # short int mcmc_output_trees;
    ctypedef __Option option

    cdef struct __Optimiz:
        bint opt_subst_param; # if opt_topo=0 and opt_subst_param=1 -> the
                              # numerical parameters of the model are optimised. if opt_topo=0 and
                              # opt_free_param=0 -> no parameter is optimised */
        # short int opt_clock_r;
        bint opt_bl_one_by_one; # =1 -> the branch lengths are optimised
        bint opt_topo;          # =1 -> the tree topology is optimised
        t_topo topo_search
        # short int opt_node_ages
        # short int opt_neff

        # phydbl init_lk; /*! initial loglikelihood value */
        # int n_it_max; /*! maximum bnumber of iteration during an optimisation step */
        # int last_opt; /*! =1 -> the numerical parameters are optimised further while
        #                 the tree topology remains fixed */
        bint random_input_tree; # boolean
        int n_rand_starts       # number of random starting points
        int brent_it_max
        int steph_spr
        int opt_five_branch
        int pars_thresh
        int hybrid_thresh
        int opt_br_len_mult
        int min_n_triple_moves
        int max_rank_triple_move
        int n_improvements
        int max_spr_depth
        int max_no_better_tree_found

        # phydbl tree_size_mult; /*! tree size multiplier */
        # phydbl min_diff_lk_local;
        # phydbl min_diff_lk_global;
        # phydbl min_diff_lk_move;
        # phydbl p_moves_to_examine;
        # int    fast_nni;
        # int    greedy;
        # int    general_pars;
        # int    quickdirty;
        # int    spr_pars;
        # int    spr_lnL;
        # int    max_depth_path;
        # int    min_depth_path;
        # int    deepest_path;
        # int    eval_list_regraft;
        # phydbl max_delta_lnL_spr;
        # phydbl max_delta_lnL_spr_current;
        # phydbl worst_lnL_spr;
        # int    br_len_in_spr;
        # int    opt_free_mixt_rates;
        # int    constrained_br_len;
        # int    opt_gamma_br_len;
        # int    first_opt_free_mixt_rates;
        # int    wim_n_rgrft;
        # int    wim_n_globl;
        # int    wim_max_dist;
        # int    wim_n_optim;
        # int    wim_n_best;
        # int    wim_inside_opt;

        # int opt_rmat_weight;
        # int opt_efrq_weight;

        # int skip_tree_traversal;
        # int serial_free_rates;

        # int curr_opt_free_rates;

        # int nni_br_len_opt;

        # int apply_spr_right_away;
        # int apply_spr;

        # phydbl l_min_spr;
    ctypedef __Optimiz t_opt

    cdef struct __Model:
        __Optimiz* s_opt # pointer to parameters to optimize
        # struct __Eigen      *eigen;
        # struct __M4         *m4mod;
        __Option*  io
        __Model*   next
        __Model*   prev
        __Model*   next_mixt
        __Model*   prev_mixt
        # struct __RateMatrix *r_mat;
        # struct __EquFreq    *e_frq;
        # struct __RAS        *ras;

        # t_string *aa_rate_mat_file;
        # FILE     *fp_aa_rate_mat;

        # t_string *modelname;
        # t_string *custom_mod_string; /*! string of characters used to define custom
        #                                 models of substitution */

        # int mod_num; /*! model number */

        # int update_eigen; /*! update_eigen=1-> eigen values/vectors need to be updated
        #                     */
        # int cv_type;      /* Type of cross-validation method */

        # int whichmodel;
        # int is_mixt_mod;
        # int augmented;
        # int ns; /*! number of states (4 for ADN, 20 for AA) */

        # int use_m4mod; /*! Use a Markov modulated Markov model ? */

        # scalar_dbl *kappa;  /*! transition/transversion rate */
        # scalar_dbl *lambda; /*! parameter used to define the ts/tv ratios in the F84
        #                         and TN93 models */
        # scalar_dbl *br_len_mult; /*! when users want to fix the relative length of
        #                             edges and simply estimate the total length of the
        #                             tree. This multiplier does the trick */
        # scalar_dbl *br_len_mult_unscaled;

        # vect_dbl   *Pij_rr; /*! matrix of change probabilities */
        # scalar_dbl *mr;     /*! mean rate = branch length/time interval  mr =
        #                         -sum(i)(vct_pi[i].mat_Q[ii]) */
        # scalar_dbl *aic;
        # scalar_dbl *bic;

        # short int
        #     log_l; /*! Edge lengths are actually log(Edge lengths) if log_l == YES !*/
        # phydbl l_min; /*! Minimum branch length !*/
        # phydbl l_max; /*! Maximum branch length !*/

        # scalar_dbl *l_var_sigma; /*! For any edge b we have b->l_var->v = l_var_sigma
        #                             * (b->l->v)^2 */
        # phydbl l_var_min; /*! Min of variance of branch lengths (used in conjunction
        #                     with gamma_mgf_bl == YES) */
        # phydbl l_var_max; /*! Max of variance of branch lengths (used in conjunction
        #                     with gamma_mgf_bl == YES) */

        # int gamma_mgf_bl; /*! P = \int_0^inf exp(QL) p(L) where L=\int_0^t R(s) ds and
        #                     p(L) is the gamma density. Set to NO by default !*/

        # int n_mixt_classes; /* Number of classes in the mixture model. */

        # scalar_dbl *r_mat_weight;
        # scalar_dbl *e_frq_weight;
        # #ifdef BEAGLE
        # int  b_inst;
        # bint optimizing_topology; /*! This is a flag that prevents the resetting of
        #                             category weights. Why? Read */
        # /*  Recall that while optimizing the topology, PhyML temporarily only uses 2
        # *  rate categories. Recall also that a BEAGLE instance is created with all
        # * the required categories, but we temporarily assign 0 weight to the other
        # * categories thus effectively using only 2 categories. However, subsequent
        # * calls to update the rates (i.e. update_beagle_ras()) will reset the
        # * weights. This flag prevents this resetting from happening */
        # #endif

    ctypedef __Model t_mod

    # void    Unroot_Tree(char **subtrees);
    # void    Set_Edge_Dirs(t_edge *b, t_node *a, t_node *d, t_tree *tree);
    # void    Restrict_To_Coding_Position(align **data, option *io);
    void      Uppercase(char *ch)
    void      Lowercase(char *ch)
    calign*   Compact_Data(align **data, option *io)
    calign*   Compact_Cdata(calign *data, option *io)
    # void    Traverse_Prefix_Tree(int site, int seqnum, int *patt_num, int *n_patt,
    #                             align **data, option *io, pnode *n);
    # pnode  *Create_Pnode(int size);
    # void    Get_Base_Freqs(calign *data);
    # void    Get_AA_Freqs(calign *data);
    # void    Swap_Nodes_On_Edges(t_edge *e1, t_edge *e2, int swap, t_tree *tree);
    # void    Connect_Edges_To_Nodes_Recur(t_node *a, t_node *d, t_tree *tree);
    # void    Connect_One_Edge_To_Two_Nodes(t_node *a, t_node *d, t_edge *b,
    #                                     t_tree *tree);
    # void    Update_Dirs(t_tree *tree);
    # void    Exit(char *message);
    # void   *mCalloc(int nb, size_t size);
    # void   *mRealloc(void *p, int nb, size_t size);
    # int     Sort_Phydbl_Decrease(const void *a, const void *b);
    # void    Qksort_Int(int *A, int *B, int ilo, int ihi);
    # void    Qksort(phydbl *A, phydbl *B, int ilo, int ihi);
    # void    Qksort_Matrix(phydbl **A, int col, int ilo, int ihi);
    # void    Order_Tree_Seq(t_tree *tree, align **data);
    char*     Add_Taxa_To_Constraint_Tree(FILE* fp, calign* cdata)
    void      Check_Constraint_Tree_Taxa_Names(t_tree* tree, calign* cdata)
    # void    Order_Tree_CSeq(t_tree *tree, calign *cdata);
    # void    Init_Mat(matrix *mat, calign *data);
    # void    Copy_Tax_Names_To_Tip_Labels(t_tree *tree, calign *data);
    # void    Share_Lk_Struct(t_tree *t_full, t_tree *t_empt);
    # void    Share_Spr_Struct(t_tree *t_full, t_tree *t_empt);
    # void    Share_Pars_Struct(t_tree *t_full, t_tree *t_empt);
    # int     Sort_Edges_NNI_Score(t_tree *tree, t_edge **sorted_edges, int n_elem);
    # int     Sort_Edges_Depth(t_tree *tree, t_edge **sorted_edges, int n_elem);
    # void    NNI(t_tree *tree, t_edge *b_fcus, int do_swap);
    # void    NNI_Pars(t_tree *tree, t_edge *b_fcus, int do_swap);
    # void    Swap(t_node *a, t_node *b, t_node *c, t_node *d, t_tree *tree);
    # void    Update_SubTree_Partial_Lk(t_edge *b_fcus, t_node *a, t_node *d,
    #                                 t_tree *tree);
    # void    Copy_Seq_Names_To_Tip_Labels(t_tree *tree, calign *data);
    # calign *Copy_Cseq(calign *ori, option *io, t_tree *tree);
    # int     Filexists(char *filename);
    # int     Is_Invar(int patt_num, int stepsize, int datatype, calign *data);
    # int     Is_Ambigu(char *state, int datatype, int stepsize);
    # void    Check_Ambiguities(calign *data, int datatype, int stepsize);
    # int     Get_State_From_Ui(int ui, int datatype);
    # int     Assign_State(char *c, int datatype, int stepsize);
    # char    Reciproc_Assign_State(int i_state, int datatype);
    # int     Assign_State_With_Ambiguity(char *c, int datatype, int stepsize);
    # void    Clean_Tree_Connections(t_tree *tree);
    # void    Bootstrap(t_tree *tree);
    void      Br_Len_Involving_Invar(t_tree* tree)
    void      Br_Len_Not_Involving_Invar(t_tree* tree)
    # void    Getstring_Stdin(char *s);
    # phydbl  Num_Derivatives_One_Param(phydbl (*func)(t_tree *tree), t_tree *tree,
    #                                 phydbl f0, phydbl *param, int which,
    #                                 int n_param, phydbl stepsize, short int logt,
    #                                 short int expt, phydbl *err, int precise,
    #                                 int is_positive);
    # phydbl  Num_Derivatives_One_Param_Nonaligned(
    #     phydbl (*func)(t_tree *tree), t_tree *tree, phydbl f0, phydbl **param,
    #     int which, int n_param, phydbl stepsize, short int logt, short int expt,
    #     phydbl *err, int precise, int is_positive);
    # int     Num_Derivative_Several_Param(t_tree *tree, phydbl *param, int n_param,
    #                                     phydbl stepsize, short int logt,
    #                                     short int expt, phydbl (*func)(t_tree *tree),
    #                                     phydbl *derivatives, int is_positive);
    # int     Num_Derivative_Several_Param_Nonaligned(t_tree *tree, phydbl **param,
    #                                                 int n_param, phydbl stepsize,
    #                                                 short int logt, short int expt,
    #                                                 phydbl (*func)(t_tree *tree),
    #                                                 phydbl *derivatives,
    #                                                 int     is_positive);
    # int     Compare_Two_States(char *state1, char *state2, int state_size);
    # void    Copy_One_State(char *from, char *to, int state_size);
    # void    Copy_Dist(phydbl **cpy, phydbl **orig, int n);
    # t_mod  *Copy_Model(t_mod *ori);
    # void    Record_Model(t_mod *ori, t_mod *cpy);
    void      Set_Defaults_Input(option *io)
    void      Set_Defaults_Model(t_mod *mod)
    void      Set_Defaults_Optimiz(t_opt *s_opt)
    # void    Test_Node_Table_Consistency(t_tree *tree);
    void      Get_Bip(t_node* a, t_node* d, t_tree* tree)
    void      Alloc_Bip(t_tree* tree)
    # int     Sort_Phydbl_Increase(const void *a, const void *b);
    # int     Sort_String(const void *a, const void *b);
    # phydbl  Compare_Bip(t_tree *tree1, t_tree *tree2, int on_existing_edges_only, int comparison_criterion);
    # void    Compare_Bip_Distance(t_tree *tree1, t_tree *tree2);
    # void    Match_Tip_Numbers(t_tree *tree1, t_tree *tree2);
    # void    Test_Multiple_Data_Set_Format(option *io);
    # int     Are_Compatible(char *statea, char *stateb, int stepsize, int datatype);
    # void    Hide_Ambiguities(calign *data);
    # void    Copy_Tree(t_tree *ori, t_tree *cpy);
    # void    Prune_Subtree(t_node *a, t_node *d, t_edge **target, t_edge **residual,
    #                     t_tree *tree);
    # void    Graft_Subtree(t_edge *target, t_node *link, t_node *link_daughter,
    #                     t_edge *residual, t_node *target_nd, t_tree *tree);
    # void    Reassign_Node_Nums(t_node *a, t_node *d, unsigned int *curr_ext_node,
    #                         unsigned int *curr_int_node, t_tree *tree);
    # void    Reassign_Edge_Nums(t_node *a, t_node *d, int *curr_br, t_tree *tree);
    # void    Find_Mutual_Direction(t_node *n1, t_node *n2, short int *dir_n1_to_n2,
    #                             short int *dir_n2_to_n1);
    # void    Update_Dir_To_Tips(t_node *a, t_node *d, t_tree *tree);
    # void    Fill_Dir_Table(t_tree *tree);
    # int     Get_Subtree_Size(t_node *a, t_node *d);
    # void    Init_Eigen_Struct(eigen *this);
    # phydbl  Triple_Dist(t_node *a, t_tree *tree);
    # phydbl  Triple_Dist_Approx(t_node *a, t_edge *b, t_tree *tree);
    # void    Make_Symmetric(phydbl **F, int size);
    # void    Divide_Mat_By_Vect(phydbl **F, phydbl *vect, int size);
    # void    Found_In_Subtree(t_node *a, t_node *d, t_node *target, int *match,
    #                         t_tree *tree);
    # void    Get_List_Of_Target_Edges(t_node *a, t_node *d, t_edge **list,
    #                                 int *list_size, t_tree *tree);
    # void    Fix_All(t_tree *tree);
    void      Record_Br_Len(t_tree *tree)
    void      Restore_Br_Len(t_tree *tree)
    # void    Get_Dist_Btw_Edges(t_node *a, t_node *d, t_tree *tree);
    # void    Detect_Polytomies(t_edge *b, phydbl l_thresh, t_tree *tree);
    # void    Get_List_Of_Nodes_In_Polytomy(t_node *a, t_node *d, t_node ***list,
    #                                     int *size_list);
    # void    Check_Path(t_node *a, t_node *d, t_node *target, t_tree *tree);
    # void    Connect_Two_Nodes(t_node *a, t_node *d);
    # void    Get_List_Of_Adjacent_Targets(t_node *a, t_node *d, t_node ***node_list,
    #                                     t_edge ***edge_list, int *list_size,
    #                                     int curr_depth, int max_depth);
    # void    Sort_List_Of_Adjacent_Targets(t_edge ***list, int list_size);
    # t_node *Common_Nodes_Btw_Two_Edges(t_edge *a, t_edge *b);
    # void    Random_Tree(t_tree *tree);
    # void    Reorganize_Edges_Given_Lk_Struct(t_tree *tree);
    # void    Random_NNI(int n_moves, t_tree *tree);
    # void    Fill_Missing_Dist(matrix *mat);
    # void    Fill_Missing_Dist_XY(int x, int y, matrix *mat);
    # phydbl  Least_Square_Missing_Dist_XY(int x, int y, phydbl dxy, matrix *mat);
    void      Check_Memory_Amount(t_tree *tree)
    # int     Get_State_From_P_Lk(phydbl *p_lk, int pos, t_tree *tree);
    # int     Get_State_From_P_Pars(short int *p_pars, int pos, t_tree *tree);
    # void    Check_Dirs(t_tree *tree);
    # void    Warn_And_Exit(const char *s);
    # void    Randomize_Sequence_Order(calign *cdata);
    # void    Update_Root_Pos(t_tree *tree);
    void      Add_Root(t_edge *target, t_tree *tree)
    # void    Update_Ancestors(t_node *a, t_node *d, t_edge *b, t_tree *tree);
    # #if (defined PHYTIME || defined SERGEII)
    # t_tree *Generate_Random_Tree_From_Scratch(int n_otu, int rooted);
    # #endif
    # void    Random_Lineage_Rates(t_node *a, t_node *d, t_edge *b, phydbl stick_prob,
    #                             phydbl *rates, int curr_rate, int n_rates,
    #                             t_tree *tree);
    # t_edge *Find_Edge_With_Label(char *label, t_tree *tree);
    # void    Site_Diversity(t_tree *tree);
    # void    Site_Diversity_Post(t_node *a, t_node *d, t_edge *b, t_tree *tree);
    # void    Site_Diversity_Pre(t_node *a, t_node *d, t_edge *b, t_tree *tree);
    # void    Subtree_Union(t_node *n, t_edge *b_fcus, t_tree *tree);
    # void    Binary_Decomposition(int value, int *bit_vect, int size);
    # void    Print_Diversity_Header(FILE *fp, t_tree *tree);
    # void    Best_Of_NNI_And_SPR(t_tree *tree);
    # int     Polint(phydbl *xa, phydbl *ya, int n, phydbl x, phydbl *y, phydbl *dy);
    t_tree*   Dist_And_BioNJ(calign *cdata, t_mod *mod, option *io)
    void      Add_BioNJ_Branch_Lengths(t_tree *tree, calign *cdata, t_mod *mod, matrix *mat)
    char*     Bootstrap_From_String(char* s_tree, calign *cdata, t_mod *mod, option *io)
    char*     aLRT_From_String(char *s_tree, calign *cdata, t_mod *mod, option *io)
    # void    Find_Common_Tips(t_tree *tree1, t_tree *tree2);
    phydbl    Get_Tree_Size(t_tree *tree)
    # void    Dist_To_Root_Pre(t_node *a, t_node *d, t_edge *b, t_tree *tree);
    # void    Dist_To_Root(t_tree *tree);
    char*     Basename(char* path)
    # t_node *Find_Lca_Pair_Of_Nodes(t_node *n1, t_node *n2, int *dist, t_tree *tree);
    # t_node *Find_Lca_Clade(t_node **node_list, int node_list_size, t_tree *tree);
    # int     Get_List_Of_Ancestors(t_node *ref_node, t_node **list, int *size, t_tree *tree);
    # int     Edge_Num_To_Node_Num(int edge_num, t_tree *tree);
    # void    Branch_Lengths_To_Rate_Lengths(t_tree *tree);
    # void    Branch_Lengths_To_Rate_Lengths_Pre(t_node *a, t_node *d, t_tree *tree);
    # int     Find_Clade(char **tax_name_list, int list_size, t_tree *tree);
    # void    Find_Clade_Pre(t_node *a, t_node *d, int *tax_num_list, int list_size, int *num, t_tree *tree);
    # t_edge *Find_Root_Edge(FILE *fp_input_tree, t_tree *tree);
    void      Copy_Tree_Topology_With_Labels(t_tree *ori, t_tree *cpy)
    void      Set_Model_Name(t_mod *mod)
    void      Adjust_Min_Diff_Lk(t_tree *tree)
    void      Translate_Tax_Names(char **tax_names, t_tree *tree)
    void      Skip_Comment(FILE *fp)
    void      Get_Best_Root_Position(t_tree *tree)
    void      Get_Best_Root_Position_Post(t_node *a, t_node *d, int *has_outgrp, t_tree *tree)
    void      Get_Best_Root_Position_Pre(t_node *a, t_node *d, t_tree *tree)
    void      Get_OutIn_Scores(t_node *a, t_node *d)
    int       Check_Sequence_Name(char *s)
    # int     Scale_Subtree_Height(t_node *a, phydbl K, phydbl floor, int *n_nodes, t_tree *tree);
    # void    Scale_Node_Heights_Post(t_node *a, t_node *d, phydbl K, phydbl floor, int *n_nodes, t_tree *tree);
    # int     Scale_Subtree_Rates(t_node *a, phydbl mult, int *n_nodes, t_tree *tree);
    void      Check_Br_Len_Bounds(t_tree *tree)
    # int     Scale_Subtree_Rates_Post(t_node *a, t_node *d, phydbl mult, int *n_nodes, t_tree *tree);
    # void    Get_Node_Ranks(t_tree *tree);
    # void    Get_Node_Ranks_Pre(t_node *a, t_node *d, t_tree *tree);
    void      Log_Br_Len(t_tree *tree)
    # phydbl  Diff_Lk_Norm_At_Given_Edge(t_edge *b, t_tree *tree);
    void      Adjust_Variances(t_tree *tree)
    # phydbl  Effective_Sample_Size(phydbl first_val, phydbl last_val, phydbl sum, phydbl sumsq, phydbl sumcurnext, int n);
    # void    Rescale_Free_Rate_Tree(t_tree *tree);
    phydbl    Rescale_Br_Len_Multiplier_Tree(t_tree *tree)
    phydbl    Unscale_Br_Len_Multiplier_Tree(t_tree *tree)
    phydbl    Reflect(phydbl x, phydbl l, phydbl u)
    int       Are_Equal(phydbl a, phydbl b, phydbl eps)
    bint      Check_Topo_Constraints(t_tree *big_tree, t_tree *small_tree)
    void      Prune_Tree(t_tree *big_tree, t_tree *small_tree)
    void      Match_Nodes_In_Small_Tree(t_tree *small_tree, t_tree *big_tree)
    void      Find_Surviving_Edges_In_Small_Tree(t_tree *small_tree, t_tree *big_tree)
    void      Find_Surviving_Edges_In_Small_Tree_Post(t_node *a, t_node *d, t_tree *small_tree, t_tree *big_tree)
    void      Set_Taxa_Id_Ranking(t_tree *tree)
    void      Get_Edge_Binary_Coding_Number(t_tree *tree)
    void      Make_Ancestral_Seq(t_tree *tree)
    void      Make_MutMap(t_tree *tree)
    int       Get_Mutmap_Val(int edge, int site, int mut, t_tree *tree)
    void      Get_Mutmap_Coord(int idx, int *edge, int *site, int *mut, t_tree *tree)
    void      Copy_Edge_Lengths(t_tree *to, t_tree *from_)
    void      Init_Scalar_Dbl(scalar_dbl *p)
    void      Init_Scalar_Int(scalar_int *p)
    void      Init_Vect_Dbl(int len, vect_dbl *p)
    void      Init_Vect_Int(int len, vect_int *p)
    char*     To_Lower_String(char *in_)
    phydbl    String_To_Dbl(char *string)
    int       String_To_Int(char *string)
    char*     To_Upper_String(char *in_)
    void      Connect_CSeqs_To_Nodes(calign *cdata, option *io, t_tree *tree)
    # void    Joint_Proba_States_Left_Right(phydbl *Pij, phydbl *p_lk_left,
    #                                     phydbl *p_lk_rght, vect_dbl *pi,
    #                                     int scale_left, int scale_rght, phydbl *F,
    #                                     int n, int site, t_tree *tree);
    void      Set_Both_Sides(bint yesno, t_tree *tree)
    # void    Set_D_States(calign *data, int datatype, int stepsize);
    # void    Path_Length(t_node *dep, t_node *arr, phydbl *len, t_tree *tree);
    # phydbl* Dist_Btw_Tips(t_tree *tree);
    # void    Best_Root_Position_IL_Model(t_tree *tree);
    void      Set_Br_Len_Var(t_edge *b, t_tree *tree)
    void      Check_Br_Lens(t_tree *tree)
    # void    Calculate_Number_Of_Diff_States_Post(t_node *a, t_node *d, t_edge *b,
    #                                             t_tree *tree);
    # void    Calculate_Number_Of_Diff_States_Pre(t_node *a, t_node *d, t_edge *b,
    #                                             t_tree *tree);
    # void    Calculate_Number_Of_Diff_States_Core(t_node *a, t_node *d, t_edge *b,
    #                                             t_tree *tree);
    # void    Calculate_Number_Of_Diff_States(t_tree *tree);
    # void    Build_Distrib_Number_Of_Diff_States_Under_Model(t_tree *tree);
    # int     Number_Of_Diff_States_One_Site(int site, t_tree *tree);
    # void    Number_Of_Diff_States_One_Site_Post(t_node *a, t_node *d, t_edge *b,
    #                                             int site, t_tree *tree);
    # int     Number_Of_Diff_States_One_Site_Core(t_node *a, t_node *d, t_edge *b,
    #                                             int site, t_tree *tree);
    # phydbl  Get_Lk(t_tree *tree);
    # phydbl  Get_d2Lk(t_tree *tree);
    # phydbl  Get_dLk(t_tree *tree);
    # void    Connect_Edges_To_Nodes_Serial(t_tree *tree);
    # phydbl  Mean_Identity(calign *data);
    # phydbl  Pairwise_Identity(int i, int j, calign *data);
    # phydbl  Fst(int i, int j, calign *data);
    # phydbl  Nucleotide_Diversity(calign *data);
    # void    Swap_Partial_Lk(t_edge *a, t_edge *b, int side_a, int side_b,
    #                         t_tree *tree);
    # scalar_dbl **Copy_Br_Len_Var(t_tree *mixt_tree);
    # scalar_dbl **Copy_Br_Len(t_tree *mixt_tree);
    # void         Transfer_Br_Len_To_Tree(scalar_dbl **bl, t_tree *tree);
    # void         Copy_Scalar_Dbl(scalar_dbl *from, scalar_dbl *to);
    # scalar_dbl  *Duplicate_Scalar_Dbl(scalar_dbl *from);
    # scalar_dbl  *Read_Weights(option *io);
    # phydbl       Scalar_Elem(int pos, scalar_dbl *scl);
    # int          Scalar_Len(scalar_dbl *scl);
    # int          Vect_Len(vect_dbl *scl);
    # void         Set_Scalar_Dbl(phydbl val, scalar_dbl *from);
    # void         Set_Scalar_Dbl_Min_Thresh(phydbl thresh, scalar_dbl *from);
    # void         Set_Scalar_Dbl_Max_Thresh(phydbl thresh, scalar_dbl *from);
    # void List_Of_Regraft_Nodes(t_node *a, t_node *d, phydbl time_thresh, t_ll *list,
    #                         t_tree *tree);
    # void Push_Bottom_Linked_List(void *what, t_ll **list, bint remove_duplicates);
    # void Remove_From_Linked_List(t_ll *elem, void *val, t_ll **list);
    # int  Linked_List_Len(t_ll *list);
    # void   *Linked_List_Elem(int pos, t_ll *ll);
    # void    Randomize_Tree(t_tree *tree, int n_prune_regraft);
    # t_ll   *Get_List_Of_Reachable_Tips(t_node *a, t_node *d, t_tree *tree);
    # void    Get_List_Of_Reachable_Tips_Post(t_node *a, t_node *d, t_ll **list,
    #                                         t_tree *tree);
    # phydbl  Length_Of_Path_Between_List_Of_Tips(t_ll *tips0, t_ll *tips1,
    #                                             matrix *mat);
    # void    Set_Update_Eigen_Lr(int yn, t_tree *tree);
    # void    Set_Use_Eigen_Lr(int yn, t_tree *tree);
    # void    Random_Walk_Along_Tree_On_Radius(t_node *a, t_node *d, t_edge *b,
    #                                         phydbl *radius, t_edge **target_edge,
    #                                         t_node **target_nd, phydbl *target_time,
    #                                         t_tree *tree);
    # void    Table_Top(unsigned int width);
    # void    Table_Row(unsigned int width);
    # void    Table_Bottom(unsigned int width);
    # t_cal  *Duplicate_Calib(t_cal *from);
    # t_clad *Duplicate_Clade(t_clad *from);
    # void    Swap_Partial_Lk_Extra(t_edge *b, t_node *d, int whichone, t_tree *tree);
    void      Remove_Duplicates(calign *data, option *io, t_tree *tree)
    # short int Are_Sequences_Identical(align *seq1, align *seq2);
    # char     *Mutation_Id(int mut_idx, t_tree *tree);
    # void      Random_Tax_Idx(t_node *a, t_node *d, int *idx, t_tree *tree);
    # void      List_Taxa_In_Clade(t_node *a, t_node *d, t_tree *tree);
    # void      Alias_Subpatt_Pre(t_node *a, t_node *d, t_tree *tree);
    # void      Alias_Subpatt_Post(t_node *a, t_node *d, t_tree *tree);
    # void      Alias_One_Subpatt(t_node *a, t_node *d, t_tree *tree);
    # void      Alias_Subpatt(t_tree *tree);
    # void   Map_Mutations(t_node *a, t_node *d, int sa, int sd, t_edge *b, int site,
    #                     int rcat, int *muttype, phydbl *muttime, int *muttax,
    #                     int *n_mut, t_tree *tree);
    void     Set_Update_Eigen(bint yesno, t_mod *mod);
    # int   *Order_Int(const int *u, const int n);
    # int   *Order_Dbl(const phydbl *u, const int n);
    # char   Integer_To_IUPAC_Code(int x);
    # void   Shuffle_Sites(const phydbl prop, align **data, const int n_otu);
    # void   Multiply_Scalar_Dbl(phydbl mult, scalar_dbl *x);
    void     Insert_Duplicates(t_tree* tree)
    # void   Get_Node_Ranks_From_Dist_To_Root(t_tree *tree);
    # void   Get_Node_Ranks_From_Times(t_tree *tree);
    # void   Get_Node_Ranks_From_Tip_Times(t_tree *tree);
    # phydbl Tree_Height(t_tree *tree);
    # void   Post_Inflate_Times_To_Get_Reasonnable_Edge_Lengths(t_node *a, t_node *d,
    #                                                         t_edge *b, phydbl min_l,
    #                                                         t_tree *tree);
    # void  Inflate_Times_To_Get_Reasonnable_Edge_Lengths(phydbl min_l, t_tree *tree);
    # void  Refactor_Tree(t_tree *tree);
    # void  Refactor_External(t_node *a, t_node *d, int *idx, t_tree *tree);
    # void  Refactor_Internal(t_node *a, t_node *d, t_edge *b, int *idx_nd,
    #                         int *idx_br, t_tree *tree);
    # int  *Integer_To_Bit(int val, const int ns);
    # char *Bit_To_Character_String(int *bit, int ns);
    # t_tree *Duplicate_Tree(t_tree *ori);
    # matrix *K80_dist(calign *data, phydbl g_shape);
    # matrix *JC69_Dist(calign *data, t_mod *mod);
    # matrix *Hamming_Dist(calign *data, t_mod *mod);
    # phydbl  Tree_Length(t_tree *tree);
    # void    Remove_Duplicates_From_Tree(calign *data, t_tree *tree);
    # void    Reset_Lk(t_tree *tree);
    # void    Set_Lk(t_tree *tree);
    # void    Reset_Prior(t_tree *tree);
    # void    Set_Prior(t_tree *tree);
    # t_edge *Get_Edge(t_node *a, t_node *d, t_tree *tree);
    # void Exchange_Nodes(t_node *a, t_node *d, t_node *w, t_node *v, t_tree *tree);
    # void Init_T_Beg(t_tree *tree);
    # void Set_Ignore_Root(bint yesno, t_tree *tree);
    # void Set_Bl_From_Rt(bint yesno, t_tree *tree);
    # void Replace_Short_With_Long_Tax_Names(t_tree *tree, option *io);
    # void Convert_Lengths_From_Calendar_To_Substitutions(t_tree *tree);
    # void Convert_Lengths_From_Calendar_To_Substitutions_Post(t_node *a, t_node *d,
    #                                                         t_tree *tree);
    # t_label *Get_Next_Label(t_label *curr_lab);
    # int      Scale_Subtree_Veloc(t_node *a, phydbl mult, int *n_nodes, int dim,
    #                             t_tree *tree);
    # int   Scale_Subtree_Veloc_Post(t_node *a, t_node *d, phydbl mult, int *n_nodes,
    #                             int dim, t_tree *tree);
    # void  Label_Edges(t_tree *tree);
    # void  Label_Nodes_With_Velocities(t_tree *tree);
    # void  Label_Nodes_With_Locations(t_tree *tree);
    # void  Edge_Labels_To_Rates(t_tree *tree);
    # void  Node_Labels_To_Velocities(t_tree *tree);
    # void  Node_Labels_To_Locations(t_tree *tree);
    # int   Add_Subtree_Veloc(t_node *a, phydbl add, int *n_nodes, int dim,
    #                         t_tree *tree);
    # int   Add_Subtree_Veloc_Post(t_node *a, t_node *d, phydbl add, int *n_nodes,
    #                             int dim, t_tree *tree);
    # int   Contmod_Start(short int datatype, short int dim_idx, t_tree *tree);
    # t_ll *Get_Velocity_Targets(t_node *a, t_node *d, t_tree *tree);
    # void Get_Velocity_Targets_Post(t_node *a, t_node *d, t_ll **list, t_tree *tree);
    # char D_State_To_Character(int d_state, t_tree *tree);
    # void ROC(phydbl *probs, short int *truth, int nclasses, int n_obs,
    #         phydbl *weights, char *tag, FILE *fp);
    # phydbl AIC(t_tree *tree);
    # phydbl BIC(t_tree *tree);
    # void   Set_Edge_Length_Optimizer(t_tree *tree);
    # int    Number_Of_Free_Params(t_tree *mixt_tree);
