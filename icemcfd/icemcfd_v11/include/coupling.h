#ifdef __cplusplus
extern "C" {
#endif

/*
 * structured couplings 
 */

int df_add_coupling(int df, int ncoarse[2], int *coarse_nodes, int nfine[2],
			int *fine_nodes);

int df_get_couplings_n(int df, int cn, int ncoarse[2], int nfine[2]);

int df_get_coupling_coarse(int df, int cn, int *coarse_nodes);

int df_get_coupling_fine(int df, int cn, int *fine_nodes);

int df_get_n_couplings(int df, int *n_couplings);

/*
 * unstructured couplings 
 */

int df_add_unstruct_coupling(int df, int nf1, int *nodes1, 
				int nf2, int *nodes2);

int df_get_n_unstruct_couplings(int df, int *n_couplings);
/* nf1       - number of faces on side 1
 * nodes1[0]... nodes1[3] nodes in face 1 side 1
 * nodes1[4] ... nodes1[7] nodes in face 2 side 1
 * etc. and same for face 2
 */
/* returns the number of faces on one side of the coupling */
int df_get_unstruct_coupling_size(int df, int cn, int side, int *nf);

/* returns the nodes for one side of the coupling */
int df_get_unstruct_coupling(int df, int cn, int side, int *nodes);

#ifdef __cplusplus
}
#endif
