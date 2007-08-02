#ifdef __cplusplus
extern "C"
{
#endif
/*
 * topolib_c.c
 */
int tp_open(char *fname);
void tp_close(void);
int tp_read(void);
int tp_init(char *fname);
void tp_quad_elements(void);
int tp_ndom(void);
int tp_nfac(void);
int tp_nsfac(void);
int tp_nar(void);
int tp_nvert(void);
int tp_nfsf(void);
int tp_nloop(void);
int tp_npc(void);
int tp_npp(void);
int tp_dom(int domain, int face);
int tp_domain_type(int domain);
int tp_tcor(int face, int domain);
int tp_ndom_nsub_ad(int idom);
int tp_ndomad(int idom, int *ier);
int tp_nsubad(int idom, int *ier);
int tp_node_position(int idom, int i, int j, int k);
int tp_domain_range(int idom, int imin[3], int imax[3]);
int tp_gsdomad(int idom, int *domface, int *subface);
int tp_domad(int idom, int dim, int *par_1, int *par_2,int *par_3,int *par_4);
int tp_domad_record(int idom, int line, int domad_record[5]);
int tp_ndim(int idom, int iface, int ndim[9]);
int tp_subface_range(int idom, int isub, int *imin, int *imax);
int tp_subface_points(int idom, int isub, int *ifirst, int *ilast);
int tp_subface_side(int idom, int isub, int *iside);
int tp_face_side(int idom, int iface, int *iside);
int tp_face_name(int idom, int isub, char *cface);
int tp_domains_of_subface(int isub, int *idomnum, int domain_list[2]);
int tp_domains_of_edge(int iedge, int *domains);
int tp_gtior(int idom, int isub, int *ior);
int tp_gcldims(int idom, int isub, int ncldim[9]);
int tp_gsfdims(int idom, int isub, int nsfdim[3][3]);
int tp_subface_3_points(int idom, int isub, int p1[3],int p2[3],int p3[3]);
int tp_edgdims(int idom, int isub, int iedge, int edgdim[3][3]);
int tp_subface_of_face(int idom, int iface, int *subfaces);
int tp_subface_of_edge(int iedge, int *subfaces);
int tp_edge_of_subface(int isub, int iedge, int *edges);
int tp_nedge_of_lin_domain(int idom);
int tp_edge_of_lin_domain(int idom, int *edges);
int tp_vtx_of_lin_domain(int idom, int *vtx);
int tp_vertex(int iedge, int vertex[2]);
int tp_vertex_of_domain(int idom, int vtx[8]);
int tp_vertex_pnt(int idom, int ivert, int kpnt[3]);
int tp_n_shared_sf(void);
void tp_shared_sf(int *subfaces);
int tp_domain_face_range(int idom, int iface, int imin[3], int imax[3]);
int tp_edge_entity_range(int idom, int iedge, int imin[3], int imax[3]);
int tp_edge_points(int idom, int iedge, int pnt1[3], int pnt2[3]);
int tp_edges_of_domain(int idom, int *edges);
int tp_entity_name(int iedge, char *name);
int tp_entity_number(char *name);
int tp_edge_in_subface(int isub, int iedge, int ndim[2][3]);
int tp_vertex_in_edge(int iedge, int ivtx, int ndim[2]);
int tp_subface_dimension(int isub, int dimension[2]);
int tp_subface_bc_code(int isub, int *bc_code);
int tp_create_edge_table(void);
void tp_free_edge_table(void);
/*
 *  global_coord.c 
 */
int tp_global_coord(int idom, int *gl_coord, int dimension);
int tp_global_coord_list(int ndomsel, int *domsel, int dim, int *koord);
/*
 *  free_subface.c
 */
int tp_get_vertex(int iv, int iedge);
int tp_fsf_edge_numbers(int ifsub);
int tp_fsf_entl(int ifsub, int ied);
int tp_fsf_loop_numbers(int ifsf, int *loarray);
int tp_fsf_data(int ifsf, int info, int *idata);
int tp_fsf_loop_status(int ifsf, int *lsarray);
int tp_fsf_curve_numbers(int ifsf, int *pcarray);
int tp_fsf_point_data(int ifsf, float *pparray);
/*
 * twin_subface.c
 */
int tp_entity_group(int entity_type, int num, int *group_number);
int tp_twin_subfaces(int subface_no, int *n_twin, int *list);
int tp_twin_edges(int edge_enity_no, int *n_twin, int *list);
/*
 * periodic.c
 */
int tp_n_periodic_pair(int entity_type);
int tp_periodic_data(int entity_type, int npair, int *list1, int *list2);

void tp_error_exit(void);
void tp_error_print(void);

#ifdef __cplusplus
}
#endif
