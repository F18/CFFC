#ifndef GLOBAL_DOMAIN_H
#define GLOBAL_DOMAIN_H
#ifdef __cplusplus
extern "C" {
#endif

int df_add_material(int file_no, char *name, int material_no);

int df_close_global();

int df_write_volume(int file_no,
        unsigned short ijk[3], unsigned short material, char boundary,
        int n_edge, int *edge_axis, double *tees, 
			double *normals, int *edge_idents,
        int n_face, int *face_axis, double *xys, int *face_idents,
	int n_p_points, double *p_point_locs, int *p_point_idents);

int
df_write_stations(int file_no, int dim, int n_stations, double *stations);


int df_global_grid_type(int fn, int *type);

int df_read_stations(int file_no, int dim, double *stations);

int df_global_nijk(int fn, int *nijk);

int df_read_volume(int file_no, int *IV,
        unsigned short ijk[3], unsigned short *material, char *boundary,
        int *n_edge, int *edge_axis, double *tees, 
			double *normals, int *edge_idents,
        int *n_face, int *face_axis, double *xys, int *face_idents,
	int *n_p_point_out, double *p_point_locs, int *p_point_idents);

int df_set_global_grid_type(int fn, int type);

int df_global_grid_type(int fn, int *type);

#ifdef __cplusplus
}
#endif

#define CARTESIAN_DOMAIN	0
#define CYLINDRICAL_DOMAIN	1

#define BOUNDARY_MATERIAL	0xffff
#endif
