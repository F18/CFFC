#ifndef DOMAIN_H
#define DOMAIN_H
/*
 * domain.h included in user's files
 */
#define DF_ERROR	-1
/*
 * modes for opening domain file
 */
#define	MODE_READ 	0
#define MODE_WRITE	1
#define MODE_CLOSED	2
#define MODE_MODIFY	3
#define MODE_SLEEP	4
/*
 * domain file types
 */

#define UNDEFINED_DOMAIN			-1
#define GLOBAL_DOMAIN                           0
#define STRUCTURED_DOMAIN                       1
#define UNSTRUCTURED_DOMAIN			2

/*
 * coodinate systems
 */
#define CARTESIAN_SYSTEM                0
#define CYLINDRICAL_SYSTEM              1
#define SPHERICAL_SYSTEM                2
#define TOROIDAL_SYSTEM                 3

/*
 * sub domain types (value of type field in general_domain)
 */
#define RANGE_DOMAIN    0
#define EXPLICIT_DOMAIN 1
#define REAL_RESULTS	2
#define DOUBLE_RESULTS	3
#define BY_PID		4
#define STRING_DOMAIN	5

#define MAX_NAME_LENGTH 500

/*
 * returns an integer ( >= 0) identifying the domain 
 * else returns -1 for an error opening the file
 */

/* create a new domain file */
#ifdef __cplusplus
extern "C" {
#endif

void df_get_date(int *date);

int
df_add_struct_domain(int file_number, const char *domain_name,
                int istart, int skip[3], int emin[3], int emax[3], int pid);

int df_set_coordinate_system(int dn, int c_system);
int df_coordinate_system(int dn, int *c_system);

const char *df_get_error(void);

int df_add_primary_domain(int file_no, const char *name, int imin[3], int imax[3]);

int DF_Add_Secondary_Domain(int file_no, const char *name, 
					int imin[3], int imax[3], int pid);

int df_add_secondary_domain(int file_no, const char *name, int imin[3], int imax[3]);

int df_open(const char *filename, int mode, int type);

int df_close(int file_number);

int df_flush(int file_number);

int df_n_domain(int file_no, int *n_domain);

const char *df_element_name(int type);

int df_element_type_id(const char *name);

int df_type(int dn, int *type);

int df_date(int dn, int *mod_date);

char *df_ascii_time(int df);
    
int df_n_elements(int fn, int *n_el);

int df_get_domain(int dom_file_number, const char *domain_name, int *sub_dom_number);

int df_domain_name(int file_number, int domain_number, const char **name);

int df_domain_dimension(int fn, int dn, int *dim);

int df_struct_domain_range(int file_number, int domain_number,
                                        int imin[3], int imax[3]);

int df_struct_read_nodes(int file_number, int domain_number,
                int imin[3], int imax[3], double *pnt);

int df_struct_rewrite_nodes(int file_number, int domain_number,
                int imin[3], int imax[3], double *pnt);

int df_struct_degree(int file_no, int *deg);

int df_struct_set_degree(int file_no, int deg);

int df_n_nodes(int file_no, int *n_nodes);

int df_unstruct_read_nodes(int file_number, int start_no, 
		int npnts, double *pnt);

int df_unstruct_update_nodes(int file_number, int start_no, int npnts,
                        double *pnt);
int df_read_elements(int file_no, int el_no, int n_el, int *data);

int df_write_elements(int file_no, int n_el, int *data, int type);

int
df_update_elements(int file_no, int el_no, int n_el, int *data);


int
df_domain_type(int file_no, int domain_no, int *type);


int
df_n_members(int file_no, int domain_no, int *nmem);

int
df_get_range(int file_no, int domain_no, int *start, int *end);

int
df_get_membership(int file_no, int domain_no,
		  short *is_member, int mstart, int size, short mark);
int
df_get_members(int file_no, int domain_no, int start, int n, int *members);

int
df_unstruct_write_nodes(int file_number, int npnts,  double *pnt);

int
df_struct_write_nodes(int file_number, int npnts,  double *pnt);

int
df_add_unstruct_domain(int file_no, const char *name, int dim, int start, int end);

int df_set_family_pid(int file_no, const char *name, int pid);

int df_set_family_pids(int file_no, int npids, char **names, int *pids);

int df_get_family_pid(int file_no, const char *name, int *pid);

int df_get_pid_family(int file_no, int pid, const char **name);

/* Grab all the family names in one pass for efficiency */
int df_get_all_family_names(int file_no, char** names, int namesize);
    
/* given a domain number, return the pid that is associated with it */
int df_get_pid_of_domain(int file_no, int domain_no, int *pid);
    
/* return the max pid of any family named in the file (not necessarily used
 * by an element
 */
int df_max_pid(int file_no);
    
int df_get_string(int file_no, int dn, const char **value);

int df_add_string_domain(int file_no, const char *name, const char *string);

/* result data stuff */

int df_add_unstruct_result_data(int file_no, const char *name, 
		int start, int end, 
		int type, void *data);

int df_add_struct_result_data(int file_no, const char *name,
		int start[3], int end[3],
		int type, void *data);
		

int df_result_data_type(int file_no, int dn, int *type);

int df_read_result_data(int file_no, int dn, void *data);



int df_nodes_per_element(int type);


void display_element(int type, float (*)[3], int solid, int color);
/* color is specified in rgb */
void display_element_rgb(int type, float (*)[3], int solid, int color);

void display_element_shrink_d(
	int type, double (*)[3], int solid, int color, int shrnk);
void display_element_shrink(
	int type, float (*)[3], int solid, int color, int shrnk);
void display_element_shrink_rgb(
	int type, float (*)[3], int solid, int color, int shrnk);

int df_n_sections(int file_no, int *n_sections);

int
df_section_info(int file_no, int section_number,
                int *start, int *end, int *type);

int
df_domain_element_type(int file_no, int domain_no, int *type);

int df_domain_material(int fn, int dn, int *mat);

int df_add_domain_members(int file_no, int domain_no,
                                     int nmem, int* members);

int df_get_node_dimension(int file_no, short *dims);

int df_set_node_dimension(int file_no, short *dims);

int df_create_connectivity(int fn);

void df_delete_connectivity(int fn);

int df_neighbors(int fn, int face_type, int max_ret,
                        int face_nodes[], int elements[]);

int df_read_elnodes(int fn, int en, int *npe, int *el_type, int *pid,
                int *ext_num, int **el_nodes);

int df_file_version(int fn, int *vers);

int df_lib_version(int *vers);

int
df_element_nodes(int fn, int en, int *npe, int *el_type, int **el_nodes);


int DF_Struct_Write_Nodes(int file_number, int npnts, double *pnt, int *pid);

int DF_Struct_Rewrite_Nodes(int file_number, int domain_number,
                int imin[3], int imax[3], double *pnt, int *pid);
int DF_Write_Elements(int file_no, int n_el,
                        int *data, int numbered, int type);
int DF_Read_Elements(int file_no, int el_no, int n_el, int *data,
                        int numbered);
int DF_Unstruct_Write_Nodes(int file_number, int npnts,  double *pnt,
                        int *nos);
int DF_Struct_Read_Nodes(int file_number, int domain_number,
                int imin[3], int imax[3], double *pnt, int *pid);

int DF_Unstruct_Read_Nodes(int file_number, int start_no, int npnts,
                        double *pnt, int *nos);
const char *df_node_names(int dim);

#ifdef __cplusplus
}
#endif

#endif
