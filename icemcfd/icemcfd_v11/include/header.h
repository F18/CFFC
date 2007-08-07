#ifndef HEADER_H
#define HEADER_H

#include "domain.h"

#include <stdio.h>
#include <time.h>


/* 
 * if the size of the file is larger than 4GB then
 * the version becomes 4
 */
#define DOMAIN_VERSION		3

typedef struct {
/*
 * structured domain in the domain file
 */
/* emin and emax are the external min and max values of the
 * struct_domain: that is, the user can request to
 * load values within these bounds */
	int emin[3];
	int emax[3];
/* an (i,j,k) triple given to this structured_domain
 * is mapped to the array of nodes by
 *
 * nodes[istart + i*skip[0] + j*skip[1] + k*skip[2]]
 */
	int istart;
	int skip[3];
/* property id of nodes in this domain */
	int pid;
} struct_sub_mesh;

/*
 * unstructured domain of the domain
 */

typedef struct {
	int dim; 
} unstruct_sub_mesh;


/*
 * material region of a global domain
 */
typedef struct {
	int material_number;
} global_sub_mesh;

typedef struct {
	int start, end;
} range;


typedef struct {
/* number of members */
	int n_members;
/* offset into file of list
this is now stored with the memory general domain
	int offset;
 */
} Explicit;

typedef struct {
	int pid;
} by_pid;


typedef union {
	range r;
	Explicit e;
	by_pid p;
	int string_len;
} list_union;
/*
 * subset of a mesh
 */
typedef union { 
	struct_sub_mesh s; /* used to be struct_ domain */
	unstruct_sub_mesh u; /* used to be unstruct_domain */
	global_sub_mesh g; /* used to be global_domain */
} mesh_subset;	 /* used to be any_domain */

typedef struct {
	int namelen;
/* information specific to each mesh type */
	mesh_subset specific;
/* how is the list stored */
	int type;
/* list of members */
	list_union list;
} file_addendum; /* used to tbe file_general_domain */

#include "dflong.h"

typedef struct {
	char *name;
/* string if a string addendum */
	char *value;
/* information specific to each mesh type */
	mesh_subset specific;
/* how is the list stored */
	int type;
/* list of members */
	list_union list;
/* offset to data immediately following addendum header */
	dflong offset;
} memory_addendum; /* used to be memory_general_domain */


/*
 * in the header file there is a section that contains information
 * specific to each mesh type. The following structures
 * define that specific information
 */
/*
 * fields only valid for UNSTRUCTURED omains
 */
typedef struct {
	int element_type;
	int n_elements;
} section;

typedef struct {
	int n_nodes;
	int n_sections;
	section *sections;
	dflong sections_offset;
} unstruct_header_info;


/*
 * fields only valid for GLOBAL meshes
 */
typedef struct {
/* cartesian or cylindrical */
	int type;
	int nijk[3];
} global_header_info;


/*
 * fields valid only for STRUCTURED meshes 
 */
typedef struct {
	int nijk[3];
/* linear or quadratic elements */
	int degree;
} struct_header_info;


/*
 * in memory copy of file header that starts every domian file
 */
typedef struct {
/* version of the domain file library that created this file */
	int version;
/* seconds since Jan. 1, 1970 */
	int modified_date; 
/* type: global, structured or unstructured */
	int type;
/* coordinate system */
	int c_system;
/* number addenda in the domain file */
	int n_addenda; /* used to be n_domains */
/* byte offset to addenda */
	dflong addenda_offset;
/* specific info for each mesh type */
	union { unstruct_header_info u;
		struct_header_info s;
		global_header_info g;
	} specific;
} domain_header;


size_t size_file_domain_header();

typedef struct {
	int ncoarse[2];
	int *coarse_nodes;
	int nfine[2];
	int *fine_nodes;
} coupling;

typedef int quad_face[4];

typedef struct {
	int n_faces[2];
	quad_face *faces[2];
} uns_coupling;

/*
 * in-memory copy of domain header plus some information
 */
typedef struct {
/* reading or writing */
	int mode;
/* file containing data */
	FILE *file;
/* name of file in case it goes into sleep mode */
	char *filename;
/* in-memory copy of file header */
	domain_header header;
/* pointers to in-memory copies of addenda */
	memory_addendum **addenda;
/* counters used in creation of various mesh types */
	int m, n;
/* number of couplings and couplings */
	int n_couplings;
	coupling *couplings;
	int coupling_sub_domain;

/* number of unstructured couplings */
	int n_uns_couplings;
	uns_coupling *uns_couplings;
	int uns_coupling_domain;

/* external identifier */
	int file_number;
} domain_file;


/* retrieve an open domain file from the list of currently open
 * domains */
domain_file *df_get_file(int domain_number);


void domain_error(char *msg, ... );

memory_addendum *df_add_domain(int file_no, const char *name);


int df_close(int file_number);

int df_struct_read(FILE *file, domain_header *header);
int df_global_read(FILE *file, domain_header *header);
int df_unstruct_read(FILE *file, domain_header *header);
void df_unstruct_init(domain_header *dh);
void df_global_init(domain_header *dh);
void df_struct_init(domain_header *dh);
int df_close_write(int file_number);
int df_close_unstruct(domain_file *df);
int df_close_global(domain_file *df);
int df_write_domain(domain_file *df, memory_addendum *orig, int do_seek);
int df_close_struct(domain_file *df);
int df_write_sections(domain_file *df);
int df_write_couplings(domain_file *df);

#include "fwrite.h"

#ifndef SEEK_SET
#define SEEK_SET	0
#define SEEK_CUR	1
#endif

int fread_header(FILE *file, domain_header *header);
int fwrite_header(FILE *file, domain_header *header);


dflong df_ftell(FILE *f);

#ifdef __alpha
int fseeko(FILE *stream, dflong val, int type);
#endif

dflong df_section_offset(unstruct_header_info *hi);

#endif

