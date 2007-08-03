
/* Interface for writing out result files in ICEM CFD format.
 *
 * Copyright (c) 1996 ICEM CFD Engineering.
 * Author:  Wayne A. Christopher
 *          ICEM CFD Engineering, Berkeley, California
 *          wayne@icemcfd.com
 */

#ifndef WRITE_RESULTS_H
#define WRITE_RESULTS_H

#ifdef __cplusplus
extern "C" {
#endif

/* There is one of these for each structured domain in the mesh. */

struct DFStructData {
    int imax[3];
    int nsubfaces;
    struct DFSubfaceData* subfaces;
};

/* Each subface needs this record. */

struct DFSubfaceData {
    int family;
    int imin[3];
    int imax[3];
};

/* There is one of these for each distinct type of cells in the mesh. */

struct DFCellData {
    int type;   /* see below */
    int num_cells;
    int* cell_info;  /* For each cell, the family and then all the nodes. */
};

/* There is one of these for each distinct type of face in the mesh --
 * i.e. TRI_3 and QUAD_4.
 */

struct DFFaceData {
    int type;   /* see below */
    int num_faces;
    int* face_info;  /* For each face, the family and then all the nodes. */
};

/* There is one of these for each family that the faces can be in.  Note that
 * families[n] is the data for family n.
 */

struct DFFamilyData {
    char* name;
};

/* One for each variable in the file. */

struct DFResultData {
    int type;      /* One of the DF_*_RESULT flags below. */
    char* name;
    int length;    /* Length of the data (only for consistency checking). */
    int prec;      /* Either 1 or 2 for single or double precision. */
    double* d_data;
    float* f_data;
};

/* Name/value pairs that can be anything you want to convey to
 * the postprocessor.
 */

struct DFPropertyData {
    char* name;
    char* value;
};

/* These are the possible values for the cell and face types, above. */

#ifndef TETRA_4
#define TETRA_4		4	/* cell */
#define HEXA_8		6	/* cell */
#define TRI_3		8	/* face */
#define PENTA_6		12	/* cell */
#define QUAD_4		14	/* face */
#define PYRA_5		18	/* cell */
#endif

#ifndef DF_NODE_RESULT
#define DF_NODE_RESULT 0
#define DF_CELL_RESULT 1
#define DF_FACE_RESULT 2
#endif

/* This function returns 1 if the file was writen and 0 if not.
 * Note that the cells must come in the order that the results for them
 * appear, if there are cell-based results.
 */

extern int
df_write_results(char* filename, char* title,
		 int nstructdata, struct DFStructData* structdata,
		 int ncelldata, struct DFCellData* celldata,
		 int nfacedata, struct DFFaceData* facedata,
		 int nnodes, double* nodes,
		 int nfamilies, struct DFFamilyData* families,
		 int nresults, struct DFResultData* results,
		 int nprops, struct DFPropertyData* props);

extern int df_add_results(char* filename, int nresults,
			  struct DFResultData* results);
    
extern int df_add_results_handle(int dom, int nresults,
			  struct DFResultData* results);

/* If you set this before calling df_write_results you will get a bunch
 * of (slow) checking done on the data you pass in.  Also it will do malloc
 * checking on the SGI using mallopt(M_DEBUG, 1) but it turns it off before
 * returning.
 */

extern int df_write_results_debug;

#ifdef __cplusplus
}
#endif

#endif
