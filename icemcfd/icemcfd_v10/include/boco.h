

/* the variables type (below) take on these values */

#define BOCO_DOMAIN     1
#define BOCO_SUBFACE    2
#define BOCO_ENTITY     3
#define BOCO_VERTEX     4
#define BOCO_PATCH      5
#define BOCO_FAMILY     6
#define BOCO_N_TYPES    6

#define BOCO_STRING 1
#define BOCO_FLOAT 2

#define NREAL_PER_BOCO  10
#ifdef   __cplusplus

int read_boco(const char *bocofile,
// function called for each grouping name 
        void (*name_func) (int type, int num, char *name, float vals[NREAL_PER_BOCO]),
// function called for each boundary condition 
        void (*bc_func) (int type, int num, char* famname,
                        float vals[NREAL_PER_BOCO], char *f1, char *f2));

int bc_get_solver_name(const char *filename, char *solvername);

extern "C" void bc_error_print();
#endif /*  __cplusplus */

