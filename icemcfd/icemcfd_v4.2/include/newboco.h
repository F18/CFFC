

#define MAX_NAME_LEN    200
#define TOPOTYPE  7
#define TOPOTYPE2 4

/* new data-structure */
extern struct family *fmy;
extern int nfam;

struct boco {
        int nstrings;
        int nnums;
        char **strings;
        float *nums;
};
struct topo {
        int mode;
        int ent;
};
struct family {
        char *name;
        int nboco;
        struct boco *bc;
        int ntopo;
        struct topo *tp;
};

#ifdef __cplusplus
extern "C"
{
#endif

/* Functions used in read_newboco.c to read family_boco and family_topo */
int bc_read_boco(void);
int bc_read_topo(void);
int bc_go_end_of_line(FILE *Bc_file);
int bc_read_line(FILE *Bc_file, char name[MAX_NAME_LEN], int *nnums, int *nstrings,
        float **nums, char ***strings);
int bc_read_string(FILE *Bc_file, int *nstrings, char ***strings);
int bc_read_number(FILE *Bc_file, int *nnums, float **nums);

/* Other internal functions */
int bc_open_both(char *family_boco, char *family_topo);
void boco_error(char *format, ...);
void bc_free_data();
void bc_type();

/* New functions */
int bc_init_both(char *family_boco, char *family_topo, int ndom, int nsf,
        int nedge, int nvtx, int npatch, int npfam, int nfreesf);
char* bc_get_error();
int bc_family_name(int mode, int num, char **family_name);
int bc_family_name_pid(int mode, int num, char *family_name);
int bc_family_pid(const char *family_name, int *pid);
int bc_pid_of_entity(int mode, int num, int *pid);
int bc_get_data(int mode, int num, int bcnum, int **types, void ***values);
int bc_oldformat(char *bocofile);

int write_fam_boco(char *filename);
int write_fam_topo(char *filename);

#ifdef __cplusplus
}
#endif

