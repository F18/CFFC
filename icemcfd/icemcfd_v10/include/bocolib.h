#ifdef __cplusplus
extern "C"
{
#endif
/*
 * bocolib_c.c
 */
int bc_init(char *bocofile);

int bc_open(char *bocofile);
void bc_close(void);
int bc_read(void);
int bc_num(int mode, int *num, int *numbc);
int bc_info(int mode, int num, int *nposib, int *nsiboca);
int bc_nboco(int mode, int num, int *ier);

int bc_get_values(int iboco, float bcvalue[10]);
int bc_get_flags(int iboco, char bcflag1[6], char bcflag2[6]);
int bc_get_data(int mode, int num, int ibc, int **types, void ***values);

int bc_value(int mode, int num, int ibc, float bcvalue[10]);
int bc_flag(int mode, int num, int ibc, char bcflag1[6], char bcflag2[6]);

int bc_n_analysis_datas(void);
int bc_flags_analysis_data(int iandat, char anflag1[6], char anflag2[6]);
int bc_values_analysis_data(int iandat, float advalue[10]);
int bc_nbtotal(void);
int bc_nbsize(void);
int bc_pcube(void);
void bc_error_exit(void);
void bc_error_print(void);
char* bc_get_error();
int bc_merge_fambocos(char **fam_bocos, int nfambocos, char *bocofile);

#ifdef __cplusplus
}
#endif


/*
 * boco.h
 */
#include "boco.h"
