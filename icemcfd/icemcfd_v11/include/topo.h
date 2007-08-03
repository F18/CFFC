#include "../util/fortran.h"

#ifdef __cplusplus
extern "C"
{
#endif

int tp_read_file(void);
int tp_parse(char *string, int num, int digit, int *value);
void topo_error(char *strng);
void tp_warning_print();
int tp_adrdom_3d(int idom, int jdom, int nvect[3]);
int tp_adrdom_2d(int idom, int jdom, int nvect[3]);
void topo_warning(char *strng);


void FMNAME(topopen, TOPOPEN) (STR_PSTR(fname), int *len, int *iunit, 
				int *ier STR_PLEN(fname));
void FMNAME(topclse, TOPCLSE) (int *iunit);
void FMNAME(topread, TOPREAD) (int *iunit, int *ier);
void FMNAME(topinit, TOPINIT) (STR_PSTR(fname), int *len, int *ier STR_PLEN(fname));
void FMNAME(tprdnow, TPRDNOW) (int *iunit, int *ier);
void FMNAME(topquad, TOPQUAD) (void);
int FMNAME(gtndom, GTNDOM) (void);
int FMNAME(gtnfac, GTNFAC) (void);
int FMNAME(gtnsfac, GTNSFAC) (void);
int FMNAME(gtnar, GTNAR) (void);
int FMNAME(gtvtx, GTVTX) (void);
int FMNAME(gtnfsf, GTNFSF) (void);
int FMNAME(gtnloop, GTNLOOP) (void);
int FMNAME(gtnpc, GTNPC) (void);
int FMNAME(gtnpp, GTNPP) (void);
int FMNAME(gtdom, GTDOM) (int *face, int *domain);
int FMNAME(ntypdeg, NTYPDEG) (int *idom);
int FMNAME(ntcor, NTCOR) (int *iface, int *idom);
int FMNAME(ndomad, NDOMAD) (int *idom, int *ier);
int FMNAME(gndomad, GNDOMAD) (int *idom, int *ier);
int FMNAME(gnsubad, GNSUBAD) (int *idom, int *ier);
int FMNAME(ndposit, NDPOSIT) (int *idom, int *i, int *j, int *k);
int FMNAME(gtvertx, GTVERTX) (int *iv, int *iedge);
int FMNAME(gtneds, GTNEDS) (int *ifsub);
int FMNAME(gtentl, GTENTL) (int *ifsub, int *ied);
void FMNAME(dmrange, DMRANGE) (int *idom, int *imin, int *imax, int *ier);
void FMNAME(gsdomad, GSDOMAD) (int *idom, int *nadsub, int *domface, int *subface, int *ier);
void FMNAME(gdomad, GDOMAD) (int *idom, int *nadsub, int *domface, int *subface,
                             int *domain, int *facead, int *ier);
void FMNAME(gsubfad, GSUBFAD) (int *idom, int *nedge, int *sfofdom, int *domadj,
                               int *sfadj, int *cedge, int *ier);
void FMNAME(domadv, DOMADV) (int *idom, int *line, int *domad_record, int *ier);
void FMNAME(gndims, GNDIMS) (int *idom, int *iface, int *ndim, int *ier);
void FMNAME(gtsreg, GTSREG) (int *idom, int *isub, int *imin, int *imax, int *ier);
void FMNAME(sregion, SREGION) (int *idom, int *isub, int *ifirst, int *ilast, int *ier);
void FMNAME(gtside, GTSIDE) (int *idom, int *isub, int *iside, int *ier);
void FMNAME(gtfside, GTFSIDE) (int *idom, int *iface, int *iside, int *ier);
void FMNAME(gtface, GTFACE) (int *idom, int *isub, STR_PSTR(cface), int *ier STR_PLEN(cface));
void FMNAME(gtdomnr, GTDOMNR) (int *iface, int *idomnum, int *idomv, int *ier);
void FMNAME(gtior, GTIOR) (int *idom, int *isub, int *ior, int *ier);
void FMNAME(gcldims, GCLDIMS) (int *idom, int *isub, int *ncldim, int *ier);
void FMNAME(gsfdims, GSFDIMS) (int *idom, int *isub, int *nsfdim, int *ier);
void FMNAME(sfspann, SFSPANN) (int *idom, int *isub, int p1[3], int p2[3], int p3[3],int *ier);
void FMNAME(sfoface, SFOFACE) (int *idom, int *iface, int *nsub, int *nsuba, int *ier);
void FMNAME(sfofed, SFOFED) (int *iedge, int *nsub, int *nsuba, int *ier);
void FMNAME(edgofsf, EDGOFSF) (int *isub, int *iedge, int *nedgee, int *nedgeea, int *ier);
void FMNAME(edofldm, EDOFLDM) (int *idom, int *nedgee, int *nedgeea, int *ier);
void FMNAME(vtofldm, VTOFLDM) (int *idom, int *vtx, int *ier);
void FMNAME(edgdims, EDGDIMS) (int *idom, int *isub, int *iedge, int *ndim, int *ier);
void FMNAME(vertex, VERTEX) (int *iedge, int *ivert);
void FMNAME(vtxpnt, VTXPNT) (int *idom, int *ivert, int *kpnt, int *ier);
void FMNAME(vtofdom, VTOFDOM) (int *idom, int vtx[8], int *ier);
void FMNAME(blintsf, BLINTSF) (int *n, int *shared_sf, int *ier);
void FMNAME(dfrange, DFRANGE) (int *idom, int *iface, int imin[3], int imax[3], int *ier);
void FMNAME(edgrnge, EDGRNGE) (int *idom, int *iedge, int imin[3], int imax[3], int *ier);
void FMNAME(edgpnts, EDGPNTS) (int *idom, int *iedge, int pnt1[3], int pnt2[3], int *ier);
void FMNAME(edgofdm, EDGOFDM) (int *idom, int *nedge, int *nedgea, int *ier);
void FMNAME(entname, ENTNAME) (int *iedge, STR_PSTR(name), int *ier STR_PLEN(name));
void FMNAME(entnr, ENTNR) (STR_PSTR(name), int *iedge STR_PLEN(name));
void FMNAME(edginsf, EDGINSF) (int *isub, int *iedge, int n_dim[6], int *ier);
void FMNAME(vtxined, VTXINED) (int *iedge, int *ivtx, int ndim[2], int *ier);
void FMNAME(gtsfacd, GTSFACD) (int *isub, int ndim[2], int *ier);
void FMNAME(gtsfcode, GTSFCODE) (int *isub, int *code, int *ier);
void FMNAME(gloco, GLOCO) (int *idom, int gl_coord[3], int *ier);
void FMNAME(gloco2d, GLOCO2D) (int *idom, int gl_coord[3], int *ier);
void FMNAME(glocox, GLOCOX) (int *ndom, int *domsel, int *dimension,int *koord,int *ier);
void FMNAME(gtfsflo, GTFSFLO) (int *ifsf, int *loarray, int *ier);
void FMNAME(gtfsfda, GTFSFDA) (int *ifsf, int *info, int *idata, int *ier);
void FMNAME(gtfsfls, GTFSFLS) (int *ifsf, int *lsarray, int *ier);
void FMNAME(gtfsfpc, GTFSFPC) (int *ifsf, int *pcarray, int *ier);
void FMNAME(gtfsfpp, GTFSFPP) (int *ifsf, float *pparray, int *ier);
void FMNAME(gtentgr, GTENTGR) (int *enttype, int *num, int *group, int *ier);
void FMNAME(gttwsf, GTTWSF) (int *isub, int *n_twin, int *list, int *ier);
void FMNAME(gttwedg, GTTWEDG) (int *iedge, int *n_twin, int *list, int *ier);
void FMNAME(gtnpair, GTNPAIR) (int *enttype, int *npair);
void FMNAME(gtperdt, GTPERDT) (int *enttype, int *npair, int *list1, int *list2, int *ier);

#ifdef __cplusplus
}
#endif
