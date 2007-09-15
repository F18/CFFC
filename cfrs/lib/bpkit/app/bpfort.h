c 
c                 BPKIT 2.0 Block Preconditioner Toolkit
c                   Authors:  E. Chow and M. A. Heroux
c     Copyright (c) 1995-1996  The Regents of the University of Minnesota

      integer BP_DENSE, BP_SPARSE
      integer BP_NONE, BP_BJACOBI, BP_BSOR, BP_BSSOR, BP_BILUK, BP_BTIF
      integer BP_LU, BP_INVERSE, BP_SVD, BP_RILUK, BP_ILUT
      integer BP_APINV_TRUNC, BP_APINV_BANDED, BP_APINV0, BP_APINVS
      integer BP_DIAG, BP_TRIDIAG, BP_SOR, BP_SSOR, BP_GMRES

      parameter (BP_DENSE = 0)
      parameter (BP_SPARSE = 1)

      parameter (BP_NONE = 0)
      parameter (BP_BJACOBI = 1)
      parameter (BP_BSOR = 2)
      parameter (BP_BSSOR = 3)
      parameter (BP_BILUK = 4)
      parameter (BP_BTIF = 5)

      parameter (BP_LU = 1)
      parameter (BP_INVERSE = 2)
      parameter (BP_SVD = 3)
      parameter (BP_RILUK = 4)
      parameter (BP_ILUT = 5)
      parameter (BP_APINV_TRUNC = 6)
      parameter (BP_APINV_BANDED = 7)
      parameter (BP_APINV0 = 8)
      parameter (BP_APINVS = 9)
      parameter (BP_DIAG = 10)
      parameter (BP_TRIDIAG = 11)
      parameter (BP_SOR = 12)
      parameter (BP_SSOR = 13)
      parameter (BP_GMRES = 14)
