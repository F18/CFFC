      program f1example
      implicit none
c-----------------------------------------------------------------------
c     Example program that calls BPKIT library from FORTRAN.
c-----------------------------------------------------------------------
c     Hardcoded for SHERMAN1 matrix from HB collection.
c     With the original .BpResource, the output should be:
c                 88      9.520099361807e-09
c-----------------------------------------------------------------------
      include  'bpfort.h'

      integer   nmax, nzmax
      parameter (nmax =  1000, nzmax = 3750)
      integer   ia(nmax+1), ja(nzmax), kvst(nmax+1)
      real*8    a(nzmax), rhs(nmax), sol(nmax), exact(nmax)
      real*8    temp(nmax*3)
      integer   i, n, nnz, nrhs
      integer*8 bmat, precon

      ! read matrix
      n = nmax
      nnz = nzmax
      nrhs = nmax*3
      call readhb('SHERMAN1', 8, n, nnz, nrhs, a, ja, ia,
     *    rhs, sol, exact, temp)

      do i = 1, n
	  rhs(i) = temp(i)   ! right-hand side
	  sol(i) = 0.d0
      enddo

      do i = 1, 11
	  kvst(i) = (i-1)*100 + 1
      enddo

      call bpinitialize()
      call blockmatrix(bmat, n, a, ja, ia, 10, kvst, BP_SPARSE)
      call preconditioner(precon, bmat, BP_BJACOBI, 0.d0, 0.d0, 
     *                                  BP_INVERSE, 0.d0, 0.d0)
      call flexgmres(bmat, sol, rhs, precon, 20, 600, 1.d-8)
      call freepreconditioner(precon)
      call freeblockmatrix(bmat)

      end
