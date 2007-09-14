      program f3example
      implicit none
c-----------------------------------------------------------------------
c     Example program that calls the BPKIT library from FORTRAN.
c
c     Tests the subroutine "blockmatrix2" for converting a matrix in 
c     block coordinate format to the BPKIT BlockMat format.  
c
c     The following matrix is tested.  The data for this matrix is 
c     given below in block coordinate format.  The calling sequence 
c     for the "blockmatrix2" function is also shown.
c
c     [  4  -1 | -1   0 |        ]
c     [ -1   4 |  0  -1 |        ]
c     [--------+--------+--------]
c     [ -1   0 |  4  -1 | -1   0 ]
c     [  0  -1 | -1   4 |  0  -1 ]
c     [--------+--------+--------]
c     [        | -1   0 |  4  -1 ]
c     [        |  0  -1 | -1   4 ]
c
c     This matrix has the following properties:
c     nrows = 3
c     nnz = 7
c     blocksize = 2
c
c-----------------------------------------------------------------------
      include  'bpfort.h'

      integer   nrows, nnz, blocksize
      real*8    rhs(6), sol(6)
      integer*8 bmat, precon

      integer   row(7), col(7)
      real*8    A(28)

      data row /1,1,2,2,2,3,3/
      data col /1,2,1,2,3,2,3/
      data A   /4.0,-1.0,-1.0,4.0,
     *          -1.0,0.0,0.0,-1.0,
     *          -1.0,0.0,0.0,-1.0,
     *          4.0,-1.0,-1.0,4.0,
     *          -1.0,0.0,0.0,-1.0,
     *          -1.0,0.0,0.0,-1.0,
     *          4.0,-1.0,-1.0,4.0/

      data rhs /1.0,1.0,1.0,1.0,1.0,1.0/
      data sol /0.0,0.0,0.0,0.0,0.0,0.0/

      nrows = 3
      nnz = 7
      blocksize = 2
      
      call bpinitialize()

      call blockmatrix2(bmat, nrows, nnz, row, col, A, blocksize)

      call preconditioner(precon, bmat, BP_BILUK,   0.d0, 0.d0, 
     *                                  BP_INVERSE, 0.d0, 0.d0)
      call flexgmres(bmat, sol, rhs, precon, 20, 600, 1.d-8)
      call freepreconditioner(precon)
      call freeblockmatrix(bmat)

      end
