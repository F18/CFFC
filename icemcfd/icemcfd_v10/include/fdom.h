C
C fortran version of domain.h
C
	integer dferr
	parameter (dferr = -1)
C
C modes for opening domain file
C
	integer dfread, dfwrit, dfmod
	parameter (dfread = 0)
	parameter (dfwrit = 1)
	parameter (dfmod = 3)
C
C domain types
C

	integer dfundef, dfglob, dfstrt, dfunst
	parameter (dfundef = -1)
	parameter (dfglob = 0)
	parameter (dfstrt = 1)
	parameter (dfunst = 2)
C
C coodinate systems
C
	integer dfcart, dfcyl, dfsph
	parameter (dfcart = 0)
	parameter (dfcyl = 1)
	parameter (dfsph = 2)
C
C type of list representing domains
C
	integer dfrang, dfexpl
	parameter (dfrang = 0)
	parameter (dfexpl = 1)
