      program test
      implicit real*8(a-h,m,o-z)
      real*8 k,k0,kratio
      data e0/78408.4D00/,r0/1.243D00/,
     &     t0/273.15D00/,k0/1.87915D-02/
      mu0=1.058D-06*16.5273D00
      open(unit=1,file='ugas_test.dat',status='unknown')
c
      write(1,'(''zone  '')')
      den=r0*1.00D-05
      dx1=(3.50D00-0.50D00)/99.00D00
      dx2=(15.00D00-0.50D00)/99.00D00
      x10=0.50D00
      x20=0.50D00
      do 100 i=1,100
         x1=(i-1)*dx1+x10
         x2=(i-1)*dx2+x20
         eint=(10.00D00**x1)*e0
         temp=1000.00D00*x2
         call ugas1(temp,den,mu,ierror)
         if (ierror.ne.0) write(6,*) 'ugas1: t=',temp,' d=',den
         call ugas2(temp,den,pr,ierror)
         if (ierror.ne.0) write(6,*) 'ugas2: t=',temp,' d=',den
         call ugas4(eint,den,k,ierror)
         if (ierror.ne.0) write(6,*) 'ugas4: e=',eint,' d=',den
         eratio=log10(eint/e0)
         tratio=temp/1000.00D00
         mratio=mu/mu0
         kratio=k/k0
         write(1,'(1pd13.5,4d13.5)') tratio,eratio,mratio,pr,kratio
  100 continue
c
      write(1,'(''zone  '')')
      den=r0*1.00D-04
      dx1=(3.50D00-0.50D00)/99.00D00
      dx2=(15.00D00-0.50D00)/99.00D00
      x10=0.50D00
      x20=0.50D00
      do 200 i=1,100
         x1=(i-1)*dx1+x10
         x2=(i-1)*dx2+x20
         eint=(10.00D00**x1)*e0
         temp=1000.00D00*x2
         call ugas1(temp,den,mu,ierror)
         if (ierror.ne.0) write(6,*) 'ugas1: t=',temp,' d=',den
         call ugas2(temp,den,pr,ierror)
         if (ierror.ne.0) write(6,*) 'ugas2: t=',temp,' d=',den
         call ugas4(eint,den,k,ierror)
         if (ierror.ne.0) write(6,*) 'ugas4: e=',eint,' d=',den
         eratio=log10(eint/e0)
         tratio=temp/1000.00D00
         mratio=mu/mu0
         kratio=k/k0
         write(1,'(1pd13.5,4d13.5)') tratio,eratio,mratio,pr,kratio
  200 continue
c
      write(1,'(''zone  '')')
      den=r0*1.00D-03
      dx1=(3.50D00-0.50D00)/99.00D00
      dx2=(15.00D00-0.50D00)/99.00D00
      x10=0.50D00
      x20=0.50D00
      do 300 i=1,100
         x1=(i-1)*dx1+x10
         x2=(i-1)*dx2+x20
         eint=(10.00D00**x1)*e0
         temp=1000.00D00*x2
         call ugas1(temp,den,mu,ierror)
         if (ierror.ne.0) write(6,*) 'ugas1: t=',temp,' d=',den
         call ugas2(temp,den,pr,ierror)
         if (ierror.ne.0) write(6,*) 'ugas2: t=',temp,' d=',den
         call ugas4(eint,den,k,ierror)
         if (ierror.ne.0) write(6,*) 'ugas4: e=',eint,' d=',den
         eratio=log10(eint/e0)
         tratio=temp/1000.00D00
         mratio=mu/mu0
         kratio=k/k0
         write(1,'(1pd13.5,4d13.5)') tratio,eratio,mratio,pr,kratio
  300 continue
c
      write(1,'(''zone  '')')
      den=r0*1.00D-02
      dx1=(3.50D00-0.50D00)/99.00D00
      dx2=(15.00D00-0.50D00)/99.00D00
      x10=0.50D00
      x20=0.50D00
      do 400 i=1,100
         x1=(i-1)*dx1+x10
         x2=(i-1)*dx2+x20
         eint=(10.00D00**x1)*e0
         temp=1000.00D00*x2
         call ugas1(temp,den,mu,ierror)
         if (ierror.ne.0) write(6,*) 'ugas1: t=',temp,' d=',den
         call ugas2(temp,den,pr,ierror)
         if (ierror.ne.0) write(6,*) 'ugas2: t=',temp,' d=',den
         call ugas4(eint,den,k,ierror)
         if (ierror.ne.0) write(6,*) 'ugas4: e=',eint,' d=',den
         eratio=log10(eint/e0)
         tratio=temp/1000.00D00
         mratio=mu/mu0
         kratio=k/k0
         write(1,'(1pd13.5,4d13.5)') tratio,eratio,mratio,pr,kratio
  400 continue
c
      write(1,'(''zone  '')')
      den=r0*1.00D-01
      dx1=(3.50D00-0.50D00)/99.00D00
      dx2=(15.00D00-0.50D00)/99.00D00
      x10=0.50D00
      x20=0.50D00
      do 500 i=1,100
         x1=(i-1)*dx1+x10
         x2=(i-1)*dx2+x20
         eint=(10.00D00**x1)*e0
         temp=1000.00D00*x2
         call ugas1(temp,den,mu,ierror)
         if (ierror.ne.0) write(6,*) 'ugas1: t=',temp,' d=',den
         call ugas2(temp,den,pr,ierror)
         if (ierror.ne.0) write(6,*) 'ugas2: t=',temp,' d=',den
         call ugas4(eint,den,k,ierror)
         if (ierror.ne.0) write(6,*) 'ugas4: e=',eint,' d=',den
         eratio=log10(eint/e0)
         tratio=temp/1000.00D00
         mratio=mu/mu0
         kratio=k/k0
         write(1,'(1pd13.5,4d13.5)') tratio,eratio,mratio,pr,kratio
  500 continue
c
      write(1,'(''zone  '')')
      den=r0*1.00D00
      dx1=(3.50D00-0.50D00)/99.00D00
      dx2=(15.00-0.50)/99.00D00
      x10=0.50D00
      x20=0.50D00
      do 600 i=1,100
         x1=(i-1)*dx1+x10
         x2=(i-1)*dx2+x20
         eint=(10.00D00**x1)*e0
         temp=1000.00D00*x2
         call ugas1(temp,den,mu,ierror)
         if (ierror.ne.0) write(6,*) 'ugas1: t=',temp,' d=',den
         call ugas2(temp,den,pr,ierror)
         if (ierror.ne.0) write(6,*) 'ugas2: t=',temp,' d=',den
         call ugas4(eint,den,k,ierror)
         if (ierror.ne.0) write(6,*) 'ugas4: e=',eint,' d=',den
         eratio=log10(eint/e0)
         tratio=temp/1000.00D00
         mratio=mu/mu0
         kratio=k/k0
         write(1,'(1pd13.5,4d13.5)') tratio,eratio,mratio,pr,kratio
  600 continue
c
      write(1,'(''zone  '')')
      den=r0*1.00D01
      dx1=(3.50D00-0.50D0)/99.00D00
      dx2=(15.00D00-0.50D00)/99.00D00
      x10=0.50D00
      x20=0.50D00
      do 700 i=1,100
         x1=(i-1)*dx1+x10
         x2=(i-1)*dx2+x20
         eint=(10.00D00**x1)*e0
         temp=1000.00D00*x2
         call ugas1(temp,den,mu,ierror)
         if (ierror.ne.0) write(6,*) 'ugas1: t=',temp,' d=',den
         call ugas2(temp,den,pr,ierror)
         if (ierror.ne.0) write(6,*) 'ugas2: t=',temp,' d=',den
         call ugas4(eint,den,k,ierror)
         if (ierror.ne.0) write(6,*) 'ugas4: e=',eint,' d=',den
         eratio=log10(eint/e0)
         tratio=temp/1000.00D00
         mratio=mu/mu0
         kratio=k/k0
         write(1,'(1pd13.5,4d13.5)') tratio,eratio,mratio,pr,kratio
  700 continue
c
      close(unit=1)
      stop
      end
