      program test
      implicit real*8(a-h,m,o-z)
      data e0/78408.4D00/,r0/1.292D00/,p0/1.0133D05/,
     &     t0/273.15D00/,gascon/287.06D00/,gam0/1.40D00/
      a0=sqrt(gam0*gascon*t0)
      open(unit=1,file='tgas_test.dat',status='unknown')
c
      write(1,'(''zone  '')')
      den=r0*1.00D-07
      dx=(3.50D00-0.50D00)/99.00D00
      x0=0.50D00
      do 100 i=1,100
         x=(i-1)*dx+x0
         eint=(10.00D00**x)*gascon*t0
         call tgas1(eint,den,pre,ssp,temp,3,ier)
         call tgas3(pre,den,temp3,ier)
         call tgas4(pre,den,h,ier)
         call tgas2(eint,den,s,ier)
         call tgas5(pre,temp,den5,eint5,ssp5,ier)
         write(6,'(i5,1pd13.5,2d13.5,i5)') i,den,den5,
     &         den-den5,ier
         eratio=log10(eint/(gascon*t0))
         pratio=log10(pre/p0)
         tratio=log10(temp/t0)
         aratio=log10(ssp/a0)
         t3ratio=log10(temp3/t0)
         hratio=log10(h/(gascon*t0))
         sratio=log10(s/gascon)
         write(1,'(1pd13.5,6d13.5)') eratio,pratio,
     &         tratio,aratio,t3ratio,hratio,sratio
  100 continue
c
      write(1,'(''zone  '')')
      den=r0*1.00D-06
      dx=(3.50D00-0.50D00)/99.00D00
      x0=0.50D00
      do 200 i=1,100
          x=(i-1)*dx+x0
         eint=(10.00**x)*gascon*t0
         call tgas1(eint,den,pre,ssp,temp,3,ier)
         call tgas3(pre,den,temp3,ier)
         call tgas4(pre,den,h,ier)
         call tgas2(eint,den,s,ier)
         call tgas5(pre,temp,den5,eint5,ssp5,ier)
         write(6,'(i5,1pd13.5,2d13.5,i5)') i,den,den5,
     &         den-den5,ier
         eratio=log10(eint/(gascon*t0))
         pratio=log10(pre/p0)
         tratio=log10(temp/t0)
         aratio=log10(ssp/a0)
         t3ratio=log10(temp3/t0)
         hratio=log10(h/(gascon*t0))
         sratio=log10(s/gascon)
         write(1,'(1pd13.5,6d13.5)') eratio,pratio,
     &         tratio,aratio,t3ratio,hratio,sratio
  200 continue
c
      write(1,'(''zone  '')')
      den=r0*1.00D-05
      dx=(3.25D00-0.50D00)/99.00D00
      x0=0.50D00
      do 300 i=1,100
         x=(i-1)*dx+x0
         eint=(10.00D00**x)*gascon*t0
         call tgas1(eint,den,pre,ssp,temp,3,ier)
         call tgas3(pre,den,temp3,ier)
         call tgas4(pre,den,h,ier)
         call tgas2(eint,den,s,ier)
         call tgas5(pre,temp,den5,eint5,ssp5,ier)
         write(6,'(i5,1pd13.5,2d13.5,i5)') i,den,den5,
     &         den-den5,ier
         eratio=log10(eint/(gascon*t0))
         pratio=log10(pre/p0)
         tratio=log10(temp/t0)
         aratio=log10(ssp/a0)
         t3ratio=log10(temp3/t0)
         hratio=log10(h/(gascon*t0))
         sratio=log10(s/gascon)
         write(1,'(1pd13.5,6d13.5)') eratio,pratio,
     &         tratio,aratio,t3ratio,hratio,sratio
  300 continue
c
      write(1,'(''zone  '')')
      den=r0*1.00D-04
      dx=(3.25D00-0.50D00)/99.00D00
      x0=0.50D00
      do 400 i=1,100
         x=(i-1)*dx+x0
         eint=(10.00D00**x)*gascon*t0
         call tgas1(eint,den,pre,ssp,temp,3,ier)
         call tgas3(pre,den,temp3,ier)
         call tgas4(pre,den,h,ier)
         call tgas2(eint,den,s,ier)
         call tgas5(pre,temp,den5,eint5,ssp5,ier)
         write(6,'(i5,1pd13.5,2d13.5,i5)') i,den,den5,
     &         den-den5,ier
         eratio=log10(eint/(gascon*t0))
         pratio=log10(pre/p0)
         tratio=log10(temp/t0)
         aratio=log10(ssp/a0)
         t3ratio=log10(temp3/t0)
         hratio=log10(h/(gascon*t0))
         sratio=log10(s/gascon)
         write(1,'(1pd13.5,6d13.5)') eratio,pratio,
     &         tratio,aratio,t3ratio,hratio,sratio
  400 continue
c
      write(1,'(''zone  '')')
      den=r0*1.00D-03
      dx=(3.25-0.50)/99.0
      x0=0.50
      do 500 i=1,100
         x=(i-1)*dx+x0
         eint=(10.00**x)*gascon*t0
         call tgas1(eint,den,pre,ssp,temp,3,ier)
         call tgas3(pre,den,temp3,ier)
         call tgas4(pre,den,h,ier)
         call tgas2(eint,den,s,ier)
         call tgas5(pre,temp,den5,eint5,ssp5,ier)
         write(6,'(i5,1pd13.5,2d13.5,i5)') i,den,den5,
     &         den-den5,ier
         eratio=log10(eint/(gascon*t0))
         pratio=log10(pre/p0)
         tratio=log10(temp/t0)
         aratio=log10(ssp/a0)
         t3ratio=log10(temp3/t0)
         hratio=log10(h/(gascon*t0))
         sratio=log10(s/gascon)
         write(1,'(1pd13.5,6d13.5)') eratio,pratio,
     &         tratio,aratio,t3ratio,hratio,sratio
  500 continue
c
      write(1,'(''zone  '')')
      den=r0*1.00D-02
      dx=(3.25D00-0.50D00)/99.00D00
      x0=0.50
      do 600 i=1,100
         x=(i-1)*dx+x0
         eint=(10.00D00**x)*gascon*t0
         call tgas1(eint,den,pre,ssp,temp,3,ier)
         call tgas3(pre,den,temp3,ier)
         call tgas4(pre,den,h,ier)
         call tgas2(eint,den,s,ier)
         call tgas5(pre,temp,den5,eint5,ssp5,ier)
         write(6,'(i5,1pd13.5,2d13.5,i5)') i,den,den5,
     &         den-den5,ier
         eratio=log10(eint/(gascon*t0))
         pratio=log10(pre/p0)
         tratio=log10(temp/t0)
         aratio=log10(ssp/a0)
         t3ratio=log10(temp3/t0)
         hratio=log10(h/(gascon*t0))
         sratio=log10(s/gascon)
         write(1,'(1pd13.5,6d13.5)') eratio,pratio,
     &         tratio,aratio,t3ratio,hratio,sratio
  600 continue
c
      write(1,'(''zone  '')')
      den=r0*1.00D-01
      dx=(3.00D00-0.50D00)/99.00D00
      x0=0.50
      do 700 i=1,100
         x=(i-1)*dx+x0
         eint=(10.00D00**x)*gascon*t0
         call tgas1(eint,den,pre,ssp,temp,3,ier)
         call tgas3(pre,den,temp3,ier)
         call tgas4(pre,den,h,ier)
         call tgas2(eint,den,s,ier)
         call tgas5(pre,temp,den5,eint5,ssp5,ier)
         write(6,'(i5,1pd13.5,2d13.5,i5)') i,den,den5,
     &         den-den5,ier
         eratio=log10(eint/(gascon*t0))
         pratio=log10(pre/p0)
         tratio=log10(temp/t0)
         aratio=log10(ssp/a0)
         t3ratio=log10(temp3/t0)
         hratio=log10(h/(gascon*t0))
         sratio=log10(s/gascon)
         write(1,'(1pd13.5,6d13.5)') eratio,pratio,
     &         tratio,aratio,t3ratio,hratio,sratio
  700 continue
c
      write(1,'(''zone  '')')
      den=r0*1.00D00
      dx=(3.00D00-0.50D00)/99.00D00
      x0=0.50
      do 800 i=1,100
         x=(i-1)*dx+x0
         eint=(10.00D00**x)*gascon*t0
         call tgas1(eint,den,pre,ssp,temp,3,ier)
         call tgas3(pre,den,temp3,ier)
         call tgas4(pre,den,h,ier)
         call tgas2(eint,den,s,ier)
         call tgas5(pre,temp,den5,eint5,ssp5,ier)
         write(6,'(i5,1pd13.5,2d13.5,i5)') i,den,den5,
     &         den-den5,ier
         eratio=log10(eint/(gascon*t0))
         pratio=log10(pre/p0)
         tratio=log10(temp/t0)
         aratio=log10(ssp/a0)
         t3ratio=log10(temp3/t0)
         hratio=log10(h/(gascon*t0))
         sratio=log10(s/gascon)
         write(1,'(1pd13.5,6d13.5)') eratio,pratio,
     &         tratio,aratio,t3ratio,hratio,sratio
  800 continue
c
      write(1,'(''zone  '')')
      den=r0*1.00D01
      dx=(3.00D00-0.50D00)/99.00D00
      x0=0.50
      do 900 i=1,100
         x=(i-1)*dx+x0
         eint=(10.00D00**x)*gascon*t0
         call tgas1(eint,den,pre,ssp,temp,3,ier)
         call tgas3(pre,den,temp3,ier)
         call tgas4(pre,den,h,ier)
         call tgas2(eint,den,s,ier)
         call tgas5(pre,temp,den5,eint5,ssp5,ier)
         write(6,'(i5,1pd13.5,2d13.5,i5)') i,den,den5,
     &         den-den5,ier
         eratio=log10(eint/(gascon*t0))
         pratio=log10(pre/p0)
         tratio=log10(temp/t0)
         aratio=log10(ssp/a0)
         t3ratio=log10(temp3/t0)
         hratio=log10(h/(gascon*t0))
         sratio=log10(s/gascon)
         write(1,'(1pd13.5,6d13.5)') eratio,pratio,
     &         tratio,aratio,t3ratio,hratio,sratio
  900 continue
c
      write(1,'(''zone  '')')
      den=r0*1.00D02
      dx=(2.80D00-0.50D00)/99.00D00
      x0=0.50
      do 1000 i=1,100
         x=(i-1)*dx+x0
         eint=(10.00D00**x)*gascon*t0
         call tgas1(eint,den,pre,ssp,temp,3,ier)
         call tgas3(pre,den,temp3,ier)
         call tgas4(pre,den,h,ier)
         call tgas2(eint,den,s,ier)
         call tgas5(pre,temp,den5,eint5,ssp5,ier)
         write(6,'(i5,1pd13.5,2d13.5,i5)') i,den,den5,
     &         den-den5,ier
         eratio=log10(eint/(gascon*t0))
         pratio=log10(pre/p0)
         tratio=log10(temp/t0)
         aratio=log10(ssp/a0)
         t3ratio=log10(temp3/t0)
         hratio=log10(h/(gascon*t0))
         sratio=log10(s/gascon)
         write(1,'(1pd13.5,6d13.5)') eratio,pratio,
     &         tratio,aratio,t3ratio,hratio,sratio
 1000 continue
c
      write(1,'(''zone  '')')
      den=r0*1.00D03
      dx=(2.80D00-0.50D00)/99.00D00
      x0=0.50
      do 1100 i=1,100
         x=(i-1)*dx+x0
         eint=(10.00D00**x)*gascon*t0
         call tgas1(eint,den,pre,ssp,temp,3,ier)
         call tgas3(pre,den,temp3,ier)
         call tgas4(pre,den,h,ier)
         call tgas2(eint,den,s,ier)
         call tgas5(pre,temp,den5,eint5,ssp5,ier)
         write(6,'(i5,1pd13.5,2d13.5,i5)') i,den,den5,
     &         den-den5,ier
         eratio=log10(eint/(gascon*t0))
         pratio=log10(pre/p0)
         tratio=log10(temp/t0)
         aratio=log10(ssp/a0)
         t3ratio=log10(temp3/t0)
         hratio=log10(h/(gascon*t0))
         sratio=log10(s/gascon)
         write(1,'(1pd13.5,6d13.5)') eratio,pratio,
     &         tratio,aratio,t3ratio,hratio,sratio
 1100 continue
c
      close(unit=1)
      stop
      end
