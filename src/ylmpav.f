c     Program for calculating Path-averages of spherical harmonics between two points
      subroutine ylmpav(infile,avylm)
      
      include 'parameters.h'

      real*4 xlat1,xlat2,xlon1,xlon2,wt
      real*8 dc
      integer mima,nsim,npath
      character*100 infile


      dimension ylmh(MXLENY),ylmt(MXLENY)
      double precision d((2*LMAX+1)**2)
      dimension sar(LMAX+1),wk1(LMAX+1),wk2(LMAX+1),wk3(LMAX+1)
      real*8, dimension(MXPATHS,MXLENY) :: avylm

     
      lenyh = (LMAX+1)**2
      npath=0

c     Read data file
      open(1,file=infile)
      open(2,file='../outputs/000S061YlmavL35',form='formatted')
   10 read(1,*,end=99) xlat1,xlon1,xlat2,xlon2,dc,wt,mima,nsim
      npath = npath+1
    
c     Calculate path average of Ylms
      call ylmav(xlat1,xlon1,xlat2,xlon2,LMAX,ylmh,ylmt
     1   ,wk1,wk2,wk3,sar,d)
      if(mima.eq.1) then
        write(2,*) xlat1,xlon1,xlat2,xlon2,mima,wt,nsim
     1     ,(ylmt(i),i=1,lenyh)
    !    avylm(npath,1:lenyh) = ylmt(1:lenyh)
    !   else
    !    ch=360./(360.-del)
    !    ct=-del/(360.-del)
    !    write(2) xlat1,xlon1,xlat2,xlon2,mima,wt,nsim
    !  1     ,(ylmh(i)*ch+ylmt(i)*ct,i=1,lenyh)

      endif
      
      
      if(mod(npath,1000).eq.0) write(6,'(i8,'' paths read'')') npath
      goto 10
   99 continue
      write(6,'(i8,'' paths read'')') npath

      end


