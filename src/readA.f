      program readA
      implicit none

      integer lmaxh,ldummy
      integer mima,nsim,numatd
      integer lenyh,npath,i
      real*4 xlat1,xlon1,xlat2,xlon2
      real*4 wt
      complex*8 yt(15395) ! each element is a path integral of a spherical harmonic (lmax=20)
      
      
      npath=0
      write(*,*) lenyh,npath
      
      open(2,file='outputs/testA'
     1 ,action='read',form='formatted',status='old')
      write(*,*) 'file opened'

      read(2,*) lmaxh,ldummy
      lenyh=(lmaxh+1)**2
      write(*,*) lmaxh,ldummy,lenyh

   10 continue 
      write(*,*) npath 
      read(2,*,end=199) xlat1,xlon1,xlat2,xlon2,mima,wt,nsim
     1 ,(yt(i),i=1,lenyh)
      npath = npath+1
      write(*,*) npath
      goto 10
   99 continue

  199 write(*,*) xlat1,xlon1,xlat2,xlon2,mima,wt,nsim,yt
      write(*,*) 'DONE'
      write(*,*) npath
      end program