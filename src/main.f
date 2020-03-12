      program main
      
      include 'parameters.h'
      character*100 infile

      real*8, dimension(MXPATHS,MXLENY) :: avylm,clm
      real*8, dimension(MXPATHS) :: pred


      infile = '../inputs/hvh.000S061.asc.rwt3_std4'
      call ylmpav(infile,avylm)
      write(*,*) "File Read. Back in main"
c     call choose clm via S2LET
      do i=1,MXPATHS
        do j=1,MXLENY
            clm(i,j) = 1
        enddo
      enddo
      call forward_modelling(clm,avylm,pred)
      end program main