      subroutine minos_start(
     &  model_path,station,str1,str2,
     &  motion,motion_number,modelin,str3,str4,str5,str6) 

      character*256 filnam
      character*200  model_path
      character*256  path,str3,str4
      integer*4      iproc
      character*6    station
      character*20   str5,str6   
      character*20   str2 
      character*3    motion
      character*1    motion_number
      character*100  modelin 
      character*2    str1


c      write(*,*) "model_path ",model_path
c      write(*,*) "station ",station
c      write(*,*) "str1 ", str1
c      write(*,*) "str2 ",str2
c      write(*,*) "motion ",motion
c      write(*,*) "motion_number ",motion_number
c      write(*,*) "modelin ",modelin
c      write(*,*) "str3 ",str3
c      write(*,*) "str4 ",str4
c      write(*,*) "str5 ",str5
c      write(*,*) "str6 ",str6
c
c      print *,'input model file:'
      read(model_path,*)
c  100 format(a256)
c MB added one line below
c      print *,model_path
c(1:lnblnk(station))
      open(7,file=model_path,status='old',form='formatted',iostat=iret)
c      print *,'output file:'
      read(str3,*)
c      print*, str3
cc MB added one line below
      open(8,file=str3,form='formatted',iostat=iret)
c      call model(iin,iout)
      call model(7,8) 
      close(7)
c      print *,'eigenfunction file (output):'
      read(str4,*)
cc MB added one line below
c      print *,str4
      ifreq=1
      if(filnam(1:4).eq.'none') ifreq=0
      open(3,file=str4,form='unformatted',iostat=iret)
      call wtable(8,3,ifreq,str5,motion_number,str6)
c      write(*,*) "----- ", ifreq,str5,motion_number," ",str6
c      write(*,*) "str5 ", str5
c      write(*,*) "motion_number ", motion_number
c      write(*,*) "str6 ", str6
      close(8)  
      close(3)
      end subroutine 
