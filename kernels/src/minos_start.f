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


      read(model_path,*)
      open(7,file=model_path,status='old',form='formatted',iostat=iret)
      read(str3,*)
      open(8,file=str3,form='formatted',iostat=iret)
      call model(7,8) 
      close(7)
      read(str4,*)
      ifreq=1
      if(filnam(1:4).eq.'none') ifreq=0
      open(3,file=str4,form='unformatted',iostat=iret)
      call wtable(8,3,ifreq,str5,motion_number,str6)
      close(8)  
      close(3)
      end subroutine 
