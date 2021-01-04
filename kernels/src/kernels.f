      program kernels

      implicit real*8(a-h,o-z)

      character*256  model_file,out_plain_file,out_bin_file
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,jcom,nord,l,kg,kount,knsw,ifanis,iback

      jcom=3
      eps=1e-7
      wgrav=10

      lmin=0
      lmax=300
      wmin=0
      wmax=50
      nmin=0
      nmax=0

      write(*,*) "building kernels!"
      model_file="/usr/local/opt/mineos/DEMO"
     1 //"/models/prem_noocean.txt"
      out_plain_file="/Users/auggiemarignier/Documents/PhD/PhaseVel/"
     1 //"kernels/outputs/properties.txt"
      out_bin_file="/Users/auggiemarignier/Documents/PhD/PhaseVel/"
     1 //"kernels/outputs/eigenfunctions"

      write(*,*) model_file
      open(7,file=model_file,status='old',form='formatted',iostat=iret)
      open(8,file=out_plain_file,form='formatted',iostat=iret)
      call model(7,8) 
      close(7)
      ifreq=1
      open(3,file=out_bin_file,form='unformatted',iostat=iret)
      call wtable(8,3,ifreq,lmin,lmax,wmin,wmax,nmin,nmax)
      close(8)  
      close(3)


      end program