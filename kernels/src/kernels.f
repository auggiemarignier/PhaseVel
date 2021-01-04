      program kernels

      implicit real*8(a-h,o-z)

      character*200  model_file,out_plain_file,out_bin_file
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
      imod=7
      out_plain_file="/Users/auggiemarignier/Documents/PhD/PhaseVel/"
     1 //"kernels/outputs/properties.txt"
      iplain=8
      out_bin_file="/Users/auggiemarignier/Documents/PhD/PhaseVel/"
     1 //"kernels/outputs/eigenfunctions"
      ibin=3

      write(*,*) model_file
      open(imod,file=model_file,status='old',form='formatted',
     1 iostat=iret)
      open(iplain,file=out_plain_file,form='formatted',iostat=iret)
      call model(imod,iplain) 
      close(imod)
      ifreq=1
      open(ibin,file=out_bin_file,form='unformatted',iostat=iret)
      call wtable(iplain,ibin,ifreq,lmin,lmax,wmin,wmax,nmin,nmax)
      close(iplain)  
      close(ibin)


      end program