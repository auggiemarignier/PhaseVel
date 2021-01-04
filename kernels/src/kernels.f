      program kernels

      implicit real*8(a-h,o-z)

      character*200  model_path
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,jcom,nord,l,kg,kount,knsw,ifanis,iback

      jcom=3
      eps=1e-7
      wgrav=10

      write(*,*) "building kernels!"
      model_path="/usr/local/opt/mineos/DEMO"
     1 //"/models/prem_noocean.txt"

      write(*,*) model_path
      call minos_start(model_path,lmin,lmax,wmin,wmax,nmin,nmax)


      end program