      program kernels

      implicit real*8(a-h,o-z)

      integer*4 mk
      parameter (mk=3000)

      character*256  model_file,out_plain_file,out_bin_file,dbase_name,
     1 eigenasc 
      real*4      rout(mk)
      real*4      U(mk),Up(mk),V(mk),Vp(mk),P(mk),Pp(mk),W(mk),Wp(mk)

      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,jcom,nord,l,kg,kount,knsw,ifanis,iback

      jcom=3
      eps=1e-7
      wgrav=10

      lmin=0
      lmax=300
      wmin=0
      wmax=100
      nmin=0
      nmax=0

      write(*,*) "building kernels!"
      model_file="/usr/local/opt/mineos/DEMO"
     1 //"/models/prem_noocean.txt"
      out_plain_file="/Users/auggiemarignier/Documents/PhD/PhaseVel/"
     1 //"kernels/outputscheck/properties.txt"
      out_bin_file="/Users/auggiemarignier/Documents/PhD/PhaseVel/"
     1 //"kernels/outputscheck/eigenfunctions"

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

      dbase_name="/Users/auggiemarignier/Documents/PhD/PhaseVel/"
     1 //"kernels/outputscheck/database"
      call eigcon(jcom,model_file,out_plain_file,out_bin_file,
     1 dbase_name,6371.0)

      eigenasc=trim(out_bin_file)//"asc"
    !   call eigen2asc(nmin,nmax,lmin,lmax,dbase_name,eigenasc)

      do n=0,nmax
        do l=0,lmax
          call read_nleigenfucntion(n,l,dbase_name,
     1                              rout,U,Up,V,Vp,P,Pp,W,Wp)
          call write_eigenfunctions_asc(eigenasc,
     1                                  rout,U,Up,V,Vp,P,Pp,W,Wp)
        enddo
      enddo
      end program