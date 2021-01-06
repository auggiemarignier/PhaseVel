      program kernels
    
      use kernelfunctions

      implicit real*8(a-h,o-z)

      integer*4 mk
      parameter (mk=3000)

      character*256  model_file,out_plain_file,out_bin_file,dbase_name,
     1 eigenasc 
      real*4      rad(mk)
      real*4      U(mk),Up(mk),V(mk),Vp(mk),P(mk),Pp(mk),W(mk),Wp(mk)
      real*4      omega,wavenum
      real*4      kkappa(mk),kmu(mk)
      real*8      alpha(mk),beta(mk)

      real*4    per_eigen,phvel_eigen,grvel_eigen,attn_eigen
      integer*4 norder_eigen,lorder_eigen,eigid_eigen,
     +          nraw_eigen,ncol_eigen,npar_eigen,foff_eigen,
     +          commid_eigen
      character*2 datatype_eigen
      character*64 dir_eigen
      character*32 dfile_eigen
      character*17 lddate_eigen
      character*1 typeo_eigen
      common/c_eigen/norder_eigen,lorder_eigen,
     +      eigid_eigen,per_eigen,phvel_eigen,grvel_eigen,
     +      attn_eigen,nraw_eigen,ncol_eigen,npar_eigen,
     +      foff_eigen,commid_eigen,typeo_eigen,
     +      datatype_eigen,dir_eigen,dfile_eigen,lddate_eigen

      common r(mk),fmu(mk),flam(mk),qshear(mk),qkappa(mk),
     + xa2(mk),xlam(mk),rho(mk),qro(3,mk),g(mk),qg(3,mk),
     + fcon(mk),fspl(3,mk),lcon(mk),lspl(3,mk),ncon(mk),
     + nspl(3,mk),ccon(mk),cspl(3,mk),acon(mk),aspl(3,mk)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,jcom,nord,l,kg,kount,knsw,ifanis,iback
      common/eifx/vpv(mk),vph(mk),vsv(mk),vsh(mk),eta(mk),wrk(mk*10)

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
      alpha=vpv
      beta=vsv
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

      do n=0,nmax
        do l=lmax,lmax
          call read_nleigenfucntion(n,l,dbase_name,
     1                              rad,U,Up,V,Vp,P,Pp,W,Wp)
          omega = 2*pi/per_eigen
          wavenum = omega/phvel_eigen
          kkappa = kernel_kappa(omega,wavenum,rad,U,Up,V)
          kmu = kernel_mu(omega,wavenum,rad,U,Up,V,Vp,W,Wp)
          call write_eigenfunctions_asc(eigenasc,
     1                                  rad,U,Up,V,Vp,P,Pp,W,Wp)
        enddo
      enddo
      end program