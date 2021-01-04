      subroutine eigcon_subs(imod,
     &dmax,modelin,str2,motion,str3,str4,str8)
      implicit none
      integer*4 mk
      parameter (mk=350)
      include "fdb_eigen.h"
c --- other variables
      real*8    r(mk),pi
      real*4    param(7),buf(7,mk)
      real*4    bigg,tau,rhobar,con,accn
      real*4    tref,rx,r1,r2,dr,dum,rn,gn,vn,wn
      real*4    dmax,rmax,ww,qq,qc,tf,fl,fl1,fl3
      integer*4 i,ieig,idat,iin,ifanis,ifdeck,k,j,jj,jcom,l,ll,lll
      integer*4 n,nic,noc,nreg,nn,nlay,nstart,nrad,ni,nrecl
      integer*4 lnblnk,ierr,iflag,ititle(20)
      character*20 str
      character*64 dir
      character*1 typeo
      character*256 fmodel,fflatin,fbinin,fout,fin
      character*200 model_path
      character*100 modelin
      character*20  str2
      character*256 str3,str4,str8
      character*3    motion
c--- my added declarations below
      integer*4 q
      integer*4   in,nmine,nmaxe,lmine,lmaxe
      real*4      rout(mk)
      real*8      Amp
      real*8      W(299),S(299)
      logical     ttf
      character*256 fdir,cmd
c --
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,jcom,nord,l,kg,kount,knsw,ifanis,iback
      common/prm/ nn,ll,ww,qq,rn,vn,accn
      equivalence (param(1),nn)
c ---
      data bigg,tau,rhobar/6.6723e-11,1000.0,5515.0/
      pi=4.0d0*datan(1.d0)
      con=pi*bigg
      nstart = 0

c read model file name
      read(model_path,'(a256)')fmodel
c read minos_bran output text file
      read(str3,'(a256)') fflatin
c read minos_bran output binary unformatted file
      read(str4,'(a256)') fbinin
      read(str8,'(a256)')fout
c
c  read in radial knots from model
c
      iin = 7
      open(iin,file=fmodel,status='old')
      read(iin,101) (ititle(i),i=1,20)
  101 format(20a4)
      read(iin,*) ifanis,tref,ifdeck
      if(ifdeck.eq.0) go to 1000
c*** card deck model ***
      read(iin,*) n,nic,noc
      read(iin,105) (r(i),i=1,n)
  105 format(f8.0)
      go to 2000
c*** polynomial model ***
 1000 read(iin,*) nreg,nic,noc,rx
      n=0
      jj=5
      if(ifanis.ne.0) jj=8
      do 10 nn=1,nreg
      read(iin,*) nlay,r1,r2
      r1=r1*tau
      r2=r2*tau
      dr=(r2-r1)/float(nlay-1)
      do 15 i=1,nlay
      n=n+1
   15 r(n)=r1+dr*float(i-1)
      do 10 j=1,jj
   10    read(iin,110) dum
  110 format(f9.5)
 2000 close(iin)
c
c  rn     : radius at surface
c  wn     : frequency normalization
c  vn     : velocity normalisation
c  accn   : acceleration normalisation
c  n      : index of surface grid point
c  nstart : index of lowest grid point
c  nrad   : # of gridpoints of reduced eigenfunctions
c
      rn=r(n)
      gn=con*rhobar*rn
      vn=sqrt(gn*rn)
      wn=vn/rn
      accn=1.e+20/(rhobar*rn**4)
c  normalize radius
      do i=1,n
         r(i)=r(i)/rn
         if(i.gt.1.and.dabs(r(i)-r(i-1)).lt.1.d-7) r(i)=r(i-1)
      enddo
c cut radius knots lower than max depth
      rmax=1.-dmax*tau/rn
      j=0
      do i=1,n
         j=j+1
         if(r(i).ge.rmax) goto 30
      enddo
   30 nstart=max0(j-1,1)
      j=0
      do i=nstart,n
         j=j+1
         buf(1,j)=r(i)
      enddo
      nrad=j
c open minos_bran plane output file and search mode part
      open(7,file=fflatin,form='formatted',status='old')
  1   read(7,'(a)',end=9) str
      ni=0
      do i = 1,20
        if(str(i:i).eq.'m') ni=i
      enddo
      if(ni.eq.0) goto 1
      if(str(ni:ni+3).eq.'mode') goto 2
      goto 1
  9   stop 'ERR004:eigcon: Wrong minos_bran output text file'
  2   read(7,'(a)') str
c open minos_bran binary unformatted file
      open(8,file=fbinin,form='unformatted',status='old')
c***************************************************************
c Create .eigen relations and binary eigenfunctions file
c***************************************************************
      ieig = 9
      idat = 10
      call null_eigen
      eigid_eigen = 1
      foff_eigen = 0
      ncol_eigen = 3
      if(jcom.eq.3) ncol_eigen = 7
      npar_eigen = 7
      nraw_eigen = nrad
      nrecl = (ncol_eigen*nraw_eigen+npar_eigen)*4
      call open_eigen(fout,ieig,idat,nrecl,dir,'w',ierr)
      dir_eigen = dir
      dfile_eigen = 'eigen'
c
c Main loop ----
c
 200  read(8,end=99) nn,ll,ww,qq,qc,((buf(l,lll),lll=1,n),
     +               l=2,ncol_eigen)
      read(7,*) norder_eigen,typeo,lorder_eigen,phvel_eigen,
     +          tf,per_eigen,grvel_eigen,attn_eigen
c      write(*,*) "eigcon l ", lorder_eigen
      if(nn.ne.norder_eigen.or.ll.ne.lorder_eigen) then
        write(*,*)
     +      'ERR001: eigcon: Input plane and binary files differ: ',
     +      nn,ll,norder_eigen,lorder_eigen
        stop
      endif
c check that jcom corresponds to minos_bran mode
      iflag = 0
      if(jcom.eq.4) then
        if(typeo.ne.'c') iflag = 1
        typeo_eigen = 'C'
      else if(jcom.eq.3)then
        if(typeo.ne.'s') iflag = 1
        typeo_eigen = 'S'
      else if(jcom.eq.2) then
        if(typeo.ne.'t') iflag = 1
        typeo_eigen = 'T'
      else if(jcom.eq.1) then
        if(typeo.ne.'s') iflag = 1
        typeo_eigen = 'S'
        grvel_eigen = -1.0
      else
        write(*,*) 'ERR002: eigcon: Unknown jcom ',jcom
        stop
      endif
      if(iflag.eq.1) then
        write(*,*) 'ERR003: eigcon: jcom=',jcom,
     +            ' does not fit mode ',typeo_eigen
        stop
      endif
c additional normalization V, V' or W, W' by sqrt(l(l+1))
      if(typeo_eigen.eq.'S'.or.typeo_eigen.eq.'T') then
        fl=ll
        fl1=fl+1.0
        fl3=sqrt(fl*fl1)
        do i=nstart,n
          do j =2,3
             if(typeo_eigen.eq.'T') then
               buf(j,i)=buf(j,i)/fl3
             else
               buf(j+2,i)=buf(j+2,i)/fl3
             endif
          enddo
        enddo
      endif
      if(qq.gt.0.01) then
         qq=0.5*ww/qq
      else
         qq = 0.0
      endif
      call write_eigen(ieig,ierr)
c
c  Form output buffer
c
      do k = 2,ncol_eigen
        j=0
        do i=nstart,n
           j=j+1
           buf(k,j)=buf(k,i)
        enddo
      enddo
      call put_eigen(idat,npar_eigen,param,nraw_eigen,ncol_eigen,buf,
     +               eigid_eigen,ierr)
      eigid_eigen = eigid_eigen+1
      foff_eigen = foff_eigen+nrecl
      goto 200
  99  close(7)
      close(8)

      call close_eigen(ieig,idat)
      end
