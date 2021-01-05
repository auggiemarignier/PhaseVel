      subroutine eigen2asc(nmine,nmaxe,lmine,lmaxe,fin,fdir)
      implicit none
      integer*4 mk
      parameter (mk=3000)
c ---  eigen relation common block
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
c --- other variables
      character*64 dir
      character*256 fin,fout,fdir,cmd
      real*4      rout(mk),buf(6,mk)
      real*4      U(mk),Up(mk),V(mk),Vp(mk),P(mk),Pp(mk),W(mk),Wp(mk)
      real*4      pi2,rn,vn,accn
      real*4      ww,qq
      integer*4   narg,iargc,nrecl,ieig,idat,ierr
      integer*4   nn,ll,lll,i,j
      integer*4   in,nmine,nmaxe,lmine,lmaxe,lnblnk
      logical tf
 
      pi2 = atan(1.0)*8.0
      in = 0
c --- if file fdir doesnt exist create it
      inquire(file=fdir,exist=tf)
      if(.not.tf) then
          write(cmd,'("mkdir -p ",a247)') fdir
          call system(cmd)
      endif
      nrecl = 2000
      ieig = 9
      idat = 10
c get record length of direct eigen file
      call open_eigen(fin,ieig,idat,nrecl,dir,'r',ierr)
      call read_eigen(ieig,ierr)
      nrecl = (ncol_eigen*nraw_eigen+npar_eigen)*4
      call close_eigen(ieig,idat)
      call open_eigen(fin,ieig,idat,nrecl,dir,'r',ierr)
c find record by indices n and l
  1   call read_eigen(ieig,ierr)
      if(ierr.ne.0) goto 99
      if(norder_eigen.lt.nmine.or.norder_eigen.gt.nmaxe.or.
     *   lorder_eigen.lt.lmine.or.lorder_eigen.gt.lmaxe) goto 1
      if(ncol_eigen.eq.3) then 
        read(idat,rec=eigid_eigen) nn,ll,ww,qq,rn,vn,accn,
     +     (rout(lll),W(lll),Wp(lll),lll=1,nraw_eigen)
      else
        read(idat,rec=eigid_eigen) nn,ll,ww,qq,rn,vn,accn,
     +     (rout(lll),U(lll),Up(lll),V(lll),Vp(lll),P(lll),Pp(lll),
     +      lll=1,nraw_eigen)
      endif      
      do i=1,nraw_eigen
      rout(i) = rout(i)*6371000.0
      enddo
c output data
c form output file name
      write(cmd,'(a1,".",i7,".",i7,".ASC")'),typeo_eigen(1:1),
     *           norder_eigen,lorder_eigen
      do i = 3,17
          if(cmd(i:i).eq.' ') cmd(i:i)='0'
      enddo
      fout = fdir(1:lnblnk(fdir))//'/'//cmd
      open(11,file=fout,status='unknown')
      if(in.eq.0)
     +write(11,1002) norder_eigen,lorder_eigen,typeo_eigen,
     +      eigid_eigen,per_eigen,phvel_eigen,grvel_eigen,
     +      attn_eigen,nraw_eigen,ncol_eigen
      do i = nraw_eigen,1,-1
      j = nraw_eigen+1-i
      if(ncol_eigen.eq.3) then
          write(11,1000) rout(i),W(i),Wp(i)
      else
          write(11,1001) rout(i),U(i),Up(i),V(i),Vp(i),P(i),Pp(i)
      endif
      enddo
      close(11)
      goto 1
  99  call close_eigen(ieig,idat)
 1000 format(f8.0,2e15.7)
 1001 format(f8.0,6e15.7)
 1002 format(i8,1x,i8,1x,a1,1x,i8,1x,4(f16.5,1x),i8,1x,2(i4,1x),
     +       a2,1x,i10,1x,a64,1x,a32,1x,i8,1x,a17)
      end
