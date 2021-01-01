c
c MINEOS version 1.0 by Guy Masters, John Woodhouse, and Freeman Gilbert
c
c This program is free software; you can redistribute it and/or modify
c it under the terms of the GNU General Public License as published by
c the Free Software Foundation; either version 2 of the License, or
c (at your option) any later version.
c
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU General Public License for more details.
c
c You should have received a copy of the GNU General Public License
c along with this program; if not, write to the Free Software
c Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
c
c**************************************************************************
c
c preigen program prints out eigen functions by mode n and closest period T.
c
c**************************************************************************
      subroutine eigen2asc(str8,str9,Amp,W,S,eigenvalue)
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
      real*4      pi2,rn,vn,accn
      real*4      ww,qq
      integer*4   narg,iargc,nrecl,ieig,idat,ierr
      integer*4   nn,ll,lll,i,j
      integer*4   in,nmine,nmaxe,lmine,lmaxe,lnblnk
      logical tf
      real*8      Amp
      real*8      W(299),S(299),Amp0(299)
      integer*4   q
      character*256   str8,str9
      character*1 eigenvalue

   
c      write(*,*) "str8 ", str8
c      write(*,*) "str9 ",str9
c      write(*,*) "Amp ",Amp
c      write(*,*) "W ", W
c      write(*,*) "S ", S
c      write(*,*) "eigenvalue ",eigenvalue
c      write(*,*) "Amp ",Amp

      nmine=0
      nmaxe=0
      lmine=0
      lmaxe=300
      pi2 = atan(1.0)*8.0
      in = 0
c get input arguments
c      call getarg(1,fin)
c      if(fin(1:2).eq.'-n') in = 1
c      narg=iargc()
c      if (narg.ne.6+in) then
c        write(*,*) ' Usage: eigen2asc nmin nmax lmin lmax db_name ',
c     *                'out_dir '
c        write(*,*) '        where,'
c        write(*,*) '        nmin, nmax - mininal and maximal values ',
c     *                 'of the radial mode number n'
c        write(*,*) '        lmin, lmax - mininal and maximal values ',
c     *                 'of the lateral mode number l'
c        write(*,*) '        db_name    - data base name, must include',
c     *                 ' .eigen relation;'
c        write(*,*) '        out_dir    - output directory for ASCII ',
c     *                 'files'
c        call exit(1)
c      endif
c ---
c      call getarg(1+in,fin)
c      read(fin,*) nmine
c      call getarg(2+in,fin)
c      read(fin,*) nmaxe
c      call getarg(3+in,fin)
c      read(fin,*) lmine
c      call getarg(4+in,fin)
c      read(fin,*) lmaxe
c      call getarg(5+in,fin)
c      call getarg(6+in,fdir)
c --- if file fdir doesn't exist create it
      inquire(file=str9,exist=tf)
      if(.not.tf) then
          write(cmd,'("mkdir -p ",a247)') str9
          call system(cmd)
      endif
      nrecl = 2000
      ieig = 9
      idat = 10
c get record length of direct eigen file
      call open_eigen(str8,ieig,idat,nrecl,dir,'r',ierr)
      call read_eigen(ieig,ierr)
      nrecl = (ncol_eigen*nraw_eigen+npar_eigen)*4
      call close_eigen(ieig,idat)
      call open_eigen(str8,ieig,idat,nrecl,dir,'r',ierr)
c find record by indices n and l
      q=1
  1   call read_eigen(ieig,ierr)
      if(ierr.ne.0) goto 99
      if(norder_eigen.lt.nmine.or.norder_eigen.gt.nmaxe.or.
     *   lorder_eigen.lt.lmine.or.lorder_eigen.gt.lmaxe) goto 1
c      write(*,*)"e2a lorder ",lorder_eigen
c      write(*,*)"e2a l,gr,U ",lorder_eigen,grvel_eigen,buf(1,nraw_eigen)
      read(idat,rec=eigid_eigen) nn,ll,ww,qq,rn,vn,accn,
     +     (rout(lll),(buf(ll,lll),ll=1,ncol_eigen-1),lll=1,nraw_eigen)
      do i=1,nraw_eigen
        rout(i) = rout(i)*6371000.0
      enddo
c output data
c form output file name
c      write(cmd,'(a1,".",i7,".",i7,".ASC")'),typeo_eigen(1:1),
c     *           norder_eigen,lorder_eigen
c      do i = 3,17
c         if(cmd(i:i).eq.' ') cmd(i:i)='0'
c      enddo
c      fout = fdir(1:lnblnk(fdir))//'/'//cmd
c      open(11,file=fout,status='unknown')
c      if(in.eq.0)
c     +write(11,1002) norder_eigen,lorder_eigen,typeo_eigen,
c     +      eigid_eigen,per_eigen,phvel_eigen,grvel_eigen,
c     +      attn_eigen,nraw_eigen,ncol_eigen
c      do i = nraw_eigen,1,-1
c        j = nraw_eigen+1-i
c        if(ncol_eigen.eq.3) then
c          write(11,1000) rout(i),(buf(lll,i),lll=1,ncol_eigen-1)
c        else
c          write(11,1001) rout(i),(buf(lll,i),lll=1,ncol_eigen-1)
c        endif
c      enddo
c      close(100)

      Amp=ABS(buf(1,nraw_eigen))/(sqrt(grvel_eigen))
c      write(*,*) "per_eigen ", per_eigen
c      write(*,*) "nraw_eigen ", buf(1,nraw_eigen)
c      write(*,*) "sqrt(grvel) ", sqrt(grvel_eigen)
c      write(*,*) "eigen2asc amp ", Amp
      W(299-q)=per_eigen
      S(299-q)=Amp
c      write(*,*) "q-------------------", q
      q=q+1
c      write(*,*) "q2-------------------", q
c      W(299-q)=per_eigen
c      S(299-q)=Amp
c      write(*,*) per_eigen, Amp


      goto 1
  99  call close_eigen(ieig,idat)

 1000 format(f8.0,2e15.7)
 1001 format(f8.0,6e15.7)
 1002 format(i8,1x,i8,1x,a1,1x,i8,1x,4(f16.5,1x),i8,1x,2(i4,1x),
     +       a2,1x,i10,1x,a64,1x,a32,1x,i8,1x,a17)
      end subroutine
