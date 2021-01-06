      module kernelfunctions

        contains
      
        function kernel_kappa(omega,wavenum,r,U,Up,V) result(kkappa)
            implicit none
            integer*4 mk
            parameter (mk=3000)
            real*4      r(mk),U(mk),Up(mk),V(mk)
            real*4      omega,wavenum
            real*4      kkappa(mk)

            kkappa = (r*Up + 2*U - wavenum*V)**2
            kkappa = (1./(2.*omega))*kkappa
        end function kernel_kappa

        function kernel_mu(omega,wavenum,r,U,Up,V,Vp,W,Wp) result(kmu)
            implicit none
            integer*4 mk
            parameter (mk=3000)
            real*4      r(mk),U(mk),Up(mk),V(mk),Vp(mk),W(mk),Wp(mk)
            real*4      omega,wavenum
            real*4      kmu(mk)

            kmu = (1./3.)*(2*r*Up -2*U +wavenum*V)**2
            kmu = kmu + (r*Vp - V + wavenum*U)**2
            kmu = kmu + (r*Wp - W)**2
            kmu = kmu + (wavenum**2 - 2)*(V**2 + W**2)
            kmu = (1./(2.*omega))*kmu
        end function kernel_mu

        function kernel_alpha(alpha,rho,kkappa) result(kalpha)
            implicit none
            integer*4 mk
            parameter (mk=3000)
            real*8      alpha(mk),rho(mk)
            real*4      kkappa(mk)
            real*4      kalpha(mk)

            kalpha = 2*rho*alpha*kkappa
        end function kernel_alpha
        
        function kernel_beta(beta,rho,kkappa,kmu) result(kbeta)
            implicit none
            integer*4 mk
            parameter (mk=3000)
            real*8      beta(mk),rho(mk)
            real*4      kkappa(mk),kmu(mk)
            real*4      kbeta(mk)

            kbeta = 2*rho*beta*(kmu - (4./3.)*kkappa)
        end function kernel_beta

      end module kernelfunctions

      subroutine write_kernels_asc(fdir,rout,kkappa,kmu,kalpha,kbeta)
        implicit none
        integer*4 mk
        parameter (mk=3000)
        real*4    per_eigen,phvel_eigen,grvel_eigen,attn_eigen
        integer*4 norder_eigen,lorder_eigen,eigid_eigen,
     +               nraw_eigen,ncol_eigen,npar_eigen,foff_eigen,
     +               commid_eigen
        character*2 datatype_eigen
        character*64 dir_eigen
        character*32 dfile_eigen
        character*17 lddate_eigen
        character*1 typeo_eigen
        common/c_eigen/norder_eigen,lorder_eigen,
     +           eigid_eigen,per_eigen,phvel_eigen,grvel_eigen,
     +           attn_eigen,nraw_eigen,ncol_eigen,npar_eigen,
     +           foff_eigen,commid_eigen,typeo_eigen,
     +           datatype_eigen,dir_eigen,dfile_eigen,lddate_eigen
        character*256 fout,fdir,cmd
        real*4      rout(mk)
        real*4      kkappa(mk),kmu(mk),kalpha(mk),kbeta(mk)
        integer i
        logical tf

        inquire(file=fdir,exist=tf)
        if(.not.tf) then
                write(cmd,'("mkdir -p ",a247)') fdir
            call system(cmd)
        endif
          
        write(cmd,'(a1,".",i7,".",i7,".ASC")'),typeo_eigen(1:1),
     *            norder_eigen,lorder_eigen
        do i = 3,17
            if(cmd(i:i).eq.' ') cmd(i:i)='0'
        enddo
        fout = fdir(1:lnblnk(fdir))//'/'//cmd
        open(11,file=fout,status='unknown')
        do i = nraw_eigen,1,-1
            write(11,1000) 
     1          rout(i),kkappa(i),kmu(i),kalpha(i),kbeta(i)
        enddo
        close(11)
 1000   format(f8.0,4e15.7)
        end subroutine