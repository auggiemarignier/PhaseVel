      module kernelfunctions

        contains

        function kernel_C(Up) result(kC)
            implicit none
            integer*4 mk
            parameter (mk=3000)
            real*4      Up(mk)
            real*4      kC(mk)

            kC = Up**2
        end function

        function kernel_A(lorder,r,U,V) result(kA)
            implicit none
            integer*4 mk
            parameter (mk=3000)
            real*4      r(mk),U(mk),V(mk)
            real*4      lorder
            real*4      kA(mk)

            kA = (1./r**2)*(2*U - lorder*(lorder+1)*V)**2 
        end function

        function kernel_F(lorder,r,U,Up,V) result(kF)
            implicit none
            integer*4 mk
            parameter (mk=3000)
            real*4      r(mk),U(mk),Up(mk),V(mk)
            real*4      lorder
            real*4      kF(mk)

            kF = (2./r)*Up*(2*U - lorder*(lorder+1)*V) 
        end function

        function kernel_L(lorder,r,U,V,Vp) result(kL)
            implicit none
            integer*4 mk
            parameter (mk=3000)
            real*4      r(mk),U(mk),V(mk),Vp(mk)
            real*4      lorder
            real*4      kL(mk)

            kL = lorder*(lorder+1)*(Vp + (U - V)/r)**2
        end function

        function kernel_N(lorder,r,U,V) result(kN)
            implicit none
            integer*4 mk
            parameter (mk=3000)
            real*4      r(mk),U(mk),V(mk),kA(mk)
            real*4      lorder
            real*4      kN(mk)

            kA = kernel_A(lorder,r,U,V)
            kN = (lorder+2)*(lorder+1)*lorder*(lorder-1)*V**2
            kN = kN - kA*r**2
            kN = (1./r**2)*kN
        end function
        
        function kernel_kappa(lorder,r,U,Up,V) result(kkappa)
            implicit none
            integer*4 mk
            parameter (mk=3000)
            real*4      r(mk),U(mk),Up(mk),V(mk)
            real*4      lorder
            real*4      kA(mk),kC(mk),kF(mk)
            real*4      kkappa(mk)

            kC = kernel_C(Up)
            kA = kernel_A(lorder,r,U,V)
            kF = kernel_F(lorder,r,U,Up,V)
            kkappa = kA + kC + kF
        end function kernel_kappa

        function kernel_mu(lorder,r,U,Up,V,Vp) result(kmu)
            implicit none
            integer*4 mk
            parameter (mk=3000)
            real*4      r(mk),U(mk),Up(mk),V(mk),Vp(mk)
            real*4      kC(mk),kA(mk),kF(mk),kL(mk),kN(mk)
            real*4      lorder
            real*4      kmu(mk)

            kL = kernel_L(lorder,r,U,V,Vp)
            kN = kernel_N(lorder,r,U,V)
            kA = kernel_A(lorder,r,U,V)
            kC = kernel_C(Up)
            kF = kernel_F(lorder,r,U,Up,V)
            kmu = kL + kN + (2./3.)*(2*kA + 2*kC - kF)
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