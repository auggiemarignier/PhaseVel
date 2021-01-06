      module kernelfunctions

        contains
      
        function kernel_kappa(omega,wavenum,r,U,Up,V) result(kkappa)
        implicit none
        integer*4 mk
        parameter (mk=3000)
        real*4      r(mk),U(mk),Up(mk),V(mk)
        real*4      omega,wavenum
        real*4      kkappa(mk)

        print*,"calculating kkappa",size(U),size(Up)
        kkappa = (r*Up + 2*U - wavenum*V)**2
        kkappa = (1./2.*omega)*kkappa
        end function

      end module
