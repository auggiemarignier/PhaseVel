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
            real*4      alpha(mk),rho(mk),kkappa(mk)
            real*4      kalpha(mk)

            kalpha = 2*rho*alpha*kkappa
        end function kernel_alpha
        
        function kernel_beta(beta,rho,kkappa,kmu) result(kbeta)
            implicit none
            integer*4 mk
            parameter (mk=3000)
            real*4      beta(mk),rho(mk),kkappa(mk),kmu(mk)
            real*4      kbeta(mk)

            kbeta = 2*rho*beta*(kmu - (4./3.)*kkappa)
        end function kernel_beta

      end module kernelfunctions
