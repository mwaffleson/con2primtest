module kastaun

use iso_fortran_env, only: dp => real64
implicit none
real(dp), parameter :: h0=1+1e-7_dp



!this module implements a conserved to primitive variable
!scheme from Kastaun et al. https://arxiv.org/pdf/2005.01821
!the code will be labled with equation and section numbers see the paper for reference

contains 

    pure function ideal_eos(rho,eps) result(p)
        real(dp), intent(in) :: rho, eps
        real(dp) ::p
        real(dp) :: gamma
        gamma=5.0_dp/3.0_dp
        p=(gamma-1)*rho*eps
    end function ideal_eos

    subroutine preprocessing(si, bi, s_sqr, b_sqr, s_dot_b)
        implicit none 
        real(dp), intent(in) :: si(3), bi(3)
        real(dp), intent(out) :: s_sqr, b_sqr, s_dot_b
        s_sqr=dot_product(si,si)
        b_sqr=dot_product(bi,bi)
        s_dot_b=dot_product(si,bi)
    end subroutine preprocessing

    function master_function(mu, d, tau, s_sqr, b_sqr, s_dot_b) result(f)
        implicit none
        real(dp), intent(in) :: mu, d, tau
        real(dp), intent(in) :: s_sqr, b_sqr, s_dot_b
        real(dp)  ::  q, x
        real(dp) :: r_sqr, myb_sqr, r_dot_b
        real(dp) :: r_bar_sqr, q_bar
        real(dp) :: v0_sqr
        real(dp) :: v_hat_sqr, w_hat
        real(dp) :: rho_hat, eps_hat
        real(dp) :: p_hat
        real(dp) :: a_hat
        real(dp) :: nu_a, nu_b
        real(dp) :: nu_hat
        real(dp) :: f

        
        r_sqr=s_sqr/d**2 !(23)

        q=tau/d !(23)

        myb_sqr=b_sqr/d !(24)

        r_dot_b=s_dot_b !necessary for (25)

        v0_sqr=r_sqr/h0**2/(1 + r_sqr/h0**2) ! (32)(33)

        x=1.0_dp/(1.0_dp + mu*myb_sqr) !(26)

        r_bar_sqr= r_sqr*x**2 + mu*x*(1+x)*r_dot_b**2/d**3 !(38)

        q_bar= q - 0.5_dp*myb_sqr - 0.5_dp*mu**2*x**2*(r_sqr*myb_sqr-r_dot_b**2/d**3) !(39)

        v_hat_sqr=min(v0_sqr, mu**2*r_bar_sqr) !(40)

        w_hat=1/sqrt(1 - v_hat_sqr) !(40)

        rho_hat=d/w_hat !(41)

        
        if (rho_hat<0) then

            rho_hat=0
        end if

        eps_hat=w_hat*(q_bar - mu*r_bar_sqr) + v_hat_sqr*w_hat**2/(1 + w_hat) !(42)

        !add a limiter for valid eps and rho
        if (eps_hat<0) then
            eps_hat=0
        end if

        p_hat=ideal_eos(rho_hat,eps_hat) !(43)


        a_hat=p_hat/rho_hat/(1 + eps_hat) !(43)

        nu_a=(1 + a_hat)*(1 + eps_hat)/w_hat !(46)

        nu_b=(1 + a_hat)*(1 + q_bar - mu*r_bar_sqr) !(47)
         
        nu_hat=max(nu_a, nu_b) !(48)

        f=mu-1.0_dp/(nu_hat + mu*r_bar_sqr)!(44)(45)


    end function master_function
        
    !bisection algorithm to find roots of f
    function bisection(mu_plus, d, tau, s_sqr, b_sqr, s_dot_b) result(mu)
        implicit none
        real(dp), intent(in) :: d, tau
        real(dp), intent(in) :: s_sqr, b_sqr, s_dot_b
        real(dp), intent(in) :: mu_plus
        real(dp) :: mu
        real(dp) :: a, b, c
        real(dp) :: fa, fb, fc
        real(dp) :: tol
        integer :: maxiter, i
        logical :: success, bound_error



        success=.false.
        bound_error=.false.
        maxiter=1000
        tol=1e-10

        a=0

        if (mu_plus==0) then
            b=1/h0
        else
            b=mu_plus
        end if

        fa=master_function(a, d, tau, s_sqr, b_sqr, s_dot_b)

        fb=master_function(b, d, tau, s_sqr, b_sqr, s_dot_b)

        if (fa>0 .eqv. fb>0) then 
            bound_error=.true.
        end if

        do i=1, maxiter
            if (bound_error) then
                exit
            end if

            c=(b+a)/2.0_dp

            fc=master_function(c, d, tau, s_sqr, b_sqr, s_dot_b)

            if (abs((b-a)/c)<tol*mu) then
                mu=c
                success = .true.
                exit
            end if

            if (fa>0 .and. fc>0) then
                a=c
                fa=fc

            else if(.not. fa>0 .and. .not. fc>0) then
                a=c
                fa=fc
            else 
                b=c
                fb=fc
            end if


        end do

        if (.not. success .and. .not. bound_error) then
            mu=c
        end if
    
    end function bisection

    function aux_f(mu, d, s_sqr, s_dot_b, b_sqr) result(f)
        implicit none
        real(dp), intent(in) :: mu
        real(dp), intent(in) :: d, s_sqr, s_dot_b, b_sqr
        real(dp)             :: mymu
        real(dp)             :: r_dot_b, myb_sqr, myd
        real(dp)             :: f
        real(dp)             :: r_bar_sqr, r_sqr, x


        mymu=mu
        r_dot_b=s_dot_b
        myd=d
        myb_sqr=b_sqr/myd
        r_sqr=s_sqr/myd**2

        x=1/(1+mu*myb_sqr)

        r_bar_sqr= r_sqr*x**2 + mu*x*(1+x)*r_dot_b**2/myd**3 !(38)

        f=mu*sqrt(h0**2+r_bar_sqr)-1
    end function aux_f

    function aux_bisection(d, s_sqr, s_dot_b, b_sqr) result(b)
        implicit none
        real(dp), intent(in) :: d
        real(dp), intent(in) :: s_sqr, b_sqr, s_dot_b
        real(dp) :: mu
        real(dp) :: a, b, c
        real(dp) :: fa, fb, fc
        real(dp) :: tol
        integer :: maxiter, i
        logical :: success, bound_error

        success=.false.
        bound_error=.false.
        maxiter=1000
        tol=1e-7

        a=0

        b=1/h0

        fa=aux_f(a, d, s_sqr, b_sqr, s_dot_b)

        fb=aux_f(b, d, s_sqr, b_sqr, s_dot_b)

        if (fa>0 .eqv. fb>0) then 
            bound_error=.true.
        end if

        do i=1, maxiter
            if (bound_error) then
                exit
            end if


            c=(b+a)/2.0_dp


            fc=aux_f(c, d, s_sqr, b_sqr, s_dot_b)

            if (abs(fc)<=tol) then
                mu=c
                success = .true.
                exit
            end if

            if (fa>0 .and. fc>0) then
                a=c
                fa=fc
            else if(.not. fa>0 .and. .not. fc>0) then
                a=c
                fa=fc
            else 
                b=c
                fb=fc
            end if
            
        end do

        if (.not. success .and. .not. bound_error) then
        end if
    
    end function aux_bisection


    subroutine con2prim(vi, lfac, d, tau, si, bi)
        implicit none
        real(dp), intent(out) :: lfac
        !real(dp), intent(out) ::rho, eps, p, h
        real(dp), intent(out) :: vi(3)
        real(dp), intent(in)  :: d, tau
        real(dp), intent(in)  :: si(3), bi(3)
        real(dp)              :: myd, mytau
        real(dp)              :: mysi(3), mybi(3)
        real(dp)              :: s_sqr, b_sqr, s_dot_b
        real(dp)              :: mu, mu_plus
        !real(dp)              :: q_bar, r_bar_sqr
        !real(dp)              :: x
        

        myd=d
        mytau=tau
        mysi=si
        mybi=bi
        call preprocessing(mysi, mybi, s_sqr, b_sqr, s_dot_b)

        if (s_sqr/myd**2<h0**2) then
            mu_plus=0.0_dp
        else
            mu_plus=aux_bisection(myd, s_sqr, s_dot_b, b_sqr)
            !mu_plus=1/h0
        end if
        


        mu= bisection(mu_plus, myd, mytau, s_sqr, b_sqr, s_dot_b)

        vi= mu/(1+mu*b_sqr/myd)*(si/myd + mu*s_dot_b*mybi/myd**2)

        lfac=1/sqrt(1 - sum(vi**2))

        !rho=d/lfac

        !x=1.0_dp/(1.0_dp + mu*b_sqr) !(26)

        !r_bar_sqr= s_sqr*x**2/d**2 + mu*x*(1+x)*s_dot_b**2 !(38)

        !q_bar=  tau/d - 0.5_dp*b_sqr - 0.5_dp*mu**2*x**2*(s_sqr*myb_sqr/d**2-s_dot_b**2/d**3) !(39)

        !eps=lfac*(q_bar - mu*r_bar_sqr) + sum(vi**2)*lfac**2/(1 + lfac) !(42)

        !p=ideal_eos(eps,rho)

        !h = 1 + eps + p/rho

    end subroutine

end module kastaun