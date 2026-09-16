program con2prim_plots
    use iso_fortran_env, only: dp => real64
    use kastaun
    implicit none
    real(dp) :: v_in(3), v_out(3)   !3-velocity of eulerian frame
    real(dp) :: lfac_in, lfac_out !lorentz factor
    real(dp) :: rho     !mass density
    real(dp) :: bi(3)    !magnetic field
    real(dp) :: si(3)     !covariant 3-momentum
    real(dp) :: tau      !energy density
    real(dp) :: d, eps, p, h
    real(dp) :: z
    integer  :: io
    real(dp) :: logsteps_x, logsteps_y
    real(dp) :: logx_min, logx_max
    real(dp) :: logy_min, logy_max
    real(dp) :: x, y
    integer  :: steps
    integer  :: i, j
    real(dp) :: p_out, p_err, eps_out



    open(newunit=io, file="data.txt")
    !write(io, :) x_min, x_max, logsteps_x, y_min, y_max, logsteps_y

    !rho=6e12
    rho=1
    bi=[20, 0, 0]

    logx_max=3
    logx_min=-2

    logy_max=1
    logy_min=-4

    steps= 1000
    
    logsteps_x=(logx_max-logx_min)/steps

    logsteps_y=(logy_max-logy_min)/steps

    x=logx_min
    y=logy_min

    do i = 1, steps
        z=10**x

        v_in=[1, 0, 0]*z/sqrt(1+z**2)

        lfac_in=1/sqrt(1-sum(v_in)**2)
        
        do j = 1, steps

            eps=10**y

            p = ideal_eos(rho,eps)

            h = 1 + eps + p/rho

            d = rho*lfac_in

            si=(rho*h* lfac_in**2+sum(bi**2))*v_in - dot_product(bi,v_in)*bi

            tau=rho * h * lfac_in**2- p+ 0.5*sum(bi**2)* (1+sum(v_in**2))- 0.5*dot_product(bi,v_in)**2 - d

            call con2prim_extras(v_out, lfac_out, d, tau, si, bi, rho, eps_out, p_out, h)

            p_err=abs(p-p_out)/p

            !p_err=abs(p-p_out)

            write(io,*) x, y, p_err

            y=y + logsteps_y

        
        end do
        y=logy_min
        print*, i*100/steps,"% done"
        
        !print*, p_out, v_in
        x=x + logsteps_x
    end do
    


    close(io)
end program con2prim_plots