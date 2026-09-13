program con2prim_test
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
!real(dp) :: eps_out, p_out, h_out
integer :: tests
real(dp), allocatable :: vs(:,:)
real(dp), allocatable :: bs(:,:)
real(dp), allocatable :: rhos(:)
real(dp), allocatable :: epss(:)
character(len=15), allocatable :: test_names(:)
integer :: i

tests=10
allocate(vs(tests,3))
allocate(bs(tests,3))
allocate(rhos(tests))
allocate(epss(tests))
allocate(test_names(tests))

rhos(:)=1
epss(:)=0.5
vs=0
bs=0

! vs(:,1)=0.1_dp
! test_names(1)="Default"


! vs(2,:)=[0.0_dp, 0.0_dp, 0.3_dp]
! test_names(2)="0.3 v"

! vs(3,:)=[-0.05_dp, 0.2_dp, 0.1_dp]
! test_names(3)="mixed v"

! bs(4,:)=[20, 0, 0]
! test_names(4)="Big B"

! vs(5,:)=[0.45_dp, 0.0_dp, 0.0_dp]
! test_names(5)="0.45 v"

! vs(6,:)=[0.0_dp, 0.0_dp, 0.0_dp]
! test_names(6)="0 v"

! rhos(7)=100
! bs(7,:)=[2, 0, 0]
! test_names(7)="big rho"

! vs(8,:)= [0.6_dp, 0.0_dp, 0.0_dp]
! test_names(8)= "big v"

do i=1, tests-1
    vs(i,:)= [i*0.1_dp, 0.0_dp, 0.0_dp]
end do



vs(tests,1)=0.9975

do i=1, tests
    print *, "Test case",i,", ", trim(test_names(i))
    rho=rhos(i)
    eps=epss(i)
    bi=bs(i,:)
    v_in=vs(i,:)
    lfac_in = 1/sqrt(1-sum(v_in**2))

    p = ideal_eos(rho,eps)

    h = 1 + eps + p/rho

    d = rho*lfac_in

    si=(rho*h* lfac_in**2+sum(bi**2))*v_in - dot_product(bi,v_in)*bi

    tau=rho * h *lfac_in**2- p+ 0.5*sum(bi**2)* (1+sum(v_in**2))- 0.5*dot_product(bi,v_in) - d


    call con2prim(v_out, lfac_out, d, tau, si, bi)
    !call bound_test(v_out, lfac_out, d, tau, si, bi)

    print *, "Rho=", rho, " Epsilon=", eps
    print *, "B=", bi
    print *, "v in:", v_in, "v out:", v_out
    print *, "Absolute v error:", abs(norm2(v_out)-norm2(v_in))
    print *, "lfac in:", lfac_in," lfac out:", lfac_out
    print *, "Absolute lfac error:", abs(lfac_in-lfac_out)
    print *, " "
    print *, " "

end do





end program con2prim_test