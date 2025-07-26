module gfn0
    use gfn0_interface
    implicit none

contains


subroutine gfn0_sp(nat, uhf, ichrg, numbers, positions, fail, energy, grad)
    integer,intent(in) :: nat
    integer,intent(in) :: uhf
    integer,intent(in) :: ichrg
    integer,intent(in) :: numbers(nat)
    double precision,intent(in) :: positions(3, nat)      ! Bohr
    !f2py depend(positions) :: nat=shape(positions,1)
    !f2py depend(numbers) :: nat=size(numbers)

    logical,intent(out) :: fail
    double precision,intent(out) :: energy          ! Hartree
    double precision,intent(out) :: grad(3, nat)    ! Hartree / Bohr
    type(gfn0_results) :: res
    type(gfn0_data) :: gdat

    call gfn0_init(nat,numbers,positions,ichrg,uhf,gdat)
    !call gfn0_init(nat,numners,positions,ichrg,uhf,gdat,solv='h2o',alpb=.false.)
    call gfn0_singlepoint(nat,numbers,positions,ichrg,uhf,gdat,energy,grad,fail,res)
end subroutine gfn0_sp


! subroutine gfn0_alpb(nat, ichrg, numbers, posotions, sol, energy, grad, io)
!     ......
! end subroutine gfn0_alpb

end module gfn0
