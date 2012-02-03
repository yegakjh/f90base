!----------------------------------------------------------------
! Category : module
! Name : ode_dverk
! Purpose : Define workspace and variables for DVERK
!----------------------------------------------------------------
! DVERK for Time Integration
! 03/07/07 DVERK

module odesolver
  use mod_const
  use mod_four2d
  use set_ivp
  use myfield

  use clock

  implicit none
  !! Required Parameters
  integer :: indicator

  integer, parameter :: nparam=28 ! original 1:24
  ! 26,27, 28 are time-info
  real(dp) :: tolerance, time_limit
  real(dp), dimension(nparam) :: tparam
  integer :: nw
  real(dp), dimension(neq,9) :: wrk

contains
  subroutine odesolver_open(tol,utim, mxstep,mxfcn)
    implicit none
    real(dp), intent(in) :: tol, utim
    integer, intent(in) :: mxstep, mxfcn
    
    tparam = zero
    tparam(1) = 1
    tparam(7) = mxfcn
    tolerance = tol
    time_limit = utim
    nw = neq
    return
  end subroutine odesolver_open

  subroutine odesolver_driver(t,tn)
    real(dp) :: t, tn
    type(myclock) :: t0, t1

    indicator = 1
    call time_check(t0)
    call dverk(neq,fcn2d,t,phik,tn,tolerance,indicator,tparam(1:24),neq,wrk)
    call time_check(t1,t0)

    if (indicator.lt.1) then
       write(*,*) "DVERK Ind = ",indicator
       write(*,'("22-24: ",3(e10.3,1x))') tparam(22:24)
    end if

    tparam(26) = t1%cpu
    tparam(27) = t1%wall
    tparam(28) = tparam(26)/tparam(27)

    return
  end subroutine odesolver_driver
end module odesolver

