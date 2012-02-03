! 2/3 Rule
module four_fftw3_2d_r2c_jacb23
  use mydef
  use four_fftw3_def
  use mod_const
  
  implicit none
  
  private 
  integer, parameter :: nomp=4
  
  public :: myjacb_open, myjacb, myjacb_close
  
  integer*8 :: fplanfwd,  fplanbwd(nomp)
  
  integer, parameter :: n1=nx, n2=ny, d1=nx+2, d2=ny
  integer, parameter :: on1=nx, on2=ny, od1=nx+2, od2=ny
  
  real(dp) :: scale
  logical :: init=.FALSE.
  character(len=*), parameter :: modn_jacb='(four_fftw3_2d_r2c_jacb)'
  
  complex(dpc) :: derx(d1/2), dery(d2)
  complex(dpc), dimension(od1/2, od2) :: fx, fy, gx, gy
  real(dp), dimension(od1, od2) :: rfx, rfy, rgx, rgy
  
  equivalence (fx, rfx), (fy, rfy), (gx, rgx), (gy, rgy)
  
contains
  
  subroutine myjacb_open(derx0, dery0)
    implicit none
    complex(dpc), intent(in) :: derx0(d1/2), dery0(d2)
    integer, dimension(2) :: ta, inemb, outemb
    integer :: i1
    
    derx = derx0
    dery = dery0
    
    write(fftout,*) " MYJACB initialization"//modn_jacb
    write(fftout,*) " ( nx,  ny)=", nx, ny
    
    ta(1) = n1; ta(2) = n2
    inemb(1) = d1; inemb(2) = d2
    outemb(1) = d1/2; outemb(2) = d2
    
    scale = one/real(n1*n2,dp)
    
    write(fftout,*) 'ta = ', ta
    write(fftout,*) 'inemb = ', inemb
    write(fftout,*) 'outemb = ', outemb
    
    call dfftw_plan_many_dft_r2c(fplanfwd,2, ta, 1, rfx, inemb, 1, 1, fx, outemb, 1,1, wisdom)
    
    do i1=1, nomp
       call dfftw_plan_many_dft_c2r(fplanbwd(i1),2, ta, 1, fx, outemb, 1, 1, rfx, inemb, 1, 1,wisdom)
    enddo
    
    write(fftout,*) "fwd plan =   ",fplanfwd
    write(fftout,*) "bckwd plan = ",fplanbwd
    
    init = .TRUE.
    return
  end subroutine myjacb_open
  
  subroutine myjacb(f,g,h)
    implicit none
    complex(dpc), dimension(d1/2, d2), intent(in) :: f,g
    complex(dpc), dimension(d1/2, d2), intent(out) :: h
    integer :: i1
    
    if (.not. init) stop modn_jacb

    !$omp parallel sections
    !$omp section
    do i1=1, d2
       fx(:,i1) = derx*f(:,i1)
    end do
    call dfftw_execute_dft_c2r(fplanbwd(1), fx, rfx)    
    
    !$omp section
    do i1=1, d2
       fy(:,i1) = dery(i1)*f(:,i1)
    end do
    call dfftw_execute_dft_c2r(fplanbwd(2), fy, rfy)    
    
    !$omp section
    do i1=1, d2
       gx(:,i1) = derx*g(:,i1)
    end do
    call dfftw_execute_dft_c2r(fplanbwd(3), gx, rgx)    
    
    !$omp section
    do i1=1, d2
       gy(:,i1) = dery(i1)*g(:,i1)
    end do
    call dfftw_execute_dft_c2r(fplanbwd(4), gy, rgy)    
    !$omp end parallel sections
    
    rgy = rfx*rgy-rfy*rgx
    
    call dfftw_execute_dft_r2c(fplanfwd, rgy, h)
    h = scale*h
    
    return
  end subroutine myjacb
  
  subroutine myjacb_close
    integer :: i1
    if (.not.init) stop modn_jacb
    call dfftw_destroy_plan(fplanfwd)
    do i1=1, nomp
       call dfftw_destroy_plan(fplanbwd(i1))
    end do
    
    return
  end subroutine myjacb_close
end module four_fftw3_2d_r2c_jacb23

module four_fftw3_2d_r2c_jacb23_2x
  use mydef
  use four_fftw3_def
  use mod_const
  
  implicit none
  
  private 
  integer, parameter :: nomp=8
  
  public :: myjacb_open, myjacb, myjacb_close, myjacb_2x
  
  integer*8 :: fplanfwd,  fplanbwd(nomp)
  
  integer, parameter :: n1=nx, n2=ny, d1=nx+2, d2=ny
  integer, parameter :: on1=nx, on2=ny, od1=nx+2, od2=ny
  
  real(dp) :: scale
  logical :: init=.FALSE.
  character(len=*), parameter :: modn_jacb='(four_fftw3_2d_r2c_jacb23_2x)'
  
  complex(dpc) :: derx(d1/2), dery(d2)
  
  complex(dpc), dimension(od1/2, od2) :: f1x, f1y, g1x, g1y
  real(dp), dimension(od1, od2) :: rf1x, rf1y, rg1x, rg1y
  equivalence (f1x, rf1x), (f1y, rf1y), (g1x, rg1x), (g1y, rg1y)

  complex(dpc), dimension(od1/2, od2) :: f2x, f2y, g2x, g2y
  real(dp), dimension(od1, od2) :: rf2x, rf2y, rg2x, rg2y
  equivalence (f2x, rf2x), (f2y, rf2y), (g2x, rg2x), (g2y, rg2y)  
  
contains
  
  subroutine myjacb_open(derx0, dery0)
    implicit none
    complex(dpc), intent(in) :: derx0(d1/2), dery0(d2)
    integer, dimension(2) :: ta, inemb, outemb
    integer :: i1
    
    derx = derx0
    dery = dery0
    
    write(fftout,*) " MYJACB initialization"//modn_jacb
    write(fftout,*) " ( nx,  ny)=", nx, ny
    
    ta(1) = n1; ta(2) = n2
    inemb(1) = d1; inemb(2) = d2
    outemb(1) = d1/2; outemb(2) = d2
    
    scale = one/real(n1*n2,dp)
    
    write(fftout,*) 'ta = ', ta
    write(fftout,*) 'inemb = ', inemb
    write(fftout,*) 'outemb = ', outemb
    
    call dfftw_plan_many_dft_r2c(fplanfwd,2, ta, 1, rf1x, inemb, 1, 1, f1x, outemb, 1,1, wisdom)
    
    do i1=1, nomp
       call dfftw_plan_many_dft_c2r(fplanbwd(i1),2, ta, 1, f1x, outemb, 1, 1, rf1x, inemb, 1, 1,wisdom)
    enddo
    
    write(fftout,*) "fwd plan =   ",fplanfwd
    write(fftout,*) "bckwd plan = ",fplanbwd
    
    init = .TRUE.
    return
  end subroutine myjacb_open
  
  subroutine myjacb(f,g,h)
    implicit none
    complex(dpc), dimension(d1/2, d2), intent(in) :: f,g
    complex(dpc), dimension(d1/2, d2), intent(out) :: h
    integer :: i1
    
    if (.not. init) stop modn_jacb

    !$omp parallel sections
    !$omp section
    do i1=1, d2
       f1x(:,i1) = derx*f(:,i1)
    end do
    call dfftw_execute_dft_c2r(fplanbwd(1), f1x, rf1x)    
    
    !$omp section
    do i1=1, d2
       f1y(:,i1) = dery(i1)*f(:,i1)
    end do
    call dfftw_execute_dft_c2r(fplanbwd(2), f1y, rf1y)    
    
    !$omp section
    do i1=1, d2
       g1x(:,i1) = derx*g(:,i1)
    end do
    call dfftw_execute_dft_c2r(fplanbwd(3), g1x, rg1x)    
    
    !$omp section
    do i1=1, d2
       g1y(:,i1) = dery(i1)*g(:,i1)
    end do
    call dfftw_execute_dft_c2r(fplanbwd(4), g1y, rg1y)    
    !$omp end parallel sections
    
    rg1y = rf1x*rg1y-rf1y*rg1x
    
    call dfftw_execute_dft_r2c(fplanfwd, rg1y, h)
    h = scale*h
    
    return
  end subroutine myjacb

  subroutine myjacb_2x(f1,g1,h1, f2, g2, h2)
    implicit none
    complex(dpc), dimension(d1/2, d2), intent(in) :: f1,g1, f2,g2
    complex(dpc), dimension(d1/2, d2), intent(out) :: h1, h2
    integer :: i1
    
    if (.not. init) stop modn_jacb
    
    !$omp parallel sections
    !$omp section
    do i1=1, d2
       f1x(:,i1) = derx*f1(:,i1)
    end do
    call dfftw_execute_dft_c2r(fplanbwd(1), f1x, rf1x)    
    
    !$omp section
    do i1=1, d2
       f1y(:,i1) = dery(i1)*f1(:,i1)
    end do
    call dfftw_execute_dft_c2r(fplanbwd(2), f1y, rf1y)    
    
    !$omp section
    do i1=1, d2
       g1x(:,i1) = derx*g1(:,i1)
    end do
    call dfftw_execute_dft_c2r(fplanbwd(3), g1x, rg1x)    
    
    !$omp section
    do i1=1, d2
       g1y(:,i1) = dery(i1)*g1(:,i1)
    end do
    call dfftw_execute_dft_c2r(fplanbwd(4), g1y, rg1y)    

    !$omp section
    do i1=1, d2
       f2x(:,i1) = derx*f2(:,i1)
    end do
    call dfftw_execute_dft_c2r(fplanbwd(5), f2x, rf2x)    
    
    !$omp section
    do i1=1, d2
       f2y(:,i1) = dery(i1)*f2(:,i1)
    end do
    call dfftw_execute_dft_c2r(fplanbwd(6), f2y, rf2y)    
    
    !$omp section
    do i1=1, d2
       g2x(:,i1) = derx*g2(:,i1)
    end do
    call dfftw_execute_dft_c2r(fplanbwd(7), g2x, rg2x)    
    
    !$omp section
    do i1=1, d2
       g2y(:,i1) = dery(i1)*g2(:,i1)
    end do
    call dfftw_execute_dft_c2r(fplanbwd(8), g2y, rg2y)    
    
    !$omp end parallel sections
    
    rg1y = rf1x*rg1y-rf1y*rg1x
    rg2y = rf2x*rg2y-rf2y*rg2x
    
    call dfftw_execute_dft_r2c(fplanfwd, rg1y, h1)
    call dfftw_execute_dft_r2c(fplanfwd, rg2y, h2)

    h1 = scale*h1
    h2 = scale*h2
    
    return
  end subroutine myjacb_2x
  
  subroutine myjacb_close
    integer :: i1
    if (.not.init) stop modn_jacb
    call dfftw_destroy_plan(fplanfwd)
    do i1=1, nomp
       call dfftw_destroy_plan(fplanbwd(i1))
    end do
    
    return
  end subroutine myjacb_close
end module four_fftw3_2d_r2c_jacb23_2x



! this is for 3/2
!     if ( n1 < on1) then
!        ta(1) = on1; ta(2) = on2
!        inemb(1) = od1; inemb(2) = od2
!        outemb(1) = od1/2; outemb(2) = od2
       
!        scale = one/real(on1*on2,dp)
       
!        write(fftout,*) 'ta = ', ta
!        write(fftout,*) 'inemb = ', inemb
!        write(fftout,*) 'outemb = ', outemb
       
!     else


!   subroutine myjacb32(derx, dery, f, g,h)
!     implicit none
!     complex(dpc), intent(in) :: derx(d1/2), dery(d2)
!     complex(dpc), dimension(d1/2, d2), intent(in) :: f,g
!     complex(dpc), dimension(d1/2, d2), intent(out) :: h
!     complex(dpc), dimension(od1/2, od2) :: fx, gx, fy,gy
!     integer :: i1
!     integer :: omp_get_num_threads, omp_get_thread_num

!     if (.not.jacb_init) stop modn_jacb
!     if (d1 > od1 .or. d2  > od2) stop 'This subroutine is to use 3/2 method'

!     !$omp parallel sections 
!     !$omp section
!     do i1=1, d2/2-1
!        fx(1:d1/2-1,i1) = derx(1:d1/2-1)*f(1:d1/2-1,i1)
!        fx(1:d1/2-1,od2+1-i1) = derx(1:d1/2-1)*f(1:d1/2-1,d2+1-i1)
!        fx(d1/2:od1/2,i1) = zero
!        fx(d1/2:od1/2,od2+1-i1) = zero
!     end do
!     fx(:,d2/2:od2-d2/2+1) =zero

!     call dfftw_execute_dft_c2r(fplanbwd(1), fx, fx)    

!     !$omp section

!     do i1=1, d2/2-1
!        fy(1:d1/2-1,i1) = dery(i1)*f(1:d1/2-1,i1)
!        fy(1:d1/2-1,od2+1-i1) = dery(d2+1-i1)*f(1:d1/2-1,d2+1-i1)
!        fy(d1/2:od1/2,i1) = zero
!        fy(d1/2:od1/2,od2+1-i1) = zero
!     end do
!     fy(:,d2/2:od2-d2/2+1) =zero

!     call dfftw_execute_dft_c2r(fplanbwd(1), fy, fy)    
    
!     !$omp section 

!     do i1=1, d2/2-1
!        gx(1:d1/2-1,i1) = derx(1:d1/2-1)*g(1:d1/2-1,i1)
!        gx(1:d1/2-1,od2+1-i1) = derx(1:d1/2-1)*g(1:d1/2-1,d2+1-i1)
!        gx(d1/2:od1/2,i1) = zero
!        gx(d1/2:od1/2,od2+1-i1) = zero
!     end do
!     gx(:,d2/2:od2-d2/2+1) =zero

!     call dfftw_execute_dft_c2r(fplanbwd(1), gx, gx)    

!     !$omp section
!     do i1=1, d2/2-1
!        gy(1:d1/2-1,i1) = dery(i1)*g(1:d1/2-1,i1)
!        gy(1:d1/2-1,od2+1-i1) = dery(d2+1-i1)*g(1:d1/2-1,d2+1-i1)
!        gy(d1/2:od1/2,i1) = zero
!        gy(d1/2:od1/2,od2+1-i1) = zero
!     end do
!     gy(:,d2/2:od2-d2/2+1) =zero

!     call dfftw_execute_dft_c2r(fplanbwd(1), gy, gy)    

!     !$omp end parallel sections

!     gy = (real(fx,dp)*real(gy,dp)-real(fy,dp)*real(gx,dp)) &
!          +ic*(dimag(fx)*dimag(gy)-dimag(fy)*dimag(gx))
    
!     call dfftw_execute_dft_r2c(fplanfwd, gy, gy)
!     h(1:d1/2, 1:d2/2) = scale*gy(1:d1/2, 1:d2/2)
!     h(1:d1/2, d2/2+1:d2) = scale*gy(1:d1/2,od2/2-d2/2:od2/2)

!     return
!   end subroutine myjacb32
  
module fcheb_fftw3_jacb_low
  use mydef
  use four_fftw3_def
  use mod_const
  implicit none
  
  private 
  public :: myjacb_open, myjacb, myjacb_close, myjacb_pre, myjacb_v2
  
  integer, parameter :: nomp=2
  integer*8 :: cplanfwd, fplanfwd, cplanbwd(nomp), fplanbwd(nomp)
  real(dp) :: cscale, fscale
  
  integer, parameter :: n1=nx, d1=nx+1, n2=ny, d2=ny+2
  integer, parameter :: onx=3*nx/2, ony=3*ny/2
  integer, parameter :: on1=onx, od1=onx+1, on2=ony, od2=ony+2
  
  complex(dpc), dimension(:), allocatable :: iky
  real(dp) :: lx, ilx

  complex(dpc), dimension(od2/2, od1) :: fx, fy, gx, gy
  real(dp), dimension(od1, od2) :: rfx, rfy, rgx, rgy
  
  equivalence (fx, rfx), (fy, rfy), (gx, rgx), (gy, rgy)

  logical :: init=.FALSE.
  character(len=*), parameter :: modn='(fcheb_fftw3_jacb_low)'
contains
  
  !---------------------------------------------
  function fdy(x) result(y)
    complex(dpc), dimension(d2/2,d1) :: x,y
    integer :: i1
    
    do i1=1,d1
       y(:,i1) = iky*x(:,i1)
    end do
    
    return
  end function fdy
!---------------------------------------------

  function cdx(x) result(y)
    complex(dpc), dimension(d2/2,d1), intent(in) :: x
    complex(dpc), dimension(d2/2,d1) :: y
    integer :: i1,i2,i3
    real(dp) :: tmp(d1)
    complex(dpc) :: ctmp
    
    tmp(1) = ilx; tmp(2:d1-1) = two*ilx; tmp(d1) = ilx
    do i1=1,d2/2
       do i2=1, d1
          ctmp = zero
          do i3=i2+1,d1,2
             ctmp = ctmp+tmp(i2)*(i3-1)*x(i1,i3)
          end do
          y(i1,i2) = ctmp
       end do
    end do
    y(:, d1) = zero
    return
  end function cdx
!----------------------------------------------

  subroutine myjacb_open(delky0, lx0)
    implicit none
    real(dp), intent(in) :: lx0, delky0

    real(dp), dimension(:), allocatable :: in1, in2, out1
    complex(dpc), dimension(:), allocatable :: out2
    
    integer  :: cta(2), cineb(2), coneb(2)
    integer  :: fta(2), fineb(2), foneb(2),i1

    integer  :: ta, inemb, outemb


    ! initialize for fdy, cdx
    allocate(iky(d2/2))

    iky = ic*delky0*(/ (i1, i1=0,d2/2) /); iky(d2/2) = zero
    lx = lx0; ilx = one/lx
    
    allocate(in1(od1), out1(od1))

    ta = od1
    inemb = od1
    call dfftw_plan_many_r2r(cplanfwd, 1, ta, 1, in1,inemb, 1, 1, out1, inemb, 1, 1, FFTW_REDFT00, wisdom)
    
    do i1=1, nomp
       call dfftw_plan_many_r2r(cplanbwd(i1), 1, ta, 1, in1,inemb, 1, 1, out1, inemb, 1, 1, FFTW_REDFT00, wisdom)
    end do
    deallocate(in1,out1)

    cscale = one/real(on1,dp)
    
    write(fftout,*) " MYCFT2D initialization"//modn
    write(fftout,*) "Cosft fwd plan =   ",cplanfwd
    write(fftout,*) "Cosft bckwd plan = ",cplanbwd
    
    allocate(in2(od2), out2(od2/2))
    ta = on2
    inemb = od2;     outemb = od2/2
    
    call dfftw_plan_many_dft_r2c(fplanfwd,1, ta, 1, in2, inemb, 1, 1, out2, outemb, 1, 1, wisdom)
    
    do i1=1,nomp
       call dfftw_plan_many_dft_c2r(fplanbwd(i1),1, ta, 1, out2, outemb, 1, 1, in2, inemb, 1, 1,wisdom)
    end do
    
    fscale = one/real(on2,dp)

    deallocate(in2, out2)
    
    write(fftout,*) "Four. fwd plan =   ",fplanfwd
    write(fftout,*) "Four. bckwd plan = ",fplanbwd
    
    init = .TRUE.
    return
  end subroutine myjacb_open

  subroutine mycft2d_fwd(in, out)
    implicit none
    real(dp), dimension(od1,od2), intent(in) :: in
    real(dp), dimension(od2,od1)  :: med
    complex(dpc), dimension(od2/2,od1), intent(out) :: out
    integer :: i1

    if (.not.init) stop modn
    do i1=1,od2
       call dfftw_execute_r2r(cplanfwd, in(:,i1), med(i1,:))
    enddo
    
    med(:,1) = half*med(:,1)
    !    med(:,2:od1-1) = med(:,2:od1-1)
    med(:,od1) = half*med(:,od1)
    
    do i1=1,od1
       call dfftw_execute_dft_r2c(fplanfwd, med(:,i1), out(:,i1))
    enddo
    out = (cscale*fscale)*out
    
    return
  end subroutine mycft2d_fwd

  subroutine mycft2d_bwd(in, out, inum)
    implicit none
     complex(dpc), dimension(od2/2,od1), intent(in) :: in
     complex(dpc), dimension(od2/2,od1) :: tmp
     real(dp), dimension(od2,od1) :: med
     real(dp), dimension(od1,od2), intent(out) :: out
     integer :: inum, i1

    if (.not.init) stop modn
     tmp = in
     do i1=1, od1
        call dfftw_execute_dft_c2r(fplanbwd, tmp(:,i1), med(:,i1))
     enddo
     med(:,2:od1-1) = half*med(:,2:od1-1)
     do i1=1, od2
        call dfftw_execute_r2r(cplanbwd(inum), med(i1,:), out(:,i1))
     enddo

     return
   end subroutine mycft2d_bwd

   subroutine myjacb_close
     integer :: i1
     
    if (.not.init) stop modn
     deallocate(iky)

     call dfftw_destroy_plan(cplanfwd)
     call dfftw_destroy_plan(fplanfwd)
     do i1=1, nomp
        call dfftw_destroy_plan(cplanbwd(i1))
        call dfftw_destroy_plan(fplanbwd(i1))
     end do
   end subroutine myjacb_close
   
  function myjacb(f,g) result(h)
    complex(dpc), dimension(d2/2, d1), intent(in) :: f,g
    complex(dpc), dimension(d2/2, d1)  :: h

    if (nomp /= 4) stop modn
    if (.not.init) stop modn
    !$omp parallel sections shared(f,g)
    !$omp section
    fx(1:d2/2,1:d1) = cdx(f)
    fx(1:d2/2,d1+1:od1) = zero; fx(d2/2+1:od2/2,:) = zero
    call mycft2d_bwd(fx, rfx, 1)

    !$omp section
    fy(1:d2/2,1:d1) = fdy(f)
    fy(1:d2/2,d1+1:od1) = zero; fy(d2/2+1:od2/2,:) = zero
    call mycft2d_bwd(fy, rfy, 2)

    !$omp section
    gx(1:d2/2,1:d1) = cdx(g)
    gx(1:d2/2,d1+1:od1) = zero; gx(d2/2+1:od2/2,:) = zero
    call mycft2d_bwd(gx, rgx, 3)
    
    !$omp section
    gy(1:d2/2,1:d1) = fdy(g)
    gy(1:d2/2,d1+1:od1) = zero; gy(d2/2+1:od2/2,:) = zero
    call mycft2d_bwd(gy, rgy, 4)
    
    !$omp end parallel sections
    
    rfx = rfx*rgy-rfy*rgx; call mycft2d_fwd(rfx,fx)
    
    h(1:d2/2-1,1:d1) = fx(1:d2/2-1, 1:d1);     h(d2/2,1:d1) = fx(d2/2, 1:d1)
    return
  end function myjacb

  function myjacb_v2(g) result(h)
    complex(dpc), dimension(d2/2, d1), intent(in) :: g
    complex(dpc), dimension(d2/2, d1)  :: h

    if (nomp /= 2) stop modn
    if (.not.init) stop modn
    !$omp parallel sections shared(g)
    !$omp section
    gx(1:d2/2,1:d1) = cdx(g)
    gx(1:d2/2,d1+1:od1) = zero; gx(d2/2+1:od2/2,:) = zero
    call mycft2d_bwd(gx, rgx, 1)
    
    !$omp section
    gy(1:d2/2,1:d1) = fdy(g)
    gy(1:d2/2,d1+1:od1) = zero; gy(d2/2+1:od2/2,:) = zero
    call mycft2d_bwd(gy, rgy, 2)
    !$omp end parallel sections
    
    rgy = rfx*rgy-rfy*rgx; call mycft2d_fwd(rgy,gy)
    
    h(1:d2/2-1,1:d1) = gy(1:d2/2-1, 1:d1);     h(d2/2,1:d1) = gy(d2/2, 1:d1)
    return
  end function myjacb_v2

  subroutine myjacb_pre(f)
    complex(dpc), dimension(d2/2, d1), intent(in) :: f

    if (.not.init) stop modn
    !$omp parallel sections shared(f)
    !$omp section
    fx(1:d2/2,1:d1) = cdx(f)
    fx(1:d2/2,d1+1:od1) = zero; fx(d2/2+1:od2/2,:) = zero
    call mycft2d_bwd(fx, rfx, 1)

    !$omp section
    fy(1:d2/2,1:d1) = fdy(f)
    fy(1:d2/2,d1+1:od1) = zero; fy(d2/2+1:od2/2,:) = zero
    call mycft2d_bwd(fy, rfy, 2)
    !$omp end parallel sections
    return
  end subroutine myjacb_pre
end module fcheb_fftw3_jacb_low

module fcheb_fftw3_jacb
  use mydef
  use four_fftw3_def
  use mod_const
  implicit none
  
  private 
  public :: myjacb_open, myjacb, myjacb_close, myjacb_pre, myjacb_v2

  integer :: nomp
  integer*8 :: cplanfwd, fplanfwd
  integer*8, dimension(:), allocatable :: cplanbwd, fplanbwd
  real(dp) :: cscale, fscale
  
  integer, parameter :: n1=nx, d1=nx+1, n2=ny, d2=ny+2
  integer, parameter :: onx=3*nx/2, ony=3*ny/2
  integer, parameter :: on1=onx, od1=onx+1, on2=ony, od2=ony+2
  
  complex(dpc), dimension(:), allocatable :: iky
  real(dp) :: lx, ilx

  complex(dpc), dimension(od2/2, od1) :: fx, fy, gx, gy
  real(dp), dimension(od1, od2) :: rfx, rfy, rgx, rgy
  
  equivalence (fx, rfx), (fy, rfy), (gx, rgx), (gy, rgy)

  logical :: init=.FALSE.
  character(len=*), parameter :: modn='(fcheb_fftw3_jacb)'
contains
  
  !---------------------------------------------
  function fdy(x) result(y)
    complex(dpc), dimension(d2/2,d1) :: x,y
    integer :: i1
    
    do i1=1,d1
       y(:,i1) = iky*x(:,i1)
    end do
    
    return
  end function fdy
!---------------------------------------------

  function cdx(x) result(y)
    complex(dpc), dimension(d2/2,d1), intent(in) :: x
    complex(dpc), dimension(d2/2,d1) :: y
    integer :: i1,i2,i3
    real(dp) :: tmp(d1)
    complex(dpc) :: ctmp
    
    tmp(1) = ilx; tmp(2:d1-1) = two*ilx; tmp(d1) = ilx
    do i1=1,d2/2
       do i2=1, d1
          ctmp = zero
          do i3=i2+1,d1,2
             ctmp = ctmp+tmp(i2)*(i3-1)*x(i1,i3)
          end do
          y(i1,i2) = ctmp
       end do
    end do
    y(:, d1) = zero
    return
  end function cdx
!----------------------------------------------

  subroutine myjacb_open(delky0, lx0, nomp0)
    implicit none
    real(dp), intent(in) :: lx0, delky0

    real(dp), dimension(:,:), allocatable :: in, med
    complex(dpc), dimension(:,:), allocatable :: out
    
    integer  :: cta(2), cineb(2), coneb(2)
    integer  :: fta(2), fineb(2), foneb(2),i1

    integer, optional :: nomp0

    ! initialize for fdy, cdx
    allocate(iky(d2/2))

    nomp = 4
    if (present(nomp0)) then
       if ((nomp0 == 2) .or. (nomp0 == 4)) nomp = nomp0
    end if

    allocate(cplanbwd(nomp), fplanbwd(nomp))

    iky = ic*delky0*(/ (i1, i1=0,d2/2) /); iky(d2/2) = zero
    lx = lx0; ilx = one/lx
    
    allocate(in(od1,od2), med(od2,od1), out(od2/2,od1))
    
    cta(1) = od1; cta(2)=1
    cineb(1) = od1; cineb(2)=1
    coneb(1) = od1; coneb(2)=1
    
    call dfftw_plan_many_r2r(cplanfwd, 2, cta, on2, in,cineb, 1, od1, med, coneb, od2,1, FFTW_REDFT00, wisdom)
    
    do i1=1, nomp
       call dfftw_plan_many_r2r(cplanbwd(i1), 2, cta, on2, med,coneb, od2, 1, in, cineb, 1, od1, FFTW_REDFT00, wisdom)
    end do
    
    cscale = one/real(on1,dp)
    
    write(fftout,*) " MYCFT2D initialization"//modn
    write(fftout,*) "Cosft fwd plan =   ",cplanfwd
    write(fftout,*) "Cosft bckwd plan = ",cplanbwd
    
    fta(1) = on2; fta(2)=1
    fineb(1) = od2; fineb(2)=1
    foneb(1) = od2/2; foneb(2)=1
    
    call dfftw_plan_many_dft_r2c(fplanfwd,2, fta, od1, med, fineb, 1, od2, out, foneb, 1, od2/2, wisdom)
    
    do i1=1,nomp
       call dfftw_plan_many_dft_c2r(fplanbwd(i1),2, fta, od1, out, foneb, 1, od2/2, med, fineb, 1, od2,wisdom)
    end do
       
    fscale = one/real(on2,dp)

    deallocate(in, med, out)
    
    write(fftout,*) "Four fwd plan =   ",fplanfwd
    write(fftout,*) "Four bckwd plan = ",fplanbwd
    
    init = .TRUE.
    return
  end subroutine myjacb_open

  subroutine mycft2d_fwd(in, out)
    implicit none
    real(dp), dimension(od1,od2), intent(in) :: in
    real(dp), dimension(od2,od1)  :: med
    complex(dpc), dimension(od2/2,od1), intent(out) :: out

    if (.not.init) call error_exit('mycft2d_fwd'//modn)
    call dfftw_execute_r2r(cplanfwd, in, med)

    med(:,1) = half*med(:,1)
    med(:,2:od1-1) = med(:,2:od1-1)
    med(:,od1) = half**med(:,od1)

    call dfftw_execute_dft_r2c(fplanfwd, med, out)
    out = cscale*fscale*out
    
    return
  end subroutine mycft2d_fwd
  
  subroutine mycft2d_bwd(in, out, inum)
    implicit none
     complex(dpc), dimension(od2/2,od1), intent(in) :: in
     complex(dpc), dimension(od2/2,od1) :: tmp
     real(dp), dimension(od2,od1) :: med
     real(dp), dimension(od1,od2), intent(out) :: out
     integer :: inum 

    if (.not.init) stop modn
     tmp = in
     call dfftw_execute_dft_c2r(fplanbwd(inum), tmp, med)

     med(:,2:od1-1) = half*med(:,2:od1-1)
     call dfftw_execute_r2r(cplanbwd(inum), med, out)
     return
   end subroutine mycft2d_bwd
   
   subroutine myjacb_close
     integer :: i1
     
    if (.not.init) call error_exit('myjacb_close'//modn)
     deallocate(iky)

     call dfftw_destroy_plan(cplanfwd)
     call dfftw_destroy_plan(fplanfwd)
     do i1=1, nomp
        call dfftw_destroy_plan(cplanbwd(i1))
        call dfftw_destroy_plan(fplanbwd(i1))
     end do

     deallocate(cplanbwd)
     deallocate(fplanbwd)
     return
   end subroutine myjacb_close
   
  function myjacb(f,g) result(h)
    complex(dpc), dimension(d2/2, d1), intent(in) :: f,g
    complex(dpc), dimension(d2/2, d1)  :: h

    if (nomp /= 4) call error_exit('myjacb omp'//modn)
    if (.not.init) call error_exit('myjacb'//modn)
    !$omp parallel sections shared(f,g)
    !$omp section
    fx(1:d2/2,1:d1) = cdx(f)
    fx(1:d2/2,d1+1:od1) = zero; fx(d2/2+1:od2/2,:) = zero
    call mycft2d_bwd(fx, rfx, 1)

    !$omp section
    fy(1:d2/2,1:d1) = fdy(f)
    fy(1:d2/2,d1+1:od1) = zero; fy(d2/2+1:od2/2,:) = zero
    call mycft2d_bwd(fy, rfy, 2)

    !$omp section
    gx(1:d2/2,1:d1) = cdx(g)
    gx(1:d2/2,d1+1:od1) = zero; gx(d2/2+1:od2/2,:) = zero
    call mycft2d_bwd(gx, rgx, 3)
    
    !$omp section
    gy(1:d2/2,1:d1) = fdy(g)
    gy(1:d2/2,d1+1:od1) = zero; gy(d2/2+1:od2/2,:) = zero
    call mycft2d_bwd(gy, rgy, 4)
    
    !$omp end parallel sections
    
    rfx = rfx*rgy-rfy*rgx; call mycft2d_fwd(rfx,fx)
    
    h(1:d2/2-1,1:d1) = fx(1:d2/2-1, 1:d1);     h(d2/2,1:d1) = fx(d2/2, 1:d1)
    return
  end function myjacb

  function myjacb_v2(g) result(h)
    complex(dpc), dimension(d2/2, d1), intent(in) :: g
    complex(dpc), dimension(d2/2, d1)  :: h

    if (nomp /= 2) call error_exit('myjacb_v2 omp'//modn)
    if (.not.init) call error_exit('myjacb_v2 '//modn)
    !$omp parallel sections shared(g)
    !$omp section
    gx(1:d2/2,1:d1) = cdx(g)
    gx(1:d2/2,d1+1:od1) = zero; gx(d2/2+1:od2/2,:) = zero
    call mycft2d_bwd(gx, rgx, 1)
    
    !$omp section
    gy(1:d2/2,1:d1) = fdy(g)
    gy(1:d2/2,d1+1:od1) = zero; gy(d2/2+1:od2/2,:) = zero
    call mycft2d_bwd(gy, rgy, 2)
    !$omp end parallel sections
    
    rgy = rfx*rgy-rfy*rgx; call mycft2d_fwd(rgy,gy)
    
    h(1:d2/2-1,1:d1) = gy(1:d2/2-1, 1:d1);     h(d2/2,1:d1) = gy(d2/2, 1:d1)
    return
  end function myjacb_v2

  subroutine myjacb_pre(f)
    complex(dpc), dimension(d2/2, d1), intent(in) :: f

    if (.not.init) stop modn
    !$omp parallel sections shared(f)
    !$omp section
    fx(1:d2/2,1:d1) = cdx(f)
    fx(1:d2/2,d1+1:od1) = zero; fx(d2/2+1:od2/2,:) = zero
    call mycft2d_bwd(fx, rfx, 1)

    !$omp section
    fy(1:d2/2,1:d1) = fdy(f)
    fy(1:d2/2,d1+1:od1) = zero; fy(d2/2+1:od2/2,:) = zero
    call mycft2d_bwd(fy, rfy, 2)
    !$omp end parallel sections
    return
  end subroutine myjacb_pre
end module fcheb_fftw3_jacb

