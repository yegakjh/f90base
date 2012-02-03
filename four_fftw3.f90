module four_fftw3_def
  use mydef
  implicit none
  include "fftw3.f"

  type pfftw3
     integer*8 :: fwd, bwd
     integer :: n
  end type pfftw3
  integer, protected :: fftout=6
  integer, parameter :: wisdom=FFTW_EXHAUSTIVE
  
contains
  subroutine set_fft_output(in) 
    integer, intent(in) :: in
    fftout = in
    return
  end subroutine set_fft_output

  subroutine error_exit(msg)
    character(len=*) :: msg
    
    write(*,*) msg
    stop
  end subroutine error_exit
end module four_fftw3_def

module four_fftw3_1d_c2c
  use four_fftw3_def
  implicit none
  
  private 
  public :: myfft1d_open, myfft1d_fwd, myfft1d_bwd, myfft1d_close
  
  type(pfftw3) :: cplan
  real(dp) :: scale
  
  integer :: n1, d1

  logical :: init=.FALSE.
  character(len=*), parameter :: modn='(four_fftw3_1d_c2c)'

contains
  subroutine myfft1d_open(nn)
    implicit none
    integer, intent(in) :: nn
    complex(dpc), dimension(nn)  :: in, out
    
    integer  :: ta, inemb, outemb

    n1 = nn; d1 = nn
    ta = n1
    inemb = d1; outemb = d1
    
    call dfftw_plan_many_dft(cplan%fwd, 1, ta, 1, in, inemb, 1, 1, out, outemb, 1, 1,FFTW_FORWARD, wisdom)
    call dfftw_plan_many_dft(cplan%bwd, 1, ta, 1, out, outemb, 1, 1, in, inemb, 1, 1,FFTW_BACKWARD, wisdom)
    
    scale = one/real(n1,dp)
    
    write(fftout,*) " MYFFT1D initialiation "//modn
    write(fftout,*) " ta, inemb, outemb = ", ta, inemb, outemb
    write(fftout,*) " fwd plan =   ",cplan%fwd
    write(fftout,*) " bckwd plan = ",cplan%bwd
    
    init=.TRUE.
    return
  end subroutine myfft1d_open
  
  subroutine myfft1d_fwd(in, out)
    implicit none
    complex(dpc), dimension(d1), intent(in) :: in
    complex(dpc), dimension(d1), intent(out) :: out
    
    if (.not.init) stop modn
    call dfftw_execute_dft(cplan%fwd, in, out)
    out = scale*out
    return
  end subroutine myfft1d_fwd
  
  subroutine myfft1d_bwd(in, out)
    implicit none
    complex(dpc), dimension(d1), intent(in) :: in
    complex(dpc), dimension(d1)  :: tmp
    complex(dpc), dimension(d1), intent(out) :: out
    
    if (.not.init) stop modn
    tmp = in
    call dfftw_execute_dft(cplan%bwd, tmp, out)
    return
  end subroutine myfft1d_bwd
  
  subroutine myfft1d_close
    if (init) then
       call dfftw_destroy_plan(cplan%fwd)
       call dfftw_destroy_plan(cplan%bwd)
    end if
    return
  end subroutine myfft1d_close
  
end module four_fftw3_1d_c2c

! Multiple 1D FFT transform on the first rank of 2D array
module four_fftw3_1d_2a_c2c
  use four_fftw3_def
  implicit none
  
  private 
  public :: myfft1d2a_open, myfft1d2a_fwd, myfft1d2a_bwd, myfft1d2a_close
  
  type(pfftw3) :: cplan
  real(dp) :: scale
  
  integer :: n1, n2, d1,d2
  
  logical :: init=.FALSE.
  character(len=*), parameter :: modn='(four_fftw3_1d_2a_c2c)'
contains
  
  subroutine myfft1d2a_open(nx, ny)
    implicit none
    integer, intent(in) :: nx, ny
    complex(dpc), dimension(nx, ny)  :: out, in
    
    integer, dimension(2) :: ta, inemb, outemb

    ! 1 : transform on 1st dimension
    ! 2 : transform on 2nd dimension
    
    n1 = nx; n2 = ny
    d1 = nx; d2 = ny
    
    ta(1) = n1; ta(2) = 1
    inemb(1) = d1; inemb(2) = 1
    outemb(1) = d1; outemb(2) = 1
    
    call dfftw_plan_many_dft(cplan%fwd,2, ta, n2, in, inemb, 1, d1, out, outemb, 1, d1, wisdom)
    call dfftw_plan_many_dft(cplan%bwd,2, ta, n2, out, outemb, 1, d1, in, inemb, 1, d1, wisdom)
    
    scale = one/real(n1,dp)
    
    write(fftout,*) " MYFFT1D2A"//modn
    write(fftout,*) "fwd plan =   ",cplan%fwd
    write(fftout,*) "bckwd plan = ",cplan%bwd
    
    init=.TRUE.
    return
  end subroutine myfft1d2a_open
  
  subroutine myfft1d2a_fwd(in, out)
    implicit none
    complex(dpc), dimension(d1, d2), intent(in) :: in
    complex(dpc), dimension(d1, d2), intent(out) :: out

    if (.not.init) stop modn    
    call dfftw_execute_dft(cplan%fwd, in, out)
    out = scale*out
    return
  end subroutine myfft1d2a_fwd
  
  subroutine myfft1d2a_bwd(in, out)
    implicit none
    complex(dpc), dimension(d1,d2), intent(in) :: in
    complex(dpc), dimension(d1,d2) :: tmp
    complex(dpc), dimension(d1,d2), intent(out) :: out
    
    if (.not.init) stop modn
    tmp = in
    call dfftw_execute_dft(cplan%bwd, tmp, out)
    return
  end subroutine myfft1d2a_bwd
  
  subroutine myfft1d2a_close
    if (init) then
       call dfftw_destroy_plan(cplan%fwd)
       call dfftw_destroy_plan(cplan%bwd)
    end if
  end subroutine myfft1d2a_close
  
end module four_fftw3_1d_2a_c2c

! Multiple 1D complex FFT transform on the second rank of 2D array
module four_fftw3_1d_2a_c2c_r2
  use four_fftw3_def
  implicit none
  
  private 
  public :: myfft1d2a_open, myfft1d2a_fwd, myfft1d2a_bwd, myfft1d2a_close
  
  type(pfftw3) :: cplan
  real(dp) :: scale
  
  integer :: n1, n2, d1,d2
  
  logical :: init=.FALSE.
  character(len=*), parameter :: modn='(four_fftw3_1d_2a_c2c_r2)'
contains
  
  subroutine myfft1d2a_open(nx, ny)
    implicit none
    integer, intent(in) :: nx, ny
    complex(dpc), dimension(nx, ny)  :: out, in
    
    integer, dimension(2) :: ta, inemb, outemb

    ! 1 : transform on 1st dimension
    ! 2 : transform on 2nd dimension
    
    n1 = nx; n2 = ny
    d1 = nx; d2 = ny
    
    ta(1) = n2; ta(2) = 1
    inemb(1) = d2; inemb(2) = 1
    outemb(1) = d2; outemb(2) = 1
    
    call dfftw_plan_many_dft(cplan%fwd,2,ta,n1,in,inemb,d1,1,out,outemb,d1,1,wisdom)
    call dfftw_plan_many_dft(cplan%bwd,2,ta,n1,out,outemb,d1,1,in,inemb,d1,1,wisdom)
    
    scale = one/real(n2,dp)
    
    write(fftout,*) " MYFFT1D2A initialiation"//modn
    write(fftout,*) "fwd plan =   ",cplan%fwd
    write(fftout,*) "bckwd plan = ",cplan%bwd

    init=.TRUE.
    return
  end subroutine myfft1d2a_open
  
  subroutine myfft1d2a_fwd(in, out)
    implicit none
    complex(dpc), dimension(d1, d2), intent(in) :: in
    complex(dpc), dimension(d1, d2), intent(out) :: out
    
    if (.not.init) stop modn
    call dfftw_execute_dft(cplan%fwd, in, out)
    out = scale*out
    
  end subroutine myfft1d2a_fwd
  
  subroutine myfft1d2a_bwd(in, out)
    implicit none
    complex(dpc), dimension(d1,d2), intent(in) :: in
    complex(dpc), dimension(d1,d2) :: tmp
    complex(dpc), dimension(d1,d2), intent(out) :: out
    
    if (.not.init) stop modn
    tmp = in
    call dfftw_execute_dft(cplan%bwd, tmp, out)
  end subroutine myfft1d2a_bwd
  
  subroutine myfft1d2a_close
    if (.not.init) stop modn
    call dfftw_destroy_plan(cplan%fwd)
    call dfftw_destroy_plan(cplan%bwd)
  end subroutine myfft1d2a_close
  
end module four_fftw3_1d_2a_c2c_r2

! Multiple 1D real-complex FFT transform on the first rank of 2D array
module four_fftw3_1d_2a_r2c
  use four_fftw3_def
  implicit none
  
  private 
  public :: myfft1d2a_open, myfft1d2a_fwd, myfft1d2a_bwd, myfft1d2a_close
  
  type(pfftw3) :: cplan
  real(dp) :: scale
  
  integer :: n1, n2, d1,d2
  
  logical :: init=.FALSE.
  character(len=*), parameter :: modn='(four_fftw3_1d_2a_r2c)'
contains
  
  subroutine myfft1d2a_open(nx, ny)
    implicit none
    integer, intent(in) :: nx, ny
    real(dp), dimension(nx+2,ny)  :: in
    complex(dpc), dimension(nx/2+1, ny)  :: out

    integer, dimension(2) :: ta, inemb, outemb

    ! 1 : transform on 1st dimension
    ! 2 : transform on 2nd dimension
    
    n1 = nx; n2 = ny
    d1 = nx+2; d2 = ny
    
    ta(1) = n1; ta(2) = 1
    inemb(1) = d1; inemb(2) = 1
    outemb(1) = d1/2; outemb(2) = 1
    
    call dfftw_plan_many_dft_r2c(cplan%fwd,2, ta, n2, in, inemb, 1, d1, out, outemb, 1, d1/2, wisdom)
    call dfftw_plan_many_dft_c2r(cplan%bwd,2, ta, n2, out, outemb, 1, d1/2, in, inemb, 1, d1,wisdom)
    
    scale = one/real(n1,dp)
    
    write(fftout,*) " MYFFT1D2A initialiation"//modn
    write(fftout,*) "fwd plan =   ",cplan%fwd
    write(fftout,*) "bckwd plan = ",cplan%bwd
    init=.TRUE.
    return
  end subroutine myfft1d2a_open
  
  subroutine myfft1d2a_fwd(in, out)
    implicit none
    real(dp), dimension(d1,d2), intent(in) :: in
    complex(dpc), dimension(d1/2, d2), intent(out) :: out

    if (.not.init) stop modn
    call dfftw_execute_dft_r2c(cplan%fwd, in, out)
    out = scale*out
    
  end subroutine myfft1d2a_fwd
  
  subroutine myfft1d2a_bwd(in, out)
    implicit none
    complex(dpc), dimension(d1/2,d2), intent(in) :: in
    complex(dpc), dimension(d1/2,d2) :: tmp
    real(dp), dimension(d1,d2), intent(out) :: out

    if (.not.init) stop modn
    tmp = in
    call dfftw_execute_dft_c2r(cplan%bwd, tmp, out)
  end subroutine myfft1d2a_bwd
  
  subroutine myfft1d2a_close
    if (init) then
       call dfftw_destroy_plan(cplan%fwd)
       call dfftw_destroy_plan(cplan%bwd)
    end if
  end subroutine myfft1d2a_close
  
end module four_fftw3_1d_2a_r2c

! Multiple 2D FFT transform on the 2D array
module four_fftw3_2d_r2c
  use four_fftw3_def
  implicit none
  
  private 
  public :: myfft2d_open, myfft2d_fwd, myfft2d_bwd, myfft2d_close
  
  type(pfftw3) :: cplan
  real(dp) :: scale
  
  integer :: n1, n2, d1,d2

  public :: myjacb_open, myjacb, myjacb32,myjacb_close
  
  type(pfftw3) :: c1p
  real(dp) :: jscale
  
  logical :: init=.FALSE., jacb_init=.FALSE.
  character(len=*), parameter :: modn='(four_fftw3_2d_r2c)'
  character(len=*), parameter :: modn_jacb='(four_fftw3_2d_jacb)'

  ! jn, jd : original size
  ! jno, jdo : JACOBIAN
  integer :: jn1,jn2,jd1,jd2, jno1, jno2, jdo1, jdo2

contains

  subroutine myfft2d_open(nx, ny)
    implicit none
    integer, intent(in) :: nx, ny
    real(dp), dimension(nx+2,ny)  :: in
    complex(dpc), dimension(nx/2+1, ny)  :: out
    
    integer, dimension(2) :: ta, inemb, outemb
    
    ! 1 : transform on 1st dimension
    ! 2 : transform on 2nd dimension
    
    n1 = nx; n2 = ny
    d1 = nx+2; d2 = ny
    
    ta(1) = n1; ta(2) = n2
    inemb(1) = d1; inemb(2) = d2
    outemb(1) = d1/2; outemb(2) = d2
    
    call dfftw_plan_many_dft_r2c(cplan%fwd,2, ta, 1, in, inemb, 1, 1, out, outemb, 1, 1, wisdom)
    call dfftw_plan_many_dft_c2r(cplan%bwd,2, ta, 1, out, outemb, 1, 1, in, inemb, 1, 1,wisdom)
    
    scale = one/real(n1*n2,dp)
    
    write(fftout,*) "MYFFT2D initialiation"//modn
    write(fftout,*) "fwd plan =   ",cplan%fwd
    write(fftout,*) "bckwd plan = ",cplan%bwd
    
    init=.TRUE.
    return
  end subroutine myfft2d_open
  
  subroutine myfft2d_fwd(in, out)
    implicit none
    real(dp), dimension(d1,d2), intent(in) :: in
    complex(dpc), dimension(d1/2, d2), intent(out) :: out
    
    if (.not.init) stop modn
    call dfftw_execute_dft_r2c(cplan%fwd, in, out)
    out = scale*out
    
  end subroutine myfft2d_fwd
  
  subroutine myfft2d_bwd(in, out)
    implicit none
    complex(dpc), dimension(d1/2,d2), intent(in) :: in
    complex(dpc), dimension(d1/2,d2) :: tmp
    real(dp), dimension(d1,d2), intent(out) :: out

    if (.not.init) stop modn    
    tmp = in
    call dfftw_execute_dft_c2r(cplan%bwd, tmp, out)
  end subroutine myfft2d_bwd
  
  subroutine myfft2d_close
    if (init) then
       call dfftw_destroy_plan(cplan%fwd)
       call dfftw_destroy_plan(cplan%bwd)
       init=.true.
    end if
  end subroutine myfft2d_close

  subroutine myjacb_open(nx, ny, nxo, nyo)
    implicit none
    integer, intent(in) :: nx, ny, nxo, nyo
    complex(dpc), dimension(nxo/2+1,nyo)  :: in
    
    integer, dimension(2) :: ta, inemb, outemb
    
    ! 1 : transform on 1st dimension
    ! 2 : transform on 2nd dimension
    
    jn1 = nx; jn2 = ny
    jd1 = nx+2; jd2 = ny

    jno1 = nxo ; jno2 = nyo
    jdo1 = nxo+2; jdo2 = nyo

    write(fftout,*) " MYJACB initialiation (four_fftw3_2d_jacb)"
    write(fftout,*) " ( nx,  ny)=", nx, ny
    write(fftout,*) " (nxo, nyo)=", nxo, nyo

    if ( jn1 < jno1) then
       ta(1) = jno1; ta(2) = jno2
       inemb(1) = jdo1; inemb(2) = jdo2
       outemb(1) = jdo1/2; outemb(2) = jdo2
       
       jscale = one/real(jno1*jno2,dp)
       
       write(fftout,*) 'ta = ', ta
       write(fftout,*) 'inemb = ', inemb
       write(fftout,*) 'outemb = ', outemb
       
    else
       ta(1) = jn1; ta(2) = jn2
       inemb(1) = jd1; inemb(2) = jd2
       outemb(1) = jd1/2; outemb(2) = jd2

       jscale = one/real(jn1*jn2,dp)
       
       write(fftout,*) 'ta = ', ta
       write(fftout,*) 'inemb = ', inemb
       write(fftout,*) 'outemb = ', outemb
       
    end if

       call dfftw_plan_many_dft_r2c(c1p%fwd,2, ta, 1, in, inemb, 1, 1, in, outemb, 1,1, wisdom)
       call dfftw_plan_many_dft_c2r(c1p%bwd,2, ta, 1, in, outemb, 1, 1, in, inemb, 1, 1,wisdom)

       write(fftout,*) "fwd plan =   ",c1p%fwd
       write(fftout,*) "bckwd plan = ",c1p%bwd
       jacb_init = .TRUE.
    return
  end subroutine myjacb_open
  
  subroutine myjacb(derx, dery, f, g,h)
    implicit none
    complex(dpc), intent(in) :: derx(jd1/2), dery(jd2)
    complex(dpc), dimension(jd1/2, jd2), intent(in) :: f,g
    complex(dpc), dimension(jd1/2, jd2), intent(out) :: h

    complex(dpc), dimension(jd1/2, jd2) :: fx, gx, fy
    real(dp), dimension(jd1, jd2) :: rfx, rgx, rfy, rgy
    integer :: i1
    integer :: omp_get_num_threads, omp_get_thread_num


    if (.not.jacb_init) stop modn_jacb

    if (jd1 < jdo1 .or. jd2 < jdo2) stop 'This subroutine is to use 2/3 method'

    !$omp parallel sections
    !$omp section
    do i1=1, jd2
       fx(:,i1) = derx*f(:,i1)
    end do
    call dfftw_execute_dft_c2r(c1p%bwd, fx, rfx)    

    !$omp section
    do i1=1, jd2
       fy(:,i1) = dery(i1)*f(:,i1)
    end do
    call dfftw_execute_dft_c2r(c1p%bwd, fy, rfy)    

    !$omp section
    do i1=1, jd2
       gx(:,i1) = derx*g(:,i1)
    end do
    call dfftw_execute_dft_c2r(c1p%bwd, gx, rgx)    

    !$omp section
    do i1=1, jd2
       h(:,i1) = dery(i1)*g(:,i1)
    end do
    call dfftw_execute_dft_c2r(c1p%bwd, h, rgy)    
    !$omp end parallel sections

    !h = (real(fx,dp)*real(h,dp)-real(fy,dp)*real(gx,dp)) &
    !+ic*(dimag(fx)*dimag(h)-dimag(fy)*dimag(gx))
    rgy = rfx*rgy-rfy*rgx

    call dfftw_execute_dft_r2c(c1p%fwd, rgy, h)
    h = jscale*h

    return
  end subroutine myjacb

  subroutine myjacb32(derx, dery, f, g,h)
    implicit none
    complex(dpc), intent(in) :: derx(jd1/2), dery(jd2)
    complex(dpc), dimension(jd1/2, jd2), intent(in) :: f,g
    complex(dpc), dimension(jd1/2, jd2), intent(out) :: h
    complex(dpc), dimension(jdo1/2, jdo2) :: fx, gx, fy,gy
    integer :: i1
    integer :: omp_get_num_threads, omp_get_thread_num

    if (.not.jacb_init) stop modn_jacb
    if (jd1 > jdo1 .or. jd2  > jdo2) stop 'This subroutine is to use 3/2 method'

    !$omp parallel sections 
    !$omp section
    do i1=1, jd2/2-1
       fx(1:jd1/2-1,i1) = derx(1:jd1/2-1)*f(1:jd1/2-1,i1)
       fx(1:jd1/2-1,jdo2+1-i1) = derx(1:jd1/2-1)*f(1:jd1/2-1,jd2+1-i1)
       fx(jd1/2:jdo1/2,i1) = zero
       fx(jd1/2:jdo1/2,jdo2+1-i1) = zero
    end do
    fx(:,jd2/2:jdo2-jd2/2+1) =zero

    call dfftw_execute_dft_c2r(c1p%bwd, fx, fx)    

    !$omp section

    do i1=1, jd2/2-1
       fy(1:jd1/2-1,i1) = dery(i1)*f(1:jd1/2-1,i1)
       fy(1:jd1/2-1,jdo2+1-i1) = dery(jd2+1-i1)*f(1:jd1/2-1,jd2+1-i1)
       fy(jd1/2:jdo1/2,i1) = zero
       fy(jd1/2:jdo1/2,jdo2+1-i1) = zero
    end do
    fy(:,jd2/2:jdo2-jd2/2+1) =zero

    call dfftw_execute_dft_c2r(c1p%bwd, fy, fy)    
    
    !$omp section 

    do i1=1, jd2/2-1
       gx(1:jd1/2-1,i1) = derx(1:jd1/2-1)*g(1:jd1/2-1,i1)
       gx(1:jd1/2-1,jdo2+1-i1) = derx(1:jd1/2-1)*g(1:jd1/2-1,jd2+1-i1)
       gx(jd1/2:jdo1/2,i1) = zero
       gx(jd1/2:jdo1/2,jdo2+1-i1) = zero
    end do
    gx(:,jd2/2:jdo2-jd2/2+1) =zero

    call dfftw_execute_dft_c2r(c1p%bwd, gx, gx)    

    !$omp section
    do i1=1, jd2/2-1
       gy(1:jd1/2-1,i1) = dery(i1)*g(1:jd1/2-1,i1)
       gy(1:jd1/2-1,jdo2+1-i1) = dery(jd2+1-i1)*g(1:jd1/2-1,jd2+1-i1)
       gy(jd1/2:jdo1/2,i1) = zero
       gy(jd1/2:jdo1/2,jdo2+1-i1) = zero
    end do
    gy(:,jd2/2:jdo2-jd2/2+1) =zero

    call dfftw_execute_dft_c2r(c1p%bwd, gy, gy)    

    !$omp end parallel sections

    gy = (real(fx,dp)*real(gy,dp)-real(fy,dp)*real(gx,dp)) &
         +ic*(dimag(fx)*dimag(gy)-dimag(fy)*dimag(gx))
    
    call dfftw_execute_dft_r2c(c1p%fwd, gy, gy)
    h(1:jd1/2, 1:jd2/2) = jscale*gy(1:jd1/2, 1:jd2/2)
    h(1:jd1/2, jd2/2+1:jd2) = jscale*gy(1:jd1/2,jdo2/2-jd2/2:jdo2/2)

    return
  end subroutine myjacb32
  
  subroutine myjacb_close
    if (.not.jacb_init) stop modn_jacb
    call dfftw_destroy_plan(c1p%fwd)
    call dfftw_destroy_plan(c1p%bwd)
  end subroutine myjacb_close
end module four_fftw3_2d_r2c

module four_fftw3_1d_r2c
  use four_fftw3_def
  implicit none
  
  private 
  public :: myfft1d_open, myfft1d_fwd, myfft1d_bwd, myfft1d_close
  
  type(pfftw3) :: cplan
  real(dp) :: scale
  
  integer :: n1, d1

  logical :: init=.FALSE.
  character(len=*), parameter :: modn='(four_fftw3_1d_r2c)'
contains
  subroutine myfft1d_open(nn)
    implicit none
    integer, intent(in) :: nn
    real(dp), dimension(nn+2)  :: in
    complex(dpc), dimension(nn/2+1)  :: out
    
    integer  :: ta, inemb, outemb
    
    n1 = nn;     d1 = nn+2
    ta = n1
    inemb = d1;     outemb = d1/2
    
    call dfftw_plan_many_dft_r2c(cplan%fwd,1, ta, 1, in, inemb, 1, 1, out, outemb, 1, 1, wisdom)
    call dfftw_plan_many_dft_c2r(cplan%bwd,1, ta, 1, out, outemb, 1, 1, in, inemb, 1, 1,wisdom)
    
    scale = one/real(n1,dp)
    
    write(fftout,*) " MYFFT1D initialiation"//modn
    write(fftout,*) "fwd plan =   ",cplan%fwd
    write(fftout,*) "bckwd plan = ",cplan%bwd
    init = .TRUE.
    return
  end subroutine myfft1d_open
  
  subroutine myfft1d_fwd(in, out)
    implicit none
    real(dp), dimension(d1), intent(in) :: in
    complex(dpc), dimension(d1/2), intent(out) :: out
    
    if (.not.init) stop modn
    call dfftw_execute_dft_r2c(cplan%fwd, in, out)
    out = scale*out
    return
  end subroutine myfft1d_fwd
  
  subroutine myfft1d_bwd(in, out)
    implicit none
    complex(dpc), dimension(d1/2), intent(in) :: in
    complex(dpc), dimension(d1/2) :: tmp
    real(dp), dimension(d1), intent(out) :: out
    
    if (.not.init) stop modn
    tmp = in
    call dfftw_execute_dft_c2r(cplan%bwd, tmp, out)
    return
  end subroutine myfft1d_bwd
  
  subroutine myfft1d_close
    if (init) then
       call dfftw_destroy_plan(cplan%fwd)
       call dfftw_destroy_plan(cplan%bwd)
       init=.false.
    end if
    return
  end subroutine myfft1d_close
  
end module four_fftw3_1d_r2c

module cheb_fftw3_1d
  use four_fftw3_def
  implicit none
  
  private 
  interface mycft1d_fwd
     module procedure mycft1d_fwd_d1, mycft1d_fwd_d2_inplace
  end interface
  
  interface mycft1d_bwd
     module procedure mycft1d_bwd_d1, mycft1d_bwd_d2_inplace
  end interface

  public :: mycft1d_open, mycft1d_fwd, mycft1d_bwd, mycft1d_close, set_fft_output
  
  type(pfftw3) :: cplan
  real(dp) :: scale
  
  integer :: n1, d1

  logical :: init=.FALSE.
  character(len=*), parameter :: modn='(cheb_fftw3_1d)'
contains
  subroutine mycft1d_open(nn)
    implicit none
    integer, intent(in) :: nn
    real(dp), dimension(nn+1)  :: in,out
    
    integer  :: ta, inemb

    n1 = nn;     d1 = nn+1;   
    ta = d1
    inemb = d1

    call dfftw_plan_many_r2r(cplan%fwd, 1, ta, 1, in,inemb, 1, 1, out, inemb, 1, 1, FFTW_REDFT00, wisdom)
    
    scale = one/real(n1,dp)
    
    write(fftout,*) " MYCFT1D initialiation"//modn
    write(fftout,*) "fwd plan =   ",cplan%fwd
    init = .TRUE.
    
    return
  end subroutine mycft1d_open
  
  subroutine mycft1d_fwd_d1(in, out)
    implicit none
    real(dp), dimension(d1), intent(in) :: in
    real(dp), dimension(d1), intent(out) :: out
    
    if (.not.init) stop modn
    call dfftw_execute_r2r(cplan%fwd, in, out)
    out = scale*out
    out(1) = half*out(1)
    out(d1) = half*out(d1)
    return
  end subroutine mycft1d_fwd_d1
  
  subroutine mycft1d_bwd_d1(in, out)
    implicit none
    real(dp), dimension(d1), intent(in) :: in
    real(dp), dimension(d1) :: tmp
    real(dp), dimension(d1), intent(out) :: out
    
    if (.not.init) stop modn
    tmp(1) = in(1)
    tmp(2:d1-1) = half*in(2:d1-1)
    tmp(d1) = in(d1)
    call dfftw_execute_r2r(cplan%fwd, tmp, out)
    return
  end subroutine mycft1d_bwd_d1
  
  subroutine mycft1d_bwd_d2_inplace(nd2,out)
    integer, intent(in) :: nd2
    real(dp), dimension(d1, nd2), intent(inout) :: out
    integer :: i1

    do i1=1, nd2
       call mycft1d_bwd_d1(out(:,i1),out(:,i1))
    end do
    return
  end subroutine mycft1d_bwd_d2_inplace
  
  subroutine mycft1d_fwd_d2_inplace(nd2,out)
    integer, intent(in) :: nd2
    real(dp), dimension(d1, nd2), intent(inout) :: out
    integer :: i1
    
    do i1=1, nd2
       call mycft1d_fwd_d1(out(:,i1),out(:,i1))
    end do
    return
  end subroutine mycft1d_fwd_d2_inplace
  
  subroutine mycft1d_close
    if (init) then
       call dfftw_destroy_plan(cplan%fwd)
       call dfftw_destroy_plan(cplan%bwd)
       init = .false.
    end if
  end subroutine mycft1d_close
end module cheb_fftw3_1d

module cheb_fftw3_1d_2a
  use four_fftw3_def
  implicit none
  
  private 
  public :: mycft1d2a_open, mycft1d2a_fwd, mycft1d2a_bwd, mycft1d2a_close
  
  type(pfftw3) :: cplan
  real(dp) :: scale
  
  integer :: n1, d1, n2, d2

  logical :: init=.FALSE.
  character(len=*), parameter :: modn='(cheb_fftw3_1d_2a)'
contains
  subroutine mycft1d2a_open(nx,ny)
    implicit none
    integer, intent(in) :: nx,ny
    real(dp), dimension(nx+1,ny)  :: in,out
    
    integer  :: ta(2), inemb(2)
    
    n1 = nx;     d1 = nx+1;   
    n2 = ny;     d2 = ny

    ta(1) = d1; ta(2)=1
    inemb(1) = d1; inemb(2)=1

    call dfftw_plan_many_r2r(cplan%fwd, 2, ta, n2, in,inemb, 1, d1, out, inemb, 1, d1, FFTW_REDFT00, wisdom)
    
    scale = one/real(n1,dp)
    
    write(fftout,*) " MYCFT1D2A initialiation"//modn
    write(fftout,*) "fwd plan =   ",cplan%fwd
    init = .TRUE.
    return
  end subroutine mycft1d2a_open
  
  subroutine mycft1d2a_fwd(in, out)
    implicit none
    real(dp), dimension(d1,d2), intent(in) :: in
    real(dp), dimension(d1,d2), intent(out) :: out
    
    if (.not.init) stop modn
    call dfftw_execute_r2r(cplan%fwd, in, out)
    out(1,:) = half*scale*out(1,:)
    out(2:d1-1,:) = scale*out(2:d1-1,:)
    out(d1,:) = half*scale*out(d1,:)
    return
  end subroutine mycft1d2a_fwd
  
  subroutine mycft1d2a_bwd(in, out)
    implicit none
    real(dp), dimension(d1,d2), intent(in) :: in
    real(dp), dimension(d1,d2) :: tmp
    real(dp), dimension(d1,d2), intent(out) :: out

    if (.not.init) stop modn    
    tmp(1,:) = in(1,:)
    tmp(2:d1-1,:) = half*in(2:d1-1,:)
    tmp(d1,:) = in(d1,:)

    call dfftw_execute_r2r(cplan%fwd, tmp, out)
    return
  end subroutine mycft1d2a_bwd
  
  subroutine mycft1d2a_close
    if (init) then
       call dfftw_destroy_plan(cplan%fwd)
       call dfftw_destroy_plan(cplan%bwd)
       init = .false.
    end if
    return
  end subroutine mycft1d2a_close
  
end module cheb_fftw3_1d_2a

! cosft on (x, y) -> (y,n) 
! x  : 1st to 2nd dimension

module cheb_fftw3_1d_2a_1t2
  use four_fftw3_def
  implicit none
  
  private 
  public :: mycft1d2a_open, mycft1d2a_fwd, mycft1d2a_bwd, mycft1d2a_close
  
  type(pfftw3) :: cplan
  real(dp) :: scale
  
  integer :: n1, d1, n2, d2

  logical :: init=.FALSE.
  character(len=*), parameter :: modn='(cheb_fftw3_1d_2a_1t2)'
contains
  subroutine mycft1d2a_open(nx,ny)
    implicit none
    integer, intent(in) :: nx,ny
    real(dp), dimension(0:nx,ny)  :: in
    real(dp), dimension(ny, 0:nx)  :: out
    
    integer  :: ta(2), iinemb(2), oinemb(2)

    n1 = nx;     d1 = nx+1;   
    n2 = ny;     d2 = ny

    ta(1) = d1; ta(2)=1
    iinemb(1) = d1; iinemb(2)=1
    oinemb(1) = d1; oinemb(2)=1

    call dfftw_plan_many_r2r(cplan%fwd, 2, ta, n2, in,iinemb, 1, d1, out, oinemb, d2,1, FFTW_REDFT00, wisdom)

    call dfftw_plan_many_r2r(cplan%bwd, 2, ta, n2, out,oinemb, d2, 1, in, iinemb, 1,d1, FFTW_REDFT00, wisdom)
    
    scale = one/real(n1,dp)
    
    write(fftout,*) " MYCFT1D2A initialiation"//modn
    write(fftout,*) "fwd plan =   ",cplan%fwd
    write(fftout,*) "bckwd plan = ",cplan%bwd
    
    init = .TRUE.
    return
  end subroutine mycft1d2a_open
  
  subroutine mycft1d2a_fwd(in, out)
    implicit none
    real(dp), dimension(d1,d2), intent(in) :: in
    real(dp), dimension(d2,d1), intent(out) :: out

    if (.not.init) stop modn

    call dfftw_execute_r2r(cplan%fwd, in, out)

    out(:,1) = half*scale*out(:,1)
    out(:,2:d1-1) = scale*out(:,2:d1-1)
    out(:,d1) = half*scale*out(:,d1)
    
    return
  end subroutine mycft1d2a_fwd
  
  subroutine mycft1d2a_bwd(in, out)
     implicit none
     real(dp), dimension(d2,d1), intent(in) :: in
     real(dp), dimension(d2,d1) :: tmp
     real(dp), dimension(d1,d2), intent(out) :: out
     
     if (.not.init) stop modn
     tmp(:,1) = in(:,1)
     tmp(:,2:d1-1) = half*in(:,2:d1-1)
     tmp(:,d1) = in(:,d1)
     
     call dfftw_execute_r2r(cplan%bwd, tmp, out)
     return
   end subroutine mycft1d2a_bwd
   
   subroutine mycft1d2a_close
     if (init) then
        call dfftw_destroy_plan(cplan%fwd)
        call dfftw_destroy_plan(cplan%bwd)
        init=.false.
     end if
     return
   end subroutine mycft1d2a_close
end module cheb_fftw3_1d_2a_1t2

module fcheb_fftw3_2d_1t2
  use four_fftw3_def
  implicit none
  
  private 
  public :: mycft2d_open, mycft2d_fwd, mycft2d_bwd, mycft2d_close
  
  type(pfftw3) :: cplan, fplan
  real(dp) :: cscale, fscale
  
  integer :: n1, d1, n2, d2

  logical :: init=.FALSE.
  character(len=*), parameter :: modn='(fcheb_fftw3_2d_1t2)'

contains
  subroutine mycft2d_open(nx,ny)
    implicit none
    integer, intent(in) :: nx,ny
    real(dp), dimension(0:nx,ny+2)  :: in
    real(dp), dimension(ny+2, 0:nx)  :: med
    complex(dpc), dimension(ny/2+1,0:nx) :: out
    
    integer  :: cta(2), cineb(2), coneb(2)
    integer  :: fta(2), fineb(2), foneb(2)

    n1 = nx;     d1 = nx+1;   
    n2 = ny;     d2 = ny+2;

    cta(1) = d1; cta(2)=1
    cineb(1) = d1; cineb(2)=1
    coneb(1) = d1; coneb(2)=1

    call dfftw_plan_many_r2r(cplan%fwd, 2, cta, n2, in,cineb, 1, d1, med, coneb, d2,1, FFTW_REDFT00, wisdom)

    call dfftw_plan_many_r2r(cplan%bwd, 2, cta, n2, med,coneb, d2, 1, in, cineb, 1,d1, FFTW_REDFT00, wisdom)

    cscale = one/real(n1,dp)

    fta(1) = n2; fta(2)=1
    fineb(1) = d2; fineb(2)=1
    foneb(1) = d2/2; foneb(2)=1
    
    call dfftw_plan_many_dft_r2c(fplan%fwd,2, fta, d1, med, fineb, 1, d2, out, foneb, 1, d2/2, wisdom)
    call dfftw_plan_many_dft_c2r(fplan%bwd,2, fta, d1, out, foneb, 1, d2/2, med, fineb, 1, d2,wisdom)
    
    fscale = one/real(n2,dp)
    
    write(fftout,*) " MYCFT2D initialiation"//modn
    write(fftout,*) "Cosft fwd plan =   ",cplan%fwd
    write(fftout,*) "Cosft bckwd plan = ",cplan%bwd
    write(fftout,*) "Four. fwd plan =   ",fplan%fwd
    write(fftout,*) "Four. bckwd plan = ",fplan%bwd
    
    init = .TRUE.
    return
  end subroutine mycft2d_open
  
  subroutine mycft2d_fwd(in, out)
    implicit none
    real(dp), dimension(d1,d2), intent(in) :: in
    real(dp), dimension(d2,d1)  :: med
    complex(dpc), dimension(d2/2,d1), intent(out) :: out

    if (.not.init) stop modn
    call dfftw_execute_r2r(cplan%fwd, in, med)

    med(:,1) = half*cscale*med(:,1)
    med(:,2:d1-1) = cscale*med(:,2:d1-1)
    med(:,d1) = half*cscale*med(:,d1)

    call dfftw_execute_dft_r2c(fplan%fwd, med, out)
    out = fscale*out

    return
  end subroutine mycft2d_fwd
  
  subroutine mycft2d_bwd(in, out)
    implicit none
     complex(dpc), dimension(d2/2,d1), intent(in) :: in
     complex(dpc), dimension(d2/2,d1) :: tmp
     real(dp), dimension(d2,d1) :: med
     real(dp), dimension(d1,d2), intent(out) :: out

    if (.not.init) stop modn
     tmp = in
     call dfftw_execute_dft_c2r(fplan%bwd, tmp, med)
     med(:,2:d1-1) = half*med(:,2:d1-1)
     call dfftw_execute_r2r(cplan%bwd, med, out)

     return
   end subroutine mycft2d_bwd
   
   subroutine mycft2d_close

    if (.not.init) stop modn
     call dfftw_destroy_plan(cplan%fwd)
     call dfftw_destroy_plan(cplan%bwd)
     
     call dfftw_destroy_plan(fplan%fwd)
     call dfftw_destroy_plan(fplan%bwd)
     return
  end subroutine mycft2d_close
end module fcheb_fftw3_2d_1t2

module fcheb_fftw3_2d_1t2_ver_1d
  use four_fftw3_def
  implicit none
  
  private 
  public :: mycft2d_open, mycft2d_fwd, mycft2d_bwd, mycft2d_close
  
  type(pfftw3) :: cplan, fplan
  real(dp) :: cscale, fscale
  
  integer :: n1, d1, n2, d2

  logical :: init=.FALSE.
  character(len=*), parameter :: modn='(fcheb_fftw3_2d_1t2_ver_1d)'

contains
  subroutine mycft2d_open(nx,ny)
    implicit none
    integer, intent(in) :: nx,ny
    real(dp), dimension(0:nx)  :: in1, out1
    real(dp), dimension(ny+2)  :: in2
    complex(dpc), dimension(ny/2+1) :: out2
    
    integer  :: ta, inemb, outemb

    n1 = nx;     d1 = nx+1;   
    ta = d1
    inemb = d1
    call dfftw_plan_many_r2r(cplan%fwd, 1, ta, 1, in1,inemb, 1, 1, out1, inemb, 1, 1, FFTW_REDFT00, wisdom)
    cscale = one/real(n1,dp)

    write(fftout,*) " MYCFT2D initialiation"//modn
    write(fftout,*) "Cosft fwd plan =   ",cplan%fwd

    n2 = ny;     d2 = ny+2;
    ta = n2
    inemb = d2;     outemb = d2/2
    
    call dfftw_plan_many_dft_r2c(fplan%fwd,1, ta, 1, in2, inemb, 1, 1, out2, outemb, 1, 1, wisdom)
    call dfftw_plan_many_dft_c2r(fplan%bwd,1, ta, 1, out2, outemb, 1, 1, in2, inemb, 1, 1,wisdom)

    fscale = one/real(n2,dp)    
    
    write(fftout,*) "Four. fwd plan =   ",fplan%fwd
    write(fftout,*) "Four. bckwd plan = ",fplan%bwd
    
    init = .TRUE.
    return
  end subroutine mycft2d_open
  
  subroutine mycft2d_fwd(in, out)
    implicit none
    real(dp), dimension(d1,d2), intent(in) :: in
    real(dp), dimension(d2,d1)  :: med
    complex(dpc), dimension(d2/2,d1), intent(out) :: out
    integer :: i1

    if (.not.init) stop modn
    do i1=1,d2
       call dfftw_execute_r2r(cplan%fwd, in(:,i1), med(i1,:))
    enddo

    med(:,1) = half*med(:,1)
!    med(:,2:d1-1) = med(:,2:d1-1)
    med(:,d1) = half*med(:,d1)
    
    do i1=1,d1
       call dfftw_execute_dft_r2c(fplan%fwd, med(:,i1), out(:,i1))
    enddo
    out = (cscale*fscale)*out

    return
  end subroutine mycft2d_fwd
  
  subroutine mycft2d_bwd(in, out)
    implicit none
    complex(dpc), dimension(d2/2,d1), intent(in) :: in
    complex(dpc), dimension(d2/2,d1) :: tmp
    real(dp), dimension(d2,d1) :: med
    real(dp), dimension(d1,d2), intent(out) :: out
    integer :: i1

    if (.not.init) stop modn
     tmp = in
     do i1=1, d1
        call dfftw_execute_dft_c2r(fplan%bwd, tmp(:,i1), med(:,i1))
     enddo
     med(:,2:d1-1) = half*med(:,2:d1-1)
     do i1=1, d2
        call dfftw_execute_r2r(cplan%fwd, med(i1,:), out(:,i1))
     enddo

     return
   end subroutine mycft2d_bwd
   
   subroutine mycft2d_close

    if (.not.init) stop modn
     call dfftw_destroy_plan(cplan%fwd)
     call dfftw_destroy_plan(cplan%bwd)
     
     call dfftw_destroy_plan(fplan%fwd)
     call dfftw_destroy_plan(fplan%bwd)
     return
  end subroutine mycft2d_close
end module fcheb_fftw3_2d_1t2_ver_1d

module fcheb_fftw3_jacb_low
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
    
    write(fftout,*) " MYCFT2D initialiation"//modn
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
    
    write(fftout,*) " MYCFT2D initialiation"//modn
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

