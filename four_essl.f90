module fcheb_essl_2d_1t2
  use mydef
  implicit none
  
  private
  real(dp) :: cscale, fscale
  
  ! used for (x,y) -> (y,n)
  integer :: ncaux1, ncaux2
  real(dp), dimension(:), allocatable :: caux2, caux1bwd, caux1fwd
  
  ! used for (y,n) -> (ky, n)
  integer :: nfaux1, nfaux2
  real(dp), dimension(:), allocatable :: faux2, faux1bwd, faux1fwd

  integer  :: n1, n2, d1, d2, d2c, n1t
  public :: mycft2d_open, mycft2d_fwd, mycft2d_bwd, mycft2d_close
  

contains
  subroutine mycft2d_open(nx, ny)
    implicit none
    integer :: nx, ny
    real(dp), dimension(nx+1,ny) :: f
    real(dp), dimension(ny+2, nx+1) :: g
    complex(dpc), dimension(ny/2+1, nx+1) :: h

    n1 = nx; d1 = nx+1; n1t = 2*nx
    n2 = ny; d2 = ny+2; d2c = ny/2+1
    
    cscale = two/nx;  fscale = one/ny    
    if (nx .gt. 16384) then
       ncaux1 = 20000+0.6*n1t
       ncaux2 = 20000+0.64*n1t
    else
       ncaux1 = 35000
       ncaux2 = 20000
    end if
    ncaux2 = ncaux2+(nx/2+257)*int(min(128.,  real(ny)))
    
    allocate(caux2(ncaux2))
    allocate(caux1fwd(ncaux1))
    allocate(caux1bwd(ncaux1))

    call dcosf(1, f, 1, d1, g, d2, 1, n1t, n2, cscale, caux1fwd, ncaux1, caux2,ncaux2)
    call dcosf(1, g, d2,1, f, 1,d1, n1t, n2, one, caux1bwd, ncaux1, caux2,ncaux2)
    
    if (ny .gt. 4096) then
       nfaux1 = 20000+1.64*ny
       nfaux2 = 20000+1.14*ny
    else
       nfaux1 = 22000
       nfaux2 = 20000
    end if
    
    allocate(faux2(nfaux2))
    allocate(faux1fwd(nfaux1))
    allocate(faux1bwd(nfaux1))
    
    call drcft(1, g, d2, h, d2c, n2,d1, 1, fscale, faux1fwd, nfaux1, faux2,nfaux2)
    call dcrft(1, h, d2c, g, d2, n2,d1, -1, one, faux1bwd, nfaux1, faux2,nfaux2)

    return
  end subroutine mycft2d_open
  
  subroutine mycft2d_fwd(f,h)
    real(dp), dimension(d1,d2), intent(in) :: f
    real(dp), dimension(d2,d1) :: g
    complex(dpc), dimension(d2c,d1) :: h
    
    call dcosf(0, f, 1, d1, g, d2,1,n1t, n2, cscale, caux1fwd, ncaux1, caux2,ncaux2)
    g(:,1) = half*g(:,1); g(:,d1) = half*g(:,d1)
    call drcft(0, g, d2, h, d2c, n2,d1, 1, fscale, faux1fwd, nfaux1, faux2,nfaux2)
    
    return
  end subroutine mycft2d_fwd
  
  subroutine mycft2d_bwd(h,g)
    complex(dpc), dimension(d2c,d1) :: h
    real(dp), dimension(d2, d1)  :: f
    real(dp), dimension(d1,d2), intent(out) :: g
    integer :: i1

    call dcrft(0, h, d2c, f, d2, n2,d1, -1, one, faux1bwd, nfaux1, faux2,nfaux2) 
    f(:,1) = two*f(:,1); f(:,d1) = two*f(:,d1)
    call dcosf(0, f, d2,1, g, 1,d1, n1t, n2, one, caux1bwd, ncaux1, caux2,ncaux2)

    return
  end subroutine mycft2d_bwd
  
  subroutine mycft2d_close
    deallocate(caux2,caux1fwd, caux1bwd)
    deallocate(faux2,faux1fwd, faux1bwd)
    return
  end subroutine mycft2d_close
  
end module fcheb_essl_2d_1t2

module fcheb_essl_jacb_32
  use mod_const
  implicit none
  
  private
  
  integer, parameter :: onx = nx*3/2
  integer, parameter :: ony = ny*3/2, hony = ony/2, onym1=ony-1,onyp1=ony+1
  
  integer, parameter :: od1 = onx+1, od2 = ony+2, on1 = onx+1, on2=ony
  real(dp), parameter :: cscale0 = two/dble(onx), fscale0 = one/real(ony,dp) 
 
  integer :: ncaux1_0, ncaux2_0
  real(dp), dimension(:), allocatable :: caux1fwd0, caux1bwd0
  real(dp), dimension(:,:), allocatable :: caux223

  integer :: nfaux123, nfaux223
  real(dp), dimension(:), allocatable :: faux1fwd0, faux1bwd0
  real(dp), dimension(:,:), allocatable :: faux223
  
  complex(dpc), dimension(0:hny) :: iky
  real(dp), dimension(0:nx) :: cb_v
  real(dp) :: lx

  public :: myjacb_open, myjacb, myjacb_close
contains
  !---------------------------------------------
  function fdy(x) result(y)
    complex(dpc), dimension(0:hny,0:nx) :: x,y
    integer :: i1

    do i1=0,nx
       y(:,i1) = iky*x(:,i1)
    end do

    return
  end function fdy
!---------------------------------------------

  function cdx(x) result(y)
    complex(dpc), dimension(0:hny,0:nx), intent(in) :: x
    complex(dpc), dimension(0:hny,0:nx) :: y
    integer :: i1,i2,i3
    real(dp) :: tmp(0:nx)
    
    tmp = two/Lx/cb_v
    do i1=0,hny
       do i2=0, nx-1
          y(i1,i2) = zero
          do i3=i2+1,nx,2
             y(i1,i2) = y(i1,i2)+tmp(i2)*i3*x(i1,i3)
          end do
       end do
    end do
    y(:, nx) = zero
    return
  end function cdx
!----------------------------------------------

  subroutine myjacb_open(iky0, lx0)
    real(dp), dimension(0:onx,0:onym1) :: f
    real(dp), dimension(0:onyp1, 0:onx) :: g
    complex(dpc), dimension(0:hony, 0:onx) :: h
    integer :: i1
    complex(dpc), dimension(0:hny) :: iky0
    real(dp) :: lx0

    ! initialize for fdy, cdx
    iky = iky0
    cb_v(0) = two; cb_v(1:nxm1) = one; cb_v(nx) = two
    lx = lx0

    if (onx .gt. 16384) then
       ncaux1_0 = 30000+0.3*onx*2
       ncaux2_0 = 20000+0.32*onx*2
    else
       ncaux1_0 = 50000
       ncaux2_0 = 25000
    end if
    ncaux2_0 = ncaux2_0+(onx/2+257)*int(min(128.,  real(ony)))
    
    allocate(caux223(ncaux2_0,4))
    allocate(caux1fwd0(ncaux1_0))
    allocate(caux1bwd0(ncaux1_0))
    
    call dcosf(1, f, 1, od1, g,od2, 1, onx*2, on2, cscale0, caux1fwd0, ncaux1_0, caux223,ncaux2_0)
    call dcosf(1, g, od2,1, f, 1,od1, onx*2, on2, one, caux1bwd0, ncaux1_0, caux223,ncaux2_0)
    
    if (ny .gt. 4096) then
       nfaux123 = 30000+1.64*ony
       nfaux223 = 20000+1.14*ony
    else
       nfaux123 = 32000
       nfaux223 = 20000
    end if
    
    allocate(faux223(nfaux223,4))
    allocate(faux1fwd0(nfaux123))
    allocate(faux1bwd0(nfaux123))
    
    call drcft(1, g, od2, h, od2/2, on2, on1, 1, fscale0, faux1fwd0, nfaux123, faux223,nfaux223)
    call dcrft(1, h, od2/2, g, od2, on2, on1, -1, one, faux1bwd0, nfaux123, faux223,nfaux223)

    return
  end subroutine myjacb_open
  
  subroutine mycft2d32_fwd(f,h)
    real(dp), dimension(0:onx,0:onyp1), intent(in) :: f
    real(dp), dimension(0:onyp1,0:onx) :: g
    complex(dpc), dimension(0:hony,0:onx) :: h

    call dcosf(0, f, 1, od1, g, od2,1,onx*2, on2, cscale0, caux1fwd0, ncaux1_0, caux223,ncaux2_0)
    g(:,0) = half*g(:,0); g(:,onx) = half*g(:,onx)
    call drcft(0, g, od2, h, od2/2, on2,on1, 1, fscale0, faux1fwd0, nfaux123, faux223,nfaux223)
    
    return
  end subroutine mycft2d32_fwd
  
  subroutine mycft2d32_bwd(h,g,n)
    complex(dpc), dimension(0:hony,0:onx) :: h
    real(dp), dimension(0:onyp1, 0:onx)  :: f
    real(dp), dimension(0:hnx,0:nyp1), intent(out) :: g
    integer :: i1,n,n0=1
    optional :: n

    if (present(n))  n0 = n
    call dcrft(0, h, od2/2, f, od2, on2,on1, -1, one, faux1bwd0, nfaux123, faux223(:,n0),nfaux223) 
    f(:,0) = two*f(:,0); f(:,onx) = two*f(:,onx)
    call dcosf(0, f, od2,1, g, 1,od1, onx*2, on2, one, caux1bwd0, ncaux1_0, caux223(:,n0),ncaux2_0)
    
    return
  end subroutine mycft2d32_bwd
  
  subroutine myjacb_close
    deallocate(caux223,caux1fwd0, caux1bwd0)
    deallocate(faux223,faux1fwd0, faux1bwd0)
    return
  end subroutine myjacb_close

  function myjacb(f,g) result(h)
    complex(dpc), dimension(0:hny, 0:nx), intent(in) :: f,g
    complex(dpc), dimension(0:hony, 0:onx) :: fx, fy, gx, gy
    real(dp), dimension(0:onx, 0:onyp1) :: rfx, rfy, rgx, rgy
    complex(dpc), dimension(0:hny, 0:nx)  :: h
    equivalence (fx, rfx), (fy, rfy), (gx, rgx), (gy, rgy)

    !$omp parallel sections shared(f,g)
    !$omp section
    fx(0:hny,0:nx) = cdx(f)
    fx(0:hny,nx+1:onx) = zero; fx(hny+1:hony,:) = zero
    call mycft2d32_bwd(fx, rfx,1)
     
    
    !$omp section
    fy(0:hny,0:nx) = fdy(f)
    fy(0:hny,nx+1:onx) = zero; fy(hny+1:hony,:) = zero
    call mycft2d32_bwd(fy, rfy,2)
    
    !$omp section
    gx(0:hny,0:nx) = cdx(g)
    gx(0:hny,nx+1:onx) = zero; gx(hny+1:hony,:) = zero
    call mycft2d32_bwd(gx, rgx,3)
    
    !$omp section
    gy(0:hny,0:nx) = fdy(g)
    gy(0:hny,nx+1:onx) = zero; gy(hny+1:hony,:) = zero
    call mycft2d32_bwd(gy, rgy,4)

    !$omp end parallel sections

    rfx = rfx*rgy-rfy*rgx; call mycft2d32_fwd(rfx, fx)
    h(0:hny-1,0:nx) = fx(0:hny-1, 0:nx);     h(hny,0:nx) = fx(hny, 0:nx)
    return
  end function myjacb
end module fcheb_essl_jacb_32


module fcheb_essl_jacb_23
  use mod_const
  use mod_fcheb, only : fdy, cdx
  use fcheb_essl_2d_1t2
  
  implicit none
  private 
  integer, parameter :: nxa=2*nx/3, nya = 2*hny/3

  public :: myjacb
contains
  
  function myjacb(f,g) result(h)
    complex(dpc), dimension(0:hny, 0:nx), intent(in) :: f,g
    complex(dpc), dimension(0:hny, 0:nx) :: fx, fy, gx, gy
    real(dp), dimension(0:nx, 0:nyp1) :: rfx, rfy, rgx, rgy
    complex(dpc), dimension(0:hny, 0:nx) :: h, tmp1, tmp
    equivalence (fx, rfx), (fy, rfy), (gx, rgx), (gy, rgy)
    
    tmp = f; tmp(:,nxa:nx) = zero; tmp(nya:hny,:)=zero
    tmp1 = g; tmp1(:,nxa:nx) = zero; tmp1(nya:hny,:)=zero
    
    !$omp sections

    !$omp section
    fx = cdx(tmp); call mycft2d_bwd(fx, rfx)

    !$omp section
    fy = fdy(tmp); call mycft2d_bwd(fy, rfy)

    !$omp section
    gx = cdx(tmp1); call mycft2d_bwd(gx, rgx)

    !$omp section
    gy = fdy(tmp1); call mycft2d_bwd(gy, rgy)

    !$omp end sections

    rfx = rfx*rgy-rfy*rgx; call mycft2d_fwd(rfx, h)
    return
  end function myjacb
end module fcheb_essl_jacb_23

module four_essl_1d_r2r
  use mydef
  implicit none
  private

  integer :: ncaux1, ncaux2
  real(dp), dimension(:), allocatable :: caux2, caux1bwd, caux1fwd
  
  integer :: n1, d1
  real(dp) :: scale

  public :: mycft1d_open, mycft1d_fwd, mycft1d_bwd, mycft1d_close
contains

  subroutine mycft1d_open(an1)
    integer, intent(in) :: an1
    real(dp), dimension(0:an1) :: f,g
    
    n1 = an1; d1=an1+1
 
    if (n1 .gt. 16384) then
       ncaux1 = 30000+0.3*2*n1
       ncaux2 = 20000+0.32*2*n1
    else
       ncaux1 = 50000
       ncaux2 = 25000
    end if
    
    allocate(caux2(ncaux2))
    allocate(caux1fwd(ncaux1))
    allocate(caux1bwd(ncaux1))

    scale = two/dble(n1)

    call dcosf(1, f, 1, 1, g, 1, 1, 2*n1, 1, scale, caux1fwd, ncaux1, caux2,ncaux2)
    call dcosf(1, g, 1, 1, f, 1, 1, 2*n1, 1, one, caux1bwd, ncaux1, caux2,ncaux2)
    
    return
  end subroutine mycft1d_open
  
  subroutine mycft1d_fwd(f,g)
    real(dp), dimension(d1), intent(in) :: f
    real(dp), dimension(d1), intent(out) :: g
    
    call dcosf(0, f, 1, 1, g, 1,1,2*n1, 1, scale, caux1fwd, ncaux1, caux2,ncaux2)
    g(1) = half*g(1); g(d1) = half*g(d1)
    
    return
  end subroutine mycft1d_fwd
  
  subroutine mycft1d_bwd(f,g)
    real(dp), dimension(d1), intent(in) :: f
    real(dp), dimension(d1), intent(out) :: g
    real(dp), dimension(d1) :: f0
    
    f0(1) = two*f(1); f0(2:d1-1) = f(2:d1-1); f0(d1) = two*f(d1)
    call dcosf(0, f0, 1,1, g, 1,1, 2*n1, 1, one, caux1bwd, ncaux1, caux2,ncaux2)
    
    return
  end subroutine mycft1d_bwd
  
  subroutine mycft1d_close
    deallocate(caux2,caux1fwd, caux1bwd)
    return
  end subroutine mycft1d_close
  
end module four_essl_1d_r2r

module four_essl_1d
  use mydef
  implicit none
  private
  
  integer :: nfaux1, nfaux2
  real(dp), dimension(:), allocatable :: faux2, faux1bwd, faux1fwd
  
  integer :: n1, d1, d1c
  real(dp) :: scale
  
  public :: myfft1d_open, myfft1d_fwd, myfft1d_bwd, myfft1d_close
contains

  subroutine myfft1d_open(an1)
    integer, intent(in) :: an1
    real(dp), dimension(an1+2) :: f
    real(dp), dimension(an1/2+1) :: g    

    n1 = an1; d1=an1+2; d1c = an1/2+1

    scale = one/real(n1,dp)
    
    if (n1 .gt. 4096) then
       nfaux1 = 20000+1.64*n1
       nfaux2 = 20000+1.14*n1
    else
       nfaux1 = 22000
       nfaux2 = 20000
    end if
    
    allocate(faux2(nfaux2))
    allocate(faux1fwd(nfaux1))
    allocate(faux1bwd(nfaux1))

    call drcft(1, f, d1, g, d1c, n1, 1, 1,scale, faux1fwd, nfaux1, faux2,nfaux2)
    call dcrft(1, g, d1c, f, d1, n1, 1, -1, one, faux1bwd, nfaux1, faux2,nfaux2)
    
    return
  end subroutine myfft1d_open
  
  subroutine myfft1d_fwd(f,g)
    real(dp), dimension(d1), intent(in) :: f
    complex(dpc), dimension(d1c), intent(out) :: g
    
    call drcft(0, f, d1, g, d1c, n1, 1, 1,scale, faux1fwd, nfaux1, faux2,nfaux2)

    return
  end subroutine myfft1d_fwd
  
  subroutine myfft1d_bwd(f,g)
    complex(dpc), dimension(d1c), intent(in) :: f
    real(dp), dimension(d1), intent(out) :: g

    call dcrft(0, f, d1c, g, d1, n1,1, -1, one, faux1bwd, nfaux1, faux2,nfaux2)    
    return
  end subroutine myfft1d_bwd
  
  subroutine myfft1d_close
    deallocate(faux2,faux1fwd, faux1bwd)
    return
  end subroutine myfft1d_close
  
end module four_essl_1d


module four_essl_1d_many
  use mydef
  implicit none
  private
  
  integer :: nfaux1, nfaux2
  real(dp), dimension(:), allocatable :: faux2, faux1bwd, faux1fwd
  
  integer :: n1, d1, d1c, n2, d2
  real(dp) :: scale
  
  public :: myfft1d_open, myfft1d_fwd, myfft1d_bwd, myfft1d_close
contains

  subroutine myfft1d_open(an1, an2)
    integer, intent(in) :: an1,an2
    real(dp), dimension(an1+2, an2) :: f
    real(dp), dimension(an1/2+1,an2) :: g    

    n1 = an1; d1=an1+2; d1c = an1/2+1
    n2 = an2; d2=an2
    scale = one/real(n1,dp)
    
    if (n1 .gt. 4096) then
       nfaux1 = 20000+1.64*n1
       nfaux2 = 20000+1.14*n1
    else
       nfaux1 = 22000
       nfaux2 = 20000
    end if
    
    allocate(faux2(nfaux2))
    allocate(faux1fwd(nfaux1))
    allocate(faux1bwd(nfaux1))

    call drcft(1, f, d1, g, d1c, n1, d2, 1,scale, faux1fwd, nfaux1, faux2,nfaux2)
    call dcrft(1, g, d1c, f, d1, n1, d2, -1, one, faux1bwd, nfaux1, faux2,nfaux2)
    
    return
  end subroutine myfft1d_open
  
  subroutine myfft1d_fwd(f,g)
    real(dp), dimension(d1,d2), intent(in) :: f
    complex(dpc), dimension(d1c,d2), intent(out) :: g
    
    call drcft(0, f, d1, g, d1c, n1, d2, 1,scale, faux1fwd, nfaux1, faux2,nfaux2)

    return
  end subroutine myfft1d_fwd
  
  subroutine myfft1d_bwd(f,g)
    complex(dpc), dimension(d1c,d2), intent(in) :: f
    real(dp), dimension(d1,d2), intent(out) :: g

    call dcrft(0, f, d1c, g, d1, n1,d2, -1, one, faux1bwd, nfaux1, faux2,nfaux2)    
    return
  end subroutine myfft1d_bwd
  
  subroutine myfft1d_close
    deallocate(faux2,faux1fwd, faux1bwd)
    return
  end subroutine myfft1d_close
  
end module four_essl_1d_many


