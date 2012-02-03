module mymath
  use mydef
  implicit none

  interface Gaussian
     module procedure Gaussian_2D_r,Gaussian_2D_c, gaussian_2d_c_phase, gaussian_1d_r
  end interface

  interface vec_prod
     module procedure vec_prod_complex, vec_prod_real
  end interface
    
  interface vec_mod
     module procedure vec_mod_r,vec_mod_c, vec_mod_c_2d
  end interface

  interface swap_xy
     module procedure swap_xy_complex_0d, swap_xy_complex_1d
  end interface

  interface inner_prod
     module procedure &
          inner_prod_real_1d, &
          inner_prod_complex_1d
  end interface
contains
  function inner_prod_real_1d(vr1, met, vr2) result(prod)
    real(dp), dimension(:), intent(in) :: vr1, vr2, met
    real(dp) :: prod

    prod = dot_product(vr1, met*vr2)
    
    return
  end function inner_prod_real_1d
  
  function inner_prod_complex_1d(vr1, met, vr2) result(prod)
    complex(dpc), dimension(:), intent(in) :: vr1, vr2
    real(dp), dimension(:), intent(in) :: met
    complex(dpc) :: prod
    
    prod = dot_product(vr1, met*vr2)
    
    return
  end function inner_prod_complex_1d
  
  function frobenius_norm(mat) result(norm)
    implicit none
    complex(dpc), dimension(:,:) :: mat
    integer :: m, n,i1
    real(dp) :: norm, tmp

    m = size(mat,dim=1)
    n = size(mat,dim=2)
    
    tmp = zero
    if (m == n) then
       do i1=1, m
          tmp = tmp+sum(conjg(mat(:,i1))*mat(:,i1))
       end do
       norm = tmp
    end if
    return
  end function frobenius_norm
          
  function vec_mod_r(n,x) result(mod)
    integer, intent(in) :: n
    real(dp), intent(in), dimension(n) :: x
    real(dp) :: mod
    integer :: i1
    
    mod=zero
    do i1 = 1, n
       mod = mod+x(i1)*x(i1)
    end do
    mod = sqrt(mod)
    return
  end function vec_mod_r
  
  function vec_mod_c(n,x) result(mod)
    integer, intent(in), optional :: n
    complex(dpc), dimension(:), intent(in) :: x
    real(dp) :: mod
    integer :: i1

    if (present(n)) then
       mod = zero
       do i1=1,n
          mod = mod+conjg(x(i1))*x(i1)
       enddo
    else
       mod=sum(conjg(x)*x)
    end if
    mod = sqrt(mod)

    return
  end function vec_mod_c

  function vec_mod_c_2d(m,n,x) result(mod)
    integer, intent(in) :: n,m
    complex(dpc), intent(in), dimension(m,n) :: x
    real(dp) :: mod
    integer :: i1,i2
    
    mod=zero
    do i1 = 1, n
       do i2= 1, m
          mod = mod+x(i2,i1)*conjg(x(i2,i1))
       end do
    end do
    mod = sqrt(mod)
    return
  end function vec_mod_c_2d

  function vec_prod_complex(a,b) result(tot)
    implicit none
    complex(dpc), dimension(:) :: a, b
    complex(dpc) :: tot
    
    tot = sum(conjg(a)*b)
    return
  end function vec_prod_complex
  
  function vec_prod_real(a,b) result(tot)
    implicit none
    real(dp), dimension(:) :: a, b
    real(dp) :: tot
    
    tot = sum(a*b)
    return
  end function vec_prod_real  
     
  subroutine Gaussian_1d_r(n, x, amp, x0, dx0, data, np_op)
    implicit none
    integer, intent(in) :: n
    real(dp), intent(in) :: x(n), amp, x0, dx0
    real(dp), dimension(n), intent(out) :: data
    integer, optional :: np_op
    integer :: np,i1
    real(dp) :: tmp

    if (present(np_op)) then
       np = np_op
    else
       np = 2
    end if

    do i1 = 1, n
       tmp = (x(i1)-x0)/dx0
       data(i1) = amp*exp(-tmp**np)
    end do
  end subroutine Gaussian_1d_r

  !---- Gaussian
  subroutine Gaussian_2D_r(nx,ny,x,y,amp,x0,y0,dx0,dy0,data)
    implicit none
    integer, intent(in) :: nx,ny
    real(dp), intent(in) :: x(nx), y(ny), x0,y0,dx0, dy0, amp
    real(dp), dimension(nx,ny) :: data
    
    !local
    integer :: i1,i2
    real(dp) :: sigx, sigy
    
    sigx = dx0*dx0
    sigy = dy0*dy0
    forall(i1=1:nx,i2=1:ny) 
       data(i1,i2) = exp(-(x(i1)-x0)**2/sigx-(y(i2)-y0)**2/sigy)*amp
    end forall
    return
  end subroutine Gaussian_2D_r

  subroutine Gaussian_2D_c(nx,ny,x,y,amp,x0,y0,dx0,dy0,data)
    implicit none
    integer, intent(in) :: nx,ny
    real(dp), intent(in) :: x(nx), y(ny), x0,y0,dx0, dy0, amp
    complex(dpc), dimension(nx,ny) :: data
    
    !local
    integer :: i1,i2
    real(dp) :: sigx, sigy

    sigx = dx0*dx0
    sigy = dy0*dy0
    forall(i1=1:nx,i2=1:ny) 
       data(i1,i2) = exp(-(x(i1)-x0)**2/sigx-(y(i2)-y0)**2/sigy)*amp
    end forall
    return
  end subroutine Gaussian_2D_c

  subroutine Gaussian_2D_c_phase(nx,ny,x,y,amp,x0,y0,dx0,dy0,phase,data)
    implicit none
    integer, intent(in) :: nx,ny
    real(dp), intent(in) :: x(nx), y(ny), x0,y0,dx0, dy0, amp
    real(dp), dimension(nx,ny), intent(in) :: phase
    complex(dpc), dimension(nx,ny) :: data
    
    !local
    integer :: i1,i2
    real(dp) :: sigx, sigy
    
    sigx = dx0*dx0
    sigy = dy0*dy0
    forall(i1=1:nx,i2=1:ny) 
       data(i1,i2) = exp(-(x(i1)-x0)**2/sigx-(y(i2)-y0)**2/sigy+ic*phase(i1,i2))*amp
    end forall
    return
  end subroutine Gaussian_2D_c_phase

  subroutine swap_xy_complex_0d(x,y)
    complex(dpc) :: x, y
    complex(dpc) :: tmp

    tmp = x
    x = y
    y = tmp
    return
  end subroutine swap_xy_complex_0d

  subroutine swap_xy_complex_1d(x,y)
    complex(dpc), dimension(:) :: x, y
    complex(dpc), dimension(size(x)):: tmp

    tmp = x
    x = y
    y = tmp
    return
  end subroutine swap_xy_complex_1d


end module mymath
