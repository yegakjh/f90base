!----------------------------------------------------------------
! Category : module
! Name : four_grid
! Purpose : Define grid size and coordinate information
!----------------------------------------------------------------
module mod_four2d
  use mod_const
  implicit none
  
  real(dp),dimension(-hnx:hnx-1) :: xx,kx

  real(dp),dimension(0:hnx) :: kxf,DerX2,kxp
  complex(dpc), dimension(0:hnx) :: DerX

  real(dp),dimension(-hny:hny-1) :: yy,ky

  real(dp),dimension(0:nym1) :: kyf,DerY2
  real(dp),dimension(0:hny) :: kyp

  complex(dpc),dimension(0:nym1) :: DerY

  real(dp) :: dx, dkx, Lx
  real(dp) :: dy, dky, Ly

  real(dp) :: efactor_k, sum_int_x,sum_int_k

  interface force_real
     module procedure force_real_2D, force_real_2D_multi
  end interface

  interface inverse
     module procedure inverse_2D
  end interface

  interface set_zero_23
     module procedure set_zero_23_xy_real, set_zero_23_xy_cmplx
  end interface
  
  interface statk2d
     module procedure statk2d_ampsq, statk2d_weight
  end interface
  
  interface autocorr
     module procedure autocorr_tot, autocorr_k
  end interface
  
  interface crosscorr_real
     module procedure crosscorr_tot_real, crosscorr_k_real
  end interface
  
  interface crosscorr_complex
     module procedure crosscorr_tot_real, crosscorr_k_real
  end interface

  interface crosscorr_ky
     module procedure crosscorr_ky_real
  end interface

contains

  subroutine four2d_open(delkx,delky)
    implicit none
    real(dp), intent(in) :: delkx, delky
    integer :: i

    dkx = delkx
    dky = delky

    ! set the box size
    Lx = twopi/dkx
    Ly = twopi/dky

    ! set x space grid
    dx = Lx/real(nx,dp)
    dy = Ly/real(ny,dp)
    
    xx  = (/(dx*real(i,dp), i = -hnx,hnx-1)/)
    yy  = (/(dy*real(i,dp), i = -hny,hny-1)/)

    ! set k space grid
    kx  = (/(dkx*real(i,dp), i = -hnx,hnx-1)/)
    ky  = (/(dky*real(i,dp), i = -hny,hny-1)/)
    

    ! set k space for fft
    kxf(0:hnx-1) = kx(0:hnx-1); kxf(hnx) = kx(-hnx)

    DerX = cmplx(zero, kxf); DerX(hnx) = zero
    DerX2 = -kxf*kxf

    kyf(0:hny-1) = ky(0:hny-1);     kyf(hny:nym1) = ky(-hny:-1)

    DerY = cmplx(zero, kyf);     DerY(hny) = zero
    DerY2 = -kyf*kyf

    sum_int_k = dkx*dky
    sum_int_x = dx*dy
    efactor_k = (Lx*Ly/twopi)**2/real(nx*ny,dp)*sum_int_k

    ! set kspace for spectrum analysis
    kxp = abs(kxf)
    kyp = abs(kyf(0:hny))
    
    return
  end subroutine four2d_open
  
  function Lapl2D(fin)
    implicit none
    complex(dpc), dimension(0:hnx,0:nym1) :: fin,lapl2d
    integer :: i
    
    do i=0,nym1
       lapl2d(:,i) = (DerX2+DerY2(i))*fin(:,i)
    end do
    return
  end function Lapl2D

  ! set_zero_23
  subroutine set_zero_23_xy_cmplx(data)
    complex(dpc),dimension(0:hnx,0:nym1), intent(inout) :: data

    data(hnx23+1:hnx,:) = zero
    data(0:hnx23,hny23+1:nym1-hny23) = zero

    return
  end subroutine set_zero_23_xy_cmplx

  subroutine set_zero_23_xy_real(data)
    real(dp),dimension(0:hnx,0:nym1), intent(inout) :: data

    data(hnx23+1:hnx,:) = zero
    data(0:hnx23,hny23+1:nym1-hny23) = zero

    return
  end subroutine set_zero_23_xy_real

  subroutine fft2d_xy_shuffle(in,out)
    implicit none
    real(dp), dimension(0:nxm1,0:nym1) :: in
    real(dp), dimension(-hnx:hnx-1,-hny:hny-1) :: out
    
    out(0:hnx-1,0:hny-1)=in(0:hnx-1,0:hny-1)
    out(0:hnx-1,hny-ny:-1)=in(0:hnx-1,hny:ny-1)
    out(hnx-nx:-1,0:hny-1)=in(hnx:nx-1,0:hny-1)
    out(hnx-nx:-1,hny-ny:-1)=in(hnx:nx-1,hny:ny-1)
    return
  end subroutine fft2d_xy_shuffle

  subroutine fft2d_kk_shuffle(in,out)
    implicit none
    complex(dp), dimension(0:hnx, 0:nym1) :: in
    complex(dp), dimension(-hnx:hnx-1,-hny:hny-1) :: out

    out(0:hnx-1, 0:hny-1) = in(0:hnx-1,0:hny-1)
    out(0:hnx-1,-hny:-1) =  in(0:hnx-1,hny:nym1)
    out(-hnx, 0:hny-1) = in(hnx,0:hny-1)
    out(-hnx, -hny:-1) = in(hnx,0:hny-1)
    out(-hnx+1:-1,-hny+1:hny-1) = conjg(out(hnx-1:1:-1,hny-1:-hny+1:-1))
    out(-hnx+1:-1,-hny) = conjg(out(hnx-1:1:-1,-hny))
    return
  end subroutine fft2d_kk_shuffle
  
  function fdy(f) result(g)
    implicit none
    complex(dpc), dimension(0:hnx,0:nym1) :: f,g
    integer :: i
    
    do i=0, nym1
       g(:,i) = DerY(i)*f(:,i)
    enddo
    return
  end function fdy
  
  function fdx(f) result(g)
    implicit none
    complex(dpc), dimension(0:hnx,0:nym1) :: f,g
    integer :: i
    
    do i=0, nym1
       g(:,i) = DerX*f(:,i)
    enddo
    return
  end function fdx
  
  function statx2d(xdata) result(stat)
    real(dp), dimension(0:nxp1,0:nym1) :: xdata
    real(dp), dimension(5) :: stat
    real(dp) :: avg, avg3rd, avg4th, normal
    
    normal = one/real(nx*ny,8)
    
    ! maximum
    stat(1) = maxval(abs(xdata(0:nxm1,0:nym1))) ! Abs Max Value
    ! rms value
    stat(2) = sqrt(sum(xdata(0:nxm1,0:nym1)**2)*normal) ! RMS 
    ! standard deviation
    avg = sum(xdata(0:nxm1,0:nym1))*normal ! Mean Value
    stat(3) = sqrt(sum((xdata(0:nxm1,0:nym1)-avg)**2)*normal) ! STD
    ! skewness
    avg3rd = sum(xdata(0:nxm1,0:nym1)**3)*normal ! 3rd order Mean Value
    stat(4) = sum((xdata(0:nxm1,0:nym1)-avg)**3)*normal/stat(3)**3 ! Skewness
    ! kurtosis
    avg4th = sum(xdata(0:nxm1,0:nym1)**4)*normal ! 4th order Mean Value
    stat(5) = sum((xdata(0:nxm1,0:nym1)-avg)**4)*normal/stat(3)**4-3.0_dp ! Kurtosis
    
    return
  end function statx2d
  
  function statk2d_ampsq(kdata) result(stat)
    ! weighting for kx, ky : Amplitude square.
    complex(dpc), dimension(0:hnx,0:nym1) :: kdata
    real(dp), dimension(6) :: stat
    real(dp) :: rmskx,rmsky,normal
    real(dp),dimension(0:hnx,0:nym1) :: sqkdata
    integer :: i, ind(1)
    
    normal = one/real((hnx+1)*ny,8)
    sqkdata = kdata*conjg(kdata)
    
    stat(1) = maxval(abs(kdata)) ! Abs Max Value
    stat(2) = sqrt(sum(sqkdata))*normal ! RMS 
    
    rmskx = zero
    rmsky = zero
    do i=0,nym1 
       rmskx = rmskx+sum((abs(kxf)*sqkdata(:,i)))
       rmsky = rmsky+sum(abs(kyf(i))*sqkdata(:,i))
    enddo
    
    stat(3) = sqrt(rmskx/stat(2))*normal ! k_x
    stat(4) = sqrt(rmsky/stat(2))*normal ! k_y
    
    ind = maxloc(abs(sqkdata(:,i)))-1
    stat(5) = abs(kxf(ind(1)))
    stat(6) = abs(kyf(ind(2)))
    
    return
  end function statk2d_ampsq
  
  function statk2d_weight(kdata,weight) result(stat)
    complex(dpc), dimension(0:hnx,0:nym1) :: kdata
    real(dp), dimension(6) :: stat
    real(dp) :: rmskx,rmsky,normal, tmp
    real(dp),dimension(0:hnx,0:nym1) :: weight, sqkdata
    integer :: i, ind(2)
    
    normal = one/real((hnx+1)*ny,8)
    sqkdata = kdata*conjg(kdata)
    
    stat(1) = maxval(abs(kdata)) ! Abs Max Value
    stat(2) = sqrt(sum(sqkdata)*normal) ! RMS 
    
    rmskx = zero
    rmsky = zero
    do i=0,nym1 
       rmskx = rmskx+sum(abs(kxf)*weight(:,i))
       rmsky = rmsky+sum(weight(:,i))*abs(kyf(i))
    enddo
    
    tmp =  sum(weight)
    stat(3) = sqrt(rmskx/tmp) ! rms k_x
    stat(4) = sqrt(rmsky/tmp) ! rms k_y
    
    ind = maxloc(abs(weight))-1
    stat(5) = abs(kxf(ind(1)))
    stat(6) = abs(kyf(ind(2)))
    return
  end function statk2d_weight
     
  subroutine get_energy_spectrum(x,y,z)
    implicit none
    real(dp), dimension(0:hnx,0:nym1) :: x
    real(dp), dimension(hnx) :: y
    integer :: i1,i2,nk
    real(dp) :: k,dk, z, multiplier
    
    y = 0;z=0
    dk = dkx
    do i2=0,ny-1
       do i1=0,hnx
          if (i1 .eq. 0 .or. i1.eq.hnx) then 
             multiplier = one
          else
             multiplier = two
             endif
             k = sqrt(-derY2(i2)-DerX2(i1))
             nk = ceiling(k/dk)
             if (nk .gt. 0 ) then 
                y(nk) = y(nk)+ x(i1,i2)*multiplier
             elseif (nk .gt. hnx) then
                z = z+x(i1,i2)*multiplier
             end if
          enddo
       enddo
     end subroutine get_energy_spectrum

  function autocorr_tot(f,korder,ek0) result(total)
    implicit none
    complex(dpc), dimension(0:hnx,0:nym1), intent(in) :: f
    integer, intent(in) :: korder
    real(dp), dimension(0:hnx,0:nym1) :: ek
    real(dp), optional :: ek0(2) ! (Zonal, Streamer)    
    real(dp) :: total
    
    integer :: i,j
    
    forall (i=0:hnx, j=0:nym1) 
       ek(i,j) = f(i,j)*conjg(f(i,j))    
    end forall

    select case (korder)
    case (0)
    case (2)
       forall (i=0:hnx, j=0:nym1) 
          ek(i,j) = ek(i,j)*(-DerX2(i)-DerY2(j))
       end forall
    case (4)
       forall (i=0:hnx, j=0:nym1) 
          ek(i,j) = ek(i,j)*(DerX2(i)+DerY2(j))*(DerX2(i)+DerY2(j))
       end forall
    case default
       write(*,*) ' korder should be 0,2, or 4. (autocorr, )'
    end select
    
    total = two*sum(ek(1:hnx-1,:))+sum(ek(0,:))+sum(ek(hnx,:))
    if (present(ek0)) then
       ek0(1) = two*sum(ek(1:hnx,0))
       ek0(2) = sum(ek(0,1:nym1))
    endif

    return
  end function autocorr_tot
  
  function autocorr_k(f,korder,ek,ek0) result(total)
    implicit none
    complex(dpc), dimension(0:hnx,0:nym1), intent(in) :: f
    integer, intent(in) :: korder
    real(dp), dimension(0:hnx,0:nym1) :: ek
    real(dp), optional :: ek0(2) ! (Zonal, Streamer)    
    real(dp) :: total
    
    integer :: i,j
    
    forall (i=0:hnx, j=0:nym1) 
       ek(i,j) = f(i,j)*conjg(f(i,j))    
    end forall

    select case (korder)
    case (0)
    case (2)
       forall (i=0:hnx, j=0:nym1) 
          ek(i,j) = ek(i,j)*(-DerX2(i)-DerY2(j))
       end forall
    case (4)
       forall (i=0:hnx, j=0:nym1) 
          ek(i,j) = ek(i,j)*(DerX2(i)+DerY2(j))*(DerX2(i)+DerY2(j))
       end forall
    case default
       write(*,*) ' korder should be 0,2, or 4. (autocorr_k, )'
    end select

    total = two*sum(ek(1:hnx-1,:))+sum(ek(0,:))+sum(ek(hnx,:))
    if (present(ek0)) then
       ek0(1) = two*sum(ek(1:hnx,0))
       ek0(2) = sum(ek(0,1:nym1))
    endif

    return
  end function autocorr_k

  function autocorr_ky(f,korder0) result(ek)
    implicit none
    complex(dpc), dimension(0:hnx,0:nym1), intent(in) :: f
    integer, intent(in), optional :: korder0
    real(dp), dimension(0:nym1)  :: ek
    
    real(dp) :: rtmpa(0:hnx)
    
    integer :: iy, korder
    integer, parameter :: hnxm1 = hnx-1
    
    korder = 0
    if (present(korder0)) korder = korder0

    select case (korder)
    case (0)
       do iy=0,nym1 
          ek(iy) = two*sum(f(1:hnxm1,iy)*conjg(f(1:hnxm1,iy))) &
               +f(0,iy)*conjg(f(0,iy))+f(hnx,iy)*conjg(f(hnx,iy))
       enddo
    case (2)
       do iy=0,nym1
          rtmpa = (-DerX2-DerY2(iy))
          ek(iy) = two*sum(rtmpa(1:hnxm1)*f(1:hnxm1,iy)*conjg(f(1:hnxm1,iy))) &
               +rtmpa(0)*f(0,iy)*conjg(f(0,iy)) &
               +rtmpa(hnx)*f(hnx,iy)*conjg(f(hnx,iy)) 
       end do
    case (4)
       do iy=0,nym1
          rtmpa = (DerX2+DerY2(iy))*(DerX2+DerY2(iy))
          ek(iy) = two*sum(rtmpa(1:hnxm1)*f(1:hnxm1,iy)*conjg(f(1:hnxm1,iy))) &
               +rtmpa(0)*f(0,iy)*conjg(f(0,iy)) &
               +rtmpa(hnx)*f(hnx,iy)*conjg(f(hnx,iy))
       end do
    case default
       write(*,*) ' korder should be 0 or 2. (autocorr_ky, )'
    end select
    
    return
  end function autocorr_ky
  
  function crosscorr_tot_real(f1,f2,korder,ek0) result(total)
    implicit none
    complex(dpc), dimension(0:hnx,0:nym1), intent(in) :: f1,f2
    integer, intent(in) :: korder
    real(dp), dimension(0:hnx,0:nym1) :: ek
    real(dp), optional :: ek0(2) ! (Zonal, Streamer)
    real(dp) :: total
    
    integer :: i,j

    ! take real (f1 f2^*)+(f1^* f2)
    forall (i=0:hnx, j=0:nym1) 
       ek(i,j) = two*f1(i,j)*conjg(f2(i,j))    
    end forall
    
    select case (korder)
    case (0)
    case (1) 
       forall (i=0:hnx, j=0:nym1) 
          ek(i,j) = ek(i,j)*sqrt((-DerX2(i)-DerY2(j)))
       end forall
    case (2)
       forall (i=0:hnx, j=0:nym1) 
          ek(i,j) = ek(i,j)*(-DerX2(i)-DerY2(j))
       end forall
    case (3)
       forall (i=0:hnx, j=0:nym1) 
          ek(i,j) = ek(i,j)*(-DerX2(i)-DerY2(j))*sqrt((-DerX2(i)-DerY2(j)))
       end forall
    case (4)
       forall (i=0:hnx, j=0:nym1) 
          ek(i,j) = ek(i,j)*(DerX2(i)+DerY2(j))*(DerX2(i)+DerY2(j))
       end forall
    case default
       write(*,*) ' korder should be between 0 and 4. (crosscorr, )'
    end select
    
    total = sum(ek)+sum(ek(1:hnx-1,:))
    
    total = two*sum(ek(1:hnx-1,:))+sum(ek(0,:))+sum(ek(hnx,:))
    if (present(ek0)) then
       ek0(1) = two*sum(ek(1:hnx,0))
       ek0(2) = sum(ek(0,1:nym1))
    endif

    return
  end function crosscorr_tot_real
  
  function crosscorr_ky_real(f1,f2, korder0) result(cor_ky)
    ! take real (f1 f2^*)+(f1^* f2)
    implicit none
    complex(dpc), dimension(0:hnx,0:nym1), intent(in) :: f1,f2
    integer, optional :: korder0
    integer :: korder
    real(dp), dimension(0:nym1) :: cor_ky

    real(dp), dimension(0:hnx) :: rtmpa
    integer :: iy
    integer, parameter :: hnxm1 = hnx-1

    !    real(dp), dimension(0:hnx,0:nym1) :: cor
    !    integer :: i,j
    
    !     forall (i=0:hnx,j=0:nym1) 
    !        cor(i,j) = two*f1(i,j)*conjg(f2(i,j))    
    !     end forall
    
    !     do i = 0, nym1
    !        cor_ky(i) = two*sum(cor(1:hnx-1,i))+cor(0,i)+cor(hnx,i)
    !     enddo

    korder = 0 !default
    if (present(korder0)) korder = korder0
    
    select case (korder)
    case (0)
       do iy=0,nym1 
          cor_ky(iy) = two*sum(f1(1:hnxm1,iy)*conjg(f2(1:hnxm1,iy))) &
               +f1(0,iy)*conjg(f2(0,iy))+f1(hnx,iy)*conjg(f2(hnx,iy))
       enddo
    case (2)
       do iy=0,nym1
          rtmpa = (-DerX2-DerY2(iy))
          cor_ky(iy) = two*sum(rtmpa(1:hnxm1)*f1(1:hnxm1,iy)*conjg(f2(1:hnxm1,iy))) &
               +rtmpa(0)*f1(0,iy)*conjg(f2(0,iy)) &
               +rtmpa(hnx)*f1(hnx,iy)*conjg(f2(hnx,iy)) 
       end do
    case (4)
       do iy=0,nym1
          rtmpa = (DerX2+DerY2(iy))*(DerX2+DerY2(iy))
          cor_ky(iy) = two*sum(rtmpa(1:hnxm1)*f1(1:hnxm1,iy)*conjg(f2(1:hnxm1,iy))) &
               +rtmpa(0)*f1(0,iy)*conjg(f2(0,iy)) &
               +rtmpa(hnx)*f1(hnx,iy)*conjg(f2(hnx,iy))
       end do
    case default
       write(*,*) ' korder should be 0 or 2. (crosscorr_ky, )'
    end select
    
    cor_ky = two*cor_ky ! since f1 f2^*+f1^* f2
    
    return
  end function crosscorr_ky_real

  function crosscorr_k_real(f1,f2,korder,ek,ek0) result(total)
    implicit none
    complex(dpc), dimension(0:hnx,0:nym1), intent(in) :: f1,f2
    integer, intent(in) :: korder
    real(dp), dimension(0:hnx,0:nym1) :: ek
    real(dp), optional :: ek0(2) ! (Zonal, Streamer)
    real(dp) :: total

    integer :: i,j

    ! take real (f1 f2^*)+(f1^* f2)
    forall (i=0:hnx, j=0:nym1) 
       ek(i,j) = two*f1(i,j)*conjg(f2(i,j))    
    end forall
    
    select case (korder)
    case (0)
    case (1) 
       forall (i=0:hnx, j=0:nym1) 
          ek(i,j) = ek(i,j)*sqrt((-DerX2(i)-DerY2(j)))
       end forall
    case (2)
       forall (i=0:hnx, j=0:nym1) 
          ek(i,j) = ek(i,j)*(-DerX2(i)-DerY2(j))
       end forall
    case (3)
       forall (i=0:hnx, j=0:nym1)
          ek(i,j) = ek(i,j)*(-DerX2(i)-DerY2(j))*sqrt((-DerX2(i)-DerY2(j)))
       end forall
    case (4)
       forall (i=0:hnx, j=0:nym1) 
          ek(i,j) = ek(i,j)*(DerX2(i)+DerY2(j))*(DerX2(i)+DerY2(j))
       end forall
    case default
       write(*,*) ' korder should be between 0 and 4. (crosscorr, )'
    end select
    
    total = two*sum(ek(1:hnx-1,:))+sum(ek(0,:))+sum(ek(hnx,:))
    if (present(ek0)) then
       ek0(1) = two*sum(ek(1:hnx,0))
       ek0(2) = sum(ek(0,1:nym1))
    endif
    return
  end function crosscorr_k_real
  
  function crosscorr_tot_complex(f1,f2,korder,ek0) result(total)
    implicit none
    complex(dpc), dimension(0:hnx,0:nym1), intent(in) :: f1,f2
    integer, intent(in) :: korder
    complex(dpc), dimension(0:hnx,0:nym1) :: ek
    complex(dpc), optional :: ek0(2) 
    ! (Zonal, Streamer)
    complex(dpc) :: total
    
    integer :: i,j

    ! take (f1 f2^*)
    forall (i=0:hnx, j=0:nym1)
       ek(i,j) = f1(i,j)*conjg(f2(i,j))    
    end forall
    select case (korder)
    case (0)
    case (1) 
       forall (i=0:hnx, j=0:nym1)
          ek(i,j) = ek(i,j)*sqrt((-DerX2(i)-DerY2(j)))
       end forall
    case (2)
       forall (i=0:hnx, j=0:nym1) 
          ek(i,j) = ek(i,j)*(-DerX2(i)-DerY2(j))
       end forall
    case (3)
       forall (i=0:hnx, j=0:nym1) 
          ek(i,j) = ek(i,j)*(-DerX2(i)-DerY2(j))*sqrt((-DerX2(i)-DerY2(j)))
       end forall
    case (4)
       forall (i=0:hnx, j=0:nym1) 
          ek(i,j) = ek(i,j)*(DerX2(i)+DerY2(j))*(DerX2(i)+DerY2(j))
       end forall
    case default
       write(*,*) ' korder should be between 0 and 4. (crosscorr, )'
    end select

    total = two*sum(ek(1:hnx-1,:))+sum(ek(0,:))+sum(ek(hnx,:))
    if (present(ek0)) then
       ek0(1) = two*sum(ek(1:hnx,0))
       ek0(2) = sum(ek(0,1:nym1))
    endif
    return
  end function crosscorr_tot_complex
  
  function crosscorr_k_complex(f1,f2,korder,ek,ek0) result(total)
    implicit none
    complex(dpc), dimension(0:hnx,0:nym1), intent(in) :: f1,f2
    integer, intent(in) :: korder
    complex(dpc), dimension(0:hnx,0:nym1) :: ek
    complex(dpc), optional :: ek0(2) ! (Zonal, Streamer)
    complex(dpc) :: total

    integer :: i,j

    ! take (f1 f2^*)
    forall (i=0:hnx, j=0:nym1) 
       ek(i,j) = f1(i,j)*conjg(f2(i,j))    
    end forall
    
    select case (korder)
    case (0)
    case (1) 
       forall (i=0:hnx, j=0:nym1) 
          ek(i,j) = ek(i,j)*sqrt((-DerX2(i)-DerY2(j)))
       end forall
    case (2)
       forall (i=0:hnx, j=0:nym1) 
          ek(i,j) = ek(i,j)*(-DerX2(i)-DerY2(j))
       end forall
    case (3)
       forall (i=0:hnx, j=0:nym1) 
          ek(i,j) = ek(i,j)*(-DerX2(i)-DerY2(j))*sqrt((-DerX2(i)-DerY2(j)))
       end forall
    case (4)
       forall (i=0:hnx, j=0:nym1) 
          ek(i,j) = ek(i,j)*(DerX2(i)+DerY2(j))*(DerX2(i)+DerY2(j))
       end forall
    case default
       write(*,*) ' korder should be between 0 and 4. (crosscorr, )'
    end select
    
    total = two*sum(ek(1:hnx-1,:))+sum(ek(0,:))+sum(ek(hnx,:))
    if (present(ek0)) then
       ek0(1) = two*sum(ek(1:hnx,0))
       ek0(2) = sum(ek(0,1:nym1))
    endif
    return
  end function crosscorr_k_complex
  
  subroutine inverse_2D(y0,y)
    complex(dpc), dimension(0:hnx, 0:nym1), intent(inout) :: y
    real(dp), intent(in) :: y0
    real(dp) :: kpsq
    integer :: iy, ix
    
    iyloop : do iy=0, nym1 
       if (iy > hny23 .and. iy <= nym1-hny23) cycle iyloop
       do ix=0, hnx
          kpsq = -(DerX2(ix)+DerY2(iy))
          y(ix,iy) = y(ix,iy)/(y0+kpsq)
       enddo
    enddo iyloop
    
    return
  end subroutine inverse_2D

  subroutine force_real_2D_multi(y)
    complex(dpc), dimension(0:hnx,0:nym1,ndim), intent(inout) :: y
    integer :: i1
    
    do i1=1, ndim
       y(0,nym1:hny+1:-1,i1)= conjg(y(0,1:hny-1,i1))    
    enddo
    
    return
  end subroutine force_real_2D_multi
  
  subroutine force_real_2D(y)
    complex(dpc), dimension(0:hnx,0:nym1), intent(inout) :: y
    integer :: i1
    
    y(0,nym1:hny+1:-1)= conjg(y(0,1:hny-1))    
    
    return
  end subroutine force_real_2D

  ! index change from regular-to-pack
  function i_r2p(nr,n1,np,n2) result(flag)
    ! translate index nr in 0:n1-1
    !    to index np with packed 0:2*n2
    integer, intent(out) :: np
    integer, intent(in) :: n1,n2,nr
    logical :: flag
    
    if (nr <= n2) then
       np = nr; flag = .true.
    else if (nr > n1-1-n2) then
       np = 2*n2-(n1-1)+nr ;  flag = .true.
    else 
       np = -1 ; flag = .false.
    end if
    return
  end function i_r2p

  ! index change from pack-to-regular
  function i_p2r(np,np0,nr,nr0) result(flag)
    ! translate index np with packed 0:2*np0
    ! index nr in 0:nr0-1
    !    to 
    integer, intent(out) :: nr
    integer, intent(in) :: np, np0, nr0
    logical :: flag
    
    if (np <= np0) then
       nr = np; flag = .true.
    else if (np > np0 ) then
       nr = nr0-1-2*np0+np ;  flag = .true.
    else 
       nr = -1 ; flag = .false.
    end if
    return
  end function i_p2r

end module mod_four2d


