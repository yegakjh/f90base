MODULE clock
  use mydef
  implicit none

  type myclock
     real :: cpu=zero
     real :: wall=zero
  end type myclock

  interface operator(-)
     module procedure subtract_myclock
  end interface

  interface operator(+)
     module procedure add_myclock
  end interface

!   function subtract_myclock(a,b) result(c)
!     type(myclock) :: c
!     type(myclock), intent(in) :: a,b

!   function add_myclock(a,b) result(c)
!     type(myclock) :: c
!     type(myclock), intent(in) :: a,b

!   subroutine time_check(clk, clk0)
!     type(myclock), intent(out) :: clk
!     type(myclock), optional, intent(in) :: clk0

contains
  function subtract_myclock(a,b) result(c)
    type(myclock) :: c
    type(myclock), intent(in) :: a,b
    
    c%cpu = a%cpu-b%cpu
    c%wall = a%wall-b%wall
    return
  end function subtract_myclock

  function add_myclock(a,b) result(c)
    type(myclock) :: c
    type(myclock), intent(in) :: a,b
    
    c%cpu = a%cpu+b%cpu
    c%wall = a%wall+b%wall
    return
  end function add_myclock

  subroutine time_check(clk, clk0)
    type(myclock), intent(out) :: clk
    type(myclock), optional, intent(in) :: clk0
    integer :: dt(8)
    
    call cpu_time(clk%cpu)
    call date_and_time(values=dt)
    clk%wall =   (60.*dt(5)+dt(6))*60.+dt(7)+0.001*dt(8)
    if (present(clk0)) then
       clk = clk-clk0
       if (clk%wall < 0.0) clk%wall= clk%wall+3600.
    end if
  end subroutine time_check
  
  subroutine time_print(clk, msg0)
    type(myclock), intent(in) :: clk
    character(len=*), intent(in), optional :: msg0

    real :: ratio
    character(len=*), parameter :: clk_fmt1= '(a20,a14,3(1x,e10.3))'
    character(len=*), parameter :: clk_fmt2= '(a14,3(1x,e10.3))'
    
    ratio = clk%cpu/clk%wall
    
    if (present(msg0)) then
       write(*,clk_fmt1) msg0,'(CPU,Wall,C/W):', clk%cpu, clk%wall, ratio
    else
       write(*,clk_fmt2) '(CPU,Wall,C/W):', clk%cpu, clk%wall, ratio
    end if
    return
  end subroutine time_print
end module clock

!
! Written by Juhyung Kim
!
module banded_matrix
  use mydef
  implicit none

  interface dgbsv
     module procedure dgbsv_1d, dgbsv_2d
  end interface
  interface zgbsv
     module procedure zgbsv_1d, zgbsv_2d
  end interface
  
contains
  subroutine cidx2_matrix_pack(aa, lda, n, lx)
    integer, intent(in) :: lda, n
    real(dp), dimension(0:lda,0:n), intent(out) :: aa
    real(dp), intent(in) :: lx

    integer :: i1, ibase
    real(dp) :: lxsq
    real(dp), dimension(0:n) :: cb
    
    cb = one; cb(0) = two; cb(n) = two
    lxsq = lx*lx*0.25_dp
    aa = zero
    
    if (lda < 6) stop 'cidx2_matrix_pack small lda'
    ibase = 4
    do i1=2, n
       aa(ibase+2,i1-2) = lxsq*cb(i1-2)/real(i1*(i1-1),dp)
       aa(ibase,i1) = -two*lxsq/real((i1*i1-1),dp)
       if (i1+2 <= n) aa(ibase-2,i1+2) = lxsq/real(i1*(i1+1),dp)
    enddo
    return
  end subroutine cidx2_matrix_pack

  subroutine dgbsv_pack(n, kl, ku, ab, ldab, aa, ldaa,info)
    ! this routine pack aa(ldb,n) into ab(ldab,n) which is to be passed
    ! to dgbsv as an argument.

    integer, intent(in) :: n, kl, ku, ldab, ldaa
    real(dp), dimension(ldaa, n), intent(in) :: aa
    real(dp), dimension(ldab, n), intent(out) :: ab
    integer, intent(out) :: info
    
    integer :: i1, i2, lb, ub, ibase
    
    ! check the size
    info = 0
    if (ldab < 2*kl+ku+1) info = 1
    if (ldaa < n) info =2
    
    ibase = kl+ku+1 
    if (info == 0) then
       do i2=1, n 
          if (i2 < ku+1) then
             lb=1
             ub = i2+kl
          else if (i2 <= n-ku) then
             lb = i2-ku
             ub = i2+kl
          else
             lb = i2-ku
             ub = n
          endif
          
          
          do i1=lb, ub
             ab(ibase+i1-i2,i2) = aa(i1,i2)
          enddo
       enddo
    endif
    return
  end subroutine dgbsv_pack
    
  SUBROUTINE DGBSV_2D( N, KL, KU, NRHS, AB, LDAB, B, LDB, INFO )
    INTEGER :: INFO, KL, KU, LDAB, LDB, N, NRHS
    INTEGER :: IPIV(n)
    real(dp) :: AB(LDAB,n), B(LDB,nrhs)
    
    ! Refer to dgbsv.f in BLAS

    ! test input parameter
    INFO = 0
    IF( N.LT.0 ) THEN
       INFO = -1
    ELSE IF( KL.LT.0 ) THEN
       INFO = -2
    ELSE IF( KU.LT.0 ) THEN
       INFO = -3
    ELSE IF( NRHS.LT.0 ) THEN
       INFO = -4
    ELSE IF( LDAB.LT.2*KL+KU+1 ) THEN
       INFO = -6
    ELSE IF( LDB.LT.MAX( N, 1 ) ) THEN
       INFO = -9
    END IF
    IF( INFO.NE.0 ) THEN
       CALL XERBLA( 'DGBSV ', -INFO )
       RETURN
    END IF

    ! Compute the LU factorization of the band matrix A.
    CALL DGBTRF( N, N, KL, KU, AB, LDAB, IPIV, INFO )

    IF( INFO.EQ.0 ) THEN
       ! Solve the system A*X = B, overwriting B with X.
       CALL DGBTRS( 'No transpose', N, KL, KU, NRHS, AB, LDAB, IPIV, &
            B, LDB, INFO )
    END IF
    RETURN
  end SUBROUTINE DGBSV_2D

  SUBROUTINE DGBSV_1D( N, KL, KU, AB, LDAB, B, LDB, INFO )
    !
    ! 1D version of DGBSV_2D (nrhs = 1)
    INTEGER :: INFO, KL, KU, LDAB, LDB, N, nrhs
    
    INTEGER :: IPIV(n)
    real(dp) :: AB(LDAB,n), B(LDB)
    
    nrhs = 1
    ! test input parameter
    INFO = 0
    IF( N.LT.0 ) THEN
       INFO = -1
    ELSE IF( KL.LT.0 ) THEN
       INFO = -2
    ELSE IF( KU.LT.0 ) THEN
       INFO = -3
    ELSE IF( NRHS.LT.0 ) THEN
       INFO = -4
    ELSE IF( LDAB.LT.2*KL+KU+1 ) THEN
       INFO = -6
    ELSE IF( LDB.LT.MAX( N, 1 ) ) THEN
       INFO = -9
    END IF
    IF( INFO.NE.0 ) THEN
       CALL XERBLA( 'DGBSV ', -INFO )
       RETURN
    END IF

    ! Compute the LU factorization of the band matrix A.
    CALL DGBTRF( N, N, KL, KU, AB, LDAB, IPIV, INFO )

    IF( INFO.EQ.0 ) THEN
       ! Solve the system A*X = B, overwriting B with X.
       CALL DGBTRS( 'No transpose', N, KL, KU, NRHS, AB, LDAB, IPIV, &
            B, LDB, INFO )
    END IF
    RETURN
  end SUBROUTINE DGBSV_1D

  subroutine zgbsv_pack(n, kl, ku, ab, ldab, aa, ldaa, info)
    ! this routine pack aa(ldb,n) into ab(ldab,n) which is to be passed
    ! to zgbsv as an argument.

    integer, intent(in) :: n, kl, ku, ldab, ldaa
    complex(dpc), dimension(ldaa, n), intent(in) :: aa
    complex(dpc), dimension(ldab, n), intent(out) :: ab
    integer, intent(out) :: info
    
    integer :: i1, i2, lb, ub, ibase
    
    ! check the size
    info = 0
    if (ldab < 2*kl+ku+1) info = 1
    if (ldaa < n) info =2
    
    ibase = kl+ku+1 
    if (info == 0) then
       do i2=1, n 
          if (i2 < ku+1) then
             lb=1
             ub = i2+kl
          else if (i2 <= n-ku) then
             lb = i2-ku
             ub = i2+kl
          else
             lb = i2-ku
             ub = n
          endif
          
          
          do i1=lb, ub
             ab(ibase+i1-i2,i2) = aa(i1,i2)
          enddo
       enddo
    endif
    return
  end subroutine zgbsv_pack
    
  SUBROUTINE ZGBSV_2D( N, KL, KU, NRHS, AB, LDAB, B, LDB, INFO )
    INTEGER :: INFO, KL, KU, LDAB, LDB, N, NRHS
    INTEGER :: IPIV(n)
    complex(dpc) :: AB(LDAB,n), B(LDB,nrhs)
    
    ! Refer to zgbsv.f in BLAS

    ! test input parameter
    INFO = 0
    IF( N.LT.0 ) THEN
       INFO = -1
    ELSE IF( KL.LT.0 ) THEN
       INFO = -2
    ELSE IF( KU.LT.0 ) THEN
       INFO = -3
    ELSE IF( NRHS.LT.0 ) THEN
       INFO = -4
    ELSE IF( LDAB.LT.2*KL+KU+1 ) THEN
       INFO = -6
    ELSE IF( LDB.LT.MAX( N, 1 ) ) THEN
       INFO = -9
    END IF
    IF( INFO.NE.0 ) THEN
       CALL XERBLA( 'ZGBSV ', -INFO )
       RETURN
    END IF

    ! Compute the LU factorization of the band matrix A.
    CALL ZGBTRF( N, N, KL, KU, AB, LDAB, IPIV, INFO )

    IF( INFO.EQ.0 ) THEN
       ! Solve the system A*X = B, overwriting B with X.
       CALL ZGBTRS( 'No transpose', N, KL, KU, NRHS, AB, LDAB, IPIV, &
            B, LDB, INFO )
    END IF
    RETURN
  end SUBROUTINE ZGBSV_2D

  SUBROUTINE ZGBSV_1D( N, KL, KU, AB, LDAB, B0, LDB, INFO )
    !
    ! 1D version of ZGBSV_2D 
    integer, parameter :: nrhs=1
    INTEGER :: INFO, KL, KU, LDAB, LDB, N
    
    INTEGER :: IPIV(n)
    complex(dpc) :: B(LDB,nrhs)
    complex(dpc), intent(inout) :: B0(LDB), ab(ldab,n)
    
    B(:,1) = B0

    ! test input parameter
    INFO = 0
    IF( N.LT.0 ) THEN
       INFO = -1
    ELSE IF( KL.LT.0 ) THEN
       INFO = -2
    ELSE IF( KU.LT.0 ) THEN
       INFO = -3
    ELSE IF( NRHS.LT.0 ) THEN
       INFO = -4
    ELSE IF( LDAB.LT.2*KL+KU+1 ) THEN
       INFO = -6
    ELSE IF( LDB.LT.MAX( N, 1 ) ) THEN
       INFO = -9
    END IF
    IF( INFO.NE.0 ) THEN
       CALL XERBLA( 'ZGBSV ', -INFO )
       RETURN
    END IF

    ! Compute the LU factorization of the band matrix A.
    CALL ZGBTRF( N, N, KL, KU, AB, LDAB, IPIV, INFO )

    IF( INFO.EQ.0 ) THEN
       ! Solve the system A*X = B, overwriting B with X.
       CALL ZGBTRS( 'No transpose', N, KL, KU, NRHS, AB, LDAB, IPIV, &
            B, LDB, INFO )
       B0 = B(:,1)
    END IF
    RETURN
  end SUBROUTINE ZGBSV_1D

end module banded_matrix


module myutil
  use mydef
  implicit none

  interface print_mxmn
     module procedure print2d_r_mxmn, print2d_c_mxmn
  end interface

  interface solve_lin_eq
     module procedure solve_lin_eq_1d_r, solve_lin_eq_1d_c, &
          solve_lin_eq_2d_r, solve_lin_eq_2d_c
  end interface

  ! function bb=cheb_int(aa,nn,lx)
  interface cheb_int
     module procedure cheb_int_real, cheb_int_complex
  end interface
     
  interface error_check
     module procedure error_check_r_sp, error_check_r_dp, error_check_r_dp_2d
     module procedure error_check_c_sp, error_check_c_dp
  end interface

  interface inverse
     module procedure inverse_real, inverse_real_v1, inverse_cmplx
  end interface

  interface sort_map
     module procedure sort_map_1d_r, sort_map_1d_c, &
          sort_map_2d_r, sort_map_2d_c, &
          sort_map_3d_r, sort_map_3d_c
  end interface

  interface chck_nan
     ! check NaN and print out the index of first NaN and return .true.
     module procedure chck_nan_1d, chck_nan_2d, chck_nan_3d
  end interface
   
  interface add_identity
     module procedure add_identity_real, add_identity_complex, &
          add_identity_real_rec, add_identity_complex_rec
  end interface

  interface add_cxmul_matrix
     module procedure add_cxmul_matrix_complex, &
          add_cxmul_matrix_complex_rec
  end interface

contains
!A
  ! add chebyshev matrix for x multiplciation
  subroutine add_cxmul_matrix_complex(xmul, n, lx)
    integer, intent(in) :: n
    complex(dpc), intent(in) :: lx
    complex(dpc), dimension(0:n, 0:n), intent(inout) :: xmul
    integer :: i1
    complex(dpc) :: hlx
    
    hlx = half*lx
    xmul(0,1) = xmul(0,1)+hlx
    xmul(1,0) = xmul(1,0)+lx
    do i1=1,n-1
       xmul(i1, i1+1) = xmul(i1,i1+1)+hlx
       xmul(i1+1, i1) = xmul(i1+1,i1)+hlx
    enddo
    return
  end subroutine add_cxmul_matrix_complex

  subroutine add_cxmul_matrix_complex_rec(xmul, n1, n2, lx)
    integer, intent(in) :: n1, n2
    complex(dpc), intent(in) :: lx
    complex(dpc), dimension(0:n1, 0:n2), intent(inout) :: xmul
    integer :: i1,nn
    complex(dpc) :: hlx
    
    hlx = half*lx
    xmul(0,1) = xmul(0,1)+hlx
    xmul(1,0) = xmul(1,0)+lx

    if (n1 < n2) then 
       nn = n1 
    else 
       nn = n2
    end if

    do i1=1,nn-1
       xmul(i1, i1+1) = xmul(i1, i1+1) + hlx
       xmul(i1+1, i1) = xmul(i1+1, i1) + hlx
    enddo
    
    if (n1 > n2) xmul(nn+1,nn) = xmul(nn+1,nn)+hlx
    if (n1 < n2) xmul(nn,nn+1) = xmul(nn, nn+1)+hlx

    return
  end subroutine add_cxmul_matrix_complex_rec
  
  subroutine add_identity_real(mtx, n, a)
    ! mtx = mtx+ a*i_m
  integer, intent(in) :: n
  real(dp), intent(in) :: a
  real(dp), dimension(n,n) :: mtx
 
  integer :: i1
  do i1=1, n
     mtx(i1,i1) = mtx(i1,i1)+a
  enddo
  return
end subroutine add_identity_real

subroutine add_identity_real_rec(mtx, n1,n2, a)
  ! mtx = mtx+ a*i_m
  integer, intent(in) :: n1,n2
  real(dp), intent(in) :: a
  real(dp), dimension(n1,n2) :: mtx
 
  integer :: i1,nn
  nn = n1
  if (n1 > n2) nn = n2
  do i1=1, nn
     mtx(i1,i1) = mtx(i1,i1)+a
  enddo
  return
end subroutine add_identity_real_rec

subroutine add_identity_complex(mtx, n, a)
  ! mtx = mtx+ a*i_m
  integer, intent(in) :: n
  complex(dpc), intent(in) :: a
  complex(dpc), dimension(n,n) :: mtx
  
  integer :: i1
  do i1=1, n
     mtx(i1,i1) = mtx(i1,i1)+a
  enddo
  return
end subroutine add_identity_complex

subroutine add_identity_complex_rec(mtx, n1,n2, a)
  ! mtx = mtx+ a*i_m
  integer, intent(in) :: n1,n2
  complex(dpc), intent(in) :: a
  complex(dpc), dimension(n1,n2) :: mtx
  
  integer :: i1,nn

  if (n1 > n2) then 
     nn = n2 
  else 
     nn = n1
  end if

  do i1=1, nn
     mtx(i1,i1) = mtx(i1,i1)+a
  enddo
  return
end subroutine add_identity_complex_rec

!B
!BI
  subroutine bin_get(nn,xmin, xmax,bin)
    implicit none
    
    !   xmin < x < xmax
    ! bin = [x0, x1, x2, ... , xnn]
    
    integer:: nn
    real(dp), intent(in) :: xmin, xmax
    real(dp), dimension(0:nn), intent(out) :: bin

    real(dp) :: dx
    integer :: i1
    
    dx = (xmax-xmin)/real(nn,dp)
    do i1=0, nn
       bin(i1) = xmin+dx*real(i1,dp)
    enddo
    return
  end subroutine bin_get
  
  subroutine binning1d(nn, bin,res, n1, x, y)
    ! if bin(n-1) < x < bin(n) then res(n) = y+res(n)
    integer, intent(in) :: nn, n1
    real(dp), dimension(0:nn), intent(in) :: bin
    real(dp), dimension(n1), intent(in) :: x
    real(dp), dimension(n1), intent(in), optional :: y
    real(dp), dimension(nn), intent(out) :: res
    
    real(dp) :: xmin, dx
    integer :: i1,i2,k
    logical :: flag = .FALSE.

    if (present(y))        flag = .TRUE.
    res = zero
    xmin = bin(0)
    dx = bin(1)-bin(0)

    if (flag) then
       ! y data is given
       do i1=1,n1
          k = floor((x(i1)-xmin)/dx)+1
          if (k >= 1 .and. k <= nn) res(k) = res(k) + y(i1)
       END do
    else
       do i1=1,n1
          k = floor((x(i1)-xmin)/dx)+1
          if (k >= 1.and.k <= nn) res(k) = res(k) + 1.0
       END do
    end if

    return
  end subroutine binning1d
  
  subroutine binning2d(nn, bin,res, n1, n2, x, y)
    ! if bin(n-1) < x < bin(n) then res(n) = y+res(n)
    integer, intent(in) :: nn, n1, n2
    real(dp), dimension(0:nn), intent(in) :: bin
    real(dp), dimension(n1, n2), intent(in) :: x
    real(dp), dimension(n1, n2), intent(in), optional :: y
    real(dp), dimension(nn), intent(out) :: res
    
    real(dp) :: xmin, dx, hdx
    integer :: i1,i2,k
    logical :: flag=.FALSE.
    
    if (present(y)) flag = .TRUE.
    
    res = zero
    xmin = bin(0)
    dx = bin(1)-bin(0)

    if (flag) then
       ! y data is given
       do i2=1,n2
          do i1=1,n1
             k = floor((x(i1,i2)-xmin)/dx)+1
             if (k >= 1.and.k <= nn) res(k) = res(k) + y(i1,i2)
          END do
       end do
    else
       do i2=1,n2
          do i1=1,n1
             k = floor((x(i1,i2)-xmin)/dx)+1
             if (k >= 1 .and.k <= nn) res(k) = res(k) + 1.0
          END do
       end do
    end if
    
    return
  end subroutine binning2d

  ! binning with two measure
  subroutine binning2d_2m(nn1, bin1,nn2, bin2,res,  n1, n2, x1, x2, y)
    ! if bin(n-1) < x < bin(n) then res(n) = y+res(n)
    integer, intent(in) :: nn1, nn2, n1, n2
    real(dp), dimension(0:nn1), intent(in) :: bin1
    real(dp), dimension(0:nn2), intent(in) :: bin2
    real(dp), dimension(n1, n2), intent(in) :: x1,x2
    real(dp), dimension(n1, n2), intent(in), optional :: y
    real(dp), dimension(nn1,nn2), intent(out) :: res
    
    real(dp) :: xmin1, xmin2, dx1, dx2
    integer :: i1,i2,k1, k2
    logical :: flag=.FALSE.

    if (present(y)) flag = .TRUE.
    
    res = zero
    xmin1 = bin1(0);  xmin2 = bin2(0)
    dx1 = bin1(1)-bin1(0)
    dx2 = bin2(1)-bin2(0)
    
    do i2=1,n2
       do i1=1,n1
          k1 = floor((x1(i1,i2)-xmin1)/dx1)+1
          k2 = floor((x2(i1,i2)-xmin2)/dx2)+1
          if (k1 >= 1 .and. k1 <= nn1 .and. k2 >= 1 .and. k2 <= nn2) then
             if (flag) then
                res(k1,k2) = res(k1,k2) + y(i1,i2)
             else
                res(k1,k2) = res(k1,k2) + 1.0
             end if
          end if
       END do
    end do
    return
  end subroutine binning2d_2m

!CH
function chck_nan_1d(n1, arr) result(flg)
  real(dp), dimension(n1), intent(in) :: arr
  integer, intent(in) :: n1
  logical :: flg
  integer :: i1
  character(*), parameter :: fmt='(a,(1x,i3), a, (1x,i3))'

  flg = .false.  

  do i1 = 1, n1
     if (isnan(arr(i1))) then
        write(*,fmt), ' NaN, at the index ', i1, ' of ', n1
        flg = .true.
        exit
     end if
  end do
  return
end function chck_nan_1d

function chck_nan_2d(n1, n2,  arr) result(flg)
  real(dp), dimension(n1,n2), intent(in) :: arr
  integer, intent(in) :: n1, n2
  logical :: flg
  integer :: i1, i2
  character(*), parameter :: fmt='(a,2(1x,i3), a, 2(1x,i3))'

  flg = .false.  
  main : do i2 = 1, n2
     do i1 = 1, n1
        if (isnan(arr(i1,i2))) then
           write(*,fmt), ' NaN, at the index ', i1,i2, ' of ', n1, n2
           flg = .true.
           exit main
        end if
     end do
  end do main
  
  return
end function chck_nan_2d

function chck_nan_3d(n1, n2, n3, arr) result(flg)
  real(dp), dimension(n1,n2,n3), intent(in) :: arr
  integer, intent(in) :: n1, n2, n3
  logical :: flg
  integer :: i1, i2, i3
  character(*), parameter :: fmt='(a,3(1x,i3), a, 3(1x,i3))'
  flg = .false.
  
  main : do i3 = 1, n3
     do i2 = 1, n2
        do i1 = 1, n1
           if (isnan(arr(i1,i2,i3))) then
              write(*,fmt), ' NaN, at the index ', i1,i2,i3, ' of ', n1,n2,n3
              flg = .true.
              exit main
           end if
        end do
     end do
  end do main
  
  return
end function chck_nan_3d

!CH
function cheb_int_real(aa,nn,lx) result(bb)
  integer, intent(in) :: nn
  real(dp), intent(in) :: lx
  real(dp), dimension(0:nn), intent(in) :: aa
  real(dp), dimension(0:nn) :: bb
  
  real(dp) :: lxsq
  integer :: i1
  
  lxsq = lx*lx*0.25_dp

  bb(0) = zero; bb(1) = zero  
  bb(2) = lxsq*(aa(0)-(two*aa(2)+half*aa(4))/3.0_dp)
  do i1=3, nn-2
     bb(i1) = aa(i1-2)/real(i1*(i1-1),dp) &
          -two*aa(i1)/real((i1*i1-1),dp) &
          +aa(i1+2)/real(i1*(i1+1),dp)
     bb(i1) = lxsq*bb(i1)
  end do
  bb(nn-1) = aa(nn-3)/real((nn-1)*(nn-2),dp) &
       -two*aa(nn-1)/real(((nn-1)*(nn-1)-1),dp)
  bb(nn-1) = lxsq*bb(nn-1)
  
  bb(nn) = aa(nn-2)/real(nn*(nn-1),dp) &
       -two*aa(nn)/real((nn*nn-1),dp)
  bb(nn) = lxsq*bb(nn)
  
  return
end function cheb_int_real

function cheb_int_complex(aa,nn,lx) result(bb)
  integer, intent(in) :: nn
  real(dp), intent(in) :: lx
  complex(dpc), dimension(0:nn), intent(in) :: aa
  complex(dpc), dimension(0:nn) :: bb
  
  real(dp) :: lxsq
  integer :: i1
  
  lxsq = lx*lx*0.25_dp
  bb(0) = zero; bb(1) = zero
  bb(2) = lxsq*(aa(0)-(two*aa(2)+half*aa(4))/3.0_dp)
  do i1=3, nn-2
     bb(i1) = aa(i1-2)/real(i1*(i1-1),dp) &
          -two*aa(i1)/real((i1*i1-1),dp) &
          +aa(i1+2)/real(i1*(i1+1),dp)
     bb(i1) = lxsq*bb(i1)
  end do
    
  bb(nn-1) = aa(nn-3)/real((nn-1)*(nn-2),dp) &
         -two*aa(nn-1)/real(((nn-1)*(nn-1)-1),dp)
  bb(nn-1) = lxsq*bb(nn-1)
  
  bb(nn) = aa(nn-2)/real(nn*(nn-1),dp) &
       -two*aa(nn)/real((nn*nn-1),dp)
  bb(nn) = lxsq*bb(nn)
  
  return
end function cheb_int_complex

!CR
! create (n,n) identity matrix
subroutine create_identity_matrix(in_m, n)
  integer, intent(in) :: n
  real(dp), dimension(n,n), intent(out) :: in_m
  integer :: i1
  
  in_m = zero
  do i1=1,n
     in_m(i1,i1) = one
  end do
end subroutine create_identity_matrix

! create (n+1, n+1) Chebyshev n-th derivative matrix
subroutine create_cdx_matrix(mtx, n, lx, nth)
  integer :: n, nth, i1,i2
  real(dp) :: lx
  real(dp), dimension(0:n,0:n), intent(out) :: mtx
  real(dp), dimension(0:n) :: cb
  intent(in) :: n, nth, lx
  
  cb = two/lx; cb(0) = one/lx; cb(n) = one/lx
  mtx = zero
  
  do i1=0,n-1
     do i2=i1+1, n, 2
        mtx(i1,i2) = cb(i1)*real(i2,dp)
     end do
  end do
  if (nth > 1) then
     do i1=2, nth
        mtx = matmul(mtx,mtx)
     end do
  end if
  return
end subroutine create_cdx_matrix

! create (n+1, n+1) Chebyshev inverse n-th derivative matrix
! or call n-th order indefinite integral matrix
subroutine create_cidx2_matrix(mtx, n,lx)
  integer :: n, i1
  real(dp) :: lx, lxsq
  real(dp), dimension(0:n,0:n), intent(out) :: mtx
  real(dp), dimension(0:n) :: cb
  intent(in) :: n, lx
  
  cb = one; cb(0) = two; cb(n) = two
  lxsq = lx*lx*0.25_dp
  mtx = zero
  
  do i1=2, n
     mtx(i1,i1-2) = lxsq*cb(i1-2)/real(i1*(i1-1),dp)
     mtx(i1,i1) = -2.0_dp*lxsq/real((i1*i1-1),dp)
     if (i1+2 <= n) mtx(i1,i1+2) = lxsq/real(i1*(i1+1),dp)
  enddo
  return
end subroutine create_cidx2_matrix

subroutine create_cabsx_matrix(absx, n,lx)
  integer, intent(in) :: n
  real(dp), intent(in) :: lx
  real(dp), dimension(0:n,0:n), intent(out) :: absx
  
  real(dp), dimension(0:n) :: cb
  integer :: i1,i2
  
  cb = -lx*two/pi
  cb(0) = cb(0)*half
  
  absx = zero
  forall (i1=0:n, i2=0:n, mod(i1+i2,2) == 0)
     absx(i1,i2) = (-1)**(mod((i1+i2)/2,2))/((i1+i2-one)*(i1+i2+one))+ &
          (-1)**(mod((i1-i2)/2,2))/((i1-i2-one)*(i1-i2+one))
     absx(i1,i2) = absx(i1,i2)*cb(i1)
  end forall
  return
end subroutine create_cabsx_matrix

! create Chebyshev matrix for x multiplication
subroutine create_cxmul_matrix(xmul, n, lx)
  integer, intent(in) :: n
  real(dp), intent(in) :: lx
  real(dp), dimension(0:n, 0:n), intent(out) :: xmul
  integer :: i1
  
  xmul = zero
  xmul(0,1) = half*lx
  xmul(1,0) = lx
  do i1=1,n-1
     xmul(i1, i1+1) = half*Lx
     xmul(i1+1, i1) = half*Lx
  enddo
  return
end subroutine create_cxmul_matrix

!D
!E
subroutine error_exit(msg)
  character(len=*) :: msg
  
  write(*,*) msg
  stop
end subroutine error_exit

subroutine error_check_r_sp(a,b,thres,msg,errout)
  real(sp), intent(in) :: a, b, thres
  real(sp), optional :: errout
  character(len=*) :: msg
  real(sp) :: large, small,err
  
  if (abs(a) > abs(b)) then
     large = a; small = b
  else
     large = b; small = a
  end if
  
  if (large == 0.0) then
     err = zero
  else if (small == 0.0) then
     err = abs(large)
  else
     err = abs((a-b)/b)
  end if
  
  if (err > thres) then
     write(*,'(a,e10.3,1x,a,1x,e10.3)') msg, err, ' > ', thres
  end if

  if (present(errout)) errout = err

  return
end subroutine error_check_r_sp

subroutine error_check_r_dp(a,b,thres,msg,errout)
  real(dp), intent(in) :: a, b, thres
  real(dp), optional :: errout
  character(len=*) :: msg
  real(dp) :: large, small,err
  
  if (abs(a) > abs(b)) then
     large = a; small = b
  else
     large = b; small = a
  end if

  if (large == 0.0) then
     err = zero
  else if (small == 0.0) then
     err = abs(large)
  else
     err = abs((a-b)/b)
  end if

  if (err > thres) then
     write(*,'(a,e10.3,1x,a,1x,e10.3)') msg, err, ' > ', thres
  end if
  if (present(errout)) errout = err
  return
end subroutine error_check_r_dp

subroutine error_check_c_sp(a,b,thres,msg,errout)
  complex(spc), intent(in) :: a, b
  real(sp), intent(in) :: thres
  real(sp), optional :: errout
  character(len=*) :: msg
  real(dp) ::err, large, small
  
  if (abs(a) > abs(b)) then
     large = abs(a); small = abs(b)
  else
     large = abs(b); small = abs(a)
  end if
  
  if (large == 0.0) then
     err = zero
  else if (small == 0.0) then
     err = abs(large)
  else
     err = abs((a-b)/b)
  end if
  
  if (err > thres) then
     write(*,'(a,e10.3,1x,a,1x,e10.3)') msg, err, ' > ', thres
  end if

  if (present(errout)) errout = err

  return
end subroutine error_check_c_sp

subroutine error_check_c_dp(a,b,thres,msg,errout)
  complex(dpc), intent(in) :: a, b
  real(dp), intent(in) :: thres
  real(dp), optional :: errout
  character(len=*) :: msg
  real(dp) ::err, large, small
  
  if (abs(a) > abs(b)) then
     large = abs(a); small = abs(b)
  else
     large = abs(b); small = abs(a)
  end if
  
  if (large == 0.0) then
     err = zero
  else if (small == 0.0) then
     err = abs(large)
  else
     err = abs((a-b)/b)
  end if
  
  if (err > thres) then
     write(*,'(a,e10.3,1x,a,1x,e10.3)') msg, err, ' > ', thres
  end if

  if (present(errout)) errout = err

  return
end subroutine error_check_c_dp


  
  subroutine error_check_r_dp_2d(n1, n2, a,b,thres,msg)
    integer, intent(in) :: n1, n2
    real(dp), dimension(n1,n2), intent(in) :: a, b
    real(dp), intent(in) :: thres
    character(len=*) :: msg

    integer :: i1, i2
    real(dp) :: large, small,err
    character(len=*), parameter :: fmt='(a,2(i3,1x),e10.3,1x,a,1x,e10.3)'
    
    outer : do i2 = 1, n2
       inner : do i1 = 1, n1
          if (abs(a(i1,i2)) > abs(b(i1,i2))) then
             large = a(i1,i2); small = b(i1,i2)
          else
             large = b(i1,i2); small = a(i1,i2)
          end if
          
          if (large == 0.0) cycle inner
          if (small == 0.0 .and. abs(large) > thres) then
             write(*,fmt) msg, i1, i2, large, ' > ', thres
             cycle inner
          end if
          err = abs((large-small)/large)
          if (err > thres) then
             write(*,fmt) msg, i1, i2, err, ' > ', thres
          end if
       enddo inner
    enddo outer
    return
  end subroutine error_check_r_dp_2d
          
!F
!G
!H
!IN
  function inv_helmholtz(na,aa,bb) result(cc)
    integer, intent(in) :: na
    real(dp), dimension(na,na),intent(in) :: aa
    complex(dpc), dimension(na),intent(in) :: bb
    complex(dpc), dimension(na)  :: cc,cc0

    cc0 = bb;     cc0(na-1:na) = zero
    cc = matmul(aa,cc0)
    return
    
    return
  end function inv_helmholtz

  subroutine inverse_real(na, A, iA)
    integer, intent(in) :: na
    real(dp), dimension(na,na), intent(in) :: a
    real(dp), dimension(na,na), intent(out) :: ia
    real(dp), dimension(na,na) :: mat
    integer, dimension(na) :: ipiv
    integer  :: lwork, info
    real(dp), dimension(na*64) :: work
    
    lwork = na*64
    ia  = a
    call dgetrf(na,na, ia, na,ipiv, info)
    if (info.ne.0)        call print_err(info, 'in dgetrf') 
  
    call dgetri(na,ia,na, ipiv, work, lwork,info)
    if (info.ne.0)        call print_err(info, 'in dgetri')

    
    mat = matmul(a,ia)
    if (maxval(mat).gt. 1.01) then
       call print_mxmn(6, mat, "|ia*a| in inverse = ")
    end if
    return
  end subroutine inverse_real

  subroutine inverse_real_v1(na, nb, A, iA)
    integer, intent(in) :: na,nb
    real(dp), dimension(na,na), intent(in) :: a
    real(dp), dimension(na,na), intent(out) :: ia
    real(dp), dimension(na,na) :: mat
    integer, dimension(na) :: ipiv
    integer  :: lwork, info
    real(dp), dimension(na*64) :: work
    
    lwork = na*64
    ia(:,1:na) = a(:,2:na+1)
    call dgetrf(na,nb, ia, na,ipiv, info)
    if (info.ne.0)        call print_err(info, 'in dgetrf') 
  
    call dgetri(na,ia,na, ipiv, work, lwork,info)
    if (info.ne.0)        call print_err(info, 'in dgetri')

    mat = matmul(a(:,2:na+1),ia(1:na,:))
    if (maxval(mat).gt. 1.01) then
       call print_mxmn(6,mat, "|ia*a| in inverse = ")
    end if
    return
  end subroutine inverse_real_v1

  subroutine inverse_cmplx(na, A, iA)
    integer, intent(in) :: na
    complex(dpc), dimension(na,na), intent(in) :: a
    complex(dpc), dimension(na,na), intent(out) :: ia
    complex(dpc), dimension(na,na) :: mat
    integer, dimension(na) :: ipiv
    integer  :: lwork, info
    complex(dpc), dimension(na*64) :: work
    
    lwork = na*64
    ia  = a
    call zgetrf(na,na, ia, na,ipiv, info)
    if (info.ne.0)        call print_err(info, 'in zgetrf') 
  
    call zgetri(na,ia,na, ipiv, work, lwork,info)
    if (info.ne.0)        call print_err(info, 'in zgetri')

    if (maxval(abs(mat)).gt. 1.01) then
       call print_mxmn(7,abs(mat), "|ia*a| in inverse = ")
    end if

    return
  end subroutine inverse_cmplx

!IO
  subroutine ios_check(filename,ios)
    character(len=*) :: filename
    integer :: ios
    
    if (ios.gt.0) then
       write(*,*) "Read Error ",filename, " ios = ", ios
       stop
    end if
    return
  end subroutine ios_check
!J
!K
!L
!M
!MO
  subroutine moment(n1, n2, data, out)
    implicit none
    integer, intent(in) :: n1, n2
    real(dp), dimension(n1,n2), intent(in) :: data
    real(dp), dimension(6), intent(out) :: out
    ! given an array of data, this routine returns 
    ! ave : mean 
    ! adev : average deviation
    ! sdev : standard deviation
    ! var : variance
    ! skew : sskewness
    ! curt : kurtosis
    
    integer :: n
    real(dp) :: ep
    real(dp), dimension(n1,n2) :: p, s
    
    n = n1*n2
    if ( n<= 1) call print_err(n, 'moment : n must be at least 2')
    
    ! average
    out(1) = sum(data)/n 

    s = data(:,:)-out(1)
    ep = sum(s)

    ! adev
    out(2) = sum(abs(s))/n

    ! var
    p = s*s;     out(4) = sum(p)

    ! skew
    p = s*p;     out(5) = sum(p)
    
    ! kurtosis
    p = p*s;     out(6) = sum(p)
    
    out(4) = (out(4)-ep**2/n)/(n-1)
    out(3) = sqrt(out(4))
    
    if (out(4) /= 0.0) then
       out(5) = out(5)/(n*out(3)**3)
       out(6) = out(6)/(n*out(4)**2)-3.0_sp
    else
       call print_err(11, 'moment : no out(5) or kurtosis when zero variance')
    end if
  end subroutine moment
!N
!O
!P
!PR
  subroutine print_err(info, lbl)
    integer :: info
    character(*) :: lbl
    character(*), parameter :: head='(a,1x,i3,1x,a)'
    write(*,head) 'Error code = ', info, lbl
    stop
  end subroutine print_err

  subroutine print2d_r_mxmn(nunit, data, str)
    integer :: nunit
    real(dp), dimension(:,:) :: data
    character(*) :: str
    !local
    real(dp) :: mx, mn
    integer :: imx(2), imn(2), id(2)

    id = shape(data)
    mx = maxval(data); imx = maxloc(data) 
    mn = minval(data); imn = minloc(data) 
    write(nunit,'(a, e10.3, a,2(i4,1x), a, 2(i4,1x),a)')  & 
         'Max = ' , mx, ' at ', imx, ' of ', id, str
    write(nunit,'(a, e10.3, a,2(i4,1x), a, 2(i4,1x),a)')  & 
         'Min = ' , mn, ' at ', imn, ' of ', id, str
    return
  end subroutine print2d_r_mxmn

  subroutine print2d_c_mxmn(nunit, data,str)
    integer :: nunit
    complex(dpc), dimension(:,:) :: data
    character(*) :: str
    !local
    real(dp) :: mx, mn
    integer :: imx(2), imn(2), id(2)

    id = shape(data)
    mx = maxval(abs(data)); imx = maxloc(abs(data)) 
    mn = minval(abs(data)); imn = minloc(abs(data)) 
    write(nunit,'(a, e10.3, a,2(i4,1x), a, 2(i4,1x),a)')  & 
         'Max = ' , mx, ' at ', imx, ' of ', id,str
    write(nunit,'(a, e10.3, a,2(i4,1x), a, 2(i4,1x),a)')  & 
         'Min = ' , mn, ' at ', imn, ' of ', id,str
    return
  end subroutine print2d_c_mxmn
!Q
!R
!S
!SO
  subroutine solve_lin_eq_1d_r(na,aa,bb,cc)
    ! aa * cc = bb, get cc
    integer, intent(in) :: na
    real(dp), dimension(na,na),intent(inout) :: aa
    real(dp), dimension(na),intent(in) :: bb
    real(dp), dimension(na),intent(out) :: cc
    integer, dimension(na) :: ipiv
    integer :: info

    !A*x = B
    call dgetrf(na,na, aa, na,ipiv, info)
    if (info.ne.0)        call print_err(info, 'in dgetrf')
    
    cc = bb
    call dgetrs('N',na,1, aa,na,ipiv,cc,na,info)
    if (info.ne.0)        call print_err(info, 'in dgetrs')
    return
  end subroutine solve_lin_eq_1d_r
  
  subroutine solve_lin_eq_1d_c(na,aa,bb,cc)
    ! aa is replaced
    integer, intent(in) :: na
    complex(dpc), dimension(na,na),intent(inout) :: aa
    complex(dpc), dimension(na),intent(in) :: bb
    complex(dpc), dimension(na),intent(out) :: cc
    integer, dimension(na) :: ipiv
    integer :: info

    !A*x = B
    call zgetrf(na,na, aa, na,ipiv, info)
    if (info.ne.0)        call print_err(info, 'in zgetrf')
 
    cc = bb
    call zgetrs('N',na,1, aa,na,ipiv,cc,na,info)
    if (info.ne.0)        call print_err(info, 'in zgetrs')
    return
  end subroutine solve_lin_eq_1d_c

  subroutine solve_lin_eq_2d_r(na,nb,aa,bb,cc)
    integer, intent(in) :: na, nb
    real(dp), dimension(na,na),intent(in) :: aa
    real(dp), dimension(na,nb),intent(in) :: bb
    real(dp), dimension(na,nb),intent(out) :: cc
    real(dp), dimension(na,na) :: aa1
    integer, dimension(na) :: ipiv
    integer :: info

    aa1 = aa
    !A*x = B
    call dgetrf(na,na, aa1, na,ipiv, info)
    if (info.ne.0)        call print_err(info, 'in dgetrf')
 
    cc = bb
    call dgetrs('N',na,nb, aa1,na,ipiv,cc,na,info)
    if (info.ne.0)        call print_err(info, 'in dgetrs')
    return
  end subroutine solve_lin_eq_2d_r

  subroutine solve_lin_eq_2d_c(na,nb,aa,bb,cc)
    integer, intent(in) :: na, nb
    complex(dpc), dimension(na,na),intent(in) :: aa
    complex(dpc), dimension(na,nb),intent(in) :: bb
    complex(dpc), dimension(na,nb),intent(out) :: cc
    complex(dpc), dimension(na,na) :: aa1
    integer, dimension(na) :: ipiv
    integer :: info

    aa1 = aa
    !A*x = B
    call zgetrf(na,na, aa1, na,ipiv, info)
    if (info.ne.0)        call print_err(info, 'in zgetrf')
 
    cc = bb
    call zgetrs('N',na,nb, aa1,na,ipiv,cc,na,info)
    if (info.ne.0)        call print_err(info, 'in zgetrs')
    return
  end subroutine solve_lin_eq_2d_c

  subroutine sort_map_1d_c(nm, x1, map)
    integer, intent(in) :: nm
    complex(dpc), dimension(nm), intent(inout) :: x1
    integer, dimension(nm), intent(in) :: map

    complex(dpc), dimension(nm) :: ctmp
    integer :: i1

    ctmp = x1
    do i1=1, nm
       x1(i1) = ctmp(map(i1))
    end do
  end subroutine sort_map_1d_c

  subroutine sort_map_1d_r(nm, x1, map)
    integer, intent(in) :: nm
    real(dp), dimension(nm), intent(inout) :: x1
    integer, dimension(nm), intent(in) :: map
    
    real(dp), dimension(nm) :: ctmp
    integer :: i1
    
    ctmp = x1
    do i1=1, nm
       x1(i1) = ctmp(map(i1))
    end do
  end subroutine sort_map_1d_r
  
  subroutine sort_map_2d_c(n1, nm, x1, map)
    integer, intent(in) :: n1, nm
    complex(dpc), dimension(n1,nm), intent(inout) :: x1
    integer, dimension(nm), intent(in) :: map
    
    complex(dpc), dimension(n1,nm) :: ctmp
    integer :: i1

    ctmp = x1
    do i1=1, nm
       x1(:,i1) = ctmp(:,map(i1))
    end do
  end subroutine sort_map_2d_c
  
  subroutine sort_map_2d_r(n1, nm, x1, map)
    integer, intent(in) :: n1, nm
    real(dp), dimension(n1, nm), intent(inout) :: x1
    integer, dimension(nm), intent(in) :: map

    real(dp), dimension(n1, nm) :: ctmp
    integer :: i1

    ctmp = x1
    do i1=1, nm
       x1(:,i1) = ctmp(:,map(i1))
    end do
  end subroutine sort_map_2d_r
  
  subroutine sort_map_3d_c(n1, n2, nm, x1, map)
    integer, intent(in) :: n1, n2, nm
    complex(dpc), dimension(n1,n2, nm), intent(inout) :: x1
    integer, dimension(nm), intent(in) :: map

    complex(dpc), dimension(n1, n2, nm) :: ctmp
    integer :: i1

    ctmp = x1
    do i1=1, nm
       x1(:,:,i1) = ctmp(:,:,map(i1))
    end do
  end subroutine sort_map_3d_c

  subroutine sort_map_3d_r(n1, n2, nm, x1, map)
    integer, intent(in) :: n1, n2, nm
    real(dp), dimension(n1,n2, nm), intent(inout) :: x1
    integer, dimension(nm), intent(in) :: map

    real(dp), dimension(n1, n2, nm) :: ctmp
    integer :: i1

    ctmp = x1
    do i1=1, nm
       x1(:,:,i1) = ctmp(:,:,map(i1))
    end do
  end subroutine sort_map_3d_r

!SU
  function sum_cf_2d(nc,nf, x, lxc0, dxf0) result(y)
    integer, intent(in)  :: nc, nf
    real(dp), dimension(0:nc,nf), intent(in)  :: x
    real(dp), intent(in), optional :: lxc0, dxf0
    real(dp) :: y
    
    real(dp) :: sinv(0:nc), cn(0:nc),lxc, dxf
    integer :: i1
    
    if (present(lxc0)) then
       lxc = lxc0
    else
       lxc = one
    end if
    
    if (present(dxf0)) then
       dxf = dxf0
    else
       dxf = twopi/nf
    end if
    sinv = (/ (sin(i1*pi/nc), i1=0,nc )/)
    
    y = (dxf*lxc*pi/nc)*sum(matmul(sinv,x))
    return
  end function sum_cf_2d

  function sum_fc_2d(nf,nc, x, dxf0, lxc0) result(y)
    integer, intent(in)  :: nc, nf
    real(dp), dimension(nf,0:nc), intent(in)  :: x
    real(dp), intent(in), optional :: lxc0, dxf0
    real(dp) :: y
    
    real(dp) :: sinv(0:nc), cn(0:nc),lxc, dxf
    integer :: i1
    
    if (present(lxc0)) then
       lxc = lxc0
    else
       lxc = one
    end if
    
    if (present(dxf0)) then
       dxf = dxf0
    else
       dxf = twopi/nf
    end if
    sinv = (/ (sin(i1*pi/nc), i1=0,nc )/)
    
    y = (dxf*lxc*pi/nc)*sum(matmul(x,sinv))
    return
  end function sum_fc_2d
  
  function sum_cb_1d(nn, x, lxc0) result(y)
    integer, intent(in)  :: nn
    real(dp), dimension(0:nn), intent(in)  :: x
    real(dp), intent(in), optional :: lxc0
    real(dp) :: y
    
    real(dp) :: sinv(0:nn), cn(0:nn),lxc
    integer :: i1
    
    if (present(lxc0)) then
       lxc = lxc0
    else
       lxc = one
    end if
    
    sinv = (/ (sin(i1*pi/nn), i1=0,nn )/)
    y = pi/nn*dot_product(x,sinv)*lxc
    
    return
  end function sum_cb_1d

!T
!TI
!U
!V
!W
!X
!Y
!Z
end module myutil
  

subroutine dverk (n, fcn, x, y, xend, tol, ind, c, nw, w)
  use mydef, only : sp, dp
  integer n, ind, nw, k
  real(dp) :: x,xend, tol, c(24), temp
  real(dp), dimension(n) :: y 
  real(dp), dimension(nw,9) :: w
  !
  !***********************************************************************
  !                                                                      *
  ! note added 11/14/85.                                                 *
  !                                                                      *
  ! if you discover any errors in this subroutine, please contact        *
  !                                                                      *
  !        kenneth r. jackson                                            *
  !        department of computer science                                *
  !        university of toronto                                         *
  !        toronto, ontario,                                             *
  !        canada   m5s 1a4                                              *
  !                                                                      *
  !        phone: 416-978-7075                                           *
  !                                                                      *
  !        electroni! mail:                                              *
  !        uucp:   {cornell,decvax,ihnp4,linus,uw-beaver}!utcsri!krj     *
  !        csnet:  krj@toronto                                           *
  !        arpa:   krj.toronto@csnet-relay                               *
  !        bitnet: krj%toronto@csnet-relay.arpa                          *
  !                                                                      *
  ! dverk is written in fortran 66.                                      *
  !                                                                      *
  ! the constants dwarf and rreb -- c(10) and c(11), respectively -- are *
  ! set for a  vax  in  double  precision.  they  should  be  reset,  as *
  ! described below, if this program is run on another machine.          *
  !                                                                      *
  ! the C array is declared in this subroutine to have one element only, *
  ! although  more  elements  are  referenced  in this subroutine.  this *
  ! causes some compilers to issue warning messages.  there is,  though, *
  ! no  error  provided  ! is declared sufficiently large in the calling *
  ! program, as described below.                                         *
  !                                                                      *
  ! the following external statement  for  fcn  was  added  to  avoid  a *
  ! warning  message  from  the  unix  f77 compiler.  the original dverk *
  ! comments and code follow it.                                         *
  !                                                                      *
  !***********************************************************************
  !
  external fcn
  !
  !***********************************************************************
  !                                                                      *
  !     purpose - this is a runge-kutta  subroutine  based  on  verner's *
  ! fifth and sixth order pair of formulas for finding approximations to *
  ! the solution of  a  system  of  first  order  ordinary  differential *
  ! equations  with  initial  conditions. it attempts to keep the global *
  ! error proportional to  a  tolerance  specified  by  the  user.  (the *
  ! proportionality  depends  on the kind of error control that is used, *
  ! as well as the differential equation and the range of integration.)  *
  !                                                                      *
  !     various options are available to the user,  including  different *
  ! kinds  of  error control, restrictions on step sizes, and interrupts *
  ! which permit the user to examine the state of the  calculation  (and *
  ! perhaps make modifications) during intermediate stages.              *
  !                                                                      *
  !     the program is efficient for non-stiff systems.  however, a good *
  ! variable-order-adams  method  will probably be more efficient if the *
  ! function evaluations are very costly.  such a method would  also  be *
  ! more suitable if one wanted to obtain a large number of intermediate *
  ! solution values by interpolation, as might be the case  for  example *
  ! with graphical output.                                               *
  !                                                                      *
  !                                    hull-enright-jackson   1/10/76    *
  !                                                                      *
  !***********************************************************************
  !                                                                      *
  !     use - the user must specify each of the following                *
  !                                                                      *
  !     n  number of equations                                           *
  !                                                                      *
  !   fcn  name of subroutine for evaluating functions - the  subroutine *
  !           itself must also be provided by the user - it should be of *
  !           the following form                                         *
  !              subroutine fcn(n, x, y, yprime)                         *
  !              integer n                                               *
  !              double precision x, y(n), yprime(n)                     *
  !                      *** et! ***                                     *
  !           and it should evaluate yprime, given n, x and y            *
  !                                                                      *
  !     x  independent variable - initial value supplied by user         *
  !                                                                      *
  !     y  dependent variable - initial values of components y(1), y(2), *
  !           ..., y(n) supplied by user                                 *
  !                                                                      *
  !  xend  value of x to which integration is to be carried out - it may *
  !           be less than the initial value of x                        *
  !                                                                      *
  !   tol  tolerance - the subroutine attempts to control a norm of  the *
  !           local  error  in  such  a  way  that  the  global error is *
  !           proportional to tol. in some problems there will be enough *
  !           damping  of  errors, as well as some cancellation, so that *
  !           the global error will be less than tol. alternatively, the *
  !           control   can   be  viewed  as  attempting  to  provide  a *
  !           calculated value of y at xend which is the exact  solution *
  !           to  the  problem y' = f(x,y) + e(x) where the norm of e(x) *
  !           is proportional to tol.  (the norm  is  a  max  norm  with *
  !           weights  that  depend on the error control strategy chosen *
  !           by the user.  the default weight for the k-th component is *
  !           1/max(1,abs(y(k))),  which therefore provides a mixture of *
  !           absolute and relative error control.)                      *
  !                                                                      *
  !   ind  indicator - on initial entry ind must be set equal to  either *
  !           1  or  2. if the user does not wish to use any options, he *
  !           should set ind to 1 - all that remains for the user to  do *
  !           then  is  to  declare c and w, and to specify nw. the user *
  !           may also  select  various  options  on  initial  entry  by *
  !           setting ind = 2 and initializing the first 9 components of *
  !           c as described in the next section.  he may also  re-enter *
  !           the  subroutine  with ind = 3 as mentioned again below. in *
  !           any event, the subroutine returns with ind equal to        *
  !              3 after a normal return                                 *
  !              4, 5, or 6 after an interrupt (see options c(8), c(9))  *
  !              -1, -2, or -3 after an error condition (see below)      *
  !                                                                      *
  !     c  communications vector - the dimension must be greater than or *
  !           equal to 24, unless option c(1) = 4 or 5 is used, in which *
  !           case the dimension must be greater than or equal to n+30   *
  !                                                                      *
  !    nw  first dimension of workspace w -  must  be  greater  than  or *
  !           equal to n                                                 *
  !                                                                      *
  !     w  workspace matrix - first dimension must be nw and second must *
  !           be greater than or equal to 9                              *
  !                                                                      *
  !     the subroutine  will  normally  return  with  ind  =  3,  having *
  ! replaced the initial values of x and y with, respectively, the value *
  ! of xend and an approximation to y at xend.  the  subroutine  can  be *
  ! called  repeatedly  with new values of xend without having to change *
  ! any other argument.  however, changes in tol, or any of the  options *
  ! described below, may also be made on such a re-entry if desired.     *
  !                                                                      *
  !     three error returns are also possible, in which  case  x  and  y *
  ! will be the most recently accepted values -                          *
  !     with ind = -3 the subroutine was unable  to  satisfy  the  error *
  !        requirement  with a particular step-size that is less than or *
  !        equal to hmin, which may mean that tol is too small           *
  !     with ind = -2 the value of hmin  is  greater  than  hmax,  which *
  !        probably  means  that the requested tol (which is used in the *
  !        calculation of hmin) is too small                             *
  !     with ind = -1 the allowed maximum number of fcn evaluations  has *
  !        been  exceeded,  but  this  can only occur if option c(7), as *
  !        described in the next section, has been used                  *
  !                                                                      *
  !     there are several circumstances that will cause the calculations *
  ! to  be  terminated,  along with output of information that will help *
  ! the user determine the cause of  the  trouble.  these  circumstances *
  ! involve  entry with illegal or inconsistent values of the arguments, *
  ! such as attempting a normal  re-entry  without  first  changing  the *
  ! value of xend, or attempting to re-enter with ind less than zero.    *
  !                                                                      *
  !***********************************************************************
  !                                                                      *
  !     options - if the subroutine is entered with ind = 1, the first 9 *
  ! components of the communications vector are initialized to zero, and *
  ! the subroutine uses only default values  for  each  option.  if  the *
  ! subroutine  is  entered  with ind = 2, the user must specify each of *
  ! these 9 components - normally he would first set them all  to  zero, *
  ! and  then  make  non-zero  those  that  correspond to the particular *
  ! options he wishes to select. in any event, options may be changed on *
  ! re-entry  to  the  subroutine  -  but if the user changes any of the *
  ! options, or tol, in the course of a calculation he should be careful *
  ! about  how  such changes affect the subroutine - it may be better to *
  ! restart with ind = 1 or 2. (components 10 to 24 of ! are used by the *
  ! program  -  the information is available to the user, but should not *
  ! normally be changed by him.)                                         *
  !                                                                      *
  !  c(1)  error control indicator - the norm of the local error is  the *
  !           max  norm  of  the  weighted  error  estimate  vector, the *
  !           weights being determined according to the value of c(1) -  *
  !              if c(1)=1 the weights are 1 (absolute error control)    *
  !              if c(1)=2 the weights are 1/abs(y(k))  (relative  error *
  !                 control)                                             *
  !              if c(1)=3 the  weights  are  1/max(abs(c(2)),abs(y(k))) *
  !                 (relative  error  control,  unless abs(y(k)) is less *
  !                 than the floor value, abs(c(2)) )                    *
  !              if c(1)=4 the weights are 1/max(abs(c(k+30)),abs(y(k))) *
  !                 (here individual floor values are used)              *
  !              if c(1)=5 the weights are 1/abs(c(k+30))                *
  !              for all other values of c(1), including  c(1) = 0,  the *
  !                 default  values  of  the  weights  are  taken  to be *
  !                 1/max(1,abs(y(k))), as mentioned earlier             *
  !           (in the two cases c(1) = 4 or 5 the user must declare  the *
  !           dimension of ! to be at least n+30 and must initialize the *
  !           components c(31), c(32), ..., c(n+30).)                    *
  !                                                                      *
  !  c(2)  floor value - used when the indicator c(1) has the value 3    *
  !                                                                      *
  !  c(3)  hmin specification - if not zero, the subroutine chooses hmin *
  !           to be abs(c(3)) - otherwise it uses the default value      *
  !              10*max(dwarf,rreb*max(weighted norm y/tol,abs(x))),     *
  !           where dwarf is a very small positive  machine  number  and *
  !           rreb is the relative roundoff error bound                  *
  !                                                                      *
  !  c(4)  hstart specification - if not zero, the subroutine  will  use *
  !           an  initial  hmag equal to abs(c(4)), except of course for *
  !           the restrictions imposed by hmin and hmax  -  otherwise it *
  !           uses the default value of hmax*(tol)**(1/6)                *
  !                                                                      *
  !  c(5)  scale specification - this is intended to be a measure of the *
  !           scale of the problem - larger values of scale tend to make *
  !           the method more reliable, first  by  possibly  restricting *
  !           hmax  (as  described  below) and second, by tightening the *
  !           acceptance requirement - if c(5) is zero, a default  value *
  !           of  1  is  used.  for  linear  homogeneous  problems  with *
  !           constant coefficients, an appropriate value for scale is a *
  !           norm  of  the  associated  matrix.  for other problems, an *
  !           approximation to  an  average  value  of  a  norm  of  the *
  !           jacobian along the trajectory may be appropriate           *
  !                                                                      *
  !  c(6)  hmax specification - four cases are possible                  *
  !           if c(6).ne.0 and c(5).ne.0, hmax is taken to be            *
  !              min(abs(c(6)),2/abs(c(5)))                              *
  !           if c(6).ne.0 and c(5).eq.0, hmax is taken to be  abs(c(6)) *
  !           if c(6).eq.0 and c(5).ne.0, hmax is taken to be            *
  !              2/abs(c(5))                                             *
  !           if c(6).eq.0 and c(5).eq.0, hmax is given a default  value *
  !              of 2                                                    *
  !                                                                      *
  !  c(7)  maximum number of function evaluations  -  if  not  zero,  an *
  !           error  return with ind = -1 will be caused when the number *
  !           of function evaluations exceeds abs(c(7))                  *
  !                                                                      *
  !  c(8)  interrupt number  1  -  if  not  zero,  the  subroutine  will *
  !           interrupt   the  calculations  after  it  has  chosen  its *
  !           preliminary value of hmag, and just before choosing htrial *
  !           and  xtrial  in  preparation for taking a step (htrial may *
  !           differ from hmag in sign, and may  require  adjustment  if *
  !           xend  is  near) - the subroutine returns with ind = 4, and *
  !           will resume calculation at the point  of  interruption  if *
  !           re-entered with ind = 4                                    *
  !                                                                      *
  !  c(9)  interrupt number  2  -  if  not  zero,  the  subroutine  will *
  !           interrupt   the  calculations  immediately  after  it  has *
  !           decided whether or not to accept the result  of  the  most *
  !           recent  trial step, with ind = 5 if it plans to accept, or *
  !           ind = 6 if it plans to reject -  y(*)  is  the  previously *
  !           accepted  result, while w(*,9) is the newly computed trial *
  !           value, and w(*,2) is the unweighted error estimate vector. *
  !           the  subroutine  will  resume calculations at the point of *
  !           interruption on re-entry with ind = 5 or 6. (the user  may *
  !           change ind in this case if he wishes, for example to force *
  !           acceptance of a step that would otherwise be rejected,  or *
  !           vice versa. he can also restart with ind = 1 or 2.)        *
  !                                                                      *
  !***********************************************************************
  !                                                                      *
  !  summary of the components of the communications vector              *
  !                                                                      *
  !     prescribed at the option       determined by the program         *
  !           of the user                                                *
  !                                                                      *
  !                                    c(10) rreb(rel roundoff err bnd)  *
  !     c(1) error control indicator   c(11) dwarf (very small mach no)  *
  !     c(2) floor value               c(12) weighted norm y             *
  !     c(3) hmin specification        c(13) hmin                        *
  !     c(4) hstart specification      c(14) hmag                        *
  !     c(5) scale specification       c(15) scale                       *
  !     c(6) hmax specification        c(16) hmax                        *
  !     c(7) max no of fcn evals       c(17) xtrial                      *
  !     c(8) interrupt no 1            c(18) htrial                      *
  !     c(9) interrupt no 2            c(19) est                         *
  !                                    c(20) previous xend               *
  !                                    c(21) flag for xend               *
  !                                    c(22) no of successful steps      *
  !                                    c(23) no of successive failures   *
  !                                    c(24) no of fcn evals             *
  !                                                                      *
  !  if c(1) = 4 or 5, c(31), c(32), ... c(n+30) are floor values        *
  !                                                                      *
  !***********************************************************************
  !                                                                      *
  !  an overview of the program                                          *
  !                                                                      *
  !     begin initialization, parameter checking, interrupt re-entries   *
  !  ......abort if ind out of range 1 to 6                              *
  !  .     cases - initial entry, normal re-entry, interrupt re-entries  *
  !  .     case 1 - initial entry (ind .eq. 1 or 2)                      *
  !  v........abort if n.gt.nw or tol.le.0                               *
  !  .        if initial entry without options (ind .eq. 1)              *
  !  .           set c(1) to c(9) equal to zero                          *
  !  .        else initial entry with options (ind .eq. 2)               *
  !  .           make c(1) to c(9) non-negative                          *
  !  .           make floor values non-negative if they are to be used   *
  !  .        end if                                                     *
  !  .        initialize rreb, dwarf, prev xend, flag, counts            *
  !  .     case 2 - normal re-entry (ind .eq. 3)                         *
  !  .........abort if xend reached, and either x changed or xend not    *
  !  .        re-initialize flag                                         *
  !  .     case 3 - re-entry following an interrupt (ind .eq. 4 to 6)    *
  !  v        transfer control to the appropriate re-entry point.......  *
  !  .     end cases                                                  .  *
  !  .  end initialization, etc.                                      .  *
  !  .                                                                v  *
  !  .  loop through the following 4 stages, once for each trial step .  *
  !  .     stage 1 - prepare                                          .  *
  !***********error return (with ind=-1) if no of fcn evals too great .  *
  !  .        calc slope (adding 1 to no of fcn evals) if ind .ne. 6  .  *
  !  .        calc hmin, scale, hmax                                  .  *
  !***********error return (with ind=-2) if hmin .gt. hmax            .  *
  !  .        calc preliminary hmag                                   .  *
  !***********interrupt no 1 (with ind=4) if requested.......re-entry.v  *
  !  .        calc hmag, xtrial and htrial                            .  *
  !  .     end stage 1                                                .  *
  !  v     stage 2 - calc ytrial (adding 7 to no of fcn evals)        .  *
  !  .     stage 3 - calc the error estimate                          .  *
  !  .     stage 4 - make decisions                                   .  *
  !  .        set ind=5 if step acceptable, else set ind=6            .  *
  !***********interrupt no 2 if requested....................re-entry.v  *
  !  .        if step accepted (ind .eq. 5)                              *
  !  .           update x, y from xtrial, ytrial                         *
  !  .           add 1 to no of successful steps                         *
  !  .           set no of successive failures to zero                   *
  !**************return(with ind=3, xend saved, flag set) if x .eq. xend *
  !  .        else step not accepted (ind .eq. 6)                        *
  !  .           add 1 to no of successive failures                      *
  !**************error return (with ind=-3) if hmag .le. hmin            *
  !  .        end if                                                     *
  !  .     end stage 4                                                   *
  !  .  end loop                                                         *
  !  .                                                                   *
  !  begin abort action                                                  *
  !     output appropriate  message  about  stopping  the  calculations, *
  !        along with values of ind, n, nw, tol, hmin,  hmax,  x,  xend, *
  !        previous xend,  no of  successful  steps,  no  of  successive *
  !        failures, no of fcn evals, and the components of y            *
  !     stop                                                             *
  !  end abort action                                                    *
  !                                                                      *
  !***********************************************************************
  !
  !     ******************************************************************
  !     * begin initialization, parameter checking, interrupt re-entries *
  !     ******************************************************************
  !
  !  ......abort if ind out of range 1 to 6

  if (ind.lt.1 .or. ind.gt.6) go to 500
  !
  !        cases - initial entry, normal re-entry, interrupt re-entries
  go to (5, 5, 45, 1111, 2222, 2222), ind
  !        case 1 - initial entry (ind .eq. 1 or 2)
  !  .........abort if n.gt.nw or tol.le.0

5 if (n.gt.nw .or. tol.le.0.d0) go to 500

  if (ind.eq. 1) then
  !              initial entry without options (ind .eq. 1)
  !              set c(1) to c(9) equal to 0

     do  k = 1, 9
        c(k) = 0.d0
     enddo
  else if (ind.eq. 2) then
     !              initial entry with options (ind .eq. 2)
     !              make c(1) to c(9) non-negative
     do k = 1, 9
        c(k) = dabs(c(k))
     enddo
     !              make floor values non-negative if they are to be used
     if (c(1).eq.4.d0 .or. c(1).eq.5.d0) then
        do k = 1, n
           c(k+30) = dabs(c(k+30))
        enddo
     endif
  endif
  
  !           initialize rreb, dwarf, prev xend, flag, counts
  c(10) = 2.d0**(-56)
  c(11) = 1.d-35
  !           set previous xend initially to initial value of x
  c(20) = x

  do k = 21, 24
     c(k) = 0.d0
  enddo

  go to 50
  !        case 2 - normal re-entry (ind .eq. 3)
  !  .........abort if xend reached, and either x changed or xend not
45 if (c(21).ne.0.d0 .and. &
        (x.ne.c(20) .or. xend.eq.c(20))) go to 500
  !           re-initialize flag
  c(21) = 0.d0
  go to 50
  !        case 3 - re-entry following an interrupt (ind .eq. 4 to 6)
  !           transfer control to the appropriate re-entry point..........
  !           this has already been handled by the computed go to        .
  !        end cases                                                     v
50 continue
  !
  !     end initialization, etc.
  !
  !     ******************************************************************
  !     * loop through the following 4 stages, once for each trial  step *
  !     * until the occurrence of one of the following                   *
  !     *    (a) the normal return (with ind .eq. 3) on reaching xend in *
  !     *        stage 4                                                 *
  !     *    (b) an error return (with ind .lt. 0) in stage 1 or stage 4 *
  !     *    (c) an interrupt return (with ind  .eq.  4,  5  or  6),  if *
  !     *        requested, in stage 1 or stage 4                        *
  !     ******************************************************************
  !
99999 continue
  !
  !        ***************************************************************
  !        * stage 1 - prepare - do calculations of  hmin,  hmax,  etc., *
  !        * and some parameter  checking,  and  end  up  with  suitable *
  !        * values of hmag, xtrial and htrial in preparation for taking *
  !        * an integration step.                                        *
  !        ***************************************************************
  !
  !***********error return (with ind=-1) if no of fcn evals too great
  if (c(7).eq.0.d0 .or. c(24).lt.c(7)) then
  else
     ind = -1
     return
  endif

  !
  !           calculate slope (adding 1 to no of fcn evals) if ind .ne. 6
  if (ind .eq. 6) then
  else
     call fcn(n, x, y, w(1,1))
     c(24) = c(24) + 1.d0
  endif

  !
  !           calculate hmin - use default unless value prescribed
  c(13) = c(3)
  if (c(3) .ne. 0.d0) go to 165
  !              calculate default value of hmin
  !              first calculate weighted norm y - c(12) - as specified
  !              by the error control indicator c(1)
  temp = 0.d0

  if (c(1) .ne. 1.d0) go to 115
  !                 absolute error control - weights are 1
  do k = 1, n
     temp = dmax1(temp, dabs(y(k)))
  enddo
  c(12) = temp
  go to 160

115 if (c(1) .ne. 2.d0) go to 120
  !                 relative error control - weights are 1/dabs(y(k)) so
  !                 weighted norm y is 1
  c(12) = 1.d0
  go to 160
120 if (c(1) .ne. 3.d0) go to 130
  !                 weights are 1/max(c(2),abs(y(k)))
  do k = 1, n
     temp = dmax1(temp, dabs(y(k))/c(2))
  enddo
  c(12) = dmin1(temp, 1.d0)
  go to 160
130 if (c(1) .ne. 4.d0) go to 140
  !                 weights are 1/max(c(k+30),abs(y(k)))
  do k = 1, n
     temp = dmax1(temp, dabs(y(k))/c(k+30))
  enddo
  c(12) = dmin1(temp, 1.d0)
  go to 160
140 if (c(1) .ne. 5.d0) go to 150
  !                 weights are 1/c(k+30)
  do k = 1, n
     temp = dmax1(temp, dabs(y(k))/c(k+30))
  enddo
  c(12) = temp
  go to 160
150 continue
  !                 default case - weights are 1/max(1,abs(y(k)))
  do k = 1, n
     temp = dmax1(temp, dabs(y(k)))
  enddo
  c(12) = dmin1(temp, 1.d0)
160 continue
  c(13) = 10.d0*dmax1(c(11),c(10)*dmax1(c(12)/tol,dabs(x)))
165 continue
  !
  !           calculate scale - use default unless value prescribed
  c(15) = c(5)
  if (c(5) .eq. 0.d0) c(15) = 1.d0
  !
  !           calculate hmax - consider 4 cases
  !           case 1 both hmax and scale prescribed
  if (c(6).ne.0.d0 .and. c(5).ne.0.d0) &
       c(16) = dmin1(c(6), 2.d0/c(5))
  !           case 2 - hmax prescribed, but scale not
  if (c(6).ne.0.d0 .and. c(5).eq.0.d0) c(16) = c(6)
  !           case 3 - hmax not prescribed, but scale is
  if (c(6).eq.0.d0 .and. c(5).ne.0.d0) c(16) = 2.d0/c(5)
  !           case 4 - neither hmax nor scale is provided
  if (c(6).eq.0.d0 .and. c(5).eq.0.d0) c(16) = 2.d0
  !
  !***********error return (with ind=-2) if hmin .gt. hmax
  if (c(13) .gt. c(16)) then 
     ind = -2
     return
  endif
  !
  !           calculate preliminary hmag - consider 3 cases
  if (ind .gt. 2) go to 175
  !           case 1 - initial entry - use prescribed value of hstart, if
  !              any, else default
  c(14) = c(4)
  if (c(4) .eq. 0.d0) c(14) = c(16)*tol**(1./6.)
  go to 185
175 if (c(23) .gt. 1.d0) go to 180
  !           case 2 - after a successful step, or at most  one  failure,
  !              use min(2, .9*(tol/est)**(1/6))*hmag, but avoid possible
  !              overflow. then avoid reduction by more than half.
  temp = 2.d0*c(14)
  if (tol .lt. (2.d0/.9d0)**6*c(19)) &
       temp = .9d0*(tol/c(19))**(1./6.)*c(14)
  c(14) = dmax1(temp, .5d0*c(14))
  go to 185
180 continue
  !           case 3 - after two or more successive failures
  c(14) = .5d0*c(14)
185 continue
  !
  !           check against hmax
  c(14) = dmin1(c(14), c(16))
  !
  !           check against hmin
  c(14) = dmax1(c(14), c(13))
  !
  !***********interrupt no 1 (with ind=4) if requested
  if (c(8) .eq. 0.d0) go to 1111
  ind = 4
  return
  !           resume here on re-entry with ind .eq. 4   ........re-entry..
1111 continue
  !
  !           calculate hmag, xtrial - depending on preliminary hmag, xend
  if (c(14) .ge. dabs(xend - x)) go to 190
  !              do not step more than half way to xend
  c(14) = dmin1(c(14), .5d0*dabs(xend - x))
  c(17) = x + dsign(c(14), xend - x)
  go to 195
190 continue
  !              hit xend exactly
  c(14) = dabs(xend - x)
  c(17) = xend
195 continue
  !
!           calculate htrial
            c(18) = c(17) - x
!
!        end stage 1
!
!        ***************************************************************
!        * stage 2 - calculate ytrial (adding 7 to no of  fcn  evals). *
!        * w(*,2), ... w(*,8)  hold  intermediate  results  needed  in *
!        * stage 3. w(*,9) is temporary storage until finally it holds *
!        * ytrial.                                                     *
!        ***************************************************************
!
            temp = c(18)/1398169080000.d0
!
            do k = 1, n
               w(k,9) = y(k) + temp*w(k,1)*233028180000.d0
            enddo
            call fcn(n, x + c(18)/6.d0, w(1,9), w(1,2))
!
            do k = 1, n
               w(k,9) = y(k) + temp*(   w(k,1)*74569017600.d0 &
                    + w(k,2)*298276070400.d0  )
            enddo
            call fcn(n, x + c(18)*(4.d0/15.d0), w(1,9), w(1,3))
!
            do k = 1, n
               w(k,9) = y(k) + temp*(   w(k,1)*1165140900000.d0 &
                    - w(k,2)*3728450880000.d0 + w(k,3)*3495422700000.d0 )
            enddo
            call fcn(n, x + c(18)*(2.d0/3.d0), w(1,9), w(1,4))
!
            do k = 1, n
               w(k,9) = y(k) + temp*( - w(k,1)*3604654659375.d0 &
                    + w(k,2)*12816549900000.d0 &
                    - w(k,3)*9284716546875.d0 &
                    + w(k,4)*1237962206250.d0 )
            enddo
            call fcn(n, x + c(18)*(5.d0/6.d0), w(1,9), w(1,5))
!
            do  k = 1, n
               w(k,9) = y(k) + temp*(   w(k,1)*3355605792000.d0 &
                    - w(k,2)*11185352640000.d0 &
                    + w(k,3)*9172628850000.d0 &
                    - w(k,4)*427218330000.d0 &
                    + w(k,5)*482505408000.d0  )
            enddo
            call fcn(n, x + c(18), w(1,9), w(1,6))
!
            do k = 1, n
               w(k,9) = y(k) + temp*( - w(k,1)*770204740536.d0 &
                    + w(k,2)*2311639545600.d0 &
                    - w(k,3)*1322092233000.d0 &
                    - w(k,4)*453006781920.d0 &
                    + w(k,5)*326875481856.d0  )
            enddo
            call fcn(n, x + c(18)/15.d0, w(1,9), w(1,7))
!
            do k = 1, n
               w(k,9) = y(k) + temp*(   w(k,1)*2845924389000.d0 &
                    - w(k,2)*9754668000000.d0 &
                    + w(k,3)*7897110375000.d0 &
                    - w(k,4)*192082660000.d0 &
                    + w(k,5)*400298976000.d0 &
                    + w(k,7)*201586000000.d0  )
            enddo
            call fcn(n, x + c(18), w(1,9), w(1,8))
            !
            !           calculate ytrial, the extrapolated approximation and store
            !              in w(*,9)
            do k = 1, n
               w(k,9) = y(k) + temp*(   w(k,1)*104862681000.d0 &
                    + w(k,3)*545186250000.d0 &
                    + w(k,4)*446637345000.d0 &
                    + w(k,5)*188806464000.d0 &
                    + w(k,7)*15076875000.d0 &
                    + w(k,8)*97599465000.d0   )
            enddo
            !
            !           add 7 to the no of fcn evals
            c(24) = c(24) + 7.d0
            !
            !        end stage 2
            !
            !        ***************************************************************
            !        * stage 3 - calculate the error estimate est. first calculate *
            !        * the  unweighted  absolute  error  estimate vector (per unit *
            !        * step) for the unextrapolated approximation and store it  in *
            !        * w(*,2).  then  calculate the weighted max norm of w(*,2) as *
            !        * specified by the error  control  indicator  c(1).  finally, *
            !        * modify  this result to produce est, the error estimate (per *
            !        * unit step) for the extrapolated approximation ytrial.       *
            !        ***************************************************************
            !
            !           calculate the unweighted absolute error estimate vector
            do k = 1, n
               w(k,2) = (   w(k,1)*8738556750.d0 &
                    + w(k,3)*9735468750.d0 &
                    - w(k,4)*9709507500.d0 &
                    + w(k,5)*8582112000.d0 &
                    + w(k,6)*95329710000.d0 &
                    - w(k,7)*15076875000.d0 &
                    - w(k,8)*97599465000.d0)/1398169080000.d0
            enddo
            !
            !           calculate the weighted max norm of w(*,2) as specified by
            !           the error control indicator c(1)
            temp = 0.d0
            if (c(1) .ne. 1.d0) go to 310
            !              absolute error control
            do k = 1, n
               temp = dmax1(temp,dabs(w(k,2)))
            enddo
               go to 360
310            if (c(1) .ne. 2.d0) go to 320
!              relative error control
               do k = 1, n
                  temp = dmax1(temp, dabs(w(k,2)/y(k)))
               enddo
               go to 360
320            if (c(1) .ne. 3.d0) go to 330
               !              weights are 1/max(c(2),abs(y(k)))
               do k = 1, n
                  temp = dmax1(temp, dabs(w(k,2))/ dmax1(c(2), dabs(y(k))) )
               enddo
               go to 360
330            if (c(1) .ne. 4.d0) go to 340
               !              weights are 1/max(c(k+30),abs(y(k)))
               do  k = 1, n
                  temp = dmax1(temp, dabs(w(k,2))/ dmax1(c(k+30), dabs(y(k))) )
               enddo
               go to 360
340            if (c(1) .ne. 5.d0) go to 350
               !              weights are 1/c(k+30)
               do k = 1, n
                  temp = dmax1(temp, dabs(w(k,2)/c(k+30)))
               enddo
               go to 360
350            continue
!              default case - weights are 1/max(1,abs(y(k)))
               do k = 1, n
                  temp = dmax1(temp, dabs(w(k,2)) / dmax1(1.d0, dabs(y(k))) )
               enddo
360            continue
!
!           calculate est - (the weighted max norm of w(*,2))*hmag*scale
!              - est is intended to be a measure of the error  per  unit
!              step in ytrial
            c(19) = temp*c(14)*c(15)
!
!        end stage 3
!
!        ***************************************************************
!        * stage 4 - make decisions.                                   *
!        ***************************************************************
!
!           set ind=5 if step acceptable, else set ind=6
            ind = 5
            if (c(19) .gt. tol) ind = 6
!
!***********interrupt no 2 if requested
            if (c(9) .eq. 0.d0) go to 2222
               return
!           resume here on re-entry with ind .eq. 5 or 6   ...re-entry..
 2222       continue
!
            if (ind .eq. 6) go to 410
!              step accepted (ind .eq. 5), so update x, y from xtrial,
!                 ytrial, add 1 to the no of successful steps, and set
!                 the no of successive failures to zero
               x = c(17)
               do k = 1, n
                  y(k) = w(k,9)
               enddo
               c(22) = c(22) + 1.d0
               c(23) = 0.d0
!**************return(with ind=3, xend saved, flag set) if x .eq. xend
               if (x .ne. xend) go to 405
                  ind = 3
                  c(20) = xend
                  c(21) = 1.d0
                  return
  405          continue
               go to 420
  410       continue
!              step not accepted (ind .eq. 6), so add 1 to the no of
!                 successive failures
               c(23) = c(23) + 1.d0
!**************error return (with ind=-3) if hmag .le. hmin
               if (c(14) .gt. c(13)) go to 415
                  ind = -3
                  return
  415          continue
  420       continue
!
!        end stage 4
!
      go to 99999
!     end loop
!
!  begin abort action
  500 continue
!
      write(6,505) ind, tol, x, n, c(13), xend, nw, c(16), c(20), &
           c(22), c(23), c(24), (y(k), k = 1, n)
  505 format( /// 1h0, 58hcomputation stopped in dverk with the following values -  &
           / 1h0, 5hind =, i4, 5x, 6htol  =, 1pd13.6, 5x, 11hx         =,  1pd22.15  &
           / 1h , 5hn   =, i4, 5x, 6hhmin =, 1pd13.6, 5x, 11hxend      =, 1pd22.15 &
           / 1h , 5hnw  =, i4, 5x, 6hhmax =, 1pd13.6, 5x, 11hprev xend =, 1pd22.15 &
           / 1h0, 14x, 27hno of successful steps    =, 0pf8.0 &
           / 1h , 14x, 27hno of successive failures =, 0pf8.0 &
           / 1h , 14x, 27hno of function evals      =, 0pf8.0 &
           / 1h0, 23hthe components of y are &
           // (1h , 1p5d24.15)                                           )
      !
      stop
!
!  end abort action
!
    end

