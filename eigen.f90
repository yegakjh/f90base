module mod_eigen
  use mydef
  use mymath
  implicit none

  private
  integer :: nm
  integer, parameter :: nb=64, nin=5, nout=6

  real(dp), parameter :: err08 = 1.0e-08_dp

  public :: eigen, eigen_qz, eigen_qz_vec, eigen_qz_vec_ch
  public :: eigen_zgeev
  interface eigen
     module procedure eigen_real, eigen_cmplx
  end interface
contains
  subroutine eigen_real(nn, amat,BMAT, egv, egvecs, out)
    integer, intent(in) :: nn
    real(dp), dimension(nn,nn), intent(in) :: amat, BMAT
    complex(dpc), dimension(nn), intent(out) :: egv
    real(dp), dimension(nn,nn), intent(out) :: egvecs
    integer, intent(in), optional :: out ! flag for display eigenvalues
    
    !     .. Local Scalars ..
    real(dp) :: small, norm
    INTEGER ::  I, INFO, J, LWKOPT, iout=0
    
    !     .. Local Arrays ..
    real(dp), dimension(nn) ::   alphar, alphai,beta
    real(dp), dimension(nn*nb) :: work
    real(dp) :: DUMMY(1,1)

    integer, dimension(nn)  :: ipiv
    complex(dpc), dimension(nn) :: tmpvec
    real(dpc), dimension(nn,nn)  :: tmpmat, tmpmat1
    
    !      local variable for dggev
    integer :: lda, ldb, ldvr, lwork
    
    !     .. External Subroutines ..
    EXTERNAL ::  DGGEV

    lda=nn;     ldb=nn;     LDVR=nn;     LWORK=nn*nb

    if (present(out)) then
       iout = out
       if (iout.eq.0) write(nout,*) 'Set to no-printing of eigenvalues at Subroutine eigen'
    else
       iout = 0
    end if

    ! since the matrix on exit is overwritten, copy to another one.
    tmpmat = amat;     tmpmat1 = bmat

    !        Solve the generalized eigenvalue problem
    CALL DGGEV('N','V',nn,tmpmat,LDA,tmpmat1,LDB,alphar,alphai,beta,DUMMY,1,EGVECS,LDVR,WORK,LWORK,INFO)

    IF (INFO.GT.0) THEN
       WRITE (NOUT,'(1X,A,I4)') 'Failure in DGGEV. INFO =', INFO
    ELSE
       SMALL = epsilon(small)
       DO J = 1, nn
          norm = sqrt(alphar(j)*alphar(j)+alphai(j)*alphai(j))
          IF (norm*SMALL.GE.ABS(BETA(J))) then
             WRITE (NOUT, *) 'Eigenvalue(', J, ')',' is numerically infinite or undetermined'
             write(nout, '(A,1x,e11.4,1x,e11.4,1x,A,e11.4)') 'ALPHA  = ', ALPHAR(J),ALPHAI(J), ', BETA = ', BETA(J)
             write(nout,'(A,1x,e11.4)') 'Small = ', small
             stop
          ELSE
             egv(j) = cmplx(alphar(j), alphai(j))/beta(j)
             if (iout.ne.0) then
                WRITE (NOUT,'(1X,A,I2,A,"(",E11.4,",",E11.4,")")') 'Eigenvalue(', J, ') = ', egv(j)
                WRITE (NOUT,'(1x,a,i2,a,/3(1x,"(",e11.4,")",:))') 'Eigenvector(', J, ')',(EGVECS(I,J),I=1,nn)

             endif
             norm = sqrt(sum(egvecs(:,j)*egvecs(:,j)))
             egvecs(:,j) = egvecs(:,j)/norm
             norm = sqrt(sum(egvecs(:,j)*egvecs(:,j)))
             if (abs(norm-one).gt. 1.01) then
                write(nout,*) norm
             endif
          END IF
       enddo
       
       LWKOPT = WORK(1)
       IF (LWORK.LT.LWKOPT) THEN
          WRITE (NOUT,*)
          WRITE (NOUT,'(1x,a,i5,/1x,a,i5)') 'Optimum workspace required = ', LWKOPT, 'Workspace provided         = ', LWORK
       END IF
    END IF

    return
  end subroutine eigen_real

subroutine eigen_cmplx(nn, amat,BMAT, egv, egvecs, out)
    implicit none
    integer, intent(in) :: nn
    complex(dpc), dimension(nn,nn), intent(in) :: amat, BMAT
    complex(dpc), dimension(nn), intent(out) :: egv
    complex(dpc), dimension(nn,nn), intent(out) :: egvecs
    integer, intent(in), optional :: out ! flag for display eigenvalues
    
    !     .. Local Scalars ..
    real(dp) :: small, norm
    INTEGER ::  I, INFO, J, LWKOPT, iout=0
    
    !     .. Local Arrays ..
    complex(dpc), dimension(nn) ::   alpha,beta
    complex(dpc), dimension(nn*nb) :: work
    real(dp), dimension(nn*8) :: rwork
    complex(dpc) :: DUMMY(1,1)

    integer, dimension(nn)  :: ipiv
    complex(dpc), dimension(nn) :: tmpvec
    complex(dpc), dimension(nn,nn)  :: tmpmat, tmpmat1
    
    !      local variable for dggev
    integer :: lda, ldb, ldvr, lwork
    
    !     .. External Subroutines ..
    EXTERNAL ::  ZGGEV

    
    lda=nn
    ldb=nn
    LDVR=nn
    LWORK=nn*nb

    if (present(out)) then
       iout = out
       if (iout.eq.0) write(*,*) 'Set to no-printing of eigenvalues at Subroutine eigen'
    else
       iout = 0
    end if

    ! since the matrix on exit is overwritten, copy to another one.
    tmpmat = amat
    tmpmat1 = bmat

    !        Solve the generalized eigenvalue problem
    CALL ZGGEV('N','V',nn,tmpmat,LDA,tmpmat1,LDB,alpha,beta,DUMMY,1,EGVECS,LDVR,WORK,LWORK,rwork, INFO)

    IF (INFO.GT.0) THEN
       WRITE (NOUT,'(1X,A,I4)') 'Failure in ZGGEV. INFO =', INFO
    ELSE
       SMALL = epsilon(small)
       DO J = 1, nn
          norm = abs(alpha(j))
          IF (norm*SMALL.GE.ABS(BETA(J))) then
             WRITE (NOUT, *) 'Eigenvalue(', J, ')',' is numerically infinite or undetermined'
             write(nout, '(A,1x,2(e11.4,1x),A,2(e11.4,1x))') 'ALPHA  = ', alpha(j), ', BETA = ', BETA(J)
             write(nout,'(A,1x,e11.4)') 'Small = ', small
             stop
          ELSE
             egv(j) = alpha(j)/beta(j)
             if (iout.ne.0) then
                WRITE (NOUT,'(1X,A,I2,A,"(",E11.4,",",E11.4,")")') 'Eigenvalue(', J, ') = ', egv(j)
                WRITE (NOUT,'(1x,a,i2,a,/3(1x,"(",e11.4,")",:))') 'Eigenvector(', J, ')',(EGVECS(I,J),I=1,nn)

             endif
             norm = sqrt(sum(egvecs(:,j)*conjg(egvecs(:,j))))
             egvecs(:,j) = egvecs(:,j)/norm
          END IF
       enddo
       
       LWKOPT = WORK(1)
       IF (LWORK.LT.LWKOPT) THEN
          WRITE (NOUT,*)
          WRITE (NOUT,'(1x,a,i5,/1x,a,i5)') 'Optimum workspace required = ', LWKOPT, 'Workspace provided         = ', LWORK
       END IF
    END IF

    return
  end subroutine eigen_cmplx

  subroutine eigen_qz(n, a, b, n0, egv)
    implicit none
    !     .. Parameters ..
    INTEGER, parameter :: NOUT=6
    INTEGER, parameter :: LDQ=1, LDZ=1
    
  !     .. Input ..
    integer :: n,n0
    complex(dpc), dimension(n,n) :: a,b
    complex(dpc), dimension(n0) :: egv
    
    !     .. Local Scalars ..
    COMPLEX(dpc) ::       E
    INTEGER      ::    I, IFAIL, IHI, ILO, INFO, IROWS, J, JWORK, lwork, lda, ldb
    integer      ::    icount=0
    CHARACTER    ::    COMPQ, COMPZ, JOB
    real(dp) :: sfmin
    
    !     .. Local Arrays ..
    COMPLEX(dpc) ::      ALPHA(n), BETA(n), Q(LDQ,LDQ),WORK(6*n),Z(LDZ,LDZ), tau(n)
    real(dp) :: LSCALE(n), RSCALE(n), RWORK(6*n)
    CHARACTER ::       CLABS(1), RLABS(1)
    
    !     .. External Subroutines ..
    EXTERNAL ::   ZGEQRF, ZGGBAL, ZGGHRD, ZHGEQZ, ZUNMQR
    !CHECK    EXTERNAL :: X0$DBF
    !     .. Intrinsic Functions ..
    INTRINSIC        DBLE, DIMAG, NINT
    
    lwork = 6*n; lda = n; ldb = n
    sfmin = tiny(sfmin)
    icount = 0
    !        Balance matrix pair (A,B)
    JOB = 'B'
    CALL ZGGBAL(JOB,N,A,LDA,B,LDB,ILO,IHI,LSCALE,RSCALE,RWORK,INFO)

    !        Reduce B to triangular form using QR
    IROWS = IHI + 1 - ILO
    CALL ZGEQRF(IROWS,IROWS,B(ILO,ILO),LDB,TAU,WORK,LWORK,INFO)
    
    !        Apply the orthogonal transformation to A
    CALL ZUNMQR('L','C',IROWS,IROWS,IROWS,B(ILO,ILO),LDB,TAU, &
         A(ILO,ILO),LDA,WORK,LWORK,INFO)
    
    !        Compute the generalized Hessenberg form of (A,B)
    COMPQ = 'N';     COMPZ = 'N'
    CALL ZGGHRD(COMPQ,COMPZ,IROWS,1,IROWS,A(ILO,ILO),LDA,B(ILO,ILO)&
         ,LDB,Q,LDQ,Z,LDZ,INFO)
    
    !        Routine ZHGEQZ
    !        Workspace query: JWORK = -1
    JWORK = -1;     JOB = 'E'
    CALL ZHGEQZ(JOB,COMPQ,COMPZ,N,ILO,IHI,A,LDA,B,LDB,ALPHA,BETA,Q,&
         LDQ,Z,LDZ,WORK,JWORK,RWORK,INFO)
  
  !        Compute the generalized Schur form
  !        if the workspace LWORK is adequate
  IF (NINT(DBLE(WORK(1))).LE.LWORK) THEN 
     CALL ZHGEQZ(JOB,COMPQ,COMPZ,N,ILO,IHI,A,LDA,B,LDB,ALPHA, &
          BETA,Q,LDQ,Z,LDZ,WORK,LWORK,RWORK,INFO)
     
     !        Print the generalized eigenvalues
     !        Note: the actual values of beta are real and non-negative
     DO  I = 1, N
        IF (abs(BETA(I)) > sfmin) THEN
           icount = icount+1
           if (icount > n0) write(*,*) 'ICOUNT > N0', icount, n0
           egv(icount) = ALPHA(I)/BETA(I)
        END IF
     end DO
  ELSE
     WRITE (NOUT,'(1X,"Insufficient workspace for array WORK in F08XSF, ZHGEQZ")')
  END IF
  n0 = icount
  return
end subroutine eigen_qz

subroutine eigen_qz_vec(nn, a, b, n0, egv, iflag, ip0)
  ! Ax = lambda Bx
  implicit none
  INTEGER :: nout=6, nin = 5
  
  !  .. Argument ..
  integer :: nn, n0
  complex(dpc), dimension(nn,nn), intent(inout) :: a,b
  complex(dpc), dimension(n0), intent(out) :: egv
  integer, optional :: iflag, ip0
  logical :: iprint
  ! ip0 defined -> iprint = .True. , not then .False.
  ! iflag = 1, lvec only, 2 = rvec only, 0 = both and default
  
  !   .. Local Scalars ..
  INTEGER :: I, ICOLS, IFAIL, IHI, ILO, INFO, IROWS, J, JWORK, M, LWORK, ICOUNT
  LOGICAL :: ILEFT, IRIGHT
  CHARACTER :: COMPQ, COMPZ, HOWMNY, JOB, SIDE
  real(dp) :: err, sfmin

  !   .. Local Arrays ..
  complex(dpc), dimension(nn) ::  ALPHA, BETA,TAU
  complex(dpc), dimension(nn,nn) :: Q, Z

  complex(dpc), dimension(:), allocatable :: WORK
  real(dp), dimension(:), allocatable :: RWORK

  real(dp),dimension(nn) :: LSCALE, RSCALE
  
  LOGICAL  :: SELECT(NN)
  CHARACTER :: CLABS(1), RLABS(1)
  
  !   .. External Subroutines ..
  EXTERNAL :: ZGEQRF, ZGGBAK, ZGGBAL,ZGGHRD, ZHGEQZ, ZTGEVC, ZUNGQR, ZUNMQR
  
  !   ILEFT  is TRUE if left  eigenvectors are required
  !   IRIGHT is TRUE if right eigenvectors are required

  if (present(iflag)) then
     select case(iflag)
     case(0)
        ileft = .TRUE. ; iright = .TRUE.
     case(1)
        ileft = .TRUE. ; iright = .FALSE.
     case(2)
        ileft = .FALSE. ; iright = .TRUE.
     end select
  else
     ileft = .TRUE. ; iright = .TRUE.
  end if
  
  if (present(ip0)) then
     select case(ip0)
     case(:-1)
        iprint = .false.
        nout = -ip0
     case(0)
        iprint = .false.
        nout = 6
     case(1:)
        iprint = .TRUE.
        nout = ip0
     end select
  else
     iprint = .FALSE.
     nout = 6
  end if

  LWORK = 6*nn; 
  allocate(work(lwork)); allocate(rwork(lwork))
  icount = 0
  sfmin = tiny(sfmin)

  ! Balance matrix pair (A,B)
  JOB = 'B'
  CALL ZGGBAL(JOB,NN,A,NN,B,NN,ILO,IHI,LSCALE,RSCALE,RWORK,INFO)
  
  ! Reduce B to triangular form using QR
  IROWS = IHI + 1 - ILO
  ICOLS = NN + 1 - ILO
  CALL ZGEQRF(IROWS,ICOLS,B(ILO,ILO),NN,TAU,WORK,LWORK,INFO)
  
  ! Apply the orthogonal transformation to A
  CALL ZUNMQR('L','C',IROWS,ICOLS,IROWS,B(ILO,ILO),NN,TAU, &
       A(ILO,ILO),NN,WORK,LWORK,INFO)
  
  ! Initialize Q (for left eigenvectors)
  IF (ILEFT) THEN
     forall (i=1:nn, j=1:nn, i/=j) q(i,j) = zero
     forall (i=1:nn) q(i,i) = one

     !CALL F06TFF('Lower',IROWS-1,IROWS-1,B(ILO+1,ILO),NN, Q(ILO+1,ILO),NN)
     forall(i=ILO+1:ILO+IROWS-1, j=ILO: ILO+IROWS-2, i>j)
        Q(i,j) = B(i,j)
     end forall
     !      do i=ILO+1, ILO+irows-1
     !         do j=ILO, ILO+irows-2
     !            if ( i > j) then
     !               Q(i,j) = B(i,j)
     !            end if
     !         end do
     !      end do

     CALL ZUNGQR(IROWS,IROWS,IROWS,Q(ILO,ILO),NN,TAU,WORK,LWORK, INFO)
  END IF
  
  ! Initialize Z for right eigenvectors
  IF (IRIGHT) then
     forall (i=1:nn, j=1:nn, i/=j) z(i,j) = zero
     forall (i=1:nn) z(i,i) = one
  end IF
  
  ! Compute the generalized Hessenberg form of (A,B)
  COMPQ = 'V';   COMPZ = 'V'
  CALL ZGGHRD(compq,compz,NN,ILO,IHI,A,NN,B,NN,Q,NN,Z,NN,INFO)
  
  ! Routine ZHGEQZ
  ! Workspace query: JWORK = -1
  JWORK = -1;   JOB = 'S'
  CALL ZHGEQZ(JOB,COMPQ,COMPZ,NN,ILO,IHI,A,NN,B,NN,ALPHA,BETA,Q,NN,Z,NN,WORK,JWORK,RWORK,INFO)
  
  if (iprint) then
     WRITE (NOUT,*)
     WRITE (NOUT,'(1X,"Minimal required LWORK = ",I6)') NINT(DBLE(WORK(1)))
     WRITE (NOUT,'(1X,"Actual value of  LWORK = ",I6)') LWORK
  end if
  
  !      Compute the generalized Schur form
  !      if the workspace LWORK is adequate
  
  IF (NINT(DBLE(WORK(1))).LE.LWORK) THEN
     CALL ZHGEQZ(JOB,COMPQ,COMPZ,NN,ILO,IHI,A,NN,B,NN,ALPHA,BETA,Q,NN,Z,NN,WORK,LWORK,RWORK,INFO)
     
     ! Print the generalized eigenvalues
     ! Note: the actual values of beta are real and non-negative
     if (iprint) WRITE (NOUT,'(1X,"Generalized eigenvalues")')
     DO I = 1, NN
        IF (abs(BETA(I)) > sfmin) THEN
           icount = icount+1
           if (icount > n0) write(*,*) 'ICOUNT > N0', icount, n0
           egv(icount) = ALPHA(I)/BETA(I)
           select(i) = .TRUE.
           
           if (iprint) WRITE (NOUT,'(a,2(1x,i3),2(1x,e10.3))') &
                'EGV = ', I,icount,  egv(icount)
        ELSE
           select(i) = .FALSE.
        END IF
     end DO
     
     ! this is the right number of well-calculated eigenvalues
     n0 = icount
     
     ! Compute left and right generalized eigenvectors of the balanced matrix
     HOWMNY = 'B'
     IF (ILEFT .AND. IRIGHT) THEN
        SIDE = 'B'
     ELSE IF (ILEFT) THEN
        SIDE = 'L'
     ELSE IF (IRIGHT) THEN
        SIDE = 'R'
     END IF
     
     CALL ZTGEVC(SIDE,HOWMNY,SELECT,NN,A,NN,B,NN,Q,NN,Z,NN,NN,M,WORK,RWORK,INFO)
     !      Compute right eigenvectors of the original matrix
     IF (IRIGHT) THEN
        JOB = 'B';         SIDE = 'R'
        CALL ZGGBAK(JOB,SIDE,NN,ILO,IHI,LSCALE,RSCALE,NN,Z,NN,INFO)
        
        icount = 0
        do i=1, nn
           if (select(i)) then
              icount = icount+1
              b(:,icount) = z(:,i)
           end if
        enddo
     END IF
     
     !      Compute left eigenvectors of the original matrix
     IF (ILEFT) THEN
        JOB = 'B';         SIDE = 'L'
        CALL ZGGBAK(JOB,SIDE,NN,ILO,IHI,LSCALE,RSCALE,NN,Q,NN,INFO)
        
        icount = 0
        do i=1, nn
           if (select(i)) then
              icount = icount+1
              a(:,icount) = q(:,i)
           end if
        enddo
     END IF
  ELSE
     WRITE (NOUT,'(1X,"Insufficient workspace for array WORK ZHGEQZ")')
  END IF
  
  deallocate(work, rwork)
  return
end subroutine eigen_qz_vec

! This verision check the left and right eigenvectors 
subroutine eigen_qz_vec_ch(nn, a, b, n0, egv, iflag, ip0)
  implicit none
  INTEGER, parameter :: nout=6, nin = 5
  
  !  .. Argument ..
  integer :: nn, n0
  complex(dpc), dimension(nn,nn), intent(inout) :: a,b
  complex(dpc), dimension(n0), intent(out) :: egv
  integer, optional :: iflag, ip0
  logical :: iprint
  ! ip0 defined -> iprint = .True. , not then .False.
  ! iflag = 1, lvec only, 2 = rvec only, 0 = both and default
  
  !   .. Local Scalars ..
  INTEGER :: I, ICOLS, IFAIL, IHI, ILO, INFO, IROWS, J, JWORK, M, LWORK, ICOUNT
  LOGICAL :: ILEFT, IRIGHT
  CHARACTER :: COMPQ, COMPZ, HOWMNY, JOB, SIDE
  real(dp) :: err, sfmin

  !   .. Local Arrays ..
  complex(dpc), dimension(nn) ::  ALPHA, BETA,TAU,vec
  complex(dpc), dimension(nn,nn) :: Q, Z, a0, b0

  complex(dpc), dimension(:), allocatable :: WORK
  real(dp), dimension(:), allocatable :: RWORK

  real(dp),dimension(nn) :: LSCALE, RSCALE

  LOGICAL  :: SELECT(NN)
  CHARACTER :: CLABS(1), RLABS(1)
  
  !   .. External Subroutines ..
  EXTERNAL :: ZGEQRF, ZGGBAK, ZGGBAL,ZGGHRD, ZHGEQZ, ZTGEVC, ZUNGQR, ZUNMQR
  
  !   ILEFT  is TRUE if left  eigenvectors are required
  !   IRIGHT is TRUE if right eigenvectors are required

  if (present(iflag)) then
     select case(iflag)
     case(0)
        ileft = .TRUE. ; iright = .TRUE.
     case(1)
        ileft = .TRUE. ; iright = .FALSE.
     case(2)
        ileft = .FALSE. ; iright = .TRUE.
     end select
  else
     ileft = .TRUE. ; iright = .TRUE.
  end if
  
  if (present(ip0)) then
     iprint = .TRUE.
  else
     iprint = .FALSE.
  end if

  a0 = a;   b0 = b

  LWORK = 6*nn; 
  allocate(work(lwork)); allocate(rwork(lwork))
  icount = 0
  sfmin = tiny(sfmin)

  !      Balance matrix pair (A,B)
  JOB = 'B'
  CALL ZGGBAL(JOB,NN,A,NN,B,NN,ILO,IHI,LSCALE,RSCALE,RWORK,INFO)
  
  !      Reduce B to triangular form using QR
  IROWS = IHI + 1 - ILO
  ICOLS = NN + 1 - ILO
  CALL ZGEQRF(IROWS,ICOLS,B(ILO,ILO),NN,TAU,WORK,LWORK,INFO)
  
  !      Apply the orthogonal transformation to A
  CALL ZUNMQR('L','C',IROWS,ICOLS,IROWS,B(ILO,ILO),NN,TAU, &
       A(ILO,ILO),NN,WORK,LWORK,INFO)
  
  !      Initialize Q (for left eigenvectors)
  IF (ILEFT) THEN
     forall (i=1:nn, j=1:nn, i/=j) q(i,j) = zero
     forall (i=1:nn) q(i,i) = one

     forall(i=ILO+1:ILO+IROWS-1, j=ILO: ILO+IROWS-2, i>j)
        Q(i,j) = B(i,j)
     end forall

     CALL ZUNGQR(IROWS,IROWS,IROWS,Q(ILO,ILO),NN,TAU,WORK,LWORK, INFO)
  END IF
  
  !      Initialize Z for right eigenvectors
  IF (IRIGHT) then
     forall (i=1:nn, j=1:nn, i/=j) z(i,j) = zero
     forall (i=1:nn) z(i,i) = one
  end IF
  
  !      Compute the generalized Hessenberg form of (A,B)
  COMPQ = 'V';   COMPZ = 'V'
  CALL ZGGHRD(compq,compz,NN,ILO,IHI,A,NN,B,NN,Q,NN,Z,NN,INFO)
  
  !      Routine ZHGEQZ
  !      Workspace query: JWORK = -1

  JWORK = -1;   JOB = 'S'
  CALL ZHGEQZ(JOB,COMPQ,COMPZ,NN,ILO,IHI,A,NN,B,NN,ALPHA,BETA,Q,NN,Z,NN,WORK,JWORK,RWORK,INFO)

  if (iprint) then
     WRITE (NOUT,*)
     WRITE (NOUT,'(1X,"Minimal required LWORK = ",I6)') NINT(DBLE(WORK(1)))
     WRITE (NOUT,'(1X,"Actual value of  LWORK = ",I6)') LWORK
  end if
     
  !      Compute the generalized Schur form
  !      if the workspace LWORK is adequate
  
  IF (NINT(DBLE(WORK(1))).LE.LWORK) THEN
     CALL ZHGEQZ(JOB,COMPQ,COMPZ,NN,ILO,IHI,A,NN,B,NN,ALPHA,BETA,Q,NN,Z,NN,WORK,LWORK,RWORK,INFO)
     
     !      Print the generalized eigenvalues
     !      Note: the actual values of beta are real and non-negative
     if (iprint) WRITE (NOUT,'(1X,"Generalized eigenvalues")')
     DO I = 1, NN
        IF (abs(BETA(I)) > sfmin) THEN
           icount = icount+1
           if (icount > n0) write(*,*) 'ICOUNT > N0', icount, n0
           egv(icount) = ALPHA(I)/BETA(I)
           select(i) = .TRUE.
           
           if (iprint) WRITE (NOUT,'(a,2(1x,i3),2(1x,e10.3))') &
                'EGV = ', I,icount,  egv(icount)
        ELSE
           select(i) = .FALSE.
        END IF
     end DO
     
     ! this is the right number of well-calculated eigenvalues
     n0 = icount
     
     !      Compute left and right generalized eigenvectors
     !      of the balanced matrix
     
     HOWMNY = 'B'
     IF (ILEFT .AND. IRIGHT) THEN
        SIDE = 'B'
     ELSE IF (ILEFT) THEN
        SIDE = 'L'
     ELSE IF (IRIGHT) THEN
        SIDE = 'R'
     END IF

     CALL ZTGEVC(SIDE,HOWMNY,SELECT,NN,A,NN,B,NN,Q,NN,Z,NN,NN,M,WORK,RWORK,INFO)
     
     !      Compute right eigenvectors of the original matrix
     IF (IRIGHT) THEN
        JOB = 'B';         SIDE = 'R'
        CALL ZGGBAK(JOB,SIDE,NN,ILO,IHI,LSCALE,RSCALE,NN,Z,NN,INFO)
        
        icount = 0
        if (iprint) then
           write(nout,*)
           write(nout,*) 'Right Eigenvector Check (eigen_qz_vec)'
        end if
        
        do i=1, nn
           if (select(i)) then
              icount = icount+1
              b(:,icount) = z(:,i)
              vec = z(:,i)
              err = vec_mod(nn, matmul(A0,vec)-egv(icount)*matmul(B0,vec))/vec_mod(nn,vec)
              if (err > err08 ) then
                 write(*,'(a,e10.3,a,e10.3,a,i3)') &
                      'Right EigenVec Err =  ',err, '>' ,err08,' at ', i
              end if
           end if
        enddo
     END IF
     
     !      Compute left eigenvectors of the original matrix
     
     IF (ILEFT) THEN
        JOB = 'B';         SIDE = 'L'
        CALL ZGGBAK(JOB,SIDE,NN,ILO,IHI,LSCALE,RSCALE,NN,Q,NN,INFO)
        
        icount = 0
        if (iprint) then
           write(nout,*) ''
           write(nout,*) 'Left Eigenvector Check (eigen_qz_vec)'
        end if
        do i=1, nn
           if (select(i)) then
              icount = icount+1
              a(:,icount) = q(:,i)
              vec = conjg(Q(:,i))
              err = vec_mod(nn,matmul(vec,A0)-egv(icount)*matmul(vec, B0))/vec_mod(nn,vec)
              if (err > err08 ) then
                 write(*,'(a,e10.3,a,e10.3,a,i3)') 'Left EigenVec Err =  ',err, '>' ,err08,' at ', i
              end if
           end if
        enddo
     END IF
     
     ! Now Check Orthogonality of left and right eigenvectors
     if (ileft .and. iright) then
        q(1:n0,:) = matmul(transpose(conjg(a(:,1:n0))),B0(:,:))
        q(1:n0,1:n0) = matmul(Q(1:n0,:), b(:,1:n0))
        
        if (iprint) then
           write(nout,*) ''
           write(nout,*) 'Left-Right Eigenvector Orthogonality Check (eigen_qz_vec)'
        end if

        do i=1,n0
           if (abs(q(i,i)) > sfmin) then
              do j=1,n0
                 err = abs(q(i,j)/q(i,i))
                 if ( err > err08 .and. i /= j) &
                      write(*,'(a,2(1x,i3),1x,e10.3)') 'Non-orthogonal:',i,j,err
              end do
           else
              write(nout,*) ' To small q(i,i) ',i
           end if
        end do
     end if
  ELSE
     WRITE (NOUT,'(1X,"Insufficient workspace for array WORK ZHGEQZ")')
  END IF
  
  deallocate(work, rwork)
  return
end subroutine eigen_qz_vec_ch

  subroutine eigen_zgeev(jobvl, jobvr, na, aa, egv, vl, vr, info)
    use mydef
    use mymath
    implicit none
    character :: jobvl, jobvr
    integer :: na, info
    complex(dpc), dimension(na,na) :: aa, vl, vr

    complex(dpc), dimension(na) :: egv

    integer :: LWORK, LWKOPT, i1
    complex(dpc) :: WORK((1+nb)*na)
    real(dp) ::  RWORK(2*na), tmpr
    complex(dpc) :: tmpc

    EXTERNAL ZGEEV

    LWORK=(1+NB)*na
    
    CALL ZGEEV(jobvl,jobvr,na,aa,na, egv,vl,na,vr,na,WORK,LWORK,RWORK,INFO)
    LWKOPT = WORK(1)

    IF (INFO.EQ.0) THEN
       if (jobvr == 'V') then
          do i1=1, na
             tmpr = dot_product(vr(:,i1),vr(:,i1))
             vr(:,i1) = vr(:,i1)/sqrt(tmpr)
          enddo
          if (jobvl == 'V') then
             do i1=1, na
                tmpc = dot_product(vl(:,i1),vr(:,i1))
                vl(:,i1) = vl(:,i1)/conjg(tmpc)
             end do
          end if
       end if
    ELSE
       WRITE (NOUT,*)
       WRITE (NOUT,99998) 'Failure in ZGEEV.  INFO = ', INFO
    END IF

    !Print workspace information
    
    IF (LWORK.LT.LWKOPT) THEN
       WRITE (NOUT,*)
       WRITE (NOUT,99997) 'Optimum workspace required = ', LWKOPT,&
            'Workspace provided         = ', LWORK
    END IF
    
99998 FORMAT (1X,A,I4)
99997 FORMAT (1X,A,I5,/1X,A,I5)
    
    return
  end subroutine eigen_zgeev

end module mod_eigen
