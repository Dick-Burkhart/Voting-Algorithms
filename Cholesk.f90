
Module Cholesk

!  Standard routines for solving matrix equations
 
   Use Util
   Use Precisn
   Use Output
   Implicit none
   
   Interface LU_solve
     Module procedure LU_solveR, LU_solveD
   End Interface LU_solve

   Interface LU_decomp
     Module procedure LU_decompR, LU_decompD
   End Interface LU_decomp

   Interface Mat_solve
     Module procedure Mat_solveR1, Mat_solveR2, Mat_solveD1, Mat_solveD2
   End Interface Mat_solve
   
   Interface RR_solve
     Module procedure RR_solveR, RR_solveD
   End Interface RR_solve
   
!  For Line_search subroutine:
   
   Type Srch_parm
     Integer    :: ipc       ! Print code
     Integer    :: dmp_opt   ! Damping option for the function evaluation 
                             !   (1 : use residual, or 2 : use functional)
     Real(Dblp) :: step      ! Initial step size for 't'.
     Real(Dblp) :: t_low     ! Lower bound for 't'.
     Real(Dblp) :: t0        ! Initial value for 't'.
     Real(Dblp) :: t_high    ! Upper bound for 't'.
     Real(Dblp) :: tol       ! Absolute convergence tolerance for steps in 't'.
     Real(Dblp) :: change    ! Fractional convergence tolerance for change in f(t).
   End Type Srch_parm
     
   Real(Dblp), Parameter :: zero= 0.0d0, one= 1.0d0, M_one= -1.0d0
   Logical, Parameter :: Test_LU= .false.

Contains
 
   Subroutine LSQ_solve (n,m, A,y, x)
   
!    Compute the least squares solution to A.x = y
!    given n >= m, using the Cholesky factorization
!    of A^t.A.

     Integer, Intent(in) :: n       ! # rows = size of y
     Integer, Intent(in) :: m       ! # columns = size of x
     Real,    Intent(in) :: A(:,:)  !(n,m) matrix
     Real,    Intent(in) :: y(:)    !(n) data vector
     Real,   Intent(out) :: x(:)    !(m) solution vector
!  Local:
     Real(Dblp), Parameter :: eps= 1.0E-6
     Real(Dblp) :: xd(m), Ad(m,n), B(m,m+1)
     Integer       :: err_col
     
     Ad= Transpose(A);  B(:,m+1)= zero
     B(:,:m)= Matmul(Ad,A);  xd= Matmul(Ad,y)
     
     Call Cholesky (.true., .true., eps, B, err_col)
     Call Substitution (.true., .true., .true., B,xd)
     
     x= xd
     
   End Subroutine LSQ_solve
   
 
   Subroutine Sym_solve (Factor, A,x)
   
!    Compute the solution to A.y = x, doing the Cholesky factorization of A.

     Logical,          Intent(in) :: Factor  ! Factor matrix 'A' if true, else assume already factored
     Real(Dblp), Intent(inout) :: A(:,:)  ! (n,n+1) The symmetric matrix
     Real(Dblp), Intent(inout) :: x(:)    ! (n)   Right hand side (input). Solution vector (output)
!  Local:
     Real(Dblp), Parameter :: eps= 1.0E-7
     Integer :: err_col
     
     If (Factor) Call Cholesky (.true., .true., eps, A, err_col)
     Call Substitution (.true., .true., .true., A,x)
     
   End Subroutine Sym_solve
 
!   The next procedure does the Cholesky factorization of the 
!   matrix equation mat, assuming that 'mat' is
!   symmetric positive definite.

!   It is a translation of the Cholesky factorization routine
!   'spofa' from LINPACK, followed by the forward and backward
!   solver 'sposl'. In the non-standard form, the factorization
!   is R.R^t instead of R^t.R, and the solver applies back 
!   substitution first.

!   Another option is to do row and column scaling by the inverse
!   square root of the diagonal.  In this case the matrix must 
!   have an extra column to store these values.
!
!   Only the upper triangle of 'mat' is used.

!   Cost 2n(n + 1) + (if factor) n[n^2 + (9*n+1)/2] / 3

!   Standard : Standard Cholesky factorization (R^t.R), else mat= R.R^t.
!   Diagonal : Do row and column scaling by the inverse square root of the diagonal.

!   If a computed diagonal element is not positive, then it is
!   replaced by 'eps' and 'err_col' is set to the column number,
!   which should trigger an error message.
!   If the computed diagonal element is positive but less than
!   'eps', then 'err_col' is set to the negative of the
!   column number, which should trigger a warning message.

!   Success if err_col = 0.

  Subroutine Cholesky (Standard, Diagonal, eps, mat, err_col)
    Logical,           Intent(in) :: Standard ! R^t.R vs. R.R^t factorization
    Logical,           Intent(in) :: Diagonal ! Do diagonal scaling
    Real(Dblp),     Intent(in) :: eps      ! Replaces negative diagonal values.
                                              ! tests positive diagonal values.
    Real(Dblp),  Intent(inout) :: mat(:,:) ! (mm,mm+1)  Only upper triangle is used. 
                                              !            Last column used only for diagonal scaling
    Integer,          Intent(out) :: err_col  ! Last col # where 'mat'
                                              ! is not or almost not pos def.
! Local:
    Integer       :: mm, j, k, dg
    Real(Dblp) :: sm, eps_sq

    mm= Size(mat,1);  dg= Size(mat,2)
    If (mm < 1 .or. dg <= mm) then
      err_col= dg;  Call out ("Error in Cholesky: Matrix dimensions too small.");  Return
    End if
    
    err_col= 0;  eps_sq= eps**2

    If (Diagonal) then  ! Do row and column scaling
      dg= mm + 1
      If (Size(mat,2) < dg) then
        err_col= dg;  Call out ("Error in Cholesky: Matrix dimensions too small.");  Return
      End if

      Do j= 1,mm
        If (mat(j,j) < eps_sq) then
          mat(j,j)= eps;  err_col= j
        End if

        mat(j,dg)= one / Sqrt(mat(j,j))
        mat(:mm,j)= mat(j,dg) * mat(:mm,j)  ! Scale the column
        mat(j,:mm)= mat(j,dg) * mat(j,:mm)  ! Scale the row
      End do
    End if

    If (Standard) then  ! Cholesky Factorization: mat= R^t.R
      If (mat(1,1) <= eps_sq) then
        err_col= 1;  mat(1,1)= eps
      Else
        mat(1,1)= Sqrt(mat(1,1))
      End if
      mat(2:mm,1)= zero

      Do j= 2,mm        ! column
        mat(1,j)= mat(1,j) / mat(1,1)
        Do k= 2,j-1     ! row
          mat(k,j)= (mat(k,j) - Dot_product(mat(1:k-1,k), &
                     mat(1:k-1,j))) / mat(k,k)
        End do
        sm= mat(j,j) - Sum(mat(1:j-1,j)**2)

        If (sm <= eps_sq) then
          err_col= j;  sm= eps
        Else
          sm= Sqrt(sm)
        End if
        mat(j,j)= sm;  mat(j+1:mm,j)= zero
      End do

    Else                ! Cholesky Factorization: mat= R.R^t
      If (mat(mm,mm) <= eps_sq) then
        err_col= mm;  mat(mm,mm)= eps
      Else
        mat(mm,mm)= Sqrt(mat(mm,mm))
      End if
      mat(mm,1:mm-1)= zero

      Do j= mm-1,1,-1   ! row
        mat(j,mm)= mat(j,mm) / mat(mm,mm)
        Do k= mm-1,j+1,-1 ! column
          mat(j,k)= (mat(j,k) - Dot_product(mat(k,k+1:mm), &
                     mat(j,k+1:mm))) / mat(k,k)
        End do
        sm= mat(j,j) - Sum(mat(j,j+1:mm)**2)

        If (sm <= eps_sq) then
          err_col= j;  sm= eps
        Else
          sm= Sqrt(sm)
        End if
        mat(j,j)= sm;  mat(j,1:j-1)= zero
      End do
    End if

    If (err_col > 0)  Call hherr (0,"Cholesky epsilon violation at column",err_col)
  End subroutine Cholesky


  Subroutine Substitution (Standard, Both, Diagonal, mat,vec)

!   Assume that mat = R^t.R if Standard is true, else mat= R.R^t.
!   Then compute vec= R^(-1).R^(-t).vec in the Standard case,
!   and vec= R^(-t).R^(-1).vec in the non-standard case.

!   However, if Both is false, then compute only R^(-t).vec
!   [forward substitution] in the Standard case, and R^(-1).vec
!   [back substitution] in the non-Standard case.

!   Special note: In the Diagonal case, the substitution is proceeded
!   by diagonal scaling for the Standard case, but not for the Non-standard
!   case. This enables you to start with forward substitution only
!   (Standard = .true. and Both = .false.), then later finish 
!   with backward substituion (Standard = .false. and Both= .false.).
!   However, if you were truly working with the non-Standard 
!   factorization, then you would have to make sure the scalings
!   are correct by hand.

    Logical,           Intent(in) :: Standard  ! 'mat' factored as R^t.R vs R.R^t
    Logical,           Intent(in) :: Both      ! Do both substitutions, else only the first one
    Logical,           Intent(in) :: Diagonal  ! Assume diagonal scaling
    Real(Dblp),     Intent(in) :: mat(:,:)  ! (>=m,>=m+1) Last column used only for diagonal scaling
    Real(Dblp),  Intent(inout) :: vec(:)    ! (m)
! Local:
    Integer       :: mm, k, dg
    Real(Dblp) :: sm

    mm= Size(vec);  If (mm < 1) Return

    If (Standard) then

      If (Diagonal) then
        dg= mm + 1;  vec= vec * mat(:mm,dg)
      End if

!     Forward solve (cost= n^2 + 2n -1)

      vec(1)= vec(1) / mat(1,1)
      Do k= 2,mm
        sm= Dot_product(mat(1:k-1,k), vec(1:k-1))
        vec(k)= (vec(k) - sm) / mat(k,k)
      End do

!     Backward solve (cost = n^2)

      If (Both) then
        Do k= mm,2,-1
          vec(k)= vec(k) / mat(k,k)
          vec(1:k-1)= vec(1:k-1) - vec(k) * mat(1:k-1,k)
        End do
        vec(1)= vec(1) / mat(1,1)
        If (Diagonal)  vec(:mm)= vec(:mm) * mat(:mm,dg)
      End if

    Else

!     Backward solve (cost = n^2)

      Do k= mm,2,-1
        vec(k)= vec(k) / mat(k,k)
        vec(1:k-1)= vec(1:k-1) - vec(k) * mat(1:k-1,k)
      End do
      vec(1)= vec(1) / mat(1,1)

!     Forward solve (cost= n^2 + 2n -1)

      If (Both) then
        vec(1)= vec(1) / mat(1,1)
        Do k= 2,mm
          sm= Dot_product(mat(1:k-1,k), vec(1:k-1))
          vec(k)= (vec(k) - sm) / mat(k,k)
        End do
      End if

      If (Diagonal) then
        dg= mm + 1
        vec(:mm)= vec(:mm) * mat(:mm,dg)
      End if
    End if
  End subroutine Substitution


  Subroutine Multiply_RR (Standard,Symmetrize, kk, mm, mat_in, mat_out)

!   Multiply to get mat_out= mat_out + R^t.R if Standard, else
!                   mat_out= mat_out + R.R^t,
!   where R = kk by mm upper triangle of mat_in (for kk <= mm).

!   In the Standard case mat_out is equivalent to taking
!   dot products between the columns of R:
!   mat_out(i,j)= Dot_product(mat_in(:,i), mat_in(:,j))

!   In the non-Standard case mat_out is equivalent to taking
!   dot products between the rows of R:
!   mat_out(i,j)= Dot_product(mat_in(i,:), mat_in(j,:))

!   Thus mat_out is mm by mm in the Standard case and
!   kk by kk in the non-Standard case.
   
    Logical,           Intent(in) :: Standard     ! R^t.R vs. R.R^t
    Logical,           Intent(in) :: Symmetrize   ! Symmetrize the output matrix
    Integer,           Intent(in) :: kk, mm       ! Dimensions of 'mat'
    Real(Dblp),     Intent(in) :: mat_in(:,:)  ! (kk,mm) Only upper triangle is used
    Real(Dblp), Intent(inout) :: mat_out(:,:) ! (kk,mm) Upper triangle only unless symmetrized
! Local:
    Integer :: i, j, jk

    If (Standard) then  ! mat_out= mat_out + mat_in^t.mat_in
      Do i= 1,mm
        Do j= 1,Min(i,kk)
          mat_out(i,i:mm)= mat_out(i,i:mm) + mat_in(j,i) * mat_in(j,i:mm)
        End do
      End do

      If (Symmetrize) then
        Do i= 1,mm-1
          mat_out(i+1:mm,i)= mat_out(i,i+1:mm)
        End do
      End if

    Else                ! mat_out= mat_out + mat_in.mat_in^t
      Do i= 1,kk
        Do j= mm,i,-1
          jk= Min(j,kk)
          mat_out(i,i:jk)= mat_out(i,i:jk) + &
                           mat_in(i,j) * mat_in(i:jk,j)
        End do
      End do

      If (Symmetrize) then
        Do i= 1,kk-1
          mat_out(i+1:kk,i)= mat_out(i,i+1:kk)
        End do
      End if
    End if
  End Subroutine Multiply_RR


  Subroutine Multiply_by_R (Standard,Both, kk,mm, mat, vec)

!   If Standard and Both are true, then multiply 'vec' by R, then R^t,
!   assuming mat= R^t.R.  That is, vec= mat.vec= R^t.(R.vec).
!   If Both is false, multiply by R only.

!   If non-Standard, then vec= mat.vec= R.(R^t.vec) for mat= R.R^t.
!   If Both is false, then multiply by R^t only.

    Logical,           Intent(in) :: Standard  ! Start by multiplying vec by R vs. R^t
    Logical,           Intent(in) :: Both      ! Multiply by R^t.R vs. R.R^t
    Integer,           Intent(in) :: kk, mm    ! Dimensions of 'mat'
    Real(Dblp),     Intent(in) :: mat(:,:)  ! (kk,mm) Factored matrix
    Real(Dblp), Intent(inout) :: vec(:)    ! (kk) Vector
! Local:
    Integer :: i, ik

    If (Standard) then

!     Multiply by R (cost= kk*(2*mm + 1 - kk))

      Do i= 1,kk
        vec(i)= Dot_product(mat(i,i:mm), vec(i:mm))
      End do

!     Multiply by R^t (cost = kk*(2*mm + 1 - kk))

      If (Both) then
        Do i= mm,2,-1
          ik= Min(i,kk)
          vec(i)= Dot_product(mat(1:ik,i), vec(1:ik))
        End do
        vec(1)= mat(1,1) * vec(1)
      End if
    Else

!     Multiply by R^t (cost = kk*(2*mm + 1 - kk))

      Do i= mm,2,-1
        ik= Min(i,kk)
        vec(i)= Dot_product(mat(1:ik,i), vec(1:ik))
      End do
      vec(1)= mat(1,1) * vec(1)

!     Multiply by R (cost= kk*(2*mm + 1 - kk))

      If (Both) then
        Do i= 1,kk
          vec(i)= Dot_product(mat(i,i:mm), vec(i:mm))
        End do
      End if
    End if
  End subroutine Multiply_by_R


  Subroutine Invert_R (R_in, R_out)

!   Invert the upper triangular matrix R_in to get the
!   upper triangular matrix R_out.

    Real(Dblp),  Intent(in) :: R_in(:,:)    ! (nn,nn) Only upper triangle is used.
    Real(Dblp), Intent(out) :: R_out(:,:)   ! (nn,nn) Only upper triangle is computed.
! Local:
    Integer :: nn, i, j

     nn= Size(R_in,1)
     R_out= zero;  R_out(1,1)= one / R_in(1,1)

     Do i= 2,nn
       R_out(i,i)= one
       Do j= i,2,-1     ! Back substitutions
         R_out(j,i)= R_out(j,i) / R_in(j,j)
         R_out(1:j-1,i)= R_out(1:j-1,i) - R_out(j,i) * R_in(1:j-1,j)
       End do
       R_out(1,i)= R_out(1,i) / R_in(1,1)
     End do
  End Subroutine Invert_R


  Subroutine Invert_sym (Factored, Diagonal, mat_in, mat_out, err_col)

!   Invert the positive definite symmetric matrix mat_in to get the matrix mat_out.

    Logical,           Intent(in) :: Factored      ! Input matrix already factored
    Logical,           Intent(in) :: Diagonal      ! Do diagonal scaling (If factored, 'mat_in' must 
                                                   !   have been factored with Diagonal)
    Real(Dblp), Intent(inout) :: mat_in(:,:)    ! (nn,nn+1)  Input matrix: Only upper triangle is used.
                                                   !            Factored on output if not already factored.
    Real(Dblp),   Intent(out) :: mat_out(:,:)   ! (nn,nn)    Inverted matrix: Fully symmetrized.
    Integer,         Intent(out) :: err_col        ! Last col # where 'mat_in'
                                                   ! is not or almost not pos def.
! Local:
    Integer       :: nn, j, dg
    Real(Dblp) :: dum, eps= 1.0E-12, tmp(Size(mat_in,1),Size(mat_in,2))

    nn= Size(mat_in,1);  err_col= 0

    If (nn == 1) then
      mat_out(1,1)= one / mat_in(1,1)
    Else
      If (.not.Factored) then
        tmp= mat_in
        Call Cholesky (.true.,Diagonal, eps, mat_in, err_col)

        If (err_col > 0) then
          Call out ("Matrix inversion error in Invert_sym: Not positive definite")
          Call out ("Make it diagonally dominant")
          Do j= 1,nn
            dum= Sum(Abs(tmp(j,j+1:nn))) + Sum(Abs(tmp(1:j-1,j)))
            tmp(j,j)= Max(tmp(j,j), 1.001*dum)
          End do
          mat_in= tmp
          Call Cholesky (.true.,Diagonal, eps, mat_in, err_col)
        End if
      End if

      Call Invert_R (mat_in, mat_out)
 
      tmp(:nn,:nn)= mat_out(:nn,:nn);  mat_out= zero
      Call Multiply_RR (.false.,.true., nn, nn, tmp, mat_out)

      If (Diagonal) then
        dg= nn + 1
        Do j= 1,nn
          mat_out(:nn,j)= mat_in(j,dg) * mat_out(:nn,j)  ! Scale the column
          mat_out(j,:nn)= mat_in(j,dg) * mat_out(j,:nn)  ! Scale the row
        End do
      End if
    End if
  End Subroutine Invert_sym

  
  Subroutine Invert_non_sym (mat_in, mat_out)

!   Invert the non-symmetric matrix mat_in to get the matrix mat_out.

    Real(Dblp), Intent(inout) :: mat_in(:,:)   ! (nn,nn)  Input matrix
    Real(Dblp),   Intent(out) :: mat_out(:,:)  ! (nn,nn)  Inverted matrix
! Local:
    Integer :: i, nn, parity, indx(Size(mat_in,2))

    nn= Size(mat_in,2);   mat_out= zero;  Forall(i=1:nn) mat_out(i,i)= one

    Call LU_decompD (mat_in, indx, parity)

    Do i= 1,nn
      Call LU_solveD (mat_in, indx, mat_out(:,i))
    End do
  End Subroutine Invert_non_sym


  Subroutine LU_decompD (mat, indx, parity)
  
!   Do the LU decomposition of the matrix 'mat' as prelude to
!   either inverting 'mat' or solving mat.y= x

    Real(Dblp), Intent(inout) :: mat(:,:)  ! (nn,nn) Do LU decompostion on this matrix
	Integer,         Intent(out) :: indx(:)   ! (nn) Pivot columns
	Integer,         Intent(out) :: parity
! Local:
	Real(Dblp) :: vv(Size(mat,1))
	Real(Dblp) :: eps= 1.0E-14
	Integer       :: j, nn, imax, im(1)

    nn= Size(mat,2);  parity= 1
	vv= Maxval(Abs(mat),dim=2)

    If (Any(vv == zero)) then 
	  Call out('Singular matrix in LU_decompD');  Stop
    End if

  	vv= one / vv

	Do j= 1,nn
	  im= Maxloc(vv(j:nn)*Abs(mat(j:nn,j)))
   	  imax= (j-1) + im(1)
	  If (j /= imax) then
  	    Call Swap (mat(imax,:), mat(j,:))
	   parity= -parity;  vv(imax)= vv(j)
	  End if

  	  indx(j)= imax
	  If (mat(j,j) == zero) mat(j,j)= eps
	    mat(j+1:nn,j)= mat(j+1:nn,j) / mat(j,j)
      mat(j+1:nn,j+1:nn)= mat(j+1:nn,j+1:nn) - Outer_product(mat(j+1:nn,j), mat(j,j+1:nn))
  	End do
  End Subroutine LU_decompD

  
  Subroutine LU_solveD (mat, indx, vec)
  
!   Solve the matrix equation mat.x= y, assuming the LU decompostion of 'mat'
!   and x = 'vec' on input, with y = 'vec' on outut.
  
    Real(Dblp),     Intent(in) :: mat(:,:)  ! (nn,nn)  Assume the LU decomposition has been done
    Integer,           Intent(in) :: indx(:)   ! (nn)     Pivot columns from LU_decompD
    Real(Dblp),  Intent(inout) :: vec(:)    ! (nn)     RHS on input, solution on output
! Local:
    Integer       :: i, nn, ii, ll
    Real(Dblp) :: summ

    nn= size(mat,2);  ii= 0

    Do i= 1,nn
	  ll= indx(i);  summ= vec(ll);  vec(ll)= vec(i)

	  If (ii /= 0) then
		  summ= summ - Dot_product(mat(i,ii:i-1), vec(ii:i-1))
	  Else if (summ /= zero) then
	    ii= i
	  End if

	  vec(i)= summ
	End do

	Do i= nn,1,-1
	  vec(i)= (vec(i) - Dot_product(mat(i,i+1:nn), vec(i+1:nn))) / mat(i,i)
	End do
  End Subroutine LU_solveD

  Subroutine Mat_solveR1 (mat, vec)
  
!   Solve the matrix equation mat.x= y for a single vectors 'y'
!   doing the LU decompostion of 'mat'
  
    Real, Intent(inout) :: mat(:,:)  ! (n,n)  Input matrix, in LU form on output
    Real, Intent(inout) :: vec(:)    ! (n)    RHS on input, solution on output
! Local:
	Integer :: i, parity, indx(Size(vec,1))  ! (n) Pivot columns

    Call LU_decomp (mat, indx, parity)
    Call LU_solve (mat, indx, vec)
    
  End Subroutine Mat_solveR1

  Subroutine Mat_solveR2 (mat, vec)
  
!   Solve the matrix equation mat.x= y for multiple vectors 'y', 
!   doing the LU decompostion of 'mat'
  
    Real, Intent(inout) :: mat(:,:)  ! (n,n)  Input matrix, in LU form on output
    Real, Intent(inout) :: vec(:,:)  ! (n,nl) RHS on input, solution on output
! Local:
	Integer :: i, parity, indx(Size(vec,1))  ! (n) Pivot columns

    Call LU_decomp (mat, indx, parity)
    Do i= 1,Size(vec,2)
      Call LU_solve (mat, indx, vec(:,i))
    End do
  End Subroutine Mat_solveR2

  
  Subroutine Mat_solveD1 (mat, vec)
  
!   Solve the matrix equation mat.x= y for a single vectors 'y'
!   doing the LU decompostion of 'mat'
  
    Real(Dblp), Intent(inout) :: mat(:,:)  ! (n,n)  Input matrix, in LU form on output
    Real(Dblp), Intent(inout) :: vec(:)    ! (n)    RHS on input, solution on output
! Local:
	Integer :: i, parity, indx(Size(vec,1))  ! (n) Pivot columns

    Call LU_decomp (mat, indx, parity)
    Call LU_solve (mat, indx, vec)
    
  End Subroutine Mat_solveD1

  Subroutine Mat_solveD2 (mat, vec)
  
!   Solve the matrix equation mat.x= y for multiple vectors 'y', 
!   doing the LU decompostion of 'mat'
  
    Real(Dblp), Intent(inout) :: mat(:,:)  ! (n,n)  Input matrix, in LU form on output
    Real(Dblp), Intent(inout) :: vec(:,:)  ! (n,nl) RHS on input, solution on output
! Local:
	Integer :: i, parity, indx(Size(vec,1))  ! (n) Pivot columns

    Call LU_decomp (mat, indx, parity)
    Do i= 1,Size(vec,2)
      Call LU_solve (mat, indx, vec(:,i))
    End do
  End Subroutine Mat_solveD2

  Subroutine LU_decompR (mat, indx, parity)
  
!   Do the LU decomposition of the matrix 'mat' as prelude to
!   either inverting 'mat' or solving mat.y= x

    Real,  Intent(inout) :: mat(:,:)  ! (nn,nn) Do LU decompostion on this matrix
	Integer, Intent(out) :: indx(:)   ! (nn) Pivot columns
    Integer, Intent(out) :: parity
! Local:
	Real    :: vv(Size(mat,1)), eps= 1.0E-7
	Integer :: j, nn, imax, im(1)

    nn= Size(mat,2);  parity= 1
	vv= Maxval(Abs(mat),dim=2)

    If (Any(vv == 0.0)) then 
      Call out('Singular matrix in LU_decompD');  Stop
    End if

  	vv= 1.0 / vv

	Do j= 1,nn
	  im= Maxloc(vv(j:nn)*Abs(mat(j:nn,j)))
   	  imax= (j-1) + im(1)
	  If (j /= imax) then
  	    Call Swap (mat(imax,:), mat(j,:))
	    parity= -parity;  vv(imax)= vv(j)
	  End if

  	  indx(j)= imax
	  If (mat(j,j) == 0.0) mat(j,j)= eps
	  mat(j+1:nn,j)= mat(j+1:nn,j) / mat(j,j)
      mat(j+1:nn,j+1:nn)= mat(j+1:nn,j+1:nn) - Outer_product(mat(j+1:nn,j), mat(j,j+1:nn))
  	End do
  End Subroutine LU_decompR

  
  Subroutine LU_solveR (mat, indx, vec)
  
!   Solve the matrix equaion mat.x= y, assuming the LU decompostion of 'mat'
!   and x = 'vec' on input, with y = 'vec' on outut.
  
    Real,     Intent(in) :: mat(:,:)  ! (nn,nn)  Assume the LU decomposition has been done
    Integer,  Intent(in) :: indx(:)   ! (nn)     Pivot columns from LU_decomp
    Real,  Intent(inout) :: vec(:)    ! (nn)     RHS on input, solution on output
! Local:
    Integer :: i, nn, ii, ll
    Real    :: summ

    nn= size(mat,2);  ii= 0

    Do i= 1,nn
	  ll= indx(i);  summ= vec(ll);  vec(ll)= vec(i)

	  If (ii /= 0) then
		summ= summ - Dot_product(mat(i,ii:i-1), vec(ii:i-1))
      Else if (summ /= 0.0) then
	    ii= i
	  End if

	  vec(i)= summ
	End do

	Do i= nn,1,-1
	  vec(i)= (vec(i) - Dot_product(mat(i,i+1:nn), vec(i+1:nn))) / mat(i,i)
	End do
  End Subroutine LU_solveR

  
   Subroutine RR_solveR (opt,RR, bb, ier)
  
!    Routine for forward and backward substitution, assuming 
!    that RR is an upper triangular matrix with non-zero 
!    diagonal elements and that the forward substitution is 
!    for the transpose of RR, which is lower triangular.

     Integer,  Intent(in)    :: opt      ! opt <= 0 for forward substituion, >= 0 for backward substitution
     Real,     Intent(in)    :: RR(:,:)  ! (kk,kk)
     Real,     Intent(inout) :: bb(:)    ! (kk) Right hand side
     Integer,  Intent(out)   :: ier
!   Local:
     Integer :: j, kk
  
     kk= size(RR,2)
     Do ier= 1,kk
       If (RR(ier,ier) == 0) then
         Call Out ("RR_solv failure due to non-zero diagonal",ier, "of value",RR(ier,ier), ln=1)
         Return
       End if
     End do
     ier= 0
  
     If (opt <= 0) then   ! Forward substitution   
       bb(1)= bb(1)/RR(1,1)
       Do j= 2,kk
         bb(j)= ( bb(j) - Dot_product(RR(:j-1,j), bb(:j-1)) ) / RR(j,j)
       End do
     End if
     
     If (opt >= 0) then  !  Backward substitution
       bb(kk)= bb(kk)/RR(kk,kk)
       Do j= kk-1,1,-1
         bb(:j)= bb(:j) - bb(j+1) * RR(:j,j+1)
         bb(j)= bb(j)/RR(j,j)
       End do
     End if
   End Subroutine RR_solveR
  
   Subroutine RR_solveD (opt,RR, bb, ier)
  
!    Routine for forward and backward substitution, assuming 
!    that RR is an upper triangular matrix with non-zero 
!    diagonal elements and that the forward substitution is 
!    for the transpose of RR, which is lower triangular.

     Integer,        Intent(in)    :: opt      ! opt <= 0 for forward substituion, >= 0 for backward substitution
     Real(Dblp),  Intent(in)    :: RR(:,:)  ! (kk,kk)
     Real(Dblp),  Intent(inout) :: bb(:)    ! (kk) Right hand side
     Integer,        Intent(out)   :: ier
!   Local:
     Real    :: dg
     Integer :: j, kk
  
     kk= size(RR,2)
     Do ier= 1,kk
       dg= RR(ier,ier)
       If (dg == 0) then
         Call Out ("RR_solv failure due to non-zero diagonal",ier, "of value",dg, ln=1)
         Return
       End if
     End do
     ier= 0
  
     If (opt <= 0) then   ! Forward substitution   
       bb(1)= bb(1)/RR(1,1)
       Do j= 2,kk
         bb(j)= ( bb(j) - Dot_product(RR(:j-1,j), bb(:j-1)) ) / RR(j,j)
       End do
     End if
     
     If (opt >= 0) then  !  Backward substitution
       bb(kk)= bb(kk)/RR(kk,kk)
       Do j= kk-1,1,-1
         bb(:j)= bb(:j) - bb(j+1) * RR(:j,j+1)
         bb(j)= bb(j)/RR(j,j)
       End do
     End if
   End Subroutine RR_solveD
  
  Subroutine Normalized_distance (Factored, Both, cov, fx, func)
    Logical,                  Intent(in) :: Factored  ! True if cov is already Cholesky factored
    Logical,                  Intent(in) :: Both      ! Do both forward and backward substitution
    Real(Dblp),         Intent(inout) :: cov(:,:)  ! (n,n+1) Covariance matrix
    Real(Dblp),         Intent(inout) :: fx(:)     ! (n) Residual (full residual fx.cov^(-1) unless func < 0)
    Real(Dblp), Optional, Intent(out) :: func      ! Squared normalized distance, based on cov.
! Local:
    Real(Dblp), Parameter :: eps= 1.0E-12
    Logical, Parameter :: Diagonal= .false.
    Logical :: Standard
    Integer :: err_col

    Standard= .true.
    If (.not.Factored) Call Cholesky (Standard, Diagonal, eps, cov, err_col)

    If (Present(func)) then
      Call Substitution (Standard, .false., Diagonal, cov, fx)
      func= Sum(fx**2)

      If (Both) then
        Standard= .false.
        Call Substitution (Standard, .false., Diagonal, cov, fx)
      End if
    Else
      Call Substitution (Standard, Both, Diagonal, cov, fx)
    End if

  End subroutine Normalized_distance 
  

  Subroutine Min_norm (n,m, A,y, x, dif,itr)
  
!   Compute the minimum norm solution of the under-determined linear
!   system A.x = y via the Kaczmarz method of iteration.  

    Integer,  Intent(in) :: n       ! # rows of A = size of y
    Integer,  Intent(in) :: m       ! # columns of A = size of x, so m > n
    Real,     Intent(in) :: A(:,:)  ! (n,m) matrix
    Real,     Intent(in) :: y(:)    ! (n) data vector
    
    Real,    Intent(out) :: x(:)    ! (m) min norm solution vector
    Real,    Intent(out) :: dif     ! final difference 
    Integer, Intent(out) :: itr     ! # outer iterations   
    
! Local:
    Real,    Parameter :: lam= 0.5, eps= 0.001
    Integer, Parameter :: itmax= 100
    Real    :: fac, x0(m), x1(m), B(n,m)
    Integer :: i, j, it, itm
    
    Do i= 1,n
      fac= lam / Sum(A(i,:)**2)
      B(i,:)= fac * A(i,:)
    End do
    
    x0= 0;  x1= x0;  i= 1;  itr= 0;  itm= itmax * n
    
    Do it= 1,itm
      fac= y(i) - Dot_product(A(i,:),x0)
      x= x1 + fac * B(i,:)
      
      If (i >= n) then
        i= 1;  itr= itr + 1
        dif= Sqrt(Sum((x - x0)**2))
        If (dif < eps) Exit
        x0= x
      Else
        i= i + 1
      End if
      
      x1= x
    End do
  End Subroutine Min_norm
  

    Subroutine QR_fac (plm, RR)
  
  !    Do the QR factorization, using Givens rotations, of the matrix
  !    RE = upper triangular matrix RR over the diagonal matrix
  !    Sqrt(plm)*I.  That is, RE is a 2k by k matrix.
  !    RR is overwritten.
  
       Real(Dblp),    Intent(in) :: plm      ! The LM parameter: added to the diagonal of the Krylov Hessian
       Real(Dblp), Intent(inout) :: RR(:,:)  ! (k,k)

  !    plm  Nonnegative constant.
  !    RR   An upper triangular matrix on input and output.  The lower
  !     triangle is used as work space.
  
  !   Local:
       Integer       :: i, j, i1, k
       Real(Dblp) :: sqlm, cs, sn, t, wr(size(RR,2))

       k= size(RR,2);  sqlm= Sqrt(plm)
  
       Do i= 1,k
  
  !      Apply Givens rotations to the ith row of the upper triangle
  !      of the lower half of RE in order to zero it out, using the
  !      ith column of the lower triangle of RR for the actual storage.
  
         i1= i + 1;  t= RR(i,i)
         
         Do j= 1,i-1
           t= Hypot(t,RR(i,j))
           cs= RR(i,i) / t;  sn= RR(i,j) / t
           RR(i,i)= t;  RR(i,j)= 0
  
           wr(i1:)  =  cs * RR(i,i1:) + sn * RR(i1:,j)
           RR(i1:,j)= -sn * RR(i,i1:) + cs * RR(i1:,j)
           RR(i,i1:)= wr(i1:)
         End do
  
  !      Likewise zero the diagonal of the lower half of RE
  
         t= Hypot(t,sqlm)
         cs= RR(i,i) / t;  sn= sqlm / t
         RR(i,i)= t
         RR(i1:,i)= -sn * RR(i,i1:)
         RR(i,i1:)=  cs * RR(i,i1:)
       End do
       
     End Subroutine QR_fac

    
   Function Cube_fit (t0,t1, f0,f1, df0,df1)

!    Fit a cubic function to the data (t,f(t),df(t)) for t= t0,t1.
!    Then solve for the minimum value in the interval (t0,t1).

     Real(Dblp), Intent(in) :: t0, t1, f0, f1, df0, df1

     Real(Dblp) :: Cube_fit
     Real(Dblp) :: rc(4), Hc(4,4), r4(4), H4(4,4), rq(3), Hq(3,3), r3(3), H3(3,3)
     Real(Dblp) :: t, dt, dft, dfm, dsc, t02, t03, t12, t13

!    Default solution = the endpoint with the smaller value
     
     Cube_fit= t0;  If (f1 < f0) Cube_fit= t1
     
     If (t1 < t0 + 0.001 .or. (df0 == zero .and. df1 <= zero) &
                         .or. (df0 >= zero .and. df1 == zero)) Return
     dfm= (f1 - f0) / (t1 - t0)
     
     If (df0 == zero .and. df1 > zero) then  ! Quadratic fit to f0, f1, & df1
       Hq(1,:)= (/one, t0, t0**2/)
       Hq(2,:)= (/one, t1, t1**2/)
       Hq(3,:)= (/zero, one, 2*t1/)
       rq= (/f0, f1, df1/)
       
       If (Test_LU) then
         H3= Hq;  Call Mat_solve (Hq,rq);  r3= Matmul(H3,rq)
       Else
         Call Mat_solve (Hq,rq)
       End if
       
       If (rq(3) > 0.0001) then
         t= -rq(2) / (2*rq(3))
         If (Test_LU) dft= rq(2) + (2*rq(3))*t
         
         If (t >= t0 .and. t <= t1) Cube_fit= t
       End if 
       Return
     Else if (df0 < zero .and. df1 == zero) then  ! Quadratic fit to f0, f1, & df0
       Hq(1,:)= (/one, t0, t0**2/)
       Hq(2,:)= (/one, t1, t1**2/)
       Hq(3,:)= (/zero, one, 2*t0/)
       rq= (/f0, f1, df0/)
       
       If (Test_LU) then
         H3= Hq;  Call Mat_solve (Hq,rq);  r3= Matmul(H3,rq)
       Else
         Call Mat_solve (Hq,rq)
       End if
       
       If (rq(3) > 0.0001) then
         t= -rq(2) / (2*rq(3))
         If (Test_LU) dft= rq(2) + (2*rq(3))*t
         
         If (t >= t0 .and. t <= t1) Cube_fit= t
       End if 
       Return
     End if
     
!    Fit a cubic polynomial f(t)= a + b*t + c*t^2 + d*t^3 to the data
!    where f'(t)= b + 2*c*t + 3*d*t^2 and f''(t)= 2*c + 6*d*t
     
     t02= t0*t0;  t03= t02*t0
     t12= t1*t1;  t13= t12*t1
     
     Hc(1,:)= (/one, t0, t02, t03/)
     Hc(2,:)= (/one, t1, t12, t13/)
     Hc(3,:)= (/zero, one, 2*t0, 3*t02/)
     Hc(4,:)= (/zero, one, 2*t1, 3*t12/)
     rc= (/f0, f1, df0, df1/)
     
     If (Test_LU) then
       H4= Hc;  Call Mat_solve (Hc,rc);  r4= Matmul(H4,rc)
     Else
       Call Mat_solve (Hc,rc)
     End if
     
     dsc= rc(3)**2 - 3*rc(2)*rc(4)  ! Discriminant of the quadratic
     
     If (dsc >= zero .and. Abs(rc(2)) > 0.0001) then
       dsc= Sqrt(dsc)
       t= (-rc(3) + dsc)/(3*rc(4))
       dt= 2*rc(3) + (6*rc(4))*t
       If (Test_LU) dft= rc(2) + (2*rc(3))*t + (3*rc(4))*t**2
       
       If (dt >= zero .and. t0 <= t .and. t <= t1) then
         Cube_fit= t; Return
       End if
       
       t= (-rc(3) - dsc)/(3*rc(4))
       dt= 2*rc(3) + (6*rc(4))*t
       If (Test_LU) dft= rc(2) + (2*rc(3))*t + (3*rc(4))*t**2
       
       If (dt >= zero .and. t0 <= t .and. t <= t1) then
         Cube_fit= t; Return
       End if
     End if
     
   End function Cube_fit
   

   Subroutine Line_search (parm, nfn,tf, ier, Func_t)

  !  Do a line search to minimize the positive function f(t) between a lower
  !  bound t_low and upper bound t_high on t, starting at t0.
  !  Subroutine 'Func_t' is called to do the function evaluation.

  !  If Move_left is true, then the function is first computed at
  !  t0 - step until the minimum has been bracketed or t_low reached
  !  or the 'change' tolerance satisfied.  Derivatives are computed at
  !  all of these steps.  If a derivative is close enough to zero according
  !  to the 'change', tolerance, then convergence is established.

  !  If any of the above procedures halt due to tolerance satisfaction,
  !  there is an immediate return.  If the halt is due to bracketing,
  !  then quadratic interpolation is done until convergence by either
  !  the change in f ('change') or change in t ('tol') tolerance.

  !  If instead Move_left is false, then the movement is to the right,
  !  letting t= t0 + step, until a minimum is bracketed or t_high
  !  is reached.  The convergence criteria are the same.

  !  If the number of function iterations reaches the limit 'mxfun'
  !  without convergence or it stops at the lower bound t_low,
  !  the routine returns with a warning message.

  !  The function is assumed to have a single local = global minimum
  !  on the interval t_low .. t0 if Move_left, otherwise on on t0..t_high.
  !  An ordered graph of the computed trajectory is output, including
  !  derivatives.

  !  It is assumed that the function Func_t is defined in terms of an
  !  underlying vector function fcn.

     Type(Srch_parm), Intent(in) :: parm

  !  dmp_opt Damping option: Minimize func(dmp_opt)
  !  eps     Finite difference epsilon
  !  step    Step size for 't'.
  !  t_low   Lower bound for 't'.
  !  t0      Initial value for 't'.
  !  t_high  Upper bound for 't'.
  !  tol     Convergence tolerance for steps in 't'.
  !  change  Convergence tolerance for fractional change in f(t).

  !  Func_t     Function evaluation subroutine.  Usage:

  !       Call Func_t (t, ft, ier)

  !       to compute the function values ft(1:,1) at 't' and
  !       their derivatives ft(1:,2), for values = functional
  !       its gradient, and its residual (= pre-conditioned gradient)
  
  !  fcn  Subroutine called by Func_t to evaluate the underlying vector

     Integer,    Intent(out) :: nfn          ! Current function evaluation # = final index for 'tf' on output.
     Real(Dblp), Intent(out) :: tf(0:,:,0:)  ! (0:3,2,0:mxfun) Trajectory function evaluation data - see below

     Integer,    Intent(out) :: ier          ! Error code. < 0 : error
                                             !   0   : warning (need to revise dimensions or parameters)
                                             !   > 0 : success
     External    :: Func_t       ! Underlying function evaluation procedure
     
  !  nfn         = final index for 'tf' on output = # function calls - 1
     
  !  tf(0,1,n)   = t
  !  tf(1,1:2,n) = functional and its derivative (opt = 1)
  !  tf(2,1:2,n) = other value such as gradient and its derivative (optional: opt = 2)
  !  tf(3,1:2,n) = other value such as residual and its derivative (optional: opt = 3)

  !  ier   Success/error code.
  !     ier < 0   Func_t internal error
  !     ier = 0   Warning: Need to change algorithm parameters or other characteristics
  !     ier = 1   Warning: Did not convergence
  !     ier = 2   Successful convergence

  ! Local:
     Real(Dblp) :: ftnew, dft, ftold, cft, told, ft(0:3,2)
     Real(Dblp) :: t, t0, t1, t1s, t2s, t1d, t2d, t_next, step, dif_t, dif_f, dftt
     Logical       :: Do_cubic
     Integer       :: k, k1, opt, mxfun

     opt= parm%dmp_opt;  mxfun= Ubound(tf,3)
     
     t1s=  parm%t_low  + parm%step
     t2s=  parm%t_high - parm%step
     t1d=  parm%t_low  + parm%tol
     t2d=  parm%t_high - parm%tol
     
     nfn= 0;  ft= 0;  t0= parm%t0
     Call Func_t (t0, ft, ier);  tf(:,:,0)= ft
     If (ier <= 0) then
       Return
     End if  
     
     ftold= ft(opt,1);  dft= ft(opt,2)
     
     cft= parm%change * ftold
     If (Abs(dft) < cft) then ! Convergence
       ier= 2;  Return
     End if
     
     Restart_loop : Do 

       Do_cubic= .false.;  step= parm%step
     
       If (dft > 0) then

!        Bracket minimum function value to the left

         told= t0;  t= Max(t0 - step, (t0 + parm%t_low)/2, t1d)

         Left_loop : Do while (nfn < mxfun)
           k1= nfn;  nfn= nfn + 1;  k= nfn

           t1= t
           Call Func_t (t1, ft, ier);  tf(:,:,nfn)= ft
           If (ier <= 0) then
             Return
           End if

           ftnew= ft(opt,1);  dft= ft(opt,2)
           dif_t= Abs(told - t);  dif_f= Abs(ftnew - ftold)
           dftt= Abs(dft)*dif_t;  cft= parm%change * ftnew
           
           If (Min(dif_f, dftt) < cft) then  ! Convergence
             ier= 2;  Return
           End if
         
           If (t1 /= t) then             ! Restart the line search
             t0= t1;  Cycle Restart_loop
           End if
         
           If (dft < 0 .and. Max(ftnew, dftt) < 5*ftold) then  ! Interpolate
             Do_cubic= .true.;  Exit Left_loop
           Else if (ftnew > ftold) then   ! Anomalous result: Backup
             t= (t + told) / 2;  Cycle Left_loop
           Else if (t <= t1d) then       ! Failed to converge
             ier= 1;  Return
           End if

           told= t;  ftold= ftnew
           
           If (t <= t1s) then  ! Decelerate the reduction in t to t_low
             t= (parm%t_low + t)/2
           Else                ! Accelerate the reduction in t to t1s
             t_next= t - step
             If (t_next > t1s) then
               t= t_next;  step= 1.5 * step
             Else
               t= t1s
             End if
           End if
         End do Left_loop

       Else  !  Bracket minimum function value to the right

         told= t0;  t= Min(t0 + step, (t0 + parm%t_high)/2, t2d)
       
         Right_loop : Do while (nfn < mxfun)
           k= nfn;  nfn= nfn + 1;  k1= nfn
           t1= t
           Call Func_t (t1, ft, ier);  tf(:,:,nfn)= ft
           If (ier <= 0) then
             Return
           End if

           ftnew= ft(opt,1);  dft= ft(opt,2)
           dif_t= Abs(t - told);  dif_f= Abs(ftnew - ftold)
           dftt= Abs(dft)*dif_t;  cft= parm%change * ftnew

           If (Min(dif_f, dftt) < cft) then  ! Convergence
             ier= 2;  Return
           End if
         
           If (t1 /= t) then              ! Restart the line search
             t0= t1;  Cycle Restart_loop
           End if
         
           If (dft > 0 .and. Max(ftnew, dftt) < 5*ftold) then  ! Interpolate
              Do_cubic= .true.;  Exit Right_loop
           Else if (ftnew > ftold) then   ! Anomalous result: Backup
             t= (t + told) / 2;  Cycle Right_loop
           Else if (t >= t2d) then        ! Failed to converge
             ier= 1;  Return
           End if

           told= t;  ftold= ftnew
           
           If (t >= t2s) then  ! Decelerate the increase in t to t_high
             t= Min((t + parm%t_high)/2, t2d)
           Else               ! Accelerate the increase in t to t2s
             t_next= t + step
             If (t_next < t2s) then
               t= t_next;  step= 1.5 * step
             Else
               t= t2s
             End if
           End if
         End do Right_loop
       End if
     
       If (.not.Do_cubic) then  ! Failed to converge
         ier= 1;  Return
       End if

!      Do cubic interpolation

       Cubic_loop : Do while (nfn < mxfun)
         told= t
         t= Cube_fit (tf(0,1,k),tf(0,1,k1), tf(opt,1,k),tf(opt,1,k1), tf(opt,2,k),tf(opt,2,k1))
         dif_t= Abs(t-told)
         
         If (dif_t < parm%tol) then
           ier= 2;  Return    ! Convergence
         End if

         nfn= nfn + 1;  t1= t
         Call Func_t (t1, ft, ier);  tf(:,:,nfn)= ft
         If (ier <= 0) then
           Return
         End if
       
         ftnew= ft(opt,1);  dft= ft(opt,2)
         dif_f= Abs(ftnew-ftold);  dftt= Abs(dft)*dif_t
         cft= parm%change * ftnew
       
         If (Min(dif_f, dftt) < cft) then  ! Convergence
           ier= 2;  Return
         End if 
       
         If (t1 /= t) then                 ! Restart the line search
           t0= t1;  ftold= ftnew
           Cycle Restart_loop
         End if
       
         If (dft > 0) then                 ! Continue
           k1= nfn
         Else
           k = nfn
         End if

         ftold= ftnew
       End do Cubic_loop
     
       Exit Restart_loop
     End do Restart_loop 
        
     ier= 1    ! Failed to converge
   End Subroutine Line_search
   
End module Cholesk
