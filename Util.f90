!             R.H. Burkhart  1-20-96

Module Util

   Use Types
   Use Output
   Use Constant
   Use Precisn
   Implicit None

   Interface Cos_rise_drv
     Module procedure Cos_rise_drvR, Cos_rise_drvD
   End Interface

   Interface Dot_product_M
     Module procedure Dot_product_MR, Dot_product_MD
   End Interface

   Interface Dot_product_M1
     Module procedure Dot_product_M1R, Dot_product_M1D
   End Interface

   Interface Dif_angle
     Module procedure Dif_angleR, Dif_angleD
   End Interface

   Interface First_true
     Module procedure First_true_in, First_true_from
   End Interface

   Interface Last_true
     Module procedure Last_true_in, Last_true_to
   End Interface

   Interface List_of_true
  	 Module procedure List_of_true0, List_of_true1
   End Interface

   Interface Binary_search
     Module procedure Binary_search_int, Binary_search_real, Binary_search_cut
   End Interface

   Interface Assignment (=)
     Module procedure Identity_vector_int, Identity_vector_real
   End Interface

   Interface Sort  ! Normally Quick_Sort_int is not 'stable' - 
                   ! does not preserve the prior orderings 
                   ! of equal items, but we give it a stable option
     Module procedure Quick_Sort_stable, Quick_Sort_real, Quick_Sort_double
   End Interface

   Interface Matrix_diagonal
     Module procedure Matrix_diagonal_op_const_int, Matrix_diagonal_op_const_real, Matrix_diagonal_op_const_double
     Module procedure Matrix_diagonal_op_vec_int, Matrix_diagonal_op_vec_real, Matrix_diagonal_op_vec_double
     Module procedure Matrix_diagonal_out_vec_int, Matrix_diagonal_out_vec_real, Matrix_diagonal_out_vec_double
   End Interface

   Interface Swap
  	 Module procedure Swap_I,  Swap_R,  Swap_D
	 Module procedure Swap_IV, Swap_RV, Swap_DV
   End Interface

   Interface Outer_product
  	 Module procedure Outer_productR, Outer_productD
   End Interface
	
   Interface Trace
  	 Module procedure Trace_R, Trace_D
   End Interface
	
   Interface Output_angle
  	 Module procedure Output_angleR, Output_angleD
   End Interface

   Interface Array
  	 Module procedure Array_int, Array_real
   End Interface

   Interface Normalize_vec
  	 Module procedure Normalize_vec1_R, Normalize_vec1_D, Normalize_vec2_R, Normalize_vec2_D
   End Interface

   Interface Subset_of
  	 Module procedure Subset_ofI, Subset_ofR
   End Interface

   Interface Sum_constraint
     Module procedure Sum_constraint1, Sum_constraint2, &
                      Sum_constraint1D, Sum_constraint2D
   End Interface

   Interface Cubic_rise
     Module procedure Cubic_riseR, Cubic_riseD, Cubic_riseVR, Cubic_riseVD
   End Interface

   Interface Cos_cube
     Module procedure Cos_cubeVR, Cos_cubeVD
   End Interface

   Interface Cube_quad
     Module procedure Cube_quadR, Cube_quadV
   End Interface

   Interface Quad_rise1
     Module procedure Quad_rise1R, Quad_rise1D
   End Interface

   Interface Cos_rise
     Module procedure Cos_riseR, Cos_riseD, Cos_riseVD, Cos_riseVR, Cos_riseMR
   End Interface
   
   Interface Cos_fall
     Module procedure Cos_fallR, Cos_fallD, Cos_fallVR
   End Interface
   
   Interface Close_ordered
     Module procedure Close_ordered1, Close_ordered2
   End Interface

   Interface Diag_to_vec
     Module procedure Diag_to_vecR, Diag_to_vecI
   End Interface

   Interface Vec_to_diag
     Module procedure Vec_to_diagR, Vec_to_diagI
   End Interface

   Interface Const_to_diag
     Module procedure Const_to_diagI, Const_to_diagR, Const_to_diagD
   End Interface

   Interface Remove_from_matrix
     Module procedure Remove_from_matrixR, Remove_from_matrixI
   End Interface

   Interface Hypot
     Module procedure HypotR, HypotD
   End Interface

   Interface Vec_mean
     Module procedure Vec_meanR, Vec_meanD
   End Interface

   Interface Median
     Module procedure MedianR, MedianD, Median2R, Median2D
   End Interface

   Interface Partial_sums
     Module procedure Partial_sumsI, Partial_sumsR, Partial_sumsD
   End Interface

   Interface Select_rnd
     Module procedure Select_rndR, Select_rndD
   End Interface

  Contains
    
!   Modified Dot_prodcut to compute a modified correlation 'corr', 
!   assuming |yr| = 1 and that entries of 'xr' where both 'xr' and 'yr' 
!   are negative will be scaled by 'fac' and renormalized before 
!   computing the correlation. If 0 <= fac < 1 this has the effect of
!   decreasing the correlation, so that it depends more on the
!   correlation of positive entries. 
    
    Subroutine Dot_product_MR (fac, xr,yr, n,corr)

      Real,     Intent(in) :: fac    ! fac < 0 or fac >= 1:  Ordinary dot product of xr and yr
                                     ! fac >= 0: First scale xr by 'fac' to get x1 where both xr & yr
                                     !           are negative. Then normalize x1, so that 'corr' becomes
                                     !           a correlation if |yr| = 1
      Real,     Intent(in) :: xr(:)  ! (nc)  Assuming |xr| = 1 for a correlation
      Real,     Intent(in) :: yr(:)  ! (nc)  Assuming |yr| = 1 for a correlation
         
      Integer, Intent(out) :: n      ! # coordinates where both xr and yr are negative, n <= nc
      Real,    Intent(out) :: corr   ! The correlation Sum(x2*yr)
!   Local:      
      Real,    Allocatable :: x2(:)  ! (nc) 'xr' as modified
      Integer, Allocatable :: lst(:) ! (nc)

      Real, Parameter :: eps= 1.0E-6
      Real    :: x2norm
      Integer :: nc;  nc= Size(yr)
      
      If (fac < 0 .or. fac >= 1) then  ! Ordinary dot product 
        n= 0;  corr= Dot_product(xr,yr)
      Else
        Allocate(lst(nc))

        Call List_of_true (xr < 0 .and. yr < 0, n,lst)
          
        If (n <= 0) then
          corr= Dot_product(xr,yr)  ! Ordinary dot product
        Else
          Allocate(x2(nc));  x2= xr

          x2(lst(:n))= fac*x2(lst(:n))
            x2norm= Sqrt(Sum(x2**2) + eps)
          corr= Dot_product(x2,yr) / x2norm  ! Correlation

          DeAllocate(x2)
        End if
        DeAllocate(lst)
      End if

    End Subroutine Dot_product_MR

    
    Subroutine Dot_product_MD (fac, xr,yr, n,corr)

      Real(Dblp),  Intent(in) :: fac    ! fac < 0 or fac >= 1:  Ordinary dot product of xr and yr
                                           ! fac >= 0: First scale xr by 'fac' to get x1 where both xr & yr
                                           !           are negative. Then normalize x1, so that 'corr' becomes
                                           !           a correlation if |yr| = 1
      Real(Dblp),  Intent(in) :: xr(:)  ! (nc)  Assuming |xr| = 1 for a correlation
      Real(Dblp),  Intent(in) :: yr(:)  ! (nc)  Assuming |yr| = 1 for a correlation
         
      Integer,       Intent(out) :: n      ! # coordinates where both xr and yr are negative, n <= nc
      Real(Dblp), Intent(out) :: corr   ! The correlation Sum(x2*yr)
!   Local:      
      Real(Dblp), Allocatable :: x2(:)   ! (nc) 'xr' as modified
      Integer,       Allocatable :: lst(:)  ! (nc)

      Real(Dblp), Parameter :: eps= 1.0E-13
      Real(Dblp) :: x2norm
      Integer       :: nc;  nc= Size(yr)
      
      If (fac < 0 .or. fac >= 1) then  ! Ordinary dot product 
        n= 0;  corr= Dot_product(xr,yr)
      Else
        Allocate(lst(nc))

        Call List_of_true (xr < 0 .and. yr < 0, n,lst)
          
        If (n <= 0) then
          corr= Dot_product(xr,yr)  ! Ordinary dot product
        Else
          Allocate(x2(nc));  x2= xr

          x2(lst(:n))= fac*x2(lst(:n))
            x2norm= Sqrt(Sum(x2**2) + eps)
          corr= Dot_product(x2,yr) / x2norm  ! Correlation

          DeAllocate(x2)
        End if
        DeAllocate(lst)
      End if

    End Subroutine Dot_product_MD

!   Modified Dot_prodcut to compute a modified correlation 'corr', 
!   assuming |yr| = 1 and that entries of 'xr' where both 'xr' and 'yr' 
!   are negative will be scaled by 'fac' and renormalized before 
!   computing the correlation. If 0 <= fac < 1 this has the effect of
!   decreasing the correlation, so that it depends more on the
!   correlation of positive entries. 
    
    Subroutine Dot_product_M1R (fac, xr,yr, x1,ls, n,corr)

      Real,     Intent(in) :: fac     ! fac < 0 or fac >= 1:  Ordinary dot product of xr and yr.
                                      ! fac >= 0: First scale xr by 'fac' to get x1, where both xr & yr
                                      !           are negative. Then normalize x1, so that 'corr' becomes
                                      !           a correlation if |yr| = 1
      Real,     Intent(in) :: xr(0:)  ! (0:nc) Assuming |xr(1:)| = 1 for a correlation,
                                      !        with xr(0)= Euclidean norm of xr(1:) prior to normalization
      Real,     Intent(in) :: yr(:)   ! (nc)   Assuming |yr| = 1 for a correlation
      
      Real,    Intent(out) :: x1(0:)  ! (0:nc) 'xr' as modified, with the norm of x1(1:) by scaled by xr(0)
                                      !        the prior norm xr(0) to make x1(0) the modified prior norm.
      Integer, Intent(out) :: ls(:)   ! (nc) List of 'n' coordinates where both xr and yr are negative

      Integer, Intent(out) :: n       ! # coordinates where both xr(1:) and yr are negative, n <= nc,
                                      ! except n = 0 if ordinary dot product
      Real,    Intent(out) :: corr    ! The correlation Sum(x1(1:)*yr) with |x1(1:)| = 1 and |yr| = 1,
!   Local:      
      Real, Parameter :: eps= 1.0E-6
      
      Call List_of_true (xr(1:) < 0 .and. yr < 0, n,ls)
          
      If (n <= 0 .or. fac < 0 .or. fac >= 1) then
        corr= Dot_product(xr(1:),yr)  ! Ordinary dot product
        x1= xr;  ls= -1
      Else
        x1(1:)= xr(1:);  x1(ls(:n))= fac*x1(ls(:n))
        x1(0)= Sqrt(Sum(x1(1:)**2) + eps)
        x1(1:)= x1(1:) / x1(0);  x1(0)= xr(0)*x1(0);  ls(n+1:)= -1
        corr= Dot_product(x1(1:),yr)  ! Correlation
      End if
      
    End Subroutine Dot_product_M1R
    

    Subroutine Dot_product_M1D (fac, xr,yr, x1,ls, n,corr)

      Real(Dblp),  Intent(in) :: fac     ! fac < 0 or fac >= 1:  Ordinary dot product of xr and yr.
                                            ! fac >= 0: First scale xr by 'fac' to get x1, where both xr & yr
                                            !           are negative. Then normalize x1, so that 'corr' becomes
                                            !           a correlation if |yr| = 1
      Real(Dblp),  Intent(in) :: xr(0:)  ! (0:nc) Assuming |xr(1:)| = 1 for a correlation,
                                            !        with xr(0)= Euclidean norm of xr(1:) prior to normalization
      Real(Dblp),  Intent(in) :: yr(:)   ! (nc)   Assuming |yr| = 1 for a correlation
      
      Real(Dblp), Intent(out) :: x1(0:)  ! (0:nc) 'xr' as modified, with the norm of x1(1:) by scaled by xr(0)
                                            !        the prior norm xr(0) to make x1(0) the modified prior norm.
      Integer,       Intent(out) :: ls(:)   ! (nc) List of 'n' coordinates where both xr and yr are negative

      Integer,       Intent(out) :: n       ! # coordinates where both xr(1:) and yr are negative, n <= nc,
                                            ! except n = 0 if ordinary dot product
      Real(Dblp), Intent(out) :: corr    ! The correlation Sum(x1(1:)*yr) with |x1(1:)| = 1 and |yr| = 1,
!   Local:      
      Real(Dblp), Parameter :: eps= 1.0E-13
      
      Call List_of_true (xr(1:) < 0 .and. yr < 0, n,ls)
          
      If (n <= 0 .or. fac < 0 .or. fac >= 1) then
        corr= Dot_product(xr(1:),yr)  ! Ordinary dot product
        x1= xr;  n= 0;  ls= -1
      Else
        x1(1:)= xr(1:);  x1(ls(:n))= fac*x1(ls(:n))
        x1(0)= Sqrt(Sum(x1(1:)**2) + eps)
        x1(1:)= x1(1:) / x1(0);  x1(0)= xr(0)*x1(0);  ls(n+1:)= -1
        corr= Dot_product(x1(1:),yr)  ! Correlation
      End if
    End Subroutine Dot_product_M1D


    Subroutine Percentiles (val, per)

!     Given percentile levels per(:,1) as 0 to 1 fractions
!     compute their corresponding values per(:,2)
!     for the data values 'vl' in increasing order.

      Real,    Intent(in) :: val(:)    ! (nv)  Increasing statistical values

      Real, Intent(inout) :: per(:,:)  ! (nl,2)  Percentile fractions (input) & values (output)
                                       !   (:,1) = increasing fractions, 0 to 1
                                       !           repesenting to percentiles
                                       !   (:,2) = increasing values = the 
                                       !           percentile values
!   Local:
      Real    :: x
      Integer :: i, j, j1, nl, nv

      nl= Size(per,1);  nv= Size(val)
      If (nl < 1 .or. nv < 1) Return

      Do i= 1,nl
        x= per(i,1) * nv;  j= Floor(x);  j1= j + 1

        If (j > 0 .and. j < nv) then
          per(i,2)= val(j) * (j1 - x) + val(j+1) * (x - j)
        Else if (j < 1) then
          per(i,2)= val(1)
        Else if (j1 > nl) then
          per(i,2)= val(nv)
        End if
      End do

    End Subroutine Percentiles

    
    Pure Real Function MedianR (rl)

!     Compute the median of an ordered set of real numbers 'rl'

      Real, Intent(in) :: rl(:) ! (n)
!   Local:
      Integer :: i, n

      n= Size(rl);  i= n/2
    
      If (n <= 0) then
        MedianR= 0
      Else if (n == 2*i) then
        MedianR= (rl(i) + rl(i+1))/2
      Else
        MedianR= rl(i+1)
      End if
    End function MedianR 
  
    Pure Double Precision Function MedianD (rl)

!     Compute the median of an ordered set of 
!     double precision numbers 'rl'

      Real(Dblp), Intent(in) :: rl(:) ! (n)
!   Local:
      Integer :: i, n

      n= Size(rl);  i= n/2
    
      If (n <= 0) then
        MedianD= 0
      Else if (n == 2*i) then
        MedianD= (rl(i) + rl(i+1))/2
      Else
        MedianD= rl(i+1)
      End if
    End function MedianD 
  
    Pure Real Function Median2R (ord, rl)

!     Compute the median of an ordered subset 'ord' of real numbers 'rl'

      Integer, Intent(in) :: ord(:)  ! (l) 1 <= ord(i), ord(i+1) <= n
      Real,    Intent(in) :: rl(:)   ! (n) l <= n
!   Local:
      Integer :: i, l

      l= Size(ord);  i= l/2
    
      If (l <= 0) then
        Median2R= 0
      Else if (l == 2*i) then
        Median2R= (rl(ord(i)) + rl(ord(i+1)))/2
      Else
        Median2R= rl(ord(i+1))
      End if
    End function Median2R 
  
    Pure Double Precision Function Median2D (ord, rl)

!     Compute the median of an ordered subset 'ord' 
!     of double precision numbers 'rl'

      Integer,       Intent(in) :: ord(:)  ! (l) 1 <= ord(i), ord(i+1) <= n
      Real(Dblp), Intent(in) :: rl(:)   ! (n) l <= n
!   Local:
      Integer :: i, l

      l= Size(ord);  i= l/2
    
      If (l <= 0) then
        Median2D= 0
      Else if (l == 2*i) then
        Median2D= (rl(ord(i)) + rl(ord(i+1)))/2
      Else
        Median2D= rl(ord(i+1))
      End if
    End function Median2D 
  
  
! Hypotenuse computation that avoides underflow or overflow    
    
  Elemental Function HypotR(a,b)
    Real, Intent(in) :: a, b
    
    Real :: HypotR  ! Output
! Local:
    Real :: x, y, z
    
    x= Abs(a);  y= Abs(b)
    
    If (x > y) then
      z= y / x;  HypotR= x * Sqrt(1.0 + z*z)
    Else if (y > x) then
      z= x / y;  HypotR= y * Sqrt(1.0 + z*z)
    Else
      HypotR= Sqrt(2.0) * x
    End if
  End Function HypotR

  Elemental Function HypotD(a,b)
    Real(Dblp), Intent(in) :: a, b
    
    Real(Dblp) :: HypotD  ! Output
! Local:
    Real(Dblp) :: x, y, z
    
    x= Abs(a);  y= Abs(b)
    
    If (x > y) then
      z= y / x;  HypotD= x * Sqrt(1.0_Dblp + z*z)
    Else if (y > x) then
      z= x / y;  HypotD= y * Sqrt(1.0_Dblp + z*z)
    Else
      HypotD= Sqrt(2.0_Dblp) * x
    End if
  End Function HypotD
  

  Pure Real Function Lin_interp (xin, delta, xa, ya)

!   Linear interpolation function.   Compute y(x), given input x value xin,
!   x array xa, and y array ya.

!   If delta > 0, assume that the x values are separated by delta, 
!   so only the value of xa(1) is actually used.  Otherwise, do a 
!   bisection search.

!   The x values may be increasing or decreasing.

    Real, Intent(in) :: xin
    Real, Intent(in) :: delta
    Real, Intent(in) :: xa(:)
    Real, Intent(in) :: ya(:)
! Local:
    Logical :: Increasing
    Integer :: j, j1, j2, nn
    Real    :: x1, x2, last_x

    Lin_interp= 0.0;  nn= Size(ya)

    If (delta <= 0) then
      nn= Min(nn, Size(xa))
      last_x= xa(nn);  Increasing= xa(1) <= last_x
    Else
      last_x= xa(1) + delta * (nn - 1)
      Increasing= .True.
    End if

    If (Increasing) then
      If (xin >= last_x)  then
        Lin_interp= ya(nn);  Return
      Else if (xin <= xa(1)) then
        Lin_interp= ya(1);  Return
      End if
    Else
      If (xin <= last_x)  then
        Lin_interp= ya(nn);  Return
      Else if (xin >= xa(1)) then
        Lin_interp= ya(1);  Return
      End if
    End if

!   Linear interpolation

    If (delta > 0.0) then    ! Direct computation
      j1= Ceiling(xin - xa(1) / delta);  j2= j1+1
      x2= xa(1) + delta * j1;  x1= x2 - delta

      Lin_interp= ((xin - x1) * ya(j2) + (x2 - xin) * ya(j1)) / delta
    Else                     ! Bisection search
      If (Increasing) then
        j1= Binary_search (xin, xa)
      Else
        j1= Binary_search (-xin, -xa)
      End if
      j2= j1 + 1

      Lin_interp= ((xin - xa(j1)) * ya(j2) + (xa(j2) - xin) * ya(j1)) / (xa(j2) - xa(j1))
    End if

  End function Lin_interp 
  
  
  Pure Real Function BiLin_interp (xin,yin, delta, xa, ya, za)

!   Bilinear interpolation of the point (xin,yin) on the 
!   xy grid specified by the x points xa and y points ya
!   with values za at the grid points.

!   If delta(1) > 0, then only the first value of xa is used,
!   with the remaining values incremented by delta(1).
!   If delta(2) > 0, then only the first value of ya is used,
!   with the remaining values incremented by delta(2).

!   The xa and ya values may be increasing or decreasing.

    Real, Intent(in) :: xin, yin
    Real, Intent(in) :: delta(2)
    Real, Intent(in) :: xa(:), ya(:)
    Real, Intent(in) :: za(:,:)
! Local:
    Logical :: Increasing(2), Err_type(6)
    Integer :: j, i1, i2, j1, j2, nn, mm
    Real    :: x1, x2, y1, y2, z1, z2, last_x, last_y
    Real    :: dx, d1, d2

    BiLin_interp= 0.0;  nn= Size(za,1);  mm= Size(za,2)

    If (delta(1) <= 0) then
      nn= Min(nn, Size(xa))
      last_x= xa(nn);  Increasing(1)= xa(1) <= last_x
    Else
      last_x= xa(1) + delta(1) * (nn - 1)
      Increasing(1)= .True.
    End if

    If (delta(2) <= 0) then
      mm= Min(mm, Size(ya))
      last_y= ya(mm);  Increasing(2)= ya(1) <= last_y
    Else
      last_y= ya(1) + delta(2) * (mm - 1)
      Increasing(2)= .True.
    End if

    If (Increasing(1)) then
      If (xin >= last_x)  then
        BiLin_interp= Lin_interp (yin, delta(2), ya, za(nn,:));  Return
      Else if (xin <= xa(1)) then
        BiLin_interp= Lin_interp (yin, delta(2), ya, za(1,:));  Return
      End if
    Else
      If (xin <= last_x)  then
        BiLin_interp= Lin_interp (yin, delta(2), ya, za(nn,:));  Return
      Else if (xin >= xa(1)) then
        BiLin_interp= Lin_interp (yin, delta(2), ya, za(1,:));  Return
      End if
    End if

    If (Increasing(2)) then
      If (yin >= last_y)  then
        BiLin_interp= Lin_interp (xin, delta(1), xa, za(:,mm));  Return
      Else if (yin <= ya(1)) then
        BiLin_interp= Lin_interp (xin, delta(1), xa, za(:,1));  Return
      End if
    Else
      If (yin <= last_y)  then
        BiLin_interp= Lin_interp (xin, delta(1), xa, za(:,mm));  Return
      Else if (yin >= ya(1)) then
        BiLin_interp= Lin_interp (xin, delta(1), xa, za(:,1));  Return
      End if
    End if

!   x interpolation

    If (delta(1) > 0.0) then    ! Direct computation
      i1= Ceiling(xin - xa(1) / delta(1));  i2= i1+1
      x2= xa(1) + delta(1) * i1;  x1= x2 - delta(1)
      dx= delta(1);  d1= xin - x1;  d2= x2 - xin
    Else                        ! Bisection search 
      If (Increasing(1)) then
        i1= Binary_search (xin, xa)
      Else
        i1= Binary_search (-xin, -xa)
      End if
      i2= i1 + 1;  dx= xa(i2) - xa(i1)
      d1= xin - xa(i1);  d2= xa(i2) - xin
    End if

!   y interpolation

    If (delta(2) > 0.0) then   ! Direct computation
      j1= Ceiling(yin - ya(1) / delta(2));  j2= j1+1
      y2= ya(1) + delta(2) * j1;  y1= y2 - delta(2)

      z1= (d1 * za(i2,j1) + d2 * za(i1,j1)) / dx
      z2= (d1 * za(i2,j2) + d2 * za(i1,j2)) / dx
      BiLin_interp= ((yin - y1) * z2 + (y2 - yin) * z1) / delta(2)
    Else                       ! Bisection search
      If (Increasing(2)) then
        j1= Binary_search (yin, ya)
      Else
        j1= Binary_search (-yin, -ya)
      End if
      j2= j1 + 1

      z1= (d1 * za(i2,j1) + d2 * za(i1,j1)) / dx
      z2= (d1 * za(i2,j2) + d2 * za(i1,j2)) / dx
      BiLin_interp= ((yin - ya(j1)) * z2 + (ya(j2) - yin) * z1) / (ya(j2) - ya(j1))
    End if

  End function BiLin_interp 

 
   Elemental Real Function Interp_angles (ang1,ang2, del,del1)  ! Linear interpolation of angles in radians
     Real, Intent(in) :: ang1,ang2
     Real, Intent(in) :: del,del1

     If (Abs(ang2-ang1) <= PI) then
       Interp_angles= del * ang2 + del1 * ang1
     Else if (ang2 > PI) then
       Interp_angles= del * (ang2-Two_PI) + del1 * ang1
     Else
       Interp_angles= del * ang2 + del1 * (ang1-Two_PI)
     End if
     Interp_angles= Modulo(Interp_angles, Two_PI)
   End Function Interp_angles

   
   Pure Subroutine Inverse_map (in_map, inverse)
   
!    Compute the 'inverse' map to 'in_map'.
!    Note: If 'in_map' is actually a many to 1 map, then
!    the inverse 'i' of 'j' will be the first value of 'i' 
!    for which in_map(i) = j.
   
     Integer,  Intent(in) :: in_map(:)   ! (ni) The map to invert, for 1 <= in_map <= mi
     Integer, Intent(out) :: inverse(:)  ! (mi) The inverse of in_map, = 0 outside 
                                         !      the range of in_map or where in_map < 1 or > mi
!  Local:     
     Integer :: i, j, ni, mi

     ni= Size(in_map);  mi= Size(inverse)
     
     inverse= 0
     Do i= ni,1,-1
       j= in_map(i);  If (j >= 1 .and. j <= mi) inverse(j)= i
     End do
   End Subroutine Inverse_map

   
   Pure Function Reverse (in_map)
   
     Integer, Intent(in) :: in_map(:)  ! (n) 
!  Result:
     Integer :: Reverse(Size(in_map))  ! (n) 'in_map' in reverse order
!  Local:     
     Integer :: i, n, n1

     n= Size(in_map);  n1= n + 1
     Forall(i=1:n) Reverse(i)= in_map(n1-i)
     
   End Function Reverse

   
   Elemental Real Function Normal_tail (x)

!    Compute the tail of the cumulative normal distribution = 
!    1 - cumulative normal distribution function.

!    This uses a rational approximation formula that is accurate
!    to less than 1.0 E-5.  Ref:  Handbook of Mathematical Formulas,
!    Tables, Functions, Graphs, Transforms, 
!    Research and Education Association, New York, 1980

     Real, Intent(in) :: x
!  Local:
     Real :: t
     t= 1.0 / (1.0 + .33267 * x) 
     Normal_tail= (t / Sqrt(Two_PI)) * Exp(-0.5 * x**2) * &
                  (.4361836 - .1201676 * (t + 0.9372980 * t))
   End Function Normal_tail
  

   Pure Integer Recursive Function LCD (n,m)

!    Compute the least common denominator LCD of two integers n < m.
 
     Integer, Intent(in) :: n, m
!  Local:
     Integer :: q, rr
     q = m / n
     rr= m - q * n

     If (rr > 1) then
        LCD= LCD (rr,n)  
     Else if (rr == 1) then
       LCD= 1
     Else
       LCD= n
     End if
   End Function LCD


   Pure Integer Function Index_in_list (xx, list)

!    Search for the index in 'list' which yields the "nearest" 
!    value to 'xx'.  We assume that 'list' is roughly linear.

!    We also assume that 'list' represents the mid-points of intervals,
!    and we are actually trying to find the nearest interval start point
!    to which 'xx' belongs.

     Real, Intent(in) :: xx
     Real, Intent(in) :: list(:)
!  Local:
     Integer :: i, j, mx

     mx= Size(list)

     If (list(1) <= list(mx)) then  ! List is increasing
       If (xx <= list(1)) then
         Index_in_list= 1
       Else if (xx > list(mx)) then
         Index_in_list= mx + 1
       Else
         j= Max(Nint(mx * (xx - list(1)) / (list(mx) - list(1))), 1) ! Initial guess
  
         If (xx <= list(j)) then  ! Search down
           Do i= j-1,1,-1
             If (xx > list(i))  Exit
           End do
           Index_in_list= i + 1
         Else
           Do i= j+1,mx           ! Search up
             If (xx <= list(i))  Exit
           End do
           Index_in_list= i
         End if
       End if
     Else  ! List is decreasing
       If (xx >= list(1)) then
         Index_in_list= 1
       Else if (xx < list(mx)) then
         Index_in_list= mx + 1
       Else
         j= Max(Nint(mx * (list(1) - xx) / (list(1) - list(mx))), 1) ! Initial guess
  
         If (xx >= list(j)) then  ! Search down
           Do i= j-1,1,-1
             If (xx < list(i))  Exit
           End do
           Index_in_list= i + 1
         Else
           Do i= j+1,mx           ! Search up
             If (xx >= list(i))  Exit
           End do
           Index_in_list= i
         End if
       End if
     End if

     Index_in_list= Min(Index_in_list, mx)
   End Function Index_in_list

   
   Elemental Real Function To_dB (x)
     Real, Intent(in) :: x
     To_dB= 10.0 * Log10(x)   ! In dB
   End Function To_dB

   
   Elemental Real Function From_dB (x)
     Real, Intent(in) :: x   ! In dB
     From_dB= 10.0**(x / 10.0)
   End Function From_dB


   Pure Subroutine Subset_ofI (set,Mask, n,sub)

!    Compute the subset 'sub' of 'set' where Mask is true.

     Integer,    Intent(in) :: set(:)    ! (m) The set
     Logical,    Intent(in) :: Mask(:)   ! (m) The true elements of the set
     Integer,   Intent(out) :: n         ! Size of the subset
     Integer, Intent(inout) :: sub(:)    ! (m to n) The subset
     
     n= Count(Mask);  sub(:n)= Pack(set, Mask);  sub(n+1:)= 0
     
   End Subroutine Subset_ofI
   
   Pure Subroutine Subset_ofR (set,Mask, n,sub)

!    Compute the subset 'sub' of 'set' where Mask is true.

     Real,       Intent(in) :: set(:)    ! (m) The set
     Logical,    Intent(in) :: Mask(:)   ! (m) The true elements of the set
     Integer,   Intent(out) :: n         ! Size of the subset
     Real,    Intent(inout) :: sub(:)    ! (m to n) The subset
     
     n= Count(Mask);  sub(:n)= Pack(set, Mask);  sub(n+1:)= 0
     
   End Subroutine Subset_ofR
   

   Pure Subroutine List_of_true0 (Mask, n,ls)

!    Compute the list of indices where Mask is true, 

     Logical,  Intent(in) :: Mask(:) ! (m)
     Integer, Intent(out) :: n       ! # true indices found out of 'm'
     Integer, Intent(out) :: ls(:)   ! (>= n) List of the indices where Mask is true
! Local
     Integer :: lid(Size(Mask))
     
     n= Count(Mask);  lid= "ID";  ls(:n)= Pack(lid,Mask)
     
   End Subroutine List_of_true0


   Pure Subroutine List_of_true1 (k,Mask, n,ls)

!    Compute the list of indices where 'Mask' is true, 
!    assuming that 'Mask' starts at 'k'

     Integer,  Intent(in) :: k         ! First index of Mask
     Logical,  Intent(in) :: Mask(k:)  ! (k:l) where l = k-1 + m and m = Size(Mask)

     Integer, Intent(out) :: n         ! # true indices found out of 'm'
     Integer, Intent(out) :: ls(:)     ! (>= n) List of the indices where Mask is true
                                       !        starting from 'k'
! Local
     Integer :: i, j, m, l, p
     
     m= Size(Mask);  p= Size(ls);  l= (k-1) + m

     n= 0
     Do i= k,l
       If (Mask(i)) then
         n= n + 1;  ls(n)= i;  If (n >= p) Exit
       End if
     End do

     If (n < p) ls(n+1:)= 0
     
   End Subroutine List_of_true1

   
   Pure Subroutine Array_int (x, n,y)
!    Array assignment over a common length
     Integer,    Intent(in) :: x(:)
     Integer,   Intent(out) :: n
     Integer, Intent(inout) :: y(:)

     n= Min(Size(x), Size(y))
     If (n > 0) y(:n)= x(:n)
   End Subroutine Array_int   

   
   Pure Subroutine Array_real (x, n,y)
!    Array assignment over a common length
     Real,     Intent(in) :: x(:)
     Integer, Intent(out) :: n
     Real,  Intent(inout) :: y(:)

     n= Min(Size(x), Size(y))
     If (n > 0) y(:n)= x(:n)
   End Subroutine Array_real   

   
   Pure Subroutine Pack_pos (list_in, n,list_out)

!    Pack the "n" positive elements of list_in into list_out.

     Integer,  Intent(in) :: list_in(:)
     Integer, Intent(out) :: n           ! # list_in values packed into list_out
     Integer, Intent(out) :: list_out(:) ! = packed positive list_in values
!  Local:
     Call Array (Pack(list_in, list_in > 0), n,list_out)
   End Subroutine Pack_pos


   Elemental Real Function Horizon (plat_alt, targ_el)
     Real, Intent(in) :: plat_alt, targ_el
!  Local:
     Real(Dblp) :: b, c, Ea, disc

     Ea= earth_rad43 + plat_alt
     b = Ea * Sin(targ_el)
     c = plat_alt * (earth_rad43 + Ea)
     disc= b**2 - c

     If (disc > 0.0 .and. targ_el < 0.0) then
       Horizon= -b - Sqrt(disc)
     Else
       Horizon= Huge(plat_alt)
     End if
   End Function Horizon

  
  Pure Subroutine Comp_ellipse (cov, ellip)

!   This routine computes an error ellipse from a covariance matrix.

    Real,  Intent(in) :: cov(:,:) ! (>=2,>=2)
    Real, Intent(out) :: ellip(:) ! (>=3) 1 = half length of longest axis 
                                  ! 2 = half length of shortest axis
                                  ! 3 = angle of longest axis wrt the x axis (radians)
                                  ! 4 = half length of z axis if applicable
!   Use trace = ellip(1) + ellip(2) = cov(1,1) + cov(2,2)
!   Use determinant = ellip(1) * ellip(2) = cov(1,1) * cov(2,2) - cov(1,2)**2

! Local:
    Real(Dblp) :: tr, det, dif, lm(2)

    ellip= 0.0;  tr= cov(1,1) + cov(2,2)

    If (Abs(cov(1,2)) > 0.0) then
      det= cov(1,1) * cov(2,2) - cov(1,2)**2
      dif= Sqrt(Abs(tr**2 - 4.0 * det))
      lm(1)= (tr + dif) / 2.0
      lm(2)= Abs(tr - dif) / 2.0
      ellip(3)  = Atan2(lm(1) - cov(1,1), cov(1,2))
      ellip(1:2)= Sqrt(lm)
    Else
      If (cov(1,1) >= cov(2,2)) then
        ellip(1)= Sqrt(cov(1,1));  ellip(2)= Sqrt(cov(2,2))
      Else
        ellip(1)= Sqrt(cov(2,2));  ellip(2)= Sqrt(cov(1,1))
      End if
    End if

    If (Size(cov,1) > 2)  ellip(4)= Sqrt(cov(3,3))
  End Subroutine Comp_ellipse


  Pure Subroutine Comp_covar (ellip, cov)

!   This routine computes a covariance matrix from an error ellipse.

    Real,  Intent(in) :: ellip(:) ! (>=3) 1 = half length of longest axis 
                                  ! 2 = half length of shortest axis
                                  ! 3 = angle of longest axis wrt the x axis (radians)
                                  ! 4 = half length of z axis if applicable
    Real, Intent(out) :: cov(:,:) ! (>=2,>=2)

! Local:
    Real :: e1, e2, e1_sq, e2_sq, var_p1, var_p2

    cov= 0.0
    e1= Cos(ellip(3));  e1_sq= e1**2
    e2= Sin(ellip(3));  e2_sq= e2**2
    var_p1= ellip(1)**2;  var_p2= ellip(2)**2
    cov(1,1)= e1_sq * var_p1 + e2_sq * var_p2
    cov(2,2)= e1_sq * var_p2 + e2_sq * var_p1
    cov(1,2)= (var_p1 - var_p2) * e1 * e2;  cov(2,1)= cov(1,2)

    If (Size(cov,1) > 2)  cov(3,3) = ellip(4)**2
  End Subroutine Comp_covar


   Pure Subroutine Project_to_ellipse (n_sig, xx, cov, xp, x_sig)

!    Project a point out side an ellipse to the nearest point in the ellipse, 
!    approximately.

     Real,  Intent(in) :: n_sig    ! Ellipse size = # sigmas 
     Real,  Intent(in) :: xx(:)    ! (2) Point to be projected, wrt ellipse center
     Real,  Intent(in) :: cov(:,:) ! (>=2,>=2) Covariance matrix that determines the ellipse
     Real, Intent(out) :: xp(:)    ! (2) Projected point, wrt to ellipse center
     Real, Intent(out) :: x_sig    ! Length of xx in ellipse sigmas.
!  Local:
     Real :: x0, xe, ye, long_dir(2), short_dir(2), ellip(3)
     Real :: a, b, c, long, short, sl, tmp, x1, y1

     Call Comp_ellipse (cov(:2,:2), ellip)

     long_dir(1)= Cos(ellip(3)); long_dir(2)= Sin(ellip(3))
     xe= Dot_product (xx, long_dir)   ! Long axis coordinate
     short_dir(1)= -long_dir(2);  short_dir(2)= long_dir(1)
     ye= Dot_product (xx, short_dir)  ! Short axis coordinate
     
     x_sig= Sqrt((xe/ellip(1))**2 + (ye/ellip(2))**2)
     If (x_sig <= n_sig) then  ! xx is inside the ellipse
       xp= xx;  Return       
     End if

     long= n_sig * ellip(1);  short= n_sig * ellip(2)
     x0= 0.75 * Sign(Min(Abs(xe), long), xe)
     sl= ye / (xe - x0);  tmp= (sl / short)**2
     a= (1.0/long)**2 + tmp;  tmp= tmp * x0
     b= 2.0*tmp;  c= tmp*x0 - 1.0

     x1= (b + Sign(Sqrt(Abs(b**2 - 4.0*a*c)), xe)) / (2.0*a)  ! Quadratic formula: Solve (x1/long)**2 + (y1/short)**2= 1
                                                              ! where y1 = sl * (x1 - x0)
     y1= sl * (x1 - x0)
     xp= x1 * long_dir + y1 * short_dir
   End subroutine Project_to_ellipse

   
   Subroutine Ran_G (nn,seed, Gaus)

!    seed<= 0 : Continue (don't re-set the seed for "Random_number").
!         > 0 : Re-set the seed for "Random_number" according to the 
!               procedure below.  "seed" may be a small integer.
     Integer,         Intent(in) :: nn
     Integer,         Intent(in) :: seed
     Real, Optional, Intent(out) :: Gaus(nn)
!   Local:
     Integer :: i, n1, sd(2)
     Real    :: ut(2), uniform(nn+1), fac, tmp, eps

     If (seed > 0) then
       Call Random_number (ut)
       sd= Nint(100000 * seed * ut)
       Call Random_seed (put=sd)
     End if
     If (.not.Present(Gaus))  Return
         
     n1= 2 * ((nn+1)/2);  eps= 1.0E-6
     Call Random_number (uniform(:n1))

     Do i= 1,nn-1,2
       fac= Sqrt(-2.0 * Log(Max(uniform(i), eps)))
       tmp= Two_PI * uniform(i+1)
       Gaus(i)  = fac * Cos(tmp)
       Gaus(i+1)= fac * Sin(tmp)
     End do

     If (n1 > nn) then
       fac= Sqrt(-2.0 * Log(Max(uniform(i), eps)))
       Gaus(i)= fac * Cos(Two_PI * uniform(i+1))
     End if
   End subroutine Ran_G

   
   Function Ran_Gaus (nn)

!    Use the law of large numbers approximation.
!    Note that a uniform [0,1] random number has variance 1/12,
!    so the sum of 12 of these has variance 1.0 with mean 6.0.

     Integer, Intent(in) :: nn
!  Result:
     Real                :: Ran_Gaus(nn)
!  Local:     
     Integer :: i
     Real    :: uniform(12)

     Do i= 1,nn
       Call Random_number (uniform)
       Ran_Gaus(i)= Sum(uniform) - 6.0
     End do
   End Function Ran_Gaus


     Elemental Function Quad_fit (t0,t1, f0,f1, df0,df1)

  !    Fit a quadratic function to the data (t,f(t),df(t)) for t= t0,t1.
  !    Then solve for the minimum value in the interval (t0,t1).

       Real(Dblp), Intent(in) :: t0, t1, f0, f1, df0, df1

       Real(Dblp) :: Quad_fit
       Real(Dblp) :: a, b, d, e, dI, d0, d1

       d= t1 - t0

       If ((f0 <= f1 .and. df0 > 0) .or. Abs(d) < 0.001) then  ! Left end
         Quad_fit= t0
       Else if (f1 <= f0 .and. df1 < 0) then  ! Right end
         Quad_fit= t1

       Else
         d0= Min(df0, 0.0_Dblp);  d1= Max(df1, 0.0_Dblp)

         If (d0 >= 0 .and. d1 <= 0) then  ! Linear interpolation
           Quad_fit= t0 + f0*d/(f0+f1)

         Else  !  Do a quadratic fit, where
               !        f(t)  = a*t^2 + b*t + f0
               !        df(t) = 2*a*t + b
           dI= 1/d;  e= (f1-f0)*dI

           If (d1 <= 0 .or. f0 <= f1) then  ! Ignore d1
             b= d0;  a= (e-b)*dI
           Else                             ! Ignore d0
             a= (d1-e)*dI
             b= e - d*a
           End if

           Quad_fit= t0 - b/(2*a)
           Quad_fit= Max(t0, Min(Quad_fit, t1))
         End if
       End if
     End function Quad_fit

   Subroutine Quad_search (Func, tol,change,t0, tl,tr, nfn,tf)

!    Do quadratic fits to minimize the function f(t) given 
!    a local minimum at t0 between tl and tr.

     External                         Func
     Real(Dblp),     Intent(in) :: tol, change
     Real(Dblp),     Intent(in) :: t0      ! Initial value for t
     Real(Dblp),  Intent(inout) :: tl, tr  ! Lower and upper bounds for t

!    tol    Stopping tolerance for change in the variable t
!    change Stopping tolerance for fraction change in f(t) or absolute value of f'(t)

!    Func   Function evaluation subroutine.  Usage:

!    Call Func (Both,t, fn,fd)

!      Both : Compute both 'fn' and 'fd'.
!      t    : Variable at which f(t)= 'fn' and f(t) = 'fd' are evaluated.
!      fn   : Value of f(t);  A value of Huge(real) indicates infeasibility.
!      fd   : Value of the derivative f'(t);  A value of Huge(real) indicates infeasibility.

     Integer,       Intent(out) :: nfn
     Real(Dblp), Intent(out) :: tf(:,0:)  ! (3,0:)

!    tf(1:3,j) = t, f(t), f'(t), values in evaluation order.

!  Local:
     Logical       :: Both= .true.
     Integer       :: i, i1, i2, lim
     Real          :: rl
     Real(Dblp) :: fn, fd, t_old, t_new, fn_chng, fn_old
     
     rl= Huge(rl) - 1.0
     tl= Min(tl,0.999*t0);  tr= Max(tr,1.001*t0)

     Call Func (Both,t0, fn,fd)

     nfn= 0;  tf(1,0)= t0;  tf(2,0)= fn;  tf(3,0)= fd
     If (Abs(fn) >= rl) Return

     If (fd > 0.0) then
       lim= Ceiling(Log((t0-tl)/tol) / Log(2.0))
       Do i= 1,lim
         Call Func (Both,tl, fn,fd)

         If (Abs(fn) < rl) Exit

         tl= (tl + t0) / 2.0
       End do

       If (i <= lim) then
         t_new= tl;  i1= 1;  i2= 0
       Else
         Call Func (Both,t0, fn,fd);  Return
       End if
     Else
       lim= Ceiling(Log((t0-tl)/tol) / Log(2.0))
       Do i= 1,lim
         Call Func (Both,tr, fn,fd)

         If (Abs(fn) < rl) Exit

         tr= (tr + t0) / 2.0
       End do

       If (i <= lim) then
         t_new= tr;  i1= 0;  i2= 1
       Else
         Call Func (Both,t0, fn,fd);  Return
       End if
     End if

     tf(1,1)= t_new;  tf(2,1)= fn;  tf(3,1)= fd
     nfn= 1;  lim= Size(tf,2) - 1

     Do while (nfn < lim)
       t_old= t_new
       t_new= Quad_fit (tf(1,i1),tf(1,i2), tf(2,i1),tf(2,i2), tf(3,i1),tf(3,i2))

       If (Abs(t_new-t_old) <= tol) Return

       nfn= nfn + 1;  fn_old= fn
       Call Func (Both,t_new, fn,fd)

       If (Abs(fn) >= rl) then  ! End at lowest value to date
         i= Minloc(tf(2,0:nfn-1), 1) - 1;  t_new= tf(1,i)
         Call Func (Both,t_new, fn,fd)
         tf(1,nfn)= t_new;  tf(2,nfn)= fn;  tf(3,nfn)= fd
         Return
       End if

       tf(1,nfn)= t_new;  tf(2,nfn)= fn;  tf(3,nfn)= fd
       fn_chng= change * (Abs(fn_old) + 0.001)

       If (Abs(fn-fn_old) <= fn_chng .or. Abs(fd) <= fn_chng)  Return

       If (fd >= 0.0) then
         If (tf(2,i1) <= tf(2,i2) .or. fn < tf(2,i2)) then
           i2= nfn
         Else
           i1= nfn
         End if
       Else
         If (tf(2,i1) >= tf(2,i2) .or. fn < tf(2,i1)) then
           i1= nfn
         Else
           i2= nfn
         End if
       End if
     End do
   End Subroutine Quad_search

   
   Elemental Real Function Normal (xx)

!    Compute the probability that a (0,1) normal random variable exceeds xx.

     Real, Intent(in) :: xx
!  Local:
     Integer       :: i
     Real(Dblp) :: s, t, parm(0:5)

     parm= (/0.2316419, 0.319381530, -0.356563782, 1.781477937, -1.821255978, 1.330274429/)
     t= 1.0 / (1.0 + parm(0)*xx);  s= parm(5)

     Do i= 4,1,-1  ! Horner's rule
       s= parm(i) + s * t
     End do

     Normal= t * s * Exp(-0.5 * xx**2) / Sqrt(Two_PI)
   End Function Normal

   
   Elemental Real Function Dif_angleR (ang)  ! -PI <= Result <= PI
     Real, Intent(in) :: ang

     Dif_angleR= Modulo(ang, Two_PI)
     If (Dif_angleR > PI)  Dif_angleR= Dif_angleR - Two_PI
   End Function Dif_angleR


   Elemental Function Dif_angleD (ang)  ! -PI <= Result <= PI
     Real(Dblp), Intent(in) :: ang
     Real(Dblp) :: Dif_angleD

     Dif_angleD= Modulo(ang, Two_PI)
     If (Dif_angleD > PI)  Dif_angleD= Dif_angleD - Two_PI
   End Function Dif_angleD

   Pure Integer Function First_true_in (Msk)  
     Logical, Intent(in) :: Msk(:)
!  Local:
     Integer, Allocatable :: ara(:)  
     Integer :: k, n  

     n= Size(Msk)
     If (n < 1 .or. All(Msk == .false.)) then
       First_true_in= 0
     Else if (n == 1) then
       First_true_in= 1
     Else
       Allocate(ara(n));  ara= "ID"
       First_true_in= Minloc(ara,1,Msk)
     End if

!     Do First_true_in= 1,n
!       If (Msk(First_true_in))  Exit
!     End do
!     If (First_true_in > n)  First_true_in= 0
   End Function First_true_in

   Pure Integer Function First_true_from (i1, Msk)  
     Integer, Intent(in) :: i1
     Logical, Intent(in) :: Msk(:)
!  Local:
     Integer :: n;  n= Size(Msk)

     Do First_true_from= i1,n
       If (Msk(First_true_from))  Exit
     End do
   End Function First_true_from


   Pure Integer Function Last_true_in (Msk)  
     Logical, Intent(in) :: Msk(:)
!  Local:
     Integer, Allocatable :: ara(:)  
     Integer :: k, n  

     n= Size(Msk)
     If (n < 1 .or. All(Msk == .false.)) then
       Last_true_in= 0
     Else if (n == 1) then
       Last_true_in= 1
     Else
       Allocate(ara(n));  ara= "ID"
       Last_true_in= Maxloc(ara,1,Msk)
     End if
   End Function Last_true_in


   Pure Integer Function Last_true_to (i1,lf)  
     Integer, Intent(in) :: i1
     Logical, Intent(in) :: lf(:)

     Do Last_true_to= Size(lf),i1,-1
       If (lf(Last_true_to))  Exit
     End do
   End Function Last_true_to


   Pure Integer Function Binary_search_int (t, list)
!
!    Assume that "list" is non-decreasing sequence of integers.
!    Then search for "t" in "list", outputting
!    the index of "t" if found, 0 if not found.

!    Approx op count = 3 * Log2(n)

     Integer, Intent(in) :: t
     Integer, Intent(in) :: list(:)
!  Local:
     Integer :: l, r, n, mid

     Binary_search_int= 0;  n= Size(list)

     If (n <= 0) then
       Return
     Else if (t < list(1) .or. t > list(n)) then
       Return
     Else if (t == list(n)) then
       Binary_search_int= n
       Return
     End if

!    Locate the interval [l,r) containing t.

     l= 1;  r= n
     Do
       mid= (l+r) / 2
       If (mid == l)  Exit

       If (t < list(mid)) then
         r= mid
       Else
         l= mid
       End if
     End do

     If (list(l) == t) Binary_search_int= l
   End Function Binary_search_int

   
   Pure Integer Function Binary_search_real (t, list)
!
!    Assume that "list" is a non-decreasing sequence of reals.
!    Then search for the interval that "t" belong to in "list", 
!    outputting index "l" if list(l) <= t < list(l+1), 
!    1 if t <= list(1), n if t >= list(n), for n = Length(list).

     Real, Intent(in) :: t
     Real, Intent(in) :: list(:)
!  Local:
     Integer :: l, r, n, mid

     n= Size(list)

     If (t <= list(1) ) then
       Binary_search_real= 1
       Return
     Else if (t >= list(n)) then
       Binary_search_real= n
       Return
     End if

     l= 1;  r= n
     Do
       mid= (l+r) / 2
       If (mid == l)  Exit

       If (t < list(mid)) then
         r= mid
       Else
         l= mid
       End if
     End do

     Binary_search_real= l
   End Function Binary_search_real


   Pure Integer Function Binary_search_cut (t, cut, list)
!
!    Assume that "list" is non-decreasing sequence of integers > cut,
!    except for occassional integers <= cut.
!    Then search for "t" in the integers > cut, outputting
!    the index of "t" if found, 0 if not found.
     Integer, Intent(in) :: t
     Integer, Intent(in) :: cut
     Integer, Intent(in) :: list(:)
!  Local:
     Integer :: i, np, map(Size(list))

     Binary_search_cut= 0;  np= 0

     Do i= 1,Size(list)
       If (list(i) > cut) then
         np= np + 1;  map(np)= i
       End if
     End do

     If (np > 0) then
       i= Binary_search_int (t, list(map(1:np)))
       If (i > 0)  Binary_search_cut= map(i)
     End if
   End Function Binary_search_cut


   Pure Subroutine Identity_vector_int (ID_vec, ID)
     Integer,     Intent(out) :: ID_vec(:)
     Character(2), Intent(in) :: ID
     Integer :: i
     Forall(i=1:Size(ID_vec)) ID_vec(i)= i
   End Subroutine Identity_vector_int
     

   Pure Subroutine Identity_vector_real (ID_vec, ID)
     Real,        Intent(out) :: ID_vec(:)
     Character(2), Intent(in) :: ID
     Integer :: i
     Forall(i=1:Size(ID_vec)) ID_vec(i)= i
   End Subroutine Identity_vector_real
   
     
    Pure Subroutine Quick_Sort_stable (Ascending, ix, key,ReOrd)

!     Apply a stable algorithm to integer data "ix"
!     by applying the real quick sort algorith

      Logical,            Intent(in) :: Ascending
      Integer,         Intent(inout) :: ix(:)  ! (nn)
      Integer, Optional, Intent(out) :: key(:) ! (nn)
      Logical, Optional, Intent(out) :: ReOrd  ! .true if 'ix' is re-ordered, otherwise
                                               !  the original ordering is correct
!   Local:
      Real, Parameter :: eps_fac= 0.000001
      Integer, Allocatable :: key1(:)
      Real    :: eps, xx(Size(ix))
      Integer :: i, nn

      nn= Size(ix);  If (nn <= 0) Return

      If (Ascending) then
        Do i= 1,nn-1
          If (ix(i) > ix(i+1)) Exit
        End do
      Else
        Do i= 1,nn-1
          If (ix(i) < ix(i+1)) Exit
        End do
      End if
      
      If (Present(ReOrd)) ReOrd= i < nn

      If (Present(key)) then
        If (Size(key) < nn) then  ! Error condition
          key= -1;  Return
        End if
        key(:nn)= 'ID';  If (nn == 1 .or. i >= nn)  Return

        eps= eps_fac * Maxval(ix);  xx= ix
        Call Quick_Sort_real (Ascending, xx, key(:nn), ReOrd, ez=eps)
        ix= ix(key(:nn))

      Else
        If (nn == 1 .or. i >= nn)  Return

        Allocate(key1(nn));  eps= eps_fac * Maxval(ix);  xx= ix
        Call Quick_Sort_real (Ascending, xx, key1, ReOrd, ez=eps)
        ix= ix(key1)
      End if

    End Subroutine Quick_Sort_stable


    Pure Subroutine Quick_Sort_real (Ascending, xx, key, ReOrd, ez)

!     Apply a quick sort algorithm to real data 'xx'.   

!     Optional arguments:    
    
!     Compute the 'key' for further sorting.  That is, to sort yy 
!     in the same way as xx, let yy= yy(key).
    
!     'ReOrd' is true if 'xx' needs to be reordered. Otherwise it is 
!     already in correct order.
    
!     'ez' is an epsilon value for determining the ordering of successive 
!     elements of 'xx' (after reordering) which are almost identical 
!     (within 'ez' of each other). For such elements the final ordering
!     will be the same as the original ordering (this 'stable' property
!     is unlikely with quick sort applied to equal elements).

      Logical,            Intent(in) :: Ascending  ! Sort into increasing order if true, else decreasing
      Real,            Intent(inout) :: xx(:)      ! (nn)  The data to be sorted
      Integer, Optional, Intent(out) :: key(:)     ! (nn)  The sorting key
      Logical, Optional, Intent(out) :: ReOrd      ! True if 'xx' is re-ordered otherwise the 
                                                   !   original ordering is correct
      Real,    Optional,  Intent(in) :: ez         ! Use this value to retain the original order for 
                                                   !   elements of 'xx' which are almost equal 
                                                   !   (within 'ez', successively). Requires 'key'.
!   Local:
      Integer, Parameter :: n_max= 10, stack_max= 200
      Integer :: stack, low(stack_max), high(stack_max)
      Integer :: i, j, k, l, o, r, it, iv, kt, l2, r1, nn
      Real    :: xv, xt

      nn= Size(xx);  If (Present(ReOrd)) ReOrd= .false.
      If (nn <= 0)  Return
      
      If (Present(key)) then
        If (Size(key) < nn) then  ! Error condition
          key= -1;  Return
        End if
        key(:nn)= 'ID'
      End if
      If (nn == 1)  Return
      
      If (Ascending) then
        If (Present(ez)) then
          Do i= 1,nn-1
            If (xx(i) > xx(i+1) - ez) Exit
          End do
        Else
          Do i= 1,nn-1
            If (xx(i) > xx(i+1)) Exit
          End do
        End if
      Else
        If (Present(ez)) then
          Do i= 1,nn-1
            If (xx(i) < xx(i+1) + ez) Exit
          End do
        Else
          Do i= 1,nn-1
            If (xx(i) < xx(i+1)) Exit
          End do
        End if
      End if
      
      If (Present(ReOrd)) ReOrd= i < nn
      If (i >= nn) Return

!     Sort as if ascending
      
      If (.not.Ascending) xx(:nn)= xx(nn:1:-1)
      
      If (Present(key) .and. Present(ez)) then
        Forall(o=1:nn) xx(o)= xx(o) + (o-1)*ez
      End if  

      stack= 0;  l= 1;  r= nn

        Do
          If (r-l < n_max) then ! Insertion sort for the interval [l=left, r=right]
            If (Present(key)) then
              Do j= l+1,r
                If (xx(j-1) <= xx(j))  Cycle
                xv= xx(j);  kt= key(j)
                Do i= j-1,l,-1
                  If (xx(i) <= xv)  Exit
                  xx(i+1)= xx(i);  key(i+1)= key(i)
                End do
                xx(i+1)= xv;  key(i+1)= kt
              End do

            Else   ! No key
              Do j= l+1,r
                If (xx(j-1) <= xx(j))  Cycle
                xv= xx(j)
                Do i= j-1,l,-1
                  If (xx(i) <= xv)  Exit
                  xx(i+1)= xx(i)
                End do
                xx(i+1)= xv
              End do
            End if

            If (stack <= 0) Exit

            r= high(stack);  l= low(stack)
            stack= stack - 1

          Else   ! Quick sort partitioning step for the interval [l=left, r=right]

            k= (l+r)/2;  i= l;  j= r

            If (Present(key)) then

!             Sort left(l), middle(k), and right(r) entries

              If (xx(l) > xx(r)) then
                xt= xx(l);  kt= key(l)
                xx(l)= xx(r);  key(l)= key(r)
                xx(r)= xt;  key(r)= kt
              End if

              If (xx(k) > xx(r)) then
                xt= xx(k);  kt= key(k)
                xx(k) = xx(r);  key(k)= key(r)
                xx(r) = xt;  key(r)= kt
              End if

              If (xx(l) > xx(k)) then
                xt= xx(l);  kt= key(l)
                xx(l)= xx(k);  key(l)= key(k)
                xx(k)= xt;  key(k)= kt
              End if

!             Swap the entries [l+1, r-1] as necessary so that if
!             xx(i) <= xx(k) < xx(j) or xx(i) < xx(k) <= xx(j)  then i < j.

              xv= xx(k)
              Do
                Do
                  i= i+1;  If (xx(i) >= xv) Exit
                End do

                Do
                  j= j-1;  If (xx(j) <= xv) Exit
                End do

                If (j <= i) Exit

                If (xx(i) > xx(j)) then
                  xt= xx(i);  kt= key(i)
                  xx(i)= xx(j);  key(i)= key(j)
                  xx(j)= xt;  key(j)= kt
                End if
              End do

            Else  ! No key

!             Sort left(l), middle(k), and right(r) entries

              If (xx(l) > xx(r)) then
                xt= xx(l);  xx(l)= xx(r);  xx(r)= xt
              End if

              If (xx(k) > xx(r)) then
                xt= xx(k);  xx(k)= xx(r);  xx(r)= xt
              End if

              If (xx(l) > xx(k)) then
                xt= xx(l);  xx(l)= xx(k);  xx(k)= xt
              End if

!             Swap the entries [l+1, r-1] as necessary so that if
!             xx(i) <= xx(k) < xx(j) or xx(i) < xx(k) <= xx(j)  then i < j.

              xv= xx(k)
              Do
                Do
                  i= i+1;  If (xx(i) >= xv) Exit
                End do

                Do
                  j= j-1;  If (xx(j) <= xv) Exit
                End do

                If (j <= i) Exit

                If (xx(i) > xx(j)) then
                  xt= xx(i);  xx(i)= xx(j);  xx(j)= xt
                End if
              End do
            End if

            stack= stack + 1

            If (stack > stack_max) then
!              Quick_sort error: stack_max too small. Result may not be fully sorted.
              stack= -1;  Exit
            End if

            If (xx(i) == xv) then
              l2= i+1
            Else
              l2= i
            End if

            If (xx(j) == xv) then
              r1= j-1
            Else
              r1= j
            End if

            If (r-l2 >= r1-l) then
              high(stack)= r;  low(stack)= l2;  r= r1
            Else
              high(stack)= r1;  low(stack)= l;  l= l2
            End if
          End if
        End do

      If (Present(key) .and. Present(ez)) then
        Forall(o=1:nn) xx(o)= xx(o) - (key(o)-1)*ez
      End if  
      
      If (.not.Ascending) then  
        xx(:nn)= xx(nn:1:-1)  
        If (Present(key)) key(:nn)= nn+1 - key(nn:1:-1)
      End if
      
    End Subroutine Quick_Sort_real


    Pure Subroutine Quick_Sort_double (Ascending, xx, key, ReOrd, ez)

!     Apply a quick sort algorithm to real data 'xx'.   

!     Optional arguments:    
    
!     Compute the 'key' for further sorting.  That is, to sort yy 
!     in the same way as xx, let yy= yy(key).
    
!     'ReOrd' is true if 'xx' needs to be reordered. Otherwise it is 
!     already in correct order.
    
!     'ez' is an epsilon value for determining the ordering of successive 
!     elements of 'xx' (after reordering) which are almost identical 
!     (within 'ez' of each other). For such elements the final ordering
!     will be the same as the original ordering (this 'stable' property
!     is unlikely with quick sort applied to equal elements).

      Logical,                 Intent(in) :: Ascending  ! Sort into increasing order if true, else decreasing
      Real(Dblp),        Intent(inout) :: xx(:)      ! (nn)  The data to be sorted
      Integer, Optional,      Intent(out) :: key(:)     ! (nn)  The sorting key
      Logical, Optional,      Intent(out) :: ReOrd      ! True if 'xx' is re-ordered otherwise the 
                                                        !   original ordering is correct
      Real(Dblp), Optional, Intent(in) :: ez         ! Use this value to retain the original order for 
                                                        !   elements of 'xx' which are almost equal 
                                                        !   (within 'ez', successively). Requires 'key'.
!   Local:
      Integer, Parameter :: n_max= 10, stack_max= 200
      Integer :: stack, low(stack_max), high(stack_max)
      Integer :: i, j, k, l, o, r, it, iv, kt, l2, r1, nn
      Real(Dblp)  :: xv, xt

      nn= Size(xx);  If (Present(ReOrd)) ReOrd= .false.
      If (nn <= 0)  Return
      
      If (Present(key)) key(:nn)= 'ID'
      If (nn == 1)  Return
      
      If (Present(key)) then
        If (Size(key) < nn) then  ! Error condition
          key= -1;  Return
        End if
      End if
      
      If (Ascending) then
        If (Present(ez)) then
          Do i= 1,nn-1
            If (xx(i) > xx(i+1) - ez) Exit
          End do
        Else
          Do i= 1,nn-1
            If (xx(i) > xx(i+1)) Exit
          End do
        End if
      Else
        If (Present(ez)) then
          Do i= 1,nn-1
            If (xx(i) < xx(i+1) + ez) Exit
          End do
        Else
          Do i= 1,nn-1
            If (xx(i) < xx(i+1)) Exit
          End do
        End if
      End if
      
      If (Present(ReOrd)) ReOrd= i < nn
      If (i >= nn) Return

!     Sort as if ascending
      
      If (.not.Ascending) xx(:nn)= xx(nn:1:-1)
      
      If (Present(key) .and. Present(ez)) then
        Forall(o=1:nn) xx(o)= xx(o) + (o-1)*ez
      End if  

      stack= 0;  l= 1;  r= nn

        Do
          If (r-l < n_max) then ! Insertion sort for the interval [l=left, r=right]
            If (Present(key)) then
              Do j= l+1,r
                If (xx(j-1) <= xx(j))  Cycle
                xv= xx(j);  kt= key(j)
                Do i= j-1,l,-1
                  If (xx(i) <= xv)  Exit
                  xx(i+1)= xx(i);  key(i+1)= key(i)
                End do
                xx(i+1)= xv;  key(i+1)= kt
              End do

            Else   ! No key
              Do j= l+1,r
                If (xx(j-1) <= xx(j))  Cycle
                xv= xx(j)
                Do i= j-1,l,-1
                  If (xx(i) <= xv)  Exit
                  xx(i+1)= xx(i)
                End do
                xx(i+1)= xv
              End do
            End if

            If (stack <= 0) Exit

            r= high(stack);  l= low(stack)
            stack= stack - 1

          Else   ! Quick sort partitioning step for the interval [l=left, r=right]

            k= (l+r)/2;  i= l;  j= r

            If (Present(key)) then

!             Sort left(l), middle(k), and right(r) entries

              If (xx(l) > xx(r)) then
                xt= xx(l);  kt= key(l)
                xx(l)= xx(r);  key(l)= key(r)
                xx(r)= xt;  key(r)= kt
              End if

              If (xx(k) > xx(r)) then
                xt= xx(k);  kt= key(k)
                xx(k) = xx(r);  key(k)= key(r)
                xx(r) = xt;  key(r)= kt
              End if

              If (xx(l) > xx(k)) then
                xt= xx(l);  kt= key(l)
                xx(l)= xx(k);  key(l)= key(k)
                xx(k)= xt;  key(k)= kt
              End if

!             Swap the entries [l+1, r-1] as necessary so that if
!             xx(i) <= xx(k) < xx(j) or xx(i) < xx(k) <= xx(j)  then i < j.

              xv= xx(k)
              Do
                Do
                  i= i+1;  If (xx(i) >= xv) Exit
                End do

                Do
                  j= j-1;  If (xx(j) <= xv) Exit
                End do

                If (j <= i) Exit

                If (xx(i) > xx(j)) then
                  xt= xx(i);  kt= key(i)
                  xx(i)= xx(j);  key(i)= key(j)
                  xx(j)= xt;  key(j)= kt
                End if
              End do

            Else  ! No key

!             Sort left(l), middle(k), and right(r) entries

              If (xx(l) > xx(r)) then
                xt= xx(l);  xx(l)= xx(r);  xx(r)= xt
              End if

              If (xx(k) > xx(r)) then
                xt= xx(k);  xx(k)= xx(r);  xx(r)= xt
              End if

              If (xx(l) > xx(k)) then
                xt= xx(l);  xx(l)= xx(k);  xx(k)= xt
              End if

!             Swap the entries [l+1, r-1] as necessary so that if
!             xx(i) <= xx(k) < xx(j) or xx(i) < xx(k) <= xx(j)  then i < j.

              xv= xx(k)
              Do
                Do
                  i= i+1;  If (xx(i) >= xv) Exit
                End do

                Do
                  j= j-1;  If (xx(j) <= xv) Exit
                End do

                If (j <= i) Exit

                If (xx(i) > xx(j)) then
                  xt= xx(i);  xx(i)= xx(j);  xx(j)= xt
                End if
              End do
            End if

            stack= stack + 1

            If (stack > stack_max) then
!              Quick_sort error: stack_max too small. Result may not be fully sorted.
              stack= -1;  Exit
            End if

            If (xx(i) == xv) then
              l2= i+1
            Else
              l2= i
            End if

            If (xx(j) == xv) then
              r1= j-1
            Else
              r1= j
            End if

            If (r-l2 >= r1-l) then
              high(stack)= r;  low(stack)= l2;  r= r1
            Else
              high(stack)= r1;  low(stack)= l;  l= l2
            End if
          End if
        End do

      If (Present(key) .and. Present(ez)) then
        Forall(o=1:nn) xx(o)= xx(o) - (key(o)-1)*ez
      End if  
      
      If (.not.Ascending) then  
        xx(:nn)= xx(nn:1:-1)  
        If (Present(key)) key(:nn)= nn+1 - key(nn:1:-1)
      End if
      
    End Subroutine Quick_Sort_double


    Pure Subroutine Matrix_diagonal_op_const_int (Op, const, mat)
      Character(3), Intent(in) :: Op
      Integer,      Intent(in) :: const
      Integer,   Intent(inout) :: mat(:,:)  ! (n,n+)
!   Local:
      Integer :: i, n

      n= Size(mat,1)

      Select case (Op)
        Case ('Add')
          Forall(i=1:n)  mat(i,i)= mat(i,i) + const
        Case ('Mul')
          Forall(i=1:n)  mat(i,i)= mat(i,i) * const
        Case default
          Forall(i=1:n)  mat(i,i)= const
      End select
    End Subroutine Matrix_diagonal_op_const_int

    
    Pure Subroutine Matrix_diagonal_op_const_real (Op, const, mat)
      Character(3), Intent(in) :: Op
      Real,         Intent(in) :: const
      Real,      Intent(inout) :: mat(:,:)  ! (n,n+)
!   Local:
      Integer :: i, n

      n= Size(mat,1)

      Select case (Op)
        Case ('Add')
          Forall(i=1:n)  mat(i,i)= mat(i,i) + const
        Case ('Mul')
          Forall(i=1:n)  mat(i,i)= mat(i,i) * const
        Case default
          Forall(i=1:n)  mat(i,i)= const
      End select
    End Subroutine Matrix_diagonal_op_const_real

   
    Pure Subroutine Matrix_diagonal_op_const_double (Op, const, mat)
      Character(3),     Intent(in) :: Op
      Real(Dblp),    Intent(in) :: const
      Real(Dblp), Intent(inout) :: mat(:,:)  ! (n,n+)
!   Local:
      Integer :: i, n

      n= Size(mat,1)

      Select case (Op)
        Case ('Add')
          Forall(i=1:n)  mat(i,i)= mat(i,i) + const
        Case ('Mul')
          Forall(i=1:n)  mat(i,i)= mat(i,i) * const
        Case default
          Forall(i=1:n)  mat(i,i)= const
      End select
    End Subroutine Matrix_diagonal_op_const_double


    Pure Subroutine Matrix_diagonal_op_vec_int (Op, vec, mat)
      Character(3), Intent(in) :: Op
      Integer,      Intent(in) :: vec(:)    ! (n)
      Integer,   Intent(inout) :: mat(:,:)  ! (n,n+)
!   Local:
      Integer :: i, n

      n= Size(vec)

      Select case (Op)
        Case ('Add')
          Forall(i=1:n)  mat(i,i)= mat(i,i) + vec(i)
        Case ('Mul')
          Forall(i=1:n)  mat(i,i)= vec(i) * mat(i,i)
        Case default
          Forall(i=1:n)  mat(i,i)= vec(i)
      End select
    End Subroutine Matrix_diagonal_op_vec_int

    
    Pure Subroutine Matrix_diagonal_op_vec_real (Op, vec, mat)
      Character(3), Intent(in) :: Op
      Real,         Intent(in) :: vec(:)    ! (n)
      Real,      Intent(inout) :: mat(:,:)  ! (n,n+)
!   Local:
      Integer :: i, n

      n= Size(vec)

      Select case (Op)
        Case ('Add')
          Forall(i=1:n)  mat(i,i)= mat(i,i) + vec(i)
        Case ('Mul')
          Forall(i=1:n)  mat(i,i)= vec(i) * mat(i,i)
        Case default
          Forall(i=1:n)  mat(i,i)= vec(i)
      End select
    End Subroutine Matrix_diagonal_op_vec_real

   
    Pure Subroutine Matrix_diagonal_op_vec_double (Op, vec, mat)
      Character(3),     Intent(in) :: Op
      Real(Dblp),    Intent(in) :: vec(:)    ! (n)
      Real(Dblp), Intent(inout) :: mat(:,:)  ! (n,n+)
!   Local:
      Integer :: i, n

      n= Size(vec)

      Select case (Op)
        Case ('Add')
          Forall(i=1:n)  mat(i,i)= mat(i,i) + vec(i)
        Case ('Mul')
          Forall(i=1:n)  mat(i,i)= mat(i,i) * vec(i)
        Case default
          Forall(i=1:n)  mat(i,i)= vec(i)
      End select
    End Subroutine Matrix_diagonal_op_vec_double


    Pure Subroutine Matrix_diagonal_out_vec_int (mat, vec)
      Integer,  Intent(in) :: mat(:,:)  ! (n,n+)
      Integer, Intent(out) :: vec(:)    ! (n)
!   Local:
      Integer :: i

      Forall(i=1:Size(vec))  vec(i)= mat(i,i)
    End Subroutine Matrix_diagonal_out_vec_int


    Pure Subroutine Matrix_diagonal_out_vec_real (mat, vec)
      Real,  Intent(in) :: mat(:,:)  ! (n,n+)
      Real, Intent(out) :: vec(:)    ! (n)
      Integer :: i

      Forall(i=1:Size(vec))  vec(i)= mat(i,i)
    End Subroutine Matrix_diagonal_out_vec_real


    Pure Subroutine Matrix_diagonal_out_vec_double (mat, vec)
      Real(Dblp),  Intent(in) :: mat(:,:)  ! (n,n+)
      Real(Dblp), Intent(out) :: vec(:)    ! (n)
      Integer :: i

      Forall(i=1:Size(vec))  vec(i)= mat(i,i)
    End Subroutine Matrix_diagonal_out_vec_double


   Elemental Subroutine Swap_I (a,b)
     Integer, Intent(inout) :: a,b
     Integer                 :: dum
     dum=a;  a=b;  b=dum
   End Subroutine Swap_I
   

   Elemental Subroutine Swap_R (a,b)
     Real, Intent(inout) :: a,b
     Real                 :: dum
     dum=a;  a=b;  b=dum
   End Subroutine Swap_R
   

   Elemental Subroutine Swap_D (a,b)
     Real(Dblp), Intent(inout) :: a,b
     Real(Dblp)                 :: dum
     dum=a;  a=b;  b=dum
   End Subroutine Swap_D
   

   Pure Subroutine Swap_IV (a,b, mask)
     Integer, Intent(inout) :: a(:), b(:)
     Logical, Optional, Intent(in) :: mask(:)

     Integer :: dum(Size(a))

     If (Present(mask)) then
       Where (mask)
         dum= a;  a= b;  b= dum
       End where
     Else
       dum= a;  a= b;  b= dum
     End if
   End Subroutine Swap_IV
   

   Pure Subroutine Swap_RV (a,b, mask)
     Real,          Intent(inout) :: a(:), b(:)
     Logical, Optional, Intent(in) :: mask(:)

     Real :: dum(Size(a))

     If (Present(mask)) then
       Where (mask)
         dum= a;  a= b;  b= dum
       End where
     Else
       dum= a;  a= b;  b= dum
     End if
   End Subroutine Swap_RV
   

   Pure Subroutine Swap_DV (a,b, mask)
     Real(Dblp), Intent(inout) :: a(:), b(:)
     Logical, Optional, Intent(in) :: mask(:)

     Real(Dblp) :: dum(Size(a))

     If (Present(mask)) then
       Where (mask)
         dum= a;  a= b;  b= dum
       End where
     Else
       dum= a;  a= b;  b= dum
     End if
   End Subroutine Swap_DV
   

   Pure Function Outer_productR (V1,V2)
     Real, Intent(in) :: V1(:)  ! (n)
     Real, Intent(in) :: V2(:)  ! (m)
     Real :: Outer_productR(Size(V1),Size(V2))
     Integer :: i
     Forall(i=1:Size(V1))  Outer_productR(i,:)= V1(i) * V2
   End Function Outer_productR

   
   Pure Function Outer_productD (V1,V2)
     Real(Dblp), Intent(in) :: V1(:)  ! (n)
     Real(Dblp), Intent(in) :: V2(:)  ! (m)
     Real(Dblp) :: Outer_productD(Size(V1),Size(V2))
     Integer :: i
     Forall(i=1:Size(V1))  Outer_productD(i,:)= V1(i) * V2
   End Function Outer_productD


   Pure Real Function Trace_R (Mat)
     Real, Intent(in) :: Mat(:,:)  ! (n,n)
!  Local:
     Real    :: ss
     Integer :: i

     ss= 0.0
     Do i= 1,Size(Mat,1)
       ss= ss + Mat(i,i)
     End do
     Trace_R= ss
   End Function Trace_R


   Pure Function Trace_D (Mat)
     Real(Dblp), Intent(in) :: Mat(:,:)  ! (n,n)
     Real(Dblp)             :: Trace_D
!  Local:
     Real(Dblp) :: ss
     Integer       :: i

     ss= 0.0
     Do i= 1,Size(Mat,1)
       ss= ss + Mat(i,i)
     End do
     Trace_D= ss
   End Function Trace_D

!  Output angles in degrees clockwise wrt north

   Elemental Real Function Output_angleR (in_ang)
     Real, Intent(in) :: in_ang
     Output_angleR= Modulo(90.000001 - rad_to_deg*in_ang, 360.0)
   End Function Output_angleR

   Elemental Real Function Output_angleD (in_ang)
     Real(Dblp), Intent(in) :: in_ang
     Output_angleD= Modulo(90.000001 - rad_to_deg*in_ang, 360.0)
   End Function Output_angleD

   Integer Function Subset (Set1, Set2, ni,Intr) 
   
!    Is Set1 a proper subset of Set2 (Subset = 1), 
!    or is Set2 a proper subset of Set1 (Subset = 2),
!    or are they equal (Subset = 0), or overlapping
!    but neither is a subset of the other (Subset = -1),
!    or non-overlapping (Subset = -2, or ni = 0?

!    Assume that both sets are represented as lists of integers 
!    in increasing order. Also compute the intersection set
!    Intr(:ni) of Set1 and Set2.

     Integer,  Intent(in) :: Set1(:)  ! (n1)
     Integer,  Intent(in) :: Set2(:)  ! (n2)
     Integer, Intent(out) :: ni       ! Size of the intersection of Set1 and Set2
     Integer, Intent(out) :: Intr(:)  ! (Min(n1,n2)=>ni) The intersection
!  Local:
     Logical :: Intersect
	 Integer :: mi, n1, n2

     n1= Size(Set1);  n2= Size(Set2);  mi= Min(n1,n2)
     
     If (Size(Intr) < mi) then
       Call Out ("Error in 'Subset'. Intesection set too small")
       Call Out ("First set size",n1, "Second set size",n2)
       Call Out ("Intesection set size", mi)
       Stop
     End if
     
     Subset= -1;  ni= 0;  Intr= -1
     
     If (n1 <= 0 .and. n2 <= 0) then
       Subset= 0
     Else if (n1 == 0) then
       Subset= 1
     Else if (n2 == 0) then
       Subset= 2
     Else if (Set1(n1) < Set2(1) .or. Set2(n2) < Set1(1)) then
       Subset= -2
     Else if (All(Set1(:mi) == Set2(:mi))) then
       If (mi == n1) then
         Subset= 1;  ni= n1;  Intr= Set1
         If (n1 == n2) Subset= 0
       Else
         Subset= 2;  ni= n2;  Intr= Set2
       End if
     Else
        Call Two_set_analysis (Set1,Set2, ni,Intr)        
        
        If (ni == 0) then
          Subset= -2
        Else if (ni == n1) then
          Subset= 1
        Else if (ni == n2) then
          Subset= 2
        End if
     End if
   End Function Subset


   Pure Subroutine Insert_in_list (v,i, n,ls, ls_in)

!    Insert "v" in the ordered integer list "ls" of length "n"
!    Insert at position "i" if i > 0, else compute i. 
!    'n' is increased by 1.

     Integer,    Intent(in) :: v
     Integer, Intent(inout) :: i, n
     Integer, Intent(inout) :: ls(:)
     Integer, Optional, Intent(in) :: ls_in(:)

     If (i <= 0) then
       If (Present(ls_in)) then
         i= First_true(1, ls_in(1:n) > v)
       Else
         i= First_true(1, ls(1:n) > v)
       End if
     End if

     If (Present(ls_in)) then
       If (i > 1)  ls(1:i-1)  = ls_in(1:i-1)
       If (i <= n) ls(i+1:n+1)= ls_in(i:n)
     Else
       If (i <= n) ls(i+1:n+1)= ls(i:n)
     End if

     ls(i)= v;  n= n + 1  ! Insert "v"
   End Subroutine Insert_in_list

   
   Pure Subroutine Remove_from_list (Search_forward, t, n,ls)

!    Remove "t" from the ordered integer list "ls" of length "n" if present

     Logical,    Intent(in) :: Search_forward
     Integer,    Intent(in) :: t
     Integer, Intent(inout) :: n, ls(:)
!  Local:
     Integer :: i

     If (Search_forward) then
       i= First_true(1, ls(1:n) == t)
       If (i <= n) then
         ls(i:n-1)= ls(i+1:n);  n= n - 1
       End if
     Else
       i= Last_true(ls(1:n) == t)
       If (i >= 1) then
         ls(i:n-1)= ls(i+1:n);  n= n - 1
       End if
     End if
   End Subroutine Remove_from_list


   Pure Subroutine Remove_from_matrixI (ind, mat)

!    Remove the row / colum 'ind' from the matrix 'mat'

     Integer, Intent(inout) :: ind       ! index to be removed. -1 on output if error
     Integer, Intent(inout) :: mat(:,:)  ! (n,n)
!  Local:
     Integer :: i, j, n, m, i1, m1, n1

     n= Size(mat,1);  m= Size(mat,2)
     If (ind < 1 .or. ind > Min(n,m)) Return
     
     i1= ind + 1;  m1= m - 1;  n1= n - 1
     If (ind < m) then
       Do i= 1,n
         mat(i,ind:m1)= mat(i,i1:)
       End do
     End if
     mat(:,m)= 0
     
     If (ind < n) then
       Do j= 1,m1
         mat(ind:n1,j)= mat(i1:,j)
       End do
     End if
     mat(n,:)= 0

   End Subroutine Remove_from_matrixI


   Pure Subroutine Remove_from_matrixR (ind, mat)

!    Remove the row / colum 'ind' from the matrix 'mat'

     Integer, Intent(inout) :: ind       ! index to be removed. -1 on output if error
     Real,    Intent(inout) :: mat(:,:)  ! (n,n)
!  Local:
     Integer :: i, j, n, m, i1, m1, n1

     n= Size(mat,1);  m= Size(mat,2)
     If (ind < 1 .or. ind > Min(n,m)) Return
     
     i1= ind + 1;  m1= m - 1;  n1= n - 1
     If (ind < m) then
       Do i= 1,n
         mat(i,ind:m1)= mat(i,i1:)
       End do
     End if
     mat(:,m)= 0
     
     If (ind < n) then
       Do j= 1,m1
         mat(ind:n1,j)= mat(i1:,j)
       End do
     End if
     mat(n,:)= 0

   End Subroutine Remove_from_matrixR

   
   Pure Subroutine Two_set_analysis (Set1,Set2, ni,Intr, Sub1,Sub2, &
                                     Unyn,Ind1,Ind2, OneMs,TwoMs)
   
!    Compute the intersection and union of two sets 'Set1' and 'Set2', plus 
!    associated lists
   
     Integer,  Intent(in) :: Set1(:)  ! (n1) The first set in increasing order
     Integer,  Intent(in) :: Set2(:)  ! (n2) The second set in increasing order
     
     Integer, Intent(out) :: ni       ! The size of the intersection set
     Integer, Intent(out) :: Intr(:)  ! (mi=>ni) The intersection set, initially of size mi= Min(n1,n2)
     
     Integer, Optional, Intent(out) :: Sub1(:)  ! (mi=>ni) Indices to get the intersection from Set1: 
                                                !          Intr(:ni)= Set1(Sub1(:ni))
     Integer, Optional, Intent(out) :: Sub2(:)  ! (mi=>ni) Indices to get the intersection from Set2: 
                                                !          Intr(:ni)= Set2(Sub2(:ni)))
     
     Integer, Optional, Intent(out) :: Unyn(:)  ! (n1+n2=>nu) The union set in increasing order
     Integer, Optional, Intent(out) :: Ind1(:)  ! (n1) Indices that yield Set1 in Unyn: Set1= Unyn(Ind1)
     Integer, Optional, Intent(out) :: Ind2(:)  ! (n2) Indices that yield Set2 in Unyn: Set2= Unyn(Ind2)
     
     Integer, Optional, Intent(out) :: OneMs(:) ! (n1=>n1-ni) Indices of Set1 - Set2 as a subset of Set1
     Integer, Optional, Intent(out) :: TwoMs(:) ! (n2=>n2-ni) Indices of Set2 - Set1 as a subset of Set2
     
!  Local:
     Integer :: i, j, k, l, i1, j1, n1, n2, nu
     
     n1= Size(Set1);  n2= Size(Set2)
     Intr= -1;  
     If (Present(Unyn)) Unyn= -1
     If (Present(Sub1)) Sub1= -1
     If (Present(Sub2)) Sub2= -1
     If (Present(Ind1)) Ind1= -1
     If (Present(Ind2)) Ind2= -1
     If (Present(OneMs)) OneMs= -1
     If (Present(TwoMs)) TwoMs= -1
     
     i= 1;  j= 1;  k= 0;  l= 0;  i1= 0;  j1= 0
     
     Do While (i <= n1 .or. j <= n2)
       If (j > n2) then
         l= l + 1;  If (Present(Unyn)) Unyn(l)= Set1(i)
         If (Present(Ind1)) Ind1(i)= l
         If (Present(OneMs)) then
           i1= i1 + 1;  OneMs(i1)= i
         End if
         i= i + 1
       Else if (i > n1) then
         l= l + 1;  If (Present(Unyn)) Unyn(l)= Set2(j)  
         If (Present(Ind2)) Ind2(j)= l
         If (Present(TwoMs)) then
           j1= j1 + 1;  TwoMs(j1)= j
         End if
         j= j + 1
       Else if (Set1(i) < Set2(j)) then
         l= l + 1;  If (Present(Unyn)) Unyn(l)= Set1(i)
         If (Present(Ind1)) Ind1(i)= l
         If (Present(OneMs)) then
           i1= i1 + 1;  OneMs(i1)= i
         End if
         i= i + 1
       Else if (Set1(i) > Set2(j)) then
         l= l + 1;  If (Present(Unyn)) Unyn(l)= Set2(j)  
         If (Present(Ind2)) Ind2(j)= l
         If (Present(TwoMs)) then
           j1= j1 + 1;  TwoMs(j1)= j
         End if
         j= j + 1
       Else
         k= k + 1;  Intr(k)= Set1(i)
         l= l + 1;  If (Present(Unyn)) Unyn(l)= Set1(i)
         If (Present(Ind1)) Ind1(i)= l
         If (Present(Ind2)) Ind2(j)= l
         If (Present(Sub1)) Sub1(k)= i
         If (Present(Sub2)) Sub2(k)= j
         i= i + 1;  j= j + 1
       End if
     End do
     
     ni= k;  nu= l
   End Subroutine Two_set_analysis
   

   Pure Function Partial_sumsI (in_dat)
     Integer,  Intent(in) :: in_dat(:)
     Integer              :: Partial_sumsI(Size(in_dat))
     Integer :: i
     Partial_sumsI(1)= in_dat(1)
     Do i= 2,Size(in_dat)
       Partial_sumsI(i)= Partial_sumsI(i-1) + in_dat(i)
     End do
   End Function Partial_sumsI

   
   Pure Function Partial_sumsR(in_dat)
     Real,  Intent(in) :: in_dat(:)
     Real              :: Partial_sumsR(Size(in_dat))
     Integer :: i
     Partial_sumsR(1)= in_dat(1)
     Do i= 2,Size(in_dat)
       Partial_sumsR(i)= Partial_sumsR(i-1) + in_dat(i)
     End do
   End Function Partial_sumsR

   Pure Function Partial_sumsD (in_dat)
     Real(Dblp), Intent(in) :: in_dat(:)
     Real(Dblp)             :: Partial_sumsD(Size(in_dat))
     Integer :: i
     Partial_sumsD(1)= in_dat(1)
     Do i= 2,Size(in_dat)
       Partial_sumsD(i)= Partial_sumsD(i-1) + in_dat(i)
     End do
   End Function Partial_sumsD

!  Compute a random permutation of 1...n
   
   Function Random_permutation(n)
     Integer,  Intent(in) :: n
     Integer :: Random_permutation(n)
!  Local:
     Integer :: i, k, l, ls(n)
     Real    :: rnd

     If (n <= 0) Return

     If (n == 1) then
       Random_permutation(1)= 1
     Else
       l= n;  ls= 'ID'
       Do i= 1,n
         Call Random_number(rnd)
         k= Ceiling(rnd * l);  Random_permutation(i)= ls(k)
         If (k < l) ls(k:l-1)= ls(k+1:l);  l= l - 1
       End do
     End if
   End Function Random_permutation

!  Compute a random subset of given size of a given set   

   Subroutine Random_subset (n,m, universe, sub)
     Integer,   Intent(in) :: n
     Integer,   Intent(in) :: m
     Integer,   Intent(in) :: universe(:)  ! (m)
     Integer,  Intent(out) :: sub(:)       ! (n)
!  Local:
     Integer :: i, k, l, ls(m)
     Real    :: rnd

     If (n >= m) then
       sub(1:n)= universe
     Else if (n > 0) then
       l= m;  ls= universe
       Do i= 1,n
         Call Random_number(rnd)
         k= Ceiling(rnd * l);  sub(i)= ls(k)
         If (k < l) ls(k:l-1)= ls(k+1:l);  l= l - 1
       End do
       Call Sort (.true.,sub)
     End if
   End Subroutine Random_subset

!  Test a mapping to see if it is a permutation
   
   Pure Logical Function Permutation (l, perm)
     Integer, Intent(in) :: l
     Integer, Intent(in) :: perm(l:)
!  Local:
     Logical :: targ(l:Ubound(perm,1))

     Permutation= All(perm >= l .and. perm <= Ubound(perm,1))
     If (.not.Permutation) Return

     targ= .false.;  targ(perm)= .true.
     Permutation= All(targ)
   End Function Permutation

  
  Subroutine Matrix_norms (A, norms, iter,error)

!   Compute the L_inf, L_1, Frobenius, and L_2 (or spectral radius) matrix norms

    Real(Dblp),            Intent(in) :: A(:,:)   ! (m,n) Possibly non-symmetric matrix
    Real(Dblp),           Intent(out) :: norms(:) ! (4) L_inf, L_1, Frobenius, and L_2 norms
    Integer, Optional,       Intent(out) :: iter     ! # iterations for the L_2 norm
    Real(Dblp), Optional, Intent(out) :: error    ! Final fractional error of the L_2 norm
! Local:
    Real(Dblp),   Parameter :: eps= 1.0E-5
    Real(Dblp), Allocatable :: x(:), y(:)
    Real(Dblp) :: err
    Integer       :: i, j, m, n, it

    norms= 0;  m= Size(A,1);  n= Size(A,2);  Allocate (x(n), y(m))

!   L_infinity norm    
    
    Do i= 1,m
      y(i)= Sum(Abs(A(i,:)))
    End do
    norms(1)= Maxval(y)
    
!    L_1 norm    

    Do j= 1,n
      x(j)= Sum(Abs(A(:,j)))
    End do
    norms(2)= Maxval(x)
    
!   Frobenius norm    

    norms(3)= Sqrt(Sum(A**2))
    
!   Spectral radius = L_2 norm
    
    Call Spectral_radius (eps, A, norms(4), it,err)
    
    If (Present(iter))  iter = it
    If (Present(error)) error= err

  End Subroutine Matrix_norms

  
  Subroutine Spectral_radius (eps, A, norm, it,err)

!   Compute the matrix two norm of A [=Sqrt(max eigenvalue of Trans(A)*A],
!   or spectral radius, by the power method. 

    Real(Dblp),  Intent(in) :: eps     ! Convergence tolerance for norm**2
    Real(Dblp),  Intent(in) :: A(:,:)  ! (m,n) Input matrix
    Real(Dblp), Intent(out) :: norm    ! The two norm = spectral radius of A
    Integer,       Intent(out) :: it      ! # iterations to convergence
    Real(Dblp), Intent(out) :: err     ! Final fractional error of norm**2
! Local:
    Integer, Parameter :: itmax= 30
    Real(Dblp), Parameter   :: del= 1.0E-12
    Real(Dblp), Allocatable :: x(:), y(:), z(:), At(:,:)
    Real(Dblp) :: eig, old_eig
    Integer       :: m, n
    
    norm= 0;  it= 0;  err= 0
    m= Size(A,1);  n= Size(A,2)
    Allocate (x(n), y(m), z(n), At(n,m))
    At= Transpose(A);  old_eig= -1
    
    Call Random_seed;  Call Random_number(z)
    x= z / Sqrt(Sum(z**2))

    Do it= 1,itmax
      y= Matmul(A,x);  eig= Sum(y**2);  If (Abs(eig) < del) Return
      z= Matmul(At,y)
      x= z / Sqrt(Sum(z**2))
      err= Abs((eig - old_eig)/old_eig)
      If (err < eps) Exit
      old_eig= eig
    End do

    norm= Sqrt(eig)
  End Subroutine Spectral_radius

  
  Pure Subroutine ReOrder_matrix (key, A)

!   Reorder a matrix 'A' according to the mapping 'key'.
!   Assume that n = Size(key) <= m where m= Min(Size(A,1), Size(A,2))
!   with the range of 'key' a subset of 1...m.

    Integer, Intent(in) :: key(:)  ! The reordering
    Real, Intent(inout) :: A(:,:)  ! The matrix to be reordered
! Local:
    Integer :: i, n

    n= Size(key)

    Do i= 1,Size(A,1)
      A(i,:n)= A(i,key)
    End do
    Do i= 1,Size(A,2)
      A(:n,i)= A(key,i)
    End do

  End Subroutine ReOrder_matrix


   Subroutine Linear_regression (n, x,y, a,b, var_a,var_b,cov_ab, var_py)

     Integer,  Intent(in) :: n      ! # data pairs (x,y)
     Real,     Intent(in) :: x(:)   ! (n) The independent data
     Real,     Intent(in) :: y(:)   ! (n) The dependent data
     Real,    Intent(out) :: a, b   ! The computed coefficients for y= ax + b
     Real,    Intent(out) :: var_a  ! Variance of a
     Real,    Intent(out) :: var_b  ! Variance of b
     Real,    Intent(out) :: cov_ab ! Covariance between a and b
     Real,    Intent(out) :: var_py ! Variance of y from its predicted value ax+b
!  Local:
     Real    :: d, sx, sx2, sxy, sy

     sx= Sum(x);  sy= Sum(y);  sx2= Sum(x * x);  sxy= Sum(x * y)

     d= sx2 * n - sx * sx
     a= (sxy * n - sy * sx) / d
     b= (sx2 * sy - sx * sxy) / d

     var_a= n / d;  var_b= sx2 / d;  cov_ab= -sx / d
     var_py= Sum((a * x + b - y)**2) / (n-2)
   End Subroutine Linear_regression

   
   Subroutine Cubic_regression00 (n, x,y, a,b, var_a,var_b,cov_ab, var_py)

     Integer,  Intent(in) :: n      ! # data pairs (x,y)
     Real,     Intent(in) :: x(:)   ! (n) The independent data
     Real,     Intent(in) :: y(:)   ! (n) The dependent data
     Real,    Intent(out) :: a, b   ! The computed coefficients for y= a x^3 + b x^2
     Real,    Intent(out) :: var_a  ! Variance of a
     Real,    Intent(out) :: var_b  ! Variance of b
     Real,    Intent(out) :: cov_ab ! Covariance between a and b
     Real,    Intent(out) :: var_py ! Variance of y from its predicted value a x^3 + b x^2
!  Local:
     Real :: x2(n), x4(n), yx2(n)
     Real :: d, sx2y, sx3y, sx4, sx5, sx6

     x2= x * x;  yx2= y * x2;  sx2y= Sum(yx2);  sx3y= Sum(yx2 * x)
     x4= x2 * x2;  sx4= Sum(x4);  sx5= Sum(x * x4);  sx6= Sum(x2 * x4)

     d= sx6 * sx4 - sx5*sx5
     a= (sx3y*sx4 - sx2y*sx5) / d
     b= (sx6*sx2y - sx5*sx3y) / d

     var_a= sx4 / d;  var_b= sx6 / d;  cov_ab= -sx5 / d
     var_py= Sum(((a * x + b) * x2 - y)**2) / (n-2)

   End Subroutine Cubic_regression00


   Subroutine Quadratic_regression (n, x, y, a, b, c, var_py, var_a, var_b, var_c)

     Integer,  Intent(in) :: n         ! # data pairs (x,y)
     Real,     Intent(in) :: x(:)      ! (n) The independent data
     Real,     Intent(in) :: y(:)      ! (n) The dependent data
     Real,    Intent(out) :: a, b, c   ! The computed coefficients for y= ax**2 + bx + c
     Real,    Intent(out) :: var_py    ! Variance of y from its predicted value ax**2 + bx + c
     Real,    Intent(out) :: var_a     ! Variance of a
     Real,    Intent(out) :: var_b     ! Variance of b
     Real,    Intent(out) :: var_c     ! Variance of c
!  Local:
     Real    :: mean_x, mean_y, mean_x2, mean_x4, var_x, var_x2
     Real    :: cor_xy, cor_x2y, cor_x2x
     Real    :: d, m1, m2, m3, r1, r2

     mean_x= Sum(x) / n;  mean_y= Sum(y) / n
     mean_x2= Sum(x**2) / n;  mean_x4= Sum(x**4) / n

     var_x = Sum((x - mean_x)**2);  var_x2= Sum((x**2 - mean_x2)**2)
     cor_xy= Sum((x - mean_x) * (y - mean_y))
     cor_x2y= Sum((x**2 - mean_x2) * (y - mean_y))
     cor_x2x= Sum((x**2 - mean_x2) * (x - mean_x))

     m1= var_x2;  m2= cor_x2x;  m3= var_x
     r1= cor_x2y;  r2= cor_xy
     d= m1 * m3 - m2**2
     a= (m3*r1 - m2*r2) / d;  b= (m1*r2 - m2*r1) / d
     c= mean_y - (a*mean_x2 + b*mean_x)

     var_py= Sum((y - (a * x**2 + b * x + c))**2) / (n-3)

     var_a= var_py * Sqrt(m3 / d)
     var_b= var_py * Sqrt(m1 / d)
     var_c= var_a * mean_x4 + var_b * mean_x2

   End Subroutine Quadratic_regression


  Pure Subroutine Sum_constraint1 (v_sum, v)
!   Constrain the vector 'v' to sum to 'v_sum' by adjusting 
!   the last coordinate of 'v'
    Integer,    Intent(in) :: v_sum
    Real,    Intent(inout) :: v(:)

    Integer :: n, n1;  n= Size(v);  n1= n - 1
    v(n)= v_sum - Sum(v(1:n1))
  End Subroutine Sum_constraint1

  Pure Subroutine Sum_constraint2(v_sum, v)
!   Constrain the vectors 'v(:,i)' to sum to 'v_sum' by adjusting 
!   the last coordinate of each 'v(:,i)'
    Integer,    Intent(in) :: v_sum
    Real,    Intent(inout) :: v(:,:)

    Integer :: i, n, n1;  n= Size(v,1);  n1= n - 1
    Do i= 1,Size(v,2)
      v(n,i)= v_sum - Sum(v(1:n1,i))
    End do
  End Subroutine Sum_constraint2

  Pure Subroutine Sum_constraint1D (v_sum, v)
!   Constrain the vector 'v' to sum to 'v_sum' by adjusting 
!   the last coordinate of 'v'
    Integer,          Intent(in) :: v_sum
    Real(Dblp), Intent(inout) :: v(:)

    Integer :: n, n1;  n= Size(v);  n1= n - 1
    v(n)= v_sum - Sum(v(1:n1))
  End Subroutine Sum_constraint1D

  Pure Subroutine Sum_constraint2D(v_sum, v)
!   Constrain the vectors 'v(:,i)' to sum to 'v_sum' by adjusting 
!   the last coordinate of each 'v(:,i)'
    Integer,          Intent(in) :: v_sum
    Real(Dblp), Intent(inout) :: v(:,:)

    Integer :: i, n, n1;  n= Size(v,1);  n1= n - 1
    Do i= 1,Size(v,2)
      v(n,i)= v_sum - Sum(v(1:n1,i))
    End do
  End Subroutine Sum_constraint2D

  Pure Function Quad_rise1R (x2, x)

!   This quadratic function rises smoothly from the constant value 1 
!   for x <= 1 to x-c for x >= x2. 
!   This is done so as to match derivatives at both ends.

    Real, Intent(in) :: x2
    Real, Intent(in) :: x(:)

    Real :: Quad_rise1R(Size(x))
! Local:
    Real :: d, c

    d= 2 * (x2 - 1);  c= (x2 - 1) / 2

    Where (x >= x2)
      Quad_rise1R= x - c
    Elsewhere (x <= 1)
      Quad_rise1R= 1
    Elsewhere
      Quad_rise1R= 1 + (x - 1)**2 / d
    End where
  End Function Quad_rise1R

  Pure Function Quad_rise1D (x2, x)

!   This quadratic function rises smoothly from the constant value 1 
!   for x <= 1 to x-c for x >= x2. 
!   This is done so as to match derivatives at both ends.

    Real,          Intent(in) :: x2
    Real(Dblp), Intent(in) :: x(:)

    Real(Dblp) :: Quad_rise1D(Size(x))
! Local:
    Real(Dblp) :: d, c

    d= 2 * (x2 - 1);  c= (x2 - 1) / 2

    Where (x >= x2)
      Quad_rise1D= x - c
    Elsewhere (x <= 1)
      Quad_rise1D= 1
    Elsewhere
      Quad_rise1D= 1 + (x - 1)**2 / d
    End where
  End Function Quad_rise1D

  
  Pure Function Cubic_riseR (x1,x2, x)

!   This cubic function rises smoothly from the constant value 'x1' 
!   for x <= x1 to become the straight line f(x)= x for x >= x2 > x1.
!   This is done so as to match derivatives at both ends:
!   0 at x1 and 1 at x2.

!   Note that for d= x2 - x1 and dx= (x-x1)/d, the function is 
!   f(x) = x1 + dx**2 * (2-dx) * d for x1 < x < x2.
!   with f'(x) = (4 - 3*dx)*dx and f''(x) = (4 - 6dx) / d.
!   Then for dx = 2/3, f''= 0, f'= 4/3, and f= x1 + (16/27)d.
!   Also f''(x1)= 4/d, f''(x2)= -2/d.

    Real, Intent(in) :: x1
    Real, Intent(in) :: x2
    Real, Intent(in) :: x

    Real :: Cubic_riseR
! Local:
    Real :: d, dx

    d= x2 - x1
    If (x >= x2) then
      Cubic_riseR= x
    Else if (x <= x1) then
      Cubic_riseR= x1
    Else
      dx= (x - x1) / d
      Cubic_riseR= x1 + (dx)**2 * (2 - dx) * d
    End if
  End Function Cubic_riseR
  
  Pure Function Cubic_riseVR (x1,x2, x)

    Real, Intent(in) :: x1
    Real, Intent(in) :: x2
    Real, Intent(in) :: x(:)

    Real :: Cubic_riseVR(Size(x))
! Local:
    Real :: d, dx(Size(x))

    d= x2 - x1
    Where (x >= x2)
      Cubic_riseVR= x
    Elsewhere (x <= x1)
      Cubic_riseVR= x1
    Elsewhere
      dx= (x - x1) / d
      Cubic_riseVR= x1 + (dx)**2 * (2 - dx) * d
    End where
  End Function Cubic_riseVR

  Pure Function Cubic_riseD (x1,x2, x)

    Real, Intent(in) :: x1
    Real, Intent(in) :: x2
    Real(Dblp), Intent(in) :: x

    Real(Dblp) :: Cubic_riseD
! Local:
    Real(Dblp) :: d, dx

    d= x2 - x1
    If (x >= x2) then
      Cubic_riseD= x
    Else if (x <= x1) then
      Cubic_riseD= x1
    Else
      dx= (x - x1) / d
      Cubic_riseD= x1 + (dx)**2 * (2 - dx) * d
    End if
  End Function Cubic_riseD

  Pure Function Cubic_riseVD (x1,x2, x)

    Real, Intent(in) :: x1
    Real, Intent(in) :: x2
    Real(Dblp), Intent(in) :: x(:)

    Real(Dblp) :: Cubic_riseVD(Size(x))
! Local:
    Real(Dblp) :: d, dx(Size(x))

    d= x2 - x1
    Where (x >= x2)
      Cubic_riseVD= x
    Elsewhere (x <= x1)
      Cubic_riseVD= x1
    Elsewhere
      dx= (x - x1) / d
      Cubic_riseVD= x1 + (dx)**2 * (2 - dx) * d
    End where
  End Function Cubic_riseVD


  Pure Real Function Cube_quadR (x, x1,x2)

!   This function is a smoothed version of Max(1,y(x)) for the quadratic 
!   y(x)= c*x**2 with y(x1)= 1. Thus c= 1/x1**2. The smoothing is a via
!   a cubic function between x1 and x2, matching in both values and slopes
!   at x1 and x2.

    Real, Intent(in) :: x
    Real, Intent(in) :: x1, x2

! Local:
    Real :: c, a, b, d, d2, dx

    c= 1 / x1**2;  d= x2 - x1;  d2= d*d
    a= (c*x2*(x2+2*x1) - 3) / d2
    b= 2*(1 - c*x2*x1) / (d*d2)

    If (x >= x2) then
      Cube_quadR= c * x**2
    Else if (x <= x1) then
      Cube_quadR= 1
    Else
      dx= x - x1
      Cube_quadR= 1 + (a + b*dx)*dx**2
    End if
  End Function Cube_quadR

  Pure Function Cube_quadV (x, x1,x2)

!   This function is a smoothed version of Max(1,y(x)) for the quadratic 
!   y(x)= c*x**2 with y(x1)= 1. Thus c= 1/x1**2. The smoothing is a via
!   a cubic function between x1 and x2, matching in both values and slopes
!   at x1 and x2.

    Real, Intent(in) :: x(:)
    Real, Intent(in) :: x1, x2

    Real :: Cube_quadV(Size(x))
! Local:
    Real :: c, a, b, d, d2, dx(Size(x))

    c= 1 / x1**2;  d= x2 - x1;  d2= d*d
    a= (c*x2*(x2+2*x1) - 3) / d2
    b= 2*(1 - c*x2*x1) / (d*d2)

    Where (x >= x2)
      Cube_quadV= c * x**2
    Elsewhere (x <= x1)
      Cube_quadV= 1
    Elsewhere
      dx= x - x1
      Cube_quadV= 1 + (a + b*dx)*dx**2
    End where
  End Function Cube_quadV


  Pure Subroutine Quad_fit3 (x,fx, xm,fxm)

!   Do quadratic interpolation among 3 pairs of coordinates 
!   (x,fx), assuming x(1) < x(2) < x(3) and fx(1) <= fx(2) >= fx(3)
!   to get the value fxm of the quadratic function at xm, 
!   forcing x(1) <= xm <= x(3) if necessary.

!   However first do bisection, if necessary, to make the
!   quadratic interpolation more reliable (= more symmetric data).

    Real,  Intent(in) :: x(:)   ! (3) Coordinate values
    Real,  Intent(in) :: fx(:)  ! (3) Corresponding function values

    Real, Intent(out) :: xm     ! Value of the interpolating quadratic at which
                                ! a max or min occurs.

    Real, Intent(out) :: fxm    ! Corresponding function value
! Local:
    Real :: a, b, c, d, det, ratio, s(3), v(2), Inv(2,2)

    a= x(1)**2 - x(3)**2;  b= x(1) - x(3);  v(1)= fx(1) - fx(3)
    c= x(2)**2 - x(3)**2;  d= x(2) - x(3);  v(2)= fx(2) - fx(3)

    Inv(1,1)= d;  Inv(1,2)= -b
    Inv(2,1)= -c;  Inv(2,2)= a;  det= a*d - b*c
    s(:2)= Matmul(Inv,v) / det;  s(3)= fx(1) - (s(1)*x(1) + s(2))*x(1)
    ratio= (x(3) - x(2)) / (x(2) - x(1))

    If (fx(3) < 0.90 * fx(1) .or. ratio >= 4) then
      xm= (x(2) + x(3)) / 2
    Else if (fx(1) < 0.90 * fx(3) .or. ratio <= 0.25) then
      xm= (x(1) + x(2)) / 2
    Else
      xm= -s(2) / (2*s(1))
      If (xm < x(1)) then
        xm= x(1)
      Else if (xm > x(3)) then
        xm= x(3)
      End if
    End if

    fxm= (s(1)*xm + s(2))*xm + s(3)
  End Subroutine Quad_fit3

  
  Subroutine Cos_rise_drvR (begin_rise,end_rise, val, drv)

!   Compute the cosine rise value corresponding to 'val' for an increasing
!   cosine where the rise value is 0 before begin_rise and 1 after end_rise,
!   and a cosine function in between.

    Real,  Intent(in) :: begin_rise  ! Begin the cosine function at this value
    Real,  Intent(in) :: end_rise    ! End the cosine function at this value (> begin_rise)
    Real,  Intent(in) :: val         ! Compute the cosine at this value
    Real, Intent(out) :: drv(0:)     ! (0:n) Cosine rise value, possibly also its first and second derivatives
                                     ! 0 = cosine rise value at 'val'
                                     ! 1 = Derivative of this with respect to 'val'
                                     ! 2 = Second derivative of this with respect to 'val'
! Local:
    Real    :: x, fac, fac2, cx, sx
    Integer :: n

    n= Ubound(drv,1);  drv= 0
    
    If (val >= end_rise) then
      drv(0)= 1
    Else if (val >= begin_rise) then
      fac= PI / (end_rise - begin_rise)
      x= fac * (end_rise - val);  cx= Cos(x)
      drv(0)= (1 + cx) / 2
      
      If (n > 0) then
        sx= Sin(x);  fac2= fac / 2
        drv(1)= fac2 * sx
        If (n > 1) drv(2)= -fac * fac2 * cx
      End if
    End if
  End Subroutine Cos_rise_drvR
  

  Subroutine Cos_rise_drvD (begin_rise,end_rise, val, drv)

!   Compute the cosine rise value corresponding to 'val' for an increasing
!   cosine where the rise value is 0 before begin_rise and 1 after end_rise,
!   and a cosine function in between.

    Real(Dblp),  Intent(in) :: begin_rise  ! Begin the cosine function at this value
    Real(Dblp),  Intent(in) :: end_rise    ! End the cosine function at this value (> begin_rise)
    Real(Dblp),  Intent(in) :: val         ! Compute the cosine at this value
    Real(Dblp), Intent(out) :: drv(0:)     ! 0 = cosine rise value at 'val'
                                              ! 1 = Derivative of this with respect to 'val'
                                              ! 2 = Second derivative of this with respect to 'val'
! Local:
    Real(Dblp) :: x, fac, fac2, cx, sx
    Integer :: n

    n= Ubound(drv,1);  drv= 0
    
    If (val >= end_rise) then
      drv(0)= 1
    Else if (val >= begin_rise) then
      fac= PI / (end_rise - begin_rise)
      x= fac * (end_rise - val);  cx= Cos(x)
      drv(0)= (1 + cx) / 2
      
      If (n > 0) then
        sx= Sin(x);  fac2= fac / 2
        drv(1)= fac2 * sx
        If (n > 1) drv(2)= -fac * fac2 * cx
      End if
    End if
  End Subroutine Cos_rise_drvD
  

  Pure Real Function Cos_riseR (begin_rise,end_rise, val)

!   Compute the cosine rise value corresponding to 'val' for an increasing
!   cosine where the rise value is 0 before begin_rise and 1 after end_rise,
!   and a cosine function in between.

    Real, Intent(in) :: begin_rise  ! Begin the cosine function at this value
    Real, Intent(in) :: end_rise    ! End the cosine function at this value (> begin_rise)
    Real, Intent(in) :: val         ! Compute the cosine at this value
! Local:
    Real :: tmp

    If (val <= begin_rise) then
      Cos_riseR= 0
    Else if (val >= end_rise) then
      Cos_riseR= 1
    Else
      tmp= (end_rise - val) * PI / (end_rise - begin_rise)
      Cos_riseR= (1 + Cos(tmp)) / 2
    End if
  End Function Cos_riseR
  

  Pure Function Cos_riseD (begin_rise,end_rise, val)

!   Compute the cosine rise value corresponding to 'val' for an increasing
!   cosine where the rise value is 0 before begin_rise and 1 after end_rise,
!   and a cosine function in between.

    Real(Dblp), Intent(in) :: begin_rise  ! Begin the cosine function at this value
    Real(Dblp), Intent(in) :: end_rise    ! End the cosine function at this value (> begin_rise)
    Real(Dblp), Intent(in) :: val         ! Compute the cosine at this value
    
    Real(Dblp) :: Cos_riseD 
! Local:
    Real(Dblp) :: tmp

    If (val <= begin_rise) then
      Cos_riseD= 0
    Else if (val >= end_rise) then
      Cos_riseD= 1
    Else
      tmp= (end_rise - val) * PI / (end_rise - begin_rise)
      Cos_riseD= (1 + Cos(tmp)) / 2
    End if
  End Function Cos_riseD

  Pure Function Cos_riseVD (begin_rise,end_rise, val)

!   Compute the cosine rise value corresponding to 'val' for an increasing
!   cosine where the rise value is 0 before begin_rise and 1 after end_rise,
!   and a cosine function in between.

    Real(Dblp), Intent(in) :: begin_rise  ! Begin the cosine function at this value
    Real(Dblp), Intent(in) :: end_rise    ! End the cosine function at this value (> begin_rise)
    Real(Dblp), Intent(in) :: val(:)      ! Compute the cosine at this value
    
    Real(Dblp) :: Cos_riseVD(Size(val))
! Local:
    Real(Dblp) :: tmp(Size(val))

    Where (val <= begin_rise) 
      Cos_riseVD= 0
    Else where (val >= end_rise) 
      Cos_riseVD= 1
    Else where 
      tmp= (end_rise - val) * PI / (end_rise - begin_rise)
      Cos_riseVD= (1 + Cos(tmp)) / 2
    End where
  End Function Cos_riseVD

  Pure Function Cos_riseVR (begin_rise, end_rise, val)

!   Compute the cosine rise value corresponding to 'val' for an increasing
!   cosine where the rise value is 0 before begin_rise and 1 after end_rise,
!   and a cosine function in between.

    Real, Intent(in) :: begin_rise  ! Begin the cosine function at this value
    Real, Intent(in) :: end_rise    ! End the cosine function at this value (> begin_rise)
    Real, Intent(in) :: val(:)      ! Compute the cosine at this value
    
    Real :: Cos_riseVR(Size(val))  
! Local:
    Real :: tmp(Size(val))

    Where (val <= begin_rise)
      Cos_riseVR= 0
    Elsewhere (val >= end_rise)
      Cos_riseVR= 1
    Elsewhere
      tmp= (end_rise - val) * PI / (end_rise - begin_rise)
      Cos_riseVR= (1 + Cos(tmp)) / 2
    Endwhere
    
  End Function Cos_riseVR
  

  Pure Function Cos_riseMR (begin_rise, end_rise, val)

!   Compute the cosine rise value corresponding to 'val' for an increasing
!   cosine where the rise value is 0 before begin_rise and 1 after end_rise,
!   and a cosine function in between.

    Real, Intent(in) :: begin_rise  ! Begin the cosine function at this value
    Real, Intent(in) :: end_rise    ! End the cosine function at this value (> begin_rise)
    Real, Intent(in) :: val(:,:)    ! Compute the cosine at this value
    
    Real  :: Cos_riseMR(Size(val,1),Size(val,2))  
! Local:
    Real :: tmp(Size(val,1),Size(val,2))

    Where (val <= begin_rise)
      Cos_riseMR= 0
    Elsewhere (val >= end_rise)
      Cos_riseMR= 1
    Elsewhere
      tmp= (end_rise - val) * PI / (end_rise - begin_rise)
      Cos_riseMR= (1 + Cos(tmp)) / 2
    Endwhere
    
  End Function Cos_riseMR
  

  Pure Real Function Cos_fallR (begin_fall,end_fall, val, base)

!   Compute the cosine fall value corresponding to 'val' for a decreasing
!   cosine where the fall value is 1 before begin_fall and 0 after end_fall
!   and a cosine function in between.

    Real, Intent(in) :: begin_fall  ! Begin the cosine function at this value
    Real, Intent(in) :: end_fall    ! End the cosine function at this value ( > begin_fall)
    Real, Intent(in) :: val         ! Compute the cosine at this value
    Real, Optional, Intent(in) :: base  ! Falls to 'base' >= 0 instead of to 0.
! Local:
    Real :: tmp

    If (Present(base)) then
      If (val <= begin_fall) then
        Cos_fallR= 1
      Else if (val >= end_fall) then
        Cos_fallR= base
      Else 
        tmp= (val - begin_fall) * PI / (end_fall - begin_fall)
        Cos_fallR= base + (1 - base) * ((1 + Cos(tmp)) / 2)
      End if
    Else
      If (val <= begin_fall) then
        Cos_fallR= 1
      Else if (val >= end_fall) then
        Cos_fallR= 0
      Else 
        tmp= (val - begin_fall) * PI / (end_fall - begin_fall)
        Cos_fallR= (1 + Cos(tmp)) / 2
      End if
    End if
  End Function Cos_fallR

  Pure Real Function Cos_fallD (begin_fall,end_fall, val, base)

!   Compute the cosine fall value corresponding to 'val' for a decreasing
!   cosine where the fall value is 1 before begin_fall and 0 after end_fall
!   and a cosine function in between.

    Real(Dblp), Intent(in) :: begin_fall  ! Begin the cosine function at this value
    Real(Dblp), Intent(in) :: end_fall    ! End the cosine function at this value ( > begin_fall)
    Real(Dblp), Intent(in) :: val         ! Compute the cosine at this value
    Real(Dblp), Optional, Intent(in) :: base  ! Falls to 'base' >= 0 instead of to 0.
! Local:
    Real(Dblp) :: tmp

    If (Present(base)) then
      If (val <= begin_fall) then
        Cos_fallD= 1
      Else if (val >= end_fall) then
        Cos_fallD= base
      Else 
        tmp= (val - begin_fall) * PI / (end_fall - begin_fall)
        Cos_fallD= base + (1 - base) * ((1 + Cos(tmp)) / 2)
      End if
    Else
      If (val <= begin_fall) then
        Cos_fallD= 1
      Else if (val >= end_fall) then
        Cos_fallD= 0
      Else 
        tmp= (val - begin_fall) * PI / (end_fall - begin_fall)
        Cos_fallD= (1 + Cos(tmp)) / 2
      End if
    End if
  End Function Cos_fallD


  Pure Function Cos_fallVR (begin_fall,end_fall, val, base)

!   Compute the cosine fall value corresponding to 'val' for a decreasing
!   cosine where the fall value is 1 before begin_fall and 0 after end_fall
!   and a cosine function in between.

    Real, Intent(in) :: begin_fall      ! Begin the cosine function at this value
    Real, Intent(in) :: end_fall        ! End the cosine function at this value ( > begin_fall)
    Real, Intent(in) :: val(:)          ! Compute the cosine fall value at these values
    Real, Optional, Intent(in) :: base  ! Falls to 'base' >=0 instead of to 0.

    Real :: Cos_fallVR(Size(val))
! Local:
    Real :: tmp(Size(val))

    If (Present(base)) then
      Where (val <= begin_fall) 
        Cos_fallVR= 1
      Elsewhere (val >= end_fall) 
        Cos_fallVR= base
      Elsewhere
        tmp= (val - begin_fall) * PI / (end_fall - begin_fall)
        Cos_fallVR= base + (1 - base) * ((1 + Cos(tmp)) / 2)
      Endwhere
    Else
      Where (val <= begin_fall) 
        Cos_fallVR= 1
      Elsewhere (val >= end_fall) 
        Cos_fallVR= 0
      Elsewhere
        tmp= (val - begin_fall) * PI / (end_fall - begin_fall)
        Cos_fallVR= (1 + Cos(tmp)) / 2
      Endwhere
    End if
  End Function Cos_fallVR

 
  Subroutine Close_ordered1 (Ascending, delta, ord, next)

!   This subroutine identifies subsets of close values in the ordered
!   arrays 'ord'. Each element of a subset is defined as close to the next one
!   if its value is within 'delta' of the next one. Thus breaks of 'delta' or more
!   are identified - by 'next(i,j)' pointing to the next index of array 'j'
!   that differs from ord(i,j) by more than 'delta'.

    Logical,  Intent(in) :: Ascending  ! Ascending or descending ordered arrays in 'ord'
    Real,     Intent(in) :: delta      ! Closeness criterion
    Real,     Intent(in) :: ord(:)     ! (n) Set of ordered arrays 
    Integer, Intent(out) :: next(:)    ! (n) next(j) points to the next non-close index
                                       !       or n + 1 if none
! Local:
    Integer :: j, j1, n

    n= Size(ord);  next= n + 1

    If (Ascending) then
      Do j= n-1,1,-1
        j1= j + 1
        If (ord(j1) < ord(j) + delta) then
          next(j)= next(j1)
        Else
          next(j)= j1
        End if
      End do
    Else
      Do j= n-1,1,-1
        j1= j + 1
        If (ord(j1) > ord(j) - delta) then
          next(j)= next(j1)
        Else
          next(j)= j1
        End if
      End do
    End if

  End Subroutine Close_ordered1


  Subroutine Close_ordered2 (Ascending, delta, ord, next)

!   This subroutine identifies subsets of close values in the ordered
!   arrays 'ord'. Each element of a subset is defined as close to the next one
!   if its value is within 'delta' of the next one. Thus breaks of 'delta' or more
!   are identified - by 'next(i,j)' pointing to the next index of array 'j'
!   that differs from ord(i,j) by more than 'delta'.

    Logical,  Intent(in) :: Ascending  ! Ascending or descending ordered arrays in 'ord'
    Real,     Intent(in) :: delta      ! Closeness criterion
    Real,     Intent(in) :: ord(:,:)   ! (n,m) Set of ordered arrays 
    Integer, Intent(out) :: next(:,:)  ! (n,m) next(j,:) points to the next non-close index
                                       !       or n + 1 if none
! Local:
    Integer :: i, j, j1, n, m

    n= Size(ord,1);  m= Size(ord,2);  next= n + 1

    Do i= 1,m
      If (Ascending) then
        Do j= n-1,1,-1
          j1= j + 1
          If (ord(j1,i) < ord(j,i) + delta) then
            next(j,i)= next(j1,i)
          Else
            next(j,i)= j1
          End if
        End do
      Else
        Do j= n-1,1,-1
          j1= j + 1
          If (ord(j1,i) > ord(j,i) - delta) then
            next(j,i)= next(j1,i)
          Else
            next(j,i)= j1
          End if
        End do
      End if
    End do

  End Subroutine Close_ordered2

  
  Pure Function Cos_cubeVR (x1,x2, x)

!   This function rises smoothly from the constant value 'x1' 
!   for x <= x1 to become the straight line f(x)= x for x >= x2 > x1
!   as in Cubic_rise, but rises back to x2 via cosine function
!   as x goes to 0.

    Real, Intent(in) :: x1, x2
    Real, Intent(in) :: x(:)

    Real :: Cos_cubeVR(Size(x))  ! Output
! Local:
    Real :: d, f, dx(Size(x))

    d= x2 - x1;  f= d / 2
    Where (x >= x2)
      Cos_cubeVR= x
    Elsewhere (x <= x1)
      Cos_cubeVR= x1 + f * (1 + Cos(x * PI / x1))
    Elsewhere
      dx= (x - x1) / d
      Cos_cubeVR= x1 + (dx)**2 * (2 - dx) * d
    End where
  End Function Cos_cubeVR

  Pure Function Cos_cubeVD (x1,x2, x)

    Real,          Intent(in) :: x1, x2
    Real(Dblp), Intent(in) :: x(:)

    Real(Dblp) :: Cos_cubeVD(Size(x))  ! Output
! Local:
    Real(Dblp) :: d, f, dx(Size(x))

    d= x2 - x1;  f= d / 2
    Where (x >= x2)
      Cos_cubeVD= x
    Elsewhere (x <= x1)
      Cos_cubeVD= x1 + f * (1 + Cos(x * PI / x1))
    Elsewhere
      dx= (x - x1) / d
      Cos_cubeVD= x1 + (dx)**2 * (2 - dx) * d
    End where
  End Function Cos_cubeVD
  
    Subroutine SymMat_indices (To1D, To2D)
    
!     Compute index conversions for a symmetric matrix,
!     between a 2D matrix representation and a 1D 
!     representation of the upper triangle of the matrix.
    
      Integer,           Intent(out) :: To1D(:,:)  ! (m,m)
      Integer, Optional, Intent(out) :: To2D(:,:)  ! (ln,2)  for ln= (m-1)m/2
!   Local:
      Integer :: i, j, m, n
      
      m= Size(To1D,1);  n= 0
      If (Present(To2D)) then
        Do i= 1,m
          To1D(i,i)= 0
          Do j= i+1,m
            n= n + 1;  To1D(i,j)= n;  To1D(j,i)= n
            To2D(n,1)= i;  To2D(n,2)= j
          End do
        End do
      Else
        Do i= 1,m
          To1D(i,i)= 0
          Do j= i+1,m
            n= n + 1;  To1D(i,j)= n;  To1D(j,i)= n
          End do
        End do
      End if
    End Subroutine SymMat_indices    
    
    
    Pure Subroutine Real_vec_to_mat (vec, mat)
    
!     Compute index conversions for a symmetric matrix,
!     between a 2D matrix representation and a 1D 
!     representation of the upper triangle of the matrix.
    
      Real,  Intent(in) :: vec(:)    ! (m^2)
      Real, Intent(out) :: mat(:,:)  ! (m,m)
!   Local:
      Integer :: i, m, n, p
      
      m= Size(mat,1);  n= 1
      Do i= 1,m
        p= n + m-i;  mat(i,:)= vec(n:p);  n= p + 1
      End do
    End Subroutine Real_vec_to_mat    

    
    Subroutine Pearson_chi_sqr (p,deg, obs,mean, chi_sqr)
    
!     Pearson chi square goodness of fit test for the p-value 'p' = 
!     level of significance and # degrees of freedom 'deg'
!     for observation vector 'obs' and expected value vector 'mean' 
!     under a null hypothesis.

!     For p = 1, acceptance of the null hypothesis has a 0.1 error rate, 
!     for p = 2, it is a 0.05 error rate, for p = 3, it's 0.01.

      Integer,  Intent(in) :: p       ! Significance level
      Integer,  Intent(in) :: deg     ! Degrees of freedom (1..25)
      Real,     Intent(in) :: obs(:)  ! (n) Observation vector
      Real,     Intent(in) :: mean(:) ! (m) Expected values under the null hypothesis:
                                      !   If m = 1, the expected value is the same for
                                      !   all observations, else m = n
      
      Real,    Intent(out) :: chi_sqr ! The chi square statistic normalized by the Q table values
                                      !   so that the test is passed if chi_sqr < 1

      Real, Parameter :: mp= 3, mdg= 25
      Real :: Q_table1(mdg) = (/1.71,4.61,6.25,7.78,9.24,10.64,12.02,13.36,14.68,15.99,17.28,  &
                                18.55,19.81,21.06,22.31,23.54,24.77,25.99,27.20,28.41,29.62,30.81,32.01,33.20,34.38/)
      Real :: Q_table2(mdg) = (/3.84,5.99,7.81,9.49,11.07,12.59,14.07,15.51,16.92,18.31,19.68, &
                                21.03,22.36,23.68,25.00,26.30,27.59,28.87,30.14,31.41,32.67,   & 
                                33.92,35.17,36.41,37.65/)
      Real :: Q_table3(mdg) = (/6.63,9.21,11.34,13.28,15.09,16.81,18.48,20.09,21.67,23.21,24.73, &
                                26.22,27.69,29.14,30.58,32.00,33.41,34.81,36.19,37.57,38.93,40.29,41.64,42.98,44.31/)
      
      If (Size(obs) <= 1) then
        chi_sqr= 0;  Return
      Else if (p < 1 .or. p > 3) then
        chi_sqr= 10;  Return
      Else if (deg < 1 .or. deg > 25) then
        chi_sqr= 100;  Return
      End if
      
      If (Size(mean) == 1) then
        chi_sqr= Sum((obs - mean(1))**2) / mean(1)
      Else
        chi_sqr= Sum((obs - mean)**2 / mean)
      End if
      
      If (p == 1) then
        chi_sqr= chi_sqr / Q_table1(deg)
      Else if (p == 2) then
        chi_sqr= chi_sqr / Q_table2(deg)
      Else
        chi_sqr= chi_sqr / Q_table3(deg)
      End if
      
    End Subroutine Pearson_chi_sqr
    
    
    Pure Real Function Percentile(frac, dat)
    
!     Compute a percentile of an ordered real array (increasing or decreasing)
    
      Real,    Intent(in) :: frac    ! The fraction corresponding to the percentile
      Real,    Intent(in) :: dat(0:) ! (0:n) The ordered array
!   Local:
      Real    :: d, x
      Integer :: m, n
      
      n= Ubound(dat,1);  x= frac * n;  m= Nint(x)
      
      If (n == Nint(m / frac)) then
        Percentile= dat(m)
      Else 
        m= x;  d= x - m 
        Percentile= (1-d) * dat(m) + d * dat(m+1)
      End if
    End Function Percentile
    
   Pure Subroutine Recover_sym_mat (Add, Sym, Mat, lst)
  
     Logical, Intent(in) :: Add      ! True  : Add 'Sym' to the 'lst' part of 'Mat'
                                     ! False : Replace the 'lst' part of 'Mat' by 'Sym'
     Real,    Intent(in) :: Sym(:)   ! (l(l+1)/2) Matrix stored in symmetric form
    
     Real, Intent(inout) :: Mat(:,:) ! (n,n) Matrix stored in full form
     Integer, Optional, Intent(in) :: lst(:) ! (l) List of indices l <= n, in any order, to update 
                                             !     in 'Mat', else update the first 'l' indices
!  Local:
     Integer :: i, j, k, l, m   
    
     j= 1
     If (Present(lst)) then
       l= Size(lst)
       If (Add) then
         j= 1;  Do i= 1,l
           m= j + l-i;    k= lst(i)
           Mat(k,lst(i:))  = Mat(k,lst(i:)) + Sym(j:m);  j= m + 1
           Mat(lst(i+1:),k)= Mat(k,lst(i+1:))
         End do
       Else
         j= 1;  Do i= 1,l
           m= j + l-i;    k= lst(i);  Mat(k,lst(i:))= Sym(j:m);  j= m + 1
           Mat(lst(i+1:),k)= Mat(k,lst(i+1:))
         End do
       End if

     Else
       l= Int(Sqrt(Real(2*Size(Sym))))
       If (Add) then
         j= 1;  Do i= 1,l
           m= j + l-i;  Mat(i,i:l)= Mat(i,i:l) + Sym(j:m);  j= m + 1
           Mat(i+1:l,i)= Mat(i,i+1:l)
         End do
       Else
         j= 1;  Do i= 1,l
           m= j + l-i;  Mat(i,i:l)= Sym(j:m);  j= m + 1
           Mat(i+1:l,i)= Mat(i,i+1:l)
         End do
       End if
     End if
   End Subroutine Recover_sym_mat

   Pure Function Recover_sym_matD (n, Sym)  Result(Mat)
  
!    Recover a matrix 'Mat' whose upper triangle is stored in the
!    linear array 'Sym'

     Integer, Intent(in) :: n        ! Array size
     Real,    Intent(in) :: Sym(:)   ! (n(n+1)/2) Matrix stored in symmetric form
    
     Real(Dblp)       :: Mat(n,n) ! (n,n) Matrix stored in full form
!  Local:
     Integer :: i, j, m   
    
     j= 1;  Do i= 1,n
       m= j + n-i;  Mat(i,i:)= Sym(j:m);  j= m + 1
       Mat(i+1:,i)= Mat(i,i+1:)
     End do
   End Function Recover_sym_matD


   Pure Subroutine Store_sym_mat (Mat, Sym)
  
     Real,  Intent(in) :: Mat(:,:)  ! (n,n) Matrix stored in full form
     Real, Intent(out) :: Sym(:)    ! (n(n+1)/2) Matrix stored in symmetric form
!  Local:
     Integer :: i, j, m, n   
    
     n= Size(Mat,1);  j= 1
    
     Do i= 1,n
       m= j + n-i;  Sym(j:m)= Mat(i,i:);  j= m + 1
     End do
   End Subroutine Store_sym_mat

   Pure Function Square_sym (n, Sym)
  
!    Square a real symmetric matrix
  
     Integer, Intent(in) :: n               ! Array size
     Real,    Intent(in) :: Sym(:,:)        ! (n,n) Symmetric matrix
    
     Real                :: Square_sym(n,n) ! (n,n) Squared symmetric matrix
!  Local:
     Integer :: i, j  
    
     Do i= 1,n
       Do j= i,n
         Square_sym(i,j)= Dot_product(Sym(:,i),Sym(:,j))
       End do
       Square_sym(i+1:,i)= Square_sym(i,i+1:)
     End do
   End Function Square_sym


   Pure Function Eig_power (n,pwr, Eigval,Eigvec)
  
!    Compute a fractional power 'pwr' of a symmetric matrix
!    via eigenvalues and eigenvectors
  
     Integer,       Intent(in) :: n              ! Array size
     Real,          Intent(in) :: pwr            ! Power of the matrix to be computed
     Real(Dblp), Intent(in) :: Eigval(:)      ! (n)   Eigenvalues of the matrix
     Real(Dblp), Intent(in) :: Eigvec(:,:)    ! (n,n) Eigenvectors of the matrix
    
     Real                      :: Eig_power(n,n) ! (n,n) Symmetric matrix to be computed
!  Local:
     Real(Dblp) :: val_pwr(n), tmp(n)
     Integer       :: i, j
    
     val_pwr= Eigval**pwr
    
     Do i= 1,n
       tmp= val_pwr * Eigvec(i,:)
       Do j= i,n
         Eig_power(i,j)= Dot_product(tmp,Eigvec(j,:))
       End do
       Eig_power(i+1:,i)= Eig_power(i,i+1:)
     End do
   End Function Eig_power
  
   Pure Function Diag_to_vecR(n, Mat)
  
!    Diagonal of a real matrix
  
     Integer, Intent(in) :: n         ! Array size
     Real,    Intent(in) :: Mat(:,:)  ! (n,n) The matrix
    
     Real                :: Diag_to_vecR(n)  ! Its diagonal
!  Local:
     Integer :: i  
     Forall(i=1:n) Diag_to_vecR(i)= Mat(i,i)
   End Function Diag_to_vecR

   Pure Function Diag_to_vecI(n, Mat)
  
!    Diagonal of an integer matrix
  
     Integer, Intent(in) :: n         ! Array size
     Integer, Intent(in) :: Mat(:,:)  ! (n,n) The matrix
    
     Integer             :: Diag_to_vecI(n)  ! Its diagonal
!  Local:
     Integer :: i  
     Forall(i=1:n) Diag_to_vecI(i)= Mat(i,i)
   End Function Diag_to_vecI

   Pure Subroutine Vec_to_diagR (vec, Mat)
  
!    Assign the diagonal of a real matrix to a vector
  
     Real,    Intent(in) :: vec(:)    ! (n)   The vector
     Real, Intent(inout) :: Mat(:,:)  ! (n,n) The matrix
!  Local:
     Integer :: i  
     Forall(i=1:Size(Mat,1)) Mat(i,i)= vec(i)
   End Subroutine Vec_to_diagR

   Pure Subroutine Vec_to_diagI (vec, Mat)
  
!    Assign the diagonal of a Integer matrix to a vector
  
     Integer,    Intent(in) :: vec(:)    ! (n)   The vector
     Integer, Intent(inout) :: Mat(:,:)  ! (n,n) The matrix
!  Local:
     Integer :: i  
     Forall(i=1:Size(Mat,1)) Mat(i,i)= vec(i)
   End Subroutine Vec_to_diagI
  
   Pure Subroutine Const_to_diagI (const, Mat)
  
!    Assign the diagonal of a integer matrix to a constant
  
     Integer,    Intent(in) :: const     ! The constant
     Integer, Intent(inout) :: Mat(:,:)  ! (n,n) The matrix
!  Local:
     Integer :: i  
     Forall(i=1:Size(Mat,1)) Mat(i,i)= const
   End Subroutine Const_to_diagI

   Pure Subroutine Const_to_diagR (const, Mat)
  
!    Assign the diagonal of a real matrix to a constant
  
     Real,    Intent(in) :: const     ! The constant
     Real, Intent(inout) :: Mat(:,:)  ! (n,n) The matrix
!  Local:
     Integer :: i  
     Forall(i=1:Size(Mat,1)) Mat(i,i)= const
   End Subroutine Const_to_diagR

   Pure Subroutine Const_to_diagD (const, Mat)
  
!    Assign the diagonal of a real matrix to a constant
  
     Real(Dblp),    Intent(in) :: const     ! The constant
     Real(Dblp), Intent(inout) :: Mat(:,:)  ! (n,n) The matrix
!  Local:
     Integer :: i  
     Forall(i=1:Size(Mat,1)) Mat(i,i)= const
   End Subroutine Const_to_diagD

   
   Integer Function Quadrant(Diag, x,y)
     Logical, Intent(in) :: Diag  ! 
     Real,    Intent(in) :: x     ! 
     Real,    Intent(in) :: y     ! 

     If (Diag) then
       If (y >= 0 .and. x > 0) then       ! 1st quadrant
         Quadrant= 1
       Else if (y > 0 .and. x <= 0) then  ! 2nd quadrant
         Quadrant= 2
       Else if (y <= 0 .and. x < 0) then  ! 3rd quadrant
         Quadrant= 3
       Else                               ! 4th quadrant
         Quadrant= 4
       End if
     Else
       If (Abs(y) <= x) then        ! Right wedge
         Quadrant= 1
       Else if (y > Abs(x)) then    ! Top wedge
         Quadrant= 2
       Else if (Abs(y) <= -x) then  ! Left wedge
         Quadrant= 3
       Else                         ! Bottom wedge
         Quadrant= 4
       End if
     End if
   End Function Quadrant
   
  Pure Integer Function Locate_set (given_set, sets)

!   Locate a 'given_set' in a list of 'sets' and output
!   its index. If not found, the default output is 0.

    Integer,  Intent(in) :: given_set(:) ! (m) The set to be located in 'sets'
    Integer,  Intent(in) :: sets(:,:)    ! (m,n) The 'n' sets of size 'm', each assumed
                                         !       to be in increasing order
    Integer :: n;  n= Size(sets,2)

    Do Locate_set= 1,n
      If (All(given_set == sets(:,Locate_set))) Exit
    End do
    If (Locate_set > n) Locate_set= 0
  End Function Locate_set

  Elemental Logical Function Even (n)
    Integer, Intent(in) :: n
    Even= 2*(n/2) == n
  End Function Even

  Elemental Logical Function Div (n,d)
    Integer, Intent(in) :: n  ! 
    Integer, Intent(in) :: d  ! Possible divisor of 'n', d > 0
    Div= d*(n/d) == n
  End Function Div

  Function Shuffle_permutation (n)
    Integer, Intent(in) :: n
!  Result:
    Integer :: Shuffle_permutation(n)
!  Local:    
    Real    :: uniform(n-1)
    Integer :: i, j, k
    
    Call Random_number (uniform)
    Shuffle_permutation= "ID"
    
    Do i= 1,n-1
      k= i - 1;  j= k + Ceiling(uniform(i) * (n-k))
      j= Max(j,i)
      Call Swap (Shuffle_permutation(i),Shuffle_permutation(j))
    End do
  End Function Shuffle_permutation
    
    
  Subroutine Move_up (sub_list, list)
  
    Integer,    Intent(in) :: sub_list(:)  ! (k < n) of indicies 1:n
    Integer, Intent(inout) :: list(:)      ! (n)
    

    Logical :: Msk(Size(list))
    Integer :: k, m, n, ls(Size(list))
    
    k= Size(sub_list);  n= size(list)
    
    If (k > n) then
      Stop
    Else if (k < n .and. k > 0) then
      ls(1:k)= sub_list
      Msk= .true.;  Msk(sub_list)= .false.
      Call List_of_true (Msk, m,ls(k+1:))
      list= list(ls)
    End if
    
  End Subroutine Move_up
    

  Elemental Integer Function Partition_bound (n, k)
    
!   Bound on the number of partitions of the integer n into at exactly k parts. 
!   Ref:  http://link.springer.com/article/10.1007/s11139-007-9022-z#page-1 (by W. Pribitkin, 2006)

    Integer, Intent(in) :: n, k

    Real, Parameter :: c= 2.565, c2= 2 / c,  Li1= 1.645
    Real :: df
  
    df= Sqrt(Real(n - k))
    Partition_bound= Ceiling((exp(c * df) / df**1.5) * exp(-c2 * df) * Li1)
  
  End Function Partition_bound

  Subroutine Ranking_pt0 (parm,Max_pt,nc,mr, rank_pt,unrank_pt)

!   Computes ranking points from ballot ranking levels. Assume that the 
!   ranking levels are from 1 (ranked first) to mr (ranked last). Then
!   the points assigned to each ranking level decline as the level increases.
    
    Real,    Intent(in) :: parm(:) ! (2) ! Parameters for converting ranking levels to points
                                   ! i = 1: Scale factor for the points of unranked candidates, 
                                   !   0 <= p <= 1, with 1 = "full Borda" = all ballots get
                                   !     the same # total points, with more going to the 
                                   !     unranked when fewer are ranked.
                                   ! i = 2: Factor which determines the points.
                                   !   0 <= p <= 1: Increase the increment in an arithmetic
                                   !     increase in points from the last level to the first.
                                   !     The initial increment is 1, corresponding to the  
                                   !     Borda Count = 1, 2, 3,...,mr for p = 0. More generally, 
                                   !     the point values are 1,2+p,3+3p,4+6p,..,
                                   !     mr + ((mr-1)mr/2)p, with = 0.5 typical.
                                   !   p > 1 signifies a geometric increase in point values, 
                                   !     so that p = 2 means doubling: 1,2,4,8,...2**(mr-1).
    Real,    Intent(in) :: Max_pt  ! Max possible point value
    Integer, Intent(in) :: nc      ! Total # candidates
    Integer, Intent(in) :: mr      ! Max # rankable candidates
    
    Real,   Intent(out) :: rank_pt(:)   ! (mr) Point values corresponding to rising ranking levels
    Real,   Intent(out) :: unrank_pt(:) ! (mr) Point values for the unranked candidates, 
                                        !      depending on # ranked
! Local:
    Real    :: fac, spread, inc, rise, sm(mr)
    Integer :: i, n
    
!   Compute ranking points

    rank_pt= 0;  unrank_pt= 0

    If (parm(2) >= 0 .and. parm(2) <= 1) then  ! Linear point descent
      rank_pt(mr)= 1;  inc= 1
      Do i= mr-1,1,-1
        inc= inc + parm(2)
        rank_pt(i)= rank_pt(i+1) + inc
      End do
      
      fac= Max_pt / rank_pt(1)
      rank_pt= fac * rank_pt
      
    Else if (parm(2) > 1) then  ! Geometric point descent
      rank_pt(1)= Max_pt
      Do i= 1,mr-1
        rank_pt(i+1)= rank_pt(i) / parm(2)
      End do
    End if
    
!   Compute points for the unranked candidates
    
    sm= Partial_sums(rank_pt)
    
    Do i= 1,mr-1
      unrank_pt(i)= (sm(mr) - sm(i)) / (nc - i)
    End do
    unrank_pt(mr)= 0

    If (parm(1) >= 0 .and. parm(1) < 1) then
      unrank_pt= parm(1) * unrank_pt
    End if
    
    If (pr_out >= 1) then
      Call Out ("Ranked points",rank_pt)
      Call Out ("Unranked points",unrank_pt)
    End if
  End Subroutine Ranking_pt0

  Subroutine Ranking_pts (parm,Max_pt,nc,mr, rank_pt,unrank_pt, pt_val)

!   Computes ranking points from ballot ranking levels. Assume that the 
!   ranking levels are from 1 (ranked first) to mr (ranked last). Then
!   the points assigned to each ranking level decline as the level increases.
    
    Real,    Intent(in) :: parm(:) ! (4) ! Parameters for converting ranking levels to points
                                   ! i = 1: Scale factor for the points of unranked candidates, 
                                   !   0 <= p <= 1, with 1 = "full Borda" = all ballots get
                                   !     the same # total points, with more going to the 
                                   !     unranked when fewer are ranked.
                                   ! i = 2: Factor which determines the points.
                                   !   0 <= p <= 1: Increase the increment in an arithmetic
                                   !     increase in points from the last level to the first.
                                   !     The initial increment is 1, corresponding to the  
                                   !     Borda Count = 1, 2, 3,...,mr for p = 0. More generally, 
                                   !     the point values are 1,2+p,3+3p,4+6p,..,
                                   !     mr + ((mr-1)mr/2)p, with = 0.5 typical.
                                   !   p > 1 signifies a geometric increase in point values, 
                                   !     so that p = 2 means doubling: 1,2,4,8,...2**(mr-1).
                                   ! i = 3: The sigma for ranking level 1 as a random 
                                   !     variable, assumed to be the most certain, hence the
                                   !     smallest sigma, over all the ranking levels. p = 1/4
                                   !     is typical. The sigma of the corresponding point 
                                   !     values is the level sigm scaled by the point 
                                   !     spread, which is the average of the gaps 
                                   !     to the adjacent point values,the level spread = 1.
                                   ! i = 4: The sigma for the unranked candidates as random 
                                   !     variables, assumed to be the most uncertain, hence the
                                   !     largest sigma, over all the ranking levels, with 
                                   !     p = 1  typical. A quadratic increase is assumed from
                                   !     level 1 through level mr+1, representing the unranked
                                   !     candidates.
                                   !        However, for rating data, it is the first and last
                                   !     levels (= most positive and most negative points)
                                   !     most certain, with the unrated candidates getting
                                   !     0 points, with the most uncertainty.
    Real,    Intent(in) :: Max_pt  ! Max possible point value
    Integer, Intent(in) :: nc      ! Total # candidates
    Integer, Intent(in) :: mr      ! Max # rankable candidates
    
    Real,   Intent(out) :: rank_pt(:,:)   ! (mr,2) Point values (:,1) and their variances (:,2)
                                          !        corresponding to rising ranking levels
    Real,   Intent(out) :: unrank_pt(:,:) ! (mr,2) Point values and variances for the 
                                          !        unranked candidates, depending on # ranked
    
    Real, Optional, Intent(in) :: pt_val(:) ! (mr) Specified point values for the ranking levels
                                            !      instead of those determined by 'parm'
! Local:
    Real    :: sm(mr), sqrt_pr(2), sig_lev(mr+1), sig_pt(mr+1)
    Real    :: fac, spread, inc, rise
    Integer :: i, n
    
!   Compute ranking points

    rank_pt= 0;  unrank_pt= 0

    If (Present(pt_val)) then
      rank_pt(:,1)= pt_val
    Else if (parm(2) <= 1 .and. parm(2) >= 0) then  ! Linear point descent
      rank_pt(mr,1)= 1;  inc= 1
      Do i= mr-1,1,-1
        inc= inc + parm(2)
        rank_pt(i,1)= rank_pt(i+1,1) + inc
      End do
      
      fac= Max_pt / rank_pt(1,1)
      rank_pt(:,1)= fac * rank_pt(:,1)
      
    Else if (parm(2) > 1) then  ! Geometric point descent
      rank_pt(1,1)= Max_pt
      Do i= 1,mr-1
        rank_pt(i+1,1)= rank_pt(i,1) / parm(2)
      End do
    End if
    
!   Compute points (= expected value) for the unranked candidates
    
    sm= Partial_sums(rank_pt(:,1))
    
    Do i= 1,mr-1
      unrank_pt(i,1)= (sm(mr) - sm(i)) / (nc - i)
    End do
    unrank_pt(mr,1)= 0

    If (parm(1) >= 0 .and. parm(1) < 1) then
      unrank_pt(:,1)= parm(1) * unrank_pt(:,1)
    End if
    
    If (parm(3) < 0.01 .or. parm(4) > 50.0 .or. parm(4) <= parm(3)) then ! Variance not valid
      If (pr_out >= 1) then
        Call Out ("Ranked points",rank_pt(:,1))
        Call Out ("Unranked points",unrank_pt(:,1))
      End if
      rank_pt(:,2)= 0;  unrank_pt(:,2)= -1;  Return
    End if
      
!   Variance values
    
    sqrt_pr= Sqrt(parm(3:4));  rise= sqrt_pr(2) - sqrt_pr(1)
    Forall(i=1:mr+1) sig_lev(i)= sqrt_pr(1) + rise*Real(i-1)/mr
    sig_lev= sig_lev**2

    Do i= 1,mr+1
      If (i > 1 .and. i < mr) then
        spread= (rank_pt(i-1,1) - rank_pt(i+1,1)) / 2.0
      Else if (i == 1) then
        spread= rank_pt(1,1) - rank_pt(2,1)
      Else if (i == mr) then
        spread= rank_pt(mr-1,1) / 2.0
      Else if (i == mr+1) then
        spread= rank_pt(mr,1)
      End if

      sig_pt(i)= spread * sig_lev(i)
    End do

    rank_pt(:mr,2)= sig_pt(:mr)**2

!   Average variances of the remaining unranked
!   for partially ranked ballots to estimate the variance 
!   for those left unranked by a partially ranked ballot

    Do i= 1,mr-1
      unrank_pt(i,2)= (Sum(rank_pt(i+1:mr,2)) + &
                      (nc-mr) * sig_pt(mr+1)**2) / (nc - i)
    End do
    unrank_pt(mr,2)= sig_pt(mr+1)**2

    If (pr_out >= 1) then
      Call Out (-1,"Ranked points & variances",rank_pt)
      Call Out (-1,"Unranked points & variances",unrank_pt)

      sm= rank_pt(:,1)/Sqrt(rank_pt(:,2))
      Call Out ("Ratio of ranked points to ranked sigmas",sm)
      sm= unrank_pt(:,1)/Sqrt(unrank_pt(:,2))
      Call Out ("Ratio of unranked points to unranked sigmas",sm)
    End if
    
  End Subroutine Ranking_pts

  Subroutine Rating_pts (parm,mt,pt_val, rt_var)
   
!   Compute the variances 'rt_var' of the 'pt_val' rating point values
    
    Real,    Intent(in) :: parm(:) ! (4) ! Parameters for converting ranking levels to points
                                   ! i = 1: Scale factor for the points of unranked candidates, 
                                   !   0 <= p <= 1, with 1 = "full Borda" = all ballots get
                                   !     the same # total points, with more going to the 
                                   !     unranked when fewer are ranked.
                                   ! i = 2: Factor which determines the points.
                                   !   0 <= p <= 1: Increase the increment in an arithmetic
                                   !     increase in points from the last level to the first.
                                   !     The initial increment is 1, corresponding to the  
                                   !     Borda Count = 1, 2, 3,...,mr for p = 0. More generally, 
                                   !     the point values are 1,2+p,3+3p,4+6p,..,
                                   !     mr + ((mr-1)mr/2)p, with = 0.5 typical.
                                   !   p > 1 signifies a geometric increase in point values, 
                                   !     so that p = 2 means doubling: 1,2,4,8,...2**(mr-1).
                                   ! i = 3: The sigma for ranking level 1 as a random 
                                   !     variable, assumed to be the most certain, hence the
                                   !     smallest sigma, over all the ranking levels. p = 1/4
                                   !     is typical. The sigma of the corresponding point 
                                   !     values is the level sigm scaled by the point 
                                   !     spread, which is the average of the gaps 
                                   !     to the adjacent point values,the level spread = 1.
                                   ! i = 4: The sigma for the unranked candidates as random 
                                   !     variables, assumed to be the most uncertain, hence the
                                   !     largest sigma, over all the ranking levels, with 
                                   !     p = 1  typical. A quadratic increase is assumed from
                                   !     level 1 through level mr+1, representing the unranked
                                   !     candidates.
                                   !        However, for rating data, it is the first and last
                                   !     levels (= most positive and most negative points)
                                   !     most certain, with the unrated candidates getting
                                   !     0 points, with the most uncertainty.
    Integer, Intent(in) :: mt         ! # rating levels
    Real,    Intent(in) :: pt_val(:)  ! (mt) Decreasing point values, from Max_pt to -Max_pt,
                                      !      of the 'mt' rating levels
    Real,   Intent(out) :: rt_var(0:) ! (0:mt) The corresponding variances, with the variance
                                      !        of unrated candidates at 0.
! Local:
    Real    :: sig_lev(mt), sig_pt(mt), sqrt_pr(2)
    Real    :: spread, rise
    Integer :: i, j, mtp, mtn
    
    rt_var= 0;  mtp= Last_true(pt_val > 0);  mtn= mt - mtp

    spread= (pt_val(mtp) - pt_val(mtp+1)) / 2.0
    rt_var(0)= (parm(4) * spread)**2

    sqrt_pr= Sqrt(parm(3:4));  rise= sqrt_pr(2) - sqrt_pr(1)
    Forall(i=1:mtp) sig_lev(i)= sqrt_pr(1) + rise*Real(i-1)/mtp
    sig_lev(:mtp)= sig_lev(:mtp)**2


    rise= parm(4) - parm(3)
    Forall(i=1:mtp) sig_lev(i)= parm(3) + rise*Real(i-1)/mtp

    Do i= 1,mtp
      If (i > 1 .and. i < mtp) then
        spread= (pt_val(i-1) - pt_val(i+1)) / 2.0
      Else if (i == 1) then
        spread= pt_val(1) - pt_val(2)
      Else if (i == mtp) then
        spread= pt_val(mtp-1) / 2
      End if
      sig_pt(i)= spread * sig_lev(i)
    End do

    Forall(i=1:mtn) sig_lev(i)= sqrt_pr(1) + rise*Real(i-1)/mtn
    sig_lev(:mtn)= sig_lev(:mtn)**2

    Do j= mtp+1,mt
      If (j > mtp+1 .and. j < mt) then
        spread= (pt_val(j-1) - pt_val(j+1)) / 2.0
      Else if (j == mtp+1) then
        spread= -pt_val(j)
      Else if (j == mt) then
        spread= pt_val(mt-1) - pt_val(mt)
      End if
      sig_pt(j)= spread * sig_lev(mt+1-j)
    End do

    rt_var(1:)= sig_pt**2

    If (pr_out >= 1) then
      Call Out ("Rating point variances",rt_var(1:))
      Call Out ("Unrated point variance",rt_var(0))
      sig_pt= Sqrt(rt_var(1:))/sig_pt
      Call Out ("Ratio of ratings to rating sigmas",sig_pt)
    End if
  End Subroutine Rating_pts
     

  Subroutine Rating0_pts (parm,mt,pt_val, rt_var)
   
!   Compute the variances 'rt_var' of the 'pt_val' ranking point values
    
    Real,    Intent(in) :: parm(:) ! (4) ! Parameters for converting ranking levels to points
                                   ! i = 1: Scale factor for the points of unranked candidates, 
                                   !   0 <= p <= 1, with 1 = "full Borda" = all ballots get
                                   !     the same # total points, with more going to the 
                                   !     unranked when fewer are ranked.
                                   ! i = 2: Factor which determines the points.
                                   !   0 <= p <= 1: Increase the increment in an arithmetic
                                   !     increase in points from the last level to the first.
                                   !     The initial increment is 1, corresponding to the  
                                   !     Borda Count = 1, 2, 3,...,mr for p = 0. More generally, 
                                   !     the point values are 1,2+p,3+3p,4+6p,..,
                                   !     mr + ((mr-1)mr/2)p, with = 0.5 typical.
                                   !   p > 1 signifies a geometric increase in point values, 
                                   !     so that p = 2 means doubling: 1,2,4,8,...2**(mr-1).
                                   ! i = 3: The sigma for ranking level 1 as a random 
                                   !     variable, assumed to be the most certain, hence the
                                   !     smallest sigma, over all the ranking levels. p = 1/4
                                   !     is typical. The sigma of the corresponding point 
                                   !     values is the level sigm scaled by the point 
                                   !     spread, which is the average of the gaps 
                                   !     to the adjacent point values,the level spread = 1.
                                   ! i = 4: The sigma for the unranked candidates as random 
                                   !     variables, assumed to be the most uncertain, hence the
                                   !     largest sigma, over all the ranking levels, with 
                                   !     p = 1  typical. A quadratic increase is assumed from
                                   !     level 1 through level mr+1, representing the unranked
                                   !     candidates.
                                   !        However, for rating data, it is the first and last
                                   !     levels (= most positive and most negative points)
                                   !     most certain, with the unrated candidates getting
                                   !     0 points, with the most uncertainty.
    Integer, Intent(in) :: mt         ! # rating levels
    Real,    Intent(in) :: pt_val(:)  ! (mt) Decreasing positve point values, down from Max_pt
    Real,   Intent(out) :: rt_var(0:) ! (0:mt) The corresponding variances, with the variance
                                      !        of unrated candidates at 0.
! Local:
    Real    :: sig_lev(mt), sig_pt(mt), sqrt_pr(2)
    Real    :: spread, rise
    Integer :: i, j
    
    rt_var= 0;  spread= pt_val(mt)
    rt_var(0)= (parm(4) * spread)**2

    sqrt_pr= Sqrt(parm(3:4));  rise= sqrt_pr(2) - sqrt_pr(1)
    Forall(i=1:mt) sig_lev(i)= sqrt_pr(1) + rise*Real(i-1)/mt
    sig_lev(:mt)= sig_lev(:mt)**2

    rise= parm(4) - parm(3)
    Forall(i=1:mt) sig_lev(i)= parm(3) + rise*Real(i-1)/mt

    Do i= 1,mt
      If (i > 1 .and. i < mt) then
        spread= (pt_val(i-1) - pt_val(i+1)) / 2.0
      Else if (i == 1) then
        spread= pt_val(1) - pt_val(2)
      Else if (i == mt) then
        spread= pt_val(mt-1) / 2
      End if
      sig_pt(i)= spread * sig_lev(i)
    End do

    rt_var(1:)= sig_pt**2

    If (pr_out >= 1) then
      Call Out ("Rating point variances",rt_var(1:))
      Call Out ("Unrated point variance",rt_var(0))
      sig_pt= Sqrt(rt_var(1:))/sig_pt
      Call Out ("Ratio of ratings to rating sigmas",sig_pt)
    End if
  End Subroutine Rating0_pts

  Subroutine Rank_pts0 (parm,Max_pt, rank_pt, pt_val)
  
!   Ranking point values, as specified by 'parm' or using prior point values 

    Real,    Intent(in) :: parm               !  > 1: Scale factor for a geometric increase in points, 
                                              !          starting from 1 point for the lowest ranked candidate
                                              !  0 <= parm <= 1: The point increment increases by this value
                                              !          starting from 1 point for the lowest ranked candidate,
                                              !          so that 0 means standard Borda with constant point increment of 1
                                              !          and 1 means the point increment increases by 1 with each rank  
    Real,    Intent(in) :: Max_pt             ! Max possible point value
    
    Real,   Intent(out) :: rank_pt(:)         ! (mr) Point values for the ranking levels, using Max_pt

    Real, Optional,  Intent(in) :: pt_val(:)  ! (mr) Prior point values for the ranking levels 1..mr,
! Local:
    Real    :: fac, inc
    Integer :: i, n, mr
    
    mr= Size(rank_pt,1);  rank_pt= 0
    
!   Compute ranking points    
    
    If (Present(pt_val)) then
      rank_pt= pt_val
    Else if (parm <= 1 .and. parm >= 0) then  ! Linear point descent
      rank_pt(mr)= 1;  inc= 1
      Do i= mr-1,1,-1
        inc= inc + parm
        rank_pt(i)= rank_pt(i+1) + inc
      End do
      
      fac= Max_pt / rank_pt(1)
      rank_pt= fac * rank_pt
      
    Else if (parm > 1) then  ! Geometric point descent
      rank_pt(1)= Max_pt
      Do i= 1,mr-1
        rank_pt(i+1)= rank_pt(i) / parm
      End do
    End if
    
  End Subroutine Rank_pts0

   Pure Subroutine Vec_meanR (wt,Mat, mean)
  
!    Compute the mean vector 'mean' across the 1st dimension of the 
!    matrix 'Mat' using the weights 'wt' for the 2nd dimension
  
     Real,  Intent(in) :: wt(:)     ! (m)   The non-negative vector weights. Need not sum to 1.
     Real,  Intent(in) :: Mat(:,:)  ! (n,m) The matrix
     Real, Intent(out) :: mean(:)   ! (n)   The mean vector
!  Local:
     Real    :: sm
     Integer :: i, m
     
     m= Size(wt)
     
     If (m < 1) then
       mean= 0
     Else if (m == 1) then
       mean= Mat(:,1)
     Else
       mean= Matmul(Mat,wt)
       sm= Sum(wt);  mean= mean / sm
     End if
   End Subroutine Vec_meanR

   Pure Subroutine Vec_meanD (wt,Mat, mean)
  
!    Compute the mean vector 'mean' across the 1st dimension of the 
!    matrix 'Mat' using the weights 'wt' for the 2nd dimension
  
     Real(Dblp),  Intent(in) :: wt(:)     ! (m)   The non-negative vector weights. Need not sum to 1.
     Real(Dblp),  Intent(in) :: Mat(:,:)  ! (n,m) The matrix
     Real(Dblp), Intent(out) :: mean(:)   ! (n)   The mean vector
!  Local:
     Real(Dblp) :: sm
     Integer :: i, m
     
     m= Size(wt)
     
     If (m < 1) then
       mean= 0
     Else if (m == 1) then
       mean= Mat(:,1)
     Else
       mean= Matmul(Mat,wt)
       sm= Sum(wt);  mean= mean / sm
     End if
   End Subroutine Vec_meanD
    
   Subroutine Normalize_vec1_R (rate)
     Real, Intent(inout) :: rate(0:)  ! (0:nc)
     
     rate(0) = Sqrt(Sum(rate(1:)**2))
     rate(1:)= rate(1:) / rate(0)
   End Subroutine Normalize_vec1_R  
   
   Subroutine Normalize_vec1_D (rate)
     Real(Dblp), Intent(inout) :: rate(0:)  ! (0:nc)
     
     rate(0) = Sqrt(Sum(rate(1:)**2))
     rate(1:)= rate(1:) / rate(0)
   End Subroutine Normalize_vec1_D  
   
   Subroutine Normalize_vec2_R (rate)
     Real, Intent(inout) :: rate(0:,:)  ! (0:nc,ns)
     Integer :: i
     
     Do i= 1,Ubound(rate,2)
       rate(0,i) = Sqrt(Sum(rate(1:,i)**2))
       rate(1:,i)= rate(1:,i) / rate(0,i)
     End do
   End Subroutine Normalize_vec2_R  
   
   Subroutine Normalize_vec2_D (rate)
     Real(Dblp), Intent(inout) :: rate(0:,:)  ! (0:nc,ns)
     Integer :: i
     
     Do i= 1,Ubound(rate,2)
       rate(0,i) = Sqrt(Sum(rate(1:,i)**2))
       rate(1:,i)= rate(1:,i) / rate(0,i)
     End do
   End Subroutine Normalize_vec2_D  
   
   
   Subroutine Select_rndR (nr,nb,prob, pr_sl,Selected)
   
!    Select a random subset of elements which have a weights or 
!    a discrete probability distribution 'prob'
   
     Integer,  Intent(in) :: nr           ! Size of the random subset
     Integer,  Intent(in) :: nb           ! Size of the set or distribution, nb > nr
     Real,     Intent(in) :: prob(:)      ! (nb)  Weights or probabilities of the distribution in any order
     
     Real,    Intent(out) :: pr_sl(:)     ! (nr)  Weights of the selected elements by decreasing weight
     Integer, Intent(out) :: Selected(:)  ! (nr)  Corresponding list of the selected elements 
!  Local:
     Real    :: rnd(nr), cum_wt(nb)
     Logical :: Sel(nb)
     Integer :: orig(nb), key(nr)
     Integer :: b, j, k, n
     
     If (nb < nr) then
       Call Out ("Error in Select_rndR. nb",nb, "but size of random subset",nr); Stop
     End if

     Call Random_number (rnd);  cum_wt= Partial_sums(prob)
     Sel= .false.;  orig= "ID";  n= nb
     
     Do j= 1,nr
       k= First_true(rnd(j) <= cum_wt(:n)/cum_wt(n))
       b= orig(k);  Sel(b)= .true.
       
!      Remove 'b' from the list 'orig' of eligible elements:       
           
       orig(k:n-1)= orig(k+1:n);  cum_wt(k:n-1)= cum_wt(k+1:n) - prob(b)
       n= n - 1
     End do

     Call List_of_true (Sel, k,Selected)
     If (k /= nr) then
       Call Out ("Error in Select_rndR. # selected",k, "but # expected",nr); Stop
     End if

     pr_sl= prob(Selected)
     Call Sort (.false.,pr_sl,key)
     Selected= Selected(key)
     
   End Subroutine Select_rndR
   
   Subroutine Select_rndD (nr,nb,prob, pr_sl,Selected)
   
!    Select a random subset of elements which have a weights or 
!    a discrete probability distribution 'prob'
   
     Integer,        Intent(in) :: nr           ! Size of the random subset
     Integer,        Intent(in) :: nb           ! Size of the set or distribution, nb > nr
     Real(Dblp),  Intent(in) :: prob(:)      ! (nb)  Weights or probabilities of the distribution in any order
     
     Real(Dblp), Intent(out) :: pr_sl(:)     ! (nr)  Weights of the selected elements by decreasing weight
     Integer,       Intent(out) :: Selected(:)  ! (nr)  Corresponding list of the selected elements 
!  Local:
     Real(Dblp) :: rnd(nr), cum_wt(nb)
     Logical :: Sel(nb)
     Integer :: orig(nb), key(nr)
     Integer :: b, j, k, n
     
     If (nb < nr) then
       Call Out ("Error in Select_rndR. nb",nb, "but size of random subset",nr); Stop
     End if

     Call Random_number (rnd);  cum_wt= Partial_sums(prob)
     Sel= .false.;  orig= "ID";  n= nb
     
     Do j= 1,nr
       k= First_true(rnd(j) <= cum_wt(:n)/cum_wt(n))
       b= orig(k);  Sel(b)= .true.
       
!      Remove 'b' from the list 'orig' of eligible elements:       
           
       orig(k:n-1)= orig(k+1:n);  cum_wt(k:n-1)= cum_wt(k+1:n) - prob(b)
       n= n - 1
     End do

     Call List_of_true (Sel, k,Selected)
     If (k /= nr) then
       Call Out ("Error in Select_rndR. # selected",k, "but # expected",nr); Stop
     End if

     pr_sl= prob(Selected)
     Call Sort (.false.,pr_sl,key)
     Selected= Selected(key)
     
   End Subroutine Select_rndD

  End Module Util


