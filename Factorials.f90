
Module Factorials

!  Routines for computing factorials, permutations, subsets,
!  and combinations. This includes listing these by size 
!  and finding these in lists.

  Use Util
  Use Output
  Use Precisn
  Implicit None

! Compute factorials and # subsets

  Interface Factorial
    Module procedure Factorial_s, Factorial_v, Factorial_p
  End Interface Factorial

  Interface N_subsets
    Module procedure N_subsets_s, N_subsets_v
  End Interface N_subsets

Contains

  Pure Integer Function Factorial_s(n)
	Integer, Intent(in) :: n
    Integer :: i

    Factorial_s= 1;  If (n <= 1) Return
    Do i= 2,n
      Factorial_s= i * Factorial_s
    End do
  End Function Factorial_s

  Pure Function Factorial_v(n)
	Integer, Intent(in) :: n(:)
	Integer  :: Factorial_v(Size(n))
    Integer :: i, j

    Factorial_v= 1
    Do i= 1,Size(n)
      Do j= 2,n(i)
        Factorial_v(i)= j * Factorial_v(i)
      End do
    End do
  End Function Factorial_v

  Pure Integer Function Factorial_p(n,m)
  
!   # Permuations of 'n' things taken 'm' at a time

	Integer, Intent(in) :: n, m  ! 0 <= m <= n required
    Integer :: i, lim

    Factorial_p= 1;  If (m <= 0) Return
    lim= Max(n-m+1,1)
    
    Do i= n,lim,-1
      Factorial_p= i * Factorial_p
    End do
  End Function Factorial_p

! # subsets of size k of a set of size n =
! # combinations of n things taken k at a time

  Pure Integer Function N_subsets_s(n, k)
    Integer, Intent(in) :: n, k
! Local:
    Integer :: i
    Real    :: x, y

    If (n < 0 .or. k < 0 .or. k > n) then
      N_subsets_s= -1
    Else if (k == 0 .or. k == n) then
      N_subsets_s= 1
    Else
      y= n - k;  x= y + 1
      Do i= 2,k
        x= x * (y + i) / i
      End do
      N_subsets_s= Nint(x)
    End if
  End Function N_subsets_s

! Vector version of 'N_subsets_s': Requires that n(:), k(:), 
! and N_subsets_v(:) all have the same length, so that 
! N_subsets_v(i)= N_subsets_s(n(I), k(i)). Error if N_subsets_v(i) < 0

  Pure Function N_subsets_v(n, k)
    Integer, Intent(in) :: n(:), k(:)
    Integer :: N_subsets_v(Size(n))
! Local:
    Integer :: i, j, nl, nk
    Real    :: x, y

    nl= Size(n);  nk= Size(k)

    If (nk < nl .or. nl < 1) then
      N_subsets_v= -1
    Else
      Do j= 1,nl
        If (n(j) < 0 .or. k(j) < 0 .or. k(j) > n(j)) then
          N_subsets_v(j)= -1
        Else if (k(j) == 0 .or. k(j) == n(j)) then
          N_subsets_v(j)= 1
        Else
          y= n(j) - k(j);  x= y + 1
          Do i= 2,k(j)
            x= x * (y + i) / i
          End do
          N_subsets_v(j)= Nint(x)
        End if
      End do
    End if
  End Function N_subsets_v


  Pure Integer Function Sub_ind(n, sub)

!   Compute the index of a subset 'sub' of size 'm' in a depth-first 
!   listing of all subsets of size 'm' of the set {1,...,n}.
!   Assume that 'sub' is in increasing order.

!   The depth-first listing starts with all subsets that start with 1,
!   then with 2, then with 3, etc.

    Integer, Intent(in) :: n       ! Size of the main set of integers: {1,...,n}.
    Integer, Intent(in) :: sub(:)  ! (m) The subset of the main set in increasing order
! Local:
    Integer :: i, j, m, c0, c1

    m= Size(sub)
    
    If (m > n .or. m < 1) then
      Sub_ind= -1  ! Error
    Else if (sub(1) < 1 .or. sub(m) > n .or. &
             Any(sub(:m-1) >= sub(2:))) then  
      Sub_ind= -1  ! Error
    Else if (m == 1) then
      Sub_ind= sub(1)
    Else
      c0= 0;  Sub_ind= 1

      Do i= 1,m
        c1= sub(i)

        Do j= c0+1,c1-1
          Sub_ind= Sub_ind + N_subsets(n-j, m-i)
        End do
        
        c0= c1
      End do
    End if
  End Function Sub_ind
  
  
  Subroutine Ind_sub (n,ind, sub, err)

!   Compute the subset 'sub' of size 'm' corresponding to the index 'ind'
!   in a depth-first listing of all subsets of size 'm' of the set {1,...,n}.
!    'sub' will be in increasing order.

!   The depth-first listing starts with all subsets that start with 1,
!   then with 2, then with 3, etc.

    Integer,  Intent(in) :: n       ! Size of the main set of integers: {1,...,n}.
    Integer,  Intent(in) :: ind     ! Index of the subset.
    Integer, Intent(out) :: sub(:)  ! (m) The subset of the main set in increasing order
    Integer, Intent(out) :: err     ! Error indicator. OK if > 0
! Local:
    Integer :: c, i, j, k, m, mx, sb

    m= Size(sub);  sub= -1;  err= 1
    
    If (m > n .or. m < 1) then
      err= -1;  Return    ! Error
    End if
    
    mx= N_subsets(n, m)
    If (ind < 1 .or. ind > mx) then  
      err= -1;  Return    ! Error
    End if
    
    If (m == 1) then
      sub(1)= ind;  Return
    End if
  
    c= 1;  k= ind
    
    Outer : Do i= 1,m
      Do j= c,n-(m-i)
        sb= N_subsets(n-j, m-i)
        If (k <= sb) then
          sub(i)= j;  c= j+1;  Cycle Outer
        Else
          k= k - sb  
        End if
      End do
    End do Outer
    
  End Subroutine Ind_sub
  
  
  Pure Integer Function Subset_ind (n,p, sub, n_sub)

!   Compute the index of a subset "sub" in a depth-first, or lexicographic,
!   listing of all non-empty subsets of size <= "p" of the set {1,...,n}.
!   Assume that "sub" is in increasing order.

!   The depth-first listing starts with all subsets that start with 1,
!   then {1,2},{1,2,3}, etc.

    Integer, Intent(in) :: n       ! Size of the main set of integers: {1,...,n}.
    Integer, Intent(in) :: p       ! Limit on subset size, 1 <= p <= n
    Integer, Intent(in) :: sub(:)  ! (<=p) The subset integers, smallest to largest
    Integer, Optional, Intent(in) :: n_sub(:,0:) ! (n-1,0:p-1) n_sub(i,j)= # subsets 
                                   !   of size <= j of a set of size i
! Local:
    Integer :: i, j, k, l, m, c0, c1, pmi, ind

    m= Size(sub)
    
    If (p > n .or. m > p .or. m < 1) then
      Subset_ind= -1  ! Error
    Else if (p == 1 .and. p == n) then
      Subset_ind= 1   ! Only one subset, so "sub" must be it.
    Else
      c0= 0;  ind= m

      Do i= 1,m
        c1= sub(i);  pmi= p - i

        If (Present(n_sub)) then
          Do j= 1,c1 - c0 - 1
            k= n-c0-j;  l= Min(k,pmi)
            ind= ind + n_sub(k,l)
          End do
        Else
          Do j= 1,c1 - c0 - 1
            ind= ind + N_sub_sum(n - c0 - j, pmi)
          End do
        End if
        c0= c1
      End do

      Subset_ind= ind
    End if
  End Function Subset_ind


  Pure Integer Function N_sub_sum(n, p)

!   # subsets of sizes <= p of a set of size n >= p

	Integer, Intent(in) :: n, p
! Local:
    Integer :: i, j, k, l
    Real    :: x, y

    k= Min(n, p)

    If (k < 0) then
      N_sub_sum= 0
    Else if (k == 0) then
      N_sub_sum= 1
    Else if (k > n/2) then
      N_sub_sum= 2**n
      If (k < n) then
        N_sub_sum= N_sub_sum - 1
        If (k < n-1) then
          N_sub_sum= N_sub_sum - n;  j= n - k - 1
          Do l= 2,j
            y= n - l;  x= y + 1
            Do i= 2,l
              x= x * (y + i) / i
            End do
            N_sub_sum= N_sub_sum - Nint(x)
          End do
        End if
      End if
    Else 
      N_sub_sum= 1 + n
      Do l= 2,k
        y= n - l;  x= y + 1
        Do i= 2,l
          x= x * (y + i) / i
        End do
        N_sub_sum= N_sub_sum + Nint(x)
      End do
    End if
  End Function N_sub_sum


  Recursive Subroutine Permutations (n,l, set, perm)

!   List all l= n! permutations of the set "set" of length n
!   in the array "perm".

    Integer,  Intent(in) :: n
    Integer,  Intent(in) :: l           ! = n!
    Integer,  Intent(in) :: set(:)      ! (n)
    Integer, Intent(out) :: perm(:,:)   ! (n,l)
! Local:
    Integer :: i, j, k, jk, n1, st(n-1)
      
    j= 1;  k= l / n

    Do i= 1,n
      Call List_of_true (set /= set(i), n1, st)

      jk= j + k;  perm(1,j:jk-1)= set(i)

      Call Permutations (n1,k, st, perm(2:,j:jk-1));  j= jk
    End do
  End Subroutine Permutations


  Recursive Subroutine List_subsets (n1,n2,p,l, k,ls)

!   List all non-empty subsets of {n1,...,n2} of size <= p
!   as extensions of a prior set of size l. Use depth first ordering.

    Integer,    Intent(in) :: n1, n2    ! Subsets of n1:n2
    Integer,    Intent(in) :: p         ! max size of set extension
    Integer,    Intent(in) :: l         ! length of the prior set

    Integer, Intent(inout) :: k         ! subset index: initially for the prior set,
                                        !   finally = index of last set
    Integer, Intent(inout) :: ls(0:,:)  ! (0:p,k) ls(1:i,k)= kth non-empty subset
                                        !         of size 'i' = ls(0,k), with ls(1:l,k)=
                                        !         the prior set.
! Local:
    Integer :: i, j, k0, l1, ln, p1

    If (p <= 0) Return
    l1= l + 1;  k0= k;  p1= p - 1;  ln= Size(ls,2)

    Do i= n1,n2
      If (k >= ln) then
        Call Out ("Ran out of space in List_subsets at k",k)
        Exit
      End if
      
      k= k + 1;  ls(0,k)= l1
      If (l > 0) ls(1:l,k)= ls(1:l,k0)
      ls(l1,k)= i;  p1= Min(p1, n2-i)
      If (p1 <= 0) Cycle
      Call List_subsets (i+1,n2,p1,l1, k,ls)
    End do
  End Subroutine List_subsets

	
  Pure Integer Function ind_subF(sub, nsb)

!   Compute the index of a subset 'sub' of size <= 'm' of the set {1,...,n)
!   using the order specified by 'nsb', assuming that 'sub' is in 
!   increasing order.

  	Implicit None
    Integer, Intent(in) :: sub(:)    ! (<=m) The ordered subset of {1,...,n}
    Integer, Intent(in) :: nsb(:,:)  ! (m,n) nsb(i,j)= # subsets 
                                     !       of size <= i that contain final element j
! Local:
    Integer :: i, m, n, p, mp

    m= Size(nsb,1);  n= Size(nsb,2)
    p= Size(sub);  mp= m - p; ind_subF= 0

    Do i= p,1,-1
      ind_subF= ind_subF + 1 + Sum(nsb(mp+i,:sub(i)-1))
    End do
  End Function ind_subF


  Pure Subroutine sub_subF (ind, m,n, nsb, p,sub)

!   Compute the ordered subset 'sub' of {1,...,n} and its size 'p' 
!   that correpsonds to the index 'ind', assuming the subset listing
!   specified by nsb(m,n).

  	Implicit None
	Integer, Intent(in) :: ind       ! The specified subset index
	Integer, Intent(in) :: m, n      ! Dimensions of the nsb listing
    Integer, Intent(in) :: nsb(:,:)  ! (m,n) nsb(i,j)= # subsets 
                                     !       of size <= i that contain final element j
    Integer, Intent(out) :: p         ! Computed size of 'sub'
    Integer, Intent(out) :: sub(:)    ! (<=m) The subset corresponding to 'ind'
! Local:
    Integer :: j, k, ir, m1, sm, sm1

    m1= m + 1;  ir= ind;  sm= 0;  k= n

    Do p= 1,m
      Do j= 1,k
        sm1= sm + nsb(m1-p,j)
        If (ir <= sm1) Exit;  sm= sm1
      End do

      ir= ir - sm;  If (ir <= 0) Exit

      sub(m1-p)= j;  sm= 1;  k= j - 1
    End do

    p= p - 1;  sub(:p)= sub(m1-p:m);  sub(p+1:)= 0
  End Subroutine sub_subF


  Pure Subroutine n_subF (nsb)

!   # subsets of size <= 'm' of a set of size 'n' 
!   containing the nth, or final, element

    Integer, Intent(out) :: nsb(:,:)  ! (m,n)
! Local:
    Integer :: i, j, m, n
    
    m= Size(nsb,1);  n= Size(nsb,2)

    nsb(1,:)= 1;  nsb(2:,1)= 1
    Forall(j=2:n) nsb(2,j)= j;  nsb(3:,2)= 2

    Do i= 3,m
      Do j= i,n
        nsb(i,j)= 1 + Sum(nsb(i-1,:j-1))
      End do
      nsb(i+1:,i)= nsb(i,i)
    End do
  End Subroutine n_subF

  Recursive Subroutine List_combo (m,p,i, set, sub)
    
!   List all combinations of 'n' integers taken 'm' at at time =
!   All subsets of size 'm' of a set of increasing integers of size 'n'

    Integer,     Intent(in) :: m         ! Size of the combinations = subsets
    Integer,     Intent(in) :: p         ! Next subset position in 'sub(:,i)'
    Integer,     Intent(in) :: i         ! Index of next subset
    Integer,     Intent(in) :: set(:)    ! (n)  Set of 'n' integers in increasing order
      
    Integer,  Intent(inout) :: sub(:,:)  ! (m,l) = (subset,index) Subsets of 'set', 
                                         !    extending a prior assignment on (1:p-1)
                                         !    so that sub(p:m,:) is new.
! Local:
    Integer :: j, k, l, n, nm
          
    n= Size(set)
    If (n == 1) then
      sub(p,i)= set(1);  Return
    End if
    nm= n - m + 1;  k= i
      
    Do j= 1,nm
      l= k + N_subsets(n-j,m-1) - 1
      sub(p,k:l)= set(j)
      
      If (m > 1 .and. j < n) then
        Call List_combo (m-1,p+1,k, set(j+1:), sub)
        k= l + 1
      Else if (j < nm) then
        k= k + 1
      End if
    End do
      
  End Subroutine List_combo

! Compute the index of a permutation
  
  Recursive Integer Function Perm_indx (ls)
    Integer, Intent(in) :: ls(:) ! (n) List of remaining 
                                 !     integers in the permutation of 1...m
! Local:
    Integer, Allocatable :: lt(:)
    Integer :: n

    n= Size(ls)
    
    If (n <= 1) then
      Perm_indx= 1
    Else if (n == 2) then
      If (ls(1) < ls(2)) then
        Perm_indx= 1
      Else
        Perm_indx= 2
      End if
    Else
      Allocate(lt(n-1))
      Where (ls(2:) > ls(1))
        lt= ls(2:) - 1
      Else where
        lt= ls(2:)
      End where
      
      Perm_indx= (ls(1)-1)*Factorial(n-1) + Perm_indx(lt)
    End if
  End Function Perm_indx

! Compute the permutation corresponding to an index
  
  Recursive Subroutine Perm_set (p,perm)
    Integer,  Intent(in) :: p       ! Index of the permutation
    Integer, Intent(out) :: perm(:) ! (n) The permutation of 1...n indexed by 'p'
! Local:
    Integer :: i, n, m1, mx, n1, p1
    
    n= Size(perm);  mx= Factorial(n)
    
    If (n <= 0 .or. p <= 0 .or. p > mx) then  ! Error
      Call Out ("Parameter Error in 'Perm_set'"); Stop
    Else if (n == 1) then
      perm= 1
    Else if (n == 2) then
      If (p == 1) then
        perm= (/1,2/)
      Else
        perm= (/2,1/)
      End if
    Else
      m1= mx / n;  n1= (p-1) / m1
      p1= p - n1*m1;  i= n1 + 1;  perm(1)= i

      Call Perm_set (p1, perm(2:))
      Where (perm(2:) >= i) perm(2:)= perm(2:) + 1
    End if
  End Subroutine Perm_set

! If Add > 0, compute the index 1...2^n of a sequence of integers 1...base 
! in base 'base', else the index 0...2^n - 1 for integers 0...base-1
  
  Pure Integer Function Base_num(Add,base,ls)
  
    Integer,  Intent(in) :: Add   ! Assume integers 1...base if > 0, else 0...base-1
    Integer,  Intent(in) :: base  ! The base integer
    Integer,  Intent(in) :: ls(:) ! (n) The corresponding sequence of integers
! Local:
    Integer :: i, n, p

    n= Size(ls);  p= base
    
    If (Add > 0) then
      Base_num= ls(n)
      Do i= n-1,1,-1
        Base_num= Base_num + (ls(i)-1)*p;  p= p * base
      End do
    Else
      Base_num= ls(n)
      Do i= n-1,1,-1
        Base_num= Base_num + ls(i)*p;  p= p * base
      End do
    End if
  End Function Base_num

! For any index 'ind', computer its coefficient sequence in base 'base',
! optionally adding 1 to each coefficient.
  
  Subroutine Base_seq (Add,base,n,ind, seq)
  
    Integer,  Intent(in) :: Add     ! Add 1 if > 0
    Integer,  Intent(in) :: base    ! The base integer
    Integer,  Intent(in) :: n       ! # integers in the sequence
    Integer,  Intent(in) :: ind     ! Index of the sequence
    Integer, Intent(out) :: seq(:)  ! (n) Sequence of 'n' integers selected from 0...base-1,
                                    !  or 1...base if Add > 0
! Local:
    Integer :: i, j, k, pw(0:n-1)
    
    pw(0)= 1
    Do i= 1,n-1
      pw(i)= base * pw(i-1) 
    End do
    
    If (Add > 0) then
      k= ind - 1
      Do i= 1,n-1
        j= k / pw(n-i);  k= k - j*pw(n-i)
        seq(i)= j+1
      End do
      seq(n)= k+1
    Else
      k= ind
      Do i= 1,n-1
        j= k / pw(n-i);  k= k - j*pw(n-i)
        seq(i)= j
      End do
      seq(n)= k
    End if
  End Subroutine Base_seq

End Module Factorials