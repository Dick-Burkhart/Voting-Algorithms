
Module Newton_operators

   Use Clust_gen
   Use Cholesk
   Use Util
   Use Output
   Use Types
   Use Precisn
   Implicit None
   
!  Data and parameters for the line search for cluster convergence, either Gauss-Newton or nonlinear GMRES
   
   Real(Dblp),   Parameter :: LS_eps= 3.605E-7  ! Finite difference epsilon for use with double precision
   Logical,      Save      :: LS_full= .false.  ! Use full Hessians in the Newton solver, versus block diagonal
   Logical,      Parameter :: Analytic= .true.  ! Use Analytic_Jacobian if true, else Diff_Jacobian

   Integer                 :: LS_ncl            ! # regular clusters, possibly decreasing
   
   Type(Multi_listD), Allocatable :: LS_Memb(:) ! (0:ns) Slate cluster membership data
   Real(Dblp), Allocatable :: LS_rate0(:,:)     ! (nc,ncl) Initial cluster mean rating vectors
   Real(Dblp), Allocatable :: LS_rate1(:,:)     ! (nc,ncl) Updated cluster mean rating vectors
                                                !          (from LS_rate0 via line search)
   Real(Dblp), Allocatable :: LS_rate2(:,:)     ! (nc,ncl) Final cluster mean rating vectors (from LS_rate1)
   Integer,    Allocatable :: LS_lst(:)         ! (ncl) Listing of the original clusters when
                                                !       LS_ncl has been reduced, without reordering

!  Data for the Gauss Newton line search for cluster convergence
   
   Real(Dblp), Allocatable :: GN_Jacob(:,:)   ! (nu,nu)     Jacobian matrix of the update U(x)
   Real(Dblp), Allocatable :: GN_JacobI(:,:)  ! (nu,nu)     I - GN_Jacob = Jacobian matrix of the update increment x - U(x)
   Real(Dblp), Allocatable :: GN_HinvJ(:,:)   ! (nu,nu)     Matrix for the preconditioned gradient = GN_Hess^(-1).GN_JacobI^t
   Real(Dblp), Allocatable :: GN_Hess(:,:)    ! (nu,nu)     Damped Gauss-Newton Hessian = GN_JacobI^t.GN_JacobI + dg.I
   Real(Dblp), Allocatable :: GN_HesB(:,:,:)  ! (nc,nc,ncl) Block diagonal version of GN_Hess
   Real(Dblp), Allocatable :: GN_Grad(:)      ! (nu)        Gradient of the update increment = GN_JacobI^t (x - U(x))
   Real(Dblp), Allocatable :: GN_up0(:)       ! (nu)        d = 0 update increment

   Real(Dblp) :: Jnorm(4)  ! Norms of the convergence operator: 
                           !   L_inf, L_1, Frobenius, L_2= spectral radius
   Real(Dblp) :: Pnorm(4)  ! Corresponding norms of the preconditioned convergence operator

!  Data for the Levenburg-Marquart line search for the nonlinear GMRES cluster convergence
   
   Integer                 :: LM_k            ! Krylov basis dimension 'kv'
   Real(Dblp)              :: LM_fac          ! LM factor for computing lambda(t)
   Real(Dblp), Allocatable :: LM_trifac(:,:)  ! (kv,-1:kv) Factored normal equations matrix, plus auxillary vectors
   Real(Dblp), Allocatable :: LM_basis(:,:)   ! (nu,kv+5) Basis vectors for the Krylov subspace,
                                              !   plus storage of additional vectors for the line search
   Integer, Parameter :: LM_x= 1, LM_r= 2, LM_t= 3, LM_u= 4, LM_p= 5, LM_v= 6
   Integer            :: LM_e                 ! = LM_v + LM_k - 1
   
   Integer            :: nGen_clust           ! Counts # calls to 'Generate_clusters' = the basic 
                                              !   centroiding operation, determining cluster membership
   
   
   Contains  
    
     Subroutine GN_Func_0 (nc,n0, ft, ier)
  
!     To be used before the Guass-Newton update line search
!       One call to Generate_clusters.
  
      Integer,     Intent(in) :: nc      ! # candidates
      Integer,     Intent(in) :: n0      ! Initial # clusters
      Real(Dblp), Intent(out) :: ft(:)   ! (2) Values to minimize
                                         !   1 = The functional
                                         !   2 = Norm of its gradient
      Integer,    Intent(out) :: ier     ! Error code: 1 = normal, 0 = warning, -1 = error
! Local:
      Integer :: n1, nu
      
!     Initialize the Hessian and compute the corresponding gradient
      
      n1= n0
      Do
        nu= nc*n1;  LS_full= n1 <= 3
        Call Init_Hess (nu,n1, ft, ier)  ! Functional and gradient

        If (ier < 0) then ! Serious error
          Return
        Else if (ier > 0) then  ! Success
          Exit
        Else  ! Start over with fewer clusters
          n1= LS_ncl;  LS_rate0(:,:n1)= LS_rate1(:,:n1)  
        End if
      End do
      
      LS_rate2(:,:LS_ncl)= LS_rate1(:,:LS_ncl)  ! then to LS_rate2
      
    End Subroutine GN_Func_0
  
  
    Subroutine GN_Func_t (t, tf, ier)
  
!     To be used with the Guass-Newton update line search.
!     2 calls to Generate_clusters, normally, but 3 initially.
  
      Real(Dblp), Intent(inout) :: t          ! Damping parameter, t >= 0, to be added to the Hessian
      Real(Dblp),   Intent(out) :: tf(0:3,2)  ! Update data at 't'
      Integer,      Intent(out) :: ier        ! Error code
! Local:
      Real(Dblp) :: t_eps, ft
      Integer    :: nc, n0, nu, ncl
    
      tf= 0d0;  ier= 1;  nc= LS_Memb(0)%n
      n0= LS_ncl;  nu= nc*n0;  LS_full= n0 <= 3
    
!     Functional at 't'
    
      Call Func_t (nc,n0,nu, t, ft, ier)

      ncl= LS_ncl;  If (ncl < 1) ier= -1
      
      If (ier < 0) then
        Call Out ("Error in 'GN_Func_t' with error code",ier, ln=1)
        Call Out ("for the functional at damping value",t);   Return
      Else if (ncl < n0) then
        Call Out ("Warning in 'GN_Func_t': # clusters reduced from",n0, "to", ncl, ln=1)
        Call Out ("Retained clusters",LS_lst(:ncl))
        ier= 0;  Return
      End if    
    
!     Derivative of the functional at 't
    
      tf(0,1)= t;  tf(1,1)= ft;  t_eps= t + LS_eps

      Call Func_t (nc,n0,nu, t_eps, ft, ier)
      ncl= LS_ncl;  If (ncl < 1) ier= -1
      
      If (ier < 0) then
        Call Out ("Error in 'GN_Func_t' with error code",ier, ln=1)
        Call Out ("for the perturbed functional at damping value",t_eps);   Return
      Else if (ncl < n0) then
        Call Out ("Warning in 'GN_Func_t': # clusters reduced from",n0, "to", ncl, ln=1)
        Call Out ("Retained clusters",LS_lst(:ncl))
        ier= 0;  Return
      End if    
    
      tf(0,2)= t_eps;  tf(1,2)= (ft - tf(1,1)) / LS_eps
    End Subroutine GN_Func_t
    
    Subroutine Func_t (nc,n0,nu, t, ft, ier)
  
!     To be used with the Newton update routine. 
!     One call to Generate_clusters.
    
      Integer,        Intent(in) :: nc      ! First coordinate size
      Integer,        Intent(in) :: n0      ! # clusters LS_ncl
      Integer,        Intent(in) :: nu      ! # unknowns = nc*n0
      
      Real(Dblp),  Intent(inout) :: t       ! Damping parameter, t >= 0
      Real(Dblp),    Intent(out) :: ft      ! Functional = update change
      Integer,       Intent(out) :: ier     ! Error code: 1 = normal, 0 = reduced # clusters (warning), -1 = error
!   Local:
      Type(Multi_listD) :: PD(0:n0)   !  Cluster membership and derivative data
      Real(Dblp)        :: up(nc,n0)  !  Mean ratings matrix - input
      Integer           :: cl, ns, ncl
      
!     Function evaluation with the damping parameter 't' added to the diagonal of the Hessian
      
      ns= Ubound(LS_Memb,1);  ft= 0d0;  LS_full= n0 <= 3
      
      Call Solve_Hess (t, up)

      Call Test_update (up)  !  LS_rate1= LS_rate0 + up
      
!     Evaluate the updated cluster mean ratings (at LS_rate1)
      
      Call Allocate_PD (nc,n0,ns, PD)
      Call Generate_clusters (nc,n0,ns, LS_Memb, LS_rate1(:,:n0), PD, &
                              LS_rate2(:,:n0),ft)
        nGen_clust= nGen_clust + 1

      ncl= PD(0)%n;  LS_ncl= ncl
      ier= PD(0)%k;  If (ncl < 1) ier= -2

      If (ier < 0) then
        Call Out ("Error in 'Func_t' error code",ier, ln=1)
      Else if (ncl < n0) then
        Call Out ("Warning in 'Func_t' # clusters reduced from",n0, "to",ncl, ln=1)
        ier= 0;  LS_lst(:ncl)= PD(0)%lt(:ncl);  LS_lst(ncl+1:)= -1
      End if
      Call DeAlloc_Multi_list_ar (PD)
      
    End Subroutine Func_t
    
    
    Subroutine Init_Hess (nu,n0, ft, ier)
    
!     Compute the Jacobian of U(x), then of F(x) = x - U(x) for computing 
!     the full damped Hessian 'GN_Hess' or the block diagonal 
!     damped Hessian 'GN_HesB', for x = LS_rate0 and U(x) = LS_rate1

!     GN_Hess = Grad(F)^t . Grad(F), where Grad(F)= I - J
!     and J = Grad(U) = Jacobian = GN_Jacob, with Grad(F) = GN_JacobI
!     and GN_grad= Grad(F)^t . F
    
      Integer,     Intent(in) :: nu     ! nc*n0 = # unknowns
      Integer,     Intent(in) :: n0     ! Initial # clusters = LS_ncl
      
      Real(Dblp), Intent(out) :: ft(:)  ! (2) The functional |F(x)|/nu and norm of its gradient

      Integer,    Intent(out) :: ier    ! Error code: 1 = normal, 0 = warning, start over with 
                                        !  fewer clusters,  < 0 = unknown error
!   Local:
      Type(Multi_listD) :: PD(0:n0) 
      Real(Dblp) :: up(nu)  ! (nu) for GN_Hess or (n0) for GN_HesB
      Real(Dblp) :: dg, change
      Integer    :: i, j, k, i0, nc, nr, ns, ncl
      
      If (Allocated(GN_Jacob)) DeAllocate(GN_Jacob);  Allocate(GN_Jacob(nu,nu))
      nc= LS_Memb(0)%n;  ns= Ubound(LS_Memb,1)
      
      If (Analytic) then
        Call Allocate_PD (nc,n0,ns, PD)
        Forall (k=1:n0) PD(k)%sx= LS_rate0(:,k)

        Call Analytic_Jacobian (nc,n0,ns, LS_Memb, PD, LS_ncl,LS_lst, &
                                GN_Jacob, LS_rate1(:,:n0),change)

        ier= PD(0)%k;  Call DeAlloc_Multi_list_ar (PD)
        If (ier <= 0) Return
      Else
        Call Diff_Jacobian (nc,n0,nu, LS_Memb, LS_rate0(:,:n0), &
                            GN_Jacob, LS_rate1(:,:n0),change)
      End if

      ft(1)= change;  ft(2)= -1

      If (Allocated(GN_JacobI)) DeAllocate(GN_JacobI, GN_Grad) 
      Allocate(GN_JacobI(nu,nu), GN_Grad(nu))

      GN_JacobI= -GN_Jacob;  Call Matrix_diagonal ("Add", 1.0d0, GN_JacobI)

      up= Reshape( LS_rate0(:,:n0) - LS_rate1(:,:n0), (/nu/))

      GN_Grad= Matmul(Transpose(GN_JacobI),up)

      ft(2)= Sqrt(Sum(GN_Grad**2)/nu)
      
      If (LS_full) then
        If (Allocated(GN_Hess)) DeAllocate(GN_Hess);  Allocate(GN_Hess(nu,nu))
          
        Do i= 1,nu
          GN_Hess(i,i)= Sum(GN_JacobI(:,i)**2)
          Do j= i+1,nu
            GN_Hess(i,j)= Dot_product(GN_JacobI(:,i),GN_JacobI(:,j))
            GN_Hess(j,i)= GN_Hess(i,j)
          End do
        End do
      Else
        If (Allocated(GN_HesB)) DeAllocate(GN_HesB);  Allocate(GN_HesB(nc,nc,n0))
          
        Do k= 1,n0
          i0= (k-1)*nc
          Do i= 1,nc
            GN_HesB(i,i,k)= Sum(GN_JacobI(:,i0+i)**2)
            Do j= i+1,nc
              GN_HesB(i,j,k)= Dot_product(GN_JacobI(:,i0+i),GN_JacobI(:,i0+j))
              GN_HesB(j,i,k)= GN_HesB(i,j,k)
            End do
          End do
        End do
      End if
    End subroutine Init_Hess


    Subroutine Solve_Hess (d, up)

!     Solve the Newon-Gauss update equations: up= -(GN_Hess + d.I)^(-1) . GN_Grad
    
      Real(Dblp),  Intent(in) :: d        ! Damping value, d >= 0
      
      Real(Dblp), Intent(out) :: up(:,:)  ! (nc,n0) Cluster update vector, from the damped gradient
                                          !   up= -(GN_Hess + d.I)^(-1) . GN_Grad
!   Local:      
      Real(Dblp), Allocatable :: Mat(:,:) ! Full normal equations matrix, from the Hessian 'GN_Hess'
      Real(Dblp), Allocatable :: up1(:)   ! One dimensional gradient
      Integer :: k, n0, nc, nu
      
      nc= Size(up,1);  n0= Size(up,2);  nu= nc*n0
      
      If (LS_full) then
        Allocate (Mat(nu,nu+1), up1(nu));  up1= -GN_Grad(:nu)
        Mat(:,:nu)= GN_Hess(:nu,:nu);  Mat(:,nu+1)= 0.0d0
        If (d > 0) Call Matrix_diagonal ("Add", d, Mat)

        Call Sym_solve (.true., Mat, up1)
        up= Reshape(up1, (/nc,n0/))
      Else
        Allocate (Mat(nc,nc+1));  up= -Reshape(GN_Grad(:nu), (/nc,n0/))

        Do k= 1,n0
          Mat(:,:nc)= GN_HesB(:,:,k);  Mat(:,nc+1)= 0.0d0
          If (d > 0) Call Matrix_diagonal ("Add", d, Mat)

          Call Sym_solve (.true., Mat, up(:,k))
        End do
      End if
    
    End subroutine Solve_Hess

    
    Subroutine Test_update (up)

!   Update the cluster mean rating vectors: 
!   LS_rate1= LS_rate0 + up, modified as necessary 
!   to impose bounds on LS_rate1

      Real(Dblp), Intent(inout) :: up(:,:)  ! (nc,n0)
!   Local:
      Integer :: n0, p

      n0= Size(up,2);  p= Size(LS_rate1,2)

      If (n0 < p) then
        LS_rate1(:,:n0)= LS_rate0(:,:n0) + up
        Call Limit_ratings (LS_rate1(:,:n0), Exceed)
        LS_rate1(:,p:)= 0;  LS_rate0(:,p:)= 0
        If (Any(Exceed)) up= LS_rate1(:,:n0) - LS_rate0(:,:n0)
        Call Out ("Warning in 'Test_update': Decrease in # clusters")
        LS_ncl= n0
      Else
        LS_rate1= LS_rate0 + up
        Call Limit_ratings (LS_rate1, Exceed)
        If (Any(Exceed)) up= LS_rate1 - LS_rate0
      End if
    End Subroutine Test_update
  

    Subroutine Analytic_Jacobian (nc,n0,ns, Memb,PD, ncl,lst, &
                                  Jacob, rate_out,change)

!     Compute the Jacobian of the convergence operator analytically.
!     One call to Generate_clusters.
  
!     Each row vector of 'Jacob' = gradient of that coordinate
!     Each column vector of 'Jacob' = all the partials wrt that coordinate
  
!     Assumes prior call to "Allocate_PD" or "Allocate_clust"

      Integer,              Intent(in) :: nc        ! # candidates
      Integer,              Intent(in) :: n0        ! # clusters
      Integer,              Intent(in) :: ns        ! # slate ballots
      Type(Multi_listD),    Intent(in) :: Memb(0:)  ! (0:ns)  Slate ballot data
    
      Type(Multi_listD), Intent(inout) :: PD(0:)  ! (0:n0) Cluster membership data
                                                  ! 0   (Cluster set data):
                                                  !   %k       : success / error code
                                                  !   %n       : # clusters (output)
                                                  !   %lt(ncl) : listing of original clusters if # clusters is reduced (output)
                                                  !   %wt(ns)  : Sum of initial memberships over all clusters, (output)
                                                  !              to determine the fuzzy factor for each slate ballot
                                                  !   %M0(ns,0:1): Fuzzy membership scale factors (output)
                                                  !               0: fuzzy scale factor
                                                  !               1: fuzzy derivative factor
                                                  ! 1:ncl (Cluster data):
                                                  !   %n       : # slate ballots members
                                                  !   %fux     : cluster width
                                                  !   %sum_wt  : total weight of the cluster = Sum(%wt)

                                                  !   %sx(nc)  : cluster mean rating vector (over slate ballot 'sl' rating vectors 
                                                  !              Memb(sl)%sx weighted by Memb(sl)%ux * %rx(sl) and normalized by %ux
                                                  !   %px(nc)  : variance of %sx
                                                  !   %tx(0:nc): 1 = centered and normalized %sx.  0 = |%sx|
                                                  !   %ux(nc)  : sum of fuzzy memberhsips of the slate ballots 'sl' using 
                                                  !              normalization slate ballot weights = 
                                                  !              Sum(Memb(sl)%ux * %rx(sl))

                                                  !   %ls(n)   : list of slate ballot members
                                                  !   %rx(n)   : unweighted memberhsips of the slate ballots
                                                  !   %wt(n)   : weighted memberhsips of the slate ballots

                                                  !   %M0(n,0:2): Slate ballot membership data
                                                  !               0: initial membership after the cosine cutoff
                                                  !               1: derivative of this membership wrt the modified correlation
                                                  !               2: fuzzy derivative factor * initial membership =
                                                  !                  PD(0)%M0(:,1) * PD(k)%M0(:,0)
                                                  !   %M2(nc,n) : Partial derivatives of the initial
                                                  !               memberships of the slate ballots in this cluster 
                                                  !               with respect to the mean vector of the cluster,
                                                  !               ordered by decreasing correlation with the slate ballots
                                                  !   %M3(nc,n): %ux weighted difference between the slate ballot mean rating 
                                                  !              vectors of the cluster and the cluster mean rating vector,
                                                  !              (Memb(sl)%sx - PD(k)%sx) * (Memb(sl)%ux / PD(k)%ux)
                                                      
      Integer,    Intent(out) :: ncl           ! Possibly reduced # clusters
      Integer,    Intent(out) :: lst(:)        ! (n0=>ncl) List retained clusters in case of a reduction
    
      Real(Dblp), Intent(out) :: Jacob(:,:)    ! (nc*n0,nc*n0) Jacobian, first nu= nc*ncl valid on output
      Real(Dblp), Intent(out) :: rate_out(:,:) ! (nc,n0) Computed cluster mean vectors
      Real(Dblp), Intent(out) :: change        ! Update norm = Sqrt(Sum(up**2) / nu) for up= rate_out - rate_in norm
! Local:
      Type(Adjacency) :: Lsb(ns)  ! Lists the clusters containing each slate ballot
      Real(Dblp) :: rate_in(nc,n0)
      Integer    :: d, i, k, l, n, ik, il, k1, k2, l0, sl, ier

      ncl= n0;  lst= "ID";  Jacob= 0;  rate_out= 0;  change= -1
    
      Forall(k=1:n0) rate_in(:,k)= PD(k)%sx
      Call Limit_ratings (rate_in, Exceed)

      Call Generate_clusters (nc,n0,ns, Memb, rate_in,PD, rate_out,change)
        nGen_clust= nGen_clust + 1
    
      ier= PD(0)%k;  ncl= PD(0)%n;  If (ncl < 1) ier= -1

      If (ier < 0) then
        Call Out ("Error in 'Analytic_Jacobian' with error code",ier, ln=1);  Return

      Else if (ncl < n0) then
        Call Out ("Warning in 'Analytic_Jacobian': # clusters reduced from",n0, "to",ncl, ln=1)
        lst(:ncl)= PD(0)%lt(:ncl);  lst(ncl+1:)= -1
        Call Out ("List of those retained",lst(:ncl))
        ier= 0;  Return
      End if
    
      Do k= 1,ncl
        If (Associated(PD(k)%M3)) DeAllocate(PD(k)%M3)
        n= PD(k)%n;  Allocate(PD(k)%M3(nc,n))
        Forall(i=1:n) PD(k)%M3(:,i)= (Memb(PD(k)%ls(i))%sx - PD(k)%sx) * (Memb(PD(k)%ls(i))%ux / PD(k)%ux)
      End do

      Lsb%nl= 0
      Do k= 1,ncl
        Lsb(PD(k)%ls)%nl= Lsb(PD(k)%ls)%nl + 1
      End do

      Do sl= 1,ns
        n= Lsb(sl)%nl;  If (n > 0) Allocate(Lsb(sl)%ls(n), Lsb(sl)%vl(n))
      End do

      Lsb%nl= 0
      Do k= 1,ncl
        Lsb(PD(k)%ls)%nl= Lsb(PD(k)%ls)%nl + 1

        Do i= 1,PD(k)%n
          sl= PD(k)%ls(i);   n= Lsb(sl)%nl 
          Lsb(sl)%ls(n)= k;  Lsb(sl)%vl(n)= i
        End do
      End do
    
!     Compute the Jacobian:  the partials of the updated cluster rating vectors with respect 
!     to the prior cluster rating vectors.
    
      Do k= 1,ncl     ! First coord cluster, function
        k1= nc*(k-1)+1;  k2= k*nc  ! k1:k2 = first coord cand limits
      
        slate_loop : Do ik= 1,PD(k)%n
          sl= PD(k)%ls(ik)

          Do n= 1,Lsb(sl)%nl  ! Clusters containing this slate. Partials wrt to these
            l= Lsb(sl)%ls(n);  l0= nc*(l-1)

            If (l == k) then
              Forall(d=1:nc) Jacob(k1:k2,l0+d)= Jacob(k1:k2,l0+d) + PD(k)%M3(:,ik) * &
                             (PD(k)%M0(ik,2) + PD(0)%M0(sl,0)) * PD(k)%M2(d,ik)
            Else
              il= Lsb(sl)%vl(n)
              Forall(d=1:nc) Jacob(k1:k2,l0+d)= Jacob(k1:k2,l0+d) + PD(k)%M3(:,ik) * &
                             PD(k)%M0(ik,2) * PD(l)%M2(d,il)
            End if
          End do
        End do slate_loop
        
      End do
    End Subroutine Analytic_Jacobian


    Subroutine Diff_Jacobian (nc,ncl,nu, Memb, rate_in, Jacob, rate_out, change)
      
!     Compute the Jacobian by finite differences of by nu+1 calls to 'Generate_clusters'.

!     Each row vector of 'Jacob' = gradient of that coordinate
!     Each column vector of 'Jacob' = all the partials wrt that coordinate

!     It is assumed that there reduction in the number of clusters.
  
      Integer,           Intent(in) :: nc        ! First coordinate size
      Integer,           Intent(in) :: ncl       ! # clusters
      Integer,           Intent(in) :: nu        ! # unknowns = nc*ncl
      
      Type(Multi_listD), Intent(in) :: Memb(0:)  ! (0:ns) Slate ballot cluster data
      
      Real(Dblp), Intent(inout) :: rate_in(:,:)  ! (nc,ncl) Input cluster mean vectors
      
      Real(Dblp),   Intent(out) :: Jacob(:,:)    ! (nu,nu) Euclidean finite difference Jacobian
      Real(Dblp),   Intent(out) :: rate_out(:,:) ! (nc,ncl) Computed cluster mean vectors
      Real(Dblp),   Intent(out) :: change        ! Update norm = Sqrt(Sum(up**2) / nu) for up= rate_out - rate_in norm
!   Local:
      Type(Multi_listD) :: PD(0:ncl)  ! Cluster set data
      Real(Dblp) :: up0(nu)  ! Initial update of the mean vectors %sx
      Real(Dblp) :: up1(nu)  ! Subsequent updates of the mean vectors %sx
      Real(Dblp) :: tmp, res
      Integer    :: i, j, k, ns

      ns= Ubound(Memb,1);  Jacob= -1;  rate_out= -1
      Call Allocate_PD (nc,ncl,ns, PD)

      Call Generate_clusters1 (nc,ncl,ns, Memb, rate_in, PD, rate_out)
        nGen_clust= nGen_clust + 1

      change= Sqrt(Sum((rate_out - rate_in)**2) / nu)
    
      up0= Reshape(rate_out, (/nu/))
        
      j= 1;  i= 0

      Do k= 1,nu
        i= i + 1  
        If (i > nc) then ! Next cluster
          j= j + 1;  i= 1
        End if
          
        tmp= rate_in(i,j);  rate_in(i,j)= tmp + LS_eps
          
        Call Generate_clusters1 (nc,ncl,ns, Memb, rate_in, PD, rate_out)

        rate_in(i,j)= tmp

        up1= Reshape(rate_out, (/nu/))
        Jacob(:,k)= (up1 - up0) / LS_eps
      End do

      nGen_clust= nGen_clust + nu
      
      Call DeAlloc_Multi_list_ar (PD)

    End Subroutine Diff_Jacobian


    Subroutine Allocate_PD (nc,ncl,ns, PD)
    
!     Allocate and initialize cluster arrays for use by 'Generate_clusters'
    
      Integer,              Intent(in) :: nc      ! First coordinate size
      Integer,              Intent(in) :: ncl     ! # clusters
      Integer,              Intent(in) :: ns      ! # slate ballots
      Type(Multi_listD), Intent(inout) :: PD(0:)  ! (0:m) Cluster data
                                                  ! For 0
                                                  !   %k       : success / error code
                                                  !   %n       : final # regular clusters = 'ncl'
                                                  !   %lt(ncl) : original numbering of clusters if # clusters is reduced 
                                                  !              (no change in order)
                                                  !   %wt(ns)  : For each slate: sum of initial memberships over all clusters
                                                  !              to determine the fuzzy reduction factor
                                                  !   %M0(ns,0:1): Fuzzy membership reduction factors
                                                  !               0: scale factor
                                                  !               1: derivative factor

                                                  ! For 1:m
                                                  !   %k  = # fuzzy reduced slate memberships
                                                  !   %l  = # slate ballots with membership > 1/3
                                                  !   %m  = # slate ballots which are full members
                                                  !   %n  = # slate ballots which are partial or full members
                                                  !   %o       : 1 for regular clusters, 2 for the cluster of independents
                                                  !   %fux     : cluster width
                                                  !   %sum_wt  : total weight of the cluster

                                                  !   %px(nc) = variance of %sx
                                                  !   %sx(nc) = mean rating vector
                                                  !   %tx(0:nc)= centered and normalized rating vector
                                                  !   %ux(nc)  : total cluster weight, incorporating slate ballot variance 
                                                  !              normalization (Memb%ux * %rx)
!   Local
      Integer :: i, k, m
      
      PD(0)%k= 1;  PD(0)%n= ncl
      Allocate (PD(0)%lt(ncl), PD(0)%wt(ns), PD(0)%M0(ns,0:1))
      PD(0)%lt= "ID";  PD(0)%wt= 0;  PD(0)%M0(:,0)= 1;  PD(0)%M0(:,1)= 0  
      
      m= Min(Ubound(PD,1),ncl+1)
      
      Do k= 1,m
        PD(k)%l= 0;  PD(k)%l= 0;  PD(k)%m= 0;  PD(k)%n= 0
        PD(k)%o= 1;   If (k > ncl) PD(k)%o= 2
        PD(k)%fux= 0;  PD(k)%sum_wt= 0

        Allocate (PD(k)%px(nc), PD(k)%sx(nc), PD(k)%tx(0:nc), PD(k)%ux(nc))
        PD(k)%px= 0;  PD(k)%sx= 0;  PD(k)%tx= 0;  PD(k)%ux= 0
      End do
    End Subroutine Allocate_PD


    Subroutine Damp_Hess (nc,ncl,dmp, Pmat)

!     'Pmat' = Inverse of the damped Hessian times the Hessian,
!     which is an approximation of the Jacobian
!     of the actual preconditioned operator used in the
!     Gauss-Newton method as computed in 'Init_Hess'.

!     This preconditioning operator is 
!     (GN_Hess + d.I)^(-1) . GN_grad, where GN_grad= Grad(F)^t . F
!     Thus the approximation is to take the gradient
!     only of the final 'F'. That is, the gradient

!     of (GN_Hess + d.I)^(-1) . Grad(F)^t . F is approx
!     (GN_Hess + d.I)^(-1) . (Grad(F)^t . Grad(F))
!     (GN_Hess + d.I)^(-1) . GN_Hess = Pmat

      Integer,     Intent(in) :: nc         ! First coordinate size
      Integer,     Intent(in) :: ncl        ! # clusters
      Real(Dblp),  Intent(in) :: dmp        ! damping parameter
      Real(Dblp), Intent(out) :: Pmat(:,:)  ! (nu,nu) Preconditioned Jacobian
    
!   Local:      
      Real(Dblp), Allocatable ::  Hess(:,:), Inv_Hess(:,:)
      Integer :: cl, i1, il, nu, nc0, nu0, err_col
      
      LS_full= ncl <= 3;  Pmat= 1001
      
      If (LS_full) then
        nu= nc*ncl;  If (.not.Allocated(GN_Hess)) Return
        nu0= Size(GN_Hess,1); If (nu0 /= nu) Return

        Allocate(Hess(nu,nu+1), Inv_Hess(nu,nu))

        Hess(:,:nu)= GN_Hess;  Hess(:,nu+1)= 0
        Call Matrix_diagonal ("Add", dmp, Hess(:,:nu))

        Call Invert_sym (.false.,.false., Hess, Inv_Hess, err_col)
        Pmat= Matmul(Inv_Hess,GN_Hess)
      Else
        If (.not.Allocated(GN_HesB)) Return
        nc0= Size(GN_HesB,1); If (nc0 /= nc) Return

        Allocate(Hess(nc,nc+1), Inv_Hess(nc,nc))
        Pmat= 0;  Hess(:,nc+1)= 0;  il= 0

        Do cl= 1,ncl
          i1= il + 1;  il= il + nc;  Hess(:,:nc)= GN_HesB(:,:,cl)
          Call Matrix_diagonal ("Add", dmp, Hess(:,:nc))

          Call Invert_sym (.false.,.false., Hess, Inv_Hess, err_col)
          Pmat(i1:il,I1:Il)= Matmul(Inv_Hess,GN_HesB(:,:,cl))
        End do
      End if
      
    End subroutine Damp_Hess

    
  Subroutine Jacobian_norms (nc,ncl,nu,dmp, Memb,rt_in, Jnorm,Pnorm)

!   Compute the analytical Jacobian of the convergence operator
!   and its norms 'Jnorm' and the corresponding 'Pnorm'
!   of the preconditioned operator. However, the matrix used
!   for the 'Pnorm' is only an approximation of the Jacobian
!   of the actual preconditioned operator.

!   Optionally, verify the Jacobian by a finite difference operator

    Integer,           Intent(in) :: nc, ncl, nu         ! nu= nc * ncl
    Real(Dblp),        Intent(in) :: dmp                 ! Damping parameter for he Hessian for 'Pnorm'
    Type(Multi_listD), Intent(in) :: Memb(0:)            ! (0:ns)  Slate cluster data
    Real(Dblp),        Intent(in) :: rt_in(:,:)          ! (nc,ncl)  Cluster mean vectors to be tested

    Real(Dblp),       Intent(out) :: Jnorm(:), Pnorm(:)  ! (4) L_inf, L_1, Frobenius, L_2= spectral radius
!  Local:
    Logical, Parameter :: Test_Jacob= .false.

    Type(Multi_listD) :: PD(0:ncl)
    Real(Dblp) :: rtD(nc,ncl), rt_outD(nc,ncl), JacobD(nu,nu)
    Real(Dblp) :: rt_out(nc,ncl), Jacob(nu,nu), Pmat(nu,nu)
    Integer    :: lst(ncl+1)
    
    Real(Dblp) :: dif, dif_GN, up, dif_up, penalty

    Logical :: Jacob_err
    Real    :: x, y
    Integer :: k, n, nl, ns, ier, iter, ls(ncl)

    If (pr_out >= 1) Call Out ("Enter 'Jacobian_norms' for # regular clusters",ncl, ln=1)

    If (Test_Jacob) then
      ns= Ubound(Memb,1);  Jacob_err= .false.
      Call Allocate_PD (nc,ncl,ns, PD)
      Forall (k=1:ncl) PD(k)%sx= rt_in(:,k) 

      Call Analytic_Jacobian (nc,ncl,ns, Memb, PD, nl,lst, Jacob, rt_out,up)
      
      Jacob_err= PD(0)%k < 0 .or. nl < ncl

      If (Jacob_err) then
        Call Out ("Error in 'Record_cluster_set': Analytic Jacobian not correct")
        Stop
      End if

      rtD= rt_in
      Call Diff_Jacobian (nc,ncl,nu, Memb, rtD, JacobD, rt_outD, dif_up)
       
      dif= Maxval(Abs(Jacob - JacobD))
      If (dif > 1.0E-5) then
        Call Out ("Finite diff Jacobian mismatch with analytic",dif)
      End if
      Jacob= JacobD
    Else
      Jacob= GN_Jacob
    End if

    Call Matrix_norms (Jacob, Jnorm, iter)

    Call Damp_Hess (nc,ncl, dmp, Pmat)

    If (Pmat(1,1) < 1000) then
      Call Matrix_norms (Pmat, Pnorm, iter)
    Else
      Pnorm= 1000
    End if

    x= Jnorm(4);  y= Pnorm(4)
    If (x > 1 .or. y > 1) then
      If (pr_out >= 1) then
        Call Out ("Jacobian spectral norm",x, "preconditioned spectral norm",y)
      End if
    End if
  End Subroutine Jacobian_norms 

End Module Newton_operators