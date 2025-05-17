Module Clust_gen
    
   Use Util
   Use Output
   Use Types
   Use Precisn
   Implicit None

   Logical    :: Exceed(2)       ! Used with Limit_ratings
   Integer    :: GN_rating       ! Rating alternative: 1 : rating data, 0 : weak rankings, -1 : strong ranking
   Real(Dblp) :: GN_dot_fac      ! Scale factor for the modified dot product: Dot_product_M1
   Real(Dblp) :: GN_rate_adj     ! Rating adjustment level. Used only if GN_rating < 1. = 0 if Gn_rating = 1
   Real(Dblp) :: GN_max_pt       ! Maximum absolute point value for the solution vector
   Real(Dblp) :: GN_min_wt       ! Minimum weight for regular clusters
   Real(Dblp) :: GN_min_mx       ! Minimum of the maximum rating, or adjusted ranking, for a regular cluster.
   
   Logical, Parameter :: GN_tri_motif= .false. ! Do triangular motif smoothing in cluster set formation.

 Contains
 
  Subroutine Generate_clusters (nc,n0,ns, Memb, rate_in,PD, rate_out, change)

!   Generate a new set of clusters from a prior set of mean rating vectors
!   or their centered and normalized versions.
  
!   Assumes prior call to Allocate_clust or Allocate_PD to set standard
!   array dimensions and initialize them to 0, etc.
  
    Integer,           Intent(in) :: nc           ! # candidates
    Integer,           Intent(in) :: n0           ! # clusters
    Integer,           Intent(in) :: ns           ! # slate ballots
    
    Type(Multi_listD), Intent(in) :: Memb(0:)     ! (0:ns)  Slate ballot data
    Real(Dblp),        Intent(in) :: rate_in(:,:) ! (nc,n0) Initial cluster mean vectors %sx,
                                                  !         assumed to be rate limited, according to 'GN_rating'
                                                  !         also non-centered for ranking (GN_rating < 1)
                                                  !         and not normalized 
    
    Type(Multi_listD), Intent(inout) :: PD(0:)  ! (0:n0) Cluster membership data to be computed
                                                ! 0     (Cluster set data):
                                                !   %k       : success / error code
                                                !   %n       : final # regular clusters = 'ncl'
                                                !   %lt(ncl) : original numbering of clusters if # clusters is reduced 
                                                !              (no change in order)
                                                !   %wt(ns)  : For each slate: sum of initial memberships over all clusters
                                                !              to determine the fuzzy reduction factor
                                                !   %M0(ns,0:1): Fuzzy membership reduction factors (output)
                                                !               0: scale factor
                                                !               1: derivative factor

                                                ! 1:ncl (Individual cluster data):
                                                !   %k  = # fuzzy reduced slate memberships
                                                !   %n  = # slate ballots which are partial or full members
                                                !   %o       : 1 for regular clusters, 2 for independents
                                                !   %fsx     : slate ballot weight averaged fuzzy reduction factor
                                                !   %fux     : cluster width
                                                !   %sum_wt  : total weight of the cluster = Sum(%wt)

                                                !   %px(nc)  : variance of %sx
                                                !   %sx(nc)  : cluster mean vector (over slate ballot 'sl' vectors 
                                                !              Memb(sl)%sx weighted by Memb(sl)%ux * %rx(sl))
                                                !   %tx(0:nc): 1...nc = centered and normalized %sx.  0 = the norm
                                                !   %ux(nc)  : sum of weighted memberhsips of the slate ballots 'sl' based on 
                                                !              slate ballot variance normalization = Sum(Memb(sl)%ux * %rx(sl))

                                                !   %ls(n)   : list of slate ballot members
                                                !   %rx(n)   : unweighted memberhsips of the slate ballots
                                                !   %wt(n)   : weighted memberhsips of the slate ballots = Memb(sl)%fsx*%rx(sl)

                                                !   %M0(n,0:2): Slate ballot membership data
                                                !               0: initial membership after the cosine cutoff
                                                !               1: derivative of this membership wrt the modified correlation
                                                !               2: initial membership * fuzzy derivative factor
                                                !   %M2(nc,n) : For Jacobian computation: Partial derivatives of the initial
                                                !               memberships of the slate ballots in this cluster 
                                                !               with respect to the mean vector of the cluster,
                                                !               ordered by decreasing correlation with the slate ballots
    Real(Dblp),  Intent(out) :: rate_out(:,:)   ! (nc,n0) Final cluster mean vectors PD(k)%sx, for nu nc*n0
    Real(Dblp),  Intent(out) :: change          ! Change in the cluster mean vectors |rate_out - rate_in| / nu
! Local:
    Real(Dblp) :: sm(nc), wt(n0), sx1(n0), tx(0:nc,n0)
    Integer    :: i, k, m, n, cl, nl, np, nu, sl, ier, ncl, lst(n0)

    m= Ubound(PD,1);  n= Size(rate_out,2);  np= Memb(0)%l
    
    If (m < n0 .or. n < n0) then
      Call Out ("Error in 'Generate_clusters': # data clusters",n0, ln=1)
      Call Out ("but 'PD' can record only",n, "and 'rate_out' only",n)
      PD(0)%k= -2;  PD(0)%n= 0;  Return
    End if

    PD(0)%lt= "ID";  rate_out= -1;  change= -1
    tx(1:,:)= rate_in
    If (GN_rating < 1) tx(1:,:)= tx(1:,:) - GN_rate_adj
    Call Normalize_vec (tx)
    
    Do k= 1,n0
      PD(k)%sx= rate_in(:,k);  PD(k)%tx= tx(:,k)
      
      If (Associated(PD(k)%ls)) DeAllocate(PD(k)%ls)
      If (Associated(PD(k)%rx)) DeAllocate(PD(k)%rx)
      If (Associated(PD(k)%wt)) DeAllocate(PD(k)%wt)
      If (Associated(PD(k)%M0)) DeAllocate(PD(k)%M0)
      If (Associated(PD(k)%M2)) DeAllocate(PD(k)%M2)
    End do

    Call Generate_clusters0 (nc,n0,ns, Memb, PD)

    ncl= PD(0)%n;  ier= PD(0)%k
    If (ncl < 1 .or. ier < 0) then
      Call Out ("Error condition in 'Generate_clusters': No more clusters")
      PD(0)%k= Min(ier,-1); Return
    End if
      
!   Eliminate any clusters whose top rating is not sufficiently positive

    Forall(k=1:ncl) sx1(k)= PD(k)%tx(0) * Maxval(PD(k)%tx(1:))
    Call List_of_true (sx1 > GN_min_mx, nl,lst)

    If (nl < ncl) then
      If (ncl < 1) then
        Call Out ("Error condition in 'Generate_clusters': No more clusters")
        PD(0)%k= -1;  Return
      End if

      Call Out ("# regular clusters reduced from",ncl, "to",nl)

      Do k= 1,nl
        cl= lst(k);  PD(0)%lt(k)= PD(0)%lt(cl)
        PD(k)%sx= PD(cl)%sx;  PD(k)%tx= PD(cl)%tx
      End do

      Call Generate_clusters0 (nc,nl,ns, Memb, PD(0:nl))

      ncl= PD(0)%n
      If (ncl < nl) then
        Call Out ("Error condition in 'Generate_clusters': # clusters reduced from",nl, &
                  "to",ncl,ln=1)
        PD(0)%k= -1;  Return
      End if
    End if
      
!   Cluster sigmas and special counts
    
    Do k= 1,ncl
      sm= 0
      Do i= 1,PD(k)%n
        sl= PD(k)%ls(i)
        sm= sm + (Memb(sl)%ux * PD(k)%rx(i)) * (Memb(sl)%px + (Memb(sl)%sx - PD(k)%sx)**2)
      End do
      PD(k)%px= sm / PD(k)%ux
      
      PD(k)%k  = Count(PD(0)%M0(PD(k)%ls,0) < 0.99)
      PD(k)%fsx= Sum(PD(0)%M0(:,0) * Memb(1:)%fsx) / np

      PD(k)%o= 1
    End do
    
!   Compute the final %sx and %tx, then change
    
    Forall(k=1:ncl) rate_out(:,k)= PD(k)%sx
    If (ncl < n0) rate_out(:,ncl+1:)= -1
    Call Limit_ratings (rate_out(:,:ncl), Exceed)

    If (Any(Exceed)) then 
      Forall(k=1:ncl) PD(k)%sx= rate_out(:,k)
      tx(1:,:ncl)= rate_out(:,:ncl)
      tx(1:,:ncl)= tx(1:,:ncl) - GN_rate_adj
      Call Normalize_vec (tx(:,:ncl))
      Forall(k=1:ncl) PD(k)%tx= tx(:,k)
    End if
    
    tx(1:,:ncl)= rate_out(:,:ncl) - rate_in(:,PD(0)%lt(:ncl))
    nu= nc*ncl
    change= Sqrt(Sum(tx(1:,:ncl)**2)/nu) ! Functional= mean Euclidean dif
                                         ! between initial & final mean
                                         ! cluster vectors
  End Subroutine Generate_clusters
  
  Subroutine Generate_clusters0 (nc,n0,ns, Memb, PD)

!   Generate a new set of clusters from a prior set 
!   represented by their  centered and normalized mean 
!   rating vectors 'tx'.
  
!   Assumes prior call to Allocate_clust or Allocate_PD 
  
    Integer,           Intent(in) :: nc  ! # candidates
    Integer,           Intent(in) :: n0  ! # clusters
    Integer,           Intent(in) :: ns  ! # slate ballots
    
    Type(Multi_listD), Intent(in) :: Memb(0:)   ! (0:ns)   Slate ballot data
    
    Type(Multi_listD), Intent(inout) :: PD(0:)  ! (0:n0) Cluster membership data to be computed
                                                ! 0     (Cluster set data):
                                                !   %k       : success / error code
                                                !   %n       : final # regular clusters = 'ncl'
                                                !   %lt(ncl) : original numbering of clusters if # clusters is reduced 
                                                !              (no change in order)
                                                !   %wt(ns)  : For each slate: sum of initial memberships over all clusters
                                                !              to determine the fuzzy reduction factor
                                                !   %M0(ns,0:1): Fuzzy membership reduction factors (output)
                                                !               0: scale factor
                                                !               1: derivative factor

                                                ! 1:ncl (Individual cluster data):
                                                !   %k  = # fuzzy reduced slate memberships
                                                !   %n  = # slate ballots which are partial or full members
                                                !   %o       : 1 for regular clusters, 2 for independents
                                                !   %fsx     : slate ballot weight averaged fuzzy reduction factor
                                                !   %fux     : cluster width
                                                !   %sum_wt  : total weight of the cluster = Sum(%wt)

                                                !   %sx(nc)  : cluster mean vector (over slate ballot 'sl' vectors 
                                                !              Memb(sl)%sx weighted by Memb(sl)%ux * %rx(sl))
                                                !   %px(nc)  : variance of %sx
                                                !   %tx(0:nc): 1...nc = centered and normalized %sx.  0 = the norm
                                                !   %ux(nc)  : sum of weighted memberhsips of the slate ballots 'sl' based on 
                                                !              slate ballot variance normalization = Sum(Memb(sl)%ux * %rx(sl))

                                                !   %ls(n)   : list of slate ballot members
                                                !   %rx(n)   : unweighted memberhsips of the slate ballots
                                                !   %wt(n)   : weighted memberhsips of the slate ballots = Memb(sl)%fsx*%rx(sl)

                                                !   %M0(n,0:2): Slate ballot membership data
                                                !               0: initial membership after the cosine cutoff
                                                !               1: derivative of this membership wrt the modified correlation
                                                !               2: initial membership * fuzzy derivative factor
                                                !   %M2(nc,n) : For Jacobian computation: Partial derivatives of the initial
                                                !               memberships of the slate ballots in this cluster 
                                                !               with respect to the mean vector of the cluster,
                                                !               ordered by decreasing correlation with the slate ballots

                                                ! Key input: PD(0)%lt, PD(1:)%sx, PD(1:)%tx
! Local:
    Real(Dblp),   Parameter :: eps0= 1.0E-7, eps1= 0.001
    Logical                 :: Test_err= .false.
    
    Real(Dblp), Allocatable :: corr(:), mbr(:,:), cor1(:), cor2(:), ql(:,:)
    Integer,    Allocatable :: key(:), key2(:)
    
    Type(AdjacencyD) :: Zer(ns)  ! Normalized mean cluster vector PD(k)%tx, modified
                                 ! by its slate ballot correlations based on the reduction
                                 ! of negative*negative terms by scaling by GN_dot_fac

    Real(Dblp) :: tmp(nc), sm(nc), sm0(0:nc), sm1(nc), cor0(ns)
    Real(Dblp) :: sum_wt, ind_size, reg_size, tot_size
    Integer    :: ls(ns), lst(n0)
    Real       :: tot_wt
    Integer    :: i, j, k, l, m, n, p, cl, sl, nl, np, ncl, ncr

    np= Memb(0)%l
    Do sl= 1,ns
      Allocate(Zer(sl)%wt(0:nc),Zer(sl)%ls(nc))  
    End do

    PD(0)%n= n0;  PD(0)%k= 1;  cl= 0;  lst= -1
    
    Gen1 : Do k= 1,n0
    
!     Compute the initial correlations 'cor0' between the prior clusters and the slate clusters
    
      Do sl= 1,ns
        Call Dot_product_M1 (GN_dot_fac, PD(k)%tx,Memb(sl)%tx(1:), &
                             Zer(sl)%wt,Zer(sl)%ls, Zer(sl)%nl,cor0(sl))
      End do
      
      Call List_of_true (cor0 > Memb(0)%qx(1), ncr,ls)

      If (ncr < 1) then
        Call Out ("Warning in 'Generate_clusters0'. No members for cluster",k)
        Cycle Gen1
      End if
      
      If (Allocated(ql)) DeAllocate(ql, corr, key, key2, mbr)
      Allocate(ql(nc,ncr),corr(ncr),key(ncr),key2(ncr),mbr(ncr,0:1))

      key= ls(:ncr);  ql= 0

!     Compute partial derivatives of the initial memberships wrt
!     the mean vector of the cluster: To be used for the Jacobian
      
      Do i= 1,ncr
        sl= key(i);  n= Zer(sl)%nl
        If (n > 0) then   ! Scale where both correlation factors are negative
          ql(:,i)= (Memb(sl)%tx(1:) - cor0(sl)*Zer(sl)%wt(1:)) / Zer(sl)%wt(0)
          If (GN_dot_fac > 0) then
            ql(Zer(sl)%ls(:n),i)= GN_dot_fac*ql(Zer(sl)%ls(:n),i) 
          Else
            ql(Zer(sl)%ls(:n),i)= 0
          End if
        Else
          ql(:,i)= (Memb(sl)%tx(1:) - cor0(sl)*PD(k)%tx(1:)) / PD(k)%tx(0)  
        End if
      End do
      
      If (Associated(PD(k)%M2)) DeAllocate(PD(k)%M2);  Allocate(PD(k)%M2(nc,ncr))
      PD(k)%M2= ql;  corr= cor0(key)

      Call Sort (.false., corr, key2, ez=eps0)
      PD(k)%M2= PD(k)%M2(:,key2)

      key= key(key2)
      
      Do i= 1,ncr
        Call Cos_rise_drv (0.0d0,Memb(0)%qx(2), corr(i), mbr(i,:))
        PD(k)%M2(:,i)= mbr(i,1) * PD(k)%M2(:,i)
      End do
      sum_wt= Sum(mbr(:,0)*Memb(key)%fsx)
      
      If (pr_out > 1.5) then
         Call Out ("For cluster",k, "with weight",Real(sum_wt),ln=1) 
         Call Out ("List of correlated slate ballots", key) 
         Call Out ("Corresponding correlations in decreasing order", Real(corr))
         Call Out (-1,"The memberships of these slate ballots", Real(mbr))
         Call Out (-1,"The partial derivatives", Real(PD(k)%M2))
      End if    
        
!     Remove clusters that are too small
      
      If (sum_wt <= GN_min_wt) then
        Call Out ("Warning in 'Generate_clusters0'. Delete cluster",k, &
                  "weight too small",Real(sum_wt))
        DeAllocate(PD(k)%M2);  Cycle Gen1
      End if
      
      cl= cl + 1;  lst(cl)= k;  PD(cl)%sum_wt= sum_wt
      If (cl < k) then
        PD(cl)%sx= PD(k)%sx;  PD(cl)%tx= PD(k)%tx
        Allocate (PD(cl)%M2(nc,ncr))
        PD(cl)%M2= PD(k)%M2;  DeAllocate(PD(k)%M2)
      End if
      
      Allocate (PD(cl)%ls(ncr), PD(cl)%M0(ncr,0:2))
      PD(cl)%n= ncr;  PD(cl)%ls= key;  PD(cl)%M0(:,:1)= mbr
    End do Gen1
    
    Call DeAlloc_Adjacency (Zer)
    
    PD(0)%n= cl
    If (cl < 1) then
      Call Out ("Error in 'Generate_clusters0': no clusters left");  
      PD(0)%k= -1;  Return
    End if
    ncl= PD(0)%n;  PD(0)%lt(:ncl)= PD(0)%lt(lst(:ncl))
    
!   Enforce the fuzzy membership requirement and compute fuzzy membership derivatives
    
    PD(0)%wt= 0.0d0
    Do k= 1,ncl
      PD(0)%wt(PD(k)%ls)= PD(0)%wt(PD(k)%ls) + PD(k)%M0(:,0)
    End do

    Do sl= 1,ns
      Call Deriv_scal (PD(0)%wt(sl), PD(0)%M0(sl,:1))
      PD(0)%wt(sl)= PD(0)%M0(sl,0) * PD(0)%wt(sl)
    End do

!   Update the cluster mean vectors
    
    Do k= 1,ncl
      n= PD(k)%n;  Allocate (PD(k)%rx(n), PD(k)%wt(n))
      PD(k)%rx= PD(0)%M0(PD(k)%ls,0) * PD(k)%M0(:,0) ! Apply the fuzzy scale factors to the 
                                                     ! cluster slate ballot memberships
      PD(k)%wt= Memb(PD(k)%ls)%fsx * PD(k)%rx   ! Resulting weights of the slate ballot memberships
      PD(k)%sum_wt= Sum(PD(k)%wt)               ! Total cluster weight
    End do

!   Check size of indepdendents

    If (Test_err) then
      ind_size= Sum(Memb(1:)%fsx * (1 - PD(0)%wt))
      reg_size= Sum(PD(:ncl)%sum_wt)
      tot_size= reg_size + ind_size
      tot_wt= Memb(0)%fsx

      If (Abs(tot_size - tot_wt) > eps1) then
        Call Out ("Cluster size error in 'Generate_clusters0'")
        Call Out ("total cluster size",Real(tot_size), &
        "vs slate cluster total",tot_wt)
      End if
    End if

    Do k= 1,ncl
      sm= 0;  sm1= 0
      Do i= 1,PD(k)%n
        sl= PD(k)%ls(i)
        tmp= PD(k)%wt(i) / Memb(sl)%px
        sm = sm  + tmp
        sm1= sm1 + tmp * Memb(sl)%sx

        PD(k)%M0(i,2)= PD(0)%M0(sl,1) * PD(k)%M0(i,0)
      End do
      
      PD(k)%ux= sm  
      PD(k)%sx= sm1 / sm  ! Mean rating vector of the cluster over its slate ballot members,                    
                          ! using variance normalized weights for the slate ballots
      PD(k)%tx(1:)= PD(k)%sx - GN_rate_adj
      Call Normalize_vec (PD(k)%tx)

      If (pr_out > 1.5) then
         Call Out ("For cluster",k, "with total weight",Real(PD(k)%sum_wt),ln=1) 
         Call Out ("'rx' fuzzy initial memberships", Real(PD(k)%rx))
         Call Out ("'wt' weighted memberships", Real(PD(k)%wt))
         Call Out ("'ux' normalization factors", Real(PD(k)%ux))
         Call Out ("'sx' mean rating vectors", Real(PD(k)%sx))
         Call Out ("'tx' adjusted 'sx'", Real(PD(k)%tx(1:)))
      End if
    End do

  End subroutine Generate_clusters0

  Subroutine Generate_clusters1 (nc,ncl,ns, Memb, rate_in, PD, rate_out)

!   Generate a new set of clusters from a prior set of mean rating vectors %sx
!   or their centered and normalized versions %tx

!   Used for a finite difference Jacobian computation: 'Diff_Jacobian'
!   For this purpose, assume no reductions in # clusters.
  
!   Assumes prior call to Allocate_clust or Allocate_PD 
  
    Integer,           Intent(in) :: nc  ! # candidates
    Integer,           Intent(in) :: ncl ! # clusters
    Integer,           Intent(in) :: ns  ! # slate ballots
    
    Type(Multi_listD), Intent(in) :: Memb(0:)     ! (0:ns) Slate ballot data

    Real(Dblp),        Intent(in) :: rate_in(:,:) ! (nc,ncl) Initial cluster mean vectors %sx
    
    Type(Multi_listD), Intent(inout) :: PD(0:)  ! (0:n0) Cluster membership data to be computed
                                                ! 0     (Cluster set data):
                                                !   %k       : success / error code
                                                !   %n       : final # regular clusters = 'ncl'
                                                !   %lt(ncl) : original numbering of clusters if # clusters is reduced 
                                                !              (no change in order)
                                                !   %wt(ns)  : For each slate: sum of initial memberships over all clusters
                                                !              to determine the fuzzy reduction factor
                                                !   %M0(ns,0:1): Fuzzy membership reduction factors (output)
                                                !               0: scale factor
                                                !               1: derivative factor

                                                ! 1:ncl (Individual cluster data):
                                                !   %k  = # fuzzy reduced slate memberships
                                                !   %n  = # slate ballots which are partial or full members
                                                !   %o       : 1 for regular clusters, 2 for independents
                                                !   %fsx     : slate ballot weight averaged fuzzy reduction factor
                                                !   %fux     : cluster width
                                                !   %sum_wt  : total weight of the cluster = Sum(%wt)

                                                !   %sx(nc)  : cluster mean vector (over slate ballot 'sl' vectors 
                                                !              Memb(sl)%sx weighted by Memb(sl)%ux * %rx(sl))
                                                !   %px(nc)  : variance of %sx
                                                !   %tx(0:nc): 1...nc = centered and normalized %sx.  0 = the norm
                                                !   %ux(nc)  : sum of weighted memberhsips of the slate ballots 'sl' based on 
                                                !              slate ballot variance normalization = Sum(Memb(sl)%ux * %rx(sl))

                                                !   %ls(n)   : list of slate ballot members
                                                !   %rx(n)   : unweighted memberhsips of the slate ballots
                                                !   %wt(n)   : weighted memberhsips of the slate ballots = Memb(sl)%fsx*%rx(sl)

                                                !   %M0(n,0:2): Slate ballot membership data
                                                !               0: initial membership after the cosine cutoff
                                                !               1: derivative of this membership wrt the modified correlation
                                                !               2: initial membership * fuzzy derivative factor
                                                !   %M2(nc,n) : For Jacobian computation: Partial derivatives of the initial
                                                !               memberships of the slate ballots in this cluster 
                                                !               with respect to the mean vector of the cluster,
                                                !               ordered by decreasing correlation with the slate ballots

                                                ! Key input: PD(0)%lt, PD(1:)%sx, PD(1:)%tx
    Real(Dblp),      Intent(out) :: rate_out(:,:) ! (nc,ncl) Final cluster mean vectors %sx
! Local:
    Real(Dblp),   Parameter :: eps0= 1.0E-7
    Real(Dblp), Allocatable :: corr(:), mbr(:,:)
    Integer,    Allocatable :: key(:), key2(:)
    
    Type(AdjacencyD) :: Zer(ns)  ! Normalized mean cluster vector PD(k)%tx, modified
                                 ! by its slate ballot correlations based on the reduction
                                 ! of negative*negative terms by scaling by GN_dot_fac

    Real(Dblp) :: sum_wt, tmp(nc), sm(nc), sm1(nc), cor0(ns), tx(0:nc,ncl)
    Integer    :: ls(ns)
    Integer    :: i, k, n, sl,  ncr

    PD(0)%k= 1;  PD(0)%n= ncl
    Do sl= 1,ns
      Allocate(Zer(sl)%wt(0:nc),Zer(sl)%ls(nc))  
    End do

    tx(1:,:)= rate_in
    If (GN_rating < 1) tx(1:,:)= tx(1:,:) - GN_rate_adj
    Call Normalize_vec (tx)
    
    Gen1 : Do k= 1,ncl

!     Compute the initial correlations 'cor0' between the prior clusters and the slate clusters
    
      Do sl= 1,ns
        Call Dot_product_M1 (GN_dot_fac, tx(:,k),Memb(sl)%tx(1:), &
                            Zer(sl)%wt,Zer(sl)%ls, Zer(sl)%nl,cor0(sl))
      End do
      
      Call List_of_true (cor0 > Memb(0)%qx(1), ncr,ls)

      If (ncr < 1) then
        Call Out ("Warning in 'Generate_clusters1'. No correlated slate ballots")
        Cycle Gen1
      End if
      
      If (Allocated(corr)) DeAllocate(corr, key, key2, mbr)
      Allocate(corr(ncr),key(ncr),key2(ncr),mbr(ncr,0:1))

      key= ls(:ncr);  corr= cor0(key)
      Call Sort (.false., corr, key2, ez=eps0)
      key= key(key2)
      
      Do i= 1,ncr
        Call Cos_rise_drv (0.0d0,Memb(0)%qx(2), corr(i), mbr(i,:))
      End do

      If (pr_out > 1.5) then
         Call Out ("For cluster",k, "with weight",Real(sum_wt),ln=1) 
         Call Out ("List of correlated slate ballots", key) 
         Call Out ("The corresponding correlations in decreasing order", Real(corr))
         Call Out (-1,"The memberships & derivatives of these slate ballots", Real(mbr))
      End if    
        
      If (Associated(PD(k)%ls)) DeAllocate(PD(k)%ls, PD(k)%M0)
      Allocate (PD(k)%ls(ncr), PD(k)%M0(ncr,0:1))
      PD(k)%n= ncr;  PD(k)%ls= key;  PD(k)%M0= mbr
      
    End do Gen1
    
    Call DeAlloc_Adjacency (Zer)

!   Enforce the fuzzy membership requirement and compute fuzzy membership derivatives
    
    PD(0)%wt= 0.0d0
    Do k= 1,ncl
      PD(0)%wt(PD(k)%ls)= PD(0)%wt(PD(k)%ls) + PD(k)%M0(:,0)
    End do

    Do sl= 1,ns
      Call Deriv_scal (PD(0)%wt(sl), PD(0)%M0(sl,:1))
      PD(0)%wt(sl)= PD(0)%M0(sl,0) * PD(0)%wt(sl)
    End do

!   Update the cluster rating vectors
    
    Do k= 1,ncl
      If (Associated(PD(k)%rx)) DeAllocate(PD(k)%rx, PD(k)%wt)
      n= PD(k)%n;  Allocate (PD(k)%rx(n), PD(k)%wt(n))
      PD(k)%rx= PD(0)%M0(PD(k)%ls,0) * PD(k)%M0(:,0)   ! Apply the fuzzy scale factors to the raw cluster slate ballot memberships
      PD(k)%wt= Memb(PD(k)%ls)%fsx * PD(k)%rx          ! Resulting weights of the slate ballot memberships
      PD(k)%sum_wt= Sum(PD(k)%wt)                      ! Total cluster weight
        
      sm= 0;  sm1= 0
      Do i= 1,PD(k)%n
        sl= PD(k)%ls(i);  tmp= Memb(sl)%ux * PD(k)%rx(i)
        sm = sm  + tmp
        sm1= sm1 + Memb(sl)%sx * tmp
      End do
      
      PD(k)%ux= sm;  PD(k)%sx= sm1 / sm  ! Mean rating vector of the cluster over its slate ballot members,                    
                                         ! using variance normalized weights for the slate ballots
      rate_out(:,k)= PD(k)%sx
    End do

  End subroutine Generate_clusters1
   
   Subroutine Limit_ratings (rate, Exceed)

!    Impose max and min limits on the mean rating / ranking vectors of the clusters
   
     Real(Dblp), Intent(inout) :: rate(:,:)   ! (nc,ncl) Cluster mean rating vectors
     Logical,      Intent(out) :: Exceed(2)   ! 1 = max exceeded, 2 = min exceeded
! Local    
     Real(Dblp) :: epsD= 1.0E-7, lim

     lim= GN_max_pt + epsD;  Exceed= .false.
     Exceed(1)= Any(rate > lim)

     If (Exceed(1)) then      ! Some ranking or rating data exceeded max. 
       Where(rate > GN_max_pt) rate= GN_max_pt

       If (pr_out > 1) then
         Call Out ("Warning in 'Limit_ratings': Some rankings or ratings exceeded max")
         Call Out (-1,"Cluster mean vectors",Real(rate))
       End if
     End if
       
     If (GN_rating < 1) then  ! Ranking data. Min= 0
       Exceed(2)= Any(rate < -epsD)

       If (Exceed(2)) then    ! Some ranking data below min
         Where(rate < 0.0d0) rate= 0.0d0

         If (pr_out > 1) then
           Call Out ("Warning in 'Limit_ratings: Some rankings were below min")
           Call Out (-1,"Cluster mean vectors",Real(rate))
         End if
       End if
     Else                      ! Rating data. Min= -GN_max_pt
       Exceed(2)= Any(rate < -lim)

       If (Exceed(2)) then     ! Some rating data below min
         Where(rate < -GN_max_pt) rate= -GN_max_pt

         If (pr_out > 1) then
           Call Out ("Warning in 'Limit_ratings: Some ratings below min")
           Call Out (-1,"Cluster mean vectors",Real(rate))
         End if
       End if
     End if
   End Subroutine Limit_ratings

  Pure Subroutine Deriv_scal (mb_sum, mb_fac)

!   Compute scale factor 'mb_fac' to scale down slate cluster memberships 
!   in a regular cluster when their sum 'mb_sum' over the regular clusters
!   is more than Sf1, also possibly its derivative and second derivative
!   with respect to the membership sum 'mb_sum'.

!   This keeps the total weight < 1 in accordance with fuzzy membership
!   but with a smooth cutoff. When mb_sum = 1, it is reduced to 0.9632
!   by the scale factor mb_fac(0), leaving 0.0368 of the slate cluster
!   weight to represent the independents.

    Real(Dblp),  Intent(in) :: mb_sum
    Real(Dblp), Intent(out) :: mb_fac(0:)  ! 0 = scale value, 1 = its derivative,
                                           ! 2 = its second derivative

    Real(Dblp), Parameter :: Sf1= 0.90, Sf2= 1 - Sf1
    Real(Dblp) :: ds
    Integer    :: n
    
    n= Ubound(mb_fac,1);  mb_fac= 0
    
    If (mb_sum <= Sf1) then
      mb_fac(0)= 1
    Else 
      ds= Exp((Sf1 - mb_sum) / Sf2)
      mb_fac(0)= (1 - Sf2 * ds) / mb_sum
      
      If (n > 0) then
        mb_fac(1)= (ds - mb_fac(0)) / mb_sum
        If (n > 1) mb_fac(2)= -(ds/Sf2 + 2 * mb_fac(1)) / mb_sum
      End if
    End if
  End Subroutine Deriv_scal
  
End Module Clust_gen