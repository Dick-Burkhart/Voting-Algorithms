    
Module Clusters_support

   Use Clusters0

   Use Newton_operators   
   Use Graph_algorithms
   Use Cholesk
   Use Util
   Use Output
   Use Types
   Use Precisn
   Use IEEE_Arithmetic
   Implicit None
  
   Interface Portions
     Module procedure PortionsR, PortionsD
   End Interface

   Interface Objective_values
     Module procedure Objective_valuesR, Objective_valuesD
   End Interface

   Interface Zrf_M
     Module procedure Zrf_MR, Zrf_MD
   End Interface
 Contains
    
   Subroutine Var_siz_record (unit,term, n,line)

!    Read from 'unit' a positive integer record 'line' of possibly unknown length 'n',
!    which is terminated by either an integer less than or equal to 'term'  
!    or by an end-of-record condition.

     Integer,     Intent(in) :: unit    ! File unit #
     Integer,     Intent(in) :: term    ! Line terminating value
     Integer,    Intent(out) :: n       ! # values read
     Integer,  Intent(inout) :: line(:) ! List of those values
! Local:
     Integer :: i, ios

     Do n= 1,Size(line)
       Read (unit,'(I6)', Advance='No', Iostat=ios, Eor=100) line(n)
       If (line(n) <= term) Exit
     End do
     
     n= n - 1;  Read (unit,*)
     Return
     
100  If (term >= 0) n= n - 1  ! Forced to do line advance: 
                              ! Nothing more to read.
   End Subroutine Var_siz_record


  Subroutine Indie_membership (nc,ncl, mb_sum,Memb, Clust)
  
!   Compute the membership data for the cluster of independent voters.
  
!   Compute: Clust%n, Clust%ls, Clust%wt, Clust%sum_wt, 
!            Clust%rx, Clust%px, Clust%sx
  
    Integer,    Intent(in) :: nc        ! # candidates
    Integer,    Intent(in) :: ncl       ! # regular clusters
    Real(Dblp), Intent(in) :: mb_sum(:) ! (ns) Fuzzified sum of 
                                        !   unweighted memberships of 
                                        !   all regular clusters, all < 1
    Type(Multi_listD),    Intent(in) :: Memb(0:) ! (0:ns)  Slate clusters
    Type(Multi_listD), Intent(inout) :: Clust    ! "Cluster" of independents
! Local:
    Real(Dblp) :: zx(nc), sm(nc), sm1(nc), sm2(nc)
    Real       :: sz
    Integer    :: sl, ns, ni
    
    ns= Ubound(Memb,1);  Clust%k= 1;  Clust%n= ns;  Clust%o= 2

    If (.not.Associated(Clust%ls)) then
      Allocate(Clust%ls(ns), Clust%wt(ns), Clust%rx(ns))
    Else 
      ni= Size(Clust%ls)
      If (ni /= ns) then
        DeAllocate(Clust%ls,Clust%wt,Clust%rx)
        Allocate(Clust%ls(ns), Clust%wt(ns), Clust%rx(ns))
      End if
    End if
    
    Clust%ls= "ID";  Clust%rx= 1 - mb_sum
    Clust%wt= Clust%rx * Memb(1:)%fsx

    Clust%l= Count(Clust%rx > 0.333)
    Clust%m= Count(Clust%rx > 0.99)
    Clust%sum_wt= Sum(Clust%wt);  sz= Clust%sum_wt

    If (sz < Min_wti) then ! Return to reduce # regular clusters
      if (pr_out > 1.5) Call Out ("Warning in 'Indie_membership': Size too small",sz)
      Clust%o= 0;  Return
    End if
    
    sm= 0;  sm1= 0;  sm2= 0
    Do sl= 1,ns
        zx= Clust%rx(sl) * Memb(sl)%ux
      sm = sm  + zx
      sm1= sm1 + zx * Memb(sl)%sx
    End do
    Clust%sx= sm1 / sm
    
    Do sl= 1,ns
        zx= Clust%rx(sl) * Memb(sl)%ux
      sm2= sm2 + zx * (Memb(sl)%px + (Memb(sl)%sx - Clust%sx)**2)
    End do
    Clust%px= sm2 / sm

  End Subroutine Indie_membership

  
  Subroutine Cluster_set_width (nc, Clust)
  
!   Compute cluster width %fux.
  
!   Use:      Clust(0)%L1, Clust%px, Clust%sx
!   Compute:  Clust%fux
  
    Integer,              Intent(in) :: nc         ! # candidates
    Type(Multi_listD), Intent(inout) :: Clust(0:)  ! (ind) Cluster data
! Local:
    Real(Dblp), Parameter :: del= 0.10
    Real(Dblp) :: wt(nc), var(nc)
    Integer    :: n, cl, key(nc)
 
!   Compute the cluster width %fux as the square root of the harmonic mean 
!   of the variances of top candidates
    
    Do cl= 1,Ubound(Clust,1)
      n= Clust(0)%L1(-2,cl) + 1 

      If (Rating < 1) then
        key(:n)= Clust(0)%L1(1:n,cl)
        wt(:n) = Clust(cl)%sx(key(:n)) / Max_ptD + del
      Else
        wt= Abs(Clust(cl)%sx);  Call Sort (.false.,wt, key)  
        wt(:n)= wt(key(:n)) / Max_ptD + del
      End if

      var(:n)= Max(Clust(cl)%px(key(:n)), del)

      Clust(cl)%fux= Sqrt( Sum(wt(:n))/ Sum(wt(:n)/ var(:n)) )  ! Weighted harmonic mean
      Clust(0)%M3(6,cl)= Clust(cl)%fux
    End do
    
  End Subroutine Cluster_set_width

    
  Subroutine ReOrd_clusters (ncl, Clust, key, ReOrd)
  
!   Reorder the regular clusters of a cluster set by decreasing cluster size
  
    Integer,              Intent(in) :: ncl       ! # regular clusters
    Type(Multi_listD), Intent(inout) :: Clust(:)  ! (ncl) Cluster set data
    Integer,             Intent(out) :: key(:)    ! (ncl) Reording key
    Logical,             Intent(out) :: ReOrd
! Local:
    Type(Multi_listD) :: Tmp(ncl)    
    Real    :: siz(ncl)
    Integer :: ind, lst(ncl)
    
    key= "ID";  ReOrd= .false.;  If (ncl <= 1) Return
    
    siz= Clust(:ncl)%sum_wt;  Call Sort (.false., siz, key,ReOrd)
    
    If (ReOrd) then
      Call Copy_multi_list (key,Clust(:ncl), Tmp)
      lst= "ID";  Call Copy_multi_list (lst, Tmp, Clust(:ncl))
      Call DeAlloc_Multi_list_ar (Tmp)
    End if
    
  End Subroutine ReOrd_clusters
          

  Subroutine Average_consolidated (n1,n2,nq, sing,pair,prc, Cliq_set, vec_in, map,vec_out)

!   Do a weighted average of 'vec_in' (initial cluster mean rating vectors)
!   over each clique to get 'vec_out' merged cluster mean rating vectors.

    Integer,  Intent(in) :: n1             ! # singleton components
    Integer,  Intent(in) :: n2             ! # pair components
    Integer,  Intent(in) :: nq             ! # cliques from components of size >= 3
    
    Integer,  Intent(in) :: sing(:)        ! (n1)   Singleton components
    Integer,  Intent(in) :: pair(:,:)      ! (2,n2) Pair components
    Real,     Intent(in) :: prc(:)         ! (n2) Connectivity of each pair component
    Type(Set_list), Intent(in) :: Cliq_set(:)  ! (nq) Disjoint cliques of correlated clusters
                                           !   %n = # clusters in the clique
                                           !   %set(n) = list of these clusters
                                           !   %val(n) = corresponding cluster weights
    
    Real,     Intent(in) :: vec_in(0:,:)   ! (0:nc,ncl) 0 = cluster size
                                           !            1:nc = Initial cluster mean vectors 
    
    Integer, Intent(out) :: map(:,:)       ! (ncl,2) 1 = Mapping from initial clusters to final clusters:
                                           !             first singletons, then pairs, then cliques
                                           !         2 = For each of the 'nl' final clusters, the largest 
                                           !             initial cluster mapped to it
    Real,    Intent(out) :: vec_out(0:,:)  ! (0:nc,nl) 1:nc = merged cluster mean vectors, 0 = cluster size
! Local:
    Real,    Allocatable :: szq(:)  ! Sizes of the clusters of cliques 
    Integer, Allocatable :: clq(:)  ! The corresponding clusters
    Real    :: sm       ! Sum of the clique sizes
    Real    :: sz(2)    ! Sizes of the 2 initial clusters of a pair
    Integer :: j, k, n, cl
    
    map= 0;  vec_out= -1

!   Singletons

    If (n1 > 0) then
      vec_out(:,:n1)= vec_in(:,sing)
      Forall(j=1:n1) map(sing(j),1)= j;  map(:n1,2)= sing
    End if
    
!   Pairs

    k= n1
    Do cl= 1,n2
      k= k + 1;  sz= vec_in(0,pair(:,cl))
      sm= Sum(sz);  sz= sz/sm
      vec_out(0,k) = sm
      vec_out(1:,k)= Matmul(vec_in(1:,pair(:,cl)), sz)

      map(pair(:,cl),1)= k;  j= Maxloc(sz,1);  map(k,2)= pair(j,cl)

      If (pr_out > 1.5) Call Out ("Merged pair of clusters", &
                                   pair(1,cl), "&", pair(2,cl))
    End do
    
!   Cliques

    Do cl= 1,nq
      k= k + 1;  n= Cliq_set(cl)%n;  Allocate(clq(n), szq(n))

      clq= Cliq_set(cl)%set;  szq= Cliq_set(cl)%val
!      szq= vec_in(0,clq)
      sm= Sum(szq);  szq= szq / sm

      vec_out(0,k) = sm
      vec_out(1:,k)= Matmul(vec_in(1:,clq), szq)

      map(clq,1)= k;  j= Maxloc(szq,1);  map(k,2)= clq(j)

      If (pr_out > 1.5) Call Out ("Merged clique of clusters", clq)
      DeAllocate(clq, szq)
    End do

    If (pr_out > 1.5) then
      Call Out ("Consolidated cluster sizes", vec_out(0,:))
      Call Out (-1, "Consolidated cluster mean vectors", vec_out(1:,:))
      Call Out (-1,"Mapping from original to merged clusters and its inverse", map)
    End if
    
  End Subroutine Average_consolidated
  
  Subroutine STV_clusters (np,nb,STVmb, Clust,key)
  
!   Compute the simplest clusters 'Clust' corresponding to the STV data structure STVmb:
!   Using raw ballot membership, compute the top 'np' clusters, in decreasing order by
!   ballot weight, plus the remaining ballots in the cluster np+1 of independents.

    Integer,        Intent(in) :: np            ! # candidates to be elected
    Integer,        Intent(in) :: nb            ! # ballots
    Type(Set_list), Intent(in) :: STVmb(0:)     ! (0:nb) Ballot ranking structure for clustering
                                                ! Ballot sets (0): 
                                                !   %n        = nc
                                                !   %p        = np
                                                !   %svl      = elected vote = Sum(%val(set(:np)))
                                                !   %smb      = independents vote = Sum(%val(set(np+1:)))
                                                !   %set(nc)  = candidates in elected ordering, with
                                                !               the first 'np' elected
                                                !   %lev(nc)  = candidate clusters by decreasing cluster vote
                                                !   %val(nc)  = candidate cluster votes, not reordered

!                                               ! Ballots (1:nb):
                                                !   %n     = # elected candidates
                                                !   %p     = # ranked candidates
                                                !   %set(p)= 'n' elected candidates, then unelected,
                                                !            in original ranking order
                                                !   %val(p)= corresponding membership weights 
                                                !   %svl   = elected membership weight= Sum(%val(:n))
                                                !   %smb   = total membership weight  = Sum(%val)

    Type(Set_list), Intent(out) :: Clust(:)    ! (np+1) The STV clusters, ordered by decreasing 
                                               !        cluster size for 1:np
                                               !  %n  = # ballot members of the cluster
                                               !  %smb= Cluster size = weight = Sum(Clust%val)
                                               !  %set(n)= List of member ballots by increasing ballot index
                                               !  %val(n)= Weighted memberships of the members
    Integer,        Intent(out) :: key(:)      ! (np) Cluster re-ordering key by size
! Local:
    Type(Set_list) :: Tmp(np)
    Real    :: vl(nb), siz(np)
    Integer :: lb(nb), lst(np)
    Logical :: ReOrd
    Real    :: wb
    Integer :: b, i, j, k, n, cl, n1, ind

    If (pr_out > 1) Call Out ("Enter 'STV_clusters'")

    Elected_clusters : Do i= 1,np
      cl= STVmb(0)%set(i);  k= 0;  lb= 0;  vl= 0

      Bal_loop1 : Do b= 1,nb
        j= First_true(STVmb(b)%set == cl);  If (j <= 0) Cycle Bal_loop1
        
        wb= STVmb(b)%val(j)
        If (wb > 0) then
          k= k + 1;  lb(k)= b;  vl(k)= wb
        End if
      End do Bal_loop1

      n= k;  Clust(i)%n= n
      Allocate(Clust(i)%set(n), Clust(i)%val(n))

      Do k= 1,n
        Clust(i)%set(k)= lb(k);  Clust(i)%val(k)= vl(k)
      End do

      Clust(i)%smb= Sum(Clust(i)%val)
      siz(i)= Clust(i)%smb
    End do Elected_clusters

    Call Sort (.false.,siz,key, ReOrd)

    If (ReOrd) then
      Call Copy_Set_list (key, Clust(:np), Tmp)
      lst= "ID";  Call Copy_Set_list (lst, Tmp, Clust(:np))
    End if

!   Determine the set of independents

    ind= np + 1;  k= 0;  lb= 0;  vl= 0

    Bal_loop2 : Do b= 1,nb
      n1= STVmb(b)%n + 1;  wb= Sum(STVmb(b)%val(n1:))
      If (wb > 0) then
        k= k + 1;  lb(k)= b;  vl(k)= wb
      End if
    End do Bal_loop2

    Clust(ind)%n= k;  Allocate(Clust(ind)%set(k), Clust(ind)%val(k))
    Clust(ind)%set= lb(:k);  Clust(ind)%val= vl(:k)
    Clust(ind)%smb= Sum(Clust(ind)%val)

    If (pr_out > 1.5) then
      Call Out ("Original candidate ordering by size",STVmb(0)%lev(:np))
      Call Out ("with cluster sizes",STVmb(0)%val(STVmb(0)%lev(:np)))
      Call Out ("The final sizes",Clust(:np)%smb)
      Call Out ("Original independents size",STVmb(0)%smb, "vs final",Clust(np+1)%smb)
    End if

  End Subroutine STV_clusters


  Subroutine Merge_STV (np,nc, Adj,rank_pt,parm, ballot,STV_mrg, &
                        ncl,cls, Clust,Borda,Corr)
  
!   Like 'Merge_n_delete' but uses the generic cluster
!   data structure STV_mrg for merging and deleting clusters.
!   This uses 'Merge_corr', based on cluster correlations,
!   derived from cluster mean ranking vectors, computed as
!   Borda Counts that emphasize the top ranked candidates
!   of each ballot. 

!   The algorithm:
  
!   Order the clusters from smallest to largest by cluster weight, 
!   keeping this ordering throughout.  Starting with the cluster 
!   of least weight, see if it can be merged with the next 
!   heavier cluster, or the next, etc., according to the correlation 
!   parm(1). 

!   After every merger, reorder and restart the merging cycle.
!   If no merger is possible, delete this cluster if it fails 
!   the min size parm(2). Continue until no more merging or
!   deletion is possible or until there are only parm(3)
!   merged clusters left.

    Integer,        Intent(in) :: np           ! # elected candidates = total ballot weight
    Integer,        Intent(in) :: nc           ! # candidates
    Real,           Intent(in) :: Adj          ! Noise adjustment to subtract from mean vectors
                                               !   before correlating
    Real,           Intent(in) :: rank_pt(:)   ! (mr) Point conversion of ranking levels
    Real,           Intent(in) :: parm(:)      ! (3) Algorithm parameters: 
                                               !   (1) = cluster correlation merge criterion
                                               !   (2) = cluster deletion criterion
                                               !   (3) = minimum number of clusters upon output
    Integer,        Intent(in) :: ballot(0:,:) ! (0:mr,nb) Ballot data
    Type(Set_list), Intent(in) :: STV_mrg(:)   ! (np+1) Membership data for the STV clusters,
                                               !        ordered by decreasing 1:np cluster size
                                               !  %n  = # ballot members of each cluster
                                               !  %smb= Cluster size = weight = Sum(Clust%val)
                                               !  %set(n)= Member ballots in increasing order
                                               !  %val(n)= Weighted memberships of the ballots

    Integer,        Intent(out) :: ncl         ! # final clusters
    Integer,        Intent(out) :: cls(:)      ! (np) Mapping from elected candidate clusters
                                               !      to the final clusters
    Type(Set_list), Intent(out) :: Clust(:)    ! (np+1=>ind) Final data for the STV clusters
                                               !  %n  = # ballot members of each cluster
                                               !  %smb= Cluster size = weight = Sum(Clust%val)
                                               !  %set(n)= Member ballots in increasing order
                                               !  %val(n)= Weighted memberships of the ballots
    Real,           Intent(out) :: Borda(:,:)  ! (nc,np+1=>ind) Final cluster Borda Counts
    Real,           Intent(out) :: Corr(:,:)   ! (np=>ncl,np=>ncl) Final cluster correlation matrix
                                               !   except = 0 on the diagonal
! Local:
    Real    :: mx_ovr(0:np), mx_cor(0:np), Overlap(np,np) 
    Real    :: wb, wt, cls_vt 
    Integer :: b, k, n, n1, ind, id(np)

    If (pr_out > 1) Call Out ("Enter 'Merge_STV'")

    ncl= np;  cls= "ID";  Clust= STV_mrg;  Borda= 0;  Corr= 0

    Call Cluster_overlap (STV_mrg(:np), Overlap,mx_ovr)

    Call Cluster_corr (nc,Adj,rank_pt, ballot,STV_mrg(:np), &
                       Borda(:,:np), Corr,mx_cor)

    If (pr_out > 1) then
      Call Out ("STV cluster ballot counts",STV_mrg%n)
      Call Out ("STV cluster sizes",STV_mrg%smb)
      Call Out ("Max cluster overlap by row",mx_ovr(1:))
      Call Out (-1,"Cluster overlap matrix",Overlap)
      Call Out ("Max STV cluster correlation by row",mx_cor(1:))
      Call Out (-1,"Cluster correlation matrix",Corr)
      Call Out (-1,"Borda Counts for the initial STV clusters",Borda(:,:np))
    End if
 
    n1= np + 1;  id= "ID"
    Call Copy_set_list (id, STV_mrg(:np), Clust(:np))

    Call Merge_corr (nc,np,Adj,rank_pt, parm,ballot, Clust(:np), &
                     Borda(:,:np),Corr, ncl,cls)
   
!   Note: This assumes that only clusters were merged in Merge_corr, none deleted.
!   Otherwise that would need to be merged with Clust(ind).

    ind= ncl+1;  Call Copy_set (STV_mrg(n1), Clust(ind))
    cls_vt= Sum(Clust(:ind)%smb) 

    wt= 0;  Borda(:,ind)= 0

    Do k= 1,Clust(ind)%n
     b= Clust(ind)%set(k);  wb= Clust(ind)%val(k)
        n= ballot(0,b)
      Borda(ballot(1:n,b),ind)= Borda(ballot(1:n,b),ind) + rank_pt(:n) * wb
      wt= wt + wb
    End Do 
    Borda(:,ind)= Borda(:,ind) / wt

    If (pr_out > 1) then
      Call Out ("After merge: cluster ballot counts", Clust(:ind)%n)
      Call Out ("After merge: cluster sizes", Clust(:ind)%smb)
      Call Out ("with summed cluster vote",cls_vt)
      Call Out ("with initial cluster to final mapping",cls)
      Call Out (-1,"Borda Counts for the merged STV clusters",Borda(:,:ind))
      Call Out (-1,"Correlations between STV clusters",Corr(:ncl,:ncl))
    End if
  End Subroutine Merge_STV

  Subroutine Merge_STV2 (np,nc,ptval,parm, ballot,ballot2, &
                         STV_mrg, ncl,cls, Clust,Borda,Corr, Adj)
  
!   Like 'Merge_n_delete' but uses the generic cluster
!   data structure STV_mrg for merging and deleting clusters.
!   This uses 'Merge_corr2', based on cluster correlations,
!   derived from cluster mean ranking vectors, computed as
!   Borda Counts that emphasize the top ranked candidates
!   of each ballot. 

!   The algorithm:
  
!   Order the clusters from smallest to largest by cluster weight, 
!   keeping this ordering throughout.  Starting with the cluster 
!   of least weight, see if it can be merged with the next 
!   heavier cluster, or the next, etc., according to the correlation 
!   parm(1). 

!   After every merger, reorder and restart the merging cycle.
!   If no merger is possible, delete this cluster if it fails 
!   the min size parm(2). Continue until no more merging or
!   deletion is possible or until there are only parm(3)
!   merged clusters left.

    Integer,  Intent(in) :: np            ! # elected candidates = total ballot weight
    Integer,  Intent(in) :: nc            ! # candidates

    Real,     Intent(in) :: ptval(:)      ! (mt) Point values of ranking or rating levels.
                                          !   If all positive (Rating <= 0) then 'Adj' is expected 
    Real,     Intent(in) :: parm(:)       ! (3) Algorithm parameters: 
                                          !   (1) = cluster correlation merge criterion
                                          !   (2) = cluster deletion criterion
                                          !   (3) = minimum number of clusters upon output
    Integer,  Intent(in) :: ballot(0:,:)  ! (0:mr,nb) Ranked or rated candidates,  with
                                          !    (0,b) = # ranked or rated
    Integer,  Intent(in) :: ballot2(0:,:) ! (0:mr,nb) Corresponding ranking or rating levels,  
                                          !    increasing, with (0,b) = # positive levels

    Type(Set_list),  Intent(in) :: STV_mrg(:)  ! (np+1) Membership data for the STV clusters,
                                               !        ordered by decreasing 1:np cluster size
                                               !  %n  = # ballot members of each cluster
                                               !  %smb= Cluster size = weight = Sum(Clust%val)
                                               !  %set(n)= Member ballots in increasing order
                                               !  %val(n)= Weighted memberships of the ballots

    Integer,        Intent(out) :: ncl         ! # final clusters
    Integer,        Intent(out) :: cls(:)      ! (np) Mapping from elected candidate clusters
                                               !      to the final clusters
    Type(Set_list), Intent(out) :: Clust(:)    ! (np+1=>ind) Final data for the STV clusters
                                               !  %n  = # ballot members of each cluster
                                               !  %smb= Cluster size = weight = Sum(Clust%val)
                                               !  %set(n)= Member ballots in increasing order
                                               !  %val(n)= Weighted memberships of the ballots
    Real,           Intent(out) :: Borda(:,:)  ! (nc,np+1=>ind) Final cluster Borda Counts
    Real,           Intent(out) :: Corr(:,:)   ! (np=>ncl,np=>ncl) Final cluster correlation matrix
                                               !   except = 0 on the diagonal

    Real,  Optional, Intent(in) :: Adj         ! Noise adjustment to subtract from mean vectors
                                               !   before correlating
! Local:
    Real    :: mx_ovr(0:np), mx_cor(0:np), Overlap(np,np) 
    Real    :: wb, wt, cls_vt 
    Integer :: b, k, n, n1, n2, ind, id(np)

    If (pr_out > 1) Call Out ("Enter 'Merge_STV2'")

    ncl= np;  cls= "ID";  Clust= STV_mrg;  Borda= 0;  Corr= -1

    Call Cluster_overlap (STV_mrg(:np), Overlap,mx_ovr)

    If (Present(Adj)) then 
      Call Cluster_corr2 (nc,ptval, ballot,ballot2,STV_mrg(:np), &
                          Borda(:,:np), Corr,mx_cor, Adj)
    Else
      Call Cluster_corr2 (nc,ptval, ballot,ballot2,STV_mrg(:np), &
                          Borda(:,:np), Corr,mx_cor)
    End if

    If (pr_out > 1) then
      Call Out ("STV cluster ballot counts",STV_mrg%n)
      Call Out ("STV cluster sizes",STV_mrg%smb)
      Call Out ("Max cluster overlap by row",mx_ovr(1:))
      Call Out (-1,"Cluster overlap matrix",Overlap)
      Call Out ("Max STV cluster correlation by row",mx_cor(1:))
      Call Out (-1,"Cluster correlation matrix",Corr)
      Call Out (-1,"Borda Counts for the initial STV clusters",Borda(:,:np))
    End if
 
    n1= np + 1;  id= "ID"
    Call Copy_set_list (id, STV_mrg(:np), Clust(:np))

    If (Present(Adj)) then
      Call Merge_corr2 (nc,np,ptval, parm,ballot,ballot2, Clust(:np), &
                        Borda(:,:np),Corr, ncl,cls, Adj)
    Else
      Call Merge_corr2 (nc,np,ptval, parm,ballot,ballot2, Clust(:np), &
                        Borda(:,:np),Corr, ncl,cls)
    End if

!   Note: This assumes that only clusters were merged in Merge_corr2, none deleted.
!   Otherwise that would need to be merged with Clust(ind).

    ind= ncl+1;  Call Copy_set (STV_mrg(n1), Clust(ind))
    cls_vt= Sum(Clust(:ind)%smb) 

    wt= 0;  Borda(:,ind)= 0

    Do k= 1,Clust(ind)%n
      b= Clust(ind)%set(k);  wb= Clust(ind)%val(k);  n= ballot(0,b)
      Borda(ballot(1:n,b),ind)= Borda(ballot(1:n,b),ind) + ptval(ballot2(1:n,b)) * wb
      wt= wt + wb
    End Do 
    Borda(:,ind)= Borda(:,ind) / wt

    If (pr_out > 1) then
      Call Out ("After merge: cluster ballot counts", Clust(:ind)%n)
      Call Out ("After merge: cluster sizes", Clust(:ind)%smb)
      Call Out ("with summed cluster vote",cls_vt)
      Call Out ("with initial cluster to final mapping",cls)
      Call Out (-1,"Borda Counts for the merged STV clusters",Borda(:,:ind))
      Call Out (-1,"Correlations between STV clusters",Corr(:ncl,:ncl))
    End if
  End Subroutine Merge_STV2

  Subroutine Merge_corr (nc,n0,Adj,rank_pt, parm,ballot, Clust,Borda,Corr, ncl,cls)
  
!   Like Merge_matrix but uses the Set_list and 'ballot' data structures.
!   This version typically uses substantiallyh less memory than Merge_matrix 
!   to represent the clusters in terms of the original ballots rather
!   than slate ballots.
  
!   At all intermediate stages 'Clust', 'Borda', and 'Corr' stay the
!   same except when clusters c1 and c2 merge into c2 or clusters are
!   deleted. Full deletion and reordering are done only at the end via 'ord'.

    Integer,  Intent(in) :: nc            ! # candidates
    Integer,  Intent(in) :: n0            ! Initial # clusters
    Real,     Intent(in) :: Adj           ! Noise adjustment to subtract from mean vectors
                                          !   before correlating
    Real,     Intent(in) :: rank_pt(:)    ! (mr) Point conversion of ranking levels
    Real,     Intent(in) :: parm(:)       ! (3) Algorithm parameters: 
                                          !   (1) = cluster correlation merge criterion
                                          !   (2) = cluster deletion criterion
                                          !   (3) = minimum number of clusters upon output
    Integer,  Intent(in) :: ballot(0:,:)  ! (0:mr,nb) Original ballot data 


    Type(Set_list), Intent(inout) :: Clust(:)  ! (n0=>ncl) Membership for the clusters, ordered
                                               !           by decreasing cluster size, with small 
                                               !           clusters deleted and substantially  
                                               !           overlapping clusters merged
                                               !  %n  = # members
                                               !  %smb= Cluster weight = Sum(Clust%val)
                                               !  %set(n)= List of members (increasing on output)
                                               !  %val(n)= Corresponding weighted memberships

    Real,  Intent(inout) :: Borda(:,:)  ! (nc,n0=>ncl) Borda Count for each candidate cluster
    Real,  Intent(inout) :: Corr(:,:)   ! (n0,n0) Correlation matrix between pairs of clusters, 
                                        !   except = 0 on the diagonal

    Integer, Intent(out) :: ncl         ! # final clusters
    Integer, Intent(out) :: cls(:)      ! (n0) Maps original clusters to final clusters, or -1 if none
! Local:
    Real, Parameter :: eps= 0.00001
    Type(Set_list) :: Bal(n0)
    Real,    Allocatable :: vl2(:)
    Integer, Allocatable :: Intr(:), Unyn(:), Ind1(:), Ind2(:)

    Real    :: siz(n0), siz0(n0)
    Integer :: cnt(n0), cnt0(n0), key(n0), inv(n0), lst(n0), ord(n0)
    
    Logical :: ReOrd
    Real    :: cor, tot_wt
    Integer :: i1, j1, k1, k2, mi, mr, mu, n1, n2, ni, nl, nu, ind, min_num, tot_mb
    Integer :: i, j, k, l, m, n, c1, c2, cl

    mr= Ubound(ballot,1);  min_num= Nint(parm(3))
    cls= "ID";  ord= cls

    cnt0= Clust%n;    tot_mb= Sum(cnt0)
    siz0= Clust%smb;  tot_wt= Sum(siz0)

    If (pr_out > 1) then
      Call Out ("Enter 'Merge_corr' with # clusters",n0, &
              "with tot_wt weight",tot_wt, ln=1)
       Call Out ("with cluster sizes",siz0)
    End if
    Call List_of_true (siz0 <= parm(2), n,lst)

    If (n >= 1 .and. n0 <= min_num) then
      If (pr_out > 1) Call Out ("Too few clusters: Delete small clusters with no merging")

      nl= n0  
      Do i= 1,n
        c1= lst(i);  Call Del_mrg (nl,c1,0, Clust%smb,cls,ord, ncl);  Clust(c1)%n= 0
      End do

      ind= ncl + 1;  key(:ncl)= "ID"

      If (Any(ord(:ncl) /= key(:ncl))) then
        Call Copy_set_list (ord(:ncl), Clust, Bal(:ncl))
        Call Copy_set_list (key(:ncl), Bal(:ncl), Clust(:ncl))

        Borda(:,:ncl)  = Borda(:,ord(:ncl))
        Corr(:ncl,:ncl)= Corr(ord(:ncl),ord(:ncl))
      End if

      Call DeAlloc_set_list_ar (Bal(:ncl))
      Call DeAlloc_set_list_ar (Clust(ind:n0))

      Clust(ind:)%n= 0;  Clust(ind:)%smb= 0
      Borda(:,ind:)= 0;  Corr(:,ind:)= 0;  Corr(ind:,:)= 0

      If (pr_out > 1) then
        Call Out ("Reduced # clusters in 'Merge_corr' from",n0, &
                  "to",ncl, ln=1)
        Call Out ("Initial to final cluster mapping", cls)
        Call Out ("Final cluster ordering", ord(:ncl))

        Call Out ("Final cluster sizes",Clust(:ncl)%smb)
        Call Out ("Final cluster membership",Clust(:ncl)%n)
      End if
      Return
    End if
    
!   Merge and delete from smallest to largest cluster

    ncl= n0;  siz= 0;  key= 0;  i= 1 
    ord= Reverse(cls);  cls= ord  ! Reverse the order

    Least_clust : Do while (i <= ncl .and. ncl > min_num)
      c1= ord(i);  i1= i + 1;  j= i1
      If (Clust(c1)%smb < eps) then
        i= i1;  Cycle Least_clust
      End if
      
      Merge_loop : Do while (j <= ncl)
        c2= ord(j)
        If (Clust(c2)%smb < eps) then
          j= j + 1;  Cycle Merge_loop
        End if

        cor= Corr(c1,c2);  n1= Clust(c1)%n;  n2= Clust(c2)%n

        If (cor > parm(1)) then !  Merge correlated cluster c1 with c2 by 
                                !  taking the union, including adding values
          mi= Min(n1,n2);  mu= n1 + n2
          Allocate(vl2(n2), Intr(mi),Unyn(mu), Ind1(n1), Ind2(n2))

          If (pr_out > 1) then
            Call Out ("Merge sets with correlation",cor, ln=1)
            Call Out ("for set c1",c1, "& set c2",c2)
            Call Out ("with c1 count",n1, "& c2 count",n2) 
            Call Out ("with c1 weight",Clust(c1)%smb, "& c2 weight",Clust(c2)%smb)
          End if
        
          Call Two_set_analysis (Clust(c1)%set,Clust(c2)%set, & 
                                 ni,Intr, Unyn=Unyn,Ind1=Ind1,Ind2=Ind2)
        
          nu= n1 + n2 - ni;  vl2= Clust(c2)%val

          DeAllocate(Clust(c2)%set,Clust(c2)%val)  
          Allocate(Clust(c2)%set(nu),Clust(c2)%val(nu))
          
          Clust(c2)%n= nu;  Clust(c2)%set= Unyn(:nu)
          
          Clust(c2)%val= 0;  Clust(c2)%val(Ind1)= Clust(c1)%val 
          Clust(c2)%val(Ind2)= Clust(c2)%val(Ind2) + vl2
          Clust(c2)%smb= Sum(Clust(c2)%val)

          DeAllocate(vl2, Intr,Unyn,Ind1,Ind2)

          nl= ncl;  Call Del_mrg (nl,c1,c2, Clust%smb,cls,ord, ncl)
          Clust(c1)%n= 0;  ind= ncl + 1

!         Test cluster sizes to see if the larger size of Clust(c2)
!         requires reordering

          siz(:ncl)= Clust(ord(:ncl))%smb;  siz(ind:)= 0
          Call Sort (.true.,siz(:ncl), key(:ncl),ReOrd)
          
          If (pr_out > 1) then
            Call Out ("Revised c2 set = union of c1 & c2 sets")
            Call Out ("with count",nu, "& weight",Clust(c2)%smb)
          End if

          If (ReOrd) then
            ord= 0
            Do cl= 1,n0
              If (cls(cl) > 0) then
                l= cls(cl);  k= key(l);  cls(cl)= k;  ord(k)= cl
              End if
            End do

            If (pr_out > 1) then
              Call Out ("Actual cluster sizes", Clust%smb)
              Call Out ("Reordered cluster sizes for reduced # clusters", siz(:ncl))
              Call Out ("Cluster size reordering key", key(:ncl))
              Call Out ("Revised initial to current cluster mapping", cls)
              Call Out ("Revised current cluster ordering", ord(:ncl))
            End if
          End if

!         Update Borda for the original clusters, deleting 'c1' and
!         using the already updated Clust(c2). Then recompute the 
!         corresponding correlations

          Call New_corr (c1,c2,nc,mr,Adj,rank_pt, Clust%smb, ballot,Clust(c2), Borda,Corr)

!         Restart the merge cycle after a merge 

          i= 1;  Cycle Least_clust  
        Else
          If (pr_out > 1) then
            Call Out ("Do not merge sets with correlation",cor, ln=1)
            Call Out ("for set c1",c1, "& set c2",c2)
          End if
        End if

        j= j + 1   ! Try to merge 'c1' to the next 'c2'
      End do Merge_loop

      If (Clust(c1)%smb <= parm(2)) then  ! Delete small cluster 'c1' after no merging
        nl= ncl;  Call Del_mrg (nl,c1,0, Clust%smb,cls,ord, ncl)
        Clust(c1)%n= 0;  ind= ncl + 1
      End if

      i= i1  ! Proceed to the merge cycle for the next smallest cluster 'c1'
    End do Least_clust

    ord= 0;  ind= ncl + 1
    Do cl= 1,n0
      If (Clust(cl)%smb > 0) then
        l= ind - cls(cl);  cls(cl)= l;  ord(l)= cl
      End if
    End do

    lst(:ncl)= "ID"
    If (Any(lst(:ncl) /= ord(:ncl))) then
      Call Copy_set_list (ord(:ncl), Clust, Bal(:ncl))
      Call Copy_set_list (lst(:ncl), Bal(:ncl), Clust(:ncl))

      Borda(:,:ncl)  = Borda(:,ord(:ncl))
      Corr(:ncl,:ncl)= Corr(ord(:ncl),ord(:ncl))
    End if

    Call DeAlloc_set_list_ar (Bal(:ncl))
    Call DeAlloc_set_list_ar (Clust(ind:n0))

    cnt(:ncl)= Clust(:ncl)%n;  siz(:ncl)= Clust(:ncl)%smb
    tot_mb= Sum(cnt(:ncl));  tot_wt= Sum(siz(:ncl))

    Borda(:,ind:)= 0;  Corr(:,ind:)= 0;  Corr(ind:,:)= 0
    cnt(ind:)= 0;  siz(ind:)= 0

    If (pr_out > 1) then
      Call Out ("Initial to final cluster mapping", cls)
      Call Out ("Final cluster ordering", ord(:ncl))

      Call Out ("Original cluster sizes",siz0)
      Call Out ("Final cluster sizes",siz(:ncl))
      Call Out ("with total weight",tot_wt)

      Call Out ("Original cluster membership counts",cnt0)
      Call Out ("Final cluster membership",cnt(:ncl))
      Call Out ("with total membership",tot_mb)

      Call Out (-1,"Final Borda Counts",Borda(:,:ncl))
      Call Out (1,"Final correlations",Corr(:ncl,:ncl))
    End if    
  End Subroutine Merge_corr

  Subroutine Merge_corr2 (nc,n0,ptval, parm,ballot,ballot2, Clust,Borda,Corr, ncl,cls, Adj)
  
!   Like Merge_matrix but uses the Set_list and 'ballot' data structures.
!   This version typically uses substantiallyh less memory than Merge_matrix 
!   to represent the clusters in terms of the original ballots rather
!   than slate ballots.
  
!   At all intermediate stages 'Clust', 'Borda', and 'Corr' stay the
!   same except when clusters c1 and c2 merge into c2 or clusters are
!   deleted. Full deletion and reordering are done only at the end via 'ord'.

    Integer,  Intent(in) :: nc            ! # candidates
    Integer,  Intent(in) :: n0            ! Initial # clusters
    Real,     Intent(in) :: ptval(:)      ! (mt) Point values of rating levels
    Real,     Intent(in) :: parm(:)       ! (3) Algorithm parameters: 
                                          !   (1) = cluster correlation merge criterion
                                          !   (2) = cluster deletion criterion
                                          !   (3) = minimum number of clusters upon output

    Integer,  Intent(in) :: ballot(0:,:)  ! (0:mr,nb) Ranked or rated candidates 
    Integer,  Intent(in) :: ballot2(0:,:) ! (0:mr,nb) Corresponding ranking or rating levels, increasing 

    Type(Set_list), Intent(inout) :: Clust(:)  ! (n0=>ncl) Membership for the clusters, ordered
                                               !           by decreasing cluster size, with small 
                                               !           clusters deleted and substantially  
                                               !           overlapping clusters merged
                                               !  %n  = # members
                                               !  %smb= Cluster weight = Sum(Clust%val)
                                               !  %set(n)= List of members (increasing on output)
                                               !  %val(n)= Corresponding weighted memberships

    Real,  Intent(inout) :: Borda(:,:)  ! (nc,n0=>ncl) Borda Count for each candidate cluster
    Real,  Intent(inout) :: Corr(:,:)   ! (n0,n0) Correlation matrix between pairs of clusters, 
                                        !   except = 0 on the diagonal

    Integer, Intent(out) :: ncl         ! # final clusters
    Integer, Intent(out) :: cls(:)      ! (n0) Maps original clusters to final clusters, or 0 if none

    Real, Optional, Intent(in) :: Adj   ! Noise adjustment to subtract from mean vectors
                                        !   before correlating
! Local:
    Real, Parameter :: eps= 0.00001
    Type(Set_list)  :: Bal(n0)
    Real,    Allocatable :: vl2(:)
    Integer, Allocatable :: Intr(:), Unyn(:), Ind1(:), Ind2(:)

    Real    :: siz(n0), siz0(n0)
    Integer :: cnt(n0), cnt0(n0), key(n0), inv(n0), lst(n0), ord(n0)
    
    Logical :: ReOrd
    Real    :: cor, tot_wt
    Integer :: i1, j1, k1, k2, mi, mu, n1, n2, ni, nl, nu, ind, min_num, tot_mb
    Integer :: i, j, k, l, m, n, c1, c2, cl

    min_num= Nint(parm(3))
    cls= "ID";  ord= cls

    cnt0= Clust%n;    tot_mb= Sum(cnt0)
    siz0= Clust%smb;  tot_wt= Sum(siz0)

    If (pr_out > 1) then
      Call Out ("Enter 'Merge_corr2' with # clusters",n0, &
                "with tot_wt weight",tot_wt, ln=1)
      Call Out ("with cluster sizes",siz0)
    End if

    Call List_of_true (siz0 <= parm(2), n,lst)

    If (n >= 1 .and. n0 <= min_num) then
      If (pr_out > 1) Call Out ("Too few clusters: Delete small clusters with no merging")

      nl= n0  
      Do i= 1,n
        c1= lst(i);  Call Del_mrg (nl,c1,0, Clust%smb,cls,ord, ncl);  Clust(c1)%n= 0
      End do

      ind= ncl + 1;  key(:ncl)= "ID"

      If (Any(ord(:ncl) /= key(:ncl))) then
        Call Copy_set_list (ord(:ncl), Clust, Bal(:ncl))
        Call Copy_set_list (key(:ncl), Bal(:ncl), Clust(:ncl))

        Borda(:,:ncl)  = Borda(:,ord(:ncl))
        Corr(:ncl,:ncl)= Corr(ord(:ncl),ord(:ncl))
      End if

      Call DeAlloc_set_list_ar (Bal(:ncl))
      Call DeAlloc_set_list_ar (Clust(ind:n0))

      Clust(ind:)%n= 0;  Clust(ind:)%smb= 0
      Borda(:,ind:)= 0;  Corr(:,ind:)= 0;  Corr(ind:,:)= 0

      If (pr_out > 1) then
        Call Out ("Reduced # clusters in 'Merge_corr2' from",n0, &
                  "to",ncl, ln=1)
        Call Out ("Initial to final cluster mapping", cls)
        Call Out ("Final cluster ordering", ord(:ncl))

        Call Out ("Final cluster sizes",Clust(:ncl)%smb)
        Call Out ("Final cluster membership",Clust(:ncl)%n)
      End if
      Return
    End if
    
!   Merge and delete from smallest to largest cluster

    ncl= n0;  siz= 0;  key= 0;  i= 1 
    ord= Reverse(cls);  cls= ord  ! Reverse the order

    Least_clust : Do while (i <= ncl .and. ncl > min_num)
      c1= ord(i);  i1= i + 1;  j= i1
      If (Clust(c1)%smb < eps) then
        i= i1;  Cycle Least_clust
      End if
      
      Merge_loop : Do while (j <= ncl)
        c2= ord(j)
        If (Clust(c2)%smb < eps) then
          j= j + 1;  Cycle Merge_loop
        End if

        cor= Corr(c1,c2);  n1= Clust(c1)%n;  n2= Clust(c2)%n

        If (cor > parm(1)) then !  Merge correlated cluster c1 with c2 by 
                                !  taking the union, including adding values
          mi= Min(n1,n2);  mu= n1 + n2
          Allocate(vl2(n2), Intr(mi),Unyn(mu), Ind1(n1), Ind2(n2))

          If (pr_out > 1) then
            Call Out ("Merge sets with correlation",cor, ln=1)
            Call Out ("for set c1",c1, "& set c2",c2)
            Call Out ("with c1 count",n1, "& c2 count",n2) 
            Call Out ("with c1 weight",Clust(c1)%smb, "& c2 weight",Clust(c2)%smb)
          End if
        
          Call Two_set_analysis (Clust(c1)%set,Clust(c2)%set, & 
                                 ni,Intr, Unyn=Unyn,Ind1=Ind1,Ind2=Ind2)
        
          nu= n1 + n2 - ni;  vl2= Clust(c2)%val

          DeAllocate(Clust(c2)%set,Clust(c2)%val)  
          Allocate(Clust(c2)%set(nu),Clust(c2)%val(nu))
          
          Clust(c2)%n= nu;  Clust(c2)%set= Unyn(:nu)
          
          Clust(c2)%val= 0;  Clust(c2)%val(Ind1)= Clust(c1)%val 
          Clust(c2)%val(Ind2)= Clust(c2)%val(Ind2) + vl2
          Clust(c2)%smb= Sum(Clust(c2)%val)

          DeAllocate(vl2, Intr,Unyn,Ind1,Ind2)

!         Change Borda and correlations, merging 'c1' with 'c2' and deleting 'c1'

          nl= ncl;  Call Del_mrg (nl,c1,c2, Clust%smb,cls,ord, ncl)
          Clust(c1)%n= 0;  ind= ncl + 1

!         Test cluster sizes to see if the larger size of Clust(c2)
!         requires reordering

          siz(:ncl)= Clust(ord(:ncl))%smb;   siz(ind:)= 0
          Call Sort (.true.,siz(:ncl), key(:ncl),ReOrd)
          
          If (pr_out > 1) then
            Call Out ("Revised c2 set = union of c1 & c2 sets")
            Call Out ("with count",nu, "& weight",Clust(c2)%smb)
          End if

          If (ReOrd) then
            ord= 0
            Do cl= 1,n0
              If (cls(cl) > 0) then
                l= cls(cl);  k= key(l);  cls(cl)= k;  ord(k)= cl
              End if
            End do

            If (pr_out > 1) then
              Call Out ("Cluster size reordering key", key(:ncl))
              Call Out ("Revised initial to current cluster mapping", cls)
              Call Out ("Revised current cluster ordering", ord(:ncl))
            End if
          End if

!         Update Borda for the original clusters, deleting 'c1' and
!         using the already updated Clust(c2). Then recompute the 
!         corresponding correlations

          If (Present(Adj)) then
            Call New_corr2 (c1,c2,nc,ptval, Clust%smb, ballot,ballot2,Clust(c2), Borda,Corr, Adj)
          Else
            Call New_corr2 (c1,c2,nc,ptval, Clust%smb, ballot,ballot2,Clust(c2), Borda,Corr)
          End if

!         Restart the merge cycle immediatley after a merge 

          i= 1;  Cycle Least_clust  
        Else
          If (pr_out > 1) then
            Call Out ("Do not merge sets with correlation",cor, ln=1)
            Call Out ("for set c1",c1, "& set c2",c2)
          End if
        End if

        j= j + 1   ! Try to merge 'c1' to the next 'c2'
      End do Merge_loop

      If (Clust(c1)%smb <= parm(2)) then  ! Delete small cluster 'c1' after no merging
        nl= ncl;  Call Del_mrg (nl,c1,0, Clust%smb,cls,ord, ncl)
        Clust(c1)%n= 0;  ind= ncl + 1
      End if

      i= i1  ! Proceed to the merge cycle for the next smallest cluster 'c1'
    End do Least_clust

    ord= 0;  ind= ncl + 1
    Do cl= 1,n0
      If (Clust(cl)%smb > 0) then
        l= ind - cls(cl);  cls(cl)= l;  ord(l)= cl
      End if
    End do

    lst(:ncl)= "ID"
    If (Any(lst(:ncl) /= ord(:ncl))) then
      Call Copy_set_list (ord(:ncl), Clust, Bal(:ncl))
      Call Copy_set_list (lst(:ncl), Bal(:ncl), Clust(:ncl))

      Borda(:,:ncl)  = Borda(:,ord(:ncl))
      Corr(:ncl,:ncl)= Corr(ord(:ncl),ord(:ncl))
    End if

    Call DeAlloc_set_list_ar (Bal(:ncl))
    Call DeAlloc_set_list_ar (Clust(ind:n0))

    cnt(:ncl)= Clust(:ncl)%n;  siz(:ncl)= Clust(:ncl)%smb
    tot_mb= Sum(cnt(:ncl));  tot_wt= Sum(siz(:ncl))

    Borda(:,ind:)= 0;  Corr(:,ind:)= 0;  Corr(ind:,:)= 0
    cnt(ind:)= 0;  siz(ind:)= 0

    If (pr_out > 1) then
      Call Out ("Initial to final cluster mapping", cls)
      Call Out ("Final cluster ordering", ord(:ncl))

      Call Out ("Original cluster sizes",siz0)
      Call Out ("Final cluster sizes",siz(:ncl))
      Call Out ("with total weight",tot_wt)

      Call Out ("Original cluster membership counts",cnt0)
      Call Out ("Final cluster membership",cnt(:ncl))
      Call Out ("with total membership",tot_mb)

      Call Out (-1,"Final Borda Counts",Borda(:,:ncl))
      Call Out (1,"Final correlations",Corr(:ncl,:ncl))
    End if    
  End Subroutine Merge_corr2

   Subroutine New_corr (c1,c2,nc,mr, Adj,rank_pt, siz, ballot,Clust2, Borda,Corr)

!    Recompute the Borda Count for updated original cluster 'c2',
!    Also its correlations with the other  clusters.
!    Delete the original cluster 'c1'.
  
     Integer,        Intent(in) :: c1           ! Original cluster to delete
     Integer,        Intent(in) :: c2           ! Original cluster c1 is merged to
     Integer,        Intent(in) :: nc           ! # candidates
     Integer,        Intent(in) :: mr           ! Max # candidates to be ranked or rated
     Real,           Intent(in) :: Adj          ! Noise adjustment to subtract from mean vectors
                                                !   before correlating
     Real,           Intent(in) :: rank_pt(:)   ! (mr) Point conversion of ranking levels
     Real,           Intent(in) :: siz(:)       ! (n0) Current cluster sizes

     Integer,        Intent(in) :: ballot(0:,:) ! (0:mr,nb) Ballot ranking structure for clustering
     Type(Set_list), Intent(in) :: Clust2       ! Merged cluster 'c2' membership
                                                !  %n     = # ballot members of the cluster
                                                !  %smb   = Cluster weight = size = Sum(val)
                                                !  %set(n)= List of ballot members
                                                !  %val(n)= Ballot weights = sizes

     Real,        Intent(inout) :: Borda(:,:)   ! (nc,n0) Borda Count for clusters, revised for 'c2'
     Real,        Intent(inout) :: Corr(:,:)    ! (n0,n0) Revised correlations with cluster 'c2' 
! Local:
     Real    :: bdl(nc), bd2(nc)
     Real    :: wb, wt, norm
     Integer :: b, k, n, cl, n0

     If (pr_out > 1) Call Out ("Enter 'New_corr'")

     n0= Size(siz);  Borda(:,c1)= 0;  Borda(:,c2)= 0
     Corr(:,c1)= 0;  Corr(c1,:)= 0
     Corr(:,c2)= 0;  Corr(c2,:)= 0

     wt= 0
     Do k= 1,Clust2%n
       b= Clust2%set(k);  wb= Clust2%val(k)
         n= ballot(0,b)
       Borda(ballot(1:n,b),c2)= Borda(ballot(1:n,b),c2) + rank_pt(:n) * wb
       wt= wt + wb
     End Do 

     Borda(:,c2)= Borda(:,c2) / wt
     bd2= Borda(:,c2) - Adj;  norm= Sqrt(Sum(bd2**2))
     bd2= bd2 / norm

     Do cl= 1,n0
       If (cl == c2 .or. siz(cl) <= 0) Cycle

       bdl= Borda(:,cl) - Adj;  norm= Sqrt(Sum(bdl**2))
       Corr(cl,c2)= Dot_product(bd2,bdl) / norm
       Corr(c2,cl)= Corr(cl,c2)
     End do

     If (pr_out > 1) then
       Call Out ("New cluster sizes",siz)
       Call Out ("New Borda Count for the merged set",Borda(:,c2))
       Call Out ("New correlations with this merged set",Corr(:,c2))
     End if
   End Subroutine New_corr

   Subroutine New_corr2 (c1,c2,nc,ptval, siz,ballot,ballot2, Clust2, Borda,Corr, Adj)

!    Recompute the Borda Count for updated original cluster 'c2',
!    Also its correlations with the other  clusters.
!    Delete the original cluster 'c1'.
  
     Integer,  Intent(in) :: c1            ! Cluster to delete
     Integer,  Intent(in) :: c2            ! Cluster c1 is merged to
     Integer,  Intent(in) :: nc            ! # candidates
     Real,     Intent(in) :: ptval(:)      ! (mt) Point values of rating levels
     Real,     Intent(in) :: siz(:)        ! (n0) Current cluster sizes

     Integer,  Intent(in) :: ballot(0:,:)  ! (0:mr,nb) Ranked or rated candidates 
     Integer,  Intent(in) :: ballot2(0:,:) ! (0:mr,nb) Corresponding ranking or rating levels, increasing 

     Type(Set_list), Intent(in) :: Clust2  ! Merged cluster 'c2' membership
                                           !  %n     = # ballot members of the cluster
                                           !  %smb   = Cluster weight = size = Sum(val)
                                           !  %set(n)= List of ballot members
                                           !  %val(n)= Ballot weights = sizes

     Real,  Intent(inout) :: Borda(:,:)    ! (nc,n0) Borda Count for clusters, revised for 'c2'
     Real,  Intent(inout) :: Corr(:,:)     ! (n0,n0) Revised correlations with cluster 'c2' 

     Real, Optional, Intent(in) :: Adj     ! Noise adjustment to subtract from mean vectors
                                           !   before correlating
! Local:
     Real    :: bdl(nc), bd2(nc)
     Real    :: wb, wt, norm
     Integer :: b, i, k, cl, mt, n0, n, nm

     If (pr_out > 1) Call Out ("Enter 'New_corr2'")

     mt= Size(ptval);  n0= Size(siz)
     Borda(:,c1)= 0;  Borda(:,c2)= 0;  Corr(:,c1)= -1;  Corr(c1,:)= -1
     Corr(:,c2)= -1;  Corr(c2,:)= -1

     wt= 0;  nm= Clust2%n
     Do i= 1,nm
       b= Clust2%set(i);  wb= Clust2%val(i);  n= ballot(0,b) 
       Borda(ballot(1:n,b),c2)= Borda(ballot(1:n,b),c2) + ptval(ballot2(1:n,b)) * wb
       wt= wt + wb
     End Do 

     Borda(:,c2)= Borda(:,c2) / Max(wt, 0.0001)

     bd2= Borda(:,c2);  If (Present(Adj)) bd2= bd2 - Adj
     norm= Sqrt(Sum(bd2**2));  bd2= bd2 / norm

     Do cl= 1,n0
       If (cl == c2 .or. siz(cl) <= 0) Cycle

       bdl= Borda(:,cl);  If (Present(Adj)) bdl= bdl - Adj
       norm= Sqrt(Sum(bdl**2));  bdl= bdl / norm

       Call Dot_product_M (Real(Dot_fac), bd2,bdl, i,Corr(c1,c2))
       Corr(c2,c1)= Corr(c1,c2)
     End do

     If (pr_out > 1) then
       Call Out ("New cluster sizes",siz)
       Call Out ("New Borda Count for the merged set",Borda(:,c2))
       Call Out ("New correlations with this merged set",Corr(:,c2))
     End if
   End Subroutine New_corr2


  Subroutine Merge_matrix (n0,parm, Memb, cnt,siz, ncl,cls)
  
!   Merge and delete using the matrix cluster structure 'Memb'

!   A simple algorithm for merging and deleting clusters 
!   defined by weighted memberships. Sub-optimmal, but good.
  
!   Let the clusters define a weighted graph, where the vertices 
!   represent clusters of elements drawn from a universe U. 
!   The membership Memb(e,c) of an element 'e' in a cluster 'c' 
!   is a non-negative real number, with zero denoting non-membership. 

!   The sum of these weighted memberships Memb(e,c) for a 
!   particular cluster 'c' defines the weight siz(c) the 
!   corresponding vertex. These weighted memberships also define 
!   the 'overlap' between any two clusters, a real number from 0 to 1 
!   which specifies the weight of the edge between the corresponding 
!   vertices, with zero denoting no edge.

!   The 'overlap',or intersection, between two clusters c1 and c2 
!   is computed by summing over edges 'e' the minimum of the membership 
!   weights Memb(e,c1) and Memb(e,c2), then dividing by the minimum of the 
!   cluster weights siz(c1) and siz(c2).

!   Now define the "merger" of vertices by taking the sum of 
!   membership weights:  Memb(e,c2)= Memb(e,c1) + Memb(e,c2) 
!   where c2 is the merger of c1 and c2.
  
!   The algorithm:
  
!   Order the vertices from smallest to largest by vertex weight, 
!   keeping this ordering throughout.  Starting with the vertex 
!   of least weight, see if it can be merged with the next 
!   heavier vertex, or the next, etc., according to the overlap 
!   parm(1). If no merger is possible, delete this vertex 
!   if it fails the min size parm(2)
  
!   Numerically, in the case of a merger, the first vertex is 
!   merged with the second vertex, removing the first vertex, 
!   followed by a reordering of the remaining vertices 
!   according to increasing size.
  
!   The algorithm is halted before reducing the number of 
!   vertices below the min number parm(3).

!   At all intermediate stages 'Memb', 'cnt', and 'siz' stay
!   the same except when clusters c1 and c2 merge into c2.
!   Deletion and reordering are done only at the end via 'ord'.
  
    Integer,     Intent(in) :: n0         ! Initial # clusters
    Real,        Intent(in) :: parm(:)    ! (3) Algorithm parameters: 
                                          !   (1) = cluster correlation merge criterion
                                          !   (2) = cluster deletion criterion
                                          !   (3) = minimum number of clusters upon output
    Real,     Intent(inout) :: Memb(:,:)  ! (ns,n0=>ncl) The membership weights for each 
                                          !   cluster 'cl' as a subset of (:,cl), from 
                                          !   initial to merged clusters, with merged 
                                          !   clusters ordered by decreasing weight = size
    Integer,  Intent(inout) :: cnt(:)     ! (n0=>ncl) Cluster membership counts
    Real,     Intent(inout) :: siz(:)     ! (n0=>ncl) Cluster weights = sizes,
                                          !   ordered by decreasing size and summing to 'np'
                                          !   (in and out)

    Integer,    Intent(out) :: ncl        ! # final clusters
    Integer,    Intent(out) :: cls(:)     ! (n0) Maps original clusters to final clusters
! Local:
    Real, Parameter :: eps= 0.00001
    Real    :: wt(n0), siz0(n0)
    Integer :: cnt0(n0), key(n0), inv(n0), ord(n0), lst(n0)
    
    Logical :: ReOrd
    Real    :: over, tot_wt
    Integer :: i, j, k, l, n, c1, c2, cl, i1, n1, nl, ind, min_num, tot_mb

    min_num= Nint(parm(3))
    cls= "ID";  ord= cls

    cnt0= cnt;  tot_mb= Sum(cnt0)
    siz0= siz;  tot_wt= Sum(siz0)

    If (pr_out > 1) then
      Call Out ("Enter 'Merge_matrix' with # clusters",n0, &
                "with total weight",tot_wt, ln=1)
      Call Out ("and total membership count",tot_mb)

      Call Out ("from cluster sizes",siz0)
      Call Out ("and counts",cnt0)
    End if

    Call List_of_true (siz0 <= parm(2), n,lst)

    If (n >= 1 .and. n0 <= min_num) then
      If (pr_out > 1) Call Out ("Too few clusters: Delete small clusters with no merging")

      ncl= n0  
      Do i= 1,n
        c1= lst(i);  nl= ncl;  Call Del_mrg (nl,c1,0, siz,cls,ord, ncl)
      End do

      ind= ncl + 1;  lst(:ncl)= "ID"

      If (Any(lst(:ncl) /= ord(:ncl))) then
        cnt(:ncl)= cnt0(ord(:ncl));  siz(:ncl)= siz0(ord(:ncl))
        Memb(:,:ncl)= Memb(:,ord(:ncl))
      End if

      cnt(ind:)= 0;  siz(ind:)= 0;  Memb(:,ind:)= 0

      If (pr_out > 1) then
        Call Out ("# clusters in 'Merge_matrix' reduced from",n0, "to",ncl, ln=1)
        Call Out ("Final cluster sizes",siz(:ncl))
        Call Out ("Final cluster memberships",cnt(:ncl))
      End if
      Return
    End if
    
!   Merge and delete from smallest to largest cluster
!   by reversing the order

    ncl= n0;  key= 0;  i= 1
    ord= Reverse(cls);  cls= ord

    Least_clust : Do while (i <= ncl .and. ncl > min_num)
      c1= ord(i);  i1= i + 1;  j= i1
      If (siz(c1) < eps) then
        i= i1;  Cycle Least_clust
      End if
      
      Merge_loop : Do while (j <= ncl)
        c2= ord(j)
        If (siz(c2) < eps) then
          j= j + 1;  Cycle Merge_loop
        End if

        over= Sum(Min(Memb(:,c1), Memb(:,c2))) / Min(siz(c1), siz(c2))
        
        If (over > parm(1)) then !  Merge overlapped cluster c1 with c2 by summing

          Memb(:,c2)= Memb(:,c1) + Memb(:,c2) 
          cnt(c2)= Count(Memb(:,c2) > 0.0001);  siz(c2)= Sum(Memb(:,c2))

!         Merge 'c1' with 'c2', deleting 'c1'

          nl= ncl;  Call Del_mrg (nl,c1,c2, siz,cls,ord, ncl)
          ind= ncl + 1

!         Test cluster sizes to see if the larger size of Clust(c2)
!         requires reordering

          wt(:ncl)= siz(ord(:ncl));  wt(ind:)= 0
          Call Sort (.true.,wt(:ncl), key(:ncl),ReOrd)
          
          If (pr_out > 1.5) then
            Call Out ("Merge sets c1 & c2 with overlap",over, ln=1)
            Call Out ("for c1",c1, "& c2",c2)
            Call Out ("Merged c2 = union of c1 & c2 with count",cnt(c2), &
                      "& size",siz(c2))
          End if

          If (ReOrd) then
            ord= 0
            Do cl= 1,n0
              If (cls(cl) > 0) then
                l= cls(cl);  k= key(l);  cls(cl)= k;  ord(k)= cl
              End if
            End do

            If (pr_out > 1.5) then
              Call Out ("Cluster size reordering key", key(:ncl))
              Call Out ("Revised initial to current cluster mapping", cls)
              Call Out ("Revised current cluster ordering", ord(:ncl))
            End if
          End if

!         Restart the merge cycle immediatley after a merge 

          i= 1;  Cycle Least_clust  
        Else
          If (pr_out > 1.5) then
            Call Out ("Do not merge sets c1 & c2 with overlap",over, ln=1)
            Call Out ("for c1",c1, "& c2",c2)
          End if
        End if

        j= j + 1   ! Try to merge 'c1' to the next 'c2'
      End do Merge_loop

      If (siz(c1) <= parm(2)) then  ! Delete small cluster 'c1' after no merging
        nl= ncl;  Call Del_mrg (nl,c1,0, siz,cls,ord, ncl);  ind= ncl + 1
      End if

      i= i1  ! Proceed to the merge cycle for the next smallest cluster 'c1'
    End do Least_clust

    ord= 0;  ind= ncl + 1
    Do cl= 1,n0
      If (siz(cl) > 0) then
        l= ind - cls(cl);  cls(cl)= l;  ord(l)= cl
      End if
    End do

    lst(:ncl)= "ID"
    If (Any(lst(:ncl) /= ord(:ncl))) then
      cnt(:ncl)= cnt(ord(:ncl));  siz(:ncl)= siz(ord(:ncl))
      Memb(:,:ncl)= Memb(:,ord(:ncl))
    End if

    tot_mb= Sum(cnt(:ncl));  cnt(ind:)= 0
    tot_wt= Sum(siz(:ncl));  siz(ind:)= 0
    Memb(:,ind:)= 0

    If (pr_out > 1) then
      Call Out ("Reduced # clusters in 'Merge_matrix' from",n0, "to",ncl, ln=1)
      Call Out ("Initial to final cluster mapping", cls)
      Call Out ("Final cluster ordering", ord(:ncl))

      Call Out ("Original cluster sizes",siz0)
      Call Out ("Final cluster sizes",siz)
      Call Out ("with total weight",tot_wt)

      Call Out ("Original cluster membership counts",cnt0)
      Call Out ("Final cluster memberships",cnt(:ncl))
      Call Out ("with total membership",tot_mb)
    End if
  End Subroutine Merge_matrix


  Subroutine Del_mrg (nl,c1,c2, siz,cls,ord, ncl)

!   Delete the original cluster c1. Also merge it with original cluster c2 if c2  > 0.

    Integer,    Intent(in) :: nl         ! Current # clusters, reduced by 1 here, with the deletion of c1
    Integer,    Intent(in) :: c1         ! Original cluster to delete
    Integer,    Intent(in) :: c2         ! Merge c1 with original cluster c2 if c2 > 0

    Real,    Intent(inout) :: siz(:)     ! (n0) The current sizes of the original clusters, with siz(c1)
                                         !      set to 0 here. If c2 > 0, assume that siz(c2) already 
                                         !      includes the merger with c1
    Integer, Intent(inout) :: cls(:)     ! (n0) Mapping from original clusters to current 
                                         !      clusters in order of increasing size
    Integer,   Intent(out) :: ord(:)     ! (n0) Inverse mapping from current clusters 
                                         !      to non-deleted original clusters
    Integer,   Intent(out) :: ncl        ! Reduced # clusters = nl - 1
! Local:
    Real, Parameter :: eps=0.000001
    Integer :: i, cls0(Size(cls))

    If (siz(c1) < eps) then
      Call Out ("Error in 'Del_mrg': cluster c1 already deleted",c1, "for c2",c2,  ln=1)
      Call Out ("Initial to current cluster mapping", cls)
      Call Out ("Its inverse mapping", ord(:ncl))
      Call Out ("Current cluster sizes of clusters", siz)
      Stop
    End if

    ncl= nl - 1;  cls0= cls

    Where (cls > cls(c1)) cls= cls - 1  ! Remapping after deletion

    If (c2 > 0) then
      cls(c1)= cls(c2)  ! Merge c1 into c2
    Else
      cls(c1)= 0
    End if

!   Inverse mapping

    siz(c1)= 0;  ord= 0
    Do i= 1,Size(cls)
      If (siz(i) > 0) ord(cls(i))= i
    End do

    If (pr_out > 1) then
      If (c2 > 0) then
        Call Out ("Delete and merge cluster",c1, "with cluster",c2)
      Else
        Call Out ("Just delete small cluster",c1)
      End if
      Call Out ("Prior initial to current cluster mapping", cls0)
      Call Out ("Updated initial to current cluster mapping", cls)
      Call Out ("Its inverse mapping", ord(:ncl))
      Call Out ("Updated clusters sizes", siz)
    End if

  End Subroutine Del_mrg
    

   Subroutine PortionsD (np,nc,cut, siz,mean_vec, Zrate0,por)

!    Determine the portions corresponding to the current cluster set.

!    All ratings "in the noise" are zeroed out, except for'
!    the independents, which represent the noise. The 'cut'
!    value is chosen to represent the noise influence.

!    The minimum rating for the independentsis is the positive 
!    value of 'eps' - to to guarantee that at least the indepdents
!    portion for each candidate is positive so that its portions 
!    sum to 1 over the clusters. This means that in certain cases
!    one of these candidates could be elected to retain
!    the proportionality of the other electeds.

!    The last step is to normalize the portions (cluster size * 
!    mean cluster rating vector) so that their sum
!    over all clusters, including independents, is 1.
   
     Integer,     Intent(in) :: np             ! # candidates to be elected (= 'tot_wt')
     Integer,     Intent(in) :: nc             ! # candidates
     Real,        Intent(in) :: cut            ! Cut level for zeroing, usually 'Noise_por'
                                               !   for ranking and 0 for rating

     Real(Dblp),  Intent(in) :: siz(:)         ! (ind) Cluster sizes
     Real(Dblp),  Intent(in) :: mean_vec(:,:)  ! (nc,ind) Mean vector of each cluster

     Real(Dblp), Intent(out) :: Zrate0(:)      ! (nc) Cluster weighted average of the
                                               !      zeroed ratings, including independents
     Real(Dblp), Intent(out) :: por(:,:,:)     ! (nc,ind,2)  Portion data
                                               ! 1 = zeroed ratings
                                               ! 2 = portion data 
!    Local:
     Real(Dblp), Parameter :: eps= 0.01, zero= 0
     Real(Dblp) :: tot_wt, sm(nc), por_avg(nc)
     Integer :: cl, cn, ind, ncl
    
     ind= Size(siz);  ncl= ind - 1

!    Noise zeroing for all regular clusters, but not for independents

     por(:,:,1)= mean_vec
     Where (Abs(mean_vec(:,:ncl)) <= cut) por(:,:ncl,1)= zero

!    Zero all negative values of portions for all clusters,
!    except keep an 'eps' for the independents

     Forall(cl=1:ncl) por(:,cl,2)= Max(siz(cl) * por(:,cl,1), zero)
     por(:,ind,2)= Max(siz(ind) * por(:,ind,1), eps)

     por_avg= Sum(por(:,:,2), 2)
     Forall(cl=1:ind) por(:,cl,2)= por(:,cl,2) / por_avg  ! Normalize portions

     tot_wt= Sum(siz);  sm= Sum(por(:,:,2), 2)  ! Check data: tot_wt = np, sm = 1

     Forall(cn=1:nc) Zrate0(cn)= Sum(siz * por(cn,:,1))
     Zrate0= Zrate0 / tot_wt
   
   End subroutine PortionsD

   Subroutine PortionsR (np,nc,cut, siz,mean_vec, Zrate0,por)

!    Determine the portions corresponding to the current cluster set.

!    All ratings "in the noise" are zeroed out, except for'
!    the independents, which represent the noise. The 'cut'
!    value is chosen to represent the noise influence.

!    The minimum rating for the independentsis is the positive 
!    value of 'eps' - to to guarantee that at least the indepdents
!    portion for each candidate is positive so that its portions 
!    sum to 1 over the clusters. This means that in certain cases
!    one of these candidates could be elected to retain
!    the proportionality of the other electeds.

!    The last step is to normalize the portions (cluster size * 
!    mean cluster rating vector) so that their sum
!    over all clusters, including independents, is 1.
   
     Integer, Intent(in) :: np             ! # candidates to be elected (= 'tot_wt')
     Integer, Intent(in) :: nc             ! # candidates
     Real,    Intent(in) :: cut            ! Cut level for zeroing, usually 'Noise_por'
                                           !   for ranking and 0 for rating

     Real,    Intent(in) :: siz(:)         ! (ind) Cluster sizes
     Real,    Intent(in) :: mean_vec(:,:)  ! (nc,ind) Mean vector of each cluster

     Real,   Intent(out) :: Zrate0(:)      ! (nc) Cluster weighted average of the
                                           !      zeroed ratings, including independents
     Real,   Intent(out) :: por(:,:,:)     ! (nc,ind,2)  Portion data
                                           ! 1 = zeroed ratings
                                           ! 2 = portion data 
!    Local:
     Real, Parameter :: eps= 0.01, zero= 0
     Real    :: tot_wt, sm(nc), por_avg(nc)
     Integer :: cl, cn, ind, ncl
    
     ind= Size(siz);  ncl= ind - 1

!    Noise zeroing for all regular clusters, but not for independents

     por(:,:,1)= mean_vec
     Where (Abs(mean_vec(:,:ncl)) <= cut) por(:,:ncl,1)= zero

!    Zero all negative values of portions for all clusters,
!    except keep an 'eps' for the independents

     Forall(cl=1:ncl) por(:,cl,2)= Max(siz(cl) * por(:,cl,1), zero)
     por(:,ind,2)= Max(siz(ind) * por(:,ind,1), eps)

     por_avg= Sum(por(:,:,2), 2)
     Forall(cl=1:ind) por(:,cl,2)= por(:,cl,2) / por_avg  ! Normalize portions

     tot_wt= Sum(siz);  sm= Sum(por(:,:,2), 2)  ! Check data: tot_wt = np, sm = 1

     Forall(cn=1:nc) Zrate0(cn)= Sum(siz * por(cn,:,1))
     Zrate0= Zrate0 / tot_wt
   
   End subroutine PortionsR
    

  Subroutine Density_n_boundary (ns,Memb,Clust, den_bdy)
  
!   Compute the density (average edge weight) within each cluster, plus the boundary
!   (average edge weight from within a cluster to its exterior).
  
    Integer,           Intent(in) :: ns         ! # slate clusters
    Type(Multi_listD), Intent(in) :: Memb(:)    ! (ns) Slate cluster membership data
                                                !   (0)%l = np
                                                !   (sl)%tx(0:nc) = normalized difference vector
    Type(Multi_listD), Intent(in) :: Clust(:)   ! (ind) Convergent cluster data
                                                !   %n  = # slate members
                                                !   %ls = list of slate members                                            
                                                !   %rx = unweighted slate memberships                                               
                                                !   %wt = weighted slate memberships

    Real(Dblp),    Intent(out) :: den_bdy(:,:)  ! (2,ind) (1,k) = Density for each cluster (should be strong)
                                                !         (2,k) = Boundary for each cluster (should be weak)
                                                !  density = edge weights between slate clusters within a cluster, 
                                                !    averaged over the product of the endpoint weighted cluster memberships
                                                !  boundary = edge weights from a slate clusters within a cluster 
                                                !    to slates ballot fully or partially outside it, averaged 
                                                !    over the product of the inside weighted cluster membership 
                                                !    and outside slate ballot weight
!   Local:
    Real(Dblp), Parameter :: eps= 0.001d0
    Logical    :: Inner(ns)
    Real(Dblp) :: fac, cor, wi, wt1, wt2, wt3, cor_in, cor_out, &
                  wt_in, wt_out
    Real    :: min_dens, max_bdy
    Integer :: i, j, k, l, m, n, s1, s2, sl, ind, ncl

    ind= Size(Clust);  ncl= ind - 1
    den_bdy(1,:)= 1;  den_bdy(2,:)= 0

    If (pr_out > 1) Call Out ("Enter 'Density_n_boundary' with # reg clusters", ncl, ln=1)

!   Compute the cluster densities = mean positive modified correlations
!   between the members of a cluster, including self correlations,
!   weighted by the memberships in the cluster

    Cluster_loop1 : Do k= 1,ind
      n= Clust(k)%n
      wt_in= Sum(Clust(k)%wt**2);  cor_in= wt_in

      Memb_1oop1 : Do i= 1,n-1
        s1= Clust(k)%ls(i);  wi= Clust(k)%wt(i)

        Memb_1oop2 : Do j= i+1,n  
          s2= Clust(k)%ls(j)
          Call Dot_product_M (Dot_fac, Memb(s1)%tx(1:),Memb(s2)%tx(1:), l,cor)

          wt1= wi * Clust(k)%wt(j);  wt_in= wt_in + wt1
          If (cor > 0) cor_in= cor_in + cor * wt1
        End do Memb_1oop2
      End do Memb_1oop1

      den_bdy(1,k)= cor_in / wt_in
    End do Cluster_loop1
      
    If (ind == 1) then
      If (pr_out > 1) Call Out ("For no reg clusters, Independents density", den_bdy(1,1))
      Return
    End if

!   Compute the cluster boundaries, correlating between members 
!   of a cluster and all others according to their non-membership
!   in the given cluster, weighted by the product of these
!   in and out memberships

    Cluster_loop2 : Do k= 1,ind
      n= Clust(k)%n;  wt_out= 0;  cor_out= 0
      Inner= .false.;  Inner(Clust(k)%ls)= .true.

      Inside_loop : Do i= 1,n
        s1= Clust(k)%ls(i);  wi= Clust(k)%wt(i)

        Outside_loop : Do s2= 1,ns
          wt2= wi * Memb(s2)%fsx

          If (Inner(s2)) then
                j= First_true(Clust(k)%ls == s2)
              fac= 1 - Clust(k)%rx(j)           ! Reduced outer membership
              If (fac < eps) Cycle Outside_loop
            wt2= wt2 * fac
          End if

          Call Dot_product_M (Dot_fac, Memb(s1)%tx(1:),Memb(s2)%tx(1:), l,cor)
          
          wt_out= wt_out + wt2  
          If (cor > 0) cor_out= cor_out + cor * wt2
        End do Outside_loop
      End do Inside_loop
      
      den_bdy(2,k)= cor_out / wt_out
    End do Cluster_loop2

    If (ncl > 0) then
      min_dens= Minval(den_bdy(1,:ncl))
      max_bdy = Maxval(den_bdy(2,:ncl))

      If (pr_out > 1) then
        Call Out ("For # regular clusters", ncl, ln=1)
        Call Out ("Minimum cluster density", min_dens, &
                  "Maximum cluster boundary", max_bdy)
        Call Out (1,"Cluster density & boundary data", den_bdy)
      End if      
    End if
    
  End Subroutine Density_n_boundary

  Subroutine Dens_n_bound0 (ns, Memb,Clust_memb, den_bdy)
  
!   Compute the density (average edge weight) within each cluster, plus the boundary
!   (average edge weight from within a cluster to its exterior).
  
    Integer,           Intent(in) :: ns       ! # slate clusters
    Type(Multi_listD), Intent(in) :: Memb(:)  ! (ns) Slate cluster membership data
                                              !   (0)%l = np
                                              !   (sl)%tx(0:nc) = normalized difference vector
    Real,    Intent(in) :: Clust_memb(:,:)  ! (ns,ind) Membership weight of each slate cluster
                                            !   in each cluster, with clusters in decreasing
                                            !   order of size, with sizes summing to 'np'

    Real,   Intent(out) :: den_bdy(:,:)     ! (2,ind) (1,k) = Density for each cluster (should be strong)
                                            !         (2,k) = Boundary for each cluster (should be weak)
                                            !  density = edge weights between slate clusters within a cluster, 
                                            !    averaged over the product of the endpoint weighted cluster memberships
                                            !  boundary = edge weights from a slate clusters within a cluster 
                                            !    to slates ballot fully or partially outside it, averaged 
                                            !    over the product of the inside weighted cluster membership 
                                            !    and outside slate ballot weight
!   Local:
    Real, Parameter :: eps= 0.001d0
    Logical :: Inner(ns)
    Integer :: lst(ns)

    Real(Dblp) :: cor
    Real    :: min_dens, max_bdy
    Real    :: fac, wi, wt1, wt2, cor_in, cor_out, wt_in, wt_out
    Integer :: i, j, k, l, n, s1, s2, sl, ind, ncl

    ind= Size(Clust_memb,2);  ncl= ind - 1
    den_bdy(1,:)= 1;  den_bdy(2,:)= 0

!   Compute the cluster densities = mean positive modified correlations
!   between the members of a cluster, including self correlations,
!   weighted by the memberships in the cluster

    Cluster_loop1 : Do k= 1,ind
      Call List_of_true (Clust_memb(:,k) > 0, n,lst)

      wt_in= Sum(Clust_memb(lst(:n),k));  cor_in= wt_in

      Memb_1oop1 : Do i= 1,n-1
        s1= lst(i);  wi= Clust_memb(s1,k)

        Memb_1oop2 : Do j= i+1,n
          s2= lst(j)
          Call Dot_product_M (Dot_fac, Memb(s1)%tx(1:),Memb(s2)%tx(1:), l,cor)

          wt1= wi * Clust_memb(s2,k);  wt_in= wt_in + wt1
          If (cor > 0) cor_in= cor_in + cor * wt1
        End do Memb_1oop2
      End do Memb_1oop1

      den_bdy(1,k)= cor_in / wt_in
    End do Cluster_loop1
      
    If (ind == 1) then
      If (pr_out > 1) Call Out ("For no reg clusters, Independents density", den_bdy(1,1))
      Return
    End if

!   Compute the cluster boundaries, correlating between members 
!   of a cluster and all others according to their non-membership
!   in the given cluster, weighted by the product of these
!   in and out memberships

    Cluster_loop2 : Do k= 1,ind
      Call List_of_true (Clust_memb(:,k) > 0, n,lst)

      wt_out= 0;  cor_out= 0
      Inner= .false.;  Inner(lst(:n))= .true.

      Inside_loop : Do i= 1,n
        s1= lst(i);  wi= Clust_memb(s1,k)

        Outside_loop : Do s2= 1,ns
          wt2= wi * Memb(s2)%fsx

          If (Inner(s2)) then
            fac=  1 - Clust_memb(s2,k) / Memb(s2)%fsx ! Reduced outer membership
              If (fac < eps) Cycle Outside_loop
            wt2= wt2 * fac
          End if

          Call Dot_product_M (Dot_fac, Memb(s1)%tx(1:),Memb(s2)%tx(1:), l,cor)
          
          wt_out= wt_out + wt2  
          If (cor > 0) cor_out= cor_out + cor * wt2
        End do Outside_loop
      End do Inside_loop
      
      den_bdy(2,k)= cor_out / wt_out
    End do Cluster_loop2

    If (ncl > 0) then
      min_dens= Minval(den_bdy(1,:ncl))
      max_bdy = Maxval(den_bdy(2,:ncl))

      If (pr_out > 1) then
        Call Out ("For # regular clusters", ncl, ln=1)
        Call Out ("Minimum cluster density", min_dens, &
                  "Maximum cluster boundary", max_bdy)
        Call Out (1,"Cluster density & boundary data", den_bdy)
      End if      
    End if
    
  End Subroutine Dens_n_bound0

  Subroutine Objective_valuesR (np,parm, siz,dat, den_bdy, &
                               db_pen, pen,obj)
                               
!   Compute the "clustering objective" value for a cluster set 
!   after its convergence by 'Iterate' from 'Form_clusters'.

!   This value 'obj(1)' is the product of the "coherent vote" 
!   'obj(2)' times the penalty 'obj(3)', which in turn  
!   is the product of the individual penalty factors 'pen'.

!   The coherent vote is the weight or size of cluster 
!   times its density and  boundary penalties, summed over
!   the regular clusters. Here the density and boundary data 
!   'den_bdy' is computed by 'Density_n_boundary' from 'Iterate'.

!   The other penalty factors 'pen' apply to the cluster set 
!   as a whole:: minimum cluster size, maximum correlation 
!   between clusters, normalized size of the set of independents, 
!   and the convergence residual (quality and accuracy)

    Integer, Intent(in) :: np            ! # candidates to be elected = sum of cluster weights
    Real,    Intent(in) :: parm(:)       ! (12) Parameters for the penalty factors 
                                         !      of the objective function
                                         !  1:2 = cluster density soft cutoff limits (cos_rise)
                                         !  3:4 = cluster boundary soft cutoff limits (cos_fall)
                                         !  5:6 = cluster size soft cutoff limits (cos_rise)
                                         !  7:8 = cluster correlation soft cutoff limits (cos_fall)
                                         !  9:10 = independents size cutoff limits (cos_fall)
                                         !  11:12 = convergence residual limits (cos_fall)

    Real,    Intent(in) :: siz(:)        ! (ncl) Cluster sizes, excluding independents

    Real,    Intent(in) :: dat(:)        ! (4) Data for cluster set penalty factors 'pen'
                                         ! (1) = minimum size of regular clusters
                                         ! (2) = maximum correlation of regular clusters
                                         ! (3) = independents size, as a fraction of 1.0
                                         ! (4) = final convergence residual

    Real,    Intent(in) :: den_bdy(:,:)  ! (2,ind) Cluster density (1) and boundary (2)

    Real,   Intent(out) :: db_pen(:,:)   ! (2,ind) Cluster density (1) and boundary (2) 
                                         !     penalty factors (1 = no penalty & 0 = full penalty)
    
    Real,   Intent(out) :: pen(:)        ! (4) Cluster set 0 to 1 penalty factors
                                         ! (1) = min reg cluster size penalty 
                                         ! (2) = max reg cluster correlation penalty
                                         ! (3) = independents size penalty
                                         ! (4) = residual penalty

    Real,   Intent(out) :: obj(:)        ! (3) Objective values
                                         ! (1): overall cluster set objective =
                                         ! (2): coherent vote = sizes * dens pen * boun pen
                                         !      (summed over the regular clusters), times
                                         ! (3): overall cluster set penalty = Product(pen)
! Local:
    Real, Parameter :: Scale= 10
    Real    :: vote, min_size, max_corr, siz_ind, final_res
    Integer :: cl, ncl

!   Coherent Vote, using regular cluster sizes times density 
!   and boundary penalties, then scaled

    Do cl= 1,Size(db_pen,2)
      db_pen(1,cl)= Cos_rise(parm(1),parm(2), Real(den_bdy(1,cl)))
      db_pen(2,cl)= Cos_fall(parm(3),parm(4), Real(den_bdy(2,cl)))
    End do

    db_pen= Max(db_pen, obj_eps) ! Don't allow 0 penalty factors

    ncl= Size(siz)
    vote= Sum(siz * db_pen(1,:ncl) * db_pen(2,:ncl)) * (Scale / np)  

!   Compute cluster set penalties

    min_size= dat(1);  max_corr= dat(2)
    siz_ind = dat(3);  final_res= dat(4)

    pen(1)= Cos_rise(parm(5),parm(6),   min_size)   ! Min cluster size penalty
    pen(2)= Cos_fall(parm(7),parm(8),   max_corr)   ! Max correlation penalty
    pen(3)= Cos_fall(parm(9),parm(10),  siz_ind)    ! Large independents size penalty
    pen(4)= Cos_fall(parm(11),parm(12), final_res)  ! Large convergence residual penalty

    pen= Max(pen, obj_eps)  ! Don't allow a 0 penalty factor
    
!   Compute the cluster set objective values, avoiding 0 values,
!   leaving at least 1 usable cluster set

    obj(3)= Max(Product(pen), obj_eps2)
    obj(2)= vote
    obj(1)= Max(obj(2) * obj(3), obj_eps3)
    
  End Subroutine Objective_valuesR

  
  Subroutine Objective_valuesD (np,parm, siz,dat, den_bdy, &
                                db_pen, pen,obj)
                               
!   Compute the "clustering objective" value for a cluster set 
!   after its convergence by 'Iterate' from 'Form_clusters'.

!   This value 'obj(1)' is the product of the "coherent vote" 
!   'obj(2)' times the penalty 'obj(3)', which in turn  
!   is the product of the individual penalty factors 'pen'.

!   The coherent vote is the weight or size of cluster 
!   times its density and  boundary penalties, summed over
!   the regular clusters. Here the density and boundary data 
!   'den_bdy' is computed by 'Density_n_boundary' from 'Iterate'.

!   The other penalty factors 'pen' apply to the cluster set 
!   as a whole:: minimum cluster size, maximum correlation 
!   between clusters, normalized size of the set of independents, 
!   and the convergence residual (quality and accuracy)

    Integer,       Intent(in) :: np            ! # candidates to be elected = sum of cluster weights
    Real,          Intent(in) :: parm(:)       ! (12) Parameters for the penalty factors 
                                               !      of the objective function
                                               !  1:2 = cluster density soft cutoff limits (cos_rise)
                                               !  3:4 = cluster boundary soft cutoff limits (cos_fall)
                                               !  5:6 = cluster size soft cutoff limits (cos_rise)
                                               !  7:8 = cluster correlation soft cutoff limits (cos_fall)
                                               !  9:10 = independents size cutoff limits (cos_fall)
                                               !  11:12 = convergence residual limits (cos_fall)

    Real(Dblp),    Intent(in) :: siz(:)        ! (ncl) Cluster sizes, excluding independents

    Real(Dblp),    Intent(in) :: dat(:)        ! (4) Data for cluster set penalty factors 'pen'
                                               ! (1) = minimum size of regular clusters
                                               ! (2) = maximum correlation of regular clusters
                                               ! (3) = independents size, as a fraction of 1.0
                                               ! (4) = final convergence residual

    Real(Dblp),    Intent(in) :: den_bdy(:,:)  ! (2,ind) Cluster density (1) and boundary (2)

    Real(Dblp),   Intent(out) :: db_pen(:,:)   ! (2,ind) Cluster density (1) and boundary (2) 
                                               !     penalty factors (1 = no penalty & 0 = full penalty)
    
    Real(Dblp),   Intent(out) :: pen(:)        ! (4) Cluster set 0 to 1 penalty factors
                                               ! (1) = min reg cluster size penalty 
                                               ! (2) = max reg cluster correlation penalty
                                               ! (3) = independents size penalty
                                               ! (4) = residual penalty

    Real(Dblp),   Intent(out) :: obj(:)        ! (3) Objective values
                                               ! (1): overall cluster set objective =
                                               ! (2): coherent vote = sizes * dens pen * boun pen
                                               !      (summed over the regular clusters), times
                                               ! (3): overall cluster set penalty = Product(pen)
! Local:
    Real, Parameter :: Scale= 10
    Real    :: vote, min_size, max_corr, siz_ind, final_res
    Integer :: cl, ncl

!   Coherent Vote, using regular cluster sizes times density 
!   and boundary penalties, then scaled

    Do cl= 1,Size(db_pen,2)
      db_pen(1,cl)= Cos_rise(parm(1),parm(2), Real(den_bdy(1,cl)))
      db_pen(2,cl)= Cos_fall(parm(3),parm(4), Real(den_bdy(2,cl)))
    End do

    db_pen= Max(db_pen, obj_eps) ! Don't allow 0 penalty factors

    ncl= Size(siz)
    vote= Sum(siz * db_pen(1,:ncl) * db_pen(2,:ncl)) * (Scale / np)  

!   Compute cluster set penalties

    min_size= dat(1);  max_corr= dat(2)
    siz_ind = dat(3);  final_res= dat(4)

    pen(1)= Cos_rise(parm(5),parm(6),   min_size)   ! Min cluster size penalty
    pen(2)= Cos_fall(parm(7),parm(8),   max_corr)   ! Max correlation penalty
    pen(3)= Cos_fall(parm(9),parm(10),  siz_ind)    ! Large independents size penalty
    pen(4)= Cos_fall(parm(11),parm(12), final_res)  ! Large convergence residual penalty

    pen= Max(pen, obj_eps)  ! Don't allow a 0 penalty factor
    
!   Compute the cluster set objective values, avoiding 0 values,
!   leaving at least 1 usable cluster set

    obj(3)= Max(Product(pen), obj_eps2)
    obj(2)= vote
    obj(1)= Max(obj(2) * obj(3), obj_eps3)
    
  End Subroutine Objective_valuesD

  Subroutine Objective_values0 (np,parm, siz,dat, pen,obj)
                               
!   Compute the simplified overall objective and penalties

    Integer, Intent(in) :: np         ! # elected candidates =  
    Real,    Intent(in) :: parm(:)    ! (6) Parameters for the penalty factors 
                                      !  1:2 = cluster size soft cutoff limits (cos_rise)
                                      !  3:4 = cluster correlation soft cutoff limits (cos_fall)
                                      !  5:6 = independents size cutoff limits (cos_fall)

    Real,    Intent(in) :: siz(:)     ! (ncl)  Cluster sizes

    Real,    Intent(in) :: dat(:)     ! (3) Cluster set penalty data   
                                      ! (1) = minimum size of regular clusters
                                      ! (2) = maximum correlation of regular clusters
                                      ! (3) = independents size, as a fraction of 'np'

    Real,   Intent(out) :: pen(:)     ! (3) Cluster set penalty factors, with
                                      !     1 = no penalty & 0 = full penalty
                                      ! (1) = min cluster size penalty 
                                      ! (2) = max cluster correlation penalty
                                      ! (3) = independents size penalty

    Real,   Intent(out) :: obj(:)     ! (3) Objective values
                                      ! (1): overall cluster set objective = obj(2) * obj(3)
                                      ! (2): cluster vote= Sum of sizes of regular clusters
                                      ! (3): overall cluster set penalty = Product(pen)
! Local:
    Real, Parameter :: Scale= 10
    Real :: vote

    vote= Sum(siz) * (Scale/np)

    pen(1)= Cos_rise(parm(1),parm(2),  dat(1))  ! Min size penalty
    pen(2)= Cos_fall(parm(3),parm(4),  dat(2))  ! Max correlation penalty
    pen(3)= Cos_fall(parm(5),parm(6),  dat(3))  ! Independents size penalty

!   Compute the cluster set objective values, avoiding 0 values

    pen= Max(pen, obj_eps)
    obj(3)= Max(Product(pen), obj_eps2)
    obj(2)= vote
    obj(1)= Max(obj(2) * obj(3), obj_eps3)
    
  End Subroutine Objective_values0

  
   Subroutine Electeds_objective (Elected, siz,avg_rate,portion, obj,dat)

!    Compute the fitness of a possible set of elected candidates. This is the value of an
!    objective function to be maximized = the cluster averaged points 'avg_rate' for each 
!    candidate, summed over the elected candidates. This objective is then subjected to 
!    penalty scale factors derived from the deviation between the actual vote for a cluster  
!    and its true size. Here the actual vote for a cluster is the sum of the 'portions'   
!    of each elected candidate for that cluster.

!    The max objective is recorded in obj(2) with product of the deviation penalties in obj(3), 
!    with a value < 1 indicating the strength of the penalty. Thus the full objective is 
!    obj(1) = obj(2)*obj(3). 

!    Note that if there is only one cluster, then the the mean point vector constitutes a
!    version of the Borda Count. At the opposite extreme, if each elected candidate is represented
!    by only one cluster, then the result is just the weighted average point value of the 
!    elected candidates, using each candidate's associated cluster size as the weight. Then 
!    each candidate is also the Borda elected candidate for its associated cluster. Thus this 
!    clustering algorithm may be viewed as a version of Borda modified to adapt to penalties for 
!    lack of proportionality, using a combination of ordinary Borda within clusters and 
!    cluster average Borda between clusters.
   
!    In addition, note that the average ranking or rating 'avg_rate' is computed from adjusted values  
!    which may be negative. This gives a preference to candidates with high rankings or ratings from
!    some clusters but without very low or negative rankings or ratings from others = less polarizing 
!    candidates.

!    A cosine penalty is used for the true vs actual cluster sizes, so that small deviations
!    have no penalty, large deviations have the full penalty, transitioned smoothly by a cosine 
!    function in between.
   
!    The product obj(3) of the penalty factors over the clusters then scales the total vote obj(2) 
!    to get the fitness objective obj(1).

     Integer, Intent(in) :: Elected(:)    ! (np)   Possible set of elected cand's, 
                                          !          in any order
     Real,    Intent(in) :: siz(:)        ! (ind)  True cluster sizes = 
                                          !          target fractional # elected
     Real,    Intent(in) :: avg_rate(:)   ! (nc)   Cluster size average for the 
                                          !          zeroed ratings
     Real,    Intent(in) :: portion(:,:)  ! (nc,ind) Candidate portions 
                                          !            by cluster

     Real,   Intent(out) :: obj(:)   ! (3) 
                                     ! 1 = Full objective = obj(2) * obj(3)
                                     ! 2 = max objective = cluster averaged ratings summed over the elected candidates
                                     ! 3 = deviation penalty = product of  
                                     !     reg cluster size deviation penalties 
                                     !     dat(:ncl,2)
     Real,   Intent(out) :: dat(:,:) ! (ind,2)  
                                     ! 1 = Cluster size deviations = 
                                     !     actual - true,  for the elected 
                                     !     candidates (by summing portions)
                                     ! 2 = Cluster size deviation penalties
!  Local:
     Real, Allocatable :: dev(:), parm(:,:)
     Real     :: ex
     Integer  :: cl, np, ncl, ind

     ind= Size(portion,2);  ncl= ind - 1;  np= Size(Elected)
     
!    Default objectives
     
     obj(2)  = Max(Sum(avg_rate(Elected)) / np, obj_eps)  ! Objective to be maximized
     dat(:,1)= Sum(portion(Elected,:),1) - siz ! Cluster size deviations = actual - true
     dat(ind,2)= 1                             ! No penalty for independents
     
     If (ind <= 1) then
       obj(3)= 1
     Else
       Allocate (dev(ncl), parm(ncl,2))
       dev= Abs(dat(:ncl,1))  ! Absolute cluster sizes deviations
       
       Forall(cl=1:ncl) dat(cl,2)= Cos_fall(Parm_el(1),Parm_el(2), dev(cl))

       dat(:ncl,2)= Max(dat(:ncl,2), obj_eps)      ! Limit large deviation penalties
       obj(3)= Max(Product(dat(:ncl,2)), obj_eps2) ! Limit their product
     End if

     obj(1)= obj(2) * obj(3)
   End Subroutine Electeds_objective

    
  Subroutine Correlate_clusters1 (Clust, mx_cor,Corr)

!   Compute a correlation matrix between the cluster mean rating vectors,
!   using noise adjusted values, so that correlations may be negative.
  
    Type(Multi_listD), Intent(in) :: Clust(:)  ! (ncl) Cluster data
    
    Real, Intent(out) :: mx_cor(:)   ! (ncl) Maximum correlation between distinct clusters
    Real, Intent(out) :: Corr(:,:)   ! (ncl,ncl) Correlation matrix between the clusters
! Local:
    Real    :: cor
    Integer :: n, c1, c2, ncl

    ncl= Size(Clust);  mx_cor= 0;  Corr= 0  
    If (ncl <= 1) Return

    Do c1= 1,ncl-1
      Do c2= c1+1,ncl
        Call Dot_product_M (Real(Dot_fac), Real(Clust(c1)%tx(1:)), Real(Clust(c2)%tx(1:)), n,cor)
        Corr(c1,c2)= cor;  Corr(c2,c1)= cor
      End do
      mx_cor(c1)= Maxval(Corr(c1,:))
    End do

  End Subroutine Correlate_clusters1
  
  Subroutine Correlate_clusters (nc,ncl,Adj,rate, mx_cor,Corr)

!   Compute a correlation matrix between the cluster mean rating vectors,
!   using noise adjusted values, so that correlations may be negative.
  
    Integer, Intent(in) :: nc         ! # candidates
    Integer, Intent(in) :: ncl        ! # clusters
    Real,    Intent(in) :: Adj        ! Noise adjustment to subtract from mean vectors
                                      !   before correlating
    Real,    Intent(in) :: rate(:,:)  ! (nc,ncl) Mean rating vector for each cluster
    Real,   Intent(out) :: mx_cor(:)  ! (ncl) Maximum correlation between distinct clusters
    Real,   Intent(out) :: Corr(:,:)  ! (ncl,ncl) Correlation matrix between the clusters
! Local:
    Real    :: cor, norm, tx(nc,ncl)
    Integer :: n, c1, c2

    mx_cor= 0;  Corr= 0;  If (ncl <= 1) Return

    Do c1= 1,ncl
      tx(:,c1)= rate(:,c1) - Adj
      norm= Sqrt(Sum(tx(:,c1)**2))
      tx(:,c1)= tx(:,c1) / norm
    End do

    Do c1= 1,ncl-1
      Do c2= c1+1,ncl
        Call Dot_product_M (Real(Dot_fac), tx(:,c1), tx(:,c2), n,cor)
        Corr(c1,c2)= cor;  Corr(c2,c1)= cor
      End do
      mx_cor(c1)= Maxval(Corr(c1,:))
    End do

  End Subroutine Correlate_clusters
  
  
  Subroutine Converge_clusters (iq,icyc, nc,ng, mx,RS, tol, Memb, &
                                ncl,Clust, cnvI,cnvR)
  
!   Iteratively centroid the cluster set 'Clust' - the key step 
!   in converging the mean rating vectors.

!   The only input from 'Clust' that is used are the cluster mean  
!   rating vectors Clust(cl)%sx, plus the slate ballot data 'Memb'.
  
    Integer,             Intent(in) :: iq          ! initial cluster set
    Integer,             Intent(in) :: icyc        ! convergence cycle
    Integer,             Intent(in) :: nc          ! # candidates
    Integer,             Intent(in) :: ng          ! Initial # regular clusters

    Integer,             Intent(in) :: mx          ! max # centroiding iterations for this restart level
    Integer,             Intent(in) :: RS          ! 1 : Level 1 restart (strongest merge and delete requirements) - initially or 
                                                   !     after a merge
                                                   ! 2 : Level 2 restart (merge or delete more easily) - after no merges but 
                                                   !      possibly deletions
                                                   ! 3 : No restart (after no reduction in # clusters at level 2 unless to ncl = 1). 
                                                   !     Do final convergence
    Real,                Intent(in) :: tol(:)      ! (3) Convergence tolerances for Restart values 1..3

    Type(Multi_listD),   Intent(in) :: Memb(0:)    ! (0:ns) Slate cluster membership data. Required input: fsx,sx, tx
                                                   !   'fsx' = slate membership
                                                   !   'sx'(nc) = mean rating vector
                                                   !   'tx'(nc) = normalized difference vector
    Integer,            Intent(out) :: ncl         ! Final # regular clusters
    Type(Multi_listD),Intent(inout) :: Clust(0:)   ! (0:ind) Convergent cluster data
                                                   ! For 0:
                                                   !   %k  = Error code. Invalid results if < 0.
                                                   !   %n  = # regular clusters = ng in, ncl out ***
                                                   !   %l  = # clusters, including independents
                                                   !   %lt(ncl) : original numbering of clusters if # clusters is reduced 

                                                   !   %rx(3)  : Cluster set objective value = coherent vote * penalty
                                                   !   %wt(ns): Sum of slate cluster memberships over all reg
                                                   !            clusters, reduced by the fuzzy set factor,
                                                   !            to be <= 1.0 

                                                   !   %L1(-2:nc,0:ind) Candidate ordering data for adjusted rankings/ratings
                                                   !     (-2:0,:) = # top, or significant, or viable candidates
                                                   !     (1:,:)   = Decreasing ordering for %T0(:,:,1:2)

                                                   !   %T0(nc,0:ind,3) Cluster point data
                                                   !     (:,1:,1) = Cluster mean vectors
                                                   !     (:,1:,2) = Cluster mean vectors with zeroing
                                                   !     (:,1:,3) = Cluster portions for each candidate
                                                   !     (:,0,1:2)= Cluster size average of (:,1:,1:2)

                                                   ! For 1:ind:
                                                   !   %k  = # fuzzy reduced slate memberships
                                                   !   %n  = # slate ballots which are partial or full members
                                                   !   %o       : 1 for regular clusters, 2 for independents
                                                   !   %fsx     : slate ballot weight averaged fuzzy reduction factor
                                                   !   %fux     : cluster width
                                                   !   %sum_wt = sum of slate cluster memberships (in and out)***

                                                   !   %ls(n) = list of slate members
                                                   !   %wt(n) = new slate memberships
                                                   !   %px(nc) = variance of %sx
                                                   !   %sx(nc) = mean rating vector (in and out)***
                                                   !   %tx(0:nc) = centered and normalized difference vector, with 0 = norm
    Integer,            Intent(out) :: cnvI(0:)    ! (0:4)  Integer convergence data
                                                   !   0 : Restart type
                                                   !   1 : Success code : 2 = fully converged,   1 = reducing functional,
                                                   !                      0 = damped reduction, -1 = convergence failure
                                                   !   2 : final # regular clusters
                                                   !   3 : # centroiding iterations
                                                   !   4 : total # line search points
    Real,               Intent(out) :: cnvR(0:)    ! (0:2) Line search data = Rel(:,itr)
                                                   !   0 = Final line search parameter
                                                   !   1 = Final functional = update change
                                                   !   2 = Line search |output - input| vectors
! Local:
    Real,  Parameter :: eps= 0.0001
    Real,  Parameter :: dmpRS(3)= (/0.50d0, 0.30d0, 0.15d0/) ! Initial damping parameter
    Integer, Parameter :: Func_opt= 2
    
    Real(Dblp) :: rates(nc,ng,0:mx)  ! Record cluster mean vectors for each iteration
    Real    :: chng(0:mx)   ! Update change (functional value) after each iteration
    Real    :: Rel(0:2,mx)  ! Iteration line search updates
                            !   0 = Final line search damping parameter
                            !   1 = Final functional = update change
                            !   2 = Net line search normalized (output - input) vectors
    Integer :: nfcn(mx)     ! # function calls per iteration
    
    Real(Dblp) :: dg
    Real    :: x, y, dmp, dmp0, dmp1, rtf(0:nc,ng+1)
    Integer :: ier, ind, itc, itr, Success, key(ng)
    Integer :: i, k, n, cl, n0, n1, n2, n3, n4, nl, np, nu
    
    If (pr_out > 1) Call Out ("Enter 'Converge_clusters' from initial set",iq, &
                               "& cycle",icyc, ln=1)

    np= Memb(0)%l
    cnvI= 0;  cnvR= 0;  Rel= 0;  nfcn= 0
    nfcn= 0;  dmp= -1;  dmp1= dmpRS(RS)

    k= 0;  Success= 1;  ier= 1;  Clust(0)%k= ier
    n0= ng;  ind= n0 + 1;  key= "ID"

    rates= 0;  chng(0)= 10
    Forall(cl=1:n0) rates(:,cl,0)= Clust(cl)%sx
    
    Iteration_loop : Do itc= 1,mx
        
      Call Update_vec (itc, Nint(pr_out),dmp, nc,n0,             &
                       rates(:,:n0,itc-1), ncl,rates(:,:n0,itc), &
                       Rel(:,itc), nfcn(itc), ier)
        
      chng(itc)= Rel(Func_opt,itc);  ier= Min(ier,ncl-1)
        
      If (ier < 0) then
        Success= -2;  Call Out ("Error in 'Converge_clusters': at iteration", itc, ln=1)
        Exit Iteration_loop

      Else if (ncl < n0) then
        key(:ncl)= key(LS_lst(:ncl));  key(ncl+1:)= 0
        Call Out ("Warning in 'Converge_clusters': # reg clusters reduced to", &
                  ncl, "at itertion",itc, ln=1)
        Call Out ("To-date retained clusters", key(:ncl))
      End if
        
      If (pr_out > 1.5) then
        Call Out ("Update_vec with # clusters",ncl, "at cycle",itc, ln=1)
        Call Out ("with convergence data",cnvI)
        If (dmp > 0) Call Out ("Damping data (parm,func,search dif)",Rel(:,itc))
      End if
        
!     Test for stagnation or divergence

      If (chng(itc) < 0.999*chng(itc-1)) then  ! Reduction success
        Success= 1;  itr= itc;  k= 0
        If (chng(itc) < tol(RS) .or. itc == mx) Exit Iteration_loop
        n0= ncl  ! Continue to next undamped iteration

        If (pr_out > 1) then
          Call Out ("undamped updates",chng(:itc))
        End if   
      Else
        Success= 0;  k= k + 1

        If (itc == mx .or. k >= 3) then
          itr= Minloc(chng(1:itc-1),1)
          If (RS == 3 .and. chng(itr) > tol(2)) Success= -1  ! Convergence shortfall
          Exit Iteration_loop
        Else
          If (dmp <= 0) rates(:,:n0,itc)= rates(:,:n0,itc-1) ! Back up & damp
          dmp0= dmp;  dmp= dmp1  
          n0= ncl  ! Continue to next damped iteration
        End if

        If (pr_out > 1) then
          Call Out (-1,"damping parm, update functional, net change",Rel(:,:itc))
        End if   
      End if
    End do Iteration_loop
    
    If (chng(itr) < tol(RS)) Success= 2
    n1= itr + Sum(nfcn(:itr))

    cnvI= (/ RS, Success, ncl, itr, n1 /);  cnvR= Rel(:,itr)

    If (dmp > 0) then
      nu= nc * ncl;  dg= cnvR(0)
      Call Jacobian_norms (nc,ncl,nu,dg, Memb,rates(:,:ncl,itr), Jnorm,Pnorm)

      If (pr_out > 1) then
        Call Out ("Jacobian norms: L_inf, L_1, Frobenius, L_2",Jnorm)
        Call Out ("preconditioned Jacobian norms: L_inf, L_1, Frobenius, L_2",Pnorm)
      End if
      x= Jnorm(4);  y= Pnorm(4)
    End if
    
    If (pr_out > 1) then
      Call Out ("For initial cluster set",iq, "at convergence cycle",icyc, ln=1)
      Call Out ("Restart phase, success code, # reg clusters, # iterations",cnvI(:4))
      Call Out ("Final damping parm, functional, net change", cnvR)

      If (dmp > 0) then
          n= Last_true(nfcn(:itr) == 0) 
        Call Out ("# stable iterations",n, "# update increases",k)
        Call Out ("Jacobian L2 norm",x, "Preconditioned L2 norm",y)
        Call Out (-1,"damping parm, functional, net change",Rel(:,:itc))

        If (Success < 0) then
          i= -1
        End if
      End if
    End if

    If (ier < 0) then
      Call Out ("Error in 'Converge_clusters': Convergence failed")
      Clust(0)%k= ier;  Clust(0)%n= ncl;  Return

    Else if (Success < 0) then
      Call Out ("Warning in 'Converge_clusters': Stagnation or divergence")
      Clust(0)%k= 0;  Clust(0)%n= ncl
      Return
    End if
    
    ind= ncl + 1
    Call Clusters_from_means (nc,ncl,-icyc, Memb, rates(:,:ncl,itr), Clust(:ind))
    
    ncl= Clust(0)%n;  Clust(0)%lt(:ncl)= key(Clust(0)%lt(:ncl)) 

  End Subroutine Converge_clusters

  
  Subroutine Update_vec (iter,ipc,dmp, nc,n0,rt_in, ncl,rt_out, Rel,nfn, ier)
  
    Integer,     Intent(in) :: iter         ! Called from this iteration
    Integer,     Intent(in) :: ipc          ! Print code
    Real,        Intent(in) :: dmp          ! Initial damping value
    
    Integer,     Intent(in) :: nc           ! # candidates
    Integer,     Intent(in) :: n0           ! Initial # clusters
    Real(Dblp),  Intent(in) :: rt_in(:,:)   ! (nc,n0) Initial cluster mean vectors
    
    Integer,    Intent(out) :: ncl          ! Final # clusters. Possibly reduced but not reordered
    Real(Dblp), Intent(out) :: rt_out(:,:)  ! (nc,n0) Final cluster mean vectors
    
    Real,       Intent(out) :: Rel(0:)      ! (0:2) Line search data
                                            !   0 = Final line search parameter
                                            !   1 = Final functional = update change
                                            !   2 = Net line search normalized (output - input) vectors
    Integer,    Intent(out) :: nfn          ! # line search steps = last index of 'tf'
                                            !   each step with 1 or 2 function calls (2 for a derivative)
                                            !   0 if no line search
    Integer,    Intent(out) :: ier          ! Success/error code
                                            !   < 0   Func_t internal error or ncl = 0
                                            !   = 0   Warning: Need to change algorithm features
                                            !   = 1   Warning: did not convergence in "Line_search
                                            !   > 1  Fully successful
! Local:
    Type(Multi_listD) :: PD(0:n0) !  Cluster membership and derivative data
    Type(Srch_parm)   :: parm
    Real(Dblp) :: tf(0:3,2,0:5)  ! (0,1,:) = damping parameter 't'
                                 ! (1,1,:) = functional = update change
                                 ! (1,2,:) = derivative of the functional
                                 ! (:,:,0) = initial evaluation
                                 ! (:,:,i) = line search evaluation 'i' for parameter 't'
    Real(Dblp) :: chng, ft(2), rt(nc,n0)
    Integer    :: lf, nl0, ng, ns, nt0, lst(n0)
    
    If (pr_out > 1) Call Out ("Enter 'Update_vec' from iteration",iter, ln=1)

    If (Allocated(LS_rate0)) DeAllocate(LS_rate0, LS_rate1, LS_rate2, LS_lst)
    Allocate(LS_rate0(nc,n0), LS_rate1(nc,n0), LS_rate2(nc,n0), LS_lst(n0))

    LS_rate0= rt_in;  LS_ncl= n0;  LS_lst= "ID"
    ns= Ubound(LS_Memb,1); ncl= n0;  Rel= 0;  nfn= 0;  tf= 0;  ier= 1
      
    If (dmp <= 0) then
      Call Allocate_PD (nc,n0,ns, PD)
      Call Generate_clusters (nc,n0,ns, LS_Memb, LS_rate0, PD, LS_rate1,chng)

        nGen_clust= nGen_clust + 1
      LS_ncl= PD(0)%n;  ncl= LS_ncl
      ier= PD(0)%k;  If (ncl /= n0) ier= -2

      rt_out= LS_rate1;  Rel(1:2)= chng 
      Return
    End if
    
    parm%ipc    = ipc
    parm%dmp_opt= 1
    parm%step   = 0.10d0
    parm%t_low  = 0.05d0
    parm%t0     = dmp
    parm%t_high = 0.50d0
    parm%tol    = 0.0025d0
    parm%change = 0.0025d0

    rt= rt_in;  lst= "ID";  nt0= nGen_clust

    Call GN_Func_0 (nc,n0, ft, ier)

    ncl= LS_ncl
    If (ier < 0) then
      nGen_clust_NT= nGen_clust_NT + (nGen_clust - nt0)
      Call Out ("Error in 'Update_vec': 'GN_Func_0' failure")
      Return
    Else if (ier == 0) then
      Call Out ("Warning in 'Update_vec': # clusters reduced by 'GN_Func_0' from", &
                n0, "to",ncl)
      lst(:ncl)= lst(LS_lst(:ncl))
    End if

    ncl= LS_ncl;  tf(1,1,0)= ft(1)
    
    If (pr_out > 1.5) then
      Call Out ("'Update_vec' after initial functional. Called at iter", &
                iter, ln=1)
      Call Out ("'Update_vec' initial result: functional",Real(ft(1)), &
                "gradient",Real(ft(2)))
    End if

    Call Line_search (parm, nfn,tf(:,:,1:), ier, GN_Func_t)

    If (ier < 0) then
      nGen_clust_NT= nGen_clust_NT + (nGen_clust - nt0)
      Call Out ("Error in 'Update_vec': 'Line_search' failure")
      Return
    Else if (ier == 0) then
      Call Out ("Warning in 'Update_vec': Reduced # clusters during 'Line_search'")
    End if
     
    ng= nGen_clust - nt0;  nGen_clust_NT= nGen_clust_NT + ng

    ncl= LS_ncl;  LS_lst(:ncl)= lst(:ncl)
    rt_out(:,:ncl)= LS_rate2(:,:ncl)

    Rel(0:1)= tf(0:1,1,nfn)
    Rel(2)  = Sqrt(Sum((LS_rate2(:,:ncl) - LS_rate0(:,:ncl))**2)/(nc*ncl))
    
    tf(2,1,:lf)= tf(1,2,:lf)
    
    If (pr_out > 1) then
      Call Out ("line search minimization of the functional completed")
      Call Out ("Final damping parm & residual, net line search change",Rel)
      Call Out (-1,"Line search: damping parm, residual & derivative",Real(tf(:2,1,:nfn)))
    End if
    
  End Subroutine Update_vec
  
  
   Subroutine Default_set (nc,np, Clust_set)
   
     Integer,           Intent(in) :: nc, np
     Type(Multi_listR), Intent(out):: Clust_set  ! Final cluster sets and associated data
                                                 ! %k  = first initial cluster set that converged to this cluster set
                                                 ! %l  = # clusters in the set, including independents = 'ind' = ncl + 1
                                                 ! %m  = # initial cluster sets that converged to this cluster set
                                                 ! %sum_wt: Significance level for ranking / rating adjustments
     
                                                 ! %px(6):  Objective data: min size, max correlation, 
                                                 !            independents size (frac), final residual
                                                 !            Jnorm, Pnorm
                                                 ! %qx(4):  Objective factors: min size penalty, max correlation penalty, 
                                                 !            independents size penalty, residual penalty
                                                 ! %rx(3):  Clustering objective value(1)= coherent vote(2) * penalty(3)
                                                 ! %sx(1):  Cluster sizes
                                                 ! %ux(6):   Fraction of voters represented to a given level
                                                 !           or higher with levels specified by Parm_rep

                                                 ! %vl(0:6): Convergence data
                                                 ! %L0(0:4,0:0) Integer convergence data per update call
                                                 ! %L1(-2:nc,0:1) Top and ordered candidate data, by cluster
                                                 ! %M2(0:4,0:0) Real convergence data
                                                 ! %M3(9,ind) Cluster data
                                                 ! %T0(nc,0:ind,3) Cluster point data
     
     Allocate (Clust_set%rx(3), Clust_set%px(6), Clust_set%qx(4), &
       Clust_set%ux(6), Clust_set%vl(0:6), Clust_set%L0(0:4,0:0), &
       Clust_set%L1(-2:nc,0:1), Clust_set%M2(0:4,0:0),            &
       Clust_set%M3(9,1), Clust_set%T0(nc,0:1,3), Clust_set%sx(1))

     Clust_set%k= -1;  Clust_set%l= 1;  Clust_set%m= 0
     Clust_set%rx= 0;  Clust_set%px= 0;  Clust_set%qx= 0
     Clust_set%ux= 0;  Clust_set%vl= 9;  Clust_set%L0= 0  
     Clust_set%M2= 0;  Clust_set%sx= np; Clust_set%M3= 0 
     Clust_set%M3(1,1)= np;  Clust_set%M3(2:3,1)= 1;  Clust_set%M3(4,1)= 0.5 
     Clust_set%M3(6,1)= 0.5;  Clust_set%M3(8:9,1)= 1

   End Subroutine Default_set
   
   
    Subroutine Allocate_clust (nc,ncl,ns, Clust)
    
!     Allocate and initialize cluster arrays
    
      Integer,            Intent(in) :: nc         ! # candidates
      Integer,            Intent(in) :: ncl        ! # regular clusters, assume > 0
      Integer,            Intent(in) :: ns         ! # slate clusters

      Type(Multi_listD), Intent(out) :: Clust(0:)  ! (0:ind) Cluster set data
                                                   ! For 0:
                                                   !   %k  = Error code
                                                   !   %m  = # candidates = nc
                                                   !   %n  = # regular clusters = ncl
                                                   !   %l  = ind = ncl + 1

                                                   !   %px(4):  Objective data 
                                                   !   %qx(4):  Objective penalty factors 
                                                   !   %rx(3):  Clustering obj value(1)= coherent vote(2) * penalty(3)
                                                   !   %sx(ind): Cluster sizes

                                                   !   %lt(n0) Cluster mapping

                                                   !   %vl(0:6): Convergence data
                                                   !   %wt(ns): Sum of slate cluster memberships, reduced by fuzzy factors

                                                   !   %M0(ns,0:1): Fuzzy set reduction factors (Generate_clusters)
                                                   !   %M1(ncl,ncl) Cluster correlation matrix
                                                   !   %M3(6,ind) Cluster data

                                                   !   %L1(-2:nc,0:ind) Candidate ordering data
                                                   !   %T0(nc,0:ind,3) Cluster point data

                                                   ! For clusters 1:ind
                                                   !   %k  = # fuzzy reduced slate memberships
                                                   !   %n  = # slate ballots which are partial or full members
                                                   !   %o   : 1 for regular clusters, 2 for independents
                                                   !   %fsx : slate ballot weight averaged fuzzy reduction factor
                                                   !   %fux : cluster width
                                                   !   %sum_wt = sum of slate memberships

                                                   !   %px(nc) = variance of %sx
                                                   !   %sx(nc) = mean rating vector
                                                   !   %tx(0:nc)= centered and normalized rating vector
                                                   !   %ux(nc)  = total cluster weight, with variance normalization
!   Local
      Integer :: i, k, n, ind
      
      n= Ubound(Clust,1) - 1
      If (n < ncl .or. ncl < 1) then
        Call Out ("Error in Allocate_clust: allocatable # reg clusters",n, &
                  "less than the required #",ncl);  Stop
      End if
      
      Clust(0)%k= 1;  Clust(0)%m= nc;  Clust(0)%n= ncl
      ind= ncl + 1;  Clust(0)%l= ind

      If (Associated(Clust(0)%px)) DeAllocate(Clust(0)%px, Clust(0)%qx,    & 
          Clust(0)%rx, Clust(0)%sx, Clust(0)%lt, Clust(0)%vl, Clust(0)%wt, &
          Clust(0)%M0, Clust(0)%M3, Clust(0)%L1, Clust(0)%T0)

      Allocate (Clust(0)%px(4), Clust(0)%qx(4), Clust(0)%rx(3),           & 
                Clust(0)%sx(ind), Clust(0)%lt(ncl), Clust(0)%vl(0:6),     &
                Clust(0)%wt(ns), Clust(0)%M0(ns,0:1), Clust(0)%M3(6,ind), &
                Clust(0)%L1(-2:nc,0:ind), Clust(0)%T0(nc,0:ind,3))

      Clust(0)%px= 0;  Clust(0)%qx= 1;      Clust(0)%rx= 0  
      Clust(0)%sx= 0;  Clust(0)%lt= "ID";   Clust(0)%vl= 0
      Clust(0)%wt= 0;  Clust(0)%M0(:,0)= 1; Clust(0)%M0(:,1)= 0
      Clust(0)%M3= 0;  Clust(0)%L1= 0;      Clust(0)%T0= 0

      If (ncl > 0) then
        Allocate(Clust(0)%M1(ncl,ncl));  Clust(0)%M1= 0
      End if 
      
      Do k= 1,ind
        Clust(k)%k= 0;  Clust(k)%m= 0;  Clust(k)%n= 0;  Clust(k)%o= 1 
        If (k == ind) Clust(k)%o= 2
        Clust(k)%fux= 0;  Clust(k)%sum_wt= 0

        If (Associated(Clust(k)%px)) DeAllocate (Clust(k)%px, Clust(k)%sx, &
                                                Clust(k)%tx, Clust(k)%ux)
        Allocate (Clust(k)%px(nc), Clust(k)%sx(nc), &
                  Clust(k)%tx(0:nc), Clust(k)%ux(nc))
        Clust(k)%px= 0;  Clust(k)%sx= 0
        Clust(k)%tx= 0;  Clust(k)%ux= 0
      End do
    End Subroutine Allocate_clust
    

    Subroutine Clusters_from_means (nc,nl,q, Memb, rate_in, Clust)

!     Compute the membership data for clusters represented by 
!     mean rating or ranking vectors and their normalized
!     difference vectors with the global mean vector (Euclidean norm
!     = 1, ready for computing correlations). Use fuzzy set
!     normalization of the memberships of the slate clusters
!     in the initial clusters.

!     Uses: Memb%fsx, Memb%sx, Memb%px, Memb%tx
!     Computes: Clust%n, Clust%sum_wt, Clust%ls, Clust%wt, Clust%sx
!               Clust%fux, Clust%lt, Clust%px

      Integer,  Intent(in) :: nc  ! # candidates
      Integer,  Intent(in) :: nl  ! Initial # regular clusters
      Integer,  Intent(in) :: q   ! For initial cluster set 'q' > 0, 
                                  ! else for convergence cycle -q
      
      Type(Multi_listD),  Intent(in) :: Memb(0:)       ! (0:ns) Slate clusters
      
      Real(Dblp),         Intent(in) :: rate_in(:,:)   ! (nc,nl) Mean vectors 
      
      Type(Multi_listD), Intent(out) :: Clust(0:)      ! (0:nl+1) Cluster data
                                                       ! For 0:
                                                       !   %k      = Error code. Invalid results if < 0.
                                                       !   %l      = # clusters, including independent = ind
                                                       !   %m      = # candidates = nc
                                                       !   %n      = # regular clusters = ncl

                                                       !   %lt(ncl) : original numbering of clusters if # clusters is reduced  ***

                                                       !   %px(6):  Objective data: min size, max correlation, min density, max boundary,
                                                       !            independents size (frac), final update
                                                       !   %qx(7):  Objective factors:  coherent vote, min size penalty, 
                                                       !            max correlation penalty, min density penalty, max boundary penalty, 
                                                       !            independents size penalty, final update penalty
                                                       !   %rx(3):  Cluster set objective value = coherent vote * penalty
                                                       !   %sx(ind): Cluster sizes  ***

                                                       !   %vl(0:6): Convergence data
                                                       !     0 = original # regular clusters
                                                       !     1 = final # regular clusters
                                                       !     2 = # convergence calls
                                                       !     3 = total # cluster merge operations
                                                       !     4 = total # cluster deletion operations
                                                       !     5 = total # centroiding iterations
                                                       !     6 = total # function calls
                                                       !   %wt(ns): Sum of slate cluster memberships over all reg
                                                       !            clusters, reduced by the fuzzy set factor,
                                                       !            to be <= 1.0 
      
                                                       !   %L0(0:4,0:n1) Integer convergence data per update call, 
                                                       !                 n1= %vl(2) - 1
                                                       !     0 : Restart type
                                                       !     1 : Success code
                                                       !     2 : # regular clusters
                                                       !     3 : # centroiding iterations
                                                       !     4 : # function calls
                                                       !   %M0(0:2,0:n1) Real convergence data
                                                       !     0 = Final damping parameter for the line search to minimize the residual
                                                       !     1 = Final residual = mean norm of output - input of 'Generate_clusters'
                                                       !     2 = Mean norm of of output - input of the line search

                                                       !   %M1(ncl,ncl) Final cluster correlation matrix  ***
                                                       !   %M3(6,ind) Cluster data ***
                                                       !     1 = size 
                                                       !     2 = density penalty
                                                       !     3 = boundary penalty
                                                       !     4 = density
                                                       !     5 = boundary
                                                       !     6 = width

                                                       !   %L1(-2:nc,0:ind) Candidate ordering data for original and zeroed ratings  ***
                                                       !     (-2:0,:) = # top, or significant, or viable candidates
                                                       !     (1:,:)   = Decreasing ordering for %T0(:,:,2)
                                                       !   %T0(nc,0:ind,3) Cluster point data  ***
                                                       !     (:,1:,1) = Cluster mean vectors
                                                       !     (:,1:,2) = Cluster mean vectors with zeroing
                                                       !     (:,1:,3) = Cluster portions for each candidate
                                                       !     (:,0,1:2)= Cluster size average of (:,1:,1:2)

                                                       ! For clusters 1:ind
                                                       !   %k  = # fuzzy reduced slate memberships
                                                       !   %n  = # slate ballots which are partial or full members
                                                       !   %o   : 1 for regular clusters, 2 for independents
                                                       !   %fsx : slate ballot weight averaged fuzzy reduction factor
                                                       !   %fux : cluster width
                                                       !   %sum_wt= sum of weighted slate memberships

                                                       !   %ls(n) = list of slate members
                                                       !   %wt(n)  = weighted slate memberships
                                                       !   %rx(n)  = unweighted slate memberships
                                                       !   %px(nc) = variance of %sx
                                                       !   %sx(nc) = mean rating vector
                                                       !   %tx(0:nc)= centered and normalized rating vector
                                                       !   %ux(nc)  = total cluster weight, incorporating slate ballot  
                                                       !              variance normalization (Memb%ux * %rx)
!   Local:
      Real(Dblp) :: por(nc,nl+1,2)  ! (nc,i0,2) Cluster set zeroed rating and and portion data
                                    !   1 = zeroed ratings,  2 = portions
      Real(Dblp) :: siz(nl+1), Zrate0(nc)
      Integer    :: key(nl)
      Logical    :: ReOrd
      Integer    :: cl, i1, n0, np, ns, ier, ind, ncl, ng0

      If (pr_out > 1) then
        If (q > 0) then
          Call Out ("Enter 'Clusters_from_means' for cluster set",q, &
                    "with # regular clusters",nl, ln=1)
        Else
          Call Out ("Enter 'Clusters_from_means' for convergence cycle",-q, &
                    "with # regular clusters",nl, ln=1)
        End if
      End if
      
      ng0= nGen_clust;  np= Memb(0)%l;  ns= Ubound(Memb,1)
      i1= nl + 1;  key= 0;  por= 0

      Call Allocate_clust (nc,nl,ns, Clust(:i1))

      Call Gen_full (np,nc,i1, Memb, rate_in(:,:nl), Clust(0:i1), &
                     Zrate0,por(:,:i1,:), ind, ier)

      ncl= ind - 1; Clust(0)%k= ier;  
      If (ncl < 1) Clust(0)%k= Min(Clust(0)%k, -1)
      ier= Clust(0)%k;  Clust(0)%n= ncl;  Clust(0)%l= ind

      If (ier < 0) then
        Call Out ("Error in 'Clusters_from_means': no clusters left")
        nGen_clust_CM= nGen_clust_CM + (nGen_clust - ng0)
        Return
      Else if (ncl < nl .and. pr_out > 1) then
        Call Out ("In 'Clusters_from_means', # clusters reduced from",nl, &
                  "to",ncl, ln= 1)
        Call Out ("Cluster reducion",Clust(0)%lt(:ncl))
      End if

      If (ncl > 1) then
        Call ReOrd_clusters (ncl, Clust(1:ncl), key(:ncl), ReOrd)

        If (ReOrd) then
          por(:,:ncl,:)= por(:,key(:ncl),:)
          Clust(0)%lt(:ncl)= Clust(0)%lt(key(:ncl))
          If (pr_out > 1) Call Out ("Cluster reducion and reordering",Clust(0)%lt(:ncl))
        End if
      End if

      Forall(cl=1:ind) siz(cl)= Clust(cl)%sum_wt

      Call Top_cand (np,nc,ind, siz(:ind), Zrate0,por(:,:ind,:), Clust(:ind))

      Call Cluster_set_width (nc, Clust(:ind))

      nGen_clust_CM= nGen_clust_CM + (nGen_clust - ng0)

      If (pr_out > 1) then
        Call Out ("Success in 'Clusters_from_means'")
        Call Out ("# top, significant, & viable cand's",Clust(0)%L1(-2:0,0))
        Call Out ("Cluster averaged candidate reordering",Clust(0)%L1(1:,0))
        Call Out (-1,"Cluster averaged ratings, full & zeroed",Real(Clust(0)%T0(:,0,1:2)))

        Call Out ("Final cluster sizes",Real(Clust(0)%sx(:ind)))
        Call Out (-1,"Zeroed ratings",Real(Clust(0)%T0(:,1:ind,2)))
        Call Out (-1,"and portions",Real(Clust(0)%T0(:,1:ind,3)))
      End if
    End Subroutine Clusters_from_means


   Subroutine Gen_full (np,nc,i0, Memb, rate_in, Clust, Zrate0,por, ind,ier)

     Integer,              Intent(in) :: np           ! # candidates to be elected
     Integer,              Intent(in) :: nc           ! # candidates
     Integer,              Intent(in) :: i0           ! initial # clusters, with independents
     Type(Multi_listD),    Intent(in) :: Memb(0:)     ! (0:ns)  Slate cluster data

     Real(Dblp),           Intent(in) :: rate_in(:,:) ! (nc,n0) Mean cluster vectors

     Type(Multi_listD), Intent(inout) :: Clust(0:)    ! (0:i0) Cluster data

     Real(Dblp),          Intent(out) :: Zrate0(:)    ! (nc) Cluster set zeroed rating 
     Real(Dblp),          Intent(out) :: por(:,:,:)   ! (nc,i0,2) Cluster set zeroed rating & portion data
                                                      !   1 = zeroed ratings, 2 = portions
     Integer,             Intent(out) :: ind          ! final # clusters,with ind's 
     Integer,             Intent(out) :: ier          ! error code

!    Local:
     Real(Dblp) :: tot_wt, change, siz(i0)
     Real(Dblp) :: rate0(nc,i0-1), rate1(nc,i0)
     Integer    :: cl, n0, nr, ns, ncl

     Zrate0= 0;  por= 0;  ind= i0;  ier= 1
     ncl= ind - 1;  ns= Ubound(Memb,1)

     If (pr_out > 1) Call Out ("Enter 'Gen_full': # reg clusters",ncl, ln=1)

     rate0= rate_in;  Clust(0)%lt= "ID"
     change= 0;  siz= 0;  rate1= 0
     
     Reduce_cycle : Do nr= 1,i0
       n0= ncl
       Call Generate_clusters (nc,n0,ns, Memb, rate0(:,:ncl), &
                               Clust(:ncl), rate1(:,:ncl),change)
         nGen_clust= nGen_clust + 1

       ncl= Clust(0)%n;  ind= ncl + 1;  ier= Clust(0)%k

       If (ncl < 1 .or. ier < 0) then  ! 'ncl' possibly reduced
         Call Out ("Error in 'Gen_full': Skip cluster set")
         Clust(0)%k= ier;  Clust(0)%n= ncl;  Return
       End if

       Call Indie_membership (nc,ncl, Clust(0)%wt,Memb, Clust(ind))

       siz(:ind)= Clust(1:ind)%sum_wt;  rate1(:,ind)= Clust(ind)%sx 
       tot_wt= Sum(siz(:ind))

       If (Clust(ind)%o == 2) Exit Reduce_cycle  ! Satisfactory cluster set

       If (ncl == 1) then
         Call Out ("Warning in 'Gen_full': Reduced to 1 reg cluster")
         Clust(0)%k= 0;  Exit Reduce_cycle
       End if
       
       If (pr_out > 1) Call Out ("Remove smallest cluster in 'Gen_full': ")
       cl= Minloc(siz(:ncl),1);  ind= ncl;  ncl= ind-1

       Clust(0)%lt(cl:ncl)= Clust(0)%lt(cl+1:ind); Clust(0)%lt(ind:)= 0
       rate0(:,cl:ncl)= rate0(:,cl+1:ind);  rate0(:,ind:)= 0
     End do Reduce_cycle

     Clust(0)%n= ncl;  Clust(0)%l= ind
     Clust(0)%sx(:ind)= siz(:ind);  Clust(0)%sx(ind+1:)= 0
     Clust(0)%sum_wt= Sum(siz(:ind))

     Call Portions (np,nc,Noise_por, siz(:ind),rate1(:,:ind), &
                    Zrate0, por(:,:ind,:))

   End Subroutine Gen_full

   Subroutine Top_cand (np,nc,ind, siz,Zrate0,por, Clust)

!    Determine the ordering of the candidates by the clusters,
!    individually and collectively (cluster aveaged mean vectors).
!    Also record data in Clust(0): %sx, %L1, %T0, and set
!    dimensions fro %M1 and %M3
   
     Integer,    Intent(in) :: np           ! # candidates to be elected
     Integer,    Intent(in) :: nc           ! # candidates
     Integer,    Intent(in) :: ind          ! # size of the cluster set
     Real(Dblp), Intent(in) :: siz(:)       ! (ind) Cluster sizes
     Real(Dblp), Intent(in) :: Zrate0(:)    ! (nc) Cluster averaged noise zeroed mean ratings
     Real(Dblp), Intent(in) :: por(:,:,:)   ! (nc,ind,2) Cluster set data
                                            ! (:,:,1) = noise zeroed ratings
                                            ! (:,:,2) = portions

     Type(Multi_listD), Intent(inout) :: Clust(0:) ! (0:ind) Cluster data
                                         ! For Clust(0) all values are output
                                         !   %k = Error code
                                         !   %m = # candidates = nc
                                         !   %n = # regular clusters = ncl
                                         !   %l = including independent = ind
                                         !   %sx(ind)  Cluster sizes  ***

                                         !   %M1(ncl,ncl) Cluster correlations  **
                                         !   %M3(6,ind)   Cluster data          **

                                         !   %L1(-2:nc,0:ind) Candidate ordering data for  ***
                                         !                    noise zeroed vectors
                                         !     (-2:0,:) = # top, significant, viable cand's
                                         !     (1:,1:)  = Decreasing ordering from %T0(:,1:,2)
                                         !     (1:,0)   = Decreasing ordering from %T0(:,0,1)
                                         !               
                                         !   %T0(nc,0:ind,3) Cluster vector data  ***
                                         !      (:,1:,1) = Cluster mean vectors
                                         !      (:,1:,2) = Cluster mean vectors with noise zeroing 
                                         !      (:,1:,3) = Cluster portion vectors
                                         !      (:,0,1)  = Zrate0 * declining regular portion factor
                                         !      (:,0,2)  = Zrate0 = cluster averaged noise zeroed mean vectors
!  Local:
     Integer :: cl, i1, n1, ncl

     ncl= ind - 1
     Clust(0)%k= 1;    Clust(0)%m= nc 
     Clust(0)%n= ncl;  Clust(0)%l= ind

     i1= Size(Clust(0)%sx);  n1= Ubound(Clust(0)%L1,1)
     If (i1 /= ind .or. nc /= n1) then
       DeAllocate (Clust(0)%sx, Clust(0)%M1, Clust(0)%M3, &
                   Clust(0)%L1, Clust(0)%T0)

       Allocate (Clust(0)%sx(ind),   Clust(0)%M1(ncl,ncl),     &
                 Clust(0)%M3(6,ind), Clust(0)%L1(-2:nc,0:ind), &
                 Clust(0)%T0(nc,0:ind,3))
     End if

     Clust(0)%sx= siz;  Clust(0)%M1= 0
     Clust(0)%M3= 0;    Clust(0)%M3(1,:)= siz

     Forall(cl=1:ind) Clust(0)%T0(:,cl,1)= Clust(cl)%sx
     Clust(0)%T0(:,1:,2:)= por  

     Call Zrf_M (nc, Zrate0, por(:,ind,2), Clust(0)%L1(:,0), Clust(0)%T0(:,0,1))

     Clust(0)%T0(:,0,2)= Zrate0;  Clust(0)%T0(:,0,3)= -1

     Call Best_cand (np,nc,ind, Clust(0)%T0(:,1:,2), Clust(0)%L1)

   End subroutine Top_cand

   Subroutine Zrf_MD (nc, Zrate0,por, Zrt1_ord, Zrt1)

!    Compute a modified Zrate0 and its ordering.

!    This is the cluster size average of the noise zeroed 
!    cluster mean vectors 'Zrate0' as scaled by (1+eps - por(:,ind)),
!    which is the same as Sum(por(:,:ncl)) + eps. Thus as
!    the independents become dominate (por(:,ind) ~ 1)
!    the scaled values Zrt1 are the Zrate0 values shrunk
!    down to 'eps'. 

!    This means that they will come in last in the 
!    Zrt1_ord ordering 'ls(1:)' of the candidates.

!    This ordering aides the pruning in the branch and bound
!    algorithm used by Fit_electeds. That is, the 0 or near 0
!    values of the portions of these candidates for the 
!    regular clusters, means they have no effect on 
!    proportionality. The 'eps' serves to keep the original
!    Zrate0 ordering among these.


     Integer,     Intent(in) :: nc            ! # candidates
     Real(Dblp),  Intent(in) :: Zrate0(:)     ! (nc) = Zrate0(:,0) = cluster averaged 
                                              !        mean vector of Zrate0(:,1:ind) 
     Real(Dblp),  Intent(in) :: por(:)        ! (nc) = por(:,ind)

     Integer,    Intent(out) :: Zrt1_ord(-2:) ! (-2:nc) Decreasing ordering of Zpor at (1:nc)
                                              !      with # top, significant, and viable at (-2:0)
     Real(Dblp), Intent(out) :: Zrt1(:)       ! (nc) Scaled Zrate0 = Zrate0 * [1+eps - por(:,ind)]
                                              !      for assignment algorithm viability.
!  Local:
     Real, Parameter :: eps= 0.005, sig_lev(3)= (/ 0.80, 0.20, eps /)  
                                    ! Fractions for top, significant, & viable 
                                    ! candidates as ordered by Zpor_ord(1:nc)
     Real    :: tmp(nc)
     Integer :: i, j
 
     Zrt1= Zrate0 * ((1+eps) - por)

     tmp= Zrt1;  Call Sort (.false., tmp, Zrt1_ord(1:), ez=eps)  ! Stable Sort
     tmp= tmp / tmp(1)

     j= -2;  Do i= 1,3
       Zrt1_ord(j)= Last_true(tmp > sig_lev(i));  j= j+1
     End do

   End Subroutine Zrf_MD

   Subroutine Zrf_MR (nc, Zrate0,por, Zrt1_ord, Zrt1)

!    Compute a modified Zrate0 and its ordering.

!    This is the cluster size average of the noise zeroed 
!    cluster mean vectors 'Zrate0' as scaled by (1+eps - por(:,ind)),
!    which is the same as Sum(por(:,:ncl)) + eps. Thus as
!    the independents become dominate (por(:,ind) ~ 1)
!    the scaled values Zrt1 are the Zrate0 values shrunk
!    down to 'eps'. 

!    This means that they will come in last in the 
!    Zrt1_ord ordering 'ls(1:)' of the candidates.

!    This ordering aides the pruning in the branch and bound
!    algorithm used by Fit_electeds. That is, the 0 or near 0
!    values of the portions of these candidates for the 
!    regular clusters, means they have no effect on 
!    proportionality. The 'eps' serves to keep the original
!    Zrate0 ordering among these.


     Integer,  Intent(in) :: nc            ! # candidates
     Real,     Intent(in) :: Zrate0(:)     ! (nc) = Zrate0(:,0) = cluster averaged 
                                           !        mean vector of Zrate0(:,1:ind) 
     Real,     Intent(in) :: por(:)        ! (nc) = por(:,ind)

     Integer, Intent(out) :: Zrt1_ord(-2:) ! (-2:nc) Decreasing ordering of Zpor at (1:nc)
                                           !      with # top, significant, and viable at (-2:0)
     Real,    Intent(out) :: Zrt1(:)       ! (nc) Scaled Zrate0 = Zrate0 * [1+eps - por(:,ind)]
                                           !      for assignment algorithm viability.
!  Local:
     Real, Parameter :: eps= 0.005, sig_lev(3)= (/ 0.80, 0.20, eps /)  
                                    ! Fractions for top, significant, & viable 
                                    ! candidates as ordered by Zpor_ord(1:nc)
     Real    :: tmp(nc)
     Integer :: i, j
 
     Zrt1= Zrate0 * ((1+eps) - por)

     tmp= Zrt1;  Call Sort (.false., tmp, Zrt1_ord(1:), ez=eps)  ! Stable Sort
     tmp= tmp / tmp(1)

     j= -2;  Do i= 1,3
       Zrt1_ord(j)= Last_true(tmp > sig_lev(i));  j= j+1
     End do

   End Subroutine Zrf_MR

   Subroutine Best_cand (np,nc,ind, Zrate, OrdZ)
   
!   Determine the top rated, significant, and viable candidates
!   for each cluster & overall. 
!   Required: 1 >  Parm_top(1) > Parm_top(2) > 0

     Integer,   Intent(in) :: np  ! # candidates to be elected
     Integer,   Intent(in) :: nc  ! # candidates
     Integer,   Intent(in) :: ind ! # clusters, including independents

     Real(Dblp), Intent(in) :: Zrate(:,:)  ! (nc,ind) Zeroed reg cluster 
                               ! mean ratings (when < noise level for Rating < 1). 
                                                   
     Integer, Intent(inout) :: OrdZ(-2:,0:) ! (-2:nc,0:ind) 
                               ! Candidate orderings at (1:,:) by declining Zrate, 
                               !   with (-2:0,:) = # top, significant, or viable
                               !   candidates for each ordering, except OrdZ(:,0) is input
!  Local:
     Real(Dblp) :: rt(nc)
     Logical :: Keep1(nc), Keep2(nc), Keep3(nc)
     Integer :: ord(nc), ko(3)
     Integer :: cl, cn, k1, k2, k3, n1, n2, n3
     
     If (Rating == 1) then
       Keep1= Pos_cand
     Else
       Keep1= Req_cand
     End if

!    Candidates to keep as required: top, significant, or viable
     
     ord= OrdZ(1:,0)
     Keep1(ord(:np))= .true.;  Keep2= Keep1;  Keep3= Keep1
     
!    Determine the top rated amd significant candidates for each cluster

     Cluster_loop : Do cl= 1,ind
       rt= Zrate(:,cl) 
       Call Sort (.false., rt, OrdZ(1:,cl))
       n3= Last_true(rt > 0)  ! Viable candidates

       n2= Last_true(rt(:n3) > Parm_top(2)*rt(1)) ! significant candidates
       n1= last_true(rt(:n2) > Parm_top(1)*rt(1)) 
       n1= Min(n1,np)          ! top candidates, so n1 <= n2 <= n3 

       OrdZ(-2,cl)= n1;  OrdZ(-1,cl)= n2;  OrdZ(0,cl)= n3
       Keep1(OrdZ(1:n1,cl))= .true.;  Keep2(OrdZ(1:n2,cl))= .true.
       Keep3(OrdZ(1:n3,cl))= .true.
     End do Cluster_loop

     ko(1)= Last_true(Keep1(ord));  ko(2)= Last_true(Keep2(ord))
     ko(3)= Last_true(Keep3(ord))

     If (pr_out > 1) then
       Call Out ("Overall prior candidate ordering",ord)
       Call Out ("with prior # top ranked, significant, & viable",OrdZ(:0,0))
       Call Out ("vs new # top ranked, significant, & viable",ko)
   
       Call Out (-1,"The candidate orderings by cluster",OrdZ(1:,1:))
       Call Out (-1,"# top ranked, significant, & viable",OrdZ(:0,1:))
     End if

   End Subroutine Best_cand


   Subroutine Cluster_corr (nc,Adj,rank_pt, ballot,Clust, Borda, Corr,mx_corr)

!    Compute a cluster correlation matrix based on top 'mr' Borda
!    Counts, using a 0.5 increase in the Borda increment.
  
     Integer,        Intent(in) :: nc           ! # candidates
     Real,           Intent(in) :: Adj          ! Noise adjustment for mean vectors before
                                                !   correlating
     Real,           Intent(in) :: rank_pt(:)   ! (mr) Point conversion of ranking levels
     Integer,        Intent(in) :: ballot(0:,:) ! (0:mr,nb) Ballot ranking structure for clustering
     Type(Set_list), Intent(in) :: Clust(:)     ! (np) Cluster membership 
                                                !  %n     = # ballot members of the cluster
                                                !  %smb   = Cluster weight = size = Sum(val)
                                                !  %set(n)= List of ballot members
                                                !  %val(n)= Ballot weights = sizes

     Real,       Intent(out)    :: Borda(:,:)   ! (nc,np) Borda Count for each candidate cluster
     Real,       Intent(out)    :: Corr(:,:)    ! (np,np) Correlation matrix between pairs of clusters, 
                                                !   except = 0 on the diagonal
     Real,       Intent(out)    :: mx_corr(0:)  ! (0:np).  1:np = max cluster correlations by row
                                                !    0 = overall max = Maxval(mx_corr(1:))   
! Local:
     Real    :: bd(nc), bd2(nc)
     Real    :: wb, wt, norm
     Integer :: b, i, n, c1, c2, cl, nm, np

     Call Out ("Enter 'Cluster_corr'")

     np= Size(Corr,2);  Borda= 0;  mx_corr= 0;  Corr= -1

     Do cl= 1,np
       wt= 0;  nm= Clust(cl)%n
       Do i= 1,nm
         b= Clust(cl)%set(i);  wb= Clust(cl)%val(i)
           n= ballot(0,b)
         Borda(ballot(1:n,b),cl)= Borda(ballot(1:n,b),cl) + rank_pt(:n) * wb
         wt= wt + wb
       End Do 

       Borda(:,cl)= Borda(:,cl) / Max(wt,0.0001)
     End do

     Do c1= 1,np-1
       bd= Borda(:,c1) - Adj;  norm= Sqrt(Sum(bd**2))
       bd= bd / norm

       Do c2= c1+1,np
         bd2= Borda(:,c2) - Adj;  norm= Sqrt(Sum(bd2**2))
         Corr(c1,c2)= Dot_product(bd,bd2) / norm
         Corr(c2,c1)= Corr(c1,c2)
       End do
     End do

     mx_corr(1:)= Maxval(Corr,2);  mx_corr(0)= Maxval(mx_corr(1:))
     Forall(c1=1:np) Corr(c1,c1)= 0
   End Subroutine Cluster_corr

   Subroutine Cluster_corr2 (nc,ptval, ballot,ballot2, Clust, Borda, Corr,mx_corr, Adj)

!    Compute a cluster correlation matrix based on top 'mr' Borda
!    Counts, using a 0.5 increase in the Borda increment.
  
     Integer,        Intent(in) :: nc            ! # candidates
     Real,           Intent(in) :: ptval(:)      ! (mt) Points for rating levels

     Integer,        Intent(in) :: ballot(0:,:)  ! (0:mr,nb) Ranked or rated candidates
                                                 !           with (0,b) = # ranked or rated
     Integer,        Intent(in) :: ballot2(0:,:) ! (0:mr,nb) Corresponding ranking or rating levels, 
                                                 !           increasing, with (0,b) = # positive levels 

     Type(Set_list), Intent(in) :: Clust(:)      ! (np) Cluster membership 
                                                 !  %n     = # ballot members of the cluster
                                                 !  %smb   = Cluster weight = size = Sum(val)
                                                 !  %set(n)= List of ballot members
                                                 !  %val(n)= Ballot weights = sizes

     Real,       Intent(out)    :: Borda(:,:)    ! (nc,np) Borda Count for each candidate cluster
     Real,       Intent(out)    :: Corr(:,:)     ! (np,np) Correlation matrix between pairs of clusters, 
                                                 !   except = 0 on the diagonal
     Real,       Intent(out)    :: mx_corr(0:)   ! (0:np).  1:np = max cluster correlations by row
                                                 !    0 = overall max = Maxval(mx_corr(1:)) 

     Real, Optional, Intent(in) :: Adj           ! Noise adjustment for mean vectors before
                                                 !   correlating if ptval is all positive
! Local:
     Real    :: bd(nc), bd2(nc)
     Real    :: wb, wt, norm
     Integer :: b, i, c1, c2, cl, n, nm, np

     If (pr_out > 1) Call Out ("Enter 'Cluster_corr2'")

     np= Size(Corr,2);  Borda= 0;  mx_corr= 0;  Corr= -1

     Do cl= 1,np
       wt= 0;  nm= Clust(cl)%n
       Do i= 1,nm
         b= Clust(cl)%set(i);  wb= Clust(cl)%val(i);  n= ballot(0,b)  
         Borda(ballot(1:n,b),cl)= Borda(ballot(1:n,b),cl) + ptval(ballot2(1:n,b)) * wb
         wt= wt + wb
       End Do 

       Borda(:,cl)= Borda(:,cl) / Max(wt,0.0001)
     End do

     Do c1= 1,np-1
       bd= Borda(:,c1);  If (Present(Adj)) bd= bd - Adj 
       norm= Sqrt(Sum(bd**2));  bd= bd / norm

       Do c2= c1+1,np
         bd2= Borda(:,c2);  If (Present(Adj)) bd2= bd2 - Adj 
         norm= Sqrt(Sum(bd2**2));  bd2= bd2 / norm

         Call Dot_product_M (Real(Dot_fac), bd,bd2, i,Corr(c1,c2))
         Corr(c2,c1)= Corr(c1,c2)
       End do
     End do

     mx_corr(1:)= Maxval(Corr,2);  mx_corr(0)= Maxval(mx_corr(1:))
     Forall(c1=1:np) Corr(c1,c1)= 0
   End Subroutine Cluster_corr2

  Subroutine Cluster_overlap (Clus, Overlap, mx_ovr)

!   Compute the cluster overlap matrix.
  
!   Measure overlap of two clusters by summing the minimum 
!   of the memberships of the ballots in the two clusters, then
!   normalizing by the minimum the two cluster weights. 
  
    Type(Set_list), Intent(in) :: Clus(:)       ! (ncl) Membership for the deleted & merged 
                                                !       clusters with sizes %smb decreasing 
                                                !  %n     = # ballot members of the cluster
                                                !  %smb   = Cluster size = Sum(val)
                                                !  %set(n)= List of ballot members (any order)
                                                !  %val(n)= Ballot weights = sizes
    Real,          Intent(out) :: Overlap(:,:)  ! (ncl,ncl) Overlap between pairs of clusters, 
                                                !           except for 0 on the diagonal
    Real,          Intent(out) :: mx_ovr(0:)    ! (0:ncl) 1:ncl = max row = column overlaps
                                                !         0 = overall max = Maxval(mx_ovr(1:))   

! Local:
    Integer, Allocatable :: Intr(:), Sub1(:), Sub2(:)
    Real    :: add
    Integer :: i, m, c1, c2, i1, i2, n1, n2, ni, ncl

    ncl= Size(Clus);  Overlap= 0;  mx_ovr= 0;  If (ncl <= 1) Return
    
    Candidate_loop1 : Do c1= 1,ncl-1
      Candidate_loop2 : Do c2= c1+1,ncl

        n1= Clus(c1)%n;  n2= Clus(c2)%n;  m= Min(n1,n2)
        Allocate(Intr(m), Sub1(m), Sub2(m))

        Call Two_set_analysis (Clus(c1)%set,Clus(c2)%set, ni,Intr, Sub1,Sub2)

        Do i= 1,ni
          i1= Sub1(i);  i2= Sub2(i)
          add= Min(Clus(c1)%val(i1), Clus(c2)%val(i2))
          Overlap(c1,c2)= Overlap(c1,c2) + add
        End do

        DeAllocate(Intr, Sub1, Sub2)
      End do Candidate_loop2
    End do Candidate_loop1

!   Overlap normalization

    Overlap_loop : Do c1= 1,ncl-1
      Do c2= c1+1,ncl
        Overlap(c1,c2)= Overlap(c1,c2) / Min(Clus(c1)%smb, Clus(c2)%smb)
        Overlap(c2,c1)= Overlap(c1,c2)
      End do
    End do Overlap_loop

    mx_ovr(1:)= Maxval(Overlap,2);  mx_ovr(0)= Maxval(mx_ovr(1:))
    
  End Subroutine Cluster_overlap


   Subroutine Pref_clust (np,portion, elect)

!    Compute the preferential order elect(:,1) of the 
!    elected candidates elect(:,0), also the
!    represented cluster sets elect(:,2)

     Integer,    Intent(in) :: np            ! # candidates to be elected
     Real,       Intent(in) :: portion(:,:)  ! (nc,ind) Portions for each cluster

     Integer, Intent(inout) :: elect(:,0:)   ! (np,0:2) 
!  Local:
     Real    :: por(np)
     Integer :: ls(np), clx(np), elc(np), cls(np), key(np)
     Logical :: ReOrd
     Integer :: j, k, n, cl, k1, kn, ind

!    Order by the elected candidates elect(:,0) first by the
!    cluster clx(j) of the top portion for each elect(j,0).

!    Where a cluster is represented by more than one candidate
!    order them by decreasing portion.

     ind= Size(portion,2)
     Forall(j=1:np) clx(j)= Maxloc(portion(elect(j,0),:),1)

     k= 0;  elc= 0;  cls= 0;  ls= 0

     Do cl= 1,ind
       Call List_of_true (clx == cl, n,ls)  ! Elected indices for cluster 'cl'
       If (n <= 0) Cycle;  k1= k + 1

       If (n > 1) then
         kn= k + n;  cls(k1:kn)= cl
         elc(k1:kn)= elect(ls(:n),0)
         por(:n)= portion(elc(k1:kn),cl)

         Call Sort (.false., por(:n),key(:n), ReOrd)
         If (ReOrd) then
           key(:n)= k + key(:n);  elc(k1:kn)= elc(key(:n))
         End if
         k= kn
       Else
         cls(k1)= cl;  elc(k1)= elect(ls(1),0);  k= k1
       End if
     End do

     elect(:,1)= elc;  elect(:,2)= cls

     If (k /= np) then
       Call Out ("Error in 'Pref_clust': Assignments",k, "but should be",np, ln=1)
       Stop
     End if

   End Subroutine Pref_clust

End Module Clusters_support
