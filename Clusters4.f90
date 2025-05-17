
!            This module contains "Evaluate_candidates" and its subroutines. 
    
!  "Evaluate_candidates" is called from Clustering_PR to evaluate different possible
!  sets of elected candidates as to how well they match the cluster sets
    
  Module Clusters4

   Use Clusters0
   Use Clusters_support
  
   Use Misc_methods
   Use Graph_algorithms
   Use Factorials
   Use Util
   Use Output
   Use Types
   Use Precisn
   Use IEEE_Arithmetic
   Implicit None

 Contains

   Subroutine Evaluate_candidates (idist, nc,np,ncs, Non_clust,Clust_set, &
                                   Clust_elect, me,Elect_dat, dpr)

!    Compute the feasible elected sets of candidates and how well they match each cluster set.
    
     Integer,              Intent(in) :: idist          ! distict #
     Integer,              Intent(in) :: nc             ! # candidate
     Integer,              Intent(in) :: np             ! # candidate to be elected
     Integer,              Intent(in) :: ncs            ! # cluster sets in 'Clust_set' and 'Clust_elect'

     Type(Multi_listR), Intent(inout) :: Non_clust  ! Data for 'nMt' selected non-clustering methods 'm'
                                         !  to elect multiple candidates.
                                         ! %n = # distinct sets of electeds found for the 
                                         !      current district, = 'nonC'
                                         ! %k = final # STV clusters, = 'nSTV'
                                         ! %l = final # DTV clusters, = 'nDTV'
                                         ! %vl(nMt)  = The distinct elected set in %L0 for each method, = 'Non_map'

                                         ! %L0(np,nonC) = Distinct non-clustering sets of electeds in 
                                         !                increasing order
                                         ! %Q0(nc,nMt,2)= For each method, all candidates (:,m,1) in ranked order
                                         !                so that (:np,m,1) = the electeds in rank order,
                                         !                with (:np,m,2) = electeds in increasing order, = 'Non_elc'
                                         ! %L1(nMt,2)  (m,1)= Index 'e' in Elect_dat of non-clustering method 'm' 
                                         !             (m,2)= The matching best cluster set 'q', = 'Non_dat'
                                         ! %L2(nonC,ncs) = Location, for each cluster set 'q', in Clust_elect(q)%Q0 
                                         !                 of each distinct non-clustering set of electeds, = 'Non_loc'

                                         ! %T0(0:nc,ind,2) = Cluster data from STV (:,:,1) and DTV (:,:,2).
                                         !    (0,:,:) = cluster weighs, (1:nc,:,:) = cluster mean rating vectors,
                                         !    = 'STV_mean' & 'DTV_mean' 
                                         ! %T1(ncl,ncl,2) = Cluster correlations for STV (:,:,1) and DTV (:,:,2),
                                         !    = 'STV_corr' & 'DTV_corr' .
     Type(Multi_listR), Intent(inout) :: Clust_set(:) ! (ncs) Final cluster sets
                                         ! %k   = best initial cluster set that converged to this cluster set
                                         ! %l   = # clusters with independents = 'ind'
                                         ! %m   = # initial cluster sets that converged to this cluster set
                                         ! %lt(m)= list of those initial sets 
                                         ! %fux = clustering objective ratio = (q)%rx(1) / (1)%rx(1) 

                                         ! %lt(ncl) Mapping of initial to converged regular clusters 

                                         ! %px(6):  Objective data: min size, max correlation, 
                                         !            independents size (frac), final residual,
                                         !            Jnorm, Pnorm
                                         ! %qx(4):  Objective factors: min size penalty, max correlation penalty, 
                                         !            independents size penalty, residual penalty
                                         ! %rx(3):  Objective value(1) = coherent vote(2) * penalty(3)
                                         ! %sx(ind): Cluster sizes

                                         ! %ux(6)  Frac of voters represented by regular clusters
                                         !         when that representation exceeds 'Parm_rep'

                                         ! %vl(0:6): Convergence data
                                         !     0 = original # regular clusters = n0
                                         !     1 = final # regular clusters
                                         !     2 = final success code
                                         !     3 = total # cluster merge and deletion operations
                                         !     4 = # convergence calls
                                         !     5 = total # centroiding iterations
                                         !     6 = total # line search points

                                         ! %L0(0:4,0:iv) Integer convergence data per update call, n1= %vl(4)
                                         !    0  = restart phase
                                         !    1  = success code:  2 = fully converged,   1 = reducing functional,
                                         !                        0 = damped reduction, -1 = convergence failure
                                         !    2 = # reg clusters
                                         !    3 = # centroiding iterations
                                         !    4 = # line search points

                                         ! %M2(0:2,0:n1) Real convergence data
                                         !    0: Final line search parameter 't'
                                         !    1: Final residual
                                         !    2: Net line search change

                                         ! %M1(ncl,ncl) cluster correlation matrix, 0 on diagonal

                                         ! %L1(-2:nc,0:ind) Candidate ordering data for
                                         !                    noise zeroed ratings
                                         !   (-2:0,:) = # top, significant, viable cand's
                                         !   (1:,1:)  = Decreasing ordering from 
                                         !              zeroed ratings = %T0(:,1:,2)
                                         !   (:,0)   = Ordering & data from %T0(1:,0,1)

                                         ! %T0(nc,0:ind,3) Cluster point data
                                         !    (:,1:,1) = Cluster mean vectors
                                         !    (:,1:,2) = Cluster mean vectors with noise zeroing 
                                         !    (:,1:,3) = Cluster portions for each candidate
                                         !    (:,0,1) = Scaled Zrate0 (for ordering the %L1(1:,0))
                                         !    (:,0,2) = Zrate0 = cluster averaged noise zeroed mean ratings

     Type(Multi_listR),   Intent(out) :: Clust_elect(:) ! (ncs) Cluster set data for best sets of electeds
                                         ! %k = 'e' where (q)%Q0(:,0,1) = Elect_dat(e)%ls
                                         ! %l = 'ind' = # clusters, including independents
                                         ! %m = 'nf'= # sets of possible electeds
                                         ! %fsx = full objective for the top set of electeds 
                                         !        %Q0(:,0,1) = Clust(q)%sx(1)
                                         ! %fux = full objective ratio = (q)%fsx / (1)%fsx 

                                         ! %sx(7)  Full objective data electeds %k, as in Elect_dat%M1 

                                         ! %ls(nonC) Locations in %Q0(:,0,:) of the non-clustering sets 
                                         ! %vl(nf) The indices in 'Elect_dat' of the 'nf' possible sets 
                                         !         of the electeds listed in 'Q0'

                                         ! %L0(np,0:2) Rem_frac electeds
                                         !   0 = set of electeds - increasing order
                                         !   1 = ordered by decreasing cluster size + ratings
                                         !   2 = their corresponding clusters

                                         ! %L1(np,0:2) DHondt electeds
                                         !   0 = set of electeds - increasing order
                                         !   1 = ordered by decreasing cluster size + ratings
                                         !   2 = their corresponding clusters

                                         ! %Q0(np,0:2,nf)= Best sets of electeds by decreasing fitness 
                                         !     (:,0,i) = electeds in increasing order 
                                         !     (:,1,i) = electeds in preferential order 
                                         !     (:,2,i) = the clusters that represent (:,1,i) 
                                         ! %M2(3,nf)= Electeds fitness(1) = average rating(2) * deviation penalty(3) 

                                         ! %T0(ind,2,nf) Objective data for each set of possible electeds
                                         !   1 = Cluster size deviations = actual - true
                                         !   2 = Cluster size deviation penalties

     Integer,             Intent(out) :: me  ! # possible sets of electeds tested

     Type(Multi_listR),   Intent(out) :: Elect_dat(:) ! (mxe=>me) Data for possible sets 'e' of electeds,
                                         ! ordered by 'best' full objective %M1(1,1)
                                         ! %k = The 'best' cluster set 'q' = %lt(1)
                                         ! %l = The matching 'eq' in Clust_elect(q)%Q0= %vl(1)
                                         ! %m = nq = # better cluster sets that list 'e' 
                                         !      or = 1 for lower cluster sets
                                         ! %n = # regular clusters for best cluster set 'q'
                                         ! %q = Non-clustering index if 'e' is a
                                         !      non-clustering elected set, else 0

                                         ! %fsx = overall satisfaction
                                         ! %fux = full objective ratio = (e)%M1(1,1) / (1)%M1(1,1) 

                                         ! %ls(np)    The set of electeds in increasing order
                                         ! %L0(np,nq) The set in preferential order for each
                                         !            of the cluster sets ordered as in %lt

                                         ! %lt(nq) Ordering of the 'nq' cluster sets 'q' which  
                                         !         list the set of electeds 'e' at a position 
                                         !         'eq' in Clust_elect(q)%Q0(:,:,eq) with eq <= nf. 
                                         !         The ordering is by full objective %M1(1,i) 
                                         !         where %lt(i) = q.
                                         ! %vl(nq) The corresponding indices 'eq'
     
                                         ! %M1(7,nq) Objective data for the clusters sets ordered %lt.
                                         !   (1,i) = full objective = (2,i)*(3,i)
                                         !   (2,i) = clustering objective = (4,i)*(5,i)
                                         !   (3,i) = electeds fitness = (6,i)*(7,i)
                                         !   (4,i) = cluster set coherent vote 
                                         !   (5,i) = cluster set penalty
                                         !   (6,i) = elected rating = average rating over electeds  
                                         !   (7,i) = deviation size penalty
     
                                         ! %M0(0:Mvt,2): Measures of success - over ballots
                                         !   (0,1) : mean # of top 3 elected in order (but not the next lower)
                                         !   (0,2) : mean # of top 3 elected regardless of order
                                         !   (l,1) : fraction with top l or more elected in order
                                         !   (l,2) : fraction with top l or more elected regardless of order

     Integer, Intent(out) :: dpr(0:)     ! (0:mfe) Elected sets whose best full objective 
                                         !         exceeds that of cluster set #1.
                                         !         (0) = # such elected sets.
! Local:
     Real    ::  fitness(3,mfe)  ! Elected set data for each cluster set

     Real    :: tmp(ncs)
     Integer :: nl(ncs), key(ncs), ord(ncs), Non_eq(Non_clust%n,2)
     Integer :: elect(np,0:2,mfe)

     Real,    Allocatable :: elc_dat(:,:,:), mean(:,:), elc_obj(:)
     Integer, Allocatable :: elc_set(:,:)

     Real    :: r1, r2
     Integer :: cl, eq, ml, nf, mq, mx, ind, ncl, nf0, nonC
     Integer :: e, i, j, k, l, m, n, q
    
     Call Out ("Enter 'Evaluate_candidates': Possible electeds by cluster set for district",idist,ln=1)
     
     nonC= Non_clust%n;  nl= Clust_set%l - 1;  ml= Maxval(nl) + 1
     
     me= 0;  dpr= 0;  Non_clust%L1= 0;  Non_clust%L2= 0

!    Compute possible sets of electeds

     Do q= 1,ncs
       ind= Clust_set(q)%l;  Clust_elect(q)%l= ind
       If (Allocated(elc_dat)) DeAllocate(elc_dat)  
       Allocate (elc_dat(ind,2,mfe))
       
       If (pr_out >= 1)  Call Out ("For cluster set",q, ln=1)

       Call Elect_best (q,nc,np,ind,Clust_set(q)%L1(:,0), Clust_set(q)%sx, &
                        Clust_set(q)%T0(:,0,2), Clust_set(q)%T0(:,1:,3),   &
                        nf0,elect, fitness, elc_dat)

       mx= Min(nf0+nonC, mfe);  nf= nf0

       Call Add_non_clust (q,np,ind,mx,nonC,Non_clust%L0,Clust_set(q)%sx,   &
                           Clust_set(q)%T0(:,0,2), Clust_set(q)%T0(:,1:,3), &
                           nf,elect(:,:,:mx), fitness(:,:mx),               &
                           elc_dat(:,:,:mx), Non_clust%L2(:,q))

       If (Associated(Clust_elect(q)%Q0)) then
         DeAllocate(Clust_elect(q)%sx, Clust_elect(q)%ls, Clust_elect(q)%vl, &
                     Clust_elect(q)%Q0,Clust_elect(q)%M2, Clust_elect(q)%T0)
       End if
       Allocate(Clust_elect(q)%sx(7),    Clust_elect(q)%ls(nonC),      &
                Clust_elect(q)%vl(nf),   Clust_elect(q)%Q0(np,0:2,nf), &
                Clust_elect(q)%M2(3,nf), Clust_elect(q)%T0(ind,2,nf))
       
       Clust_elect(q)%m  = nf
       Clust_elect(q)%fsx= Clust_set(q)%rx(1) * fitness(1,1)
       Clust_elect(q)%fux= Clust_elect(q)%fsx / Clust_elect(1)%fsx
       Clust_elect(q)%ls = Non_clust%L2(:,q)
       Clust_elect(q)%vl = 0

       Clust_elect(q)%M2= fitness(:,:nf)

       Clust_elect(q)%sx(1)  = Clust_elect(q)%fsx       
       Clust_elect(q)%sx(2)  = Clust_set(q)%rx(1)       
       Clust_elect(q)%sx(3)  = fitness(1,1)    
       Clust_elect(q)%sx(4:5)= Clust_set(q)%rx(2:3)       
       Clust_elect(q)%sx(6:7)= fitness(2:3,1)

       Clust_elect(q)%Q0= elect(:,:,:nf)
       Clust_elect(q)%T0= elc_dat(:,:,:nf)
     End do

!    Compile a list of all 'significant' sets of electeds 'elc_set'
!    - those from the non-clustering algorithms and all those 
!    which scored well for some cluster set.

     Allocate(elc_set(-1:np,mxe), elc_obj(mxe))
     
     Call Merge_electeds (ncs,Clust_elect, me,elc_set,elc_obj)

!    Non-clustering data

     Elect_dat(:me)%q= 0;  Non_eq= 0

     Do m= 1,nMt
       i= Non_clust%vl(m);  If (i < 1 .or. i > nonC) Cycle

       If (Non_eq(i,1) < 1) then
         eq= Clust_elect(1)%ls(i);  e= Clust_elect(1)%vl(eq)
         Elect_dat(e)%q= i;  q= elc_set(-1,e)
         Non_eq(i,1)= e;  Non_eq(i,2)= q
       End if

       Non_clust%L1(m,:)= Non_eq(i,:)
     End do

!    Compute objective data Elect_dat%M1 for each significant 
!    set of electeds: for the best overall cluster set, the top 
!    cluster set for the set of electeds, and the top cluster 
!    sets for each number of regular clusters.

     Do e= 1,me
       Allocate(Elect_dat(e)%ls(np), Elect_dat(e)%M0(0:Mvt,2))
       q= elc_set(-1,e);  eq= elc_set(0,e)
       Elect_dat(e)%k= q;  Elect_dat(e)%l= eq
       Elect_dat(e)%n= Clust_elect(q)%l - 1

       Elect_dat(e)%fux= elc_obj(e) / elc_obj(1)
       Elect_dat(e)%ls = elc_set(1:,e);  Elect_dat(e)%M0= 0
     End do
     
     If (pr_out > 1) then
       Call Out ("Cluster set full objectives",Clust_elect%fsx)
       Call Out ("Corresponding # regular clusters", nl)

         Forall(q=1:ncs) tmp(q)= Clust_set(q)%rx(1)
       Call Out ("Corresponding cluster objective",tmp)
         Forall(q=1:ncs) tmp(q)= Clust_elect(q)%M2(1,1)
       Call Out ("Corresponding best fitness",tmp)
         Forall(q=1:ncs) ord(q)= Clust_elect(q)%vl(1)
       Call Out ("Corresponding best elected set",ord)

       Call Out ("Total # elected sets",me)
       Call Out ("# elected sets associated to each cluster set", &
                 Clust_elect%m)
       Call Out (-1,"Elected set & best cluster set, for non-clustering methods", Non_eq)
       Call Out ("Ratio non-clustering to best clustering full objectives", &
                 Elect_dat(Non_eq(:,1))%fux)
     End if

     Call Organize_electeds (nc,np,ncs,me, Clust_elect, Elect_dat(:me), dpr)

!    Candidates elected by largest remaining fraction method & D'Hondt
     
     Do q= 1,ncs
       If (Associated(Clust_elect(q)%L0)) &
           DeAllocate(Clust_elect(q)%L0, Clust_elect(q)%L1)
       Allocate(Clust_elect(q)%L0(np,0:2), Clust_elect(q)%L1(np,0:2))
     End do
     Allocate (mean(0:nc,ml));  mean= 0

     Do q= 1,ncs
         ind= nl(q) + 1
       mean(0,:ind) = Clust_set(q)%sx(:ind)       ! cluster sizes
       mean(1:,:ind)= Clust_set(q)%T0(:,1:ind,1)  ! cluster mean rating vectors
       
       Call Rem_frac (q,nc,ind, mean(:,:ind), Clust_elect(q)%L0)
       Call DHondt   (q,nc,ind, mean(:,:ind), Clust_elect(q)%L1)
     End do

   End Subroutine Evaluate_candidates

                          
   Subroutine Elect_best (q,nc,np,ind,ord, siz,Zrate0, &
                          portion, nf,elect,fitness,elc_dat) 

!    Determine the top sets of 'np' candidates which are the 
!    best match for the cluster set characterized by Zrate0, 
!    siz, and ord. These will be up to the top 'mfe' 
!    such sets as measured by fitness.

!    These have the highest overall cluster rankings, 
!    where these rankings are weighted by portions 
!    which reflect candidate memberships in the clusters.
     
     Integer,    Intent(in) :: q              ! Cluster set
     Integer,    Intent(in) :: nc             ! # candidates
     Integer,    Intent(in) :: np             ! # candidates to be elected
     Integer,    Intent(in) :: ind            ! # clusters, including independents

     Integer,    Intent(in) :: ord(-2:)       ! (-2:nc) (1:nc) = Ordering of Zrate0 scaled
                                              !     by Sum(portion(:,:ncl),2)
     Real,       Intent(in) :: siz(:)         ! (ind) Cluster sizes
     Real,       Intent(in) :: Zrate0(:)      ! (nc) Cluster-averaged mean vectors with zeroing
     Real,       Intent(in) :: portion(:,:)   ! (nc,ind) Portion of each cluster for each candidate
     
     Integer,   Intent(out) :: nf             ! # matching elected sets found
     Integer,   Intent(out) :: elect(:,0:,:)  ! (np,0:2,mfe) Corresponding sets of electeds, 
                                              !            using the overall fitness ordering
                                              !   1 : with electeds in increasing order
                                              !   2 : with electeds in preferential order
                                              !   3 : represented clusters represented - by top portion
     Real,      Intent(out) :: fitness(:,:)   ! (3,mfe) Fitness values for possible electeds
                                              !   (1,e) :  overall fitness (declining ordere)
                                              !   (2,e) :  elected rating = average Zrate0  
                                              !            over elect(:,0,e)
                                              !   (3,e) :  mismatch penalty = Product(elc_dat(:ncl,2,e))
     Real,      Intent(out) :: elc_dat(:,:,:) ! (ind,2,mfe) Objective data values for possible electeds
                                              !   (1) :  Cluster size deviations = actual - true
                                              !   (2) :  Corresponding penalties
!  Local:
     Integer, Parameter :: ntf= 3
     Real, Parameter :: Top_fac(ntf)= (/ 0.85, 0.50, 0.15 /)
     Real    :: avg_rate(nc), por_ord(nc,ind)
     Real    :: fit_ratio(Size(fitness,2))
     Integer :: ls(np), elc(np), cls(np), key(np), fit_cut(ntf)
     Integer :: e, i, j, k, n, cl, ce, cn, k1, nv
     
     nf= 0;  elc= 0;   elect= 0;  fitness= 0;  elc_dat= 0

!    The 'ord(1:)' reordering of the candidates, based on Zrate0 
!    plus the corresponding portions 'por_ord' 
     
     nv= ord(0);  avg_rate= Zrate0(ord(1:))
     Forall(cl=1:ind) por_ord(:,cl)= portion(ord(1:),cl)

     If (pr_out > 1) then
       Call Out ("Enter 'Elect_best' from cluster set",q, ln=1)
       Call Out ("Compute viable sets of electeds for # reg clusters",ind-1, ln=1)
       Call Out ("with cluster sizes",siz)
       Call Out ("Using the reordering",ord(1:))
       Call Out ("from the cluster-averaged mean vector",avg_rate)
       Call Out (-1,"and ordered portions", por_ord)
     End if

!    The all-independents cluster set

     If (ind == 1) then ! No proportionality, so choose the 'np' top candidates
       nf= 1;  elect(:,0,1)= ord(1:np);  Call Sort (.true.,elect(:,0,1)) 

       Call Electeds_objective (elect(:,0,1), siz,avg_rate,por_ord, &
                                fitness(:,1),elc_dat(:,:,1))
       Call Pref_clust (np,portion,elect(:,:,1))

       If (pr_out > 1) then
         Call Out (-1,"The all-independents electeds",elect(:,:,1))
         Call Out ("Fitness for the all-independents elected set",fitness(:,1))
         Call Out (-1,"Corresponding data",elc_dat(:,:,1))
       End if
       Return
     End if

!    Compute the objective values for possible sets of electeds,
!    using candidates in the ord(1:) order
     
     If (Test_range) then
       Call Fit_range (nv,obj_eps2,siz,avg_rate,por_ord, nf,elect(:,0,:),fitness,elc_dat)
     Else
       Call Fit_all (obj_eps2,siz,avg_rate,por_ord, nf,elect(:,0,:),fitness,elc_dat)
     End if

     If (nf >= mfe .and. pr_out > 1) then
       Call Out ("Warning in 'Elect_best': Data base full for set",q, ln=1)
     End if
    
!    Restore the standard candidate numbering for the possible 
!    elected sets elect 
     
     Do i= 1,nf
       elect(:,0,i)= ord(elect(:,0,i));  Call Sort (.true.,elect(:,0,i))
       Call Pref_clust (np,portion, elect(:,:,i))
     End do
     
!    Identify # possible electeds sets in each fitness category

     fit_ratio= fitness(1,:nf) / fitness(1,1)

     Do k= 1,ntf
       fit_cut(k)= last_true(fit_ratio >= Top_fac(k))
     End do

     If (pr_out > 1) then
       Call Out ("In 'Elect_best' for cluster set",q, &
                "# ranked sets of electeds",nf, ln=1)
       Call Out ("# in each fitness ratio category (.85,.50,.15)",fit_cut)
     End if
    
     If (pr_out > 1) then
       Call Out ("For fitness ratios",fit_ratio(:nf))
       Do i= 1,nf
         Call Out ("For possible elected set",i,ln=1)
         Call Out (-1,"Increasing order, preferential order & clusters",elect(:,:,i))
         Call Out ("Fitness: overall = avg rating * proportion penalty",fitness(:,i))
         Call Out (-1,"cluster size deviations & penalties",elc_dat(:,:,i))
       End do
     End if

   End Subroutine Elect_best

   Subroutine Merge_electeds (ncs,Clust_elect, me, elc_set,elc_obj)
 
!    Merge the possible sets of electeds for each cluster set,
!    with the non-clustering sets first, second the 
!    elected sets listed by the top  cluster set categories, 
!    then the remaining elected sets - those listed by
!    the lower level cluster sets.

     Integer,              Intent(in) :: ncs             ! # cluster sets
     Type(Multi_listR), Intent(inout) :: Clust_elect(:)  ! (ncs) Final cluster set data for 
                                                         !   determining associated sets of electeds
                                                         ! %k = index of the best set of electeds in 
                                                         !      Elect_dat = %vl(1)
                                                         ! %l = 'ind' = # clusters, including independents
                                                         ! %m = 'nf'= # sets of possible electeds [nf <= mfe]
                                                         ! %fsx = full objective for the top set of electeds 
                                                         !        %Q0(:,0,1) = Clust(q)%sx(1)
                                                         ! %fux = full objective ratio = (q)%fsx / (1)%fsx

                                                         ! %sx(7) Full objective data, as in Elect_dat%M1 

                                                         ! %ls(nonC) Location in %Q0(:,0,:) of non-clustering sets 
                                                         ! %vl(nf) = The indices in 'Elect_dat' of the sets 
                                                         !           in %Q0(:,0,eq) as reordered on output

                                                         ! %Q0(np,0:2,nf)= Possible sets of electeds. 
                                                         !   1 = increasing order, 2 = preferential order,
                                                         !   3 = corresponding clusters
                                                         ! %M2(3,nf)= Electeds fitness(1) = average rating(2) * &
                                                         !            deviation penalty(3) 
     Integer,             Intent(out) :: me              ! # significant sets of electeds recorded, me <= mxe
     Integer,             Intent(out) :: elc_set(-1:,:)  ! (-1:np,mxe) All recorded sets of electeds, 
                                                         !    each set in increasing order, with the
                                                         !    best listing cluster set 'q' at -1
                                                         !    at position 'eq' in Clust_elect(q)%Q0,
                                                         !    with 'eq' at 0. 
     Real,                Intent(out) :: elc_obj(:)      ! (mxe) Full objectives for the elected sets
!  Local:
     Integer, Allocatable :: Elc(:), key1(:), key2(:), inv(:)
     Real    :: obj

     Integer :: mf, m0, m1, m2, mq, n0, n2, nf, q1
     Integer :: e, i, j, k, l, m, n, q, eq, e1, e2

     Call Out ("Enter 'Merge_electeds'")

     mq= cls_n(nct);  m1= 0;  m2= 0;  me= 0
     elc_set= 0;  elc_obj= 0

!    First add electeds from the better cluster sets
     
     Cluster_set_loop : Do q= 1,ncs
       nf= Clust_elect(q)%m

       Electeds_loop : Do eq= 1,nf
         obj= Clust_elect(q)%sx(2) * Clust_elect(q)%M2(1,eq)

         e= Locate_set(Clust_elect(q)%Q0(:,0,eq), elc_set(1:,:me))

         If (e <= 0 .and. me < mxe) then  ! Add current set to 'elc_set'
           me= me + 1;  e= me
           elc_set(1:,e)= Clust_elect(q)%Q0(:,0,eq)
           elc_set(-1,e)= q;  elc_set(0,e)= eq
           elc_obj(e)= obj
         Else if (e > 0) then
           If (obj > elc_obj(e)) then
             elc_set(-1,e)= q;  elc_set(0,e)= eq
             elc_obj(e)= obj
           End if
         Else
           Call Out ("Warning in 'Merge_electeds': 'Elect_dat' full")
           Exit Electeds_loop
         End if

         Clust_elect(q)%vl(eq)= e
       End do Electeds_loop

     End do Cluster_set_loop

!    Reorder the possible sets of electeds 'elc_set'
!    according to full objective data. Then modify  
!    the ordering using the 'clust_cut' categories

     Allocate (Elc(me), key1(me), key2(me), inv(me))

!    First reordering

     Call Sort (.false., elc_obj(:me), key1)
     elc_set(:,:me)= elc_set(:,key1)

     Call Inverse_map (key1, inv)

     Do q= 1,ncs
       Clust_elect(q)%vl= inv(Clust_elect(q)%vl)
     End do

!    Top 2 tiers of elected sets

     Elc= 0;  key2= 0
     
     Do q= 1,cls_n(1)
       n= Last_true(Clust_elect(q)%M2(1,:) >= cls_cut(1) * Clust_elect(q)%M2(1,1))
       Elc(Clust_elect(q)%vl(:n))= 1
     End do

     Do q= cls_n(1)+1, cls_n(2)
       n= Last_true(Clust_elect(q)%M2(1,:) >= cls_cut(2) * Clust_elect(q)%M2(1,1))
       Where(Elc(Clust_elect(q)%vl(:n)) == 0) Elc(Clust_elect(q)%vl(:n))= 2
     End do

!    Bottom tier cluster set electeds

     Call List_of_true (Elc == 1, m1,key2)
     Call List_of_true (Elc == 2, m2,key2(m1+1:))

     m2= m1 + m2;  elc_n(1)= m1;  elc_n(2)= m2
     Call List_of_true (Elc == 0, m0,key2(m2+1:))

!    Final reordering, based on top tier, then 2nd tier, then bottom

     elc_set(:,:me)= elc_set(:,key2)
     elc_obj(:me)  = elc_obj(key2)

     Call Inverse_map (key2, inv)

     Do q= 1,ncs
       Clust_elect(q)%vl= inv(Clust_elect(q)%vl)
       Clust_elect(q)%k = Clust_elect(q)%vl(1)
     End do

     If (pr_out > 1) then
       Call Out ("Elected set categories",Elc)
       Call Out ("Rordering key",key2)
       Call Out ("Full objectives for final reordering",elc_obj(:me))
       Call Out (1,"Matching cluster set # and its electeds index",elc_set(-1:0,:me))
       Call Out (-1,"The corresponding possible elected sets",elc_set(1:,:me))
     End if

   End Subroutine Merge_electeds

    
   Subroutine Organize_electeds (nc,np,ncs,me, Clust_elect, Elect_dat, dpr)

!    For each significant possible set of electeds, find the corresponding 
!    objective data for different cluster sets. Also sort the sets of electeds 
!    by top full objective.

     Integer,             Intent(in) :: nc              ! # candidates
     Integer,             Intent(in) :: np              ! # candidates to be elected
     Integer,             Intent(in) :: ncs             ! # cluster sets
                                                        !   by the best cluster sets 
     Integer,             Intent(in) :: me              ! # feasible sets of electeds
                                                        !           in increasing order
     Type(Multi_listR),   Intent(in) :: Clust_elect(:)  ! (ncs) Final cluster set data with sets of electeds
                                                        ! %k = index of the best set of electeds in 
                                                        !      Elect_dat = %vl(1)
                                                        ! %l = 'ind' = # clusters, including independents
                                                        ! %m = 'nf'= # sets of possible electeds [nf <= mfe]
                                                        ! %fsx = full objective for the top set of electeds 
                                                        !        %Q0(:,0,1) = Clust(q)%sx(1)
                                                        ! %fux = full objective ratio = (q)%fsx / (1)%fsx

                                                        ! %sx(7) =  Full objective data, as in Elect_dat%M1 

                                                        ! %ls(nonC) Location in %Q0(:,0,:) of non-clustering sets 
                                                        ! %vl(nf) = The indices in 'Elect_dat' of the sets 
                                                        !           in %Q0(:,0,eq) as reordered on output

                                                        ! %Q0(np,0:2,nf)= Possible sets of electeds. 
                                                        !   0 = increasing order, 1 = preferential order,
                                                        !   2 = corresponding clusters
                                                        ! %M2(3,nf)= Fitness objective(1)= 
                                                        !            electeds rating(2) * penalty(3) 

     Type(Multi_listR), Intent(inout) :: Elect_dat(:) ! (me) Cluster set and elected set data for
                                                      !      possible elected sets 'e' and 
                                                      !      associated cluster sets 'q', with 'e' 
                                                      !      ordered by decreasing full objective 
                                                      !      within cluster set categories
                                                      ! %k = The top cluster set 'q' for 'e' = %lt(1) [input]
                                                      ! %l = The matching 'eq' in Clust_elect(q)%Q0 = %vl(1) [input]
                                                      ! %m = nq = # recorded cluster sets that list 'e' 
                                                      ! %n = # regular clusters for best cluster set 'q'
                                                      ! %q = Non-clustering index if 'e' is a
                                                      !      non-clustering elected set, else 0

                                                      ! %ls(np)    The set of electeds in increasing order [input]
                                                      ! %L0(np,nq) The set in preferential order for each
                                                      !            of the cluster sets ordered as in %lt

                                                      ! %lt(nq) Ordering of the 'nq' cluster sets 'q' which  
                                                      !         list the set of electeds 'e' at a position 
                                                      !         'eq' in Clust_elect(q)%Q0(:,:,eq) with eq <= nf. 
                                                      !         The ordering is by full objective %M1(1,i) 
                                                      !         where %lt(i) = q.
                                                      ! %vl(nq)    The corresponding indices 'eq'
     
                                                      ! %M1(7,nq) Objective data for the clusters sets ordered %lt.
                                                      !   (1,i) = full objective = (2,i)*(3,i)
                                                      !   (2,i) = clustering objective = (4,i)*(5,i)
                                                      !   (3,i) = electeds fitness = (6,i)*(7,i)
                                                      !   (4,i) = cluster set coherent vote 
                                                      !   (5,i) = cluster set penalty
                                                      !   (6,i) = elected rating = average rating over electeds  
                                                      !   (7,i) = deviation size penalty

     Integer,            Intent(out) :: dpr(0:)       ! (0:mfe) Elected sets whose best full objective 
                                                      !         exceeds that of cluster set #1.
                                                      !         (0) = # such elected sets.
!  Local:
     Real    :: full1, tmp(ncs)
     Integer :: lse(ncs), lsq(ncs), ord(ncs)
     Integer :: e, i, j, n, q, eq, nd, nf, nq, nonC

     Call Out ("Enter 'Organize_electeds'")

!    Elected sets listed by the best 'ncs' cluster sets

     Electeds_loop : Do e= 1,me

       nq= 0;  lsq= 0;  lse= 0

       Do q= 1,ncs
         eq= First_true(Clust_elect(q)%vl == e)
         If (eq > 0) then
           nq= nq + 1;  lsq(nq)= q;  lse(nq)= eq 
         End if
       End do

       If (nq <= 0) then  ! Error condition
         Call Out ("Error in 'Organize_electeds': No cluster set for elected set",e, ln=1)
         Stop
       End if
       Elect_dat(e)%m= nq

       Allocate(Elect_dat(e)%lt(nq), Elect_dat(e)%vl(nq), &
                Elect_dat(e)%l0(np,nq), Elect_dat(e)%M1(7,nq))
       tmp= 0;  ord= 0

       Do i= 1,nq
         q= lsq(i);  eq= lse(i)
         tmp(i)= Clust_elect(q)%sx(2) * Clust_elect(q)%M2(1,eq) ! full objective
       End do

       Call Sort (.false., tmp(:nq), ord(:nq))

       Elect_dat(e)%lt= lsq(ord(:nq))
       Elect_dat(e)%vl= lse(ord(:nq))

       Do i= 1,nq
         q= Elect_dat(e)%lt(i);  eq= Elect_dat(e)%vl(i)

         Elect_dat(e)%L0(:,i)  = Clust_elect(q)%Q0(:,1,eq)

         Elect_dat(e)%M1(1,i)  = tmp(i)  
         Elect_dat(e)%M1(2,i)  = Clust_elect(q)%sx(2)  
         Elect_dat(e)%M1(3,i)  = Clust_elect(q)%M2(1,eq)
         Elect_dat(e)%M1(4:5,i)= Clust_elect(q)%sx(4:5)  
         Elect_dat(e)%M1(6:7,i)= Clust_elect(q)%M2(2:3,eq) 
       End do
     End do Electeds_loop

!    List those elected sets 'e' whose best full objective 
!    exceeds the #1 full objective for the top cluster set

     nd= 0;  dpr= 0;  full1= Clust_elect(1)%fsx
     Do e= 1,me
       If (Elect_dat(e)%M1(1,1) > full1) then
         nd= nd + 1;  dpr(nd)= e;  if (nd == Ubound(dpr,1)) Exit
       End if
     End do
     dpr(0)= nd

     If (pr_out > 1) then
       If (nd > 0) Call Out ("Elected sets whose best full objective > cluster set #1",dpr(1:nd))
       Do e= 1,me
         nq= Elect_dat(e)%m  
         Call Out ("For elected set",e, "# associated cluster sets",nq, ln=1)
         Call Out ("The associated cluster sets", Elect_dat(e)%lt)
         Call Out ("The elected set location in the cluster set list", Elect_dat(e)%vl)
         Call Out (-1,"Objective data for each cluster set",Elect_dat(e)%M1)
       End do
     End if
   End Subroutine Organize_electeds

                                        
   Subroutine Analyze_elected (nc,vwt, bwt,cand, Elect_dat)

!    Analyze the success of the clustering.
  
!    "Overall satisfaction measure" of voter success:  Each voter regards the set of electeds
!    as successful to the extent that his or her top candidates are in the set of electeds.
!    the top candidate gets a weight of vwt(1), the voter's second ranked candidate
!    gets a weight of vwt(2), etc. The weights must sum to 1 and the number of top
!    candidates used is 'Mvt'. Averaging this measure of success for each voter, over all
!    voters, gives the overall satisfaction. Example: vwt= (/ 4/7., 2/7., 1/7. /)

     Integer,             Intent(in) :: nc            ! # candidates
     Real,                Intent(in) :: vwt(:)        ! (Mvt) Weights corresponding to a voter's satisfaction rankings, summing to 1
     Real,                Intent(in) :: bwt(:)        ! (nb) Ballot weights, summing to 'np'
     Integer,             Intent(in) :: cand(0:,:)    ! (0:,nb) 0 = # candidates 'n' ranked or rated, 
                                                      ! 1:n = those candidates ordered by ranking or rating
     Type(Multi_listR),Intent(inout) :: Elect_dat(:)  ! (me) Data associated with the possible sets of electeds, 
                                                      !      ordered by decreasing full objective.
                                                      ! %k = The top cluster set for the set of electeds = %lt(1)
                                                      ! %l   = nq = # top cluster sets that list this elected set
                                                      ! %fsx = overall satisfaction
                                                      ! %fux = full objective ratio = (e)%M1(1,1) / (1)%M1(1,1)

                                                      ! %ls(np) = the set of electeds in increasing order

                                                      ! %lt(nq) Ordering of the 'nq' cluster sets 'q' which list the 
                                                      !         set of electeds 'e' at a position 'eq' in 
                                                      !         Clust_elect(q)%Q0(:,eq,:) with eq <= nf. The ordering 
                                                      !         is by decreasing full objective Clust_elect(q)%M2(0,eq).
                                                      ! %vl(nq)    The corresponding indices 'eq'
                                                      ! %L0(np,nq) This elected set in preferential order for each of the
                                                      !            the cluster sets ordered as in %lt
     
                                                      ! %M0(0:Mvt,2): Measures of success - over ballots
                                                      !   (0,1) : mean # of top 3 elected in order (but not the next lower)
                                                      !   (0,2) : mean # of top 3 elected regardless of order
                                                      !   (n,1) : fraction of voters with top n or more elected in order
                                                      !   (n,2) : fraction of voters with top n or more elected regardless of order

!  Local:
     Logical :: Elc(nc), Nelc(nc)
     Real    :: wt, tot_wt
     Integer :: b, e, l1, l2, n, me, Mvt, nb
    
     If (pr_out > 1)  Call Out ("Enter 'Analyze_elected'")

     me= Size(Elect_dat);  Mvt= Size(vwt)
     nb= Size(cand,2);  tot_wt= Sum(bwt)

     Do e= 1,me
       Elect_dat(e)%fsx= 0;  Elect_dat(e)%M0= 0
       Elc= .false.;  Elc(Elect_dat(e)%ls)= .true.;  Nelc= .not.Elc
      
       Do b= 1,nb
         n = Min(cand(0,b),Mvt)
         l1= First_true(1, Nelc(cand(1:n,b))) - 1  ! = # top consecutive candidates elected
         l2= Count(Elc(cand(1:n,b)))               ! = # top candidates elected, without regard to order
         
         If (l1 > 0) then
           Elect_dat(e)%M0(0,1)   = Elect_dat(e)%M0(0,1) + l1 * bwt(b)
           Elect_dat(e)%M0(1:l1,1)= Elect_dat(e)%M0(1:l1,1) + bwt(b)
         End if
        
         If (l2 > 0) then
           Elect_dat(e)%M0(0,2)   = Elect_dat(e)%M0(0,2) + l2 * bwt(b)
           Elect_dat(e)%M0(1:l2,2)= Elect_dat(e)%M0(1:l2,2) + bwt(b)

           wt= Sum(vwt(:n), Elc(cand(1:n,b)))
           Elect_dat(e)%fsx= Elect_dat(e)%fsx + wt * bwt(b)
         End if
       End do
     End do
    
     Elect_dat%fsx= Elect_dat%fsx / tot_wt
     Do e= 1,me
       Elect_dat(e)%M0= Elect_dat(e)%M0 / tot_wt
     End do
   End Subroutine Analyze_elected


   Subroutine Add_to_fit (fit0,obj,dat,cur_sub, nf,best_sub, best_fit,best_dat,loc)

!    Test for a newly evaluated set of electeds 'cur_sub' and store it in 'best_sub'
!    along with its fitness data 'obj' and 'dat' in 'best_fit' in decreasing order 
!    of fitness. The higher the fitness the better.

     Real,         Intent(in) :: fit0             ! Lower bound for fitness feasibility
     Real,         Intent(in) :: obj(:)           ! (3)  Objective data for current set
     Real,         Intent(in) :: dat(:,:)         ! (ind,2) Cluster deviations and penalties
                                                  !         for the current set
     Integer,      Intent(in) :: cur_sub(:)       ! (np) Current set, in increasing order, 
                                                  !      to be added if possible
       
     Integer,   Intent(inout) :: nf               ! # sets recorded to date, nf <= mx
     Integer,   Intent(inout) :: best_sub(:,:)    ! (np,mx) Records good subsets, ordered
                                                  !   by decreasing fitness   
     Real,      Intent(inout) :: best_fit(:,:)    ! (3,mx) Decreasing fitness values of
                                                  !  subsets recorded in best_sub
     Real,      Intent(inout) :: best_dat(:,:,:)  ! (ind,2,mx) Fitness data values of 
                                                  !   subsets recorded in best_sub
     Integer,     Intent(out) :: loc              ! Location in 'best_sub' of the new 
                                                  !   set of electeds 'cur_sub',
                                                  !   or minus this location if already
                                                  !   present in 'best_sub',
                                                  !   or mx+1 if not recorded
 ! Local:
     Real    :: cur_fit
     Integer :: j, k, j1, k1, l1, l2, mx, n1

     mx= Size(best_fit,2);  n1= Min(nf + 1,mx)

     cur_fit= obj(1);  loc= mx + 1
     If (cur_fit <= fit0) then
       If (pr_out > 1)  Call Out("Bad set not recorded");  Return
     End if
 
     If (nf <= 0) then  ! Initial call
       loc= 1
     Else  ! Subsequent call
       loc= Minloc(Abs(cur_fit - best_fit(1,:nf)),1)
       
       If (Abs(cur_fit - best_fit(1,loc)) < 0.01*fit0) then
         If (All(cur_sub == best_sub(:,loc))) then
           loc= -loc;  If (pr_out > 1)  Call Out("Repeated set not recorded")
           Return
         End if
       End if
       
       If (cur_fit <= best_fit(1,loc)) then
         l1= loc + 1;  k= n1
         
         If (k >= l1) then  ! Push right and prepare to insert at loc+1
           l2= loc + 2;  k1= k - 1
           If (k >= l2) then  
             best_fit(:,l2:k)  = best_fit(:,l1:k1)
             best_dat(:,:,l2:k)= best_dat(:,:,l1:k1)
             best_sub(:,l2:k)  = best_sub(:,l1:k1)
           End if 
           loc= l1;  ! Insert cur_sub at loc+1
         Else
           loc= mx+1;  If (pr_out > 1) Call Out("Data full: set not recorded")
           Return
         End if

       Else ! For cur_fit > best_fit(1,loc)
         k= n1
         
         If (k >= loc) then  ! Push right and prepare to insert at loc
           l1= loc + 1;  k1= k - 1  
           If (k >= l1) then  
             best_fit(:,l1:k)  = best_fit(:,loc:k1)
             best_dat(:,:,l1:k)= best_dat(:,:,loc:k1)
             best_sub(:,l1:k)  = best_sub(:,loc:k1)
           End if
         End if
       End if
     End if
     
     nf= n1;  best_fit(:,loc)= obj;  best_dat(:,:,loc)= dat;  best_sub(:,loc)= cur_sub

     If (pr_out > 1.5) then
       Call Out ("New set at position in the ordered search date",loc)
       Call Out ("The new set",cur_sub)
       Call Out ("Its fitness values",obj)
       Call Out ("Database of overall fit values",best_fit(1,:nf))
     End if
   End Subroutine Add_to_fit

   Subroutine Add_to_fit1 (fit0,obj,cur_sub, nf,best_sub, best_fit,loc)

!    Test for a newly evaluated set of electeds 'cur_sub' and store it in 'best_sub'
!    along with its fitness data 'obj' and 'dat' in 'best_fit' in decreasing order 
!    of fitness. The higher the fitness the better.

     Real,         Intent(in) :: fit0             ! Lower bound for fitness feasibility
     Real,         Intent(in) :: obj(:)           ! (3)  Objective data for current set
     Integer,      Intent(in) :: cur_sub(:)       ! (np) Current set, in increasing order, 
                                                  !      to be added if possible
       
     Integer,   Intent(inout) :: nf               ! # sets recorded to date, nf <= mx
     Integer,   Intent(inout) :: best_sub(:,:)    ! (np,mx) Records good subsets, ordered
                                                  !   by decreasing fitness   
     Real,      Intent(inout) :: best_fit(:)      ! (mx) Decreasing fitness values of
                                                  !  subsets recorded in best_sub
     Integer,     Intent(out) :: loc              ! Location in 'best_sub' of the new 
                                                  !   set of electeds 'cur_sub',
                                                  !   or minus this location if already
                                                  !   present in 'best_sub',
                                                  !   or mx+1 if not recorded
 ! Local:
     Real    :: cur_fit
     Integer :: j, k, j1, k1, l1, l2, mx, n1

     mx= Size(best_fit);  n1= Min(nf + 1,mx)

     cur_fit= obj(1);  loc= mx + 1
     If (cur_fit <= fit0) then
       If (pr_out > 1)  Call Out("Bad set not recorded");  Return
     End if
 
     If (nf <= 0) then  ! Initial call
       loc= 1
     Else  ! Subsequent call
       loc= Minloc(Abs(cur_fit - best_fit(:nf)),1)
       
       If (Abs(cur_fit - best_fit(loc)) < 0.01*fit0) then
         If (All(cur_sub == best_sub(:,loc))) then
           loc= -loc;  If (pr_out > 1)  Call Out("Repeated set not recorded")
           Return
         End if
       End if
       
       If (cur_fit <= best_fit(loc)) then
         l1= loc + 1;  k= n1
         
         If (k >= l1) then  ! Push right and prepare to insert at loc+1
           l2= loc + 2;  k1= k - 1
           If (k >= l2) then  
             best_fit(l2:k)= best_fit(l1:k1)
             best_sub(:,l2:k)= best_sub(:,l1:k1)
           End if 
           loc= l1;  ! Insert cur_sub at loc+1
         Else
           loc= mx+1;  If (pr_out > 1) Call Out("Data full: set not recorded")
           Return
         End if

       Else ! For cur_fit > best_fit(loc)
         k= n1
         
         If (k >= loc) then  ! Push right and prepare to insert at loc
           l1= loc + 1;  k1= k - 1  
           If (k >= l1) then  
             best_fit(l1:k)= best_fit(loc:k1)
             best_sub(:,l1:k)= best_sub(:,loc:k1)
           End if
         End if
       End if
     End if
     
     nf= n1;  best_fit(loc)= obj(1);  best_sub(:,loc)= cur_sub

     If (pr_out > 1.5) then
       Call Out ("New set at position in the ordered search date",loc)
       Call Out ("The new set",cur_sub)
       Call Out ("Its fitness values",obj)
       Call Out ("Database of overall fit values",best_fit(:nf))
     End if
   End Subroutine Add_to_fit1

   Subroutine Add_non_clust (q,np,ind,mx,nonC,Non_clust, siz,Zrate0,portion, &
                             nf,elect,fitness,elc_dat, non_loc) 

!    Add on the non-clustering sets 'Non_clust', as needed, to the electeds data
!    'elect', 'fitness', 'elc_dat' for the current cluster set 'q'.

     Integer,    Intent(in) :: q               ! Current cluster set
     Integer,    Intent(in) :: np              ! # candidates to be elected
     Integer,    Intent(in) :: ind             ! # clusters in this set
     Integer,    Intent(in) :: mx              ! Max # possible sets of electeds
     Integer,    Intent(in) :: nonC            ! # Non-clustering sets of electeds
     Integer,    Intent(in) :: Non_clust(:,:)  ! (np,nonC) Non-clustering sets of electeds

     Real,       Intent(in) :: siz(:)          ! (ind) Cluster sizes (= weights)
     Real,       Intent(in) :: Zrate0(:)       ! (nc) Cluster-averaged mean vectors with zeroing
     Real,       Intent(in) :: portion(:,:)    ! (nc,ind) Portion of each cluster for each candidate
     
     Integer, Intent(inout) :: nf              ! # matching elected sets found
     Integer, Intent(inout) :: elect(:,0:,:)   ! (np,0:2,mx) Corresponding orderings of the electeds 
                                               !     (:,0,i) = electeds in increasing order 
                                               !     (:,1,i) = electeds in preferential order 
                                               !     (:,2,i) = the clusters that represent (:,1,i) 
     Real,    Intent(inout) :: fitness(:,:)    ! (3,mx) Fitness values for possible electeds
                                               !   (1,e) :  overall fitness (declining ordere)
                                               !   (2,e) :  elected rating = average Zrate0  
                                               !            over elect(:,0,e)
                                               !   (3,e) :  mismatch penalty = Product(elc_dat(:ncl,2,e))
     Real,    Intent(inout) :: elc_dat(:,:,:)  ! (ind,2,mx) Objective data values for possible electeds
                                               !   (1) :  Cluster size deviations = actual - true
                                               !   (2) :  Corresponding penalties
     Integer,   Intent(out) :: Non_loc(:)      ! (nonC) Location of each non-clustering set in elect
!  Local
     Real    :: obj(nonC), fit(3,nonC), dat(ind,2,nonC), tmp(mx)
     Integer :: loc(nonC), key(mx), Non(mx), inv(mx), miss(nonC), elc(np,0:2,nonC)
     Logical :: ReOrd
     Integer :: e, i, j, k, l, n, e1, n1, nl, ns
     
     ns= 0;  Non_loc= 0;  obj= 0;  fit= 0;  dat= 0;  key= 0;  Non= 0;  inv= 0;  miss= 0

!    First: Locate the non-clustering elected sets 'Non_clust' 
!    in 'elect'. Identify in 'Non_loc', plus any missing in 'miss'

     Do j= 1,nonC
       e= Locate_set(Non_clust(:,j), elect(:,0,:nf))

       If (e >= 1 .and. e <= nf) then
         Non(e)= j;  Non_loc(j)= e
       Else
         ns= ns + 1;  miss(ns)= j
       End if
     End do

     If (ns <= 0) Return  ! None missing

!    Reduce 'nf' if necessary to make room for the rest of the non-clustering
!    elected sets. Do this by removing the least fit regular sets.

     k= ns - (mx - nf)

     Do j= 1,k
       e= Last_true(Non(:nf) == 0);  n1= nf - 1 ! Remove the last clustering elected
                                                ! to make space for a missing 
                                                ! non-clustering elected set
       If (e < nf) then
         e1= e + 1;  Non(e:n1)= Non(e1:nf)
         elect(:,:,e:n1)  = elect(:,:,e1:nf)

         fitness(:,e:n1)  = fitness(:,e1:nf)
         elc_dat(:,:,e:n1)= elc_dat(:,:,e1:nf)
       End if

       Non(nf)= 0;  elect(:,:,nf)= 0;  fitness(:,nf)= 0;  elc_dat(:,:,nf)= 0
       nf= n1
     End do

!    Evaluate the missing non-clustering sets

     Do i= 1,ns
       e= nf + i;  elect(:,0,e)= Non_clust(:,miss(i))

       Call Electeds_objective (elect(:,0,e), siz,Zrate0,portion, &
                                fitness(:,e),elc_dat(:,:,e))

       Call Pref_clust (np,portion, elect(:,:,e))
       Non(e)= miss(i)
     End do

     nf= e;  Call Inverse_map (Non,Non_loc)

!    Reorder the elected sets according to decreasinng fitness

     tmp= fitness(1,:nf);  Call Sort (.false., tmp(:nf), key(:nf), ReOrd)

     If (ReOrd) then
       fitness(:,:nf)  = fitness(:,key(:nf))
       elect(:,:,:nf)  = elect(:,:,key(:nf))
       elc_dat(:,:,:nf)= elc_dat(:,:,key(:nf))

       Call Inverse_map (key(:nf),inv(:nf))
       Non_loc= inv(Non_loc)
     End if

     If (pr_out > 1) then
       Call Out ("Add_non_clust': from cluster set",q, ln=1)
       Call Out ("# non-clusterting sets missing",ns)
       Call Out ("Non_clustering electeds locations",Non_loc)
       Call Out ("and fitnesses",fitness(1,Non_loc))
       Call Out ("Full list of fitnesses",fitness(1,:nf))
       Call Out ("ReOrd key",key(:nf))

       Call Out (-1,"Non_clustering elected sets",elect(:,0,Non_loc))
       Call Out (-1,"vs prior listing",Non_clust)
     End if

   End Subroutine Add_non_clust

   
   Subroutine Elect_fit (iq,np,nc,ind, Memb,Clust_vec, &
                         elc2, elect1, Por,Corr,fit) 

!    Determine the sets of 'np' candidates which are the 
!    best match ("optimal") for the cluster set characterized
!    by Zrate0, siz, and OrdZ. These will be up to the top 'mf' 
!    such sets as measured by fitness.

!    These have the highest overall cluster rankings, 
!    where these rankings are weighted by portions 
!    which reflect candidate memberships in the clusters.
     
     Integer,    Intent(in) :: iq             ! Cluster set index
     Integer,    Intent(in) :: np             ! # candidates to be elected
     Integer,    Intent(in) :: nc             ! # candidates
     Integer,    Intent(in) :: ind            ! # clusters, including independents

     Type(Multi_listD), Intent(in) :: Memb(:)    ! (ns) Final slate cluster data
     Real,       Intent(in) :: Clust_vec(0:,0:)  ! (0:nc,0:ind) Mean rating vector for each cluster
                                                 !   with cluster size at (0,:) and the cluster
                                                 !   averaged mean vector at (:,0)
     
     Integer,    Intent(in) :: elc2(:,:)      ! (np,2) Two possible sets of electeds
                                              !   0 : with electeds in increasing order
                                              !   1 : with electeds in preferential order
                                              !   2 : clusters represented - by top portion
     Integer,   Intent(out) :: elect1(:,0:)   ! (np,0:2) Computed optimal set of electeds 
                                              !   0 : with electeds in increasing order
                                              !   1 : with electeds in preferential order
                                              !   2 : clusters represented - by top portion

     Real,      Intent(out) :: Por(:,0:)      ! (nc,0:ind) Cluster portions (:,1:) with averaged 
                                              !            noise zeroed ratings Zrate0 at (:,0)
     Real,      Intent(out) :: Corr(:,:)      ! (ncl,ncl) Correlation matrix

     Real,      Intent(out) :: fit(:,:)       ! (7,3) Cluster objective and electeds fitness data
                                              !       for the 2 input sets of electeds and the 
                                              !       output optimal set
                                              !   1 :  full objective
                                              !   2 :  clustering objective
                                              !   3 :  overall fitness
                                              !   4 :  coherent vote
                                              !   5 :  cluster penalty
                                              !   6 :  average elected rating
                                              !   7 :  cluster size deviation penalty
!  Local:
     Integer, Parameter :: mf= 1
     Real    :: Zrt1(nc)          ! Scaled Zrate0: cluster-averaged noise-zeroed  
                                  !      mean ratings, scaled by 1+eps - portion(:,ind)

     Real    :: Zrate0(nc)        !  Cluster averaged mean vector of the 
                                  !  noise-zeroed rating matrix
     Real    :: por2(nc,ind,2)    ! (1) Noise zeroed ratins, and (2) portions of the clusters 
                                  !     that represent each candidate

     Integer :: Zrt1_ord(-2:nc)   ! Decreasing ordering of Zrt1 at (1:nc) 
                                  !   for Zrt1 = Zrate0 scaled by [1+eps - portion(:,ind)],
                                  !   with # top, significant, and viable at (-2:0)
     Integer :: elect(np,0:2,mf)  ! Corresponding for possible electeds,
                                  !   using the overall fitness ordering
                                  !   0 : with electeds in increasing order
                                  !   1 : with electeds in preferential order
                                  !   2 : represented clusters represented - by top portion
     Real    :: fitness(3,mf)     ! Fitness values for possible electeds
                                  !   1 :  overall fitness (declining order)
                                  !   2 :  elected rating = average Zrate0 over elect(:,0,e)
                                  !   3 :  mismatch penalty = Product(elc_dat(:ncl,2,e))
     Real    :: elc_dat(ind,2,mf) ! Proportionality data for possible electeds
                                  !   1 :  overall fitness (declining order)
                                  !   2 :  elected rating = average Zrate0 over elect(:,0,e)
                                  !   3 :  mismatch penalty = Product(elc_dat(:ncl,2,e))

     Real    :: siz(ind)
     Real    :: dat(3), pen(3), obj(3), mx_cor(ind-1)
     Integer :: i, k, nf, ns, ncl

     ncl= ind - 1;  ns= Size(Memb);  siz= Clust_vec(0,1:)
     elect1= 0;  Por= 0;  Corr= 0;  fit= 0

     Call Portions (np,nc,Noise_por, siz,Clust_vec(1:,1:), Zrate0,por2)
     Por(:,0)= Zrate0;  Por(:,1:)= por2(:,:,2)
     Call Zrf_M (nc, Zrate0,por2(:,ind,2), Zrt1_ord, Zrt1)

     Call Correlate_clusters (nc,ncl,Noise_cor,Clust_vec(1:,1:ncl), mx_cor,Corr)

     Call Elect_best (iq,nc,np,ind,Zrt1_ord, siz,Zrate0, &
                      por2(:,:,2), nf,elect,fitness,elc_dat) 
     elect1= elect(:,:,1)

     dat(1)= Minval(siz(:ncl));  dat(2)= Maxval(mx_cor);  dat(3)= siz(ind) / np
     Call Objective_values0 (np,Parm_obj(5:10), siz(:ncl),dat, pen,obj)

     fit(2,:)= obj(1);  fit(4,:)= obj(2);  fit(5,:)= obj(3)
     fit(3,3)= fitness(1,1);  fit(6:7,3)= fitness(2:3,1)

     Do i= 1,2
       Call Electeds_objective (elc2(:,i), siz,Zrate0,por2(:,:,2), &
                                fitness(:,1),elc_dat(:,:,1))
       fit(3,i)= fitness(1,1);  fit(6:7,i)= fitness(2:3,1)
     End do

     fit(1,:)= fit(2,:) * fit(3,:)

     If (pr_out > 1) then
       Call Out ("Data computed in 'Elect_fit' for cluster set",iq, ln=1)
       Call Out ("Cluster averaged, noise zeroed, mean vector",Zrate0)
       Call Out ("The portion modified Zrate0",Zrt1)
       Call Out ("It's candidate ordering",Zrt1_ord(1:))
       Call Out ("It's # top, significant, & viable candidates",Zrt1_ord(:0))
       Call Out (-1,"Portions",Por(:,1:))
       Call Out (1,"Cluster set correlation matrix",Corr)
       Call Out ("Best fit elected set",elect1(:,0))
       Call Out (-1,"Cluster and best electeds fitness data for 3 elected sets",fit)
     End if

   End Subroutine Elect_fit


   Subroutine Possible_electeds (np,nc,ind,ne,full, domain, Clust_set, nf,Clust_elc)

!    Compute the possible sets of elected sets of candidates and their elected fitness
!    for the specified cluster set. Also compute the dominance factors between
!    the different sets of possible electeds: 

!    The idea is that when a set 'l' has lower fitness than a set 'k', then 'k' 
!    is effectively completely dominant over 'l' when rt = 1 - fitness(l)/fitness(k) >= 'full',
!    a parameter such as 0.20. That is, when comparing whether or not 'l' is an overall
!    winner over 'k', it's important to look at 2 factors. 

!    One factor is the number of domains where one set is a clear winner over the other, 
!    and the second factor is to distinguish how competitive they are on the remaining domains, 
!    paying much more attention to small differences. This is accomplished by a cosine rise 
!    function to compute a dominance factor between 0 (no preference) and 1 (full dominance).
    
     Integer,              Intent(in) :: np   ! # candidate to be elected
     Integer,              Intent(in) :: nc   ! # candidates
     Integer,              Intent(in) :: ind  ! # clusters, including independents
     Integer,              Intent(in) :: ne   ! # subsets of 1...nc of size 'np'
     Real,                 Intent(in) :: full ! Limit parameter for the rising cosine dominance function

     Integer,              Intent(in) :: domain(:) ! (nc) Original candidates

     Type(Multi_listR),    Intent(in) :: Clust_set ! Top cluster set
                                                   ! %sx(ind): Cluster sizes
                                                   ! %T0(nc,0:ind,3) Cluster point data
                                                   !    (:,1:,1) = Cluster mean vectors
                                                   !    (:,1:,2) = Cluster mean vectors with noise zeroing 
                                                   !    (:,1:,3) = Cluster portions for each candidate
                                                   !    (:,0,1)  = Scaled Zrate0 (for ordering the %L1(1:,0))
                                                   !    (:,0,2)  = Zrate0 = cluster averaged noise zeroed mean ratings

     Integer,             Intent(out) :: nf        ! # possible sets of electeds recorded
     Type(Multi_listR), Intent(inout) :: Clust_elc ! Cluster set data for a given domain
                                                   ! %n= nc = # candidates in the domain
                                                   ! %l= ind = # clusters, including independents
                                                   ! %m= nf  = # possible elected sets

                                                   ! %fsx      = Top fitness value
                                                   ! %rx(nf)   = Fitness ratios: fitness(i) / fitness(1)
                                                   ! %L1(np,nf)= Corresponding sets of electeds (original candidates)
                                                   ! %M1(nf,nf)= Dominance factors between the %L1 sets.
                                                   !    0 = no dominance. Full dominance of set 'k' over 'l' when
                                                   !    %M1(l,k) = 1, with k < l, i.e., fitness(l) < fitness(k)
                                                 
                                                   ! %ls(nc) = Domain candidates  ***
                                                   ! %sx(nc) = Domain average candidate scores  ***
                                                   ! %vl(nc) = Domain reordering by average candidate scores  ***
! Local:
     Real    :: rt, fitness(ne)
     Integer :: i, k, l, sets(np,ne)

!    Generate 'nf' feasible elected sets 'sets', ordered by their 'fitness'

     Call Elect_ord (nc,np,ne, Clust_set%sx, Clust_set%T0(:,0,2), &
                     Clust_set%T0(:,1:,3), nf,sets, fitness)

     Clust_elc%n= nc;  Clust_elc%l= ind
     Clust_elc%m= nf;  Clust_elc%fsx= fitness(1)
     Allocate(Clust_elc%rx(nf), Clust_elc%L1(np,nf), Clust_elc%M1(nf,nf))

     Clust_elc%rx= fitness(:nf) / fitness(1)
     Forall(i=1:np) Clust_elc%L1(i,:)= domain(sets(i,:nf))

     Clust_elc%M1= 0
     Do k= 1,nf-1
       Do l= k+1,nf
         rt= 1 - fitness(l) / fitness(k)               ! 'rt' >= 0 since fitness(l) <= fitness(k)
         Clust_elc%M1(l,k)= Cos_rise(-full,full, rt)   ! From 0 at rt = 0 to 1 for rt >= 'full',
                                                       ! rising rapidly at first, then leveling off,
                                                       ! as full dominance is approached. 
       End do
     End do

   End Subroutine Possible_electeds


   Subroutine Elect_ord (nc,np,ne, siz,Zrate0, portion, nf,elect,fitness) 

     Integer,    Intent(in) :: nc             ! # candidates
     Integer,    Intent(in) :: np             ! # candidates to be elected
     Integer,    Intent(in) :: ne             ! max # elected sets to be recorded

     Real,       Intent(in) :: siz(:)         ! (ind) Cluster sizes
     Real,       Intent(in) :: Zrate0(:)      ! (nc) Cluster-averaged mean vectors with zeroing
     Real,       Intent(in) :: portion(:,:)   ! (nc,ind) Portion of each cluster for each candidate
     
     Integer,   Intent(out) :: nf             ! # top elected sets recorded
     Integer,   Intent(out) :: elect(:,:)     ! (np,ne) Corresponding sets of ordered electeds, 
                                              !         using the overall fitness ordering
     Real,      Intent(out) :: fitness(:)     ! (ne) Fitness values for possible electeds
!  Local:
     Integer :: i, k
     
     nf= 0;  elect= 0;  fitness= 0

     If (pr_out > 1) then
       Call Out ("Compute viable sets of electeds with cluster sizes",siz)
       Call Out ("from the cluster-averaged mean vector",Zrate0)
       Call Out (-1,"and ordered portions", portion)
     End if

     Call Fit_all1 (obj_eps2,siz,Zrate0,portion, nf,elect,fitness)

!    Identify # possible electeds sets in each fitness category

     If (pr_out > 1) then
       Call Out ("In 'Elect_ord' # ranked sets of electeds",nf, ln=1)

       If (pr_out > 1) then
         Call Out ("From fitness values",fitness(:nf))
         Call Out (-1,"and possible elected sets",elect(:,:nf))
       End if
     End if

   End Subroutine Elect_ord

   Subroutine Fit_all (fit0,siz,avg_rate,portion, nf,best_elect,best_fit,best_dat)
  
!    Find all possible sets of electeds and order them by fitness   
    
     Real,       Intent(in) :: fit0            ! Lower bound for fitness feasibility
     Real,       Intent(in) :: siz(:)          ! (ind) Cluster sizes = target fractional # elected

     Real,       Intent(in) :: avg_rate(:)     ! (nc)  Candidate average rating, in decreasing order
                                               !       for equal portion(:,ind)
     Real,       Intent(in) :: portion(:,:)    ! (nc,ind)  Candidate portions 

     Integer,   Intent(out) :: nf              ! updated # sets recorded
     Integer, Intent(inout) :: best_elect(:,:) ! (np,ne) Records the top 'nf'  (up to mx) best sets of electeds,
                                               !         with the sets ordered by decreasing fitness objective
     Real,    Intent(inout) :: best_fit(:,:)   ! (3,ne)  Fitness objective values of the best_elect sets
     Real,    Intent(inout) :: best_dat(:,:,:) ! (ind,2,ne) Actual - true cluster sizes(1) and corresponding penalties(2)
!  Local:
     Integer, Allocatable :: sub(:,:)
     Real    :: obj(3), dat(Size(siz),2)
     Integer :: i, j, k, n, nc, np, nsub, loc

     nc= Size(avg_rate);  np= Size(best_elect,1)

     nsub= N_sub_sum(nc,np) - 1;  Allocate (sub(0:np,nsub))
     k= 0;  sub= 0

     Call List_subsets (1,nc,np,0, k,sub)  ! k = nsub

     nf= 0
     Do i= 1,k
       n= sub(0,i)  ! # candidates in subset
       If (n /= np) Cycle

       Call Electeds_objective (sub(1:,i), siz,avg_rate,portion, obj,dat) 

       Call Add_to_fit (fit0,obj,dat,sub(1:,i), nf,best_elect, &
                        best_fit,best_dat,loc)
     End do
    
   End Subroutine Fit_all

   Subroutine Fit_all1 (fit0,siz,avg_rate,portion, nf,best_elect,best_fit)
  
!    Find all possible sets of electeds and order them by fitness   
    
     Real,       Intent(in) :: fit0            ! Lower bound for fitness feasibility
     Real,       Intent(in) :: siz(:)          ! (ind) Cluster sizes = target fractional # elected

     Real,       Intent(in) :: avg_rate(:)     ! (nc)  Candidate average rating, in decreasing order
                                               !       for equal portion(:,ind)
     Real,       Intent(in) :: portion(:,:)    ! (nc,ind)  Candidate portions 

     Integer,   Intent(out) :: nf              ! updated # sets recorded
     Integer, Intent(inout) :: best_elect(:,:) ! (np,ne) Records the top 'nf'  (up to mx) best sets of electeds,
                                               !          with the sets ordered by decreasing fitness objective
     Real,    Intent(inout) :: best_fit(:)     ! (ne)  Fitness objective values of the best_elect sets
!  Local:
     Integer, Allocatable :: sub(:,:)
     Real    :: obj(3), dat(Size(siz),2)
     Integer :: i, j, k, n, nc, np, nsub, loc

     nc= Size(avg_rate);  np= Size(best_elect,1)

     nsub= N_sub_sum(nc,np) - 1;  Allocate (sub(0:np,nsub))
     k= 0;  sub= 0

     Call List_subsets (1,nc,np,0, k,sub)  ! k = nsub

     nf= 0
     Do i= 1,k
       n= sub(0,i)  ! # candidates in subset
       If (n /= np) Cycle

       Call Electeds_objective (sub(1:,i), siz,avg_rate,portion, obj,dat) 

       Call Add_to_fit1 (fit0,obj,sub(1:,i), nf,best_elect, best_fit,loc)
     End do
    
   End Subroutine Fit_all1


   Subroutine Fit_range (nv,fit0,siz,avg_rate,portion, nf,best_elect, &
                         best_fit,best_dat)
  
!    Find a range possible sets of electeds and order them by fitness.
!    These sets are predetermined beyond candidate 'nv', to be the
!    the candidtes nv+1,...nv+j if there are np - j candidates 
!    elected through 'nv'.
    
     Integer,    Intent(in) :: nv              ! Last variable candidate
     Real,       Intent(in) :: fit0            ! Lower bound for fitness feasibility
     Real,       Intent(in) :: siz(:)          ! (ind) Cluster sizes = target fractional # elected

     Real,       Intent(in) :: avg_rate(:)     ! (nc)  Candidate average rating, in decreasing order
                                               !       for equal portion(:,ind)
     Real,       Intent(in) :: portion(:,:)    ! (nc,ind)  Candidate portions 

     Integer,   Intent(out) :: nf              ! # sets recorded
     Integer, Intent(inout) :: best_elect(:,:) ! (np,ne) Records the top 'nf'  (up to mx) best sets of electeds,
                                               !          with the sets ordered by decreasing fitness objective
     Real,    Intent(inout) :: best_fit(:,:)   ! (3,ne)  Fitness objective values of the best_elect sets
     Real,    Intent(inout) :: best_dat(:,:,:) ! (ind,2,ne) Actual - true cluster sizes(1) and corresponding penalties(2)
!  Local:
     Integer, Allocatable :: sub(:,:)
     Real    :: obj(3), dat(Size(siz),2)
     Integer :: i, j, k, n, n1, n2, nc, np, nsub, loc

     nc= Size(avg_rate);  np= Size(best_elect,1);  nf= 0
     If (nv < 1) Return

     n1= Max(np - (nc-nv),1)  ! Min # candidates <= nv
     n2= Min(nv,np)           ! Max # candidates <= nv

     nsub= N_sub_sum(nv,n2) - 1;  Allocate (sub(0:np,nsub))
     k= 0;  sub= 0

     Call List_subsets (1,nv,n2,0, k,sub)  ! k = nsub

     Do i= 1,nsub
       n= sub(0,i)  ! # candidates  <= nv in subset
       If (n < n1 .or. n > n2) Cycle

       Forall(j=1:np-n) sub(n+j,i)= nv + j

       Call Electeds_objective (sub(1:,i), siz,avg_rate,portion, obj,dat)

       Call Add_to_fit (fit0,obj,dat,sub(1:,i), nf,best_elect, &
                        best_fit,best_dat,loc)
     End do
    
   End Subroutine Fit_range


   Subroutine List_elect (nc,ne, wtb,ballot, nu,sets, top_obj,objs) 

     Integer,  Intent(in) :: nc            ! # candidates in the domain
     Integer,  Intent(in) :: ne            ! Max # 'sets' to be ranked and recorded

     Real,     Intent(in) :: wtb(:)        ! (nb) Ballot weights, sum normalized to 'np'
     Integer,  Intent(in) :: ballot(0:,:)  ! (0:mr,nb)  (1:n,b) = candidates in rank order
                                           !              (0,b) = 'n' = # ranked 
     
     Integer, Intent(out) :: nu            ! # sets recorded, which together with the top set 
                                           !   cover the complete domain
     Integer, Intent(out) :: sets(:,0:)    ! (np,ne>=nu) The sets, ordered by decreasing objective

     Real,    Intent(out) :: top_obj       ! The top objective value
     Real,    Intent(out) :: objs(0:)      ! (ne>=nu) Decreasing ratios of objective values for the sets
!  Local:
     Logical :: Cover0(nc), Cover(nc)
     Real    :: siz(nc)
     Integer :: lst(ne)
     Integer :: b, i, j, k, cn, nf, np

     np= Size(sets,1);  siz= 0

     Do b= 1,Size(wtb)
       cn= ballot(1,b);  siz(cn)= siz(cn) + wtb(b)
     End do

     sets= 0;  objs= 0

     Call Best_sets (siz, nf,sets,objs)  ! List all sets & order by objective

     top_obj= objs(0);  objs(:nf)= objs(:nf) / top_obj

!    Only keep the set pairs (1st & 'i'th) which cover the domain

     j= 0;  lst= 0;  Cover0= .false.;  Cover0(sets(:,0))= .true.  

     Do i= 1,nf
       Cover= Cover0;  Cover(sets(:,i))= .true.

       If (All(Cover)) then
         j= j + 1;  lst(j)= i
       End if
     End do
     nu= j

     sets(:,1:nu)= sets(:,lst(:nu))
     objs(1:nu)  = objs(lst(:nu))

     If (pr_out > 1) then
       Call Out ("In 'List_elect' # ranked sets",nf, "# covering sets",nu, ln=1)
       Call Out ("Top objective",top_obj)

       If (pr_out > 1.5) then
         Call Out ("Covering objective ratios",objs(:nu))
         Call Out (-1,"and sets",sets(:,:nu))
       End if
     End if

   End Subroutine List_elect


   Subroutine List_elect1 (avgW, nu,sets, top_obj,objs) 

     Real,     Intent(in) :: avgW(:)       ! (nc) Average candidate weights on the 'domain'
     
     Integer, Intent(out) :: nu            ! # sets recorded, which together with the top set 
                                           !   cover the complete domain
     Integer, Intent(out) :: sets(:,0:)    ! (np,0:ne=>nu) The sets, ordered by decreasing objective

     Real,    Intent(out) :: top_obj       ! The top objective value
     Real,    Intent(out) :: objs(0:)      ! (0:ne>=nu) Decreasing ratios of objective values for the sets
!  Local:
     Integer, Allocatable :: lst(:)
     Logical :: Cover0(Size(avgW)), Cover(Size(avgW))
     Integer :: i, j, nf

     sets= 0;  objs= 0

     Call Best_sets (avgW, nf,sets,objs)  ! List all sets & order by objective

     top_obj= objs(0);  objs(:nf)= objs(:nf) / top_obj

!    Only keep the set pairs (1st & 'i'th) which cover the domain

     Allocate(lst(nf))

     j= 0;  lst= 0;  Cover0= .false.;  Cover0(sets(:,0))= .true.  

     Do i= 1,nf
       Cover= Cover0;  Cover(sets(:,i))= .true.

       If (All(Cover)) then
         j= j + 1;  lst(j)= i
       End if
     End do
     nu= j

     sets(:,1:nu)= sets(:,lst(:nu))
     objs(1:nu)  = objs(lst(:nu))

     If (pr_out > 1) then
       Call Out ("In 'List_elect1' # ranked sets",nf, "# covering sets",nu, ln=1)
       Call Out ("Top objective",top_obj)

       If (pr_out > 1.5) then
         Call Out ("Covering objective ratios",objs(:nu))
         Call Out (-1,"and sets",sets(:,:nu))
       End if
     End if

   End Subroutine List_elect1


   Subroutine Best_sets (obj, nf,sets,objs)
  
!    Find all sets and order them by objective up to 'ne'  
    
     Real,     Intent(in) :: obj(:)     ! (nc) Objective for each candidate

     Integer, Intent(out) :: nf         ! # sets recorded beyond the top set
     Integer, Intent(out) :: sets(:,0:) ! (np,0:ne>=nf) Records the sets corresponding to objs
     Real,    Intent(out) :: objs(0:)   ! (0:ne>=nf)    Decreasing objective values of the sets
!  Local:
     Integer, Allocatable :: sub(:,:)
     Integer :: i, k, n, nc, np, nsub, loc

     nc= Size(obj);  np= Size(sets,1)

     nsub= N_sub_sum(nc,np) - 1;  Allocate (sub(0:np,nsub))
     k= 0;  sub= 0

     Call List_subsets (1,nc,np,0, k,sub)  ! k = nsub

     nf= -1
     Do i= 1,k
       n= sub(0,i)  ! # candidates in subset
       If (n /= np) Cycle

       Call Add_to_list (obj,sub(1:,i), nf,sets,objs,loc)
     End do
    
   End Subroutine Best_sets

   Subroutine Add_to_list (obj,new_set, nf,best_set, best_obj,loc)

!    Compute the objective 'obj' of a 'new_set' and store it in 'best_set'
!    along with its objective in 'best_obj' in decreasing order 
!    of objective. The higher the objective the better.

     Real,         Intent(in) :: obj(:)           ! (nc) Objective for each candidate
     Integer,      Intent(in) :: new_set(:)       ! (np) New set of candidates, in increasing
                                                  !      order, to be added to 'best_set' if feasible
       
     Integer,   Intent(inout) :: nf               ! # sets recorded beyond the top set, nf <= ne
     Integer,   Intent(inout) :: best_set(:,0:)   ! (np,0:ne=>nf) Records good sets, ordered
                                                  !         by decreasing objective   
     Real,      Intent(inout) :: best_obj(0:)     ! (0:ne=>nf) Decreasing objective values of
                                                  !      sets recorded in best_set
     Integer,     Intent(out) :: loc              ! Location in 'best_set' of the new 
                                                  !   set of electeds 'new_set',
                                                  !   or minus this location if already
                                                  !   present in 'best_set',
                                                  !   or mx+1 if not recorded
 ! Local:
     Real    :: ob1
     Integer :: j, k, j1, k1, l1, l2, ne, n1

     ne= Ubound(best_obj,1);  n1= Min(nf+1,ne)

     ob1= Sum(obj(new_set));  loc= ne
 
     If (nf < 0) then  ! Initial call
       loc= 0
     Else  ! Subsequent call
       loc= Minloc(Abs(ob1 - best_obj(:nf)),1) - 1
       
       If (Abs(ob1 - best_obj(loc)) < 0.001*ob1) then
         If (All(new_set == best_set(:,loc))) then
           loc= -loc;  If (pr_out > 1)  Call Out("Repeated set not recorded")
           Return
         End if
       End if
       
       If (ob1 <= best_obj(loc)) then
         l1= loc + 1;  k= n1
         
         If (k >= l1) then  ! Push right and prepare to insert at loc+1
           l2= loc + 2;  k1= k - 1
           If (k >= l2) then  
             best_obj(l2:k)  = best_obj(l1:k1)
             best_set(:,l2:k)= best_set(:,l1:k1)
           End if 
           loc= l1;  ! Insert new_set at loc+1
         Else
           loc= ne;  If (pr_out > 1) Call Out("Data full: set not recorded")
           Return
         End if

       Else ! For ob1 > best_obj(loc)
         k= n1
         
         If (k >= loc) then  ! Push right and prepare to insert at loc
           l1= loc + 1;  k1= k - 1  
           If (k >= l1) then  
             best_obj(l1:k)  = best_obj(loc:k1)
             best_set(:,l1:k)= best_set(:,loc:k1)
           End if
         End if
       End if
     End if
     
     nf= n1;  best_obj(loc)= ob1;  best_set(:,loc)= new_set

     If (pr_out > 1.5) then
       Call Out ("New set at database position",loc, "with objective",ob1)
       Call Out ("The new set",new_set)
       Call Out ("Database of overall objective values",best_obj(:nf))
       Call Out (-1,"for possible elected sets",best_set(:,:nf))
     End if
   End Subroutine Add_to_list

End Module Clusters4

