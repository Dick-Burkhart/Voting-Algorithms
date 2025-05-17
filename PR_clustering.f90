
!  This program implements a clustering algorithm for Proportional
!  Representation (PR), created and tested by Dick Burkhart
    
!  The clustering identifies and compares configurations of voting blocks to see which best match the ballots.
!  Then sets of candidates are compared to see which best match the sizes of the voting blocks according to their
!  popularity in those blocks. 

!  The program assumes either preferential ballots (ranking) or rating-type ballots where the ratings are converted
!  to points with a computer-implemented restriction on the number of points per ballot.
    
!  The program processes either one or many sets of ballots, which may be in different common formats. Ballots are
!  pre-processed, the consoliated prior to clustering.
    
Program PR_Clustering

   Use Clusters0  ! Contains parameters
   Use Clusters1  ! Contains "Read_ballots" and its subroutines
   Use Clusters2  ! Contains "Consolidate_ballots" and its subroutines
   Use Clusters3  ! Contains "Form_clusters" and its subroutines
   Use Clusters4  ! Contains "Evaluate_candidates" and its subroutines
   Use Clusters5  ! Contains "Write_summary" and its subroutines
   Use Clusters6  ! Contains "Run_input" and subroutines for run statistics
                  ! like "Non_clust_stats" and "Sensitivity_parm"
   Use Newton_operators
   Use Types
   Use Precisn
   Implicit None

!  Variables for "Run_input":
   
   Integer                :: dst1, dst2            ! First and last districts to process
   Integer                :: nDst                  ! # districts to process = dst2 - dst1 + 1
   Character(18), Pointer :: District(:)=>Null()   ! (mx_Dist) List of voting districts to be processed (input)
   Character(3),  Pointer :: Party(:,:)=>Null()    ! (mxCand,mx_Dist) Party affiliation of each candidate for each district, if any

!  For parameter sensitivity statistics:
   
   Integer :: sens_parm    ! Sensitivity parameter to test, from 1 though Nsens
   Integer :: cntr         ! Sensitivity index from 0 to ndev, with 0 = standard run
   
!  Data structures for processed ballots: Read_ballots

   Integer, Pointer :: orig_cand(:)=>Null()  ! (nc) Borda ordering of the original candidates

   Real,    Pointer :: pt_val(:)=>Null()     ! (mt) Decreasing point values of the rating levels

   Real,    Pointer :: wtb(:)=>Null()        ! (nb) Ballot weights summing to 'np'
   Integer, Pointer :: ballot(:,:)=>Null()   ! (0:mr,nb)  (1:nr,b) = candidates in preferential order
                                             !   (0,b) = 'nr' = # ranked or rated  
   Integer, Pointer :: ballot2(:,:)=>Null()  ! (0:mr,nb) (1:nr,b) =  increasing ranking or rating levels
                                             !   (0,b) = # positive ratings 

   Real,    Pointer :: global(:,:)=>Null()   ! (nc,5) Global ballot mean (1) & sigma (2) vectors
                                             !    with noise mean (3) and its sigma (4), plus
                                             !   normalized mean rating above the noise (5)
   
!  Data structures for Non-clustering methods

     Type(Multi_listR) :: Non_clust  ! Data for 'nMt' selected non-clustering methods 'm'
                                     !                to elect multiple candidates.
                                     ! %n = # distinct sets of electeds found for the 
                                     !      current district, = 'nonC'
                                     ! %k = final # STV clusters, = 'nSTV'
                                     ! %l = final # DTV clusters, = 'nDTV'
                                     ! %vl(nMt)     = Mapping to the distinct elected sets in %L0
                                     ! %L0(np,nonC) = Distinct non-clustering sets of electeds

                                     ! %Q0(np,0:2,nMt)= For each method 'm', the elected candidates (:,0,m) 
                                     !                in increasing order, in ranked order (:,1,m),
                                     !                with the mapping to merged and reordered candidate clusters
                                     !                in (:,2,j) for j = 1 for STV, = 2 for DTV, else 0
                                     ! %L1(nMt,2)  (m,1)= Index 'e' in Elect_dat of non-clustering method 'm' 
                                     !             (m,2)= The matching best cluster set 'q'
                                     ! %L2(nonC,ncs) = Location, for each cluster set 'q', in Clust_elect(q)%Q0 
                                     !                of each distinct non-clustering set of electeds

                                     ! %T0(0:nc,0:ind,2) = Cluster data from STV (:,:,1) and DTV (:,:,2).
                                     !    (0,:,:) = cluster sizes, (1:,:,:) = cluster mean rating vectors
                                     !    (1:,0,:) = cluster averaged mean vector
                                     ! %T1(ncl,ncl,2) = Cluster correlations for STV (:,:,1) and DTV (:,:,2),
                                     !    = 'STV_corr' & 'DTV_corr' .
!  Data structcures for ballot consolidation into slate clusters

   Type(Multi_listD), Pointer :: Memb(:)=>Null() ! (0:nsl) Final slate cluster data
                                                 ! 0 Case: 
                                                 !   %l = np = # candidates to be elected = 
                                                 !             total weight in position units
                                                 !   %n = nc = # candidates
                                                 !   %qx(2)  = correlation parameters
                                                 !  %vl(mlv). (lv) = # slate clusters with
                                                 !      # ranked or rated candidates <= lv

                                                 ! 1:nsl Case:
                                                 !  %l   = top slate level    = # top candidates
                                                 !  %m   = bottom slate level = # bottom candidates (for Rating = 1)
                                                 !  %p   = # positively rated candidates (for %tx)
                                                 !  %n   = # negatively rated candidates (for %tx & Rating = 1)
                                                 !  %k = # slate ballots positively correlated with and subsequent 
                                                 !       to this slate ballot
                                                 !  %vl = listing of these slate ballots by decreasing correlation
                                                 !  %qx = the corresponding correlations

                                                 !  %fsx= slate membership
                                                 !  %fux= slate width
                                                 !  %lt(nc) = candidates in preferential order (from %sx)
                                                 !  %px(nc) = variance vector for %sx
                                                 !  %rx(nc) = %sx * %ux
                                                 !  %sx(nc) = mean rating vector
                                                 !  %tx(0:nc)= centered and normalized mean vector
                                                 !  %ux(nc) = variance normalized weight %fsx / %px

   Real,   Allocatable :: Mean_rnd(:,:,:)  ! (0:nc,nr,N_cand) 
                                           ! Initial cluster mean vectors, standard ranking/rating
                                           ! bounds assumed, also non-centered for ranking and 
                                           ! not normalized. (0,:,:) = cluster weight 

!  Additional Variables for "Form_clusters":
     
     Type(Multi_listR) :: Init(N_init)     ! Initial cluster set data
                                           ! %k = cluster set converged to
                                           ! %n = ncl = # regular clusters
                                           ! %lt(ncl) Mapping of initial to converged clusters 

                                           ! %Q0(np,0:2,4) Elected candidates, by 4 methods 
                                           !      (:,:,:1) = Rem_frac electeds (original for STV/DTV)
                                           !      (:,:,:2) = DHondt electeds
                                           !      (:,:,:3) = Optimal electeds 
                                           !      (:,:,:4) = Converged optimal electeds
                                           !     (:,0,:) = electeds in increasing order
                                           !     (:,1,:) = electeds in ranked order
                                           !     (:,2,:) = corresponding clusters

                                           ! %M1(7,4) Objective data for the electeds %Q0(:,0,:)
                                           ! %M2(ncl,ncl) Cluster correlation matrix, 0 on diagonal

                                           ! %T0(0:nc,0:ind,2) (1:,1:,1) = Candidate mean ratings  
                                           !      for each cluster, with cluster size at (0,1:,1)
                                           !      (1:,1:,2) = Cluster portions with averaged 
                                           !                  noise zeroed ratings at (1:,0,2)

     Type(Multi_listR) :: Clust_set(N_init)  ! (ncs) Final cluster sets and associated convergence data
                                             ! %k   = best initial cluster set that converged to this set
                                             ! %l   = # clusters in this set, including independents = 'ind' = ncl + 1
                                             ! %m   = # initial cluster sets that converged to this cluster set
                                             ! %ls(m)= list of those initial sets 
                                             ! %fux = clustering objective ratio = (q)%rx(1) / (1)%rx(1) 

                                             ! %lt(ncl) Mapping of initial to converged regular clusters 

                                             ! %px(6)  Objective data: min size, max correlation, 
                                             !         independents size (frac), final residual,
                                             !         Jnorm, Pnorm
                                             ! %qx(4)  Objective factors: min size penalty, max correlation penalty, 
                                             !         independents size penalty, residual penalty
                                             ! %sx(ind): Cluster weights = sizes

                                             ! %ux(6)  Frac of voters represented by regular clusters
                                             !         when that representation exceeds 'Parm_rep'

                                             !   %vl(0:6): Convergence data
                                             !     0 = original # regular clusters = n0
                                             !     1 = final # regular clusters
                                             !     2 = final success code
                                             !     3 = total # cluster merge and deletion operations
                                             !     4 = # convergence calls
                                             !     5 = total # centroiding iterations
                                             !     6 = total # line search points

                                             !   %L0(0:4,0:iv) Integer convergence data per update call, n1= %vl(4)
                                             !   0  = restart phase
                                             !   1  = success code:  2 = fully converged,   1 = reducing functional,
                                             !                       0 = damped reduction, -1 = convergence failure
                                             !   2 = # reg clusters
                                             !   3 = # centroiding iterations
                                             !   4 = # line search points

                                             ! %L2(n0,2) Mapping data between initial and final clusters 
                                             !   1 : Maps initial clusters onto intermediate clusters, reordered
                                             !   2 : For each final cluster, the largest initial cluster mapped to it

                                             ! %M2(0:2,0:ncyc) Real convergence data
                                             !     0 = Final line search damping parameter
                                             !     1 = Residual for the updated mean rating vectors
                                             !     2 = Net line search change

                                             ! %M1(ncl,ncl) Cluster correlation matrix

                                             ! %M3(8,ind) Cluster data
                                             !   1 = size 
                                             !   2 = density penalty
                                             !   3 = boundary penalty
                                             !   4 = density
                                             !   5 = boundary
                                             !   6 = width
                                             !   7 = Frac slate cluster fuzzy reduced memberships
                                             !   8 = Average reduction factor

                                             !  %L1(-2:nc,0:ind) Candidate ordering data for
                                             !                    noise zeroed vectors
                                             !     (-2:0,:) = # top, significant, viable cand's
                                             !     (1:,1:)  = Decreasing ordering from %T0(:,1:,2)
                                             !     (1:,0)   = Decreasing ordering from %T0(:,0,1)

                                             !  %T0(nc,0:ind,3) Cluster vector data
                                             !      (:,1:,1) = Cluster mean vectors
                                             !      (:,1:,2) = Cluster mean vectors with noise zeroing 
                                             !      (:,1:,3) = Cluster portion vectors
                                             !      (:,0,1)  = Zrate0 * declining regular portion factor
                                             !      (:,0,2)  = Zrate0 = cluster averaged noise zeroed mean vectors

!  Data structures for sets of electeds that could match cluster sets:

   Type(Multi_listR), Allocatable :: Clust_elect(:)  ! (ncs) Data for sets of electeds associated to each cluster set
                                                     ! %k = index of the best set of electeds %Q0(np,0,1) in Elect_dat
                                                     ! %l = 'ind' = # clusters, including independents
                                                     ! %m = 'nf'= # sets of possible electeds
                                                     ! %fsx = full objective for the top set of electeds
                                                     ! %fux = full objective ratio = (q)%fsx / (1)%fsx 

                                                     ! %sx(7)  Full objective data
                                                     !   (1) full objective = (2) clustering objective * (3) fitness
                                                     !   (4) cluster set coherent vote * (5) cluster set penalty
                                                     !   (6) avg elected rating * (7) cluster deviation size penalty
                                                     ! %vl(nf) The indices in "Elect_dat' of the possible sets of electeds

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
                                                     ! %M2(3,nf)= Fitness objectives(1) = electeds ratings(2) * electeds penalties(3) 

                                                     ! %T0(ind,2,nf) Objective data for each set of possible electeds
                                                     !   1 = Cluster size deviations = actual - true
                                                     !   2 = Cluster size deviation penalties
   

   Type(Multi_listR) :: Elect_dat(mxe)   ! Data for possible sets 'e' of electeds,
                                         !   ordered by 'best' full objective %M1(1,1)
                                         ! %k = The 'best' cluster set 'q' for 'e' = %lt(1)
                                         ! %l = The matching 'eq' in Clust_elect(q)%Q0 = %vl(1)
                                         ! %m = nq = # cluster sets that list 'e' 
                                         !      or = 1 for lower cluster sets
                                         ! %n = # regular clusters for best cluster set 'q'
                                         ! %q = Non-clustering index if 'e' is a
                                         !      non-clustering elected set, else 0
                                         ! %fsx = overall voter satisfaction
                                         ! %fux = full objective ratio = (e)%M1(1,1) / (1)%M1(1,1) 

                                         ! %ls(np) = the set of electeds in increasing order
                                         ! %lt(nq) Ordering of the 'nq' cluster sets 'q' which  
                                         !         list the set of electeds 'e' at a position 
                                         !         'eq' in Clust_elect(q)%Q0(:,:,eq) with eq <= nf. 
                                         !         The ordering is by full objective 
                                         !         Clust_elect(q)%sx(2) * Clust_elect(q)%M2(1,eq)
                                         ! %vl(nq)    The corresponding indices 'eq'
                                         ! %L0(np,nq) This elected set in preferential order for each
                                         !            of the cluster sets ordered as in %lt
     
                                         ! %M0(0:Mvt,2): Measures of success - over ballots
                                         !   (0,1) : mean # of top 3 elected in order (but not the next lower)
                                         !   (0,2) : mean # of top 3 elected regardless of order
                                         !   (l,1) : fraction with top l or more elected in order
                                         !   (l,2) : fraction with top l or more elected regardless of order

                                         ! %M1(7,nq) Objective data for the clusters sets ordered %lt.
                                         !   (1,q) = full objective = (2,q)*(3,q)
                                         !   (2,q) = clustering objective = (4,q)*(5,q)
                                         !   (3,q) = electeds fitness = (6,q)*(7,q)
                                         !   (4,q) = cluster set coherent vote 
                                         !   (5,q) = cluster set penalty
                                         !   (6,q) = elected rating = average rating over electeds  
                                         !   (7,q) = deviation size penalty
   Integer :: dpr(0:mxe)  ! (0:mxf) Elected sets whose best full objective exceeds that of cluster set #1.
                          !         (0) = # such elected sets.

   Real,    Allocatable :: Noise_dist(:,:)  ! (-Noise_opt:Noise_opt,nDst) Deviation ranking noise levels by district = Noise_cor
   Real,    Allocatable :: Noise_obj(:,:)   ! (-Noise_opt:Noise_opt,nDst) Corresponding clustering objectives
   Integer, Allocatable :: Noise_ncl(:,:)   ! (-Noise_opt:Noise_opt,nDst) Corresponding # clusters
   Integer, Allocatable :: Noise_top(:)     ! (nDst) Deviation index of the top clustering objective
   Real,    Allocatable :: Noise_cor1(:)    ! (nDst) Top noise correlation level by district = Noise_dist(Noise_top(id),id)
   Real,    Allocatable :: Noise_dif(:,:)   ! (nDst,3)  (id,1) Noise_pr(idist) - Noise_cor, (id,2) Noise_cor1(id) - Noise_cor
                                            !           (id,3) Noise_cor1(id) - Noise_pr(idist)
   Integer, Allocatable :: lst(:)           ! (2*Noise_opt+1)
   Real,        Pointer :: Noise_pr(:)=>Null()  ! (mx_Dist) Prior computed values of Noise_cor
   Integer, Save        :: Noise_opt        ! Noise_cor option: # +- deviations from the prior value if > 0.
                                            ! If = 0, use the prior value of Noise_cor
                                            ! If < 0, use the currently computed value of Noise_cor
   Integer, Save        :: STV_opt          ! If > 0, choose the top deviation to favor STV.
                                            ! If < 0, favor DTV. If = 0, favor neither.
   Real,   Parameter    :: del= 0.025       ! Deviation delta for Noise_cor
   Real,   Parameter    :: eps= 0.01        ! 

   Integer :: nc_lim 

!  For each district:

   Logical, Parameter   :: Set_rnd_seeds= .false.
   Integer :: Rnd_seed(2,mx_Dist)  ! Random seed for each district, to be use for 'Mean_rnd'

   Integer :: STV_eq(2,2) ! (:,1) = elected sets 'e' in Elect_dat for STV & DTV = Non_clust%L1(1:2,1)
                          ! (:,2) = converged cluster sets 'q' for STV & DTV    = Init(19:20)%k

   Integer :: ncs    ! Total # cluster sets converged to 
   Integer :: nb     ! # initial ballots
   Integer :: nc0    ! # original candidates
   Integer :: nc     ! # candidates
   Integer :: np     ! # candidates to be elected
   Integer :: mrp    ! Max # positively ranked or rated candidates
   Integer :: mrn    ! Max # negatively rated candidates
   Integer :: mr     ! Max # ranked or rated candidates
   Integer :: mt     ! # rating levels in 'pt_val'
   Integer :: nsl    ! # slate ballot clusters
   Integer :: tot_wt ! Total # original ballots in a district

   Integer :: me     ! # sets of possible electeds as computed
   Integer :: max_me ! Max value of 'me' over all runs
   Integer :: ncl    ! ncl = # regular clusters in a cluster set
   Integer :: ind    ! = ncl + 1 = total # clusters, including independents, 
                     ! = last 'cluster in the set
   Integer :: ncyc   ! # convergence updates = index of the final cycle
   Integer :: irun   ! Index for runs of statistical trials
   Integer :: idist  ! Index for district

   Real    :: xm, fac, adj, xns, rnd(3)
   Integer :: i1, i2, id, iq, ml, nd, nr, ns, pl, sl, ios, no1, ns1, nonC, nslat
   Integer :: e, i, j, k, l, m, n, p, q
 
   If (Set_rnd_seeds) then
     Call Random_number (rnd)
     Do id= 1,mx_Dist
       Call Random_number (rnd)
       Call Random_seed (Get = Rnd_seed(:,id))
     End do

     Open(7, File='RndSeed.txt', IOstat=ios, Action='Write')
       Write(7,*) "Random seeds for 'Mean_rnd' for each district"
       n= mx_Dist/4;  j= 1

       Do i= 1,n
         k= Min(j+3,mx_Dist);  Write(7,'(4(2X,2I11))') (Rnd_seed(:,id), id=j,k)
         If (k == mx_Dist) Exit;  j= k + 1
       End do
     Close(7)
   Else
     Open(7, File='RndSeed.txt', IOstat=ios, Status='Old', Action='Read')
       Read(7,*);  n= mx_Dist/4;  j= 1

       Do i= 1,n
         k= Min(j+3,mx_Dist);  Read(7,'(4(2X,2I11))') (Rnd_seed(:,id), id=j,k)
         If (k == mx_Dist) Exit;  j= k + 1
       End do
     Close(7)
   End if
    
   Allocate(Non_clust%vl(nMt),Non_clust%L1(nMt,2))
   Non_clust%vl= 0;  Non_clust%L1= 0
   GN_max_pt= Max_ptD;  GN_min_wt= Min_wt;  GN_min_mx= Min_mx  

   Call PR_input (District,Party,Noise_pr, pr_out, dst1,dst2, &
                   sens_parm,ndev, Noise_opt,STV_opt)

   If (Rating < 1) then
     Dot_fac= Dot_fac_rank
   Else
     Dot_fac= Dot_fac_rate
   End if

   GN_rating= Rating;  GN_dot_fac= Dot_fac
   nDst= dst2 - dst1 + 1;  max_me= 0

   If (Rating == 1) then
     Noise_cor= 0;  GN_rate_adj= 0
   Else 
     If (Noise_opt > 0) then
       n= Noise_opt;  no1= n+1
     
       Allocate (Noise_dist(-n:n,nDst), Noise_obj(-n:n,nDst), &
                 Noise_ncl(-n:n,nDst), Noise_top(nDst),       &
                 Noise_cor1(nDst), Noise_dif(nDst,3), lst(n+no1))
       Noise_dist= 0;  Noise_obj= 0;  Noise_ncl= 0
       Noise_top= 0;  Noise_cor1= 0;  Noise_dif= 0;  lst= 0
     Else if (Noise_opt < 0) then
       Allocate (Noise_cor1(nDst));  Noise_cor1= 0
     End if
   End if

!  Standard statistics and parameter perturbation tests 
   
   Call Sensitivity_parm (-1,sens_parm,ndev, dst1,nDst, 0)

   Cntr_parm_loop : Do cntr= 0,ndev

     Standard= cntr == 0
     nGen_clust= 0;  nGen_clust_CM= 0;  nGen_clust_NT= 0

     Call General_stats (0,dst1,nDst, 0,0,0)

     If (nDst == mx_Dist) Call General_data (0,nDst)

     Call Sensitivity_parm (0,sens_parm,ndev, dst1,nDst, cntr)

     District_loop : Do idist= dst1,dst2
       
       ngc(1)= nGen_clust;  ngc(2)= nGen_clust_CM  
       ngc(3)= nGen_clust_NT
       id= idist - (dst1 - 1)

       If (idist == 12) then
         nc_lim= 8
       Else
         nc_lim= nc0
       End if

       Call Read_ballots (Trim(District(idist)),idist,nc_lim, nc0,nc,np,mr,mrp, nb,pt_val, &
                          orig_cand,global, wtb,ballot,ballot2, Non_clust)

       mt= Size(pt_val);  mrn= mr - mrp
       nr= Min(np+2,7);  Allocate(Mean_rnd(0:nc,nr,N_cand))

       If (Rating < 1) then
         If (Noise_opt > 0) then ! Use Noise_cor = top deviation from best prior

!          Record deviations from the current Noise_cor

           Noise_dif(id,1)= Noise_pr(idist) - Noise_cor
           Noise_dif(id,2)= Noise_cor

           Do i= -Noise_opt,Noise_opt
             xns= Noise_pr(idist) + i*del
             Noise_dist(i,id)= Max(Min(xns,Noise_max), Noise_min)

             Noise_cor= Noise_dist(i,id)
             Noise_por= Por_fac * Noise_cor;  GN_rate_adj= Noise_cor

             Call Random_seed (Put = Rnd_seed(:,idist))

             Call Consolidate_ballots (idist, nc,np,mr,nb, pt_val,       &
                                       wtb,ballot(:,:nb),ballot2(:,:nb), &
                                       nr,nsl,Memb, Mean_rnd)

             Call Form_clusters (idist,nc,np, Memb, Mean_rnd(:,:nr,:),  &
                                 Non_clust, Init, ncs,Clust_set)

             Call DeAlloc_Multi_list_ar (Memb);     DeAllocate(Memb)
             Call DeAlloc_Multi_list_ar (LS_Memb);  DeAllocate(LS_Memb)

             nonC= Size(Non_clust%L0,2);  Allocate (Clust_elect(ncs), Non_clust%L2(nonC,ncs))

             Call Evaluate_candidates (idist, nc,np,ncs, Non_clust,Clust_set(:ncs), &
                                       Clust_elect, me,Elect_dat, dpr)

             If (STV_opt > 0) then       ! Let Noise_top(id) tend to favor STV
               e= Non_clust%L1(1,1)
             Else if (STV_opt < 0) then  ! Let it favor DTV
               e= Non_clust%L1(2,1)
             Else                        ! Let it favor neither
               e= 1
             End if

             q= Elect_dat(e)%k;  ncl= Clust_set(q)%l-1 
             Noise_obj(i,id)= Elect_dat(e)%M1(1,1);  Noise_ncl(i,id)= ncl

             Call DeAlloc_Multi_list_ar (Init)
             Call DeAlloc_Multi_list_ar (Clust_set)
             Call DeAlloc_Multi_list_ar (Clust_elect)
             Call DeAlloc_Multi_list_ar (Elect_dat)
             DeAllocate (Clust_elect);  DeAllocate(Non_clust%L2)
           End do

           i= Maxloc(Noise_obj(:,id),1) - no1
           xm= Noise_obj(i,id) - eps;  Call List_of_true (-Noise_opt, Noise_obj(:,id) > xm, n,lst)

           If (n > 1) then
             j= Minloc(Abs(lst(:n)),1);  i= lst(j)
           End if

           Noise_top(id)= i;  Noise_cor= Noise_dist(i,id)
           Noise_por= Por_fac * Noise_cor;  GN_rate_adj= Noise_cor

         Else if (Noise_opt == 0) then ! Use best prior 'Noise_cor'
           Noise_cor= Noise_pr(idist)
           Noise_por= Por_fac * Noise_cor;  GN_rate_adj= Noise_cor
     
         Else   ! Use the statistically computed current 'Noise_cor'
           Noise_cor1(id)= Noise_cor
           Noise_por= Por_fac * Noise_cor;  GN_rate_adj= Noise_cor
         End if
       End if
       
       Call Random_seed (Put = Rnd_seed(:,idist))

       Call Consolidate_ballots (idist, nc,np,mr,nb, pt_val,       &
                                 wtb,ballot(:,:nb),ballot2(:,:nb), &
                                 nr,nsl,Memb, Mean_rnd)

       Call Form_clusters (idist,nc,np, Memb, Mean_rnd(:,:nr,:), &
                           Non_clust,Init, ncs,Clust_set)

       DeAllocate(Mean_rnd);  Call DeAlloc_Multi_list_ar (Memb);  DeAllocate(Memb)
       Call DeAlloc_Multi_list_ar (LS_Memb);  DeAllocate(LS_Memb)

       nonC= Size(Non_clust%L0,2);  Allocate (Clust_elect(ncs), Non_clust%L2(nonC,ncs))

       Call Evaluate_candidates (idist, nc,np,ncs, Non_clust,Clust_set(:ncs), &
                                 Clust_elect, me,Elect_dat, dpr)
       max_me= Max(max_me,me)
       Call Analyze_elected (nc,Votr_wt, wtb,ballot(:,:nb), Elect_dat(:me))

       Do iq= 1,N_init
         q= Init(iq)%k;  If (q < 1) Cycle
         Init(iq)%Q0(:,:,4)= Clust_elect(q)%Q0(:,:,1)
         Init(iq)%M1(:,4)= Clust_elect(q)%sx
       End do

       STV_eq(:,1)= Non_clust%L1(:2,1);  STV_eq(:,2)= Init(19:20)%k

       Call Out ("District",idist, "with Noise_cor",Noise_cor, ln=1)
       q= Elect_dat(1)%k;  ncl= Clust_set(q)%l-1
       Call Out ("with best cluster set",q, "# clusters",ncl)
       Call Out ("Top elected set objectives",Elect_dat(1)%M1(:,1))
       Call Out (-1,"STV_eq",STV_eq)

       ml= Maxval(Clust_set%l);  nd= dpr(0)

       Call Write_summary (Trim(District(idist)),idist, np,nc,mr,mrn, ml,ncs,me, &
                           orig_cand,Party(:nc0,idist), Clust_set(:ncs),         &
                           Clust_elect,Elect_dat(:me), Non_clust,Init, dpr(0:nd))
         
!      Accumulate the statistics for the current run

       If (nDst == mx_Dist) Call General_data (1,nDst, Non_clust,Elect_dat)

       Call General_stats (1,dst1,nDst,idist, np,ncs, Non_clust, Init, &
                           Clust_set(:ncs), Clust_elect,Elect_dat(:me))

       e= Non_clust%L1(1,1)
       Call Sensitivity_parm (1,sens_parm,ndev, dst1,nDst, cntr, idist, &
                              STV_eq, Clust_set(:ncs),Elect_dat(:me))
       
       Call DeAlloc_Multi_list_ar (Init)
       Call DeAlloc_Multi_list_ar (Clust_set)
       Call DeAlloc_Multi_list_ar (Clust_elect)
       Call DeAlloc_Multi_list_ar (Elect_dat)
       DeAllocate (Clust_elect);  DeAllocate(Non_clust%L2)

     End do District_loop

!    Compute the final statistics for the current run

     Call Out ("Max # elected sets over all runs",max_me, ln=1)

     If (nDst == mx_Dist) Call General_data (2,nDst)

     Call General_stats (2,dst1,nDst, 0,0,0)

     Call Sensitivity_parm (2,sens_parm,ndev, dst1,nDst, cntr)

     If (Rating < 1) then
       If (Noise_opt > 0) then
         Call Out (-1,"Deviation ranking noise signficance levels by district",Noise_dist)
         Call Out (-1,"Their clustering full objectives",Noise_obj)
         Call Out (-1,"Their # regular clusters",Noise_ncl)
         Call Out ("The top deviation",Noise_top)

         Do id= 1,nDst
           idist= (dst1-1) + id;  i= Noise_top(id)
           Noise_cor1(id) = Noise_dist(i,id)
           Noise_dif(id,2)= Noise_cor1(id) - Noise_dif(id,2)
           Noise_dif(id,3)= Noise_cor1(id) - Noise_pr(idist)
         End do

         Call Out ("Top noise correlation level by district",Noise_cor1)
         Call Out (-1,"Differences of noise correlation levels by district",Noise_dif)
         Call Out ("Dif1: Prior - current. Dif2: Top - current. Dif3: Top - prior")

       Else if (Noise_opt == 0) then
         Call Out ("Prior best noise correlation level by district",Noise_pr(dst1:dst2))
       Else
         Call Out ("Current statistical noise correlation level by district",Noise_cor1)
       End if
     End if

   End do Cntr_parm_loop


   Close(8)
   
End Program PR_Clustering