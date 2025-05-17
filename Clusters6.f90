!         This module contains "Run_input" and "Correlation_stats" and their subroutines. 
    
!  "Run_input" is called from Clustering_PR at the beginning to read in data about the election 
!   ballots to be tested and to read in a few run parameters    
    
!  "Non_clust_stats" is called from Clustering_PR at the end to compile general statistics, that
!  including the comparison of the non-clustering methods to the clustering method
    
!  "Sensitivity_parm" is called to test the statistics of perturbations of certain parameters 
!   to standard parameters.  
    
Module Clusters6

   Use Clusters4
   Use Clusters_support
   Use Clusters0
   
   Use Util
   Use Output
   Use Types
   Use Precisn
   
   Implicit None
   Integer, Parameter :: fm2(3)= (/ 8, 2, 4 /), fm3(3)= (/ 8, 3, 4 /)
    
 Contains
    
   Subroutine PR_input (District,Party,Noise_pr, pr_out, Dist_1,Dist_e, sens_parm,ndev, Noise_opt,STV_opt)
   
     Character(18),     Pointer :: District(:)  ! (mxDist) List of voting districts to be processed (file read)
     Character(3),      Pointer :: Party(:,:)   ! (mxCand,mxDist) Party affiliation of each candidate 
                                                !   for each district, if any (file read)
     Real,              Pointer :: Noise_pr(:)  ! (mxDist) Prior computed values of Noise_cor by district
     
     Real,          Intent(out) :: pr_out       ! Specifies level of diagnostic print out to 'out_file'
     
     Integer,       Intent(out) :: Dist_1       ! First district to process
     Integer,       Intent(out) :: Dist_e       ! Last district to process

     Integer,       Intent(out) :: sens_parm    ! Sensitivity parameter to test (0, or 1...Nsens)
     Integer,       Intent(out) :: ndev         ! # deviation deltas for 'sens_parm' (0, or 2...6)

     Integer,       Intent(out) :: Noise_opt    ! Noise_cor option: # +- deviations from the prior value if > 0.
                                                ! If = 0, use the prior value of Noise_core
                                                ! If < 0, use the currently computed value of Noise_cor
     Integer,       Intent(out) :: STV_opt      ! If > 0, choose the top deviation to favor STV.
                                                ! If < 0, favor DTV. If = 0, favor neither.

! Local:
     Integer,  Allocatable :: D_data(:,:) ! (4,mxDist) For each district
!              (1)= file type, (2)= # candidates, (3)= # to be elected (4) Party affiliations known (1) or not (0).
!      file type:  
!        1  : for candidates listed in decreasing order of ranking
!        2  : for ranking or rating level, or 0, for each of nc0 candidates
   
     Character(50) :: label1, label2, label3
     Character(12) :: out_file
     Integer :: i, j, k, n, ios, mxCand, mxDist, nDist
     
     Open(7, File='Vote_files.txt', IOstat=ios, Status='Old', Action='Read')
       Read(7,*) mxDist, mxCand
       Allocate (District(mxDist), Party(mxCand,mxDist), D_data(4,mxDist), Noise_pr(mxDist))
       Party= '   '
       Do i= 1,mxDist
         Read(7,*) District(i), D_data(:,i)
         If (D_data(4,i) > 0) then
           n= D_data(2,i);  Read(7,*) Party(1:n,i)
         Else
           Read(7,*)  
         End if
       End do
     Close(7)
     
     Open(7, File='Data.txt', IOstat=ios, Status='Old', Action='Read')
       Read(7,*);  n= Ceiling(mxDist / 8.0);  j= 1
       Do i= 1,n-1
         k= j + 7;  Read(7,*) Noise_pr(j:k);  j= k + 1
       End do
       Read(7,*) Noise_pr(j:)
     Close(7)

     Open(7, File='PR_opt.txt', IOstat=ios, Status='Old', Action='Read')
       Read(7,*) label1, out_file, pr_out, Rating, Noise_opt, STV_opt
       Read(7,*) label3, Dist_1, Dist_e, sens_parm, ndev
     Close(7)

     Dist_1= Max(Min(Dist_1,mxDist), 1)
     Dist_e= Min(Max(Dist_e,Dist_1), mxDist)

     sens_parm= Max(Min(sens_parm, Nsens), 0)
     If (ndev < Mndev .or. ndev > Mxdev) ndev= 0
     If (ndev == 0) sens_parm= 0

     nDist= Dist_e - Dist_1 + 1
     
     Open(8, File=Outdat//out_file, IOstat=ios, Status='Replace', Action='Write')
     Call Out ("")
     Write (8,'(A, 3X,A, F4.1, I4, 2X,2I3)') Trim(label1), out_file, pr_out, Rating,  Noise_opt,STV_opt
     Write (8,'(A, 2I4, 4X,2I4)') Trim(label3), Dist_1,Dist_e,  sens_parm,ndev

     Call Out ("Use Fit_range in Elect_best",Test_range,ln=1)

     Call Out ("Standard Run Parameters")
     Call Out ("Parameters for 'Read_ballots'")
     Call Out ("Ballot point parameters",Parm_bal)
     
     Call Out ("Parameters for 'Consolidate_ballots'")
     Call Out ("Parameters for 'Initial cluster sets'")
     Call Out ("Parameters for Core_clustering", Core_parm)
     Call Out ("Correlation bounds", Parm_cor)
     
     Call Out ("Parameters for 'Cluster_slates'")
     Call Out ("Iteration limits",Mxitr)
     Call Out ("Convergence tolerances",Real(Parm_cnv(:,1)))
     Call Out ("Min merge criterion for cluster correlation", &
                Real(Parm_cnv(:,2)))
     Call Out ("Min size criterion for cluster deletion",&
                Real(Parm_cnv(:,3)))
     Call Out ("max correlation criterion for cluster deletion", &
               Real(Parm_cnv(:,4)))
     Call Out ("More convergence parameters",Real(Parm_cnv(:,5)))
     Call Out ("Objective parameters",Real(Parm_obj))
     
     Call Out ("Parameters for 'Evaluate_candidates'")
     Call Out ("Top cand ratio, secondary ratio, cutoff ratio", Parm_top)
     Call Out ("Deviation penalty parm for elected size - true", Parm_el)
     Call Out ("")
     
   End Subroutine PR_input


   Subroutine Domain_input (District,Party,Noise_pr, pr_out, dst1,dst2, Tourn_union,full)

!    Parameter input for 'Tourn_rem'
   
     Character(18),  Pointer :: District(:)   ! (mxDist) List of voting districts to be processed (file read)
     Character(3),   Pointer :: Party(:,:)    ! (mxCand,mxDist) Party affiliation of each candidate 
                                              !   for each district, if any (file read)
     Real,           Pointer :: Noise_pr(:)   ! (mxDist) Prior computed values of Noise_cor by district
     
     Real,       Intent(out) :: pr_out        ! Specifies level of diagnostic print out to 'out_file'
     Integer,    Intent(out) :: dst1          ! First district to process
     Integer,    Intent(out) :: dst2          ! Last  district to process
     Integer,    Intent(out) :: Tourn_union   ! 1 = Make all domains = union of 2 possible elected sets of size 'np'
                                              ! < 0 Domains = all subsets of size n0-2 of a domain of size n0,
                                              !     or other domain or processing options 
     Real,       Intent(out) :: full          ! Limit parameter for the rising cosine dominance function
! Local:
     Integer,  Allocatable :: D_data(:,:)     ! (4,mxDist) For each district
!              (1)= file type, (2)= # candidates, (3)= # to be elected (4) Party affiliations known (1) or not (0).
   
     Character(50) :: label1, label2, label3
     Character(12) :: out_file
     Real          :: fac, inc
     Integer       :: i, j, k, n, mr, ios, mxCand, mxDist
     
     Open(7, File='Vote_files.txt', IOstat=ios, Status='Old', Action='Read')
       Read(7,*) mxDist, mxCand
       Allocate (District(mxDist), Party(mxCand,mxDist), D_data(4,mxDist), Noise_pr(mxDist))
       Party= '   '
       Do i= 1,mxDist
         Read(7,*) District(i), D_data(:,i)
         If (D_data(4,i) > 0) then
           n= D_data(2,i);  Read(7,*) Party(1:n,i)
         Else
           Read(7,*)  
         End if
       End do
     Close(7)

     Open(7, File='Domain_opt.txt', IOstat=ios, Status='Old', Action='Read')
       Read(7,*) label1, out_file, pr_out, dst1,dst2, Tourn_union, full
     Close(7)
     
     Open(7, File='Data.txt', IOstat=ios, Status='Old', Action='Read')
       Read(7,*);  n= Ceiling(mxDist / 8.0);  j= 1
       Do i= 1,n-1
         k= j + 7;  Read(7,*) Noise_pr(j:k);  j= k + 1
       End do
       Read(7,*) Noise_pr(j:)
     Close(7)

     Open(8, File=Outdat//out_file, IOstat=ios, Status='Replace', Action='Write')
     Call Out ("")
     Write (8,'(A, 2X,A, F5.1, 2I4, 2X,I4,F7.2)') label1, out_file, pr_out, dst1,dst2, Tourn_union, full

     Call Out ("Standard Run Parameters")
     Call Out ("Parameters for 'Read_ballots'")
     Call Out ("Ballot point parameters",Parm_bal)
     
     Call Out ("Parameters for 'Consolidate_ballots'")
     Call Out ("Parameters for 'Initial cluster sets'")
     Call Out ("Parameters for Core_clustering", Core_parm)
     Call Out ("Correlation bounds", Parm_cor)
     
     Call Out ("Parameters for 'Cluster_slates'")
     Call Out ("Iteration limits",Mxitr)
     Call Out ("Convergence tolerances",Real(Parm_cnv(:,1)))
     Call Out ("Min merge criterion for cluster correlation", &
                Real(Parm_cnv(:,2)))
     Call Out ("Min size criterion for cluster deletion",&
                Real(Parm_cnv(:,3)))
     Call Out ("max correlation criterion for cluster deletion", &
               Real(Parm_cnv(:,4)))
     Call Out ("More convergence parameters",Real(Parm_cnv(:,5)))
     Call Out ("Objective parameters",Real(Parm_obj))
     
     Call Out ("Parameters for 'Evaluate_candidates'")
     Call Out ("Top cand ratio, secondary ratio, cutoff ratio", Parm_top)
     Call Out ("Deviation penalty parm for elected size - true", Parm_el)
     Call Out ("")
     
   End Subroutine Domain_input

                         
   Subroutine Sensitivity_parm (Status,sens_parm,ndev, Dist_1,nDst,cntr, &
                                idist,STV_eq, Clust_set,Elect_dat)
   
!    Test the sensitivity of the clustering results to variations in the values of selected parameters of the algorithm.
   
   
     Integer, Intent(in) :: Status        ! Do initial allocations if < 0, 
                                          !   initialize the statistics if = 0, 
                                          !   update if = 1, finalize if = 2 
     Integer, Intent(in) :: sens_parm     ! Sensitivity parameter to test if
                                          !   1 <= sens_parm <= N_sens, else the standard run
     Integer, Intent(in) :: ndev          ! # sens_parm deviations to test, else 
                                          !  the standard run if ndev = 0
     Integer, Intent(in) :: Dist_1        ! First electoral district
     Integer, Intent(in) :: nDst          ! Total # electoral districts
     
     
     Integer, Intent(in) :: cntr          ! Sensitivity counter for sens_parm: 0 to ndev 

     Integer, Optional, Intent(in) :: idist       ! District index, when updating
     Integer, Optional, Intent(in) :: STV_eq(2,2) ! (:,1) = elected sets 'e' in Elect_dat for STV & DTV = Non_clust%L1(1:2,1)
                                                  ! (:,2) = converged cluster sets 'q'    for STV & DTV = Init(19:20)%k

     Type(Multi_listR), Optional, Intent(in) :: Clust_set(:) ! (ncs) Final cluster sets
                                                ! %k   = best initial cluster set that converged to this cluster set
                                                ! %l   = # clusters in the set = 'ind'
                                                ! %m   = # initial cluster sets that converged to this cluster set
                                                ! %ls(m)= list of those initial sets 
                                                ! %fux = clustering objective ratio = (q)%rx(1) / (1)%rx(1) 

                                                ! %lt(ncl) Mapping of initial to converged regular clusters 

                                                ! %px(4):  Objective data: min size, max correlation, 
                                                !            independents size (frac), final residual
                                                ! %qx(4):  Objective factors: min size penalty, max correlation penalty, 
                                                !            independents size penalty, residual penalty
                                                ! %rx(3):  Objective value(1) = coherent vote(2) * penalty(3)
                                                ! %sx(ind):  Cluster sizes

                                                !   %tx(nc)  %T0(:,0,2) scaled by 1 - %T0(:,ind,3)
                                                !   %ls(-2:nc) Decreasing ordering of %tx at (1:nc)
                                                !              with # top, significant, and viable at (-1:0)

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

                                                ! %M1(ncl,ncl) cluster correlation matrix

                                                ! %M3(8,ind) Cluster data
                                                !   1 = size 
                                                !   2 = density penalty
                                                !   3 = boundary penalty
                                                !   4 = density
                                                !   5 = boundary
                                                !   6 = width
                                                !   7 = Frac slate cluster fuzzy reduced memberships
                                                !   8 = Average reduction factor

                                                !   %L1(-2:nc,0:ind) Candidate ordering data for
                                                !                    noise zeroed ratings
                                                !     (-2:0,:) = # top, significant, viable cand's
                                                !     (1:,1:)  = Decreasing ordering from 
                                                !                zeroed ratings = %T0(:,1:,2)
                                                !     (:,0)   = Ordering & data from %T0(1:,0,1)

                                                !   %T0(nc,0:ind,3) Cluster point data
                                                !      (:,1:,1) = Cluster mean vectors
                                                !      (:,1:,2) = Cluster mean vectors with noise zeroing 
                                                !      (:,1:,3) = Cluster portions for each candidate
                                                !      (:,0,1) = Scaled Zrate0 (for ordering the %L1(1:,0))
                                                !      (:,0,2) = Zrate0 = cluster averaged noise zeroed mean ratings

     Type(Multi_listR), Optional, Intent(in) :: Elect_dat(:) ! (me) Data for possible sets 'e' of electeds,
                                         ! ordered by 'best' full objective %M1(1,1)
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
     
                                         ! %M0(0:Mvt,2): Measures of success
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
!  Local:
     Type(Multi_listR), Save :: Stand_run(mx_Dist)  ! Save data from the standard run
                                                    ! %n = nc;  %p = np;  %k = ncl
                                                    ! %qx(7) = objective data for the top elected set
                                                    !   1 = Full objective
                                                    !   2 = Overall voter satisfaction
                                                    !   3 = STV to clustering ratio
                                                    !   4 = Coherent vote (from clustering objective)
                                                    !   5 = Clustering penalty (from clustering objective)
                                                    !   6 = Elected rating (from electeds objective)
                                                    !   7 = Deviation penalty (from electeds objective)

                                                    !   1 = Full objective
                                                    !   2 = Clustering objective
                                                    !   3 = Fitness
                                                    !   4 = Coherent vote
                                                    !   5 = Clustering penalty
                                                    !   6 = Elected rating
                                                    !   7 = Proportionality penalty
                                                    !   8 = DTV best full objective
                                                    !   9 = STV converged full objective
                                                    ! %ls(np) Top set of possible electeds (from Elect_dat%ls)
                                                    ! %M0(0:Mvt,2) Measures of success (from Elect_dat%M0)
                                                    !   (:,1) is for the top 'l' elected:  if cand(1:l) are elected 
                                                    !         but no lower ranked candidate for l <= n= Min(Mvt,cand(0)).  
                                                    !   (:,2) is for any of the top candidates cand(1:n) being elected.
                                                    !   (0,1) : mean fraction of "top elected" 'l' out of 'n'
                                                    !   (1,1) : fraction with top choice or more elected
                                                    !   (2,1) : fraction with top two choices or more elected
                                                    !   (3,1) : fraction with top three choices or more elected
                                                    !   (0,2) : mean fraction of the top 'n' elected
                                                    !   (1,2) : fraction with at least of 1 of the top 'n' choices elected
                                                    !   (2,2) : fraction with at least of 2 of the top 'n' choices elected
                                                    !   (3,2) : fraction with at least of 3 of the top 'n' choices elected
     Type(Multi_listR), Save :: Perturbed(mx_Dist)  !  Save data from a perturbed run
                                                    ! %n = nc;  %p = np;  %k = ncl
                                                    ! %qx(7) = cluster set obj data (from Elect_dat%M1, %fsx, STV_ratio)
                                                    ! %ls(np) Top set of possible electeds (from Elect_dat%ls)
                                                    ! %M0(0:Mvt,2) Measures of success (from Elect_dat%M0)
!    Data structures for the statistics
     
     Real, Save :: Stat1(7,4,0:1)        ! Stats for changes to objective data for perturbed vs standard parameters
                                         ! (:,1:3,:) clustering objective, fitness, DTV ratio, converged STV ratio
                                         ! (1:7,:,:) The 7 percentiles
                                         ! (:,:,0:1) Standard or perturbed run
     Real, Save :: Stat2(7,0:Mvt,2,0:1)  ! Stats for changes to the top electeds
                                         ! (:,0:Mvt,1,:) Same top choices elected, as in Stand_run%M0(0:Mvt,1)
                                         ! (:,0:Mvt,2,:) Some top ranked elected, as in Stand_run%M0(0:Mvt,2)
                                         ! (1:7,:,:,:) The 7 percentiles
                                         ! (:,:,:,0:1) Standard or perturbed run

     Real, Save :: Stat3(8,0:2)          ! Objective comparisons for perturbed - standard parameter values
                                         ! (:,0) for standard parameter, (:,1) for perturbed parameter, 
                                         !       (:,2) for perturbed - standard (mean only)
                                         ! (1:2,:) = mean & sigma/mean for clustering objective
                                         ! (3:4,:) = mean & sigma/mean for fitness
                                         ! (5:6,:) = mean & sigma/mean for DTV ratio
                                         ! (7:8,:) = mean & sigma/mean for converged STV ratio

     Real, Save :: Stat4(2,0:Mvt,2,0:1)  ! Stats for changes to the top electeds
                                         ! (:,0:Mvt,1) Same top choices elected, as in Stand_run%M0(0:Mvt,1)
                                         ! (:,0:Mvt,2) Some top ranked elected, as in Stand_run%M0(0:Mvt,2)
                                         ! (1:2,:,:,:) Mean, sigma
                                         ! (:,:,:,0:1) Standard vs Perturbed
     Real, Save :: Stat5(5,0:2,0:1)      ! Statistics for # cluster set matches, from initial to  
                                         !   converged cluster sets, as in Cluster_set(q)%m. 
                                         !   (i,0,:)   Counter for # districts with an ith best cluster set
                                         !   (i,1,:)   Sum the # matches over these districts 
                                         !   (i,2,:)   Sum the squares of the # matches
                                         !   (:,:,0:1) Standard vs Perturbed
     
     Real, Save :: Stat_mat              ! Fraction of districts whose perturbed electeds are the same 
                                         ! as the standard electeds
     Real, Save :: Stat_obj(4,3)         ! Accumulates for the mean(:,1) and standard deviation(:,2) of the abs difference
                                         !   between the perturbed and standard runs of the clustering objective(1,:), 
                                         !   fitness(2,:), DTV ratio(3,:), and converged STV ratio(4,:)
                                         !   Also the final perturbed - standard means in (:,3)

     Integer, Save :: N_reg(0:mxcl,2,0:1) ! Distribution # regular clusters 'n' (n=1:mxcl, row 1 for the
                                          ! top cluster set, row 2 for all cluster sets with n>=2). 
                                          ! (:,:,0) = standard run, (:,:,1) = perturbed run
                                          ! N_reg(0,:,:)= Sum(N_reg(1:,:,:), 1)
     
     Real, Save    :: Out_obj(0:4,0:Mxdev,2)  ! Output objectives (1:4 = perturbed - standard differences): 
                                          ! (0,:,:)= frac match, (1,:,:) clustering obj, (2,:,:)= fitness
                                          ! (3,:,:)= DTV ratio, (4,:,:)= converged STV ratio
                                          ! (:,cntr,:) = parameter deviation index
                                          ! (:,cntr,1) = obj values. cntr=0 for standard, else perturbed
                                          ! (:,cntr,2) = obj differences = perturbed - standard
     Real, Save    :: frac2p(0:2,2)      ! Fractions of cluster sets with # reg clusters >= 2.
                                         ! (0,:)= standard run. (1,:)= perturbed run. (2,:)= difference.
                                         ! (:,1)= top cluster set. (:,2)= all cluster sets. 
     
     Real, Save    :: avg2p(0:2,2)       ! Average # reg clusters in sets with # reg clusters >= 2.
                                         ! (0,:)= standard run. (1,:)= perturbed run. (2,:)= difference.
                                         ! (:,1)= top cluster set. (:,2)= all cluster sets.

     Real, Parameter :: Perc(7)= (/0.01, 0.10, 0.25, 0.50, 0.75, 0.90, 0.99/)  ! Percentiles for histograms
     
     Real, Parameter :: dig= 1000.0      ! Scale factor to round off to a certain # digits
     Real, Parameter :: match_val= 3.0/8 ! Value for a good match between a general clustering method and 
                                         ! the best cluster set when it converges to the latter
     Real    :: obp(mx_Dist,4)           ! Objective district data for perturbations 
     Real    :: tpp(mx_Dist,0:Mvt,2)     ! Top ranking district data for perturbations 
     Real    :: obd0(mx_Dist,4)          ! Objective district data for standard parameters
     Real    :: tpd0(mx_Dist,0:Mvt,2)    ! Top ranking district data for standard parameters

     
     Real    :: x, x1, x2, dev, dif, top_full, obj3(0:2)
     Integer :: ind, ncl, ncs, ns2, top, last
     Integer :: i, j, k, l, m, n, q, e1, e2, id, me, n4, n5, nc, np, nt, q1, q2
     
     If (Status < 0) then  ! Initialization of data structures

       If (sens_parm /= 0 .and. (sens_parm < 1 .or. sens_parm > Nsens)) then
         Call Out ("Error in 'Sensitivity_parm' for sensitivity parameter", sens_parm,ln=1);  Stop
       End if

       If (ndev /= 0 .and. (ndev < Mndev .or. ndev > Mxdev)) then
         Call Out ("Error in 'Sensitivity_parm': for bad # deviations", ndev,ln=1);  Stop
       End if

       If (Dist_1 < 1 .or. nDst < 1 .or. (Dist_1-1) + nDst > mx_Dist) then
         Call Out ("Error in 'Sensitivity_parm': bad 1st district",Dist_1, &
                   "or # districts", nDst, ln=1);  Stop
       End if

!      Standard parameter values

       Parm_bal= Parm_bal1

       If (sens_parm >= 1 .and. sens_parm <= Nsens .and. &
           ndev >= Mndev .and. ndev <= Mxdev) then
         SenP= -1;  SenP(0)= Std_sens(sens_parm)

         Do i= 1,ndev
           Select case (ndev)
             Case (2)
               SenP(i)= SenP(0) + SenV2(i,sens_parm)
             Case (3)
               SenP(i)= SenP(0) + SenV3(i,sens_parm)
             Case (4)
               SenP(i)= SenP(0) + SenV4(i,sens_parm)
             Case (5)
               SenP(i)= SenP(0) + SenV5(i,sens_parm)
             Case (6)
               SenP(i)= SenP(0) + SenV6(i,sens_parm)
           End Select
         End do

         Call Out ("For sensitivity parameter", sens_parm, ln=1)
         Call Out ("Perturbed sensitivity values",SenP(0:ndev))
       End if

       Stat1= 0;  Stat2= 0;  Stat3= 0;  Stat4= 0;  Stat5= 0;  N_reg= 0

     Else if (Status == 0) then   ! Initialize statistics and perturb the parameter
                                  ! based on 'cntr'

       If (cntr < 0 .or. cntr > ndev) then
         Call Out ("Error in 'Sensitivity_parm': bad deviation counter", cntr, ln=1);  Stop
       End if

       Out_obj(:,cntr,:)= 0;  Stat_mat= 0;  Stat_obj= 0

       If (sens_parm < 1 .or. sens_parm > Nsens) Return

       x= SenP(cntr)
       Select case (sens_parm)
         Case (1)
           Parm_bal(1)= x
         Case (2)
           Parm_bal(2)= x
         Case (3)
           Parm_bal(3)= x
         Case (4)
           Parm_bal(4)= Parm_bal(3) * x
       End Select

       If (cntr == 0) then
         Call Out ("Standard Run for parameter",sens_parm, ln= 1)
       Else
         Call Out ("Perturbed Run for parameter",sens_parm, &
                   "with deviation index",cntr, ln= 1)

         Stat1(:,:,1)  = 0;  Stat2(:,:,:,1)= 0;  Stat3(:,1)= 0
         Stat4(:,:,:,1)= 0;  Stat5(:,:,1)  = 0;  N_reg(:,:,1)= 0
       End if

       Call Out ("Parameter value",SenP(cntr))
       
       Call Out ("Ballot parameters",Parm_bal)
       Call Out ("versus standard parameters",Std_sens)

     Else if (Status == 1) then  ! Update data for statistics, using the top set of electeds

       ncs= Size(Clust_set);  top= Elect_dat(1)%lt(1)
       nc= Size(Clust_set(top)%T0,1);  n5= Min(ncs,5)
       np= Size(Elect_dat(1)%ls);  me= Size(Elect_dat)
       ind= Clust_set(top)%l;  ncl= ind - 1
       id= idist+1 - Dist_1

       If (top /= 1) then
         Call Out ("Warning in 'Sensitivity_parm': top cluster from top elected set",top, &
                   "differs from top full objective ordering",1)
       End if
       
       If (cntr == 0) then  ! Standard run stat updates
         Call Out ("Update stats for the Standard Run for parameter",sens_parm, &
                   "for district",idist, ln= 1)

         Allocate (Stand_run(id)%qx(9), Stand_run(id)%ls(np), Stand_run(id)%M0(0:Mvt,2))
         
         Stand_run(id)%n= nc;  Stand_run(id)%p= np;  Stand_run(id)%k= ncl
         top_full= Elect_dat(1)%M1(1,1)

         Stand_run(id)%qx(:7)= Elect_dat(1)%M1(:,1);  Stand_run(id)%qx(8:)= 0
         e2= STV_eq(2,1)                      ! DTV elected index
         Stand_run(id)%qx(8)= Elect_dat(e2)%M1(1,1) / top_full
         e1= STV_eq(1,1);  q1= STV_eq(1,2)    ! STV cluster set 
         i= First_true(Elect_dat(e1)%lt == q1)
         If (i > 0) Stand_run(id)%qx(9)= Elect_dat(e1)%M1(1,i) / top_full

         Stand_run(id)%ls= Elect_dat(1)%ls
         Stand_run(id)%M0= 0;  Stand_run(id)%M0= Elect_dat(1)%M0
         
         Do q= 1,n5
           m= Clust_set(q)%m
           If (m > 0) then
             Stat5(q,0,0)= Stat5(q,0,0) + 1
             Stat5(q,1,0)= Stat5(q,1,0) + m
             Stat5(q,2,0)= Stat5(q,2,0) + m**2
           End if
         End do
         
!        Count # top (:,1) or all (:,2) cluster sets of sizes 1,2,3,4,5
       
         N_reg(ncl,1,0)= N_reg(ncl,1,0) + 1   ! sum for # reg clusters for top cluster set
         ns2= Last_true(Clust_set%l > 2)
         N_reg(0,2,0)= N_reg(0,2,0) + ns2     ! sum # cluster sets with # reg clusters >= 2

         Do j= 1,ncs
             n= Clust_set(j)%l - 1
           N_reg(n,2,0)= N_reg(n,2,0) + 1      ! # cluster sets for each # reg clusters
         End do

       Else  ! Update stats for a perturbed parameter, using the top set of electeds

         Call Out ("Update stats for the Perturbed Run for parameter",sens_parm, "& deviation",cntr, ln= 1)
         Call Out ("for district",idist)
         If (sens_parm == 5) then
           Call Out ("Noise signficance adjustment",Noise_cor)
         End if

         If (All(Elect_dat(1)%ls == Stand_run(id)%ls)) then
           Stat_mat= Stat_mat + 1
         End if

         Allocate (Perturbed(id)%qx(9), Perturbed(id)%ls(np), Perturbed(id)%M0(0:Mvt,2))
         
         Perturbed(id)%n= nc;  Perturbed(id)%p= np;  Perturbed(id)%k= ncl
         top_full= Elect_dat(1)%M1(1,1)

         Perturbed(id)%qx(:7)= Elect_dat(1)%M1(:,1);  Perturbed(id)%qx(8:9)= 0
         e2= STV_eq(2,1)                      ! DTV elected index
         Perturbed(id)%qx(8)= Elect_dat(e2)%M1(1,1) / top_full
         e1= STV_eq(1,1);  q1= STV_eq(1,2)    ! STV cluster set 
         i= First_true(Elect_dat(e1)%lt == q1)
         If (i > 0) Perturbed(id)%qx(9)= Elect_dat(e1)%M1(1,i) / top_full

         Perturbed(id)%ls= Elect_dat(1)%ls
         Perturbed(id)%M0= Elect_dat(1)%M0

         Do i= 1,4
             dif= Abs(Perturbed(id)%qx(i) - Stand_run(id)%qx(i))
           Stat_obj(i,1)= Stat_obj(i,1) + dif 
           Stat_obj(i,2)= Stat_obj(i,2) + dif**2
         End do
         
         Do q= 1,n5
           m= Clust_set(q)%m 
           If (m > 0) then
             Stat5(q,0,1)= Stat5(q,0,1) + 1
             Stat5(q,1,1)= Stat5(q,1,1) + m
             Stat5(q,2,1)= Stat5(q,2,1) + m**2
           End if
         End do
       
!        Count # top (:,1) or all (:,2) cluster sets of sizes 1,2,3,4+
       
         N_reg(ncl,1,1)= N_reg(ncl,1,1) + 1   ! sum # reg clusters for top cluster set
         ns2= Last_true(Clust_set%l > 2)
         N_reg(0,2,1)= N_reg(0,2,1) + ns2     ! sum # cluster sets with # reg clusters >= 2

         Do j= 1,ncs
             n= Clust_set(j)%l - 1
           N_reg(n,2,1)= N_reg(n,2,1) + 1     ! # cluster sets for each # reg clusters
         End do
       End if
       
     Else if (Status == 2) then  ! Finalize statistics   
           
       If (cntr == 0) then  ! Initial run = Standard run
         Call Out ("Final stats for the Standard Run for parameter",sens_parm, ln=2)
         Call Out ("Parameter value",SenP(cntr))

         obd0= -1;  tpd0= -1
         
         Do i= 1,nDst
           obd0(i,1:2)= Stand_run(i)%qx(2:3)
           obd0(i,3:4)= Stand_run(i)%qx(8:9)
           tpd0(i,:,:)= Stand_run(i)%M0
         End do   
         
         Do i= 1,4
           Call Sort (.true.,obd0(:,i))
         End do
         
         Do j= 1,2
           Do i= 0,Mvt
             Call Sort (.true.,tpd0(:,i,j))
           End do
         End do
         
         Do k= 1,7
           x= Perc(k) * nDst;  l= Floor(x)
           
           If (l > 0 .and. l < nDst) then
             If (obd0(l,1) > 0 .and. obd0(l+1,1) > 0) then
               x1= l + 1 - x;  x2= x - l
               Stat1(k,:,0)  = x1 * obd0(l,:)   + x2 * obd0(l+1,:)
               Stat2(k,:,:,0)= x1 * tpd0(l,:,:) + x2 * tpd0(l+1,:,:)
             End if
           Else if (l == 0) then
             If (obd0(1,1) > 0) then
               Stat1(k,:,0)  = obd0(1,:)
               Stat2(k,:,:,0)= tpd0(1,:,:)
             End if
           Else if (obd0(nDst,1) > 0) then
             Stat1(k,:,0)  = obd0(nDst,:)
             Stat2(k,:,:,0)= tpd0(nDst,:,:)
           End if
         End do
         
         Stat3(:,0)= 0

         Do i= 1,4
           j= 2*(i - 1);  k= Count(obd0(:,i) > 0);  If (k < 1) Cycle

           Stat3(j+1,0)= Sum(obd0(:,i), obd0(:,i) > 0) / k
           Stat3(j+2,0)= Sqrt(Sum((obd0(:,i) - Stat3(j+1,0))**2, obd0(:,i) > 0) / k)
           Stat3(j+2,0)= Stat3(j+2,0) / Stat3(j+1,0)
         End do
         
         Do j= 1,2
           Do i= 0,Mvt
             k= Count(tpd0(:,i,j) > 0)
             If (k > 0) then
               Stat4(1,i,j,0)= Sum(tpd0(:,i,j), tpd0(:,i,j) > 0) / k
               Stat4(2,i,j,0)= Sqrt(Sum((tpd0(:,i,j) - Stat4(1,i,j,0))**2, tpd0(:,i,j) > 0) / k)
             End if
           End do
         End do
         
         Do i= 1,5
           n= Nint(Stat5(i,0,0))
           If (n > 0) then
             x= Stat5(i,1,0) / n;  Stat5(i,1,0)= x
             Stat5(i,2,0)= Sqrt(Abs(Stat5(i,2,0) / n - x**2))
           End if    
         End do
         
       Else  !      Perturbed run     
         Call Out ("Final stats for the Perturbed Run for parameter",sens_parm, &
                   "& deviation",cntr, ln= 1)
         Call Out ("Parameter value",SenP(cntr))

         Stat_mat= Stat_mat / nDst
     
         Stat_obj(:,1)= Stat_obj(:,1) / nDst
         Stat_obj(:,2)= Sqrt(Stat_obj(:,2)/(nDst-1) - Stat_obj(:,1)**2)

         obp= -1;  tpp= -1
         Do i= 1,nDst
           obp(i,1:2)= Perturbed(i)%qx(2:3)
           obp(i,3:4)= Perturbed(i)%qx(8:9)
           tpp(i,:,:)= Perturbed(i)%M0
         End do
         
         Do i= 1,4
           Call Sort (.true.,obp(:,i))
         End do
         
         Do j= 1,2
           Do i= 0,Mvt
             Call Sort (.true.,tpp(:,i,j))
           End do
         End do
         
         Do k= 1,7
           x= Perc(k) * nDst;  l= Floor(x)
           
           If (l > 0 .and. l < nDst) then
             If (obp(l,1) > 0 .and. obp(l+1,1) > 0) then  
               x1= l + 1 - x;  x2= x - l
               Stat1(k,:,1)  = x1 * obp(l,:)   + x2 * obp(l+1,:)
               Stat2(k,:,:,1)= x1 * tpp(l,:,:) + x2 * tpp(l+1,:,:)
             End if
           Else if (l == 0) then
             If (obp(1,1) > 0) then  
               Stat1(k,:,1)  = obp(1,:)
               Stat2(k,:,:,1)= tpp(1,:,:)
             End if
           Else if (obp(nDst,1) > 0) then
             Stat1(k,:,1)  = obp(nDst,:)
             Stat2(k,:,:,1)= tpp(nDst,:,:)
           End if
         End do
         
         Stat3(:,1)= 0

         Do i= 1,4
           j= 2*(i - 1);  k= Count(obp(:,i) > 0);  If (k < 1) Cycle

           Stat3(j+1,1)= Sum(obp(:,i), obp(:,i) > 0) / k
           Stat3(j+2,1)= Sqrt(Sum((obp(:,i) - Stat3(j+1,1))**2, obp(:,i) > 0) / k)
           Stat3(j+2,1)= Stat3(j+2,1) / Stat3(j+1,1)
         End do
         
         Do j= 1,2
           Do i= 0,Mvt
             k= Count(tpp(:,i,j) > 0)
             If (k > 0) then
               Stat4(1,i,j,1)= Sum(tpp(:,i,j)) / k
               Stat4(2,i,j,1)= Sqrt(Sum((tpp(:,i,j) - Stat4(1,i,j,1))**2) / k)
             End if
           End do
         End do
         
         Do i= 1,5
           n= Nint(Stat5(i,0,1))
           If (n > 0) then
             x= Stat5(i,1,1) / n;  Stat5(i,1,1)= x
             Stat5(i,2,1)= Sqrt(Abs(Stat5(i,2,1) / n - x**2))
           End if    
         End do
       End if
        
!      Output 
       
       If (cntr == 0) then  ! Standard run
         N_reg(0,:,0)= Sum(N_reg(1:,:,0), 1) ! Sum stats over all # reg clusters

         Do i= 1,2
           frac2p(0,i)= Sum(N_reg(2:,i,0)) / Real(N_reg(0,i,0))

           avg2p(0,i)= 0
           Do n= 2,5
             avg2p(0,i)= avg2p(0,i) + n*N_reg(n,i,0)
           End do
           avg2p(0,i)= avg2p(0,i) / Real(N_reg(0,i,0))
         End do

         Call Out ("Frac of top cluster sets with ncl >= 2", frac2p(0,1), &
                   "Frac of all cluster sets", frac2p(0,2))
         Call Out ("Avg # clusters in top cluster sets with ncl >= 2", avg2p(0,1), &
                   "Avg # clusters of all cluster sets", avg2p(0,2))

         Call Out ("Clustering objective: mean, sigma/mean",Stat3(1:2,0))
         Call Out ("Fitness: mean, sigma/mean",Stat3(3:4,0))
         Call Out ("DTV ratio: mean, sigma/mean",Stat3(5:6,0))
         Call Out ("Converged STV ratio: mean, sigma/mean",Stat3(7:8,0))

         Call Out ("For percentile values", Nint(100 * Perc))
         Call Out ("Clust obj, fitness, DTV ratio, converged STV ratio:")
         Call Out (-1,"Percentiles for 4 objectives", Stat1(:,:,0))
         
         Call Out ("Statistics over the top electeds: how well they match voter rankings")
         Call Out ("'same top': mean # elected of top ranked, then fraction with top>=1 elected, top>=2, etc")
         Call Out ("'some top': mean # elected regardless of rank, then fraction with >=1 elected, >=2, etc")

         Call Out (-1,"Percentiles for the same top ranked elected", Stat2(:,:,1,0))
         Call Out (-1,"Percentiles for some top ranked elected", Stat2(:,:,2,0))
         
         Call Out (-1,"Mean & sigma for # initial to final cluster set matches for their top 5 sets", Stat5(:,:,0))
         
         Call Out ("Standard data by district:")
         Write(8,'(A)') "Full objective & clustering objective & fitness"
         Write(8,'(A)') "Cluster coherent vote & penalty, elected rating & penalty"
         Write(8,'(A)') "DTV ratio & converged STV ratio"
         
         Do id= 1,nDst
           Call Out (" ");  i= (Dist_1-1) + id
           Write (8,'(A,I3,2X, A,3I3)') "District",i, "  Standard nc,np,ncl:", &
                                         Stand_run(id)%n, Stand_run(id)%p, Stand_run(id)%k

           Write (8,'(A)') "   Obj Data: Full objective = clustering objective & fitness"
           Write (8,'(F8.3, 2X,2F8.3)') Stand_run(id)%qx(:3)
           Write (8,'(A)') "   Cluster coherent vote & penalty,  elected rating & penalty"
           Write (8,'(2F8.3, 3X, 2F8.3)') Stand_run(id)%qx(4:7)
           Write (8,'(A)') "   DTV ratio & converged STV ratio"
           Write (8,'(2F8.3)') Stand_run(id)%qx(8:9)
         End do

!        Key statistics

         Call Out ("Key stats for Standard Run for parameter",sens_parm, &
                   "with parameter value",SenP(0), ln=1)

         Call Out ("Frac of top cluster sets with ncl >= 2", frac2p(0,1), ln=1)
         Call Out ("Avg # clusters in top sets with ncl >= 2", avg2p(0,1))

         Out_obj(0,0,1)= 1.0;  Out_obj(1:4,0,1)= Stat3(1:7:2,0)
         Call Out ("Frac match, Clust obj, Fitness, DTV ratio, converged STV ratio", Out_obj(:,0,1))

       Else  !  Perturbed run

         Stat_obj(1:4,3)= Stat3(1:7:2,1) - Stat3(1:7:2,0)
         N_reg(0,:,1)= Sum(N_reg(1:,:,1), 1) ! Sum stats over all # reg clusters

         Do i= 1,2
           frac2p(1,i)= Sum(N_reg(2:,i,1)) / Real(N_reg(0,i,1))

           avg2p(1,i)= 0
           Do n= 2,5
             avg2p(1,i)= avg2p(1,i) + n*N_reg(n,i,1)
           End do
           avg2p(1,i)= avg2p(1,i) / Real(N_reg(0,i,1))
         End do

         frac2p(2,:)= frac2p(1,:) - frac2p(0,:)
         Call Out ("Std:Per:dif. Frac of top cluster sets with ncl >= 2", frac2p(:,1))
         Call Out ("Std:Per:dif. Frac of all cluster sets with ncl >= 2", frac2p(:,2))

         avg2p(2,:)= avg2p(1,:) - frac2p(0,:)
         Call Out ("Std:Per:dif. Avg # top cluster sets with ncl >= 2", avg2p(:,1))
         Call Out ("Std:Per:dif. Avg # all cluster sets with ncl >= 2", avg2p(:,2))

         Call Out (-1,"Mean & std of |pert - stand| for 4 obj's, + final pert-std means", Stat_obj) 

         Call Out ("The 4: Clustering obj, fitness, DTV ratio, converged STV ratio:")

         Call Out ("Perturbed:   Clustering objective: mean, sigma/mean",Stat3(1:2,1))
         Call Out ("Vs standard: Clustering objective: mean, sigma/mean",Stat3(1:2,0))

         Call Out ("Perturbed:   Fitness: mean, sigma/mean",Stat3(3:4,1))
         Call Out ("Vs standard: Fitness: mean, sigma/mean",Stat3(3:4,0))

         Call Out ("Perturbed:   DTV ratio: mean, sigma/mean",Stat3(5:6,1))
         Call Out ("Vs standard: DTV ratio: mean, sigma/mean",Stat3(5:6,0))

         Call Out ("Perturbed:   Converged STV ratio: mean, sigma/mean",Stat3(7:8,1))
         Call Out ("Vs standard: Converged STV ratio: mean, sigma/mean",Stat3(7:8,0))

         Call Out ("Statistics over top sets of electeds: how well they match voter rankings")
         Call Out ("'same top': mean # elected of top ranked, then fraction with top >=1 elected, top >=2, etc")
         Call Out ("'some top': mean # elected regardless of rank, then fraction with >=1 elected, >=2, etc")

         Call Out (-1,"Overall perturbed mean & sigma for the same top ranked elected", Stat4(:,:,1,1))
         Call Out (-1,"Versus standard mean & sigma for the same top ranked elected", Stat4(:,:,1,0))
         
         Call Out (-1,"Overall perturbed mean & sigma for some top ranked elected", Stat4(:,:,2,1))
         Call Out (-1,"Versus standard mean & sigma for some top ranked elected", Stat4(:,:,2,0))
         
         Call Out ("For percentile values", Nint(100 * Perc))
         Call Out ("Clustering obj, fitness, DTV ratio, converged STV ratio:")

         Call Out (-1,"Percentiles for 4 perturbed objectives", Stat1(:,:,1))
           Call Out (-1,"versus the standard values", Stat1(:,:,0))
         Call Out (-1,"Percentiles for the same top ranked elected", Stat2(:,:,1,1))
           Call Out (-1,"versus the standard values", Stat2(:,:,1,0))
         Call Out (-1,"Percentiles for some top ranked elected", Stat2(:,:,2,1))
           Call Out (-1,"versus the standard values", Stat2(:,:,2,0))
         
         Call Out (-1,"Mean & sigma for # initial to final cluster matches for top 5 sets", Stat5(:,:,1))
         
         Call Out ("Perturbed vs standard data by district:")
         
         Do id= 1,nDst
           Call Out (" ");  i= Dist_1-1 + id
           Write (8,'(A,I3,4X, A,3I3)') "District",i, "Perturbed nc,np,ncl:", &
                                         Perturbed(id)%n, Perturbed(id)%p, Perturbed(id)%k
           Write (8,'(3X,A,I3)') "vs Standard # regular clusters:", Stand_run(id)%k
           
           Write (8,'(3X,A,6I3)') "Standard  electeds:", Stand_run(id)%ls
           Write (8,'(3X,A,6I3)') "Perturbed electeds:", Perturbed(id)%ls
           
           Write (8,'(A)') "       Standard vs Perturbed:"
           Write (8,'(A)') "   Full objective & clustering objective & fitness"
           Write (8,'(F8.3,1X,2F8.3, 3X, F8.3,1X,2F8.3)') Stand_run(id)%qx(:3), Perturbed(id)%qx(:3)
           Write (8,'(A)') "   Cluster coherent vote & penalty"
           Write (8,'(2F8.3, 3X, 2F8.3)') Stand_run(id)%qx(4:5), Perturbed(id)%qx(4:5)
           Write (8,'(A)') "   Elected rating & penalty"
           Write (8,'(2F8.3, 3X, 2F8.3)') Stand_run(id)%qx(6:7), Perturbed(id)%qx(6:7)
           Write (8,'(A)') "   DTV ratio & converged STV ratio"
           Write (8,'(2F8.3, 3X, 2F8.3)') Stand_run(id)%qx(8:9), Perturbed(id)%qx(8:9)
         End do

!        Key statistics

         Call Out ("Key stats for Perturbed Run for parameter",sens_parm, &
                   "deviation index",cntr, ln= 1)
         Call Out ("with parameter value",SenP(cntr))
 
         Call Out ("Frac of top cluster sets with ncl >= 2", frac2p(0,1), ln=1)
         Call Out ("Avg # clusters in top sets with ncl >= 2", avg2p(0,1))

         Out_obj(0,cntr,1)= Stat_mat;  Out_obj(1:4,cntr,1)= Stat3(1:7:2,1)
         Call Out ("Frac match, Clust obj, Fitness, DTV ratio, converged STV ratio", Out_obj(:,cntr,1))

         Out_obj(:,cntr,2)= Out_obj(:,cntr,1) - Out_obj(:,0,1)
         Call Out ("Dif obj", Out_obj(:,cntr,2))
       End if
     End if
     
   End Subroutine Sensitivity_parm 


   Subroutine General_stats (Status,Dist_1,nDst, idist,np,ncs, Non_clust, &
                             Init, Clust_set,Clust_elect,Elect_dat)
   
!    Initialize, update, or finalize statistics, across the districts. Include the performance 
!    of the non-clustering methods in relation to the clustering method, also STV and DTV..
   
     Integer,           Intent(in) :: Status          ! Initialize the statistics if = 0, update if = 1, finalize if = 2 
     Integer,           Intent(in) :: Dist_1          ! First district index
     Integer,           Intent(in) :: nDst            ! Total # districts used for the statistics

     Integer,           Intent(in) :: idist           ! Current district index
     Integer,           Intent(in) :: np              ! # candidates to be elected
     Integer,           Intent(in) :: ncs             ! # cluster sets in 'Clust_set'

     Type(Multi_listR), Optional, Intent(in) :: Non_clust  ! Data for 'nMt' selected non-clustering methods 'm'
                                     !                       to elect multiple candidates.
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

     Type(Multi_listR), Optional, Intent(in) :: Init(:)  ! (N_init) Initial cluster set data
                                                ! %k = cluster set converged to, 0 for convergence failure
                                                !      -1 for severe error
                                                ! %n = ncl = # regular clusters

                                                ! %Q0(np,0:2,4) Elected candidates, by 4 methods 
                                                !      (:,:,:1) = Rem_frac electeds (original for STV/DTV)
                                                !      (:,:,:2) = DHondt electeds
                                                !      (:,:,:3) = Optimal electeds 
                                                !      (:,:,:4) = Converged optimal electeds
                                                !     (:,0,:) = electeds in increasing order
                                                !     (:,1,:) = electeds in ranked order
                                                !     (:,2,:) = corresponding clusters

                                                ! %M1(7,4)     Objective data for the electeds $Q0(:,0,:)
                                                ! %M2(ncl,ncl) Cluster correlation matrix, 0 on diagonal

                                                ! %T0(0:nc,0:ind,2) (1:,1:,1) = Candidate mean ratings  
                                                !      for each cluster, with cluster size at (0,1:,1)
                                                !      (1:,1:,2) = Cluster portions with averaged 
                                                !                  noise zeroed ratings at (1:,0,2)

     Type(Multi_listR), Optional, Intent(in) :: Clust_set(:) ! (ncs) Final cluster sets
                                                ! %k   = best initial cluster set that converged to this cluster set
                                                ! %l   = ind = ncl + 1 = # clusters, including independents
                                                ! %m   = # initial cluster sets that converged to this cluster set
                                                ! %ls(m)= list of those initial sets 
                                                ! %fux = clustering objective ratio = (q)%rx(1) / (1)%rx(1) 

                                                ! %lt(ncl) Mapping of initial to converged regular clusters 

                                                ! %px(4):  Objective data: min size, max correlation, 
                                                !            independents size (frac), final residual
                                                ! %qx(4):  Objective factors: min size penalty, max correlation penalty, 
                                                !            independents size penalty, residual penalty
                                                ! %rx(3):  Objective value(1) = coherent vote(2) * penalty(3)
;
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

                                                ! %L0(0:4,0:iv) Integer convergence data per update call
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

                                                !   %L1(-2:nc,0:ind) Candidate ordering data for
                                                !                    noise zeroed ratings
                                                !     (-2:0,:) = # top, significant, viable cand's
                                                !     (1:,1:)  = Decreasing ordering from 
                                                !                zeroed ratings = %T0(:,1:,2)
                                                !     (:,0)   = Ordering & data from %T0(1:,0,1)

                                                !   %T0(nc,0:ind,3) Cluster point data
                                                !      (:,1:,1) = Cluster mean vectors
                                                !      (:,1:,2) = Cluster mean vectors with noise zeroing 
                                                !      (:,1:,3) = Cluster portions for each candidate
                                                !      (:,0,1) = Scaled Zrate0 (for ordering the %L1(1:,0))
                                                !      (:,0,2) = Zrate0 = cluster averaged noise zeroed mean ratings

     Type(Multi_listR), Optional, Intent(in) :: Clust_elect(:) ! (ncs) Cluster set data for best sets of electeds
                                         ! %k = index of the best set of electeds %Q0(:,0,1) in Elect_dat
                                         ! %l = 'ind' = # clusters, including independents
                                         ! %m = 'nf'= # sets of possible electeds
                                         ! %fsx = full objective for the top set of electeds 
                                         !        %Q0(:,0,1) = Clust(q)%sx(1)
                                         ! %fux = full objective ratio = (q)%fsx / (1)%fsx

                                         ! %sx(7)  Full objective data for %k, or from Clust_set%rx & %M2 
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
     
                                         ! %M2(3,nf)= Electeds fitness (1) = average electeds rating (2) * 
                                         !            cluster size deviation penalty (3) 

                                         ! %T0(ind,2,nf) Objective data for each set of possible electeds
                                         !   1 = Cluster size deviations = electeds size - true
                                         !   2 = Cluster size deviation penalties

     Type(Multi_listR), Optional, Intent(in) :: Elect_dat(:) ! (me) Data for possible sets 'e' of electeds,
                                         ! ordered by 'best' full objective %M1(1,1)
                                         ! %k = The 'best' cluster set 'q' for 'e' = %lt(1)
                                         ! %l = The matching 'eq' in Clust_elect(q)%Q0 = %vl(1)
                                         ! %m = nq = # cluster sets that list 'e' 
                                         !      or = 1 for lower cluster sets
                                         ! %n = # regular clusters for best cluster set 'q'
                                         ! %q = Unique non-clustering index if 'e' is a
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
                                         ! %L0(np,nq) These elected sets are  in preferential order
                                         !            for the cluster sets ordered as in %lt
     
                                         ! %M0(0:Mvt,2): Measures of success
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

!  Local:
     Real,    Parameter :: eps= 0.00001  ! Non-zero quantity
     Integer, Parameter :: Msc= 150      ! Max # convergence cases to be recorded: all cluster sets 
                                         !  with ncl >= 2, + the best cluster set if if has ncl = 1
     Integer, Parameter :: Mp= 5         ! Max value of 'np' for statistica purposes: 
                                         !   Combine this value and above.
     Character(1) :: c1(9)= (/"1","2","3","4","5","6","7","8","9"/)

     Integer, Save :: Elc1(0:Mp,mx_Dist)  ! Top elected set for each district, (0,:)= np
     Integer, Save :: Elc2(Mp,Msc)        ! Top elected set for each cluster set with ncl >= 2, 
                                          ! or also for the overall top elected set if its
                                          ! cluster set has ncl = 1

     Integer, Parameter :: NonR= N_init-N_cand + 2  ! # non-random initial cluster sets + 2 for 
                                                    ! STV & DTC with their original electeds added.
     Integer, Save :: nonRnd_elc(0:Mp,NonR,2,mx_Dist)  ! Optimal elected set (1:np,..) for each non-random 
                                                       !  initial cluster set(1), converged set(2)
                                                       !  where (0,iq,...) = cluster set converged to
                                                       !  from initial set 'iq'
     Integer, Save :: nonRnd_ratio(NonR,mx_Dist)    ! Ratio * 1000 of each non-random converged cluster set full objective
                                                    !   for their top elected sets, to the optimal full objective.
     Integer       :: mat_nonEl(NonR,3)   ! Count # matches for non-random initial sets of elected sets
                                          ! (:,1) = # matches of initial top elected to converged top elected
                                          ! (:,2) = # matches of initial top elected to overall top elected
                                          ! (:,3) = # matches of converged top elected to overall top elected
     Integer       :: mat_nonOb(NonR)     ! Count # matches for non-random converged to optimal full objectives
     Integer       :: mat_id(nDst,NonR,2) ! List of district indices for each match for each method and elected set

     Integer       :: m2              ! # current cluster sets with # reg clusters >= 2
     Integer       :: ncl             ! # regular clusters in a cluster set
     Integer       :: ncyc            ! Index of the # convergence cycles
     Integer, Save :: ns              ! # cluster sets, summed over districts, with ncl >= 2
                                      !   or with ncl = 1 if it's the top cluster set
     Integer, Save :: cnv_dat(6,Msc)  ! Convergence data
                                      ! 1 = district index
                                      ! 2 = # regular clusters
                                      ! 3 = # convergence cycles
                                      ! 4 = # cluster merge and deletion operations
                                      ! 5 = # centroiding iterations
                                      ! 6 = # function calls
     Real,    Save :: obj_dat(4,Msc)  ! Objective data for converged cluster sets from Clust_set%px
     Real,    Save :: obj_cmp(3,Msc)  ! Objective components from Clust_set%px
                                      ! overall objective(1)= coherent vote(2) * overall penalty(3)
     Integer, Save :: Nelc(2:5)       ! # districts with each possible # elected, 
                                      !   2 thru 5 (or more)

     Integer, Save :: MismatchI(5,20) ! (1) district #, (2) best cluster set, (3) most probable set
                                      ! (4) # clusters (best), (5) # clusters (top)
     Integer, Save :: nMis            ! # mismatches = when best /= top
     Real,    Save :: MismatchR(8,20) ! (1) top full obj, (2) top clust obj, (3) top elc fitness 
                                      ! (4) best full obj, (5) best clust obj, (6) best elc fitness 
                                      ! For  q = Elect_dat(1)%k = best cluster set:
                                      ! (7) clustering objective ratio Clust_set(q)%fux
                                      ! (8) full objective ratio Clust_elect(q)%fux
     Real,    Save :: Frac_top(4,2,nMt) ! Fraction of the districts where the top 'n' <= 4 
                                        !   non-clustering electeds 'nMt' match the clustering 
                                        !   electeds as a set (:,1,:) or in rank order (:,2,:) 
     Integer, Save :: SD_cnt(2)         ! (1) = # districts where STV & DTV electeds differ
                                        ! (2) = # districts where STV & DTV converged cluster sets differ
     Integer, Save :: STV_dif(5,5)      ! # times STV elected sets differ. Initial cluster set (Init(19)):  
                                        ! 1 = original elected set, 2:3 = DHondt & optimal from Init(19) 
                                        ! Converged cluster set (q1 = Init(19)%k): 4 = DHondt, 5 = optimal 

     Integer, Save :: ntop              ! Counter for clusters from cluster sets with >= 2 regular clusters
     Real,    Save :: Top_ratios(3,2)   ! For accumlating means (1) and variances (2) for the 3 ratios 
                                        ! of top candidates from Clust_set%L1: n1/np, n2/np, n3/np

     Integer, Parameter :: ni= 5        ! # intervals - 1 in the histogram
     Real,   Save :: STV_hist(0:ni,0:2) ! Histogram of STV best full objective to clustering ratios. 
                                        !   0 = starting value of each interval
                                        !   1 = count STV to clustering ratios => output fraction
                                        !   2 = count DTV to clustering ratios => output fraction
     Real,   Save :: Cnv_hist(0:ni,0:5) ! Histogram for convergence data (from cnv_dat)
                                        !   the final column = max convergence #s for rows 2:6
                                        !   0 = fractional interval starting values, which scale the final column
                                        !   1 = # regular clusters
                                        !   2 = # convergence calls
                                        !   3 = # cluster merge and deletion operations
                                        !   4 = # centroiding iterations
                                        !   5 = # function calls

     Integer, Save :: ncl_STV(3,mx_Dist)        ! # regular clusters for each STV
                                                !  cluster set: initial, extended, converged
     Integer, Save :: ncl_DTV(3,mx_Dist)        ! # regular clusters for each DTV
                                                !  cluster set: initial, extended, converged
     Real,    Save :: obj_STV(3,mx_Dist)        !  STV electeds full objective for the 
                                                !  cluster sets: extended, converged, best (over all cluster sets)
     Real,    Save :: obj_DTV(3,mx_Dist)        ! DTV electeds full objective for the 
                                                !  cluster sets: extended, converged, best (over all cluster sets)
     Real,    Save :: Non_ratio(-1:mx_Dist,nMt) ! Ratio of each non-clustering best full
                                                ! objective to the best clustering full objective
                                                ! (-1,:) = mean ratio overall districts
                                                ! (0,:)  = mean ratio overall districts with ratio < 1
     Real,    Save :: SD_ratio(-1:mx_Dist,2)    !  (:,1) Ratio of the converged STV (Init(19)%k)
                                                !    full objective to the best clustering full objective
                                                !  (:,2) Same ratio for the converged (Init(20)%k) 
                                                !    DTV full objective

     Real, Save :: DtoS_ratio(mx_Dist,2)        ! Ratios of DTV to STV full objectives:
                                                ! (1) Using DTV & STV best cluster sets, with dif electeds
                                                ! (2) Using DTV & STV converged cluster sets, with dif objective 

     Type(Multi_listR), Save :: cls_dat(mx_Dist)  ! Data for each cluster set 
                                                  ! %k= best cluster set 'q'
                                                  ! %q= # cluster sets = 'ncs'
                                                  ! %m= total # clusters = Sum(ls)
                                                  ! %ls(ncs)= # regular clusters 'ncl' in each cluster set
                                                  ! %M1(5,m) = size, dens, bound, dens pen, bound pen (by cluster)

     Real,    Allocatable :: db(:,:) ! (ns,4) Density & boundary data for cluster sets
     Integer, Allocatable :: ky(:)   ! (ns)

     Integer :: STVset(np,5)        ! STV elected sets. Initial cluster set (Init(19)):  
                                    ! 1 = original, 2 = DHondt, 3 = optimal 
                                    ! Converged cluster set (q1 = Init(19)%k): 4 = DHondt, 5 = optimal 
     Real    :: Exact_match(nMt)    ! Fraction of districts for which a non-clustering 
                                    ! method matches clustering exactly
     Real    :: Good_match(nMt)     ! Fraction of districts for which a non-clustering 
                                    ! method either matches clustering exactly or  
                                    ! is very close in full objective
     Real    :: Exact_converged(2)  ! Fraction of districts for which STV & DTV actuals
                                    ! match clustering exactly
     Real    :: Good_converged(2)   ! Fraction of districts for which STV & DTV actuals
                                    ! either match clustering exactly or are very close 
                                    ! in full objective
     Integer, Save :: ncf, conv_fail(Msc)  ! List of districts with a convergence failure

     Integer, Save :: ncnvg(6)      ! 1 = # partially converged cluster sets, without damping
                                    ! 2 = # fully converged cluster sets
                                    ! 3 = # with damping at 1 or more phases
                                    ! 4 = # with damping in restart phase 1
                                    ! 5 = # with damping in restart phase 2
                                    ! 6 = # with damping in restart phase 3

     Integer :: DSTV_dif(2), lsi(N_init)
     Integer :: rnk_cl(np), rnk_non(np), set_cl(4,4), set_non(4,4)
     Integer :: phs(3), Nel(4), lq2(ncs), nl(ncs), ndif(nMt), ic(5), jc(5), Lsdt(nDst,5)
     Real    :: y(3), tmp(Msc), frac(nMt), mean_SDTV(2), cnvM(5)  

     Logical :: Msk(nMt), Msk0(nDst), Msk1(nDst)

     Integer :: non(4,nMt), non_ls(nDst,nMt), non_ky(nDst,nMt)
     Integer :: nona(4,2), non_lsa(nDst,2), non_kya(nDst,2)
     Real    :: non_rd(nDst,nMT), non_rda(nDst,2)

     Real    :: mean_STV1, mean_DTV1, mean0_STV1, mean0_DTV1, med_STV1, med_DTV1
     Real    :: mean_STV2, mean_DTV2, mean0_STV2, mean0_DTV2, med_STV2, med_DTV2
     Real    :: x1, x2, dif, fac, max1, mean2, sig1, sig2
     
     Integer :: STV_ind, DTV_ind, Borda_ind
     Integer :: id1, ios, np5, code, nonC, nmd, top_ncl, nitr, nline, best
     Integer :: ml, mn, n0, n1, n2, n3, n4, ne, nq, nt, q0, q1, q2, qj
     Integer :: cl, e1, e2, ej, eq, i2, id, iq, j1, j2, k0, k1, l1, l2
     Integer :: e, i, j, k, l, m, n, q, p
     
     If (Status == 0) then  ! Initialize
       Elc1= 0;  Elc2= 0;  nonRnd_elc= 0;  nonRnd_ratio= 0
       ncl_STV= 0;  ncl_DTV= 0;  obj_STV= 0;  obj_DTV= 0
       SD_ratio= 0;  DtoS_ratio= 0;  Non_ratio= 0
       Exact_match= 0;      Good_match= 0
       Exact_converged= 0;  Good_converged= 0

       ns= 0;  cnv_dat= 0;  obj_dat= 0;  obj_cmp= 0
       Nel= 0;  Nelc= 0;  Frac_top= 0;  ntop= 0;  Top_ratios= 0
       nMis= 0;  MismatchI= 0;  MismatchR= 0

       SD_cnt= 0;  STV_dif= 0;  DSTV_dif= 0;  mean_SDTV= 1

       STV_hist(:,1:)= 0;  Cnv_hist(:,1:)= 0
       STV_hist(:,0)= (/0.0, 0.50, 0.75, 0.90, 0.95, 0.9999/)
       Cnv_hist(:,0)= (/0.0, 0.20, 0.40, 0.60, 0.80, -1.0/)
       
       ntop= 0;  top_ratios= 0;  ncnvg= 0;  ncf= 0;  conv_fail= 0

     Else if (Status == 1) then  ! Update
       Call Out ("Enter 'General_stats' update for district",idist, ln=1)

       Call List_of_true (Init%k <= 0, k,lsi)
       If (k > 0) then
         n1= ncf+1;  ncf= ncf + k;  conv_fail(n1:ncf)= idist
         Call Out ("For district #",idist, ln=1)
         Call Out ("Convergence failures from initial sets",lsi(:k))
       End if

       id= idist - Dist_1 + 1;  n0= Min(np,Mp);  fac= 1000.0 / Elect_dat(1)%M1(1,1)
       Elc1(0,id)= n0;  Elc1(1:n0,id)= Elect_dat(1)%ls(:n0)

       Do i= 1,NonR-2
         iq= N_cand + i;  q= Init(iq)%k;  nonRnd_elc(0,i,:,id)= q
         nonRnd_elc(1:n0,i,1,id)= Init(iq)%Q0(:n0,0,3)             ! initial set optimal electeds

         If (q > 0) then
           nonRnd_elc(1:n0,i,2,id)= Init(iq)%Q0(:n0,0,4)           ! Clust_elect(q)%Q0(:,0,1) = converged optimal
           nonRnd_ratio(i,id)   = Nint(fac * Init(iq)%M1(1,4))   ! Clust_elect(q)%sx(1) 
         End if
       End do

       j= 0
       Do i= NonR-1,NonR
         i2= i - 2;  iq= N_cand + i2;  q= Init(iq)%k;  j= j + 1
         nonRnd_elc(0,i,:,id)= q
         nonRnd_elc(1:n0,i,1,id)= Init(iq)%Q0(:n0,0,1)             ! original STV or DTV electeds
         nonRnd_elc(1:n0,i,2,id)= Init(iq)%Q0(:n0,0,1)             ! original STV or DTV electeds

         If (q > 0) then
           ej= Non_clust%L1(j,1);  qj= First_true(Elect_dat(ej)%lt == q)
           If (qj > 0) then
             nonRnd_ratio(i,id)= Nint(fac * Elect_dat(ej)%M1(1,qj))  ! full objective for the converged cluster set from
                                                                     !   extended STV or DTV, using original STV or DTV electeds 
           End if
         End if
       End do

!      Count convergence 'Success' codes for each converged cluster set

       Do q= 1,ncs

!        final success code: partial or full convergence

         i= clust_set(q)%vl(2);  If (i >= 1) ncnvg(i)= ncnvg(i) + 1

!        Record damping by restart phase

         Do j= 1,3
           phs(j)= Count(clust_set(q)%L0(0,1:) == j .and. clust_set(q)%L0(1,1:) == 0)

           If (phs(j) > 0) ncnvg(3+j)= ncnvg(3+j) + 1
         End do

         If (Any(phs(:3) > 0)) then
           ncnvg(3)= ncnvg(3) + 1
         End if
       End do

!      STV and DTV stat updates

       STV_ind= 1;  DTV_ind= Non_clust%vl(2)
       q1= Init(19)%k;  q2= Init(20)%k
       e1= Non_clust%L1(1,1);  e2= Non_clust%L1(2,1)

       If (q1 > 0) then
         ncl_STV(1,id)= np
         ncl_STV(2,id)= Init(19)%n 
         ncl_STV(3,id)= Clust_set(q1)%l-1

         obj_STV(1,id)= Init(19)%M1(1,1)
           j1= First_true(Elect_dat(e1)%lt == q1)
         obj_STV(2,id)= Elect_dat(e1)%M1(1,j1)
         obj_STV(3,id)= Elect_dat(e1)%M1(1,1)

         SD_ratio(id,1)= obj_STV(2,id)/ Elect_dat(1)%M1(1,1)

         STVset(:,:3)= Init(19)%Q0(:,0,:3)
         STVset(:,4) = Clust_elect(q1)%L1(:,0)
         STVset(:,5) = Clust_elect(q1)%Q0(:,0,1)

         Do i= 1,4
           Do j= i+1,5
             If (Any(STVset(:,i) /= STVset(:,j))) STV_dif(i,j)= STV_dif(i,j) + 1
           End do
         End do
       Else
         Call Out ("Warning in 'General_stats': Bad STV cluster set",q1)
         STVset(:,:3)= Init(19)%Q0(:,0,:3)

         Do i= 1,2
           Do j= i+1,3
             If (Any(STVset(:,i) /= STVset(:,j))) STV_dif(i,j)= STV_dif(i,j) + 1
           End do
         End do
       End if

       If (q2 > 0) then
         ncl_DTV(1,id)= np;  ncl_DTV(2,id)= Init(20)%n 
         ncl_DTV(3,id)= Clust_set(q2)%l-1
         
         obj_DTV(1,id)= Init(20)%M1(1,1)
           j2= First_true(Elect_dat(e2)%lt == q2)
         obj_DTV(2,id)= Elect_dat(e2)%M1(1,j2)
         obj_DTV(3,id)= Elect_dat(e2)%M1(1,1)
         SD_ratio(id,2)= obj_DTV(2,id)/ Elect_dat(1)%M1(1,1)
       Else
         Call Out ("Warning in 'General_stats': Bad DTV cluster set",q2)
       End if

!      Objective ratios, noting that same elected set +
!      same converged cluster set => same full objective.
!      The 'best' full objective for an elected set 'e'
!      comes from the best cluster 'q' for 'e' = Elect_dat(e)%k,
!      which is the cluster set yielding the highest full
!      objective.

!      Compute the best full objective ratio for the elected set 
!      of each non-clustering method vs the top elected set ('Non_ratio')

       Non_ratio(id,:)= Elect_dat(Non_clust%L1(:,1))%fux

!      Record the 'DtoS_ratio' of the DTV best full objective to the 
!      STV best full objective, when their elected sets differ (1). 
!      Also on the ratio for the converged DTV & STV cluster sets,
!      if their elected sets are the same but not their converged
!      cluster sets (2).

       If (All(obj_STV(:,id) > 0) .and. All(obj_DTV(:,id) > 0)) then
         DtoS_ratio(id,1)= obj_DTV(3,id) / obj_STV(3,id)
         DtoS_ratio(id,2)= obj_DTV(2,id) / obj_STV(2,id)
       End if

       If (STV_ind /= DTV_ind) SD_cnt(1)= SD_cnt(1) + 1
       If (q1 /= q2) SD_cnt(2)= SD_cnt(2) + 1

       If (pr_out > 1) then
         Call Out ("For district",idist, ln=1)
         Call Out (-1,"STV & DTV elected indices & top cluster sets",Non_clust%L1(:2,:)) 

         Call Out ("STV best full objective ratio to best clustering",Non_ratio(id,1), &
                   "DTV best full objective ratio to best clustering",Non_ratio(id,2), ln=1)
         Call Out ("STV converged full objective ratio to best clustering",SD_ratio(id,1), &
                   "DTV converged full objective ratio to best clustering",SD_ratio(id,2))
         Call Out ("DTV to STV best full objective ratio if electeds differ",              &
                    DtoS_ratio(id,1), "DTV to STV ratio if converged cluster sets differ", &
                    DtoS_ratio(id,2))
       End if

!      Count # districts by # positions to be elected

       np5= Min(np,5);  Nelc(np5)= Nelc(np5) + 1

!      Count # districts by 'best' vs 'top' cluster set: 
!      top = 1. best = Elect_dat(1)%k, which has the top 
!      full objective subject to clustering objective >= 
!      0.90 * top objective.

       nt= Clust_set(1)%l - 1;  best= Elect_dat(1)%k;  n2= Clust_set(best)%l - 1

       If (best /= 1 .and. nMis < 20) then
         nMis= nMis + 1

         MismatchI(1,nMis)= idist;  MismatchI(2,nMis)= best
         MismatchI(3,nMis)= Maxloc(Clust_set%m,1)
         MismatchI(4,nMis)= n2;  MismatchI(5,nMis)= nt

         Forall(i=1:3) 
           MismatchR(i,nMis)  = Clust_elect(1)%sx(i) 
           MismatchR(3+i,nMis)= Clust_elect(best)%sx(i)
         End forall

           q= Elect_dat(1)%k
         MismatchR(7,nMis)= Clust_set(q)%fux
         MismatchR(8,nMis)= Clust_elect(q)%fux
       End if

!      Count all cluster sets by # reg clusters

       nl= Clust_set%l - 1
       Call List_of_true (nl > 1, m2,lq2)  ! Cluster sets with # ncl >= 2

       q= Elect_dat(1)%k  

       If (Clust_set(q)%l == 2) then  ! Add on the best cluster set if ncl = 1
         n= m2 + 1;  lq2(2:n)= lq2(:m2);  lq2(1)= q;  m2= n
       End if

       Do q0= 1,m2
         q= lq2(q0);  ncl= nl(q);  ncyc= Clust_set(q)%vl(4)
         nmd = Clust_set(q)%vl(3)
         nitr= Clust_set(q)%vl(5);  nline= Clust_set(q)%vl(6)
         
         If (ns < Msc) then     ! Record convergence data for ncl >= 2
           ns= ns + 1           ! convergence counter

           Elc2(:n0,ns)= Clust_elect(q)%Q0(:n0,0,1)

           cnv_dat(1,ns)= idist  ! district index
           cnv_dat(2,ns)= ncl   ! final # regular clusters
           cnv_dat(3,ns)= ncyc  ! # convergence cycles
           cnv_dat(4,ns)= nmd   ! # mergers and deletions
           cnv_dat(5,ns)= nitr  ! total # centroiding iterations
           cnv_dat(6,ns)= nline ! total # line search points
           
           obj_dat(:,ns)= Clust_set(q)%px  ! objective factors
           obj_cmp(:,ns)= Clust_set(q)%rx  ! objective components
         End if
         
         ntop= ntop + ncl  ! accumulate # clusters over all cluster sets with ncl >= 2

!        Stats for # top, significant, viable candidates 
!        for each clusteras fractions of # elected

         Do cl= 1,ncl
           y= Clust_set(q)%L1(-2:0,cl) / Real(np)  
           Top_ratios(:,1)= Top_ratios(:,1) + y
           Top_ratios(:,2)= Top_ratios(:,2) + y**2
         End do
       End do
         
!      Fraction of the districts where the top non-clustering electeds match
!      the clustering electeds either as a set or also with the same ranking

       set_cl= -1;  set_non= -1;  rnk_non= -1
       n4= Min(np,4);  Nel(:n4)= Nel(:n4) + 1

       rnk_cl= Elect_dat(1)%L0(:,1);  set_cl(1,1)= rnk_cl(1)

       Do k= 2,n4
         set_cl(:k,k)= rnk_cl(:k);  Call Sort (.true.,set_cl(:k,k))
       End do

!      Set matching vs rank matching of non-clustering electeds to the clustering electeds

       Method_loop2 : Do m= 1,nMt
         rnk_non= Non_clust%Q0(:np,1,m);  set_non(1,1)= rnk_non(1)

         Do k= 2,n4
           set_non(:k,k)= rnk_non(:k);  Call Sort (.true.,set_non(:k,k))
         End do

         Do k= 1,n4
           If (All(set_non(:k,k) == set_cl(:k,k))) Frac_top(k,1,m)= Frac_top(k,1,m) + 1
           If (All(rnk_non(:k)   == rnk_cl(:k)))   Frac_top(k,2,m)= Frac_top(k,2,m) + 1
         End do
       End do Method_loop2

!      Store data in cls_dat

       cls_dat(id)%k= Elect_dat(1)%k;  cls_dat(id)%q= ncs
       ml= Sum(nl);  cls_dat(id)%m= ml

       Allocate(cls_dat(id)%ls(ncs), cls_dat(id)%M1(5,ml))
       
       i= 0
       Do q= 1,ncs
         ncl= nl(q);  cls_dat(id)%ls(q)= ncl

         Do cl= 1,ncl
           i= i + 1;  cls_dat(id)%M1(1,i)= Clust_set(q)%M3(1,cl)
           cls_dat(id)%M1(2:3,i)= Clust_set(q)%M3(4:5,cl)
           cls_dat(id)%M1(4:5,i)= Clust_set(q)%M3(2:3,cl)
         End do
       End do

     Else  ! Finalize and output statistics

       Call Out ("Enter 'General_stats' end-of-run output")

!      Convergence history histogram
       
       If (ns > 0) then
         Do j= 1,5
           cnvM(j)= Maxval(cnv_dat(1+j,:ns)) 
           tmp(:ns)= cnv_dat(1+j,:ns)/(cnvM(j) + 0.0001)
         
           Do k= 1,ns
             i= Last_true(tmp(k) >= Cnv_hist(:4,0)) - 1
             Cnv_hist(i,j)= Cnv_hist(i,j) + 1
           End do
         End do
       
         Cnv_hist(:4,1:)= Cnv_hist(:4,1:) / ns
         Cnv_hist(5,1:) = cnvM
       End if
       
!      # districts electing 'n' or more candidates
       
       Do i= 1,4
         n= Nel(i);  If (n > 0) Frac_top(i,:,:)= Frac_top(i,:,:) / n
       End do

       If (ntop > 0) then
         Top_ratios= Top_ratios / ntop
         Top_ratios(:,2)= Sqrt(Top_ratios(:,2) - Top_ratios(:,1)**2)
       End if
       

       Call Out (" ")
       Call Out ("Basic statistics")

       Call Out ("# districts for each # positions to be elected (2 - 5)",Nelc)

       jc= 0
       District_loop1 : Do id= 1,nDst
         Do q= 1,cls_dat(id)%q
           ncl= Min(cls_dat(id)%ls(q), 5)
           jc(ncl)= jc(ncl) + 1
         End do
       End do District_loop1

       ic= 0
       District_loop2 : Do id= 1,nDst
         q= cls_dat(id)%k;  ncl= Min(cls_dat(id)%ls(q), 5)
         i= (Dist_1-1) + id

         n= ic(ncl) + 1;  Lsdt(n,ncl)= i;  ic(ncl)= n
       End do District_loop2

       Call Out ("# cluster sets over all districts, counted by # clusters in each set",jc)
       n= Sum(jc(2:));  Call Out ("# cluster sets with 2 or more regular clusters",n)

       If (ncf > 0) Call Out ("List of districts with convergences that failed",conv_fail(:ncf))
       Call Out ("# partly converged cluster sets, # fully, # with damping, # damped by phase",ncnvg)

       Call Out ("First, # districts, then lists, by # clusters in district best cluster sets",ic)
       Do i= 1,5
         Write (8,"(A,I2, A,30I3)", Iostat=ios) "# clusters",i, "  districts",Lsdt(:ic(i),i)
       End do
       
       Call Out (-1,"Mean & sig for top cand ratios to np (# top, # sig, # viable)", Top_ratios)

       Call Out (-1,"Histogram of the convergence data", Cnv_hist)
       Call Out (  "Row 1 = histogram intervals. Final column = max data values (for scale)")
       Write (8,'(A)')  "  Row 2 = # regular clusters.   Row 3 = # convergence calls"
       Write (8,'(A)')  "  Row 4 = # merge / delete ops. Row 5 = # centroiding iterations"
       Write (8,'(A)')  "  Row 6 = # function calls"
       
       Call Out ("Average over districts of # cluster set centroiding ops",nGen_clust/nDst, ln=1)
       Call Out ("Count only operations in Clusters_from_means",nGen_clust_CM/nDst)
       Call Out ("Count only operations in Newton_update",nGen_clust_NT/nDst)

       ns= Min(ns,Msc);  Allocate(ky(ns))
       Call List_of_true (obj_dat(4,:ns) >= 0.05, n1,ky)
       Call Out ("Convergences with poor residuals",ky(:n1))
       Write(8,'(A,20I3)') "with # clusters", cnv_dat(2,ky(:n1))
       Write(8,'(A,20F7.3)') "& clustering objectives",obj_cmp(1,ky(:n1))

       Call List_of_true (cnv_dat(5,:ns) >= 30, n2,ky)
       Call Out ("Slow convergences",ky(:n2))

       If (nMis > 0) then
         Call Out ("# districts where the 'best' cluster set /= top set",nMis, ln=1)

         Call Out (-1, "District, best clust set, most prob set, # clust(best), # clust(top)", &
                        MismatchI(:,:nMis))
         Call Out (-1, "Top full obj, clust obj, fitness,  best full obj, clust obj, fitness", &
                        MismatchR(:6,:nMis))

         Call Out (-1, "For best cluster set: clustering obj ratio & full obj ratio", &
                       MismatchR(7:8,:nMis))
           x1= Sum(MismatchR(7,:nMis)) / nMis;  x2= Sum(MismatchR(8,:nMis)) / nMis
         Call Out ("Mean of the clustering ratios",x1, "Mean of the full ratios",x2, ln=1)
       End if

         x1= Sum(cnv_dat(3,:ns)) / (ns+eps);  x2= Sum(cnv_dat(4,:ns)) / (ns+eps)
       Call Out ("Avg # convergence cycles per cluster set",x1, ln=1)
       Call Out ("with avg # merger & deletion calls",x2)

       x1= Sum(cnv_dat(5,:ns)) / (ns+eps);  x2= Sum(cnv_dat(6,:ns)) / (ns+eps)
       Call Out ("Avg # centroiding iterations per cluster set",x1)
       Call Out ("with avg # function calls",x2)

       If (pr_out >= 1.0) then
         Call Out ("Data by district, then cluster set")
         l1= 1

         Do id= 1,nDst
           n0= Elc1(0,id);  i= (Dist_1-1) + id
           Write(8,'(A)') " ";  Write(8,'(A, I4, A,I3)') "For district", i, "   with # elected",n0

           l2= Last_true(cnv_dat(1,:ns) == i)  ! district cluster sets from 'l1' to 'l2'
           If (l2 < 1) then
             k= -1;  Cycle  ! No cluster sets recorded
           End if
           
           Write(8,'(2X,A)') "1: top elected sets,  2: #_clust, #_cycle, #_reduc, #_iter, #_func"
           Write(8,'(2X,A)') "3: min_size, max_corr, ind_size, final_res,  4: obj = coherent vote * penalty"
             
           Write (8,"(10("//c1(n0)//"I3, 2X))") (Elc1(j,id),j=1,n0), ((Elc2(j,i),j=1,n0),i=l1,l2)
           Write(8,'(A)') " "

           Do i= l1,l2
             Write (8,'(5I4)')      cnv_dat(2:,i)
             Write (8,'(3X,4F7.3)') obj_dat(:,i)
             Write (8,'(3X,3F7.3)') obj_cmp(:,i)
           End do
           l1= l2+1
         End do
       End if

       Write(8,'(A)') " ";  
       Write(8,'(2X,A)') "Non-random initial cluster set & its converged set + their optimal elected sets"
       Write(8,'(2X,A)') "Next, their full objective ratios"

       mat_nonEl= 0;  mat_nonOb= 0;  mat_id= 0

       Do id= 1,nDst
         n0= Elc1(0,id)

         Do i= 1,NonR
           q= nonRnd_elc(0,i,1,id);  iq= N_cand + i

           If (q > 0) then  ! Initial set converged
             If (pr_out > 1) then
               Write (8,"(2I3,2X, 2("//c1(n0)//"I3, 2X))")  iq, q, (nonRnd_elc(j,i,1,id), j=1,n0), &
                                                                   (nonRnd_elc(j,i,2,id), j=1,n0)
               Write (8,'(9X,F7.3)') nonRnd_ratio(i,id)
             End if

             If (All(nonRnd_elc(1:n0,i,1,id) == nonRnd_elc(1:n0,i,2,id))) then
               n= mat_nonEl(i,1) + 1;  mat_nonEl(i,1)= n;  mat_id(n,i,1)= id
             End if
             If (All(nonRnd_elc(1:n0,i,1,id) == Elc1(1:n0,id))) then
               n= mat_nonEl(i,2) + 1;  mat_nonEl(i,2)= n;  mat_id(n,i,2)= id
             End if
           Else             ! Did not converge
             If (pr_out > 1) then
               Write (8,"(2I3,2X, "//c1(n0)//"I3)")  iq, 0, (nonRnd_elc(j,i,1,id), j=1,n0)
             End if
           End if

           If (All(nonRnd_elc(1:n0,i,2,id) == Elc1(1:n0,id))) mat_nonEl(i,3)= mat_nonEl(i,3) + 1
           If (nonRnd_ratio(i,id) == 1000) mat_nonOb(i)= mat_nonOb(i) + 1
         End do
       End do

       If (pr_out > 1) then
         Call Out ("For elected set matches for non-random methods")
         Do i= 1,NonR
           Call Out ("For non-random method",i,ln=1)
           n= mat_nonEl(i,1);  Call Out ("Initial to converged districts",mat_id(:n,i,1))
           n= mat_nonEl(i,2);  Call Out ("Initial to optimal districts",mat_id(:n,i,2))
         End do
       End if

       Call Out ("# matches of non-random initial set electeds to converged set electeds",mat_nonEl(:,1))
       Call Out ("# matches of non-random initial electeds to optimal electeds",mat_nonEl(:,2))
       Call Out ("# matches of non-random converged electeds to optimal electeds",mat_nonEl(:,3))

       Call Out ("# matches of non-random converged set to optimal full objective",mat_nonOb)

       If (pr_out >= 1.5) then
         Call Out ("For all clusters: size, dens, boun")
         Allocate(db(ns,4))
        
         j= 0;  k= 0 
         District_loop3 : Do id= 1,nDst
           m= cls_dat(id)%m
           Call Out ("For district #",id,"with total # clusters",m, ln=1)

           i= 0;  nq= cls_dat(id)%q
           Clust_set_loop : Do q= 1,nq
             ncl= cls_dat(id)%ls(q)
             Do cl= 1,ncl
               i= i+1;  Write (8,"(2I3,2X, 3F7.3)", Iostat=ios) &
                        q, ncl, cls_dat(id)%M1(:3,i)
             End do
             If (ncl > 1) then
               j= j + 1;  Write (8,"(5X, F7.3,1X,2F7.3, 4X,2F7.3,2X,2F7.3)", Iostat=ios) &
                          obj_cmp(:,j), obj_dat(:4,j)
               k= k + 1;  db(k,:)= cls_dat(id)%M1(2:5,i)  
             End if
           End do Clust_set_loop
         End do District_loop3

         Call Sort (.true.,db(:,1),ky);  db(:,3)= db(ky,3)
         Call Sort (.true.,db(:,2),ky);  db(:,4)= db(ky,4)

         Call Out ("For cluster sets with >= 2 cluster")
         Call Out ("Ordered cluster densities",db(:,1))
         Call Out ("Their penalties",db(:,3))
         Call Out ("Ordered cluster boundaries",db(:,2))
         Call Out ("Their penalties",db(:,4))
       End if


       Call Out (" ")
       Call Out ("Non-clustering statistics")

       Call Out (-1,"Frac with >= top 'n' non-clust electeds match top 'n' clust electeds", &
                     Frac_top(:,1,:), fm=fm3)
       Call Out (-1,"Frac with >= top 'n' non-clust electeds match in order", &
                     Frac_top(:,2,:), fm=fm3)

!      # districts where each non-clustering elected set differs from 
!      the top clustering elected set, also where the differences 
!      of the full objective are small or large

       non_ls= -1;  non_rd= -1;  non_ky= -1

       Do i= 1,nMt
         Msk0= Non_ratio(1:,i) > 0;                  k0= Count(Msk0)  ! # districts with valid ratio
         Msk1= Msk0 .and. Non_ratio(1:,i) <= 0.999;  k1= Count(Msk1)  ! # districts with ratio < 1)

         If (k0 > 0) Non_ratio(-1,i)= Sum(Non_ratio(1:,i)) / k0        ! Mean ratio over valid districts
         If (k1 > 0) Non_ratio(0,i) = Sum(Non_ratio(1:,i), Msk1) / k1  ! Mean ratio with ratio < 1

         Call List_of_true (Msk1, k1,non_ls(:,i));  non(1,i)= k1  ! List of districts with ratio < 1
           non_rd(:k1,i)= Non_ratio(non_ls(:k1,i),i)              ! List of their ratios
         Call Sort (.true.,non_rd(:k1,i), non_ky(:k1,i))          ! These ratios in increasing order

         n2= Last_true(non_rd(:k1,i) < 0.90);  non(2,i)= n2       ! # districts with ratio < 0.90 
         n3= Last_true(non_rd(:k1,i) < 0.50);  non(3,i)= n3       ! # districts with ratio < 0.50 
         n4= Last_true(non_rd(:k1,i) < 0.15);  non(4,i)= n4       ! # districts with ratio < 0.15 

         non_ky(:k1,i)= non_ls(non_ky(:k1,i),i)                   ! List of all districts with ratio < 1,
                                                                  ! sorted by increasing ratio
         If (k0 > 0) then
           Exact_match(i)= (k0 - k1) / Real(k0)                     ! Fraction with exact match
           Good_match(i) = (k0 - n2) / Real(k0)                     ! Fraction with ratio >= 0.90
         End if
       End do

       Call Out ("Non-clustering labels");  Write (8, "(7A)") Method_lab(1:)

       Call Out (1,"# districts where non-clust differs from clust: exact, 0.90 obj, 0.50, 0.15",non)
       m= Maxval(non(1,:))
       Call Out (-1,"Based on non-clustering to clustering full obj ratios",non_rd(:m,:))
       Call Out (-1,"For corresponding districts",non_ky(:m,:))

       Call Out ("Frac of each non-clustering that is an exact match to clustering",Exact_match)
       Call Out ("Frac of each non-clustering that is a very good match to clustering",Good_match)

       Call Out ("Mean non-clustering to clustering ratio of full objectives, all districts",Non_ratio(-1,:))
       Call Out ("Mean non-clustering to clustering ratio of full objectives, when ratio < 1",Non_ratio(0,:))


       Call Out (" ")
       Call Out ("STV / DTV statistics")

       DSTV_dif(1)= Count(DtoS_ratio(:,1) /= 1.0)
       DSTV_dif(2)= Count(DtoS_ratio(:,2) /= 1.0)

       If (DSTV_dif(1) > 0) mean_SDTV(1)= Sum(DtoS_ratio(:,1), DtoS_ratio(:,1) /= 1.0) / DSTV_dif(1)
       If (DSTV_dif(2) > 0) mean_SDTV(2)= Sum(DtoS_ratio(:,2), DtoS_ratio(:,2) /= 1.0) / DSTV_dif(2)


       Call Out (-1, "# STV clusters for init, extended, converged",ncl_STV(:,:nDst))
       Call Out (-1, "# DTV clusters for init, extended, converged",ncl_DTV(:,:nDst))

       Call Out (-1, "STV full objectives extended, converged, best",obj_STV(:,:nDst))
       Call Out (-1, "DTV full objectives extended, converged, best",obj_DTV(:,:nDst))

!      STV / DTV histograms

       Do id= 1,nDst
         x1= Non_ratio(id,1);  i= Last_true(x1 > STV_hist(:,0)) - 1

         If (i >= 0 .and. i <= ni) then
           STV_hist(i,1)= STV_hist(i,1) + 1   ! Count the STV ratios
         End if

         x2= Non_ratio(id,2);  j= Last_true(x2 > STV_hist(:,0)) - 1

         If (j >= 0 .and. j <= ni) then
           STV_hist(j,2)= STV_hist(j,2) + 1   ! Count the DTV ratios
         End if
       End do

       STV_hist(:,1)= STV_hist(:,1) / nDst
       STV_hist(:,2)= STV_hist(:,2) / nDst
       Call Out (-1,"Histogram of the STV & DTV best full objective to clustering ratios", STV_hist)

!      STV & DTV full objectives for their converged cluster sets (#19, #20)

       non_lsa= -1;  non_rda= -1;  non_kya= -1

       Do i= 1,2
         Msk0= SD_ratio(1:,i) > 0;  k0= Count(Msk0)
         Msk1= Msk0 .and. SD_ratio(1:,i) <= 0.999;  k1= Count(Msk1)

         SD_ratio(-1,i)= Sum(SD_ratio(1:,i)) / k0
         SD_ratio(0,i) = Sum(SD_ratio(1:,i), Msk1) / k1

         Call List_of_true (Msk1, k1,non_lsa(:,i));  nona(1,i)= k1
           non_rda(:k1,i)= SD_ratio(non_lsa(:k1,i),i)
         Call Sort (.true.,non_rda(:k1,i), non_kya(:k1,i)) ! Sorted list of ratios > 0 .and. < 1

         n2= Last_true(non_rda(:k1,i) < 0.90);  nona(2,i)= n2
         n3= Last_true(non_rda(:k1,i) < 0.50);  nona(3,i)= n3
         n4= Last_true(non_rda(:k1,i) < 0.15);  nona(4,i)= n4

         non_kya(:k1,i)= non_lsa(non_kya(:k1,i),i)          ! Corresponding districts

         Exact_converged(i)= (k0 - k1) / Real(k0)
         Good_converged(i) = (k0 - n2) / Real(k0)
       End do
        
       Call Out ("Frac of STV & DTV converged objectives = exact match to top clustering",Exact_converged)
       Call Out ("Frac of STV & DTV converged objectives = good match to top clustering",Good_converged)

       Call Out ("Mean converged STV & DTV to clustering ratio of full objectives, over all districts",SD_ratio(-1,:))
       Call Out ("Mean converged STV & DTV to clustering ratio of full objectives, when ratio < 1",SD_ratio(0,:))

       Call Out (1,"STV & DTV converged vs clust: # elc dif, # < .90, .50, .15 obj",nona)
         m= Maxval(nona(1,:))
       Call Out (-1,"STV & DTV converged to full obj clustering ratios",non_rda(:m,:))
       Call Out (-1,"For corr. districts",non_kya(:m,:))

       Do id= 1,nDst
         x1= SD_ratio(id,1);  i= Last_true(x1 > STV_hist(:,0)) - 1

         If (i >= 0 .and. i <= ni) then
           STV_hist(i,1)= STV_hist(i,1) + 1   ! Count the STV ratios
         End if

         x2= SD_ratio(id,2);  j= Last_true(x2 > STV_hist(:,0)) - 1

         If (j >= 0 .and. j <= ni) then
           STV_hist(j,2)= STV_hist(j,2) + 1   ! Count the DTV ratios
         End if
       End do

       STV_hist(:,1)= STV_hist(:,1) / nDst
       STV_hist(:,2)= STV_hist(:,2) / nDst

       Call Out (-1,"Histogram of the STV & DTV converged full objective to clustering ratios", STV_hist)

!      STV to DTV comparisons

       Call Out ("Case 0: # districts where STV & DTV elected sets differ", &
                 SD_cnt(1), ln=1)
       Call Out ("Case 1: # districts where STV & DTV converged cluster sets differ", &
                 SD_cnt(2))

       Call Out ("Mean ratio of DTV to STV best full objectives ", mean_SDTV(1), ln=1)
       Call Out ("vs mean for the converged ",mean_SDTV(2))
       Call Out ("# differ: DTV to STV full objectives for best cluster sets", DSTV_dif(1))
       Call Out ("vs # for the converged sets",DSTV_dif(2))
       Call Out (-1,"Ratios of DTV to STV full objective: best sets then converged", &
                 DtoS_ratio)

       Call Out (1,"# districts where successive STV electeds differ",STV_dif)
       Call Out ("Initial 1:3 = orig, DHondt, optimal. Converged 4:5 = DHondt, Opt")
     End if
   End Subroutine General_stats
    

   Subroutine General_data (Status,nDst, Non_clust,Elect_dat)

     Integer, Intent(in) :: Status  ! Initialize the statistics if = 0, update if = 1, output if = 2 
     Integer, Intent(in) :: nDst    ! Total # districts

     Type(Multi_listR), Optional, Intent(in) :: Non_clust    ! Data for non-clustering methods
     Type(Multi_listR), Optional, Intent(in) :: Elect_dat(:) ! (me) Data for possible sets of electeds

!  Local:
     Integer, Parameter :: Mp= 5
     Character(12) :: dat_file = "Top_dat1.txt"

     Integer, Allocatable, Save :: nt(:,:)     ! (3,nDst)
     Integer, Allocatable, Save :: Non(:,:)    ! (nMt,nDst)
     Integer, Allocatable, Save :: Elc(:,:,:)  ! (Mp,0:3,nDst)
     Integer, Allocatable, Save :: Obj(:,:)    ! (7,nDst)

     Integer, Save :: id
     Integer :: i, j, k, n, n1, n2, np, ncl, nonC, ios

     If (Status < 1) then  ! Initialize
       If (Allocated(nt)) DeAllocate(nt, Non, Elc, Obj)
       Allocate (nt(3,nDst), Non(nMt,nDst), Elc(Mp,0:3,nDst), Obj(7,nDst))
       nt= 0;  Non= 0;  Elc= 0;  Obj= 0;  id= 0

     Else if (Status == 1) then ! Update
       id= id + 1;  np= Size(Elect_dat(1)%ls);  n1= Min(np,Mp);  ncl= Elect_dat(1)%n 
       nonC= Size(Non_clust%L0,2);  n2= Min(nonC,3)

       nt(1,id)= np;  nt(2,id)= ncl;  nt(3,id)= n2
       Non(:,id)= Non_clust%vl

       Elc(:n1,0,id)= Elect_dat(1)%ls(:n1)
       Elc(:n1,1:n2,id)= Non_clust%L0(:n1,:n2)
       Obj(:,id)= Nint(100*Elect_dat(1)%M1(:,1))

     Else                       ! Output
       Open(7, File=dat_file, IOstat=ios, Status='Replace', Action='Write')

       Write(7,'(A)') "id,np,ncl, non_map(:). top elected sets (clustering,non-clustering) + objective data"
       Do k= 1,nDst
         Write(7,'(3I3,2X,7I2, 2X,5I3, 2X,3(2X,5I3))')  k, nt(:2,k), Non(:,k), (Elc(:,j,k),j=0,nt(3,k))
         Write(7,'(10X,I5,(2X,2I4), 2X,2(2X,2I4))') Obj(:,k)
       End do
       Close(7)
     End if

   End Subroutine General_data

End Module Clusters6