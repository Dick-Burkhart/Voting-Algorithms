
!    This module contains parameters for Clustering_sng    

 Module Clusters0

   Use Types
   Implicit None
   
!  Values that differ only by district
   
   Integer, Parameter :: mx_Dist= 40 ! Max # voting districts

   Integer,      Save :: Rating      ! = 1 for ratings 
                                     ! = 0 for weak ranking (may be equal)
                                     ! = -1 for strong ranking
   Integer            :: mr_spec= 8  ! Limit # ranked or rated candidates to be used
   
   Real(Dblp), Parameter :: Dot_fac_rank= 0.0d0  ! Scale factor for the modified dot product: 
                                                 !   Dot_product_M1, for ranking data.
   Real(Dblp), Parameter :: Dot_fac_rate= 0.50d0 ! Scale factor for the modified dot product: 
                                                 !   Dot_product_M1, for rating data.
   Real(Dblp), Save      :: Dot_fac  ! Ordinary dot product if < 0 or >=1
                                     ! else scales component products where
                                     ! both factors are negative - used for ranking data.

   Real, Parameter :: Noise_max= 4.00  ! Max value of Noise_cor
   Real, Parameter :: Noise_min= 2.50  ! Min value of Noise_cor
   Real, Save      :: Noise_cor        ! Noise significance value for correlations

   Real, Parameter :: Noise_fac= 1.30  ! Use to compute 'x2' for the statistical 'Noise_cor' in Clusters1
   Real, Parameter :: Por_fac= 0.50    ! Noise_por= Por_fac * Noise_cor
   Real, Save      :: Noise_por        ! Noise level for portions


!  Parameters for Read_ballots:  ***
   
   Integer,  Parameter :: nMt= 7   ! # non-clustering methods computed
   
   Character(11), Save :: Method_lab(0:nMt)= &
                          (/"Clustering ", "STV        ", "DTV        ",  &
                            "IRV        ", "Top Vote   ", "Borda      ",  &
                            "Condorcet  ", "Approval   "/)
   Logical, Save    :: Method(nMt)  !  Computed non-clustering methods
   Logical, Pointer :: Req_cand(:)  ! (nc) Candidates required to be electable, based on non-clustering methods, etc
   Logical, Pointer :: Pos_cand(:)  ! (nc) The positive candidates for rating data (electeds should be a subset)
    
   Real(Dblp), Parameter :: Max_ptD= 10.0d0  ! Max and min point values for ballot ratings.
   Real, Parameter :: Max_pt= Max_ptD
   
   Real, Parameter :: Bpm1= 0.9999, Bpm2= 0.50, Bpm3= 0.35, Bpm4= 4.00 
   Real, Parameter :: Parm_bal1(4)= (/ Bpm1, Bpm2, Bpm3, Bpm4*Bpm3 /)  ! Full Borda + alternative increment

   Real,      Save :: Parm_bal(4)= Parm_bal1
                      ! Parameters for converting ranking levels to points
                      ! i = 1: Scale factor for the points of unranked candidates, 
                      !   0 <= p <= 1, with 1 = "full Borda" = all ballots get
                      !     the same # total points, with more going to the 
                      !     unranked when fewer are ranked.
                      ! i = 2: Parameter determining the increase in point values.
                      !   0 <= p <= 1: Point increase in the increment for an 
                      !     arithmetic increase. The initial increment is 1, 
                      !     corresponding to the standard Borda Count = 
                      !     1, 2, 3,...,mr for p = 0.  More generally, 
                      !     the point values are 1,2+p,3+3p,4+6p,..,
                      !     mr + ((mr-1)mr/2)p, with = 0.5 typical.
                      !   p > 1 signifies a geometric increase in point values, 
                      !     so that p = 2 means doubling: 1,2,4,8,...2**(mr-1).
                      ! i = 3: The sigma for ranking level 1 as a random 
                      !     variable, assumed to be the most certain, hence the
                      !     smallest sigma, over all the ranking levels. p = 1/3
                      !     is typical. The sigma of the corresponding point 
                      !     values is the level sigma scaled by the point 
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
   
!  Parameters for Consolidate_ballots  *****
   
   Integer, Parameter :: Max_sl = 150  ! Max # slate ballots for Rating < 1
   Integer, Parameter :: Max_sl1= 100  ! Partial # slate ballots for Rating = 1
   Real,    Save      :: Min_slate     ! Minimum weight for slate ballots
   Integer, Parameter :: Msl= 4        ! # top candidates recorded in 'Init'
   
!  Parameters used for Initialize_clusters in Form_clusters  *****
   
   Real(Dblp), Parameter :: Min_wt = 0.05d0  ! Minimum size for provisional regular clusters
   Real      , Parameter :: Min_wti= 0.01d0  ! Minimum size for the independents
   Real(Dblp), Parameter :: Min_mx = 0.50d0  ! Minimum of the maximum rating for independents

   Integer, Parameter :: N_cand1= 10     ! Reduced # initial cluster sets starting from random subsets
   Integer, Parameter :: N_init1= 12     ! Reduced total # initial cluster sets
   Integer, Parameter :: N_cand= 15      ! Standard # initial cluster sets starting from random subsets
   Integer, Parameter :: N_init= 20      ! Standard total # initial cluster sets

   Integer, Parameter :: Clust_opt= 1    ! Option for how to compute the degree and triangle graph clustering values
   Real,    Parameter :: Corr_cut= 0.05  ! Lower correlation cutoff for the slate graphs

!   Real,    Parameter :: Core_parm(4)= (/ 0.15, 0.40, 0.60, 0.50 /) ! Parameters for 'Core_clustering'
   Real,    Parameter :: Core_parm(4)= (/ 0.10, 0.45, 0.65, 0.45 /) ! Parameters for 'Core_clustering'
                                         !  1 = correlation cutoff for the correlation graph
                                         !  2 = correlation cutoff for cluster separation
                                         !  3 = correlation cutoff for a core set                                              
                                         !  4 = total membership cutoff for small clusters
   Integer, Parameter :: Simple_opt= 1   ! Simple options for electing the strongest candidates from 
                                         !   cluster weights and strength vectors with no proportionality 
                                         !   constraints. Used by 'Simple_elect'.
                                         ! 1 : For each candidate 'c', compute 'value(c)' =  sum of the 
                                         !      reduced ratings over all clusters.
                                         ! 2 : For each candidate 'c', compute 'value(c)' =  sum of the 
                                         !     squares of the reduced ratings over all clusters.
                                         ! 3 : For each candidate 'c', compute 'value(c)' =  the maximum 
                                         !      reduced ratings over all clusters.
                                         ! In all cases, sort 'value' into decreasing order and 
                                         ! elect the top 'np' candidates.
   
   Real, Parameter :: STVmd_parm(3)= (/0.50, 0.25, 2.0/)
                      ! Merge & delete parameters for STV clustering
                      !   (1) = cluster overlap merge criterion
                      !   (2) = cluster size deletion criterion
                      !   (3) = minimum number of clusters upon output
   
   Real, Parameter :: Cnv_mrg(3)= (/0.55, 0.25, 5.0/)
                      ! Merge and delete parameters for initial clustering
                      !   (1) = cluster overlap merge criterion
                      !   (2) = cluster deletion criterion
                      !   (3) = minimum number of clusters upon output; 
                      !         convergence clustering beyond this
   
!  Parameters used for Clustering in Form_clusters  *****

   Logical, Parameter :: Test_range= .true.  ! Use Fit_range instead of Fit_all

   Integer, Parameter :: mxcl= 5   ! Max # regular clusters in a converged cluster set
   Integer, Parameter :: mxcl0= 8  ! Maximum of regular clusters in an initial cluster set

   Real,    Parameter :: Parm_cor(2)= (/0.035d0, 0.70d0/) 
                         ! Lower and upper Correlation bounds for cluster
                         ! memberships
   Real,    Parameter :: Parm_top(2)= (/0.67d0,0.25d0/) 
                         ! Top & secondary candidate ratio in 'Best_cand'
   
   Integer, Parameter :: Mxitr(3)=  (/ 6, 10, 15 /) 
                               ! Iteration limits for Restart levels to measure
                               ! progress at each iteration, also convergence, 
                               ! using Parm_cnv(:,1)
   Real,   Parameter :: Parm_cnv(3,5)= &
                          (/ (/ 0.10d0, 0.05d0, 0.01d0 /), &
                             (/ 0.80d0, 0.65d0, 0.50d0 /), &
                             (/ 0.20d0, 0.35d0, 0.50d0 /), &
                             (/ 0.15d0, 0.30d0, 0.50d0 /), &
                             (/ 0.10d0, 0.10d0, 0.75d0 /) /)
                               ! For convergence Restart phase = 1, 2, 3:
                               ! (RS,1) = convergence criteria
                               ! (RS,2) = min cluster correlation criteria for merging
                               !          of correlated clusters
                               ! (RS,3) = max cluster size criterion for deletion of 
                               !          small isolated clusters
                               ! (RS,4) = max cluster correlation criterion for deletion 
                               !          of small isolated clusters
                               ! (1,5) = maximum norm permitted for an update differential
                               ! (2,5) = match tolerance for measures of cluster similarity
                               ! (3,5) = min separation criterion for the separation 
                               !         graph for correlation clique of clusters 
   
   Real, Parameter :: Parm_obj(12)= (/ 0.15d0,0.45d0,  0.10d0,0.35d0,  &
                                       0.25d0,0.75d0,  0.25d0,0.75d0,  &
                                       0.40d0,0.80d0,  0.025d0,0.250d0 /) 
                               ! Parameters for the objective penalty factors
                               !  1:2 = cluster density  soft cutoff limits (cos_rise)
                               !  3:4 = cluster boundary soft cutoff limits (cos_fall)
                               !  5:6 = min cluster size soft cutoff limits  (cos_rise)
                               !  7:8 = max cluster correlation soft cutoff limits (cos_fall) 
                               !  9:10 = independents size (normalized) cutoff limits (cos_fall)
                               !  11:12= convergence residual limits (cos_fall)
   
!  Parameters for Evaluate_candidates  *****
   
   Integer, Parameter :: mfe= 50  ! Max # sets of electeds per cluster set
   Integer, Parameter :: mft= 100 ! Max # sets of electeds per cluster set
   Integer, Parameter :: mfi= 1   ! Max # sets of electeds per initial 
                                  !   or special cluster set
   Integer, Parameter :: mxe= 150 ! Max # sets of electeds per district

   Real, Parameter :: Parm_rep(6)= (/0.01, 0.25, 0.50, 0.75, 0.90, 0.95/)
   
   Real, Parameter :: Parm_el(2)= (/0.20, 0.70/)  !  (/0.1875, 0.625/) for stronger proportionality penalties
                      ! Parameters for penalties for deviation between true 
                      ! cluster elected size (via 'portions') in position units.
                      ! No penalty for deviations <= Parm_el(1), max penalty for >= Parm_el(2)

   Real, Parameter :: obj_eps= 0.05   ! Lower limit for cluster or elected set 
                                      ! individual penalties
   Real, Parameter :: obj_eps2= 0.01  ! Lower limit for the product of these
                                      ! cluster or elected set penalties

   Real, Parameter :: obj_eps3= 0.005 ! Lower limit for a cluster set objective

!  Parameters for Testing and Output  *****
   
   Real, Parameter :: fac_full= 10000.0  ! Normalization level for the top full
                                         ! objective (cluster set + elected set)
   Integer, Parameter :: nct= 2
   Integer,      Save :: cls_n(nct), elc_n(nct)
   Real,    Parameter :: cls_cut(nct)= (/ 0.85, 0.50 /)  
                                         ! Clustering objective cutoff fractions for
                                         ! the top 'nct' cluster set categories

   Character(13) :: Indat=  "Elections123/"
   Character(8)  :: Outdat= "Results/"

   Integer, Parameter :: Mvt= 3 ! Max # candidates to be regarded as 
                                ! a voter's top choices
   Real,    Parameter :: Votr_wt(Mvt)= (/ 4./7., 2./7., 1./7. /) 
                         ! Weights for a voter's top choice, 2nd choice,...,
                         ! as used to compute overall voter satisfaction
       
   Real,    Parameter :: Opt_wt(3)= (/ 0.70, 0.20, 0.10 /) 
                         ! Weights for combining 3 measures of the success of 
                         ! proportionalrepresentation by clustering:
                         ! (1) full objective, 
                         ! (2) voter satisfaction, and 
                         ! (3) STV ratio

!  Random seed parameters  ***

   Integer, Parameter  :: Init_clust_seed(1)= 875615 
   Integer, Parameter  :: Rnd_trial_seed(1) = 395261

!  Cluster set centroiding/membership counts

   Integer, Save :: nGen_clust_CM  ! # calls to 'Generate_clusters' from 
                                   !   Clusters_from_means
   Integer, Save :: nGen_clust_NT  ! # calls to 'Generate_clusters' from 
                                   !   Newton_update
   
   Integer, Save :: ngc(3), ngcR(3)  ! Used for Clust_gen counts:  
                    ! (1)= nGen_clust [all calls]
                    ! (2)= nGen_clust_CM [all calls in Clusters_from_means]
                    ! (3)= nGen_clust_NT [all calls in Newton_update]

!  Parameters for sensitivity studies and other run statistics

!  For use in Rnd_stats: To compare the clustering from randomized slate 
!  clusters to the standard clustering, when the slate ballots are 
!  algorithmically determined from the original ballots. For each 
!  district there is a standard run, followed by the randomized runs, 
!  accumulating statistics.

   Integer, Parameter :: nt7= 3  ! Max # top elected sets from each 
            ! randomized run to be matched to the top initial elected 
            ! sets, both sets required to have distinct best cluster sets 
   Integer, Parameter :: ni7= 3  ! Mas # top initial elected sets, 
            ! with distnct best cluster sets
   Integer, Parameter :: ki7= 3  ! Max # top elected sets from each 
            ! of the associated cluster sets 
   Integer, Parameter :: nx7= Max(ni7,nt7)
   
!  For each sensitivity parameter "sens_parm": a standard value 
!  = Std_sens(sens_parm) plus deviations. 
   
   Integer, Parameter :: Nsens= 4  ! # sensitivity parameters
   Integer, Parameter :: Mndev= 2  ! Min # deviations
   Integer, Parameter :: Mxdev= 6  ! Max # deviations

!  Sensitivity parameters to test:

   Character(40) :: Sens(Nsens)= & 
                 (/ "Unranked candidate point scale factor",    &
                    "Ranked candidate point ascent increment",  &
                    "Top candidate ballot point sigma",         &
                    "Multiplier for the unranked ballot sigma" /)

!  Standard values:
   Real :: Std_sens(Nsens)= (/ Bpm1, Bpm2, Bpm3, Bpm4/)
   

   Real, Save :: SenV2(2,Nsens)= (/ (/-0.20, -0.10/), (/-0.10, +0.10/), & 
                                    (/-0.05, +0.05/), (/-0.25, +0.25/)/)

   Real, Save :: SenV3(3,Nsens)= (/ (/-0.30, -0.15, -0.05/), (/-0.10, +0.10, +0.20/),  & 
                                    (/-0.10, -0.05, +0.05/), (/-0.25, +0.20, +0.40/) /)

   Real, Save :: SenV4(4,Nsens)= (/ (/-0.40,-0.20, -0.10, -0.05/),   &
                                    (/-0.20, -0.10, +0.10, +0.20/),  &
                                    (/-0.12, -0.08, -0.04, +0.05/),  &
                                    (/-0.40, -0.20, +0.20, +0.40/) /)

   Real, Save :: SenV5(5,Nsens)= (/ (/-0.40, -0.20, -0.10, -0.05, -0.025/), &
                                    (/-0.20, -0.10, +0.10, +0.20, +0.40/),  &
                                    (/-0.12, -0.08, -0.04, +0.04, +0.08/),  &
                                    (/-0.40, -0.20, +0.20, +0.40, +0.60/) /)

   Real, Save :: SenV6(6,Nsens)= (/ (/-0.50, -0.25, -0.10, -0.05, -0.025, -0.0125/), &
                                    (/-0.30, -0.10, -0.05, +0.05, +0.10, +0.30/),    &
                                    (/-0.10, -0.04, -0.01, +0.01, +0.04, +0.10/),    &
                                    (/-0.40, -0.10, -0.05, +0.05, +0.10, +0.40/) /)

   Real,    Save :: SenP(0:Mxdev)  ! Sensitivity parameter deviation values for the
                                   ! current parameter 'sens_parm'
                                   ! with SenP(0) = the standard value Std_sens(sens_parm)
   Integer, Save :: ndev           ! # deviation deltas, 1 to Mxdev, for sens_parm

   Logical, Save :: Standard       ! Standard vs Perturbed run

 End Module Clusters0

                     
