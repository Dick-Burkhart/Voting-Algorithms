!            This module contains "Form_clusters" and its subroutines. 
    
!  "Form_clusters" is called from Clustering_PR to compute the best 
!   sets of clusters from the slate ballots 'Memb'. First
!  'Initialize_clusters'is called to compute a variety of initial 
!  cluster sets. Then 'Clustering' is called to converge these to
!  to final cluster sets 'clust_set', merging strongly overlapping
!  clusters and deleting very small, isolated clusters. Duplicate
!  cluster sets are deleted.
    
Module Clusters3

   Use Clusters0
   Use Clusters4
   Use Clusters_support

   Use Newton_operators   
   Use Cholesk
   Use Graph_algorithms
   Use Factorials
   Use Util
   Use Output
   Use Types
   Use Precisn
   Use IEEE_Arithmetic
   Implicit None

 Contains
 
   Subroutine Form_clusters (idist, nc,np, Memb,Mean_rnd, Non_clust,Init, ncs,Clust_set)
   
!    Identify a number of initial cluster sets, based both on random selections of slate ballots 
!    and on heuristic or non-optimal algortihms.
   
!    Then converge each initial cluster set to a final cluster set, if possible, from 1 to 5 clusters per set.
   
!    Retain a new converged cluster set if it differs significantly from the prior converged set but only the best 3 
!    cluster sets for each possible # of clusters per set.

   
     Integer,           Intent(in) :: idist     ! Voting district index
     Integer,           Intent(in) :: nc        ! # candidates
     Integer,           Intent(in) :: np        ! # candidates to be elected
     Type(Multi_listD), Intent(in) :: Memb(0:)  ! (0:ns) Final slate cluster data
                                                ! 0 Case: 
                                                !   %n = nc = # candidates 
                                                !   %l = np = # candidates to be elected = total weight
                                                !             in position units
                                                !   %qx(2)  = correlation parameters
                                                !   %vl(mlv). (lv) = # slate clusters with
                                                !      # ranked or rated candidates <= lv

                                                ! 1:nsl Case:
                                                !  %l = # top slate candidates
                                                !  %m = # slate candidates (for Rating = 1)
                                                !  %p = # positively rated candidates (for %tx)
                                                !  %n = # negatively rated candidates (for %tx)
                                                !  %k = # slate ballots positively correlated with and 
                                                !       subsequent to this slate ballot
                                                !  %vl = listing of these by decreasing correlation
                                                !  %qx = the corresponding correlations

                                                !  %fsx= slate membership
                                                !  %fux= slate width
                                                !  %lt(nc) = candidates in preferential order (from %sx)
                                                !  %px(nc) = variance vector for %sx
                                                !  %rx(nc) = mean rating vector with variance normalized weight
                                                !          = %sx * %ux
                                                !  %sx(nc) = mean rating vector
                                                !  %tx(0:nc)= centered (subtract Noise_cor if Rating < 1) 
                                                !             and normalized mean vector
                                                !  %ux(nc) = slate ballot weight as variance normalized 
                                                !          = %fsx / %px

     Real,    Intent(in) :: Mean_rnd(0:,:,:)  ! (0:nc,nr,N_cand1) For randomized cluster set initialization.
                                              ! 'nr' randomly selected mean vectors for each of 'N_cand1"
                                              ! initial cluster sets, with (0,:,:) = cluster weights 
     Type(Multi_listR), Intent(in) :: Non_clust  ! Data for 'nMt' selected non-clustering methods 'm'
                                     !             to elect multiple candidates.
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
     
     Type(Multi_listR), Intent(out) :: Init(:)  ! (N_init) Initial cluster set data
                                                ! %k = cluster set converged to, 0 if none
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

                                                ! %M1(7,4)    Objective data for the electeds %Q0(:,0,:)
                                                ! %M2(ncl,ncl) Cluster correlation matrix, 0 on diagonal

                                                ! %T0(0:nc,0:ind,2) (1:,1:,1) = Candidate mean ratings  
                                                !      for each cluster, with cluster size at (0,1:,1)
                                                !      and cluster average mean vector at (1:,0,1)
                                                !      (1:,1:,2) = Cluster portions with averaged 
                                                !                  noise zeroed ratings at (1:,0,2)

     Integer,          Intent(out) :: ncs          ! Final # cluster sets
     Type(Multi_listR),Intent(out) :: Clust_set(:) ! (N_init) Final cluster sets and associated convergence data
                                                   !        where the clusters in a set are ordered by decreasing weight
                                                   ! %k   = best initial cluster set that converged to this cluster set
                                                   ! %l   = 'ind' = # clusters, including independents
                                                   ! %m   = # initial cluster sets that converged to this cluster set
                                                   ! %ls(m)= list of those initial sets 
                                                   ! %fux = clustering objective ratio = (q)%rx(1) / (1)%rx(1) 

                                                   ! %lt(ncl) Mapping of initial to converged clusters 

                                                   ! %px(6):  Objective data: min size, max correlation, 
                                                   !            independents size (frac), final residual,
                                                   !            Jnorm, Pnorm
                                                   ! %qx(4):  Objective factors: min size penalty, max correlation penalty, 
                                                   !            independents size penalty, residual penalty
                                                   ! %rx(3):   Clustering objective value(1)= coherent vote(2) * penalty(3)
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
                                                   !               
                                                   !  %T0(nc,0:ind,3) Cluster vector data
                                                   !      (:,1:,1) = Cluster mean vectors
                                                   !      (:,1:,2) = Cluster mean vectors with noise zeroing 
                                                   !      (:,1:,3) = Cluster portion vectors
                                                   !      (:,0,1)  = Zrate0 * declining regular portion factor
                                                   !      (:,0,2)  = Zrate0 = cluster averaged noise zeroed mean vectors
! Local:
     Real,    Allocatable :: tmp(:,:)
     Integer, Allocatable :: nl(:), key(:), map(:,:), Prob(:,:)
     Integer :: ind, ncl, ncyc
     Integer :: i, m, q, cl, nr, ns
     
     Call Out ("Enter 'Form_clusters' for district #",idist, ln=1)

     ns= Ubound(Memb,1);  nr= Ubound(Mean_rnd,2)

     Call Initialize_clusters (idist,nc,np,ns,nr, Mean_rnd, Memb, Non_clust%Q0(:,:,:2), &
                               Non_clust%T0(:,1:,:), Init)

     Call Clustering (idist,N_init,nc, Memb, Init,ncs,Clust_set)

     Call Out ("Clustering output for district #",idist, ln=1)
     Call Out ("Total # cluster sets",ncs, "# top cluster sets",cls_n(1), ln=1)
     Call Out ("# viable cluster sets",cls_n(2))

     Allocate (nl(ncs), key(ncs), tmp(3,ncs), map(N_init,2))

     Do q= 1,ncs
       nl(q)= Clust_set(q)%l - 1
       Forall(i=1:3) tmp(i,q)= Clust_set(q)%rx(i)
     End do

     Allocate (Prob(4,ncs))
     Prob(1,:)= Clust_set%m
     Call Sort (.false., Prob(1,:), Prob(2,:))
     Call Inverse_map (Prob(2,:), Prob(3,:))
     Prob(4,:)= nl(Prob(2,:))

     map= 0
     Do i= 1,N_init
       map(i,1)= Init(i)%n  
       q= Init(i)%k;  If (q > 0) map(i,2)= nl(q)
     End do

     If (pr_out > 1) then
       Call Out (1,"# convergences, to clust set #, inverse, # clusters",Prob)

       Call Out ("Initial to final cluster set mapping",Init%k)
       Call Out (-1,"Initial vs final # clusters",map)

       Call Out ("Fractional clustering objectives",Clust_set(:ncs)%fux)
       Call Out (-1,"Cluster set objectives: clust obj= coherent vote * clust penalty",tmp)
     
       Do q= 1,ncs
         ncl= nl(q)
         Call Out ("For cluster set",q, "with # regular clusters ncl",ncl, ln=1)
         Call Out ("Cluster sizes",Clust_set(q)%sx)
         If (ncl > 1) Call Out (1,"Cluster correlation matrix", Clust_set(q)%M1)
         Call Out ("Obj data (min size, max corr, ind size, resid, J & P norms)",Clust_set(q)%px)
         Call Out ("Objective factors",Clust_set(q)%qx)
         Call Out ("Objective = coherent vote * penalties",Clust_set(q)%rx)

         If (Rating < 1) Call Out ("Ranking noise level",Noise_cor)
         Call Out ("Cluster size averaged zeroed ratings ('Zrate')", Clust_set(q)%T0(:,0,2))
         Call Out ("Zrate * (1- por(:,ind)) for ordering", Clust_set(q)%T0(:,0,1))
         Call Out ("Corresponding ordering", Clust_set(q)%L1(1:,0))
       
         If (pr_out > 1.5) then
           Call Out (-1,"Cluster size, dens pen, bound pen, dens, boun, width",Clust_set(q)%M3(:6,:ncl))
           Call Out (-1,"Frac slate cluster fuzzy reduced memb, avg reduction factor",Clust_set(q)%M3(7:,:ncl))
           Call Out (-1,"Original cluster mean vectors", Clust_set(q)%T0(:,1:,1))
           Call Out (-1,"The corresponding portions", Clust_set(q)%T0(:,1:,3))
           Call Out (-1,"Corresponding candidate ordering data by cluster", Clust_set(q)%L1(1:,1:))
           Call Out (-1,"Top #, significant #, viable # candidates by cluster", Clust_set(q)%L1(:0,1:))

           Call Out ("Cluster set convergence data:")
           Call Out ("Initial & final # clusters, # convergence cycles", Clust_set(q)%vl(0:2))
           Call Out ("# merges, deletes, # centroidings, # function calls", &
                    Clust_set(q)%vl(3:6))
         
           Call Out (-1,"Convergence: restart, success, ncl, # itr, # fcn calls", &
                       Clust_set(q)%L0)
           Call Out (-1,"Convergence: final t, residual, net line search change", &
                       Clust_set(q)%M2)
         End if
       End do
     End if

   End Subroutine Form_clusters


   Subroutine Initialize_clusters (idist, nc,np,ns,nr, Mean_rnd, Memb,elc_SD,mean_SD, Init)

!    Compute initial cluster sets for the convergence process. Use a selection of slate cluster sets that are randomly gnerated
!    from slate ballots based on their top 1 or 2 ranked/rated candidates, called candidate cluster sets, alonng with several
!    cluster sets computed by heuristic or non-optimal algortihms.
  
!    These include a "core set" algorithm, a "spectral" algorithm, a "clustering degree" algorithm, 
!    a "triangular coefficient" algorithm, and the STV algorithm.

     
!    Parameters from Cluster0 used: N_init
  
     Integer,           Intent(in) :: idist     ! Index the voting district
     Integer,           Intent(in) :: nc        ! # candidates
     Integer,           Intent(in) :: np        ! # candidates to be elected
     Integer,           Intent(in) :: ns        ! # consolidated slate ballots
     Integer,           Intent(in) :: nr        ! max # regular clusters in an initial cluster set
    
     Real,              Intent(in) :: Mean_rnd(0:,:,:)  ! (0:nc,nr,N_cand) Mean rating vectors for random initial clusters,
                                                !   ordered by decreasing cluster size, with sizes summing to 'np'

     Type(Multi_listD), Intent(in) :: Memb(0:)  ! (0:ns) Slate ballot data
     Integer,  Intent(in) :: elc_SD(:,0:,:)     ! (np,0:2,2)  Elected set: increasing (:,0,:), ranked (:,1,:), 
                                                !             clusters (:,2,:). STV (:,:,1) and DTV (:,:,2)         
     Real,     Intent(in) :: mean_SD(0:,:,:)    ! (0:nc,ind,2) = Cluster data from STV (:,:,1) and DTV (:,:,2).
                                                !    (0,1:,:) = cluster weighs, (1:nc,1:,:) = cluster mean vectors
                                                !     with cluster average mean vector at (1:,0,:)
     
!    Output:
     Type(Multi_listR), Intent(out) :: Init(:)  ! (N_init) Initial cluster set data
                                                ! %k = cluster set converged to, 0 if none
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

                                                ! %M1(7,4)    Objective data for the electeds %Q0(:,0,:)
                                                ! %M2(ncl,ncl) Cluster correlation matrix, 0 on diagonal

                                                ! %T0(0:nc,0:ind,2) (1:,1:,1) = Candidate mean ratings  
                                                !      for each cluster, with cluster size at (0,1:,1)
                                                !      and cluster average mean vector at (1:,0,1)
                                                !      (1:,1:,2) = Cluster portions with averaged 
                                                !                  noise zeroed ratings at (1:,0,2)
      
 ! Local:
     Real, Parameter :: cut(4)= (/0.20, 0.40, 0.60, 0.80/)
     Type(Adjacency) :: Slate_G(ns), Slate_G2(ns)
     Type(Set_list)  :: Core(nr+1) ! The core sets of the clusters for 'Core_clustering'
                                   !  %n       = # vertices in the core set
                                   !  %set(n)  = list of those vertices
                                   !  %mbr(nv) = 0..1 memberships of all vertices in the core
                                   !  %svl     = sum of the core vertex weights
                                   !  %smb     = total weighted membership of the cluster
     Real    :: clust_vec(0:nc,nr+1), clust_memb(ns,nr+1), Coef(ns)
     Real    :: mrg2(3), siz(nr+1), pk_val(2*nr)
     Integer :: cls(nr), cnt(nr+1), key(nr), elc(np,2), peak(2*nr)

     Logical :: ReOrd
     Real    :: max_deg, tot_wt, wt_fac
     Integer :: mc, n1, nf, sl, ng0, ind, ncl, nr1, sum_deg
     Integer :: c, j, n, q, cl, iq

     Call Out ("Enter 'Initialize_clusters' from district",idist, ln=1)
     
     ng0= nGen_clust;  mrg2= Cnv_mrg;  mrg2(3)= np + 1

!    Compute the initial 'candidate' cluster sets
       
     Cand_loop : Do q= 1,N_cand

       cnt= 0;  siz=0;  clust_vec(:,:nr)= Mean_rnd(:,:,q);  nr1= nr + 1

       Call Clust_matrix (q,nc,np,ns,nr,Noise_cor,Memb(1:), clust_vec, &
                          cnt,siz, clust_memb)

       Call Merge_matrix (nr,mrg2, clust_memb(:,:nr), &
                          cnt(:nr),siz(:nr), n1,cls(:nr))
       ind= n1 + 1 
       
       If (n1 < nr) then
         clust_memb(:,ind)= clust_memb(:,nr1)
         cnt(ind)= cnt(nr1);  siz(ind)= siz(nr1)
         clust_memb(:,nr1)= 0;  cnt(nr1)= 0;  siz(nr1)= 0
       End if

       Call Clust_mean (nc,np, Memb(1:),clust_memb(:,:ind), &
                        siz(:ind),clust_vec(:,:ind))

       Allocate(Init(q)%T0(0:nc,0:ind,2));  Init(q)%T0= 0

       Init(q)%n= n1;  Init(q)%T0(:,1:,1)= clust_vec(:,:ind)
       Forall(c=1:nc) Init(q)%T0(c,0,1)= Sum(clust_vec(0,:ind) * clust_vec(c,:ind)) / np

       If (pr_out > 1) then
         Call Out ("For random initial clustering set",q, "# regular clusters",n1, ln=1)
         Call Out ("with cluster weights",clust_vec(0,:ind))
         If (pr_out > 1.5) Call Out (-1,"Their mean vectors",clust_vec(1:,:ind))
       End if
     End do Cand_loop

!    Special methods. First the core clustering algorithm
    
     Slate_G%sum_wt= Memb(1:)%fsx;  mc= Size(Core) - 1
     Call Slate_graph (ns,Core_parm(1), Memb, Slate_G, sum_deg)

     Call Core_clustering (Core_parm(2:), mc,ns, Slate_G, n1,Core)

     ind= n1 + 1;  key= 0;  siz= 0;  clust_memb= 0;  clust_vec= 0

     Do cl= 1,ind
       Where (Core(cl)%mbr > 0.001) clust_memb(:,cl)= Core(cl)%mbr * Memb(1:)%fsx
       siz(cl)= Sum(clust_memb(:,cl))
     End do

     Call Clust_mean (nc,np, Memb(1:),clust_memb(:,:ind), &
                      siz(:ind),clust_vec(:,:ind))

     Allocate(Init(q)%T0(0:nc,0:ind,2));  Init(q)%T0= 0

     Init(q)%n= n1;  Init(q)%T0(:,1:,1)= clust_vec(:,:ind)
     Forall(c=1:nc) Init(q)%T0(c,0,1)= Sum(clust_vec(0,:ind) * clust_vec(c,:ind)) / np

     tot_wt= Sum(siz(:ind)) 

     If (pr_out > 1) then
       Call Out ("For the core clustering set",q, "# regular clusters",n1, ln=1)
       Call Out ("with total cluster weight",tot_wt)
       Call Out ("from cluster weights",siz(:ind))
       If (pr_out > 1.5) Call Out (-1,"and mean vectors",clust_vec(1:,:ind))
     End if
      
!    Degree graph peak clustering

     q= q + 1;  n1= 0
     Call DeAlloc_Adjacency (Slate_G);  Slate_G%sum_wt= Memb(1:)%fsx
     Do j= 1,4
       Call Slate_graph (ns,cut(j), Memb,Slate_G, sum_deg)
       If (sum_deg < 1) Cycle
       Call Degree_graph (Clust_opt, Slate_G, Coef, max_deg)
       If (max_deg < 0.001) Cycle
       Call Peak_clustering (ns,Slate_G, 0.05,Coef, n1,peak,pk_val)
       If (n1 > np) Exit
     End do

     If (n1 > 0) then
       n1= Min(n1,nr);  ind= n1 + 1;  cnt= 0;  siz=0

       Call Clust_matrix (q,nc,np,ns,n1,Noise_cor,Memb(1:), clust_vec(:,:ind), &
                          cnt(:ind),siz(:ind),clust_memb(:,:ind), peak(:n1))

       Allocate(Init(q)%T0(0:nc,0:ind,2));  Init(q)%T0= 0

       Init(q)%T0(:,1:,1)= clust_vec(:,:ind)
       Forall(c=1:nc) Init(q)%T0(c,0,1)= Sum(clust_vec(0,:ind) * clust_vec(c,:ind)) / np

       tot_wt= Sum(siz(:ind)) 

       If (pr_out > 1) then
         Call Out ("For the degree graph peak clustering set",q, &
                   "# regular clusters",n1, ln=1)
         Call Out ("with total cluster weight",tot_wt)
         Call Out ("from cluster weights",siz(:ind))
         Call Out ("Peak coefficient values",pk_val(:n1))
         If (pr_out > 1.5) Call Out (-1,"and mean vectors",clust_vec(1:,:ind))
       End if
     End if

     Init(q)%n= n1
     
!    Compute and store the trianglar coefficient initial cluster set

     q= q + 1;  n1= 0;  peak= 0
     Call DeAlloc_Adjacency (Slate_G);  Slate_G%sum_wt= Memb(1:)%fsx

     Do j= 1,4
       Call Slate_graph (ns,cut(j), Memb,Slate_G, sum_deg)
       If (sum_deg < 1) Cycle
       Call Triangle_graph (Clust_opt, Slate_G, Coef, max_deg)
       If (max_deg < 0.001) Cycle
       Call Peak_clustering (ns,Slate_G, 0.05,Coef, n1,peak,pk_val)
       If (n1 > np+1) Exit
     End do  
     
     If (n1 > 0) then
       n1= Min(n1,nr);  ind= n1 + 1;  cnt= 0;  siz=0

       Call Clust_matrix (q,nc,np,ns,n1,Noise_cor,Memb(1:), clust_vec(:,:ind), &
                          cnt(:ind),siz(:ind), clust_memb(:,:ind), peak(:n1))

       Allocate(Init(q)%T0(0:nc,0:ind,2));  Init(q)%T0= 0

       Init(q)%T0(:,1:,1)= clust_vec(:,:ind)
       Forall(c=1:nc) Init(q)%T0(c,0,1)= Sum(clust_vec(0,:ind) * clust_vec(c,:ind)) / np

       tot_wt= Sum(siz(:ind))

       If (pr_out > 1) then
         Call Out ("For the triangular graph peak clustering set",q, "# clusters",n1, ln=1)
         Call Out ("with total cluster weight",tot_wt)
         Call Out ("from cluster weights",siz(:ind))
         Call Out ("Peak coefficient values",pk_val(:n1))
         If (pr_out > 1.5) Call Out (-1,"and mean vectors",clust_vec(1:,:ind))
       End if
     End if

     Init(q)%n= n1


!    STV cluster set
    
     q= q + 1;  ind= Last_true(mean_SD(0,:,1) > 0);  n1= ind - 1

     Allocate(Init(q)%T0(0:nc,0:ind,2));  Init(q)%T0= 0

     cnt= 0;  siz=0;  Init(q)%n= n1
     Init(q)%T0(:,1:,1)= mean_SD(:,:ind,1)

     Call Clust_matrix (q,nc,np,ns,n1,Noise_cor,Memb(1:), &
                        Init(q)%T0(:,1:,1), cnt(:ind),siz(:ind), clust_memb)

     Forall(c=1:nc) Init(q)%T0(c,0,1)= Sum(siz(:ind) * Init(q)%T0(c,1:,1)) / np
     tot_wt= Sum(siz(:ind))

!    DTV cluster set
    
     q= q + 1;  ind= Last_true(mean_SD(0,:,2) > 0);  n1= ind - 1

     Allocate(Init(q)%T0(0:nc,0:ind,2));  Init(q)%T0= 0

     cnt= 0;  siz=0;  Init(q)%n= n1
     Init(q)%T0(:,1:,1)= mean_SD(:,:ind,2)

     Call Clust_matrix (q,nc,np,ns,n1,Noise_cor,Memb(1:), &
                        Init(q)%T0(:,1:,1), cnt(:ind),siz(:ind), clust_memb)

     Forall(c=1:nc) Init(q)%T0(c,0,1)= Sum(siz(:ind) * Init(q)%T0(c,1:,1)) / np
     tot_wt= Sum(siz(:ind))

     If (pr_out > 1) Call Out ("Initial cluster sets, elect candidates by Rem_frac, DHondt, & optimal methods")

     Init_loop : Do iq= 1,N_init
       ncl= Init(iq)%n;  ind= ncl+1;  If (ncl < 1) Cycle Init_loop

       Allocate(Init(iq)%lt(ncl), Init(iq)%Q0(np,0:2,4), &
                Init(iq)%M1(7,4), Init(iq)%M2(ncl,ncl))
       Init(iq)%lt= "ID";  Init(iq)%Q0= 0;  Init(iq)%M1= 0;  Init(iq)%M2= 0  

       If (iq >= 19) then
         Init(iq)%Q0(:,:,1)= elc_SD(:,:,iq-18)
       Else
         Call Rem_frac (iq,nc,ind, Init(iq)%T0(:,1:,1), Init(iq)%Q0(:,:,1))
       End if
       Call DHondt (iq,nc,ind, Init(iq)%T0(:,1:,1), Init(iq)%Q0(:,:,2))

       Call Elect_fit (iq,np,nc,ind, Memb(1:),Init(iq)%T0(:,:,1), &
                       Init(iq)%Q0(:,0,:2), Init(iq)%Q0(:,:,3),   &
                       Init(iq)%T0(1:,:,2),Init(iq)%M2,Init(iq)%M1(:,:3)) 

     End do Init_loop
    
   End Subroutine Initialize_clusters


  Subroutine Clustering (idist,Ninit,nc,Memb, Init, nd,Clust_set)

!   Clustering for implicit proportional representation, given that each
!    voter ranks or rates one or more of 'nc' candidates for np positions.

! Input:
    Integer,              Intent(in) :: idist       ! Index of the voting district
    Integer,              Intent(in) :: Ninit       ! # initial sets
    Integer,              Intent(in) :: nc          ! # candidates
    Type(Multi_listD),    Intent(in) :: Memb(0:)    ! (0:ns) Final slate cluster data
                                                    ! 0 Case: 
                                                    !   %n = nc = # candidates 
                                                    !   %l = np = # candidates to be elected = total weight
                                                    !             in position units
                                                    !   %qx(2)  = correlation parameters
                                                    !   %vl(mlv). (lv) = # slate clusters with
                                                    !      # ranked or rated candidates <= lv

                                                    ! 1:nsl Case:
                                                    !  %l = # top slate candidates
                                                    !  %m = # slate candidates (for Rating = 1)
                                                    !  %p = # positively rated candidates (for %tx)
                                                    !  %n = # negatively rated candidates (for %tx)
                                                    !  %k = # slate ballots positively correlated with and 
                                                    !       subsequent to this slate ballot
                                                    !  %vl = listing of these by decreasing correlation
                                                    !  %qx = the corresponding correlations

                                                    !  %fsx= slate membership
                                                    !  %fux= slate width
                                                    !  %lt(nc) = candidates in preferential order (from %sx)
                                                    !  %px(nc) = variance vector for %sx
                                                    !  %rx(nc) = mean rating vector with variance normalized weight
                                                    !          = %sx * %ux
                                                    !  %sx(nc) = mean rating vector
                                                    !  %tx(0:nc)= centered (subtract Noise_cor if Rating < 1) 
                                                    !             and normalized mean vector
                                                    !  %ux(nc) = slate ballot weight as variance normalized 
                                                    !          = %fsx / %px

    Type(Multi_listR), Intent(inout) :: Init(:)  ! (Ninit) Initial cluster set data
                                                 ! %k = cluster set converged to
                                                 ! %n = ncl = # regular clusters

                                                 ! %Q0(np,0:2,4) Elected candidates, by 4 methods 
                                                 !      (:,:,:1) = Rem_frac electeds (original for STV/DTV)
                                                 !      (:,:,:2) = DHondt electeds
                                                 !      (:,:,:3) = Optimal electeds 
                                                 !      (:,:,:4) = Converged optimal electeds
                                                 !     (:,0,:) = electeds in increasing order
                                                 !     (:,1,:) = electeds in ranked order
                                                 !     (:,2,:) = corresponding clusters

                                                 ! %M0(ns,ind) Cluster slate cluster membership weights
                                                 ! %M1(7,4) Objective data for the electeds $Q0(:,0,:)
                                                 ! %M2(ncl,ncl) Cluster correlation matrix, 0 on diagonal
 
                                                 ! %T0(0:nc,0:ind,2) (1:,1:,1) = Candidate mean ratings  
                                                 !      for each cluster, with cluster size at (0,1:,1)
                                                 !      (1:,1:,2) = Cluster portions with averaged 
                                                 !                  noise zeroed ratings at (1:,0,2)

    Integer,           Intent(out) :: nd           ! # final cluster sets in Clust_set

    Type(Multi_listR), Intent(out) :: Clust_set(:) ! (Ninit) Final cluster sets and associated convergence data
                                                   ! %k   = best initial cluster set that converged to this cluster set
                                                   ! %l   = # clusters in the set, including independents = 'ind' = ncl + 1
                                                   ! %m   = # initial cluster sets that converged to this cluster set
                                                   ! %ls(m)= list of those initial sets 
                                                   ! %fux = clustering objective ratio = (q)%rx(1) / (1)%rx(1)

                                                   ! %lt(ncl) Mapping of initial to converged regular clusters 

                                                   ! %px(6):  Objective data: min size, max correlation, 
                                                   !            independents size (frac), final residual,
                                                   !            Jnorm, Pnorm
                                                   ! %qx(4):  Objective factors: min size penalty, max correlation penalty, 
                                                   !            independents size penalty, residual penalty
                                                   ! %rx(3):  Clustering objective value(1)= coherent vote(2) * penalty(3)
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

                                                   ! %L2(n0,2) Mapping data between initial and final clusters 
                                                   !     1 : Maps initial clusters onto intermediate clusters, reordered
                                                   !     2 : For each final cluster, the largest initial cluster mapped to it

                                                   ! %M2(0:2,0:n1) Real convergence data
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

                                                   !   %L1(-2:nc,0:ind) Candidate ordering data for
                                                   !                    noise zeroed vectors
                                                   !     (-2:0,:) = # top, significant, viable cand's
                                                   !     (1:,1:)  = Decreasing ordering from %T0(:,1:,2)
                                                   !     (1:,0)   = Decreasing ordering from %T0(:,0,1)
                                                   !               
                                                   !   %T0(nc,0:ind,3) Cluster vector data
                                                   !      (:,1:,1) = Cluster mean vectors
                                                   !      (:,1:,2) = Cluster mean vectors with noise zeroing 
                                                   !      (:,1:,3) = Cluster portion vectors
                                                   !      (:,0,1)  = Zrate0 * declining regular portion factor
                                                   !      (:,0,2)  = Zrate0 = cluster averaged noise zeroed mean vectors

! Local:
    Integer, Parameter :: nrd= 4   ! # data components for pr_dat
    Integer :: pr_ncl(3,Ninit)     ! 1 = # regular clusters in each 
                                   !     recorded set
                                   ! 2 = # matches from all converged 
                                   !     cluster sets to each recorded set
                                   ! 3 = average objective * 100
    Real    :: pr_dat(nrd,Ninit)   ! Data for each for each recorded set
                                   !   1 = coherent vote
                                   !   2 = min cluster size
                                   !   3 = max cluster correlation
                                   !   4 = independents size
    
    Type(Multi_listR), Allocatable :: Clust_tmp(:) 
    Real    :: tmp(Ninit)
    Integer :: lst(Size(Memb)), lq(Ninit), ls(Ninit), cls(Ninit), key(Ninit), inv(0:Ninit)
    Logical :: ReOrd
    Integer :: i, m, q, iq, ns, q0,  nd0

    Call Out ("Enter 'Clustering' to converge clusters from district",idist, ln=1)

    ns= Ubound(Memb,1)  

    If (Allocated(LS_Memb)) then
      Call DeAlloc_Multi_list_ar (LS_Memb);  DeAllocate(LS_Memb)
    End if
    Allocate(LS_Memb(0:ns));  lst= "ID"
    Call Copy_Multi_List (lst,Memb, LS_Memb)

    Call Converge_multiple (idist,Ninit,Parm_cnv,Parm_obj, Memb,Init, &
                            nd,Clust_set, pr_ncl,pr_dat)
    If (nd < 1) then
      Call Out ("Error in 'Clustering': No cluster sets")
      Stop
    End if

!   Remove any final cluster sets which converge from only 1 initial cluster set

    nd0= nd;  iq= 0;  lq= 0

    Do q= 1,nd
      m= Count(Init%k == q)
      If (m > 1) then  ! Keep this cluster set
        iq= iq + 1;  lq(iq)= q
      End if
    End do

    nd= iq;   ReOrd= .false.

!   Reorder Clust_set by decreasing clustering objective, with any removals

    If (nd > 1) then
      Do iq= 1,nd
        q= lq(iq);  tmp(iq)= Clust_set(q)%rx(1)
      End do
      tmp(nd+1:)= 0;  key(nd+1:)= 0

      Call Sort (.false., tmp(:nd), key(:nd),ReOrd)

      If (ReOrd) then
        key(:nd)= lq(key(:nd));  Allocate(Clust_tmp(nd));  ls(:nd)= "ID"

        Call Copy_multi_list (key(:nd),Clust_set(:nd0), Clust_tmp)
        Call Copy_multi_list (ls(:nd),Clust_tmp, Clust_set(:nd))

        pr_ncl(:,:nd)= pr_ncl(:,key(:nd))
        pr_dat(:,:nd)= pr_dat(:,key(:nd))

        inv= 0;  Call Inverse_map (key(:nd),inv(1:))
        cls= Init%k;  Init%k= inv(cls)

        If (pr_out >= 1) then
          Call Out ("Cluster sets to keep",lq(:nd))
          Call Out ("Cluster sets reordered using the key",key(:nd))
          Call Out ("Inverse key",inv(:nd0))
          Call Out ("Original mapping of initial to final cluster sets",cls)
          cls= Init%k;  Call Out ("Revised mapping",cls)
        End if

        Call DeAlloc_multi_list_ar (Clust_tmp)
      End if
    End if

    If (nd < nd0 .and..not.ReOrd) then
      cls= 0
      Do iq= 1,nd
        q= lq(iq);  Where(Init%k == q) cls= iq
        If (q == iq) Cycle

        Call Copy_multi_list (Clust_set(q), Clust_set(iq))
        pr_ncl(:,iq)= pr_ncl(:,q)
        pr_dat(:,iq)= pr_dat(:,q)
      End do

      Init%k= cls;  pr_ncl(:,nd+1:)= 0;  pr_dat(:,nd+1:)= 0

      If (pr_out >= 1) then
        Call Out ("Cluster sets to keep",lq(:nd))
        Call Out ("Revised mapping of initial to final cluster sets",cls)
      End if
    End if

    Do q= 1,nd
      Clust_set(q)%fux= Clust_set(q)%rx(1) / Clust_set(1)%rx(1) 

      Call List_of_true (Init%k == q, m,ls)
      Allocate (Clust_set(q)%ls(m))  
      Clust_set(q)%m= m;  Clust_set(q)%ls= ls(:m)

      Forall(i=1:m) tmp(i)= Init(ls(i))%M1(1,3)
      iq= Maxloc(tmp(:m),1);  Clust_set(q)%k= iq
    End do

    Do i= 1,nct
      cls_n(i)= Last_true(Clust_set(:nd)%fux >= cls_cut(i))
    End do

    If (pr_out > 1) then
        m= Sum(pr_ncl(2,:nd))
      Call Out ("After reordering: # cluster sets matched",nd, "vs total # matches",m)
      Call Out (-1,"# clusters, # matches, objective", pr_ncl(:,:nd))
      Call Out (-1,"normalized vote, min size, max corr, ind size", pr_dat(:,:nd))
    End if

  End Subroutine Clustering


  Subroutine Converge_multiple (idist,Ninit,Parm_cnv,Parm_obj, Memb,Init, &
                                nd,Clust_set, pr_ncl,pr_dat)

!   Converge from initial cluster sets 'Init' to final clusters.
  
    Integer, Intent(in) :: idist  ! Index of the voting district
    Integer, Intent(in) :: Ninit  ! # initial cluster sets

    Real,    Intent(in) :: Parm_cnv(:,:) ! Convergence parameters
                                 ! For convergence phase RS = 1, 2, 3:
                                 ! (RS,1) = convergence tolerances
                                 ! (RS,2) = min merge criterion for 
                                 !   cluster correlation
                                 ! (RS,3) = min size criterion for cluster
                                 !   deletion
                                 ! (RS,4) = max correlation criterion for
                                 !   cluster deletion
                                 ! (1,5) = maximum norm permitted for an
                                 !   update differential
                                 ! (2,5) = match tolerance for measures of
                                 !   cluster set similarity
                                 ! (3,5) = min separation criterion for the
                                 !   separation graph for the clique of
                                 !   correlated clusters 
    Real,    Intent(in) :: Parm_obj(:)  ! (12) Parameters for the penalty factors 
                                     !        of the objective function
                                     !  1:2 = cluster density soft cutoff limits (cos_rise)
                                     !  3:4 = cluster boundary soft cutoff limits (cos_fall)
                                     !  5:6 = cluster size soft cutoff limits (cos_rise)
                                     !  7:8 = cluster correlation soft cutoff limits (cos_fall)
                                     !  9:10 = independents size cutoff limits (cos_fall)
                                     !  11:12 = convergence residual limits (cos_fall)
    
    Type(Multi_listD),    Intent(in) :: Memb(0:) ! (0:ns) Slate cluster data

    Type(Multi_listR), Intent(inout) :: Init(:)  ! (Ninit) Initial cluster set data
                                                 ! %k = cluster set converged to (output):
                                                 !      0 if convergence failed, -1 if severe error
                                                 ! %n = ncl = # regular clusters
                                                 ! %ls(-2:nc) Ordering of %T0(1:,0,2), with 
                                                 !      # top, significant, and viable at (-2:0)

                                                 ! %Q0(np,0:2,4) Elected candidates, by 4 methods 
                                                 !      (:,:,:1) = Rem_frac electeds (original for STV/DTV)
                                                 !      (:,:,:2) = DHondt electeds
                                                 !      (:,:,:3) = Optimal electeds 
                                                 !      (:,:,:4) = Converged optimal electeds
                                                 !     (:,0,:) = electeds in increasing order
                                                 !     (:,1,:) = electeds in ranked order
                                                 !     (:,2,:) = corresponding clusters

                                                 ! %M0(ns,ind) Cluster slate cluster membership weights
                                                 ! %M1(7,4) Objective data for the electeds %Q0(:,0,:)
                                                 ! %M2(ncl,ncl) Cluster correlation matrix, 0 on diagonal
 
                                                 ! %T0(0:nc,0:ind,2) (1:,1:,1) = Candidate mean ratings  
                                                 !      for each cluster, with cluster size at (0,1:,1)
                                                 !      (1:,1:,2) = Cluster portions with averaged 
                                                 !                  noise zeroed ratings at (1:,0,2)

    Integer,           Intent(out) :: nd           ! Current # distinct, converged cluster sets
    Type(Multi_listR), Intent(out) :: Clust_set(:) ! (Ninit => ncs) Final cluster sets and associated convergence data
                                                   ! %k   = best initial cluster set that converged to this cluster set
                                                   ! %l   = # clusters in the set, including independents = 'ind' = ncl + 1
                                                   ! %m   = # initial cluster sets that converged to this cluster set
                                                   ! %ls(m)= list of those initial sets 
                                                   ! %fux = clustering objective ratio = (q)%rx(1) / (1)%rx(1) 

                                                   ! %lt(ncl) Mapping of initial to converged regular clusters 

                                                   ! %px(6):  Objective data: min size, max correlation, 
                                                   !            independents size (frac), final residual,
                                                   !            Jnorm, Pnorm
                                                   ! %qx(4):  Objective factors: min size penalty, max correlation penalty, 
                                                   !            independents size penalty, residual penalty
                                                   ! %rx(3):  Clustering objective value(1)= coherent vote(2) * penalty(3)
                                                   ! %sx(ind): Cluster sizes

                                                   ! %ux(6)  Fraction of voters represented to a given level
                                                   !         or higher with levels specified by Parm_rep

                                                   !  %vl(0:6): Convergence data
                                                   !     0 = original # regular clusters = n0
                                                   !     1 = final # regular clusters
                                                   !     2 = final success code
                                                   !     3 = total # cluster merge and deletion operations
                                                   !     4 = # convergence calls
                                                   !     5 = total # centroiding iterations
                                                   !     6 = total # line search points

                                                   !  %L0(0:4,0:iv) Integer convergence data per update call, n1= %vl(4)
                                                   !    0  = restart phase
                                                   !    1  = success code:  2 = fully converged,   1 = reducing functional,
                                                   !                        0 = damped reduction, -1 = convergence failure
                                                   !    2 = # reg clusters
                                                   !    3 = # centroiding iterations
                                                   !    4 = # line search points

                                                   ! %M2(0:2,0:n1) Real convergence data
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

                                                   !   %L1(-2:nc,0:ind) Candidate ordering data for
                                                   !                    noise zeroed vectors
                                                   !     (-2:0,:) = # top, significant, viable cand's
                                                   !     (1:,1:)  = Decreasing ordering from %T0(:,1:,2)
                                                   !     (1:,0)   = Decreasing ordering from %T0(:,0,1)
                                                   !               
                                                   !   %T0(nc,0:ind,3) Cluster vector data
                                                   !      (:,1:,1) = Cluster mean vectors
                                                   !      (:,1:,2) = Cluster mean vectors with noise zeroing 
                                                   !      (:,1:,3) = Cluster portion vectors
                                                   !      (:,0,1)  = Zrate0 * declining regular portion factor
                                                   !      (:,0,2)  = Zrate0 = cluster averaged noise zeroed mean vectors
      
    Integer, Intent(out) :: pr_ncl(:,:)    ! (3,Ninit) 1 = # regular clusters for each converged cluster set
                                           !            2 = # matches, itself included, to it from subsequent 
                                           !                converged cluster sets.
                                           !            3 = mean objective * 100
    Real,    Intent(out) :: pr_dat(:,:)    ! (nrd,Ninit) Data for converged cluster sets (up to nd)
                                           !   1 = coherent vote
                                           !   2 = min cluster size
                                           !   3 = max cluster correlation
                                           !   6 = independents size
! Local:    
    Type(Multi_listD), Allocatable :: Clust(:)
    
    Real    :: px(4)
    Integer :: ier, ind, ncl, rep, ncyc, match, iobj, iset, lst(Ninit)
    Integer :: k, n, q, iq, cl, mx, n0, nc, np, ns, nw
    
    Call Out ("Enter 'Converge_multiple' for district #",idist, ln=1)

    ns= Ubound(Memb,1);  nc= Memb(0)%n;  np= Memb(0)%l
    mx= Maxval(Init%n);  ind= mx + 1
    Allocate (Clust(0:ind))

    nd= 0;  pr_ncl= 0;  pr_dat= 0;  Init%k= 0
    
    Initial_clusters : Do iq= 1,Ninit
      n0= Init(iq)%n;  ind= n0 + 1;  If (n0 < 1) Cycle Initial_clusters

      Call Clusters_from_means (nc,n0,iq, Memb, Dble(Init(iq)%T0(1:,1:n0,1)), Clust(:ind))
      ncl= Clust(0)%n

      Call Iterate (idist,iq,nc,n0, Mxitr,Parm_cnv,Parm_obj, &
                    Memb,Clust(:ind), ncl)
      
      ier= Clust(0)%k  

      If (ier < 0) then
        Call Out ("Error in 'Converge_multiple' for district",idist, "& initial cluster set",iq, ln=1)
        Init(iq)%k= 0;  Cycle Initial_clusters
      Else 
         If (ier == 0) Call Out ("Warning in 'Converge_multiple': partial convergence for district",idist, &
                                 "& initial cluster set",iq, ln=1)
         iobj= Floor(100*Clust(0)%rx(1))
         If (iobj < 1) then
           Call Out ("Warning in 'Converge_multiple': bad objective for district",idist, &
                                 "& initial cluster set",iq, ln=1)
           Init(iq)%k= 0;  Cycle Initial_clusters
         End if
      End if

      Init(iq)%lt= Clust(0)%lt

      If (pr_out > 1) then
        Call Out ("For for district #",idist, "initial set",iq, ln=1)
        Call Out ("# clusters reduced from",n0, "to",ncl)
        Call Out ("with scaled full objective",iobj)
        Call Out ("Final clusters combined from initial clusters",Init(iq)%lt)
        Call Out ("of sizes",Real(Clust(0)%sx(:ncl)))
        Call Out ("with convergence data",Clust(0)%vl)
      End if  
      
!     Test for matching to a prior set, first normalizing the data

      px(1)= Clust(0)%rx(2) / GN_max_pt;  px(2:)= Clust(0)%px(:3)
      px(3)= Max(px(3),0.0) + 0.3  ! Max correlation = special case
      px= px + 0.20
      
      Call Cluster_match1 (iq,nc,ncl, iobj,px, nd,pr_ncl,pr_dat, match)

!     Record a new cluster set in Clust_set if feasible.
!     'match' = prior matching cluster set if > 0, or was added
!      to the end of the list at position 'nd = -match' if < 0
      
      If (match < 0) then  ! New converged cluster set
        Call Record_clusters (nd,ncl, Memb, Clust(:ncl+1), Clust_set)

        Init(iq)%k= nd
      Else if (match > 0) then  ! Matches prior cluster set 'match'
        Init(iq)%k= match
      End if
    End do Initial_clusters

    Call DeAlloc_Multi_list_ar (Clust);  DeAllocate (Clust)

    lst= Init%k
    Call List_of_true (Init%k < 1, n,lst);  lst(n+1:)= 0

    If (n > 0) then
      Call Out ("For district",idist, "# convergence failures",n, ln=1)
      If (pr_out >= 1) then
        Call Out ("convergence failures",lst(:n));  lst= Init%k
        Call Out ("Final cluster sets converged to from initial",lst)
      End if
    End if

    If (pr_out > 1) then
      Call Out ("# distinct cluster sets recorded by 'Converge_multiple'",nd, ln=1)
      Call Out (-1,"Before reordering: # clust, # matches, obj",pr_ncl(:,:nd))
      Call Out (-1,"Normalized data: coh vote, min size, max corr, ind size",pr_dat(:,:nd))
      Do q= 1,nd
        Call Out ("Cluster sizes for each set",Real(Clust_set(q)%sx))
        Call Out ("Corresponding cluster objectives",Real(Clust_set(q)%rx))
      End do
    End if
    
  End Subroutine Converge_multiple
                                

    Subroutine Clust_matrix (q,nc,np,ns,ncl,Noise_cor,Memb, clust_vec, &
                             cnt,siz, clust_memb, slate_peaks)

!     Given the mean vectors of the regular clusters 'clust_vec(1:,:ncl)',
!     or, alternatively, 'slate_peaks' whose slate cluster data from 'Memb'
!     substitues for 'clust_vec', compute the slate_cluster memberships
!     'clust_memb' of each cluster, include independents.

!     Next, reduce the 'clust_memb' values, as necessary, to meet the 
!     fuzzy set requirement: For each slate cluster from 'Memb', the  
!     sum 'wtsl' of its memberships in the regular clusters must be <= 1.0,
!     with any excess 1.0 - wtsl going assigned to the independents.

!     Finally, recompute the cluster mean vectors 'clust_vec(1:,:)' 
!     and the clust weights 'clust_vec(0,:)' from the scaled memberships,
!     reordering them by decreasing weight.

      Integer,            Intent(in) :: q                ! Called from this initial cluster set
      Integer,            Intent(in) :: nc               ! # candidates
      Integer,            Intent(in) :: np               ! # candidates to be elected = total cluster set weight
      Integer,            Intent(in) :: ns               ! # slate ballots
      Integer,            Intent(in) :: ncl              ! # regular clusters

      Real,               Intent(in) :: Noise_cor        ! Noise adjustment level for the ratings
      Type(Multi_listD),  Intent(in) :: Memb(:)          ! (ns) Slate cluster data
    
      Real,            Intent(inout) :: clust_vec(0:,:)  ! (0:nc,ind) Mean rating vectors of the clusters,
                                                         !   with sizes at (0,:), ordered by decreasing
                                                         !   size, with sizes summing to 'np' (in and out)
      
      Integer,           Intent(out) :: cnt(:)           ! (ind) Slate ballot membership count for each cluster
      Real,              Intent(out) :: siz(:)           ! (ind) Cluster sizes = clust_vec(0,:)
      Real,              Intent(out) :: clust_memb(:,:)  ! (ns,ind) Slate ballot memberships for each cluster
      Integer,  Optional, Intent(in) :: slate_peaks(:)   ! (ncl) List of slate ballots whose mean vectors 
!   Local:
      Real, Parameter :: corr_cut= 0.01  ! Cluster membership correlation cutoff level
      Real, Parameter :: eps= 0.000001
      Real     :: siz0(ncl), tx(nc,ncl), mbr0(ns), mbr(ns)
      Integer  :: key(ncl)
      Logical  :: ReOrd
      Real     :: wt, mbi, norm, corr, tot_wt
      Integer  :: n, cl, sl, ind
      
      Call Out ("Enter 'Clust_matrix' from initial cluster set",q, ln=1)
      ind= ncl + 1

!     Initialize 'clust_vec' with 'slate_peaks' if present.
!     Otherwise use 'clust_vec as input. Also compute 
!     the noise adjusted and normalized version 'tx'
!     of 'clust_vec'

      If (Present(slate_peaks)) then
        clust_vec= 0
        Do cl= 1,ncl
            sl= slate_peaks(cl)
          clust_vec(0,cl) = Memb(sl)%fsx
          clust_vec(1:,cl)= Memb(sl)%sx
          tx(:,cl)= Memb(sl)%tx(1:)
        End do

        Call Sort (.false.,clust_vec(0,:ncl),key)
        clust_vec(1:,:ncl)= clust_vec(1:,key)
        tx= tx(:,key)
      Else
        clust_vec(:,ind)= 0;  siz0= clust_vec(0,:ncl)
        If (Rating < 1) then
          tx= clust_vec(1:,:ncl) - Noise_cor
        Else
          tx= clust_vec(1:,:ncl)
        End if

        Do cl= 1,ncl 
          norm= Sqrt(Sum(tx(:,cl)**2))
          tx(:,cl)= tx(:,cl) / norm
        End do
      End if

!     The fuzzy set computation

      clust_memb= 0;  mbr0= 0

      Do sl= 1,ns
        Do cl= 1,ncl
          Call Dot_product_M (Real(Dot_fac), Real(Memb(sl)%tx(1:)), tx(:,cl), n,corr)

          If (corr > corr_cut) then
            mbr0(sl)= mbr0(sl) + corr
            clust_memb(sl,cl)= corr * Memb(sl)%fsx  
          End if
        End do

        If (mbr0(sl) > 1) clust_memb(sl,:ncl)= clust_memb(sl,:ncl) / mbr0(sl)
      End do

      Where (mbr0 > 1)
        mbr= 1 
      Else where
        mbr= mbr0
      End where

!     Recompute the regular cluster mean vectors 'clust_vec'

      siz= 0;  clust_vec= 0

      Do cl= 1,ncl
        Do sl= 1,ns
          wt= clust_memb(sl,cl)
          If (wt > 0) then
            siz(cl)= siz(cl) + wt
            clust_vec(1:,cl)= clust_vec(1:,cl) + Memb(sl)%sx * wt
          End if
        End do
        clust_vec(1:,cl)= clust_vec(1:,cl) / Max(siz(cl),eps)
      End do

!     Compute the new independents 'clust_vec(:,ind)'

      Do sl= 1,ns
        mbi= 1 - mbr(sl);  If (mbi <= 0) Cycle

        wt= mbi * Memb(sl)%fsx;  clust_memb(sl,ind)= wt
        siz(ind)= siz(ind) + wt
        clust_vec(1:,ind)= clust_vec(1:,ind) + Memb(sl)%sx * wt
      End do

      clust_vec(1:,ind)= clust_vec(1:,ind) / Max(siz(ind), eps)
      clust_vec(0,:)= siz;  tot_wt= Sum(siz)

!     Reorder the regular clusters as necessary

      Call Sort (.false.,siz(:ncl),key, ReOrd)

      If (ReOrd) then
        clust_memb(:,:ncl)= clust_memb(:,key)
        clust_vec(:,:ncl) = clust_vec(:,key)
      End if

      Forall(cl=1:ind) cnt(cl)= Count(clust_memb(:,cl) > 0)

      If (pr_out > 1) then
        Call Out ("In 'Clust_matrix': total weight",tot_wt, ln=1)
        Call Out ("Original cluster sizes",siz0)
        Call Out ("Revised and reordered cluster sizes",siz)
        Call Out ("by key",key)
        Call Out ("Cluster member counts",cnt)
        Call Out ("Original slate cluster membership weights",mbr0)
        Call Out (-1,"Revised cluster mean vectors",clust_vec(1:,:))
      End if

    End Subroutine Clust_matrix

    
    Subroutine Clust_mean (nc,np, Memb,clust_memb, siz,clust_vec)

!     Compute the cluster mean vectors 'clust_vec' from their 
!     slate ballot memberships 'clust_memb', assuming that the
!     cluster memberships sum to 'np'.

!     If necessary, reorder the regular clusters by decreasing size
    
      Integer,           Intent(in) :: nc               ! # candidates
      Integer,           Intent(in) :: np               ! # candidates to be elected = total cluster set weight
      Type(Multi_listD), Intent(in) :: Memb(:)          ! (ns) Slate ballot data

      Real,           Intent(inout) :: clust_memb(:,:)  ! (ns,ind) Membership weight of each slate ballot
                                                        !   in each cluster, with regular clusters by decreasing size
      Real,             Intent(out) :: siz(:)           ! (ind) Cluster sizes = clust_vec(0,:)
      Real,             Intent(out) :: clust_vec(0:,:)  ! (0:nc,ind) Mean rating vectors of the clusters,
                                                        !   with clusters 'cl' in decreasing order of size 
      
!   Local:
      Real, Parameter :: eps0= 0.000001
      Integer, Allocatable :: key(:)
      Logical :: ReOrd
      Real    :: wt, fac, tot_wt
      Integer :: cl, sl, ns, ind, ncl

      ind= Size(clust_memb,2);  ncl= ind - 1;  ns= Ubound(clust_memb,1)
      siz= 0;  clust_vec= 0

      Do cl= 1,ind
        Do sl= 1,ns
          wt= clust_memb(sl,cl)

          If (wt > eps0) then 
            clust_vec(1:,cl)= clust_vec(1:,cl) + Memb(sl)%sx * wt
            siz(cl)= siz(cl) + wt
          End if
        End do

        clust_vec(1:,cl)= clust_vec(1:,cl) / Max(siz(cl),eps0)
        clust_vec(0,cl)= siz(cl) 
      End do

!     Check total cluster weight

      tot_wt= Sum(siz);  fac= np / tot_wt
      If (fac /= 1) then
        siz= fac * siz;  clust_vec(0,:)= siz
        clust_memb= fac * clust_memb
      End if

!     Reorder by size as needed

      Allocate(key(ncl));  Call Sort (.false.,siz(:ncl),key,ReOrd)
      If (ReOrd) then
        clust_vec(:,:ncl) = clust_vec(:,key)
        clust_memb(:,:ncl)= clust_memb(:,key)
      End if

      If (pr_out > 1) then
        Call Out ("In 'Clust_mean': total cluster weight",tot_wt, ln=1)
        Call Out ("from cluster sizes",siz)
        Call Out (-1,"Cluster mean vectors",clust_vec(1:,:))
      End if
    End Subroutine Clust_mean

    
    Subroutine Core_clustering (Core_parm, mc,nv, G, ncl,Core)
  
!     Cluster the vertices of a graph 'G' by finding 'ncl' sets of
!     non-intersecting vertices 'Core', each of which forms the core 
!     of a cluster. Non-core vertices may be members of more than 
!     one cluster, as determined by the graph 'G', subject to the 
!     fuzzy membership rule: The sum of the memberships of each 
!     such vertex must be <= 1. Not fully assigned vertices will
!     have their remainder treated as independents.
    
      Real,            Intent(in) :: Core_parm(:) ! (3): Core_clustering parameters
                                                  !  1 = correlation cutoff for cluster separation
                                                  !  2 = correlation cutoff for a core set                                              
                                                  !  3 = total membership cutoff for small clusters
      Integer,         Intent(in) :: mc           ! max # core clusters < Size(Core) to allow for independents
      Integer,         Intent(in) :: nv           ! # vertices
      Type(Adjacency), Intent(in) :: G(:)         ! (nv) Weighted graph, with vertices ordered by decreasing weight
                                                  !  %nl      = # adjacent vertices
                                                  !  %sum_wt  = vertex weight
                                                  !  %ls(nl)  = List of adjacent vertices
                                                  !  %wt(nl)  = Their edge weights (0 to 1) = similarities, 
                                                  !             ordered by decreasing similarity
      
      Integer,        Intent(out) :: ncl          ! # clusters found = # core sets
      Type(Set_list), Intent(out) :: Core(:)      ! (mc+1) The core sets of the clusters
                                                  !  %n      = # vertices in the core set
                                                  !  %set(n) = list of those vertices
                                                  !  %mbr(nv)= 0..1 memberships of all vertices in the core
                                                  !  %svl    = sum of the core vertex weights
                                                  !  %smb    = total weighted membership of the cluster
!   Local:
      Type(Set_list), Allocatable :: CoreP(:)  ! (ncl)
      Integer, Allocatable :: set(:,:)
      Integer :: Clustered(nv)  ! Core members: Clustered(v) = j if v is a core member 
                                ! of cluster j, else Clustered(v) = -1 
      Integer :: ncore(mc)      ! # current members of the core clusters
      Real    :: vol(mc)        ! Sum of the 'G' vertex weights of the core members

      Logical :: ReOrd, Bound_v(nv), Bound_i(nv)
      Real    :: sm, sim(mc), mbr(nv), smbr(nv)
      Integer :: ls(nv)
      Integer :: i, j, k, l, m, n, v

      Call Out ("Enter 'Core_clustering'")

      ls= -1;  sim= -1
      Clustered= -1;  ncore= 0;  vol= 0
      Clustered(1)= 1;  ncl= 1;  ncore(1)= 1;  vol(1)= G(1)%sum_wt
      
      Main_loop : Do v= 2,nv
        k= 0;  Bound_v= .false.
        Where (G(v)%wt >= Core_parm(2)) Bound_v(G(v)%ls)= .true.
        
!       Test vertex 'v' for addition to a prior cluster
        
        Do j= 1,ncl
          m= Count(Clustered == j .and. Bound_v)
          
          If (m == ncore(j)) then ! Must be bound to all members of cluster 'j'
            k= k + 1;  ls(k)= j
          End if
        End do
        
!       Add vertex 'v' to the prior cluster 'j' to which it is most closely bound
        
        j= -1
        If (k == 1) then
          j= ls(1)
        Else if (k > 1) then
          Do i= 1,k
            j= ls(i);  sm= 0
            Do l= 1,G(v)%nl
              m= G(v)%ls(l)
              If (Clustered(m) == j) sm= sm + G(v)%wt(l) * G(m)%sum_wt
            End do
            sim(i)= sm / vol(j)
          End do
          
          j= Maxloc(sim,1)
        End if
    
!       Add vertex 'v' to cluster 'j'
          
        If (j > 0) then
          Clustered(v)= j;  ncore(j)= ncore(j) + 1;  vol(j)= vol(j) + G(v)%sum_wt
          Cycle Main_loop
        End if
        
!       Keep 'v' unclustered or start a new cluster
        
        If (ncl >= mc .or. Any(Clustered(G(v)%ls) > 0 .and. G(v)%wt > Core_parm(1))) Cycle Main_loop
          
!       Start a new cluster with vertex 'v'
    
        j= ncl + 1;  ncl= j;  Clustered(v)= j
        ncore(j)= 1;  vol(j)= G(v)%sum_wt

!       Test prior unclustered vertices to see if they can be added to new cluster 'j'
    
        Do i= 1,v-1
          If (Clustered(i)) Cycle
        
          If (ncore(j) == 1) then
            m= 0;  If (Bound_v(i)) m= 1
          Else
            Bound_i= .false.;  Where (G(i)%wt >= Core_parm(2)) Bound_i(G(i)%ls)= .true.
            m= Count(Clustered == j .and. Bound_i)
          End if
          
          If (m == ncore(j)) then  ! Add vertex 'i' to cluster 'j'
            Clustered(i)= j;  ncore(j)= ncore(j) + 1;  vol(j)= vol(j) + G(i)%sum_wt
          End if
        End do
      End do Main_loop
    
!     Compute vertex memberships in each cluster,
!     starting with the core members
        
      Do j= 1,ncl
        n= ncore(j);  Core(j)%n= n;  Core(j)%svl= vol(j)
        Allocate(Core(j)%set(n), Core(j)%mbr(nv))
        m= Count(Clustered == j)
        Call List_of_true (Clustered == j, n,Core(j)%set)
        Core(j)%mbr= 0;  Core(j)%mbr(Core(j)%set)= 1
      End do
        
!     Next determine non-core memberships in each cluster, 
!     and the fuzzy set adjustment, removing clusters     
!     that are too small
      
      Remove_loop : Do n= ncl,2,-1
      
        Call Core_memb (nv,ncl, G,Clustered, Core(:ncl))
    
        j= Minloc(Core(:n)%smb,1);  sm= Core(j)%smb
        If (sm >= Core_parm(3) .or. n <= 2) Exit ! No more removals

        Clustered(Core(j)%set)= -1;  ! Remove cluster 'j'
        
        Do k= j+1,n
          l= k - 1;  m= Core(k)%n;  Core(l)%n= m
          Core(l)%svl= Core(k)%svl;  Core(l)%mbr= Core(k)%mbr

          DeAllocate(Core(l)%set);  Allocate(Core(l)%set(m))
          Core(l)%set= Core(k)%set
        End do
      End do Remove_loop

!     Sort clusters by decreasing weighted membership
      
      ncl= n;  sim(:n)= Core(:n)%smb
      Call Sort (.false., sim(:n), ls(:n),ReOrd)
      
      If (ReOrd) then
        Allocate(CoreP(n))
        Call Copy_Set_list (ls(:n), Core(:n), CoreP(:n))
        ls(:n)= "ID"
        Call Copy_Set_list (ls(:n), CoreP(:n), Core(:n))
        Call DeAlloc_Set_list_ar (CoreP)
      End if

!     Compute "independents", with no core
      
      i= ncl + 1;  Core(i)%n= 0;  Core(i)%svl= 0
      If (Associated(Core(i)%set)) DeAllocate (Core(i)%set)
      Allocate(Core(i)%mbr(nv))

      Do v= 1,nv
        Forall(k=1:ncl) sim(k)= Core(k)%mbr(v)
        Core(i)%mbr(v)= Max(1.0 - Sum(sim(:ncl)), 0.0)
      End do

      Core(i)%smb= Sum(Core(i)%mbr * G%sum_wt)
      sim(:ncl)= Core(:ncl)%smb
      
      If (pr_out > 1) then 
        Call Out ("# core clusters found",ncl)
        Do j= 1,ncl
          Call Out ("For cluster",j, "# in the core set",Core(j)%n, ln=1)
          Call Out ("Core set weight",Core(j)%svl, "total cluster weight", Core(j)%smb)
          If (pr_out > 1.5) then 
            Call Out ("Vertices in the core set", Core(j)%set)
            Call Out ("0 to 1 memberships", Core(j)%mbr)
          End if
        End do

        Call Out ("Independents weight",Core(i)%smb,ln=1)
        If (pr_out > 1.5) Call Out ("0 to 1 independents memberships", Core(i)%mbr)
      End if
      
    Contains   

      Subroutine Core_memb (nv,ncr, G,Clustered, Core)

!       Compute the memberships 'mbr' of non-core vertices in each 
!       'Core' cluster. These membership values are determined by 
!       the graph 'G'. These vertices may be members of more than 
!       one vertex of one cluster, also of more than one cluster.

!       Also duo fuzzy set reduction on the non-core vertices
      
        Integer,           Intent(in) :: nv           ! # vertics
        Integer,           Intent(in) :: ncr          ! # clusters
        Type(Adjacency),   Intent(in) :: G(:)         ! (nv) Weighted graph, with vertices ordered by decreasing weight
                                                      !  %nl      = # adjacent vertices
                                                      !  %sum_wt  = vertex weight
                                                      !  %ls(nl)  = List of adjacent vertices
                                                      !  %wt(nl)  = Their decreasing edge weights (0 to 1)
                                                      !             correlations
        Integer,           Intent(in) :: Clustered(:) ! (nv) Clustered(v) = j if v is a core member of cluster j
                                                      !      else Clustered(v) = -1 
        Type(Set_list), Intent(inout) :: Core(:)      ! (ncl) The core sets of the clusters
                                                      !  %n       = # vertices 'n' in the core set
                                                      !  %set(n)  = list of those vertices
                                                      !  %mbr(nv) = 0..1 memberships of all vertices in the core set
                                                      !  %svl     = sum of the core vertex weights
                                                      !  %smb     = total weighted membership of the cluster
!     Local:
        Real    :: mbr(nv)  ! Non-core memberships for each cluster
        Real    :: smb(nv)  ! Non-core memberships summed over all clusters
        Integer :: mb(ncr)
        Integer :: i, j, k, n, v, vcr

!       Non-core vertex memberships in the clusters
        
        Cluster_loop : Do j= 1,ncr
          mbr= 0
          Do k= 1,Core(j)%n
            vcr= Core(j)%set(k);  mbr(vcr)= 1  ! vcr = kth core member of cluster 'j'
            
            Do i= 1,G(vcr)%nl                  ! ith member of 'vcr'
              v= G(vcr)%ls(i)  
              If (Clustered(v) < 1) mbr(v)= mbr(v) + G(vcr)%wt(i)  ! 'v' = non-core member of 'j'
            End do    
          End do
          
          Core(j)%mbr= mbr
        End do Cluster_loop

!       Sum non-core memberships for fuzzy set reduction

        smb= 1
        Vertex_loop : Do v= 1,nv
          If (Clustered(v) < 1) then ! Non-core member
            Forall(j=1:ncr) mb(j)= Core(j)%mbr(v)
            smb(v)= Sum(mb(:ncr))
          End if
        End do Vertex_loop
      
!       Do fuzzy set reduction

        Cluster_loop2 : Do j= 1,ncr
          Where (smb > 1) Core(j)%mbr= Core(j)%mbr / smb  ! Fuzzy reduction
          Core(j)%smb= Sum(Core(j)%mbr * G%sum_wt)  ! Total weight of cluster 'j'
        End do Cluster_loop2

      End subroutine Core_memb
      
    End subroutine Core_clustering
  

    Subroutine Slate_graph (ns,cut, Memb, G, sum_deg)
    
!     Form the graph 'G' of slate cluster correlations above a cutoff value 'cut'
    
      Integer,             Intent(in) :: ns        ! # vertices = # slate ballots
      Real,                Intent(in) :: cut       ! Correlation cut value, > 0, to determine edge connectivity
      Type(Multi_listD),   Intent(in) :: Memb(0:)  ! (0:ns) Slate cluster data
      
      Type(Adjacency),  Intent(inout) :: G(:)      ! (ns) Weighted graph
                                                   !  %sum_wt  = vertex weight
                                                   !  %nl      = # adjacent vertices = degree
                                                   !  %ls(nl)  = List of adjacent vertices, ordered by decreasing %wt
                                                   !  %wt(nl)  = Their 0.. 1 edge weights = slate cluster modified 
                                                   !             correlations ordered by decreasing correlation
      Integer,            Intent(out) :: sum_deg   ! Sum of the degrees of the graph = Sum(G%nl)
!   Local:
      Real(Dblp) :: cor, wt(ns)
      Integer :: i, j, n, s1, sl, lst(ns), key(ns)
      
      If (cut <= 0) then
        Call Out ("Error: correlation cutoff must be positive"); Stop
      End if
      
      wt= 0 ; lst= 0;  key= 0
      Do sl= 1,ns
        If (Associated(G(sl)%ls)) then  ! Assume a prior graph with a lower value of 'cut'
          n= 0
          Do j= 1,G(sl)%nl
            If (G(sl)%wt(j) > cut) then
              n= n + 1;  lst(n)= G(sl)%ls(j);  wt(n)= G(sl)%wt(j)
            End if
          End do
          DeAllocate(G(sl)%ls,G(sl)%wt)
        Else
          n= 0
          Do s1= 1,ns
            If (sl == s1) Cycle
            Call Dot_product_M (Dot_fac, Memb(sl)%tx(1:),Memb(s1)%tx(1:), i,cor)
            If (cor > cut) then
              n= n + 1;  lst(n)= s1;  wt(n)= cor
            End if
          End do
        End if
        
        G(sl)%nl= n
        
        If (n > 0) then
          Call Sort (.false., wt(:n), key(:n))
          Allocate(G(sl)%ls(n), G(sl)%wt(n))
          G(sl)%ls= lst(key(:n));  G(sl)%wt= wt(:n)
        End if
      End do
      
      sum_deg= Sum(G%nl)  
      
      If (pr_out > 1.5) then
        Call Out ("Sum of the degrees in'Slate_graph'", sum_deg, ln=1)
        Call Out ("List of the degrees", G%nl)
      End if   
    End Subroutine Slate_graph
  
    Subroutine Degree_graph (Clust_opt, G, Coef, max_deg)
  
!     Compute degree clustering coefficients 'Coef' for the vertex and edge 
!     weighted graph 'G', where 'Clust_opt' specifies the method of computation. 
    
!     Clust_opt = 1: The degree clustering value of a vertex is defined as the 'product' 
!                    of the weights of the two end vertices times the edge weight between them, 
!                    over all edges containing the given vertex with edge weight of sufficient size.
    
!     Clust_opt = 2: The degree clustering value of a vertex is defined as the 'average' 
!                    of the weights of the two end vertices times the edge weight between them, 
!                    over all edges containing the given vertex with edge weight of sufficient size.
    
      Integer,         Intent(in) :: Clust_opt ! Option for how to compute the clustering values

      Type(Adjacency), Intent(in) :: G(:)      ! (nv) Vertex and edge weighted graph
                                               !  %sum_wt  = vertex weight
                                               !  %ls(nl)  = List of adjacent vertices
                                               !  %wt(nl)  = Their edge weights (Init_cut to 1), ordered 
                                               !             by decreasing weight
      
      Real,           Intent(out) :: Coef(:)   ! (nv) Degree clustering coefficients, scaled from 0 to 100
      Real,           Intent(out) :: max_deg   ! Maximum of the degree clustering coefficients before normalization
!   Local:
      Real    :: sm
      Integer :: nv, v

      nv= Size(G);  Coef= 0
    
      Select Case (Clust_opt)
        Case (1)
          Do v= 1,nv
            If (G(v)%nl > 0) then
              sm= Sum(G(v)%wt * G(G(v)%ls)%sum_wt)
              Coef(v)= G(v)%sum_wt * sm
            End if
          End do
          
        Case (2)
          Do v= 1,nv
            If (G(v)%nl > 0) then
              sm= Sum((G(v)%sum_wt + G(v)%wt) * G(G(v)%ls)%sum_wt)
              Coef(v)= G(v)%sum_wt * sm / 2
            End if
          End do
      End Select
      
      max_deg= Maxval(Coef);  If (max_deg > 0) Coef= Coef * (100.0 / max_deg)
      
      If (pr_out > 1.5) then
        Call Out ("For degree clustering option",Clust_opt , &
                  "Max coefficient in 'Degree_graph'", max_deg, ln=1)
        If (max_deg > 0) then
          Call Out ("The degree clustering coefficients, normalized to 100", Coef)
        End if   
      End if   

    End subroutine Degree_graph
    
    Subroutine Triangle_graph (Clust_opt, G, Coef, max_deg)
  
!     Compute triangular clustering coefficients 'Coef' for the vertex and edge 
!     weighted graph 'G', where 'Clust_opt' specifies the method of computation.  
    
!     Clust_opt = 1: The triagular clustering coefficient of a vertex is defined as the 
!                    'product' of the weights of the 3 vertices of a triangle times the 
!                    'geometric mean' of the 3 edge weights of a triangle, summed over all triangles 
!                    containing the given vertex.
    
!     Clust_opt = 2: The triagular clustering coefficient of a vertex is defined as the 
!                    'arithmetic average' of the weights of the 3 vertices of a triangle times
!                    the 'geometric mean' of the 3 edge weights of a triangle, summed over all  
!                    triangles containing the given vertex.
    
!     Then locate the peaks of these coefficients, among closely related vertices.
    
      Integer,            Intent(in) :: Clust_opt ! Option for how to compute the clustering values
      Type(Adjacency), Intent(inout) :: G(:)      ! (nv) Vertex and edge weighted graph
                                                  !  %nl      = # adjacent vertices
                                                  !  %sum_wt  = vertex weight
                                                  !  %ls(nl)  = List of adjacent vertices
                                                  !  %wt(nl)  = Their edge weights (0 to 1), 
                                                  !             ordered by decreasing weight
                                                  !  %wl(nl)  = Cube root of %wt
      
      Real,              Intent(out) :: Coef(:)   ! (nv) Triangular clustering coefficients, scaled from 0 to 100
      Real,              Intent(out) :: max_deg   ! Maximum of the degree clustering coefficients, before normalization
!   Local:
      Real     :: sm, tmp
      Integer  :: i, j, k, n, nv, v, v1, v2

      nv= Size(G);  Coef= 0
      Do v= 1,nv
        n= G(v)%nl
        If (n > 0) then
          If (Associated(G(v)%wl)) DeAllocate(G(v)%wl)
          Allocate(G(v)%wl(n))
          G(v)%wl= G(v)%wt**(1.0/3.0)
        End if
      End do
      
      Select Case (Clust_opt)
        Case (1)
          Do v= 1,nv
            n= G(v)%nl;  If (n < 2) Cycle
            
            Do i= 1,n-1
              v1= G(v)%ls(i);  If (G(v1)%nl < 1) Cycle
              sm= 0
          
              Do j= i+1,n
                v2= G(v)%ls(j)
                k= First_true(G(v1)%ls == v2);  If (k < 1) Cycle ! No triangle
          
                sm= sm + (G(v)%wl(j) * G(v1)%wl(k)) * G(v2)%sum_wt
              End do
              
              Coef(v)= Coef(v) + (G(v)%wl(i) * G(v1)%sum_wt) * sm
            End do
      
            Coef(v)= G(v)%sum_wt * Coef(v)
          End do
          
        Case (2)
          Do v= 1,nv
            n= G(v)%nl;  If (n < 2) Cycle
            
            Do i= 1,n-1
              v1= G(v)%ls(i);  If (G(v1)%nl < 1) Cycle
              
              sm= 0;  tmp= G(v)%sum_wt + G(v1)%sum_wt
          
              Do j= i+1,G(v)%nl
                v2= G(v)%ls(j)
                k= First_true(G(v1)%ls == v2);  If (k < 1) Cycle
          
                sm= sm + (G(v)%wl(j) * G(v1)%wl(k)) * (tmp + G(v2)%sum_wt)
              End do
              Coef(v)= Coef(v) + G(v)%wl(i) * sm / 3
            End do
          End do
      End Select
      
      max_deg= Maxval(Coef);  If (max_deg > 0) Coef= Coef * (100 / max_deg)
      
      If (pr_out >= 1.5) then
        Call Out ("For triangle clustering option",Clust_opt ,&
                  "Max coefficient in 'Triangle_graph'", max_deg, ln=1)
        If (max_deg > 0) then
          Call Out ("The triangle clustering coefficients, normalized to 100", Coef)
        End if   
      End if   
    End subroutine Triangle_graph
    
    Subroutine Peak_clustering (nv,G, coef_cut,Coef, npk,peak,pk_val)
  
!     Compute the peaks of the vertex clustering values of a graph, except
!     for very small peaks, as specified by 'coef_cut'
    
      Integer,         Intent(in) :: nv        ! # vertices
      Type(Adjacency), Intent(in) :: G(:)      ! (nv) Weighted graph
                                               !  %nl      = # adjacent vertices
                                               !  %sum_wt  = vertex weight
                                               !  %ls(nl)  = List of adjacent vertices
                                               !  %wt(nl)  = Their edge weights or correlations (0 to 1), ordered by decreasing weight
      
      Real,            Intent(in) :: coef_cut  ! Minimum coefficient value of a peak, as a fraction of Maxval(Coef)
      Real,            Intent(in) :: Coef(:)   ! (nv) Vertex clustering values

      Integer,        Intent(out) :: npk       ! # peak vertices found
      Integer,        Intent(out) :: peak(:)   ! (nv) List of peak vertices found
      Real,           Intent(out) :: pk_val(:) ! (nv) List of peak values relative to Maxval(Coef)
!   Local:
      Logical :: ReOrd
      Real    :: mx, min_pk, mx_coef
      Integer :: v, key(nv)

      npk= 0;  peak= 0;  pk_val= 0
      mx_coef= Maxval(Coef);  min_pk= coef_cut * mx_coef
      
      Do v= 1,nv
        If (.not.Associated(G(v)%ls) .or. Coef(v) <= min_pk) Cycle
        
        mx= Maxval(Coef(G(v)%ls))
        If (Coef(v) > mx) then   ! If mx = 0, then it's an isolated peak
          npk= npk + 1;  peak(npk)= v;  If (npk >= Size(peak)) Exit
        End if
      End do
      If (npk < 1) Return
    
      pk_val(:npk)= Coef(peak(:npk)) / mx_coef
      
      If (npk > 1) then
        Call Sort (.false., pk_val(:npk), key, ReOrd)
        If (ReOrd) peak(:npk)= peak(key)
      End if
      
      If (pr_out > 1) Then 
        Call Out ("Peaks of the vertex clustering values", peak(:npk))
        Call Out ("Their peak values", pk_val(:npk))
        Call Out ("# vertices in each peak region", G(peak(:npk))%nl + 1)
      End if

    End subroutine Peak_clustering

  
    Subroutine Iterate (idist,iq,nc,n0, Mxitr,Parm_cnv,Parm_obj, Memb,Clust, ng)

!     Iterate to convergence an initial cluster set of size greater than 1 
    
      Integer,       Intent(in) :: idist     ! Index of the voting district
      Integer,       Intent(in) :: iq        ! Initial cluster set
      Integer,       Intent(in) :: nc        ! # candidates
      Integer,       Intent(in) :: n0        ! Initial # clusters in the cluster set
      Integer,       Intent(in) :: Mxitr(:)  ! (3) Max # iterations at each restart level

      Real, Intent(in) :: Parm_cnv(:,:) ! Convergence parameters
                                   ! For convergence phase RS = 1, 2, 3:
                                   ! (RS,1) = convergence tolerances
                                   ! (RS,2) = min merge criterion for 
                                   !   cluster correlation
                                   ! (RS,3) = min size criterion for cluster
                                   !   deletion
                                   ! (RS,4) = max correlation criterion for
                                   !   cluster deletion
                                   ! (1,5) = maximum norm permitted for an
                                   !   update differential
                                   ! (2,5) = match tolerance for measures of
                                   !   cluster set similarity
                                   ! (3,5) = min separation criterion for the
                                   !   separation graph for the clique of
                                   !   correlated clusters 
      Real, Intent(in) :: Parm_obj(:) ! (12) Parameters for the penalty factors 
                                   !        of the objective function
                                   !  1:2 = cluster density soft cutoff limits (cos_rise)
                                   !  3:4 = cluster boundary soft cutoff limits (cos_fall)
                                   !  5:6 = cluster size soft cutoff limits (cos_rise)
                                   !  7:8 = cluster correlation soft cutoff limits (cos_fall)
                                   !  9:10 = independents size cutoff limits (cos_fall)
                                   !  11:12 = convergence residual limits (cos_fall)
    
      Type(Multi_listD),    Intent(in) :: Memb(0:)   ! (0:ns) Slate cluster data
      
      Type(Multi_listD), Intent(inout) :: Clust(0:)  ! (0:ind) Cluster set data
                                                     ! For 0:
                                                     !   %k      = Error code. Invalid results if < 0.
                                                     !   %m      = # candidates = nc
                                                     !   %n      = # regular clusters = n0
                                                     !   %l      = # clusters, including independent = ind

                                                     !   %px(4):  Objective data: min size, max correlation, 
                                                     !              independents size (frac), final residual
                                                     !   %qx(4):  Objective factors: min size penalty, max correlation penalty, 
                                                     !              independents size penalty, residual penalty
                                                     !   %rx(3):  Clustering objective value(1)= coherent vote(2) * penalty(3)
                                                     !   %sx(ind): Cluster sizes

                                                     !   %lt(n0)  Mapping from initial to converged clusters
                                                     !   %wt(ns): Sum of slate cluster memberships over all reg clusters,
                                                     !            reduced by the fuzzy set factor to be <= 1.0 

                                                     !   %vl(0:6): Convergence data
                                                     !     0 = original # regular clusters = n0
                                                     !     1 = final # regular clusters
                                                     !     2 = final success code
                                                     !     3 = total # cluster merge and deletion operations
                                                     !     4 = # convergence calls
                                                     !     5 = total # centroiding iterations
                                                     !     6 = total # line search points

                                                     !   %L0(0:4,0:iv) Integer convergence data per update call
                                                     !     0  = restart phase
                                                     !     1  = success code:  2 = fully converged,   1 = reducing functional,
                                                     !                         0 = damped reduction, -1 = convergence failure
                                                     !     2 = # reg clusters
                                                     !     3 = # centroiding iterations
                                                     !     4 = # line search points

                                                     !   %M2(0:2,0:iv) Real convergence data
                                                     !     0 = final line search damping parameter
                                                     !     1 = Final min. functional
                                                     !     2 = net (output - input) difference

                                                     !   %M1(ng,ng) Cluster correlation matrix, 0 on diagonal
      
                                                     !   %M3(6,ind) Cluster data
                                                     !     1 = size 
                                                     !     2 = density penalty
                                                     !     3 = boundary penalty
                                                     !     4 = density
                                                     !     5 = boundary
                                                     !     6 = width

                                                     !  %L1(-2:nc,0:ind) Candidate ordering data for
                                                     !                    noise zeroed vectors
                                                     !     (-2:0,:) = # top, significant, viable cand's
                                                     !     (1:,1:)  = Decreasing ordering from %T0(:,1:,2)
                                                     !     (1:,0)   = Decreasing ordering from %T0(:,0,1)
                                                     !               
                                                     !  %T0(nc,0:ind,3) Cluster vector data
                                                     !      (:,1:,1) = Cluster mean vectors
                                                     !      (:,1:,2) = Cluster mean vectors with noise zeroing 
                                                     !      (:,1:,3) = Cluster portion vectors
                                                     !      (:,0,1)  = Zrate0 * declining regular portion factor
                                                     !      (:,0,2)  = Zrate0 = cluster averaged noise zeroed mean vectors

                                                     ! For clusters 1:ind
                                                     !   %k  = # fuzzy reduced slate memberships
                                                     !   %n  = # slate ballots which are partial or full members
                                                     !   %o   : 1 for regular clusters, 2 for independents
                                                     !   %fsx : slate ballot weight averaged fuzzy reduction factor
                                                     !   %fux : cluster width
                                                     !   %sum_wt = sum of slate memberships

                                                     !   %ls(n)  = list of slate members
                                                     !   %wt(n)  = weighted slate memberships
                                                     !   %rx(n)  = unweighted slate memberships
                                                     !   %px(nc) = variance of %sx
                                                     !   %sx(nc) = mean rating vector
                                                     !   %tx(0:nc)= centered and normalized rating vector
                                                     !   %ux(nc)  = total cluster weight, incorporating slate ballot variance 
                                                     !              normalization (Memb%ux * %rx)
      Integer,  Intent(out) :: ng  ! # clusters, initial to final
!   Local:
      Real,    Allocatable :: Corr(:,:)   ! (ng,ng) Correlation matrix
      Real,    Allocatable :: cnvR(:,:)   ! (0:2,0:md) Final 'Converge_clusters' values:
                                          !   0 = damping parm. 1 = min. functional. 2 = net (out - in) difference
      Integer, Allocatable :: cnvI(:,:)   ! (0:4,0:md) Final 'Converge_clusters' data:
                                          !   0  = restart phase
                                          !   1  = success code:  2 = fully converged,   1 = reducing functional,
                                          !                       0 = damped reduction, -1 = convergence failure
                                          !   2 = # reg clusters
                                          !   3 = # centroiding iterations
                                          !   4 = # line search points

      Integer :: Restart  ! 1 : Level 1 restart (strongest merge and delete requirements) - initially or after a merge
                          ! 2 : Level 2 restart (merge or delete more easily) - after no merges but possibly deletions
                          ! 3 : No restart (after no reduction in # clusters at level 2 unless to n0 = 1). Do final convergence

      Integer, Allocatable :: map(:,:,:)  ! (0:n0,2,md)  Mapping data between initial and intermediate clusters 
                          ! (1:n0,1,ic) Maps initial clusters onto intermediate clusters, reordered, 
                          !             for clusters change counter 'ic'
                          ! (1:ng,2,ic) For each intermediate cluster, the largest initial cluster mapped to it
                          ! (0,1,ic) = 'ng'.  (0,2,ic) = 'icyc'

      Type(Set_list), Pointer :: Cliq_set(:)=>Null()
      Real    :: siz(n0+1), vec(0:nc,n0+1), avg(0:nc,n0), prc(n0), mx_cor(n0)
      Integer :: ls(n0), key(n0), inv(n0), clq_map(n0,2), sing(n0), pair(2,n0)

      Real    :: mxc, obj, min_size
      Logical :: Err, ReOrd, Keep_mask(n0)
      Integer :: md, ind, ier, itc, itr, ndl, nmg, icyc, ncyc, Success
      Integer :: i, j, k, l, n, q, r, cl, ic, iv, mq, mx, n1, n2, nl, np, nq, ns, r1, RS

      np= Memb(0)%l;  ns= Ubound(Memb,1)
      ng= n0;  ind= ng + 1;  Clust(0)%k= 1
      
      Call Out ("Enter 'Iterate' for district",idist, "and initial cluster set",iq, ln=1)
      Call Out ("Start with # regular clusters",ng)

      k= Ubound(Clust,1)
      If (ind > k) then
        Call Out ("Error in 'Iterate': 'Clust' too small",k, "but needs",ind, ln=1)
        Clust(0)%k= -1;  Return
      End if

      vec(0,:)= Clust(1:)%sum_wt;  Forall(cl=1:ind) vec(1:,cl)= Clust(cl)%sx

      If (pr_out > 1) then
        Call Out ("with cluster sizes",vec(0,:))
        If (pr_out > 1.5) Call Out (-1,"Mean rating vectors",vec(1:,:))
      End if

      ic= 0;  ndl= 0;  nmg= 0;  Restart= 1;  Clust(0)%vl= 0
      md= 3*(ng-1) + 1;  mx= Mxitr(Restart)
      
      Allocate(cnvI(0:4,0:md), cnvR(0:2,0:md), map(0:n0,2,md))

      cnvI= 0;  cnvR= 0;  map= 0;  vec= -1;  avg= -1
      siz= -1;  ls= -1;  clq_map= -1;  sing= -1;  pair= -1
      icyc= 0;  ncyc= icyc;  iv= 0

      Call Converge_clusters (iq,0, nc,n0, mx,Restart, Parm_cnv(:,1),Memb, &
                              ng,Clust(:ind), cnvI(:,iv),cnvR(:,iv))
      ier= Clust(0)%k

      If (ier < 0) then
        Call Out ("Error in 'Iterate': final # regular clusters",ng, &
                  "error code",ier, ln=1)
        Call Out ("At initial convergence with initial # reg. clusters",n0)
        Clust(0)%k= -1;  Return  ! Skip this cluster set

      Else if (pr_out > 1) then
        If (ng < n0) then 
          Call Out ("First convergence in 'Iterate': # reg. clusters reduced from",n0, &
                     "to",ng, ln=1)
        Else
          Call Out ("First convergence in 'Iterate': Success for # reg. clusters",ng, ln=1)
        End if
      End if

      ind= ng + 1
      If (ng < n0) then
        ic= ic + 1;  map(0,1,ic)= ng;  map(0,2,ic)= 0
        map(1:ng,2,ic)= Clust(0)%lt(:ng);  Forall(cl=1:ng) map(Clust(0)%lt(cl),1,ic)= cl
      End if
      
      Merge_cycle : Do icyc= 1,md
          
        If (pr_out > 1) Call Out ("For merge/delete test cycle",icyc, ln=1)
       
!       Compute the cluster correlations and identify weighted cliques of clusters to average

        If (ng > 1) then
          If (Allocated(Corr) .and. ng /= Size(Corr,1)) DeAllocate (Corr)
          If (.not.Allocated(Corr)) Allocate(Corr(ng,ng));  mx_cor= -1
        
          vec(0,:ng)= Clust(1:)%sum_wt;  Forall(cl=1:ng) vec(1:,cl)= Clust(cl)%sx
          vec(:,ng+1:)= 0

          Call Correlate_clusters1 (Clust(1:ng), mx_cor(:ng),Corr)
          mxc= Maxval(mx_cor(:ng));  siz(:ng)= Real(Clust(1:ng)%sum_wt)

          Call Consolidate_clusters (ng, Parm_cnv(Restart,2),mxc, siz(:ng), &
                                     Corr, n1,sing,n2,pair, nq,Cliq_set)
          mq= n2 + nq;  nl= n1 + mq
          
          If (nl < ng) then
            nmg= nmg + 1
            If (pr_out > 1.5) then
              Call Out ("# singleton cluster cliques",n1, "# pair cliques",n2)
              Call Out ("# general cliques",nq)
              Call Out ("Prior # clusters",ng, "to be reduced to",nl)
            End if  

            If (nl < 1) then
              Call Out ("Error in 'Iterate': Consolidation failed for # regular clusters",&
                         nl, "at cycle",icyc, ln=1)
              ng= nl;  Clust(0)%k= -1;  Return  ! Skip this cluster set
            End if
          End if
        Else if (ng == 1) then
          n1= 1;  n2= 0;  nq= 0;  mq= 0;  nl= 1;  mxc= 0
        Else  ! Error condition (ng <= 0)
          Call Out ("Serious Error in 'Iterate' at cycle",icyc, ln=1)
          Stop   ! This should be impossible
        End if
          
        If (nl < ng) then   ! Consolidate the cliques into 'nl' clusters
          Forall(i=1:n2) prc(i)= Corr(pair(1,i),pair(2,i))
          
          Call Average_consolidated (n1,n2,nq, sing(:n1),pair(:,:n2),prc(:n2), &
                                     Cliq_set, vec(:,:ng), clq_map(:ng,:),avg(:,:nl))
          
          siz(:nl)= avg(0,:nl) 
          Call Sort (.false., siz(:nl), key(:nl), ReOrd= ReOrd)
          Clust(0)%n= nl;  Clust(1:nl)%sum_wt= siz(:nl)
          Forall(cl=1:nl) Clust(cl)%sx= avg(1:,key(cl))

          If (ReOrd) then
            Call Inverse_map (key(:nl),inv)
            Forall(cl=1:nl) inv(key(cl))= cl
            clq_map(:ng,1)= inv(clq_map(:ng,1))
            clq_map(:nl,2)= clq_map(key(:nl),2)
          End if 

          ic= ic + 1;  map(0,1,ic)= nl;  map(0,2,ic)= icyc
          map(1:ng,1,ic)= clq_map(:ng,1)
          map(1:nl,2,ic)= clq_map(:nl,2)

          If (nq > 0) then
            Call DeAlloc_set_list_ar (Cliq_set);  DeAllocate(Cliq_set)
          End if
          
          If (pr_out > 1) then
            Call Out ("Consolidation in 'Iterate': from # regular clusters",ng, &
                      "to",nl, ln=1)
            If (pr_out > 1.5) then
              vec(0,:nl)= siz(:nl);  Forall(cl=1:nl) vec(1:,cl)= Clust(cl)%sx
              Call Out ("For cluster sizes",vec(0,:nl))
              Call Out (-1,"Cluster mean vectors",vec(1:,:nl))
            End if
          End if
          
        Else if (ng >= 2) then !  No merging of clusters, but check for 
                               !  small, isolated clusters to delete
          siz(:ng)= Clust(1:ng)%sum_wt
          
          Call List_of_true (siz(:ng)    >= Parm_cnv(Restart,3) .or. &
                             mx_cor(:ng) >= Parm_cnv(Restart,4), nl,ls) 

          If (nl < 1) then
            nl= 1;  ls(1)= 1
            If (pr_out > 1) Call Out ("Keep the largest cluster, even if small and well separated")
          End if
          
          If (nl < ng) then ! Keep clusters which are not small and isolated and continue convergence
            ic= ic + 1;  map(0,1,ic)= nl;  map(0,2,ic)= icyc
            map(1:nl,2,ic)= ls(:nl)
            Forall(cl=1:nl) map(ls(cl),1,ic)= cl

            ndl= ndl + (ng-nl);  Clust(0)%n= nl
            Clust(1:nl)%sum_wt= Clust(ls(:nl))%sum_wt;  
            siz(:nl)= Clust(1:nl)%sum_wt
            Forall(cl=1:nl) Clust(cl)%sx= Clust(ls(cl))%sx
            
            If (pr_out > 1) then
              Call Out ("Keep the clusters",ls(:nl))
              If (pr_out > 1.5) then
                Call Out ("Based on cluster sizes",siz(:ng))
                Call Out ("and max correlations",mx_cor(:ng))
                Call Out ("for restart level",Restart)
              End if
            End if
          Else if (nl < 1) then
            Call Out ("Warning in 'Iterate': all clusters deleted at cycle", &
                       icyc, "for restart level",Restart, ln=1)
            Call Out ("Based on cluster sizes",siz(:ng))
            Call Out ("and maximum correlatins",mx_cor(:ng))
            ng= nl;  Clust(0)%k= -1;  Return  ! Skip this cluster set
          End if
        End if
          
        If (nl == ng) then  ! No change in # clusters, with ng >= 2
          If (Restart <= 1) then
              Restart= 2
          Else if (Restart == 2) then
            Restart= 3
          Else
            Exit Merge_cycle ! Convergence done
          End if
        Else   ! nl < ng due to merging or deletion
          
          ng= nl
          
          If (ng < 1) then  ! Error condition
            Call Out ("Serious Error in 'Iterate' at cycle",icyc, ln=1)
            Stop   ! This should be impossible
          Else if (ng == 1) then  ! Refine only
            Restart= 3
          Else if (mq > 0) then   ! Merging case
            Restart= 1
          Else                    ! Deletion case
            Restart= 2
          End if
        End if
        
        If (pr_out > 1) then
          Call Out ("Completed convergence cycle with # clusters",ng, ln=1)
          Call Out ("Next restart status",Restart)
        End if
        
        If (Restart > 1 .and. ng > mxcl) then
          ng= mxcl;  Restart= 1
        End if  
      
        nl= ng;  ind= nl + 1;  mx= Mxitr(Restart);  Clust(0)%n= nl;  
        vec(0,:nl)= Clust(1:n0)%sum_wt;  Forall(cl=1:nl) vec(1:,cl)= Clust(cl)%sx
        vec(:,ind:)= 0;  iv= iv + 1

        Call Converge_clusters (iq,icyc, nc,nl, mx,Restart, Parm_cnv(:,1),Memb, &
                                ng,Clust(:ind), cnvI(:,iv),cnvR(:,iv))

        ier= Clust(0)%k;  ng= Clust(0)%n;  ind= ng + 1;  ncyc= icyc;  Success= cnvI(1,iv)
        vec(0,:ng)= Clust(1:ng)%sum_wt;  Forall(cl=1:ng) vec(1:,cl)= Clust(cl)%sx

        If (Success < 0 .and. Restart == 3) ier= Min(ier,0)

        If (ier < 0 ) then
          Call Out ("Convergence failure in 'Iterate': for initial set",iq, &
                    "at cycle",icyc, ln=1)
          Call Out ("Restart phase, Success code, # clusters, # centroidings",cnvI(:3,iv))
          Return
        Else if (ier == 0 ) then
          Call Out ("Convergence problem in 'Iterate': for initial set",iq, &
                    "at cycle",icyc, ln=1)
          Call Out ("Restart phase, Success code, # clusters, # centroidings",cnvI(:3,iv))
          Exit Merge_cycle
        Else if (ng < nl) then 
          ic= ic + 1;  map(0,1,ic)= ng;  map(0,2,ic)= icyc
          map(1:ng,2,ic)= Clust(0)%lt(:ng);  
          Forall(cl=1:ng) map(Clust(0)%lt(cl),1,ic)= cl

          If (pr_out > 1) then
            Call Out ("Successful Cycle",icyc, "completed in 'Iterate': # clusters reduced to",ng, ln=1)
            Call Out ("Restart phase, Success code, # clusters, # centroidings",cnvI(:3,iv))
          End if
          If (ng < 2) Exit Merge_cycle

          Restart= Min(Restart,2)  ! Backup restart phase after a merge or delete
        Else
          If (pr_out > 1) then
            Call Out ("Successful Cycle",icyc, "completed  in 'Iterate' for # clusters", ng, ln=1)
            Call Out ("Restart phase, Success code, # clusters, # centroidings",cnvI(:3,iv))
          End if
          If (ng < 2) Exit Merge_cycle
        End if
        
      End do Merge_cycle
      
      ind= ng + 1
      Clust(0)%vl(0)= n0
      Clust(0)%vl(1)= ng  
      Clust(0)%vl(2)= cnvI(1,iv)
      Clust(0)%vl(3)= nmg + ndl
      Clust(0)%vl(4)= iv
      Clust(0)%vl(5)= Sum(cnvI(3,:iv))  
      Clust(0)%vl(6)= Sum(cnvI(4,:iv))

      If (Associated(Clust(0)%L0)) DeAllocate(Clust(0)%L0, Clust(0)%M2)
      Allocate(Clust(0)%L0(0:4,0:iv), Clust(0)%M2(0:2,0:iv))
      Clust(0)%L0= cnvI(:,:iv);  Clust(0)%M2= cnvR(:,:iv)

      siz(:ind)= Clust(0)%sx(:ind)

      If (Associated(Clust(0)%sx)) DeAllocate(Clust(0)%sx)
      Allocate(Clust(0)%sx(ind));  Clust(0)%sx= Clust(1:ind)%sum_wt

      If (Any(Real(siz(:ind)) /= Real(Clust(0)%sx))) then
        Call Out ("Warning in 'Iterate': cluster size mismatches")
      End if

      If (Associated(Clust(0)%M1)) DeAllocate(Clust(0)%M1)

      If (ng > 1) then
        Call Correlate_clusters1 (Clust(1:ng), mx_cor(:ng),Corr(:ng,:ng))
        Allocate(Clust(0)%M1(ng,ng));  Clust(0)%M1= Corr(:ng,:ng)
        If (pr_out > 1) Call Out (-1,"Final correlation matrix",Corr(:ng,:ng))
      Else
        mx_cor= 0
      End if

!     Mapping of initial to final clusters

      ls= "ID"
      Do i= 1,ic
        Where(ls > 0) ls= map(ls,1,i)
      End do

      If (Size(Clust(0)%lt) /= n0) then
        DeAllocate(Clust(0)%lt);  Allocate(Clust(0)%lt(n0)) 
      End if
      Clust(0)%lt= ls
      
!     Compute objective data for the converged cluster set
      
      If (pr_out > 1) then
        Call Out ("Convergence data for 'Iterate' from initial cluster set",iq, ln=1)
        If (ic > 0) then
          Call Out ("Mapping of initial clusters to final clusters", Clust(0)%lt)
          If (ic > 1) then
            Do i= 1,ic
              Call Out ("For cluster change",i, "at cycle", map(0,2,i));  nl= map(0,1,i)
              Call Out ("Mapping of initial clusters to changed clusters", map(1:n0,1,i))
              Call Out ("The changed clusters mapped to their largest original clusters", map(1:nl,2,i))
            End do
          End if
        End if
        Call Out ("With final cluster sizes", Clust(0)%sx)

        Call Out ("# convergence calls, # merger & deletion operations",&
                  Clust(0)%vl(2:4))
        Call Out ("# centroiding iterations",Clust(0)%vl(5), &
                 "# line search points", Clust(0)%vl(6))
        
        Call Out (-1, "By cycle: Restart, success, # clusters, # iterations, # fcn calls",cnvI(:,:ncyc))
        Call Out (-1, "By cycle: Final line search parm 't', residual, net change",cnvR(:,:ncyc))
      End if
      
!     Compute the final cluster set objective data

      nl= Size(Clust(0)%M3,2)
      If (nl /= ind) then
        DeAllocate(Clust(0)%M3);  Allocate(Clust(0)%M3(6,ind))
      End if

      Clust(0)%M3(1,:)= Clust(0)%sx  
      Clust(0)%M3(6,:)= Clust(1:ind)%fux

      Call Density_n_boundary (ns,Memb(1:),Clust(1:ind), Clust(0)%M3(4:5,:ind))

      If (ng > 1) then 
        Clust(0)%px(1)= Minval(siz(:ng));  Clust(0)%px(2)= Maxval(mx_cor(1:ng)) 
      Else
        Clust(0)%px(1)= siz(1);  Clust(0)%px(2)= -1
      End if
      Clust(0)%px(3)= siz(ind) / np;  Clust(0)%px(4)= cnvR(2,ncyc)   

      Call Objective_values (np,Parm_obj, Clust(0)%sx(:ng),Clust(0)%px,    &
                             Clust(0)%M3(4:5,:ind), Clust(0)%M3(2:3,:ind), &
                             Clust(0)%qx, Clust(0)%rx)
      If (pr_out > 1) then
        Call Out ("From 'Iterate': objective, coherent vote, penalty", Clust(0)%rx)  
        Call Out ("Pen data: min_siz, max_cor, ind_siz, last_res", Clust(0)%px)
        Call Out ("Corresponding penalties", Clust(0)%qx)
        Call Out (-1,"Cluster data: size, density & boundary + their penalties, width", Clust(0)%M3)
      End if

      Clust(0)%M3(2:3,ind)= 1
      
      If (Any(Clust(0)%rx <= obj_eps2)) then
        Clust(0)%k= 0
        Call Out ("Warning in Iterate. Bad cluster set objective", Clust(0)%rx)
        Call Out ("For initial # clusters",n0, "vs final",ng)
        Call Out ("For district",idist, "converged from initial cluster set",iq)
      End if
      
    End Subroutine Iterate

  Subroutine Record_clusters (nd,ncl,Memb, Clust, Clust_set)
  
!   Record the new cluster set unless the 'Clust_set' of cluster sets is already full
!   of those of the same size and the new's objective doesn't improve on any of the prior ones.

    Integer,              Intent(in) :: nd          ! new cluster set index
    Integer,              Intent(in) :: ncl         ! # regular clusters in 'Clust'
    
    Type(Multi_listD),    Intent(in) :: Memb(0:)    ! (0:ns) Final slate cluster data
                                                    ! 0 Case: 
                                                    !   %n = nc = # candidates 
                                                    !   %l = np = # candidates to be elected = total weight
                                                    !             in position units
                                                    !   %qx(2)  = correlation parameters
                                                    !   %vl(mlv). (lv) = # slate clusters with
                                                    !      # ranked or rated candidates <= lv

                                                    ! 1:nsl Case:
                                                    !  %l = # top slate candidates
                                                    !  %m = # slate candidates (for Rating = 1)
                                                    !  %p = # positively rated candidates (for %tx)
                                                    !  %n = # negatively rated candidates (for %tx)
                                                    !  %k = # slate ballots positively correlated with and 
                                                    !       subsequent to this slate ballot
                                                    !  %vl = listing of these by decreasing correlation
                                                    !  %qx = the corresponding correlations

                                                    !  %fsx= slate membership
                                                    !  %fux= slate width
                                                    !  %lt(nc) = candidates in preferential order (from %sx)
                                                    !  %px(nc) = variance vector for %sx
                                                    !  %rx(nc) = mean rating vector with variance normalized weight
                                                    !          = %sx * %ux
                                                    !  %sx(nc) = mean rating vector
                                                    !  %tx(0:nc)= centered (subtract Noise_cor if Rating < 1) 
                                                    !             and normalized mean vector
                                                    !  %ux(nc) = slate ballot weight as variance normalized 
                                                    !          = %fsx / %px

    Type(Multi_listD), Intent(inout) :: Clust(0:)    ! (0:ind) Convergent cluster data, in decreasing order of cluster rating
                                                     ! For 0:
                                                     !   %k = Error code. Invalid results if < 0.
                                                     !   %m = # candidates = nc
                                                     !   %n = # regular clusters = ncl
                                                     !   %l = # clusters, including independent = ind

                                                     !   %px(4):  Objective data: min size, max correlation, 
                                                     !              independents size (frac), final residual
                                                     !   %qx(4):  Objective factors: min size penalty, max correlation penalty, 
                                                     !              independents size penalty, residual penalty
                                                     !   %rx(3):  Clustering objective value(1)= coherent vote(2) * penalty(3)
                                                     !   %sx(ind): Cluster sizes

                                                     !   %lt(ncl)  Mapping of initial clusters to converged clusters

                                                     !   %vl(0:6): Convergence data
                                                     !     0 = original # regular clusters = n0
                                                     !     1 = final # regular clusters
                                                     !     2 = last cycle = # convergence calls - 1
                                                     !     3 = total # cluster merge operations
                                                     !     4 = total # cluster deletion operations
                                                     !     5 = total # centroiding iterations
                                                     !     6 = total # line search points                        
                                                     !   %wt(ns): Sum of slate cluster memberships over all reg
                                                     !            clusters, reduced by the fuzzy set factor,
                                                     !            to be <= 1.0 

                                                     !  %L0(0:4,0:n1) Integer convergence data per update call
                                                     !     0 : Restart type
                                                     !     1 : Success code
                                                     !     2 : # regular clusters
                                                     !     3 : # centroiding iterations
                                                     !     4 : # line search points   
                                                     !  %M2(0:2,0:n1) Real convergence data
                                                     !     0 = Final line search damping parameter
                                                     !     1 = Residual for the updated mean rating vectors
                                                     !     2 = Net line search change
                                                     !  %M1(ncl,ncl) Final cluster correlation matrix, 0 on diagonal
    
                                                     !  %M3(6,ind) Cluster data
                                                     !    1 = size 
                                                     !    2 = density penalty
                                                     !    3 = boundary penalty
                                                     !    4 = density
                                                     !    5 = boundary
                                                     !    6 = width

                                                     !  %L1(-2:nc,0:ind) Candidate ordering data for
                                                     !                    noise zeroed vectors
                                                     !     (-2:0,:) = # top, significant, viable cand's
                                                     !     (1:,1:)  = Decreasing ordering from %T0(:,1:,2)
                                                     !     (1:,0)   = Decreasing ordering from %T0(:,0,1)
                                                     !               
                                                     !  %T0(nc,0:ind,3) Cluster vector data
                                                     !      (:,1:,1) = Cluster mean vectors
                                                     !      (:,1:,2) = Cluster mean vectors with noise zeroing 
                                                     !      (:,1:,3) = Cluster portion vectors
                                                     !      (:,0,1)  = Zrate0 * declining regular portion factor
                                                     !      (:,0,2)  = Zrate0 = cluster averaged noise zeroed mean vectors

                                                     ! For 1:ind:
                                                     !   %k  = # fuzzy reduced slate memberships
                                                     !   %n  = # slate ballots which are partial or full members
                                                     !   %o   : 1 for regular clusters, 2 for independents
                                                     !   %fsx : slate ballot weight averaged fuzzy reduction factor
                                                     !   %fux : cluster width
                                                     !   %sum_wt = new sum of slate memberships

                                                     !   %ls(n)  = list of slate members
                                                     !   %wt(n)  = weighted slate memberships
                                                     !   %rx(n)  = unweighted slate memberships
                                                     !   %px(nc) = variance of %sx
                                                     !   %sx(nc) = mean rating vector
                                                     !   %tx(0:nc)= centered and normalized rating vector
                                                     !   %ux(nc)  : total cluster weight, incorporating slate ballot variance 
                                                     !              normalization (Memb%ux * %rx)
    
    Type(Multi_listR), Intent(inout) :: Clust_set(:)  ! (Ninit) Final cluster sets and associated convergence data
                                                      ! %k   = best initial cluster set that converged to this cluster set
                                                      ! %l   = # clusters in the set, including independents = 'ind' = ncl + 1
                                                      ! %m   = # initial cluster sets that converged to this cluster set
                                                      ! %ls(m)= list of those initial sets 
                                                      ! %fux = clustering objective ratio = (q)%rx(1) / (1)%rx(1) 

                                                      ! %lt(ncl) Mapping of initial to converged regular clusters 

                                                      ! %px(6):  Objective data: min size, max correlation, 
                                                      !            independents size (frac), final residual,
                                                      !            Jnorm, Pnorm
                                                      ! %qx(4):  Objective factors: min size penalty, max correlation penalty, 
                                                      !            independents size penalty, residual penalty
                                                      ! %rx(3):  Clustering objective value(1)= coherent vote(2) * penalty(3)
                                                      ! %sx(ind): Cluster sizes

                                                      ! %ux(6): Frac of voters represented by regular clusters
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
                                                      !   0  = restart phase
                                                      !   1  = success code:  2 = fully converged,   1 = reducing functional,
                                                      !                       0 = damped reduction, -1 = convergence failure
                                                      !   2 = # reg clusters
                                                      !   3 = # centroiding iterations
                                                      !   4 = # line search points

                                                      ! %M2(0:2,0:n1) Real convergence data
                                                      !   0 = Final line search damping parameter
                                                      !   1 = Residual for the updated mean rating vectors
                                                      !   2 = Net line search change
    
                                                      ! %M1(ncl,ncl) Cluster correlation matrix, 0 on diagonal

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
                                                      !               
                                                      !  %T0(nc,0:ind,3) Cluster vector data
                                                      !      (:,1:,1) = Cluster mean vectors
                                                      !      (:,1:,2) = Cluster mean vectors with noise zeroing 
                                                      !      (:,1:,3) = Cluster portion vectors
                                                      !      (:,0,1)  = Zrate0 * declining regular portion factor
                                                      !      (:,0,2)  = Zrate0 = cluster averaged noise zeroed mean vectors
!  Local:
    Real(Dblp), Allocatable :: mean_vec(:,:)
    Real(Dblp) :: Jnorm(4), Pnorm(4)
    Integer    :: cl, nc, n0, n1, nl, np, nu, ind

    Call Out ("Enter 'Record_clusters' for # regular clusters",ncl, ln=1)

    nc= Memb(0)%n;  np= Memb(0)%l
    nl= Size(Clust(0)%sx);  n1= Ubound(Clust(0)%L0,2)
    ind= ncl + 1;  nu= nc*ncl

    If (ind /= nl .or. ind < 2) then
      Call Out ("Error in 'Record_clusters': 1st # clusters",ind, &
                "vs 2nd",nl, ln=1)
      Stop
    End if

    Allocate(mean_vec(nc,ncl));  Forall(cl=1:ncl) mean_vec(:,cl)= Clust(cl)%sx
    
!   Assign 'Clust' to 'Clust_set' for export
    
    Clust_set(nd)%l= ind

    If (.not.Associated(Clust_set(nd)%rx)) then
      Allocate(Clust_set(nd)%rx(3), Clust_set(nd)%vl(0:6), &
               Clust_set(nd)%px(6), Clust_set(nd)%qx(4), Clust_set(nd)%ux(6))
    End if

    If (Associated(Clust_set(nd)%sx)) then
      DeAllocate(Clust_set(nd)%sx, Clust_set(nd)%lt,        &
      Clust_set(nd)%L0, Clust_set(nd)%L1, Clust_set(nd)%M1, &
      Clust_set(nd)%M2, Clust_set(nd)%M3, Clust_set(nd)%T0)
    End if 
        
    n0= Size(Clust(0)%lt)

    Allocate(Clust_set(nd)%sx(ind), Clust_set(nd)%lt(n0),        &
      Clust_set(nd)%L0(0:4,0:n1), Clust_set(nd)%L1(-2:nc,0:ind), &
      Clust_set(nd)%M2(0:2,0:n1), Clust_set(nd)%M3(8,ind),       &
      Clust_set(nd)%T0(nc,0:ind,3))

    Call Representation_levels (np,parm_rep, Real(Memb(1:)%fsx), &
                                Real(Clust(0)%wt), Clust_set(nd)%ux)

    Clust_set(nd)%lt= Clust(0)%lt
    Clust_set(nd)%rx= Real(Clust(0)%rx)
    Clust_set(nd)%vl= Clust(0)%vl
    
    Clust_set(nd)%M3(7,:)= Clust(1:ncl)%k / Real(Clust(1:ncl)%n)
    Clust_set(nd)%M3(8,:)= Clust(1:ncl)%fsx
         
    Clust_set(nd)%px(:4)= Real(Clust(0)%px)
    Clust_set(nd)%px(5) = Real(Jnorm(4)) 
    Clust_set(nd)%px(6) = Real(Pnorm(4))
    Clust_set(nd)%qx    = Real(Clust(0)%qx)

    Clust_set(nd)%L0= Clust(0)%L0
    Clust_set(nd)%L1= Clust(0)%L1
    
    Clust_set(nd)%M2= Real(Clust(0)%M2)
    Clust_set(nd)%M3(:6,:)= Real(Clust(0)%M3)
    Clust_set(nd)%sx= Clust(0)%sx
    Clust_set(nd)%T0= Real(Clust(0)%T0)

    If (ncl > 1) then
      Allocate(Clust_set(nd)%M1(ncl,ncl)); Clust_set(nd)%M1= Real(Clust(0)%M1)
    End if
    
  End Subroutine Record_clusters
        
  Subroutine Consolidate_clusters (ncl, corr_min,max_corr, Value,Corr, &
                                   n1,sing,n2,pair, nq,Cliq_set)
  
!   Lists maximal sets of disjoint cliques of well-correlated clusters 
!   = cliques of clusters wrt corr_min, either singletons, pairs, or
!   maximal cliques.

!   These are to be subsequently consolidated by a weighted averaging of 
!   cluster mean vectors in 'Average_consolidated'.
  
    Integer,     Intent(in) :: ncl          ! # elemental clusters
    Real,        Intent(in) :: corr_min     ! Minimum correlation permitted 
                                            !   among the clusters of a clique
    Real,        Intent(in) :: max_corr     ! Maximum off-diagonal correlation
    Real,        Intent(in) :: Value(:)     ! (ncl)     Cluster values (= weights)
    Real,        Intent(in) :: Corr(:,:)    ! (ncl,ncl) Cluster correlation matrix, 
                                            !           diagonal not used
    
    Integer,    Intent(out) :: n1           ! # singleton component cliques
    Integer,    Intent(out) :: sing(:)      ! (n1<=ncl) List of these singletons 
    Integer,    Intent(out) :: n2           ! # pair component cliques
    Integer,    Intent(out) :: pair(:,:)    ! (2,n2) List of these pairs
    
    Integer,    Intent(out) :: nq           ! # cliques from components of size >= 3
    Type(Set_list), Pointer :: Cliq_set(:)  ! (nq) Disjoint cliques of clusters of strong correlation
                                            !   n = # clusters in the clique
                                            !   p = graph component of the clique 
                                            !   set(n) = list of the clusters in the clique
                                            !   val(n) = membership weighted values of these
                                            !   mbr(n) = memberships of these in the clique, 
                                            !            summing to 1
                                            !   svl    = sum(val)
! Local:
    Type(Set_list),  Pointer :: Cliq(:)  ! (mq) Lists maximal cliques 'q' of clusters 
                                         !      from components of size >= 3

    Type(Set_list),  Pointer :: Cp3(:)   ! (n3) Components of size >= 3
    Integer,  Allocatable :: sg(:), pr(:,:), intr(:), stq(:), lsq(:)
    
    Logical :: rms(ncl)
    Real    :: svl(ncl), pvl(ncl)
    Integer :: key(ncl)
    Integer :: i, j, m, n, q, c1, cp, i1, i2, i3, mq, ms, ni, n3

    If (pr_out > 1) Call Out ("Enter 'Consolidate_clusters'")
    
!   Initialization
    
    n1= 0;  n2= 0;  nq= 0;  sing= 0;  pair= 0
    
!   Test for at least one clique
    
    If (max_corr < corr_min) then
      n1= ncl;  sing(:n1)= "ID";  Return
    End if
    
    If (pr_out > 1.5) then
      Call Out ("Cluster values",Value)
      Call Out (0,"Cluster correlation matrix",Corr)
    End if
    
!   List the 'n3' components Cp3. For each component, list
!   the top maximval cliques 'Cliq' of well-correlated
!   clusters (correlation > corr_min), ordered by decreasing
!   connectivity weighted sum of cluster sizes.

!   These cliques may overlap, so the following 'Clique_loop' 
!   selects a subset of disjoint cliques
    
    Call Best_cliques (5,corr_min,Value,Corr, n1,sing,svl, &
                       n2,pair,pvl, n3,Cp3, Cliq)
    If (n3 <= 0) Return

    mq= Size(Cliq);  ms= Maxval(Cp3%n);  Allocate (lsq(mq), stq(ms), Intr(ms))
    c1= 0;  nq= 0;  lsq= 0;  stq= 0;  key= 0
    
    Clique_loop : Do q= 1,mq  ! List of cliques in decreasing order in each
                              ! component - by connectivity weighted value
      cp= Cliq(q)%p ! Component of clique 'q'

      If (cp > c1) then
        c1= cp;  nq= nq + 1;  lsq(nq)= q
        n= Cliq(q)%n;  stq= 0;  stq(:n)= Cliq(q)%set
      Else  ! Record clique in the same component only if disjoint
        i= Subset(Cliq(q)%set, stq(:n), ni,Intr)

        If (ni == 0) then ! Disjoint: Add the new clique Cliq(q)%set to the 
                          ! union 'stq' of the prior cliques in this component
          nq= nq + 1;  lsq(nq)= q;  m= n + 1  
          n= n + Cliq(q)%n;  stq(m:n)= Cliq(q)%set
          Call Sort (.true.,stq(:n),key(:n)) 
        End if
      End if
    End do Clique_loop

    Allocate (Cliq_set(nq));  Call Copy_Set_list (lsq(:nq), Cliq, Cliq_set)
    
  ! Identify any remaining clusters and add them as singletons
  
    rms= .true.;  rms(sing(:n1))= .false.
    rms(pair(1,:n2))= .false.;  rms(pair(2,:n2))= .false.
    Forall(q=1:nq) rms(Cliq_set(q)%set)= .false.

    Call List_of_true (rms, n,sing(n1+1:));  n1= n1 + n
    
    If (pr_out > 1.5) then
      If (n1 > 0) Call Out ("Singleton components",sing(:n1))
      If (n2 > 0) Call Out (-1,"Pair components",pair(:,:n2))
      If (n3 > 0) then
        Call Out ("General components",n3, "Disjoint cliques",nq)
        Do q= 1,nq
          Call Out ("General clique",Cliq_set(q)%set)
          Call Out ("of component",Cliq_set(q)%p, "with value",Cliq_set(q)%svl)
          Call Out ("Membership weighted values",Cliq_set(q)%val)
          Call Out ("Memberships",Cliq_set(q)%mbr)
        End do
      End if
    End if
    
    Call DeAlloc_set_list_ar (Cp3);   DeAllocate (Cp3)
    Call DeAlloc_set_list_ar (Cliq);  DeAllocate (Cliq)

  End Subroutine Consolidate_clusters

  Subroutine Representation_levels (np,parm,sl_wt, sum_mbr, frac_rep)

!   Analyze the success of the representation produced 
!   by the clustering by computing the fraction 
!   of voters represented by regular clusters, 
!   when that representation exceeds 'parm'

    Integer,  Intent(in) :: np           ! # candidates to be elected
    Real,     Intent(in) :: parm(:)      ! (6) 0 to 1 fractions
    Real,     Intent(in) :: sl_wt(:)     ! (ns) Slate cluster weights
    Real,     Intent(in) :: sum_mbr(:)   ! (ns) Sum of slate cluster memberships over all reg
                                         !      clusters, reduced by the fuzzy set factor,
                                         !      to be <= 1.0 

    Real,    Intent(out) :: frac_rep(:)  ! (6) Frac of voters represented by regular clusters
                                         !     when that representation exceeds 'parm'
! Local:
    Integer :: i, nf

    nf= Size(parm)
    Do i= 1,nf
      frac_rep(i)= Sum(sl_wt, sum_mbr > parm(i)) / np
    End do

    If (pr_out > 1) then
      Call Out ("Frac of voters represented by reg clusters",frac_rep)
      Call Out ("when that representation exceeds",parm)
      Call Out ("Frac represented fully by independents", 1.0 - frac_rep(1))
    End if
  End Subroutine Representation_levels
  
  
   Subroutine Cluster_match1 (iq,nc,ncl, iobj,cur_dat,nd, pr_ncl,pr_dat, match)

!     Test the current cluster set to see if it matches a prior set. 
!     If not, add it to the list of prior sets. Must match the candidate
!     and cluster ratings.

      Integer,    Intent(in) :: iq            ! initial cluster set index
      Integer,    Intent(in) :: nc            ! # candidates
      Integer,    Intent(in) :: ncl           ! # clusters in the current cluster set

      Integer,    Intent(in) :: iobj          ! Current cluster set objective * 100
      Real,       Intent(in) :: cur_dat(:)    ! (nd) Current cluster set objective data:
                                              !      coherent vote / 10, min cluster size, max  
                                              !      cluster correlation, independents size / np
      
      Integer, Intent(inout) :: nd            ! Current # stored cluster sets
      Integer, Intent(inout) :: pr_ncl(:,:)   ! (3,Ninit) Stored # regular clusters in cluster set. 
                                              !   1 = # regular clusters, 2 = current # matches
                                              !   3 = cluster set objective * 100
      Real,    Intent(inout) :: pr_dat(:,:)   ! (nd,Ninit) Stored cluster set objective data:
                                              !     coherent vote / 10, min cluster size, 
                                              !     max cluster correlation, independents size / np
      Integer,   Intent(out) :: match         ! The current set matches prior set 'match' if > 0,
                                              !   else next stored data is at nd= nd+1, 
                                              !   with match= -nd
!   Local:
      Real,   Parameter :: wt(4)= (/ .4, .3, .2, .1 /) ! Weight the highest discrepancy the most
      Real,   Parameter :: required= 0.95

      Real,    Allocatable :: ratio_dat(:), tmp(:), mxr(:)
      Integer, Allocatable :: ls(:)

      Real     :: ratio, max_ratio
      Integer  :: i, k, n, n1, nrd
      
      If (pr_out > 1) Call Out ("Enter 'Cluster_match1' from cluster set index",iq, &
                                "with # recorded sets",nd, ln=1)

      match= -(nd+1) !  Default: No match, with -match = newly recorded cluster set
      
!     Try to match the current cluster set to a prior cluster set
      
      If (nd > 0) then
        nrd= Size(pr_dat,1)
        Allocate (ratio_dat(nrd), tmp(nrd), mxr(nd), ls(nd))
        max_ratio= 0;  tmp= 0;  mxr= 0;  ls= 0;  k= 0

        Prior_loop : Do i= 1,nd
          If (ncl /= pr_ncl(1,i)) Cycle Prior_loop

          ratio_dat= Min(cur_dat, pr_dat(:,i)) / Max(cur_dat, pr_dat(:,i))
          tmp= ratio_dat;  Call Sort(.true.,tmp)
          ratio= Dot_product(wt, tmp(:4));  mxr(i)= ratio

          If (ratio > required) then  ! Match
            k= k + 1;  ls(k)= i
  
            If (ratio > max_ratio) then
              match= i;  max_ratio= ratio  ! Best match so far
            End if
          End if
        End do Prior_loop
      End if

!     Output

      If (pr_out > 1) then
        If (match > 0) then
          Call Out ("Initial cluster set",iq, "matched to prior set",match, ln=1)
          Call Out ("with data ratios to all prior sets",mxr)
        Else 
          Call Out ("Current cluster set",iq, "is a new match, with index",-match, ln=1)
          If (nd > 0) Call Out ("with data ratios to all prior sets",mxr)
        End if
      End if

      If (match > 0) then  ! Match
        n= pr_ncl(2,match);  n1= n + 1;  pr_ncl(2,match)= n1 

        pr_ncl(3,match)= Nint((n * pr_ncl(3,match) + iobj) / Real(n1))
        pr_dat(:,match)= (n * pr_dat(:,match) + cur_dat) / n1

        If (k > 1 .and. pr_out > 1) then
          Call Out ("List of matching cluster sets",ls(:k))
          Call Out ("and data ratios",mxr(ls(:k)))
        End if
      Else                 ! No prior match, record as new
        nd= nd + 1;  pr_ncl(1,nd)= ncl
        pr_ncl(2,nd)= 1;  pr_ncl(3,nd)= iobj
        pr_dat(:,nd)= cur_dat
      End if
    End Subroutine Cluster_match1


   Subroutine Form_clust1 (idist,idm,nc,np, Memb,Mean_rnd, Init,ncs,Clust_set)
   
!    Identify a number of initial cluster sets, based both on random selections of slate ballots 
!    and on heuristic or non-optimal algortihms.
   
!    Then converge each initial cluster set to a final cluster set, if possible, from 1 to 5 clusters per set.
   
!    Retain a new converged cluster set if it differs significantly from the prior converged set but only the best 3 
!    cluster sets for each possible # of clusters per set.

   
     Integer,            Intent(in) :: idist         ! # of the voting district
     Integer,            Intent(in) :: idm           ! index of the domain
     Integer,            Intent(in) :: nc            ! # candidates
     Integer,            Intent(in) :: np            ! # candidates to be elected
     Type(Multi_listD),  Intent(in) :: Memb(0:)      ! (0:ns) Final slate cluster data

     Real,               Intent(in) :: Mean_rnd(0:,:,:)  ! (0:nc,nr,N_cand1) For randomized cluster set initialization.
     
     Type(Multi_listR), Intent(out) :: Init(:)      ! (Ninit) Initial cluster set data
     Integer,           Intent(out) :: ncs          ! # cluster sets
     Type(Multi_listR), Intent(out) :: Clust_set(:) ! (Ninit=>ncs) Final cluster sets and associated 
                                                    !              convergence data
! Local:
     Real,    Allocatable :: tmp(:,:)
     Integer, Allocatable :: nl(:), key(:), map(:,:)
     Integer :: ind, ncl, ncyc, Ninit
     Integer :: i, m, q, cl, nr, ns
     
     Call Out ("Enter 'Form_clust1' for district",idist, "and domain index",idm, ln=1)

     ns= Ubound(Memb,1);  nr= Ubound(Mean_rnd,2);  Ninit= Size(Init)

     Call Init1_clusters (nc,np,ns,nr,Ninit, Mean_rnd, Memb, Init)

     Call Clustering (idist,Ninit,nc, Memb, Init,ncs,Clust_set)

     Call Out ("Clustering output for district",idist, "and domain index",idm, ln=1)
     Call Out ("Total # cluster sets",ncs, "# top cluster sets",cls_n(1), ln=1)
     Call Out ("# viable cluster sets",cls_n(2))

     Allocate (nl(ncs), key(ncs), tmp(3,ncs), map(Ninit,2))

     Do q= 1,ncs
       nl(q)= Clust_set(q)%l - 1
       Forall(i=1:3) tmp(i,q)= Clust_set(q)%rx(i)
     End do

     map= 0
     Do i= 1,Ninit
       map(i,1)= Init(i)%n  
       q= Init(i)%k;  If (q > 0) map(i,2)= nl(q)
     End do

     If (pr_out > 1) then
       Call Out ("Initial to final cluster set mapping",Init%k)
       Call Out (-1,"Initial vs final # clusters",map)

       Call Out ("Fractional clustering objectives",Clust_set(:ncs)%fux)
       Call Out (-1,"Cluster set objectives: clust obj= coherent vote * clust penalty",tmp)
     
       Do q= 1,ncs
         ncl= nl(q)
         Call Out ("For cluster set",q, "with # regular clusters ncl",ncl, ln=1)
         Call Out ("Cluster sizes",Clust_set(q)%sx)
         If (ncl > 1) Call Out (1,"Cluster correlation matrix", Clust_set(q)%M1)
         Call Out ("Obj data (min size, max corr, ind size, resid, J & P norms)",Clust_set(q)%px)
         Call Out ("Objective factors",Clust_set(q)%qx)
         Call Out ("Objective = coherent vote * penalties",Clust_set(q)%rx)

         If (Rating < 1) Call Out ("Ranking noise level",Noise_cor)
         Call Out ("Cluster size averaged zeroed ratings ('Zrate')", Clust_set(q)%T0(:,0,2))
         Call Out ("Zrate * (1- por(:,ind)) for ordering", Clust_set(q)%T0(:,0,1))
         Call Out ("Corresponding ordering", Clust_set(q)%L1(1:,0))
       
         If (pr_out > 1) then
           Call Out (-1,"Cluster size, dens pen, bound pen, dens, boun, width",Clust_set(q)%M3(:6,:ncl))
           Call Out (-1,"Frac slate cluster fuzzy reduced memb, avg reduction factor",Clust_set(q)%M3(7:,:ncl))
           Call Out (-1,"Original cluster mean vectors", Clust_set(q)%T0(:,1:,1))
           Call Out (-1,"The corresponding portions", Clust_set(q)%T0(:,1:,3))
           Call Out (-1,"Corresponding candidate ordering data by cluster", Clust_set(q)%L1(1:,1:))
           Call Out (-1,"Top #, significant #, viable # candidates by cluster", Clust_set(q)%L1(:0,1:))

           Call Out ("Cluster set convergence data:")
           Call Out ("Initial & final # clusters, # convergence cycles", Clust_set(q)%vl(0:2))
           Call Out ("# merges, deletes, # centroidings, # function calls", &
                      Clust_set(q)%vl(3:6))
         
           Call Out (-1,"Convergence: restart, success, ncl, # itr, # fcn calls", &
                       Clust_set(q)%L0)
           Call Out (-1,"Convergence: final t, residual, net line search change", &
                       Clust_set(q)%M2)
         End if
       End do
     End if

   End Subroutine Form_clust1


   Subroutine Init1_clusters (nc,np,ns,nr,Ninit, Mean_rnd, Memb, Init)

!    Compute initial cluster sets for the convergence process. Use a selection of slate cluster sets that are randomly gnerated
!    from slate ballots based on their top 1 or 2 ranked/rated candidates, called candidate cluster sets, alonng with several
!    cluster sets computed by heuristic or non-optimal algortihms.
  
!    These include a "core set" algorithm, a "spectral" algorithm, a "clustering degree" algorithm, 
!    a "triangular coefficient" algorithm, and the STV algorithm.

     
  
     Integer,            Intent(in) :: nc                ! # candidates
     Integer,            Intent(in) :: np                ! # candidates to be elected
     Integer,            Intent(in) :: ns                ! # consolidated slate ballots
     Integer,            Intent(in) :: nr                ! max # regular clusters in an initial cluster set
     Integer,            Intent(in) :: Ninit             ! # initial cluster sets
     
     Real,               Intent(in) :: Mean_rnd(0:,:,:)  ! (0:nc,nr,N_cand1) Mean rating vectors for random initial clusters,
                                                         !   ordered by decreasing cluster size, with sizes summing to 'np'

     Type(Multi_listD),  Intent(in) :: Memb(0:)  ! (0:ns) Slate ballot data
     
!    Output:
     Type(Multi_listR), Intent(out) :: Init(:)   ! (Ninit) Initial cluster set data
      
 ! Local:
     Real, Parameter :: cut(4)= (/0.20, 0.40, 0.60, 0.80/)
     Type(Adjacency) :: Slate_G(ns), Slate_G2(ns)
     Real    :: clust_vec(0:nc,nr+1), clust_memb(ns,nr+1), Coef(ns)
     Real    :: mrg2(3), siz(nr+1), pk_val(2*nr)
     Integer :: cls(nr), cnt(nr+1), key(nr), elc(np,2), peak(2*nr)

     Logical :: ReOrd
     Real    :: max_deg, tot_wt, wt_fac
     Integer :: mc, n1, nf, sl, ng0, ind, ncl, nr1, sum_deg
     Integer :: c, j, n, q, cl, iq

     Call Out ("Enter 'Init1_clusters' to determine the initial cluster sets")
     
     ng0= nGen_clust;  mrg2= Cnv_mrg;  mrg2(3)= np + 1

!    Compute the initial 'candidate' cluster sets
       
     Cand_loop : Do q= 1,N_cand1

       cnt= 0;  siz=0;  clust_vec(:,:nr)= Mean_rnd(:,:,q);  nr1= nr + 1

       Call Clust_matrix (q,nc,np,ns,nr,Noise_cor,Memb(1:), clust_vec, &
                          cnt,siz, clust_memb)

       Call Merge_matrix (nr,mrg2, clust_memb(:,:nr), &
                          cnt(:nr),siz(:nr), n1,cls(:nr))
       ind= n1 + 1 
       
       If (n1 < nr) then
         clust_memb(:,ind)= clust_memb(:,nr1)
         cnt(ind)= cnt(nr1);  siz(ind)= siz(nr1)
         clust_memb(:,nr1)= 0;  cnt(nr1)= 0;  siz(nr1)= 0
       End if

       Call Clust_mean (nc,np, Memb(1:),clust_memb(:,:ind), &
                        siz(:ind),clust_vec(:,:ind))

       Allocate(Init(q)%T0(0:nc,0:ind,2));  Init(q)%T0= 0

       Init(q)%n= n1;  Init(q)%T0(:,1:,1)= clust_vec(:,:ind)
       Forall(c=1:nc) Init(q)%T0(c,0,1)= Sum(clust_vec(0,:ind) * clust_vec(c,:ind)) / np

       If (pr_out > 1) then
         Call Out ("For random initial clustering set",q, "final # regular clusters",n1, ln=1)
         Call Out ("with cluster weights",clust_vec(0,:ind))
         If (pr_out > 1) Call Out (-1,"Their mean vectors",clust_vec(1:,:ind))
       End if
     End do Cand_loop

!    Degree graph peak clustering

     n1= 0;  peak= 0
     Call DeAlloc_Adjacency (Slate_G);  Slate_G%sum_wt= Memb(1:)%fsx

     Do j= 1,4
       Call Slate_graph (ns,cut(j), Memb,Slate_G, sum_deg)
       If (sum_deg < 1) Cycle
       Call Degree_graph (Clust_opt, Slate_G, Coef, max_deg)
       If (max_deg < 0.001) Cycle
       Call Peak_clustering (ns,Slate_G, 0.05,Coef, n1,peak,pk_val)
       If (n1 > np) Exit
     End do

     If (n1 > 0) then
       n1= Min(n1,nr);  ind= n1 + 1;  cnt= 0;  siz=0

       Call Clust_matrix (q,nc,np,ns,n1,Noise_cor,Memb(1:), clust_vec(:,:ind), &
                          cnt(:ind),siz(:ind),clust_memb(:,:ind), peak(:n1))

       Allocate(Init(q)%T0(0:nc,0:ind,2));  Init(q)%T0= 0

       Init(q)%T0(:,1:,1)= clust_vec(:,:ind)
       Forall(c=1:nc) Init(q)%T0(c,0,1)= Sum(clust_vec(0,:ind) * clust_vec(c,:ind)) / np

       tot_wt= Sum(siz(:ind)) 

       If (pr_out > 1) then
         Call Out ("For the degree graph peak clustering set",q, &
                   "# regular clusters",n1, ln=1)
         Call Out ("with total cluster weight",tot_wt)
         Call Out ("from cluster weights",siz(:ind))
         Call Out ("Peak coefficient values",pk_val(:n1))
         If (pr_out > 1.5) Call Out (-1,"and mean vectors",clust_vec(1:,:ind))
       End if
     End if

     Init(q)%n= n1
     
!    Compute and store the trianglar coefficient initial cluster set

     q= q + 1;  n1= 0;  peak= 0
     Call DeAlloc_Adjacency (Slate_G);  Slate_G%sum_wt= Memb(1:)%fsx

     Do j= 1,4
       Call Slate_graph (ns,cut(j), Memb,Slate_G, sum_deg)
       If (sum_deg < 1) Cycle
       Call Triangle_graph (Clust_opt, Slate_G, Coef, max_deg)
       If (max_deg < 0.001) Cycle
       Call Peak_clustering (ns,Slate_G, 0.05,Coef, n1,peak,pk_val)
       If (n1 > np+1) Exit
     End do  
     
     If (n1 > 0) then
       n1= Min(n1,nr);  ind= n1 + 1;  cnt= 0;  siz=0

       Call Clust_matrix (q,nc,np,ns,n1,Noise_cor,Memb(1:), clust_vec(:,:ind), &
                          cnt(:ind),siz(:ind), clust_memb(:,:ind), peak(:n1))

       Allocate(Init(q)%T0(0:nc,0:ind,2));  Init(q)%T0= 0

       Init(q)%T0(:,1:,1)= clust_vec(:,:ind)
       Forall(c=1:nc) Init(q)%T0(c,0,1)= Sum(clust_vec(0,:ind) * clust_vec(c,:ind)) / np

       tot_wt= Sum(siz(:ind))

       If (pr_out > 1) then
         Call Out ("For the triangular graph peak clustering set",q, "# clusters",n1, ln=1)
         Call Out ("with total cluster weight",tot_wt)
         Call Out ("from cluster weights",siz(:ind))
         Call Out ("Peak coefficient values",pk_val(:n1))
         If (pr_out > 1.5) Call Out (-1,"and mean vectors",clust_vec(1:,:ind))
       End if
     End if

     Init(q)%n= n1

     If (pr_out > 1) Call Out ("Initial cluster sets, elect candidates by Rem_frac, DHondt, & optimal methods")

     Init_loop : Do iq= 1,N_init1
       ncl= Init(iq)%n;  ind= ncl+1;  If (ncl < 1) Cycle Init_loop

       Allocate(Init(iq)%lt(ncl), Init(iq)%Q0(np,0:2,4), &
                Init(iq)%M1(7,4), Init(iq)%M2(ncl,ncl))
       Init(iq)%lt= "ID";  Init(iq)%Q0= 0;  Init(iq)%M1= 0;  Init(iq)%M2= 0

       Call Rem_frac (iq,nc,ind, Init(iq)%T0(:,1:,1), Init(iq)%Q0(:,:,1))
       Call DHondt (iq,nc,ind, Init(iq)%T0(:,1:,1), Init(iq)%Q0(:,:,2))

       Call Elect_fit (iq,np,nc,ind, Memb(1:),Init(iq)%T0(:,:,1), &
                       Init(iq)%Q0(:,0,:2), Init(iq)%Q0(:,:,3),   &
                       Init(iq)%T0(1:,:,2),Init(iq)%M2,Init(iq)%M1(:,:3)) 
     End do Init_loop
    
   End Subroutine Init1_clusters

End Module Clusters3
 