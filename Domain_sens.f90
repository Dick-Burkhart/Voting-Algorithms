Program Domain_sens   

!  This program implements uses PR_clustering to test the sensitivity of clustering 
!  and non-clustering methods to perturbations or changes in the domain as a subset 
!  of the full set of candidates 1...nc. This is related  to Nic Tideman's concept 
!  for a pairwise tournament between competing sets of candidates, with winning set 
!  to be elected. Except that my clustering full objective functions are used to compare 
!  the pairs of sets instead of STV concepts. 

!  The domains may be either subsets of size nc – 2 or subsets of np+1 through 2*np, 
!  where np candidates are to be elected. The second method is like Tideman’s Tournament 
!  algorithm, where 2 competing sets of size np are compared on the domain consisting of their union. 

!  The point of a sensitivity test is that the best elected set may differ when it is computed 
!  for a domain of candidates that is a proper subset of the full set because the proportionality 
!  may be affected. To demonstrate this I compute the full objective function for the elected 
!  set on the domain ‘D’ in question. 
    
!  To compare 2 subsets A’ and ‘B’ of ‘D’ I use the ratio Rd(A,B) = Fd(A) / Fd(B) of their 
!  full objective functions Fd for the top cluster set of ‘D’. Then the overall objective 
!  function for ‘A’ becomes a weighted average Sd(A) of a function of Rd(A,B) over all 
!  distinct sets ‘B’ over all domains ‘D’ containing both ‘A’ and B’, weighted by the 
!  cluster set objective for the top cluster set of ‘D’.

!  However, it is surprisingly effective to just use N(A) = number of domains for which ‘A’ 
!  is the top elected set. In fact, to determine the “winner”, we  restrict our search 
!  to sets ‘A’ which are the top elected set (= the set with the best fitness) for at 
!  least one domain. Then Sd(A) could be used as a tie breaker.
    
!  Note that computational load may limit the # districts that can be processed on one run    
    
   Use Clusters0  ! Contains run parameters
   Use Clusters1  ! Contains "Read_ballots0" and "Domain_ballots" and thier subroutines
   Use Clusters2  ! Contains "Consol1_ballots" and its subroutines
   Use Clusters3  ! Contains "Form_clust1" and its subroutines
   Use Clusters4  ! Contains "Possible_electeds" and related subroutines
   Use Clusters6  ! Contains subroutines for input and statistics

   Use Newton_operators
   Use Factorials
   Use Types
   Use Precisn
   Implicit None

!  Data structcures for Domain_input and Read_ballots0 (ballot input)
   
   Character(18), Pointer :: District(:)=>Null()   ! (mx_Dist) List of voting districts to be processed (input)
   Character(3),  Pointer :: Party(:,:)=>Null()    ! (mxCand,mx_Dist) Party affiliation of each candidate for each district, if any
   Integer,       Save    :: Rnd_seed(2,mx_Dist)   ! Random seed for each district, to be use for 'Mean_rnd'
   Real,          Pointer :: Noise_pr(:)=>Null()   ! Prior computed values of Noise_cor by district

   Integer              :: nb0           ! Initial # ballots for a district
   Real,    Allocatable :: wtb0(:)       ! (nb0) Initial ballot weights, summing to np
   Integer, Allocatable :: ballot0(:,:)  ! (0:mrc,nb0) Initial ballots. (1:n,b) = candidates in preferential order
                                         !                              (0,b)   = # candidates ranked or rated     
   Integer          :: nb                    ! # ballots in a subdomain
   Real,    Pointer :: wtb(:)=>Null()        ! (nb) ballot weights for the subdomain
   Integer, Pointer :: ballot(:,:)=>Null()   ! (0:mr,nb) ballots for the subdomain

!  Data structcures for Consol1_ballots (consolidation to slate ballots)
   
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

!  Data structcures for "Form_clust1" (converge to cluster sets)
     
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

                                           ! %M0(nr,ind) Cluster slate cluster membership weights
                                           ! %M1(7,4) Objective data for the electeds %Q0(:,0,:)
                                           ! %M2(ncl,ncl) Cluster correlation matrix, 0 on diagonal

                                           ! %T0(0:nc,0:ind,2) (1:,1:,1) = Candidate mean ratings  
                                           !      for each cluster, with cluster size at (0,1:,1)
                                           !      (1:,1:,2) = Cluster portions with averaged 
                                           !                  noise zeroed ratings at (1:,0,2)

     Type(Multi_listR) :: Clust_set(N_init)  ! (ncs) Final cluster sets and associated convergence data
                                             ! %k   = best initial cluster set that converged to this 
                                             !        set, based on the optimal elected before convergence
                                             ! %l   = # clusters in the set, including independents = 'ind' = ncl + 1
                                             ! %m   = # initial cluster sets that converged to this cluster set
                                             ! %ls(m)= list of those initial sets 
                                             ! %fux = clustering objective ratio = (q)%rx(1) / (1)%rx(1) 

                                             ! %lt(ncl) Mapping of initial to converged regular clusters 

                                             ! %px(6)  Objective data: min size, max correlation, 
                                             !         independents size (frac), final residual,
                                             !         Jnorm, Pnorm
                                             ! %qx(4)  Objective factors: min size penalty, max correlation penalty, 
                                             !         independents size penalty, residual penalty
                                             ! %rx(3)   Clustering objective value(1)= coherent vote(2) * penalty(3)
                                             ! %sx(ind) Cluster sizes

                                             ! %ux(6)  Frac of voters represented by regular clusters
                                             !         when that representation exceeds 'Parm_rep'

                                             ! %vl(0:6): Convergence data
                                             !   0 = original # regular clusters = n0
                                             !   1 = final # regular clusters
                                             !   2 = final success code
                                             !   3 = total # cluster merge and deletion operations
                                             !   4 = # convergence calls
                                             !   5 = total # centroiding iterations
                                             !   6 = total # line search points

                                             ! %L0(0:4,0:iv) Integer convergence data per update call
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

!  Data structcures for "Possible_electeds" (list possible elected sets for each domain)
     
   Integer, Allocatable :: ncp(:)        ! (n1:n2) N_subsets(nc,np)
   Integer, Allocatable :: Subi(:)       ! (msub) Subset index in 'subU' for each tested subset
   Integer, Allocatable :: subU(:,:)     ! (0:n2,nsub) List of subsets of (1...n0) of size <= n2
                                         !   where (0,i)= size of subset 'i' & n2 = size of A U B
   Real,    Allocatable :: cls_obj(:)    ! (msub) Top cluster set objective
   Real,    Allocatable :: tot_wt(:)     ! (msub) Reduced total weight for each domain

   Type(Multi_listR), Allocatable :: Clust_elc(:) ! (msub) Data for sets of electeds associated to the
                                                  !        top cluster set of each domain
                                                  ! %n= nc = # candidates in the domain
                                                  ! %l= ind = # clusters, including independents
                                                  ! %m= nf  = # possible elected sets

                                                  ! %fsx      = Top fitness value
                                                  ! %rx(nf)   = Fitness ratios
                                                  ! %L1(np,nf)= Sets of electeds (original candidates)
                                                  ! %M1(nf,nf)= Dominance factors between the %L1 sets:
                                                  !    0 = no dominance. Full dominance of set 'k' over 'l' when
                                                  !    %M1(l,k) = 1, with k < l, i.e., fitness(l) < fitness(k)
                                                  
                                                  ! %ls(nc) = Domain candidates
                                                  ! %sx(nc) = Domain average candidate scores
                                                  ! %vl(nc) = Domain reordering by average candidate scores

!  Data structcures for 'Domain_sens':

   Real, Parameter :: eps= 0.001   
   Character(1)    :: c1(9)= (/"1","2","3","4","5","6","7","8","9"/)

   Integer            :: nc_lim= 9    ! Max # candidates - to reduce the computational load
   Real               :: full         ! Input: Limit parameter for the rising cosine dominance function

   Integer,      Save :: Sens_opt     ! Input: 1 : Domains = union of 2 possible elected sets of size 'np'
                                      ! < 0 : Domains = all subsets of size n0-2 of a domain of size n0
                                      !  -1 : Process by comparing dominance values of top sets (dmv => Sdv)
                                      !  -2 : Process by summing dominance values for each top set (Sdv)

   Integer, Save :: Tr_mat(-1:nMt,2)  ! Adds up, over districts, matches of the top elected set by 
                                      !   (-1,1) domain count top set = dominance formula top set
                                      !   (0,1)  clustering top set   = dominance formula top set
                                      !   (>0,1) non-clustering top set = dominance formula top set
                                      !   (0,2)  clustering top set   = domain count top set
                                      !   (>0,2) non-clustering top set = domain count top set

   Integer, Save :: Elc_pr(5,0:nMt,40)   ! Previously computed top elected sets, read from 'Top_dat1.txt',
                                         ! with (0) = clustering and (1:nMt) = non-clustering 

   Integer, Allocatable :: elc1(:,:)     ! (tp,2)  (:,1) = top_elc(:,dm1), (:,2) = top_elc(:,tp1)
   Logical, Allocatable :: dom0(:,:)     ! (n0,md) True on a domain
   Logical, Allocatable :: dom_k(:)      ! (n0) True on a subset of size 'np'
   Logical, Allocatable :: dom_lk(:)     ! (n0) True on a subset of size 'np'
   Real,    Allocatable :: dmv(:,:)      ! (tp,tp) Dominance value between 2 top sets
   Real,    Allocatable :: Sdv(:)        ! (tp) Summed dominance values for each top set
   Real,    Allocatable :: tpN(:)        ! (tp) top_elc(0,:) + eps * Sdv. To compute keyT
   Integer, Allocatable :: keyD(:)       ! (tp) Reordering key for Sdv
   Integer, Allocatable :: keyT(:)       ! (tp) Reordering key for tpN

   Integer, Allocatable :: top_elc(:,:)  ! (0:np,md) Top sets of domains (1:np,:) and their # (0,:)
   Integer, Allocatable :: mp(:)         ! (md)  Mapping from domains to their top sets top_elc(1:,:)
   Integer, Allocatable :: domain(:)     ! (nc)  

   Integer :: nf    ! # possible elected sets for which fitness is computed
   Integer :: tp1   ! Top set for domain count 'top_elc'
   Integer :: dm1   ! Top set for dominance values 'Sdv'

   Integer :: ncs    ! Total # cluster sets converged to 
   Integer :: nc     ! # candidates
   Integer :: np     ! # candidates to be elected
   Integer :: mr     ! Max # ranked or rated candidates
   Integer :: mrc    ! Possibly reduced 'mr' so that mrc= Min(mr,nc)
   Integer :: nsl    ! # slate ballot clusters

   Integer :: me     ! # sets of possible electeds as computed
   Integer :: ncl    ! ncl = # regular clusters in a cluster set
   Integer :: ind    ! = ncl + 1 = total # clusters, including independents, 
                     ! = last 'cluster in the set
   Integer :: ncyc   ! # convergence updates = index of the final cycle
   Integer :: idist  ! Current district index

   Integer :: dst1, dst2 ! First and last district to process
   Integer :: n12(2)     ! (1) = np, (2) = nonC
   Integer :: Non(nMt)   ! The index of the non-clustering elected set 
                         !   for each method, <= nonC

   Real,    Allocatable :: tmp(:)
   Integer, Allocatable :: jtp(:)

   Integer :: El0(5,0:3), El1(5,0:2)

   Real    :: sm
   Integer :: msub, nonC, nsub, n1sub
   Integer :: m1, m2, md, ml, n0, n1, n2, n5, nd, nl, kj, lj, ne, nr, tp, ios, nr0
   Integer :: e, i, j, k, l, m, n, p, q, id, iq, j1, k1
    
   Rating= -1;  GN_rating= -1;  Dot_fac= Dot_fac_rank 
   GN_max_pt= Max_ptD;  GN_min_wt= Min_wt;  GN_min_mx= Min_mx
   ngc(1)= nGen_clust;  ngc(2)= nGen_clust_CM  
   ngc(3)= nGen_clust_NT

   Call Domain_input (District,Party,Noise_pr, pr_out, dst1,dst2, Sens_opt,full)

   Open(7, File='RndSeed.txt', IOstat=ios, Status='Old', Action='Read')
     Read(7,*);  n= mx_Dist/4;  j= 1

     Do i= 1,n
       k= Min(j+3,mx_Dist);  Read(7,'(4(2X,2I11))') (Rnd_seed(:,id), id=j,k)
       If (k == mx_Dist) Exit;  j= k + 1
     End do
   Close(7)

   Open(7, File='Top_dat1.txt', IOstat=ios, Status='Old', Action='Read')
     Read(7,*);  Elc_pr= 0
     Do id= 1,mx_Dist
       Read(7,'(3I3,2X,7I2, 2X,5I3, 2X,3(2X,5I3))')  k, n12, Non, (El0(:,j),j=0,3)
       Read(7,*)
       n5= Min(n12(1),5);  Elc_pr(:n5,0,id)= El0(:n5,0)  
       Forall(i=1:nMt) Elc_pr(:n5,i,id)= El0(:n5,Non(i))
     End do
   Close(7)

!  Loop over the specified range of districts

   District_loop : Do idist= dst1,dst2

     Call Read_ballots0 (Trim(District(idist)),idist,nc_lim, nc,np,mr,nb, wtb,ballot)

     Noise_cor= Noise_pr(idist);  GN_rate_adj= Noise_cor
     Noise_por= Por_fac * Noise_cor;  nb0= nb;  n0= nc

     If (Allocated(wtb0)) DeAllocate(wtb0, ballot0);  Allocate(wtb0(nb), ballot0(0:mr,nb))
     wtb0= wtb(:nb);  ballot0= ballot(:,:nb)

!    Do the comparisons for each subset of (1...n0) of size of n1 to n2

     If (Sens_opt > 0) then
       n1= np + 1;  n2= Min(2*np,n0)
     Else
       n1= Max(n0-2,np+1);  n2= n1
     End if 

     nsub = N_sub_sum(n0,n2) - 1 
     n1sub= N_sub_sum(n0,n1-1);  msub= nsub+1 - n1sub

     If (Allocated(ncp)) DeAllocate(ncp, subU, Subi, tot_wt, cls_obj, Clust_elc)
     Allocate (ncp(n1:n2), subU(0:n2,nsub), Subi(msub), tot_wt(msub), &
               cls_obj(msub), Clust_elc(msub))

     Forall(nc=n1:n2) ncp(nc)= N_subsets(nc,np)
     j= 0;  k= 0;  subU= 0;  Subi= 0;  cls_obj= 0

     Call List_subsets (1,n0,n2,0, k,subU)  ! All non-empty subsets of size <= n2. 
                                            ! Note: k= nsub

     Call Out ("For district",idist, "with total # candidates",n0, ln=1)
     Call Out ("For min domain size",n1,"to max domain size",n2)
     Call Out ("For # electeds",np, "total # domains",msub)
     Call Out ("# domains of each size, min to max",ncp)

!    Loop over all domain subsets of 1...nc, computing objective data
!    for each of these domains

     Domain_loop0 : Do i= 1,nsub
       nc= subU(0,i);  If (nc < n1) Cycle Domain_loop0

       j= j + 1;  ne= ncp(nc);  mrc= Min(mr,nc)

       Call Domain_ballots (np,mrc,nc,n0,subU(1:nc,i), nb0,wtb0,ballot0, &
                            nb,tot_wt(j),wtb,ballot(:mrc,:))

       nr0= Min(np+2,7);  Allocate(Mean_rnd(0:nc,nr0,N_cand1))
       Call Random_seed (Put = Rnd_seed(:,idist))

       Call Consol1_ballots (nc,np,mrc,nb, wtb(:nb),ballot(:mrc,:nb), &
                             nr,nsl,Memb, Mean_rnd)
       Call Form_clust1 (idist,j,nc,np, Memb,Mean_rnd, Init(:N_init1), ncs,Clust_set)

       DeAllocate(Mean_rnd);  Call DeAlloc_Multi_list_ar (Init);  
       Call DeAlloc_Multi_list_ar (Memb);  DeAllocate(Memb)

!      For each domain 'j' of size 'nc', list all the possible elected sets
!      by decreasing fitness

       ind= Clust_set(1)%l;  Subi(j)= i;  cls_obj(j)= 0.10 * Clust_set(1)%rx(1)

       Call Possible_electeds (np,nc,ind,ne,full, subU(1:nc,i), Clust_set(1), nf,Clust_elc(j))
       Call DeAlloc_Multi_list_ar (Clust_set)
     End do Domain_loop0

    md= j  ! Compare to 'msub'

    If (Allocated(top_elc)) DeAllocate (top_elc, dom0, mp, jtp)

    Allocate(top_elc(0:np,md), dom0(n0,md), mp(md), jtp(md))
    top_elc= 0;  k= 0;  dom0= .false.;  mp= 0;  jtp= 0

!   For each of these domains, output the possible elected subsets 
!   and their fitness values, recording in top_elc how many times
!   a subset tops for a domain = domain count top set
   
    Domain_loop1 : Do j= 1,md
      nc= Clust_elc(j)%n;  n= Min(Clust_elc(j)%m,10);  i= Subi(j)
      dom0(subU(1:nc,i),j)= .true.

      Do m= 1,k
        If (All(Clust_elc(j)%L1(:,1) == top_elc(1:,m))) then
          top_elc(0,m)= top_elc(0,m) + 1;  mp(j)= m;  Exit
        End if
      End do

      If (m > k) then  ! New top set
        k= m;  mp(j)= k;  jtp(k)= j
        top_elc(0,k)= 1;  top_elc(1:,k)= Clust_elc(j)%L1(:,1)

        If (pr_out >= 1) then
          Write (8,*);  Write (8,'(A,I4, A,I4, A,F7.3)')  "New top set",k, "  from domain",j, "  with objective",cls_obj(j)
          Write (8,'(A,10I3)')  "  domain ", (subU(l,i), l=1,nc)
          Write (8,"(A, 10("//c1(np)//"I3, 2X))") "  ord elc sets  : ",  ((Clust_elc(j)%L1(l,p), l=1,np), p=1,n)
          Write (8,'(A, 10F7.3)')                 "  fitness ratios: ",  (Clust_elc(j)%rx(p), p=1,n)
        End if
      Else
        If (pr_out >= 1) then
          Write (8,*);  Write (8,'(A,I4, A,I4, A,F7.3)')  "Top set for domain",j,"  matches prior top set",m, "  with objective",cls_obj(j)
          Write (8,"(A,"//c1(np)//"I3, A,10I3)")  "  top set ", (top_elc(l,m), l=1,np), "  domain ", (subU(l,i), l=1,nc)
        End if
      End if

    End do Domain_loop1

    tp= k;  If (Allocated(dom_k)) DeAllocate(dom_k, dom_lk, keyD, keyT, Sdv, tpN, tmp, elc1, dmv)
    Allocate(dom_k(n0), dom_lk(n0), keyD(tp), keyT(tp), Sdv(tp), tpN(tp), tmp(tp), elc1(np,2), dmv(tp,tp))

    Call Out ("# top sets",tp, "vs # domains",md,ln=1)
    Call Out ("# domains for each top elected set",top_elc(0,:tp))
    Call Out ("Mapping of all domains to the top elected sets",mp)

!   Next compute the dominance formula values Sdv (objective function)
!   for each domain count top set.

    If (Sens_opt > 0) then  ! Domains = unions of sets to be compared by dominance values
      Sdv= 0
      Top_set_loop : Do k= 1,tp
        dom_k= .false.;  dom_k(top_elc(1:,k))= .true.

        Domain_loop2 : Do j= 1,md    ! 'dom_k' must be a subset of dom0(:,j)
          If (Any(dom_k .and. .not.dom0(:,j))) Cycle Domain_loop2

          sm= 0;  nl= 0;  nf= Clust_elc(j)%m;  kj= Locate_set(top_elc(1:,k), Clust_elc(j)%L1)

          Subset_loop : Do l= 1,nf  ! the union 'dom_lk' must equal dom0(:,j)
            dom_lk= dom_k;  dom_lk(Clust_elc(j)%L1(:,l))= .true.
            If (Any(dom_lk /= dom0(:,j))) Cycle Subset_loop

            If (l > kj) then
              sm= sm + Clust_elc(j)%M1(l,kj) * cls_obj(j)
            Else
              sm= sm - Clust_elc(j)%M1(kj,l) * cls_obj(j)
            End if

            nl= nl + 1
          End do Subset_loop

          If (nl > 0) Sdv(k)= Sdv(k) + sm / nl
        End do Domain_loop2
      End do Top_set_loop

    Else if (Sens_opt == -1) then ! Domains of size n0-2. Process by comparing dominance values for top sets
      dmv= 0
      Top_set_loop1 : Do k= 1,tp-1
        dom_k= .false.;  dom_k(top_elc(1:,k))= .true.
        Top_set_loop2 : Do l= k+1,tp
          dom_lk= dom_k;  dom_lk(top_elc(1:,l))= .true.

          Domain_loop3: Do j= 1,md
            If (Any(dom_lk .and. .not.dom0(:,j))) Cycle Domain_loop3

            kj= Locate_set(top_elc(1:,k), Clust_elc(j)%L1)
            lj= Locate_set(top_elc(1:,l), Clust_elc(j)%L1)
            If (kj < 1 .or. lj < 1) Cycle Domain_loop3

            If (lj > kj) then
              dmv(l,k)= dmv(l,k) + Clust_elc(j)%M1(lj,kj) * cls_obj(j)
            Else
              dmv(l,k)= dmv(l,k) - Clust_elc(j)%M1(kj,lj) * cls_obj(j)
            End if
          End do Domain_loop3
        End do Top_set_loop2
      End do Top_set_loop1

      Do k= 1,tp  !   Summed dominance values
        Sdv(k)= Sum(dmv(k+1:,k)) - Sum(dmv(k,:k-1))
      End do

    Else  ! Domains of size n0-2. Compute and sum dominance values for each top set
      Sdv= 0
      Top_set_loop3 : Do k= 1,tp
        dom_k= .false.;  dom_k(top_elc(1:,k))= .true.

        Domain_loop4 : Do j= 1,md    ! 'dom_k' must be a subset of dom0(:,j)
          If (Any(dom_k .and. .not.dom0(:,j))) Cycle Domain_loop4

          sm= 0;  nf= Clust_elc(j)%m;  kj= Locate_set(top_elc(1:,k), Clust_elc(j)%L1)

          Subset_loop1 : Do l= 1,nf
            If (l > kj) then
              sm= sm + Clust_elc(j)%M1(l,kj) * cls_obj(j)
            Else if (l < kj) then
              sm= sm - Clust_elc(j)%M1(kj,l) * cls_obj(j)
            End if
          End do Subset_loop1

          Sdv(k)= Sdv(k) + sm/(nf-1)

        End do Domain_loop4
      End do Top_set_loop3
    End if
    
    If (idist == 12) then
      Call Out ("Domain indices for top sets (1-7)",jtp(:7))
      Call Out ("# domains with same top set",top_elc(0,:7))
      Call Out ("Corresponding dominance values",Sdv(:7))
    End if

  ! Top set by dominance formula vs domain count

    tmp= Sdv;  Call Sort (.false.,tmp,keyD)

    tpN= top_elc(0,:tp) + eps * Sdv;  tmp= tpN
    Call Sort (.false.,tmp,keyT)

    dm1= keyD(1);  tp1= keyT(1)
    elc1(:,1)= top_elc(1:,dm1);  elc1(:,2)= top_elc(1:,tp1)

    Call Out ("For district",idist, ln=1)
    Call Out ("The summed dominance values for each top set",Nint(100*Sdv))
    Call Out ("Their ordering",keyD)
    Call Out ("Domain counts for each top set",Nint(100*tpN))
    Call Out ("Their ordering",keyT)

    Call Out ("Top set index by dominance values",dm1, "by domain count",tp1)
    Call Out (-1,"Top set by dominance values & domain count",elc1)

    n5= Min(np,5)
    Call Out ("Top overall clustering elected set",Elc_pr(:n5,0,idist))
    Call Out (-1,"Non-clustering elected sets",Elc_pr(:n5,1:,idist))

    Call Out ("# domains for each top set by dominance ordering",top_elc(0,keyD))
    Call Out ("# domains for each top set by domain count ordering",top_elc(0,keyT))

    If (Sens_opt == -1 .and. pr_out > 1) then
      Call Out ("Dominance values for original over subsequent top sets")
      Do k= 1,tp-1
        Call Out ("For original top set",k,ln=1)
        Call Out ("Dominance values * 100",Nint(100*dmv(k+1:,k)))
      End do
    End if

!   Tr_mat adds up matches of the top elected set by 
!     (-1,1) domain count top set = dominance formula top set
!     (0,1)  clustering top set   = dominance formula top set
!     (>0,1) non-clustering top set = dominance formula top set

!     (0,2)  clustering top set   = domain count top set
!     (>0,2) non-clustering top set = domain count top set

    If (tp1 == dm1) Tr_mat(-1,1)= Tr_mat(-1,1) + 1  ! Top set by dominance formula = top set by domain count

    Do k= 1,2
      Do i= 0,nMt                                   ! Matches to clustering & non-clustering elected sets 
        If (All(elc1(:n5,k) == Elc_pr(:n5,i,idist))) Tr_mat(i,k)= Tr_mat(i,k) + 1
      End do
    End do

  End Do District_loop

!  Output statistics over districts

  Call Out ("For domain option",Sens_opt,ln=1)
  Call Out ("From district",dst1, "to district",dst2)
  Call Out ("# matches of 'dominance formula top set' to 'domain count top set'",Tr_mat(-1,1))

  Call Out ("# matches of the 'dominance formula top set' to clustering top set",Tr_mat(0,1),ln=1)
  Call Out ("# matches of the 'dominance formula top set' to non-clustering electeds",Tr_mat(1:,1))

  Call Out ("# matches of the 'domain count top set' to clustering top set",Tr_mat(0,2),ln=1)
  Call Out ("# matches of the 'domain count top set' to non-clustering electeds",Tr_mat(1:,2))
  Close(8)

End Program Domain_sens