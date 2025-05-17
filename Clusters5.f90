!            This module contains "Write_summary" and its subroutines. 
    
!  "Write_summary" is called from Clustering_PR to write a summary 
!  of the clustering and the candidate evaluation to the output file
    
Module Clusters5

   Use Clusters0
   Use Clusters_support

   Use Util
   Use Output
   Use Types
   Use Precisn
   Implicit None
   
 Contains

   Subroutine Write_summary (District,idist, np,nc,mr,mrn, ml,ncs,me, & 
                             orig_cand,Party, Clust_set,Clust_elect,  &
                             Elect_dat, Non_clust,Init, dpr)
   
     Character(*), Intent(in) :: District  ! Name of the voting district
     Integer,      Intent(in) :: idist     ! District #
     Integer,      Intent(in) :: np        ! # electeds
     Integer,      Intent(in) :: nc        ! # candidates
     Integer,      Intent(in) :: mr        ! Max # ranked candidates
     Integer,      Intent(in) :: mrn       ! Max # negatively rated candidates

     Integer,      Intent(in) :: ml        ! max # clusters + 1
     Integer,      Intent(in) :: ncs       ! # cluster sets
     Integer,      Intent(in) :: me        ! # possible sets of electeds
     Integer,      Intent(in) :: orig_cand(:) ! (nc) Original candidates
     
     Character(3), Intent(in) :: Party(:)  ! (nc0) Party affiliation of each candidate 
                                           !   for each district, if any (file read)

     Type(Multi_listR), Intent(in) :: Clust_set(:) ! (ncs) Final cluster sets
                                         ! %k   = best initial cluster set that converged to this cluster set
                                         ! %l   = 'ind' = # clusters, including independents
                                         ! %m   = # initial cluster sets that converged to this cluster set
                                         ! %ls(m)= list of those initial sets 
                                         ! %fux = clustering objective ratio = (q)%rx(1) / (1)%rx(1) 

                                         ! %lt(ncl) Mapping of initial to converged regular clusters 

                                         ! %px(6):  Objective data: min size, max correlation, 
                                         !          independents size (frac), final residual,
                                         !          Jnorm, Pnorm
                                         ! %qx(4):  Objective factors: min size penalty, max correlation penalty, 
                                         !              independents size penalty, residual penalty
                                         ! %rx(3):  Objective value(1) = coherent vote(2) * penalty(3)
                                         ! %sx(ind): Cluster weights = sizes

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

                                         ! %M2(0:2,0:iv) Real convergence data
                                         !    0: line search parameter, if > 0
                                         !    1: functional
                                         !    2: Net line search change

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

     Type(Multi_listR), Intent(in) :: Clust_elect(:) ! (ncs) Cluster set data, with possible sets of electeds,
                                                     !       ordered by decreasing fitness 
                                         ! %k = index 'e' = %vl(1) in Elect_dat of the first electeds %Q0(:,0,1) 
                                         ! %l = 'ind' = # clusters, including independents
                                         ! %m = 'nf'= # sets of possible electeds
                                         ! %fsx = full objective for the top set of electeds 
                                         !        %Q0(:,0,1) = Clust(q)%sx(1)
                                         ! %fux = full objective ratio = (q)%fsx / (1)%fsx 

                                         ! %sx(7)  Full objective data, as in Elect_dat(1)%M1 
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

     Type(Multi_listR), Intent(in) :: Elect_dat(:) ! (me) Data for possible sets 'e' of electeds,
                                         !      ordered by 'best' full objective %M1(1,1)
                                         ! %k = The 'best' cluster set 'q' for 'e' = %lt(1)
                                         ! %l = The matching 'eq' in Clust_elect(q)%Q0 for q = %vl(1)
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
     Type(Multi_listR), Intent(in) :: Non_clust  ! Data for 'nMt' selected non-clustering methods 'm'
                                     !           to elect multiple candidates.
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
                                     ! %L1(nMt,2)  (m,1)= Index 'e' in Elect_dat of non-clustering set 'm' 
                                     !             (m,2)= The matching best cluster set 'q'
                                     ! %L2(nonC,ncs) = Location, for each cluster set 'q', in Clust_elect(q)%Q0 
                                     !                of each distinct non-clustering set of electeds

                                     ! %T0(0:nc,0:ind,2) = Cluster data from STV (:,:,1) and DTV (:,:,2).
                                     !    (0,:,:) = cluster sizes, (1:,:,:) = cluster mean rating vectors
                                     !    (1:,0,:) = cluster averaged mean vector
                                     ! %T1(ncl,ncl,2) = Cluster correlations for STV (:,:,1) and DTV (:,:,2),
                                     !    = 'STV_corr' & 'DTV_corr' .
     Type(Multi_listR), Intent(in) :: Init(:)  ! (N_init) Initial cluster set data
                                               ! %k = cluster set converged to, or 0 if none
                                               ! %n = ncl = # regular clusters
                                               ! %lt(ncl) Mapping from initial to converged clusters

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
                                               !      (1:,1:,2) = Cluster portions with averaged 
                                               !                  noise zeroed ratings at (1:,0,2)
     Integer,  Intent(in) :: dpr(0:)           ! (0:mfe) Elected sets whose best full objective 
                                               !         exceeds that of cluster set #1.
                                               !         (0) = # such elected sets.
!  Local
     Character(2) :: c2(25)= (/" 1"," 2"," 3"," 4"," 5"," 6"," 7"," 8"," 9","10", &
                               "11","12","13","14","15","16","17","18","19","20", &
                               "21","22","23","24","25"/)
     Real, Parameter :: D2= 100.0,  D3= 1000.0
     Integer, Allocatable :: nl(:), set_ord(:), orig_ord0(:), elc_sets(:,:,:)
     Real,    Allocatable :: rates(:,:,:), por(:,:)

     Logical :: Covered(me)
     Integer :: elc(np,2), ls(nc), Non_cls_sets(np,nMt,2), cs(4), ob(4)
     Integer :: electeds(np,2)

     Real    :: wt1, wt2, fit, obj_fac, full_obj, cls_obj(3), elc_obj(3) 
     Integer :: ncyc, nSTV, nDTV, N_init, best
     Integer :: ios, in0, ind, me2, mrp, nc0, ncl
     Integer :: iq, k1, l1, n0, n1, n2, n3, n4, nd, nf, nq, nt, p1, p2, qa, qb, qo1, qo2
     Integer :: e, i, j, k, l, m, n, o, p, q, cl
          
     Call Out (" ");  Call Out ("Enter 'Write_summary' for district "//District)

     Call Out ("# candidates",nc, "# to be elected",np)
     Call Out ("# cluster sets",ncs, "# possible elected sets",me)
     mrp= mr - mrn;  N_init= Size(Init)
     nc0= Maxval(orig_cand);  n4= Min(ncs,4)
     best= Elect_dat(1)%k  ! Best overall cluster set

     If (Rating < 1) then
       Call Out ("For max slate size",mr, ln=1)
     Else
       Call Out ("For max positive slate size",mrp, "and negative",mrn, ln=1)
     End if

     Allocate(nl(ncs));  nl= Clust_set%l - 1
     Call Out ("Fractional cluster set objectives",Clust_set%fux)
     Call Out ("# regular clusters per cluster set",nl)


     Call Out (" ")
     Call Out ("Initial cluster set data for district "//District)

     Init_loop1 : Do iq= 1,N_init

       n0= Init(iq)%n;  q= Init(iq)%k
       If (n0 < 1 .or. q < 1) Cycle Init_loop1
       
       Call Out ("Initial cluster set",iq, "converged to cluster set",q,ln=1)
       n1= nl(q);  Call Out ("With initial # regular clusters",n0, "and final #",n1)
       Call Out ("via the mapping",Init(iq)%lt)
       
       Call Out ("Initial cluster weights", Init(iq)%T0(0,1:,1))
       Call Out ("Converged cluster weights", Clust_set(q)%sx)
       
       Call Out ("Initial optimal cluster objective & fitness data", Init(iq)%M1(:,3))
       Call Out ("Converged optimal cluster objective & fitness data", Clust_elect(q)%sx)
       
       If (n0 >= 2) Call Out (1,"Initial cluster correlations", Init(iq)%M2)
       If (n1 >= 2) Call Out (1,"Converged cluster correlations", Clust_set(q)%M1)

       Call Out ("Overall optimal electeds",Elect_dat(1)%ls)
       If (iq < 19) then
         Call Out (-1,"Initial electeds (Rem_frac,DHondt,optimal,converged)",Init(iq)%Q0(:,0,:))
       Else
         Call Out (-1,"Initial electeds (original,DHondt,init method,converged)",Init(iq)%Q0(:,0,:))
         Call Out ("Original electeds cluster objective & fitness data", Init(iq)%M1(:,1))
       End if
     End do Init_loop1
     
     
     Call Out (" ")
     Call Out ("STV cluster data for district "//District)

     nSTV= Non_clust%k;  ind= nSTV + 1;  n0= Init(19)%n;  in0= n0 + 1
     e= Non_clust%L1(1,1);  qo1= Elect_dat(e)%k;  qa= Init(19)%k

     Call Out ("STV elected set with final elected index",e, "has optimal converged cluster set",qo1, ln=1)
     Call Out ("vs converged set from extended STV",qa, "with original # clusters",nSTV)
     Call Out ("and initial set # clusters",n0)

     Call Out ("Mapping of original STV clusters to extended STV clusters",Non_clust%Q0(:,2,1))

     Call Out ("Extended STV cluster sizes", Non_clust%T0(0,1:ind,1))
     Call Out ("Initial STV cluster set sizes", Init(19)%T0(0,1:in0,1))

     wt1= Sum(Non_clust%T0(0,1:ind,1));  wt2= Sum(Init(19)%T0(0,1:in0,1))
     Call Out ("with extended weight sum",wt1, "and initial weight sum",wt2)

     Call Out ("Extended STV cluster averaged mean vector", Non_clust%T0(1:,0,1))
     Call Out ("Initial STV cluster averaged mean vector", Init(19)%T0(1:,0,1))
     Call Out ("Initial STV noise zeroed mean vector", Init(19)%T0(1:,0,2))

     If (nSTV > 1) Call Out (-1,"Extended STV correlations", Non_clust%T1(:nSTV,:nSTV,1))
     If (n0 > 1) Call Out (-1,"Initial STV cluster set correlations", Init(19)%M2)

     Call Out (-1,"Initial STV cluster set portions", Init(19)%T0(1:,1:in0,2))
     Call Out ("Initial STV cluster set objectives & fitness", Init(19)%M1(:,1))

     If (qa > 0) then
       Call Out ("Converged # clusters",nl(qa), ln=1)
       Call Out ("and elected set",Elect_dat(e)%ls)
         i= First_true(Elect_dat(e)%lt == qa)
       Call Out ("Converged STV clustering objectives & fitness", Elect_dat(e)%M1(:,i))
       If (n > 1) Call Out (-1,"Converged STV correlation matrix", Clust_set(qa)%M1)
     End if

     Call Out (" ")
     Call Out ("DTV cluster data for district "//District)

     nDTV= Non_clust%l;  ind= nDTV + 1;  n0= Init(20)%n;  in0= n0 + 1
     e= Non_clust%L1(2,1);  qo2= Elect_dat(e)%k;  qb= Init(20)%k

     Call Out ("DTV elected set with final elected index",e, "has optimal converged cluster set",qo2, ln=1)
     Call Out ("vs converged set from extended DTV",qb, "with original # clusters",nDTV)
     Call Out ("and initial set # clusters",n0)

     Call Out ("Mapping of original DTV clusters to extended DTV clusters",Non_clust%Q0(:,2,2))

     Call Out ("Extended DTV cluster sizes", Non_clust%T0(0,1:ind,2))
     Call Out ("Initial DTV cluster set sizes", Init(20)%T0(0,1:in0,1))

     wt1= Sum(Non_clust%T0(0,1:ind,2));  wt2= Sum(Init(20)%T0(0,1:in0,1))
     Call Out ("with extended weight sum",wt1, "and initial weight sum",wt2)

     Call Out ("Extended DTV cluster averaged mean vector", Non_clust%T0(1:,0,2))
     Call Out ("Initial DTV cluster averaged mean vector", Init(20)%T0(1:,0,1))
     Call Out ("Initial DTV noise zeroed mean vector", Init(20)%T0(1:,0,2))

     If (nDTV > 1) Call Out (-1,"Extended DTV correlations", Non_clust%T1(:nDTV,:nDTV,2))
     If (n0 > 1) Call Out (-1,"Initial DTV cluster set correlations", Init(20)%M2)

     Call Out (-1,"Initial DTV cluster set portions", Init(20)%T0(1:,1:in0,2))
     Call Out ("Initial DTV cluster set objectives & fitness", Init(20)%M1(:,1))

     If (qb > 0) then
       Call Out ("Converged # clusters",nl(qb), ln=1)
       Call Out ("and elected set",Elect_dat(e)%ls)
         i= First_true(Elect_dat(e)%lt == qb)
       Call Out ("Converged DTV clustering objectives & fitness", Elect_dat(e)%M1(:,i))
       If (n > 1) Call Out (-1,"Converged DTV correlation matrix", Clust_set(qb)%M1)
     End if


     Call Out (" ") 
     Call Out ("Elected set data for district "//District)

     Call Out ("Objective data for each elected set / cluster set combination:")
     Call Out ("full obj= clust obj * fitness,  coh vote * clust pen,  avg rate * elc pen")

     me2= Min(Max(Non_clust%L2(1,1), 10), me)

     Elc_loop1 : Do e= 1,me2
       Call Out ("Elected set",Elect_dat(e)%ls);  iq= Elect_dat(e)%k

       If (e <= elc_n(1)) then
         Call Out ("First tier elected set",e, "with best cluster set",iq, ln=1)
       Else if (e <= elc_n(2)) then
         Call Out ("Second tier elected set",e, "with best cluster set",iq, ln=1)
       Else
         Call Out ("General elected set",e, "with best cluster set",iq, ln=1)
       End if

       Call Out ("Objective ratio",Elect_dat(e)%fux, "overall satisfaction",Elect_dat(e)%fsx)
       Call Out ("Mean # electeds ranked by ballots, in order",Elect_dat(e)%M0(0,1), &
                 "or not", Elect_dat(e)%M0(0,2))

       Call Out ("All associated cluster sets (set, # clusters, objective data):")
       Do i= 1,Elect_dat(e)%m
         q= Elect_dat(e)%lt(i);  ncl= Clust_set(q)%l - 1
         Write (8,'(2I3, F8.3,2X, 2F7.3, 4X,2F7.3, 3X,2F7.3)') q, ncl, Elect_dat(e)%M1(:,i)
       End do
     End do Elc_loop1

     nd= dpr(0) 
     If (nd > 0) then
       Call Out ("Elected sets whose best full objective exceeds that of cluster set #1",dpr(1:nd))
     End if

    
     Call Out (" ")  
     Call Out ("Cluster set data, by clustering objective, for district "//District)

     obj_fac= fac_full / Clust_elect(best)%fsx
     n2= Last_true(nl > 1);  nq= Min(Max(3,n2), ncs)

     Set_loop : Do q= 1,nq

       ncl= nl(q);  ind= ncl + 1
       Call Out ("Cluster set",q,"with # regular clusters",ncl, ln=1)
       Call Out ("and cluster sizes",Clust_set(q)%sx)

       Call Out ("Top elected set",Clust_elect(q)%Q0(:,0,1))

       If (q == best) Call Out ("This is the best cluster set: max full objective")

       If (q == cls_n(1)) Call Out ("This is the last first tier cluster set")
       If (q == cls_n(2)) Call Out ("This is the last second tier cluster set")

       If (q == qo1) Call Out ("This is the optimal STV cluster set")
       If (q == qo2) Call Out ("This is the optimal DTV cluster set")

       If (q == qa) Call Out ("This is the actual STV cluster set")
       If (q == qb) Call Out ("This is the actual DTV cluster set")

       Call Out ("Initial cluster sets converging to this cluster set",Clust_set(q)%ls)

       Call Out ("Full objective and fitness data:")
       Write (8,'(F8.3,2X, 2F7.3, 4X,2F7.3, 3X,2F7.3)') Clust_elect(q)%sx
       
       Write(8,*);  cls_obj= Clust_set(q)%rx
       Write (8,'(2(A,F7.3))')  "Cluster set objective :",  cls_obj(1), &
                                "   Coherent vote      :",  cls_obj(2)
       Write (8,'(2(A,F7.3))')  "Cluster size penalty  :",  Clust_set(q)%qx(1), &
                                "   Correlation penalty:",  Clust_set(q)%qx(2)
       Write (8,'(2(A, F7.3))') "Independents penalty  :",  Clust_set(q)%qx(3), &
                                "   Residual penalty   :",  Clust_set(q)%qx(4)
       Write(8,*)
       Write (8,'(2(A, F7.3))') "Min cluster size   :",     Clust_set(q)%px(1),  &
                                   "   Max correlation :",  Clust_set(q)%px(2)
       Write (8,'(A,F7.3,A,F9.5)') "Indep fraction     :",  Clust_set(q)%px(3),  &
                                   "   Final residual  :",  Clust_set(q)%px(4)
       Write (8,'(2(A,F7.3))')   "Jacobian norm      :",    Clust_set(q)%px(5),  &
                                 "   Precond. norm   :",    Clust_set(q)%px(6)
  
       Write(8,*)
       Write (8,'(A, 6F7.3)') "Frac of voters represented by reg clusters", Clust_set(q)%ux
       Write (8,'(A, 6F7.3)') "          when that representation exceeds", Parm_rep
       Call Out ("Frac of voters represented fully by independents", 1.0 - Clust_set(q)%ux(1))

       If (ncl > 1) then
         Write(8,*)
         Call Out ("Cluster set correlation matrix:")
         Do cl= 1,ncl
           Write (8,'(A,I2, A, (5F7.3))') "  Cluster:", cl, &
                  "  with correlations:", Clust_set(q)%M1(:,cl)
         End do
       End if

!      Print out mean rating vector data 

       If (Allocated(set_ord)) DeAllocate(set_ord, orig_ord0, rates, por)
       
       Allocate(orig_ord0(nc), set_ord(nc), rates(nc,0:ind,2), por(nc,ind))
       
       Call Out ("The original candidates, as Borda ordered");  
       Write (8,"(20I3)") orig_cand

       set_ord= Clust_set(q)%L1(1:nc,0);  orig_ord0= orig_cand(set_ord)
       
       rates= Nint(D3*Clust_set(q)%T0(:,:,1:2)) / D3
       por  = Nint(D3*Clust_set(q)%T0(:,1:,3))  / D3
       
       Call Out ("Cluster-averaged zeroed mean vector")
         Write (8,"(10X,(10F7.3))") rates(:,0,2)
       Call Out ("using noise zeroing rating cutoff value",Noise_por)
         
       Call Out ("Downscaled by regular portion")
         Write (8,"(10X,(10F7.3))") rates(:,0,1)
         
       Call Out ("Corresponding candidate reordering")
         Write (8,"(10X,(10I7))") Clust_set(q)%L1(1:,0)
         
       Call Out ("The cluster sizes")
         Write (8,"(10X,(10F7.3))") Clust_set(q)%sx
       
       If (ncl > 0) then
         n= Min(10,nc);  Call Out ("Full cluster mean vectors:")

         Do cl= 1,ind
           Write (8,"(A,I2, A,10F7.3)")  "Cluster",cl, ":", rates(:n,cl,1)
         End do

         Call Out ("Zeored cluster mean vectors:")
         Do cl= 1,ind
           Write (8,"(A,I2, A,10F7.3)")  "Cluster",cl, ":", rates(:n,cl,2)
         End do

         Call Out ("Cluster portions for each candidate:")
         Do cl= 1,ind
           Write (8,"(A,I2, A,10F7.3)")  "Cluster",cl, ":", por(:n,cl)
         End do
       End if

       Call Out ("# cand top rated, signif, & viable, overall", Clust_set(q)%L1(-2:0,0))
       Call Out ("Viable candidates listed in rank order by cluster:")

       Do cl= 1,ind
         n= Min(Clust_set(q)%L1(0,cl), 10)
         Write (8,"(A,I2, A,2X, 10I3)")  "Cluster",cl,":", Clust_set(q)%L1(1:n,cl)
       End do
       
       Write (8,'(A)') " "
       Write (8,'(A, (10F7.3))') "Cluster sizes       : ", Clust_set(q)%sx
       Write (8,'(A, (10F7.3))') "Size deviations     : ", Clust_elect(q)%T0(:,1,1)
       Write (8,'(A, (10F7.3))') "Deviation penalties : ", Clust_elect(q)%T0(:,2,1)

       Write (8,'(A, (10F7.3))') "Cluster densities   : ", Clust_set(q)%M3(4,:)
       Write (8,'(A, (10F7.3))') "Density penalties   : ", Clust_set(q)%M3(2,:ncl)
       Write (8,'(A, (10F7.3))') "Cluster boundaries  : ", Clust_set(q)%M3(5,:)
       Write (8,'(A, (10F7.3))') "Boundary penalties  : ", Clust_set(q)%M3(3,:ncl)

       Write (8,'(A, (10F7.3))') "Frac fuzzy reduced  : ", Clust_set(q)%M3(7,:ncl)
       Write (8,'(A, (10F7.3))') "Average reduction   : ", Clust_set(q)%M3(8,:ncl)

       Call Out (-1,"DHondt electeds & clusters", Clust_elect(q)%L1)
       Call Out ("Indices of non-clustering elected sets",Non_clust%L2(Non_clust%vl,q))

       Call Out ("For district "//District//" and cluster set",q,ln=1)
       Call Out ("List data for its sets of electeds in order of fitness:")
       
       nf= Clust_elect(q)%m
       
       Do j= 1,nf
         elc_obj= Clust_elect(q)%M2(:,j);  full_obj=  cls_obj(1) * elc_obj(1)
         
         Write (8, "(I2, A,F8.3, A,3(F7.3), A,10I3)" )  j, "  Full obj:",full_obj, &
                "  fitness = avg rating * size penalty:", elc_obj, "  electeds:",  Clust_elect(q)%Q0(:,0,j)
         Write (8, "(9X, A, "//c2(ind)//"F7.3, A, "//c2(ind)//"F6.3)" ) "dev:", &
                Clust_elect(q)%T0(:,1,j), "    pen:", Clust_elect(q)%T0(:,2,j)
       End do
       
       Call Out ("Cluster set convergence data, excluding matching sets")
       Call Out ("initial & final # clusters, # convergence cycles", &
                  Clust_set(q)%vl(0:2))
       Call Out ("# merge & delete ops, final # centroidings, final # search points", &
                  Clust_set(q)%vl(3:6))
       Call Out ("Convergence history, by update cycle")
       If (Associated(Clust_set(q)%L0) .and. Associated(Clust_set(q)%M2)) then
         Call Out (-1,"   Restart, success code, # clusters, # centroidings, # search points", &
                          Clust_set(q)%L0)
         Call Out (-1,"   Damp parameter, functional, & net change", Clust_set(q)%M2)
       Else
         Call Out ("Error in Write_summary: No convergence history")
       End if
     End do Set_loop
    
    
     Call Out (" ")
     Call Out ("Comparative elected set / cluster set output for district "//District)

     Call Out ("# candidates",nc, "# to be elected",np)
     Call Out ("# cluster sets",ncs, "# possible elected sets",me)
     Call Out ("# regular clusters in each set", Clust_set%l-1)
     Call Out ("Full objective ordering for the top cluster sets")
     Write (8,*) ""
     
     If (np >= 3) then
       n3= np - 3;  p1= 9 + n3*4;  p2= 7 + n3*6
     Else
       p1= 3;  p2= 3
     End if

     Write (8, "(3X,A, "//c2(p1)//"X,A,  "//c2(p2)//"X,A, 4X,A)") &
            "Method", "Electeds-Party (order)", "Cluster sets", "& full objectives"
     
!    Best set of electeds for clustering (the match to the best cluster set)

     Allocate (elc_sets(np,me2,3))
     
     Do e= 1,me2
       elc_sets(:,e,1)= Elect_dat(e)%ls
       elc_sets(:,e,2)= orig_cand(elc_sets(:,e,1))
       elc_sets(:,e,3)= Elect_dat(e)%L0(:,1)
     End do
       
!    Write out the best clustering results

     e= 1
     Write (8,"(A, "//c2(np)//"(I3,A), A,I2, "//c2(np-1)//"I3, A)", Advance="NO")   &
            Method_lab(0),  (elc_sets(i,e,1), "-"//Party(elc_sets(i,e,2)), i=1,np), &
            "  (",(elc_sets(i,e,3),i=1,np),") "

     Covered= .false.;  e= 1;  Covered(e)= .true.;  ob= 0
     nq= Min(Elect_dat(e)%m,4)

     cs(:nq)= Elect_dat(e)%lt(:n);  cs(nq+1:)= -1
     ob(:nq)= Ceiling(obj_fac* Elect_dat(e)%M1(1,:nq))
     Write (8,"(1X, 4I3, 3X, 4I6)") cs, ob(:nq)
     
     Non_cls_sets= -1
     Do i= 1,nMt
       Non_cls_sets(:,i,1)= Non_clust%Q0(:,0,i)
       Non_cls_sets(:,i,2)= orig_cand(Non_cls_sets(:,i,1))
     End do
     
!    Write out the matching objective results for the non-clustering methods

     Method_loop : Do j= 1,nMt
       Write (8,"(A, "//c2(np)//"(I3,A), A,I2, "//c2(np-1)//"I3, A)", Advance="NO")          &
              Method_lab(j), (Non_cls_sets(i,j,1), "-"//Party(Non_cls_sets(i,j,2)), i=1,np), &
              "  (",(Non_clust%Q0(i,1,j),i=1,np),") "

       e= Non_clust%L1(j,1);  Covered(e)= .true.
       nq= Min(Elect_dat(e)%m,4)

       cs(:nq)= Elect_dat(e)%lt(:nq);  cs(nq+1:)= -1  
       ob(:nq)= Ceiling(obj_fac* Elect_dat(e)%M1(1,:nq))

       Write (8,"(1X, 4I3, 3X, 4I6)") cs, ob(:nq)
     End do Method_loop

     Call Out (" ")
       
     Elc_loop2 : Do e= 1,me2
       If (Covered(e)) Cycle Elc_loop2

       nq= Min(Elect_dat(e)%m,4)

       Write (8,"(A, "//c2(np)//"(I3,A), A,I2,"//c2(np-1)//"I3, A)", Advance="NO")   &
              Method_lab(0), (elc_sets(i,e,1), "-"//Party(elc_sets(i,e,2)), i=1,np), &
              "  (",(elc_sets(i,e,3),i=1,np),") "

       cs(:nq)= Elect_dat(e)%lt(:nq);  cs(nq+1:)= -1  
       ob(:nq)= Ceiling(obj_fac* Elect_dat(e)%M1(1,:nq))
       Write (8,"(1X, 4I3, 3X, 4I6)") cs, ob(:nq)
     End do Elc_loop2
     

     Call Out ("Best cluster set summary for district "//District)

     ind= Clust_set(best)%l
     Call Out ("Portion data for best cluster set",best, ln=1)
     Call Out ("")

     Do cl= 1,ind
       Write (8,'(A,I2,2(A,F6.3),A)', Iostat=ios) 'Clust #',cl,           &
                '  Cluster size:',Clust_set(best)%sx(cl),"  Deviation: ", &
                                  Clust_elect(best)%T0(cl,1,1),'  Portions: '
       Write (8,"(12F6.3)", Iostat=ios) Clust_set(best)%T0(1:,cl,3)
     End do


     Call Out ("Computational stats")

     n= nGen_clust - ngc(1);  i= nGen_clust_CM - ngc(2);  j= nGen_clust_NT - ngc(3)
     Call Out ("For district "//District//" # cluster set centroiding operations",n, ln=1)
     Call Out ("Count only operations via Cluster_from_means",i)
     Call Out ("Count only operations via Newton_update",j)

     Call Out ("Exit Write_summary")

   End Subroutine Write_summary
                             
End Module Clusters5 