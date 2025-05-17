
!        This module contains "Read_ballots" and its subroutines. 
    
!  "Read_ballots" is called from Clustering_PR to read the raw ballot files and 
!   pre-process the ballots, doing basic ballot consolidation and reordering.
    
 Module Clusters1

   Use Clusters0    
   Use Clusters_support

   Use Misc_methods
   Use Factorials
   Use Util
   Use Output
   Use Types
   Use Precisn
   Use IEEE_Arithmetic
   Implicit None

 Contains
    
   Subroutine Read_ballots (District,idist,nc_lim, nc0,nc,np,mr,mrp,nb, pt_val,orig0, &
                            global, wtb,ballot,ballot2, Non_clust)
   
!    Read the ballots, depending on the format. Either preprocess (consolidate) and store 
!    the preprocessed ballots or read previously preprocessed ballots.

     Character(*), Intent(in) :: District  ! Voting district Name
     Integer,  Intent(in) :: idist         ! Voting district index
     Integer , Intent(in) :: nc_lim        ! limit on # candidates to be used

     Integer, Intent(out) :: nc0           ! original # candidates
     Integer, Intent(out) :: nc            ! # candidates
     Integer, Intent(out) :: np            ! # candidates to be elected
     Integer, Intent(out) :: mr            ! max # candidates ranked
     Integer, Intent(out) :: mrp           ! max # postive candidates ranked
     Integer, Intent(out) :: nb            ! # ballots read

     Real,        Pointer :: pt_val(:)     ! (mt) Decreasing point values of the rating levels, must be pos or neg

     Integer,     Pointer :: orig0(:)      ! (nc) Mean vector ordering of the original 'nc0' candidates
     
     Real,        Pointer :: global(:,:)   ! (nc,5)  1 = global mean rating vector, 2 = its sigma vector
                                           ! For rankings which are not top ranked, or ratings which are not top
                                           !   or bottom rated, representing mostly 'noise':
                                           ! 3 = noise mean rating vector, 4 = corresponding sigma vector
                                           ! 5 = mean rating above the noise, normalized
     Real,        Pointer :: wtb(:)        ! (nb) Ballot weights summing to 'np'
     Integer,     Pointer :: ballot(:,:)   ! (0:mr,nb)  (1:nr,b) = candidates in preferential order
                                           !   (0,b) = 'nr' = # ranked or rated  
     Integer,     Pointer :: ballot2(:,:)  ! (0:mr,nb) (1:nr,b) =  increasing ranking or rating levels
                                           !   (0,b) = # positive ratings 
     Type(Multi_listR), Intent(inout) :: Non_clust  ! Data for 'nMt' selected non-clustering methods 'm'
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
                                     ! %L1(nMt,2)  (m,1)= Index 'e' in Elect_dat of non-clustering set 'm' 
                                     !             (m,2)= The matching best cluster set 'q'
                                     ! %L2(nonC,ncs) = Location, for each cluster set 'q', in Clust_elect(q)%Q0 
                                     !                of each distinct non-clustering set of electeds

                                     ! %T0(0:nc,0:ind,2) = Cluster data from STV (:,:,1) and DTV (:,:,2).
                                     !    (0,:,:) = cluster sizes, (1:,:,:) = cluster mean rating vectors
                                     !    (1:,0,:) = cluster averaged mean vector
                                     ! %T1(ncl,ncl,2) = Cluster correlations for STV (:,:,1) and DTV (:,:,2),
                                     !    = 'STV_corr' & 'DTV_corr' .
!  Local:
     Type(Set_list), Allocatable :: STVmb(:)   ! (0:nb) Ballot ranking structure for clustering
                                               ! Ballot sets (0): 
                                               !   %n        = nc
                                               !   %p        = np
                                               !   %svl      = elected vote = Sum(%val(set(:np)))
                                               !   %smb      = independents vote = Sum(%val(set(np+1:)))
                                               !   %set(nc)  = candidates in elected ordering, with
                                               !               the first 'np' elected
                                               !   %lev(nc)  = candidate clusters by decreasing cluster vote
                                               !   %val(nc)  = candidate cluster votes, not reordered

!                                              ! Ballots (1:nb):
                                               !   %n     = # elected candidates
                                               !   %p     = # ranked candidates
                                               !   %set(p)= 'n' elected candidates, then unelected,
                                               !            in original ranking order
                                               !   %val(p)= corresponding membership weights 
                                               !   %svl   = elected membership weight= Sum(%val(:n))
                                               !   %smb   = total membership weight  = Sum(%val)
     Type(Set_list), Allocatable :: DTVmb(:)   ! (0:nb) Ballot ranking structure for clustering
                                               ! Ballot sets (0): 
                                               !   %n        = nc
                                               !   %p        = np
                                               !   %svl      = elected vote = Sum(%val(set(:np)))
                                               !   %smb      = independents vote = Sum(%val(set(np+1:)))
                                               !   %set(nc)  = candidates in elected ordering, with
                                               !               the first 'np' elected
                                               !   %lev(nc)  = candidate clusters by decreasing cluster vote
                                               !   %val(nc)  = candidate cluster votes, not reordered

!                                              ! Ballots (1:nb):
                                               !   %n     = # elected candidates
                                               !   %p     = # ranked candidates
                                               !   %set(p)= 'n' elected candidates, then unelected,
                                               !            in original ranking order
                                               !   %val(p)= corresponding membership weights 
                                               !   %svl   = elected membership weight= Sum(%val(:n))
                                               !   %smb   = total membership weight  = Sum(%val)
     Type(Set_list), Allocatable :: STV_mrg(:) ! (np+1) Data for the deleted & merged STV clusters
                                               !        with cluster sizes %smb decreasing. 
                                               !  %n  = # ballot members of the cluster
                                               !  %smb= Cluster weight = Sum(STV_mrg%val)
                                               !  %set(n)= Ballot  members in increasing order
                                               !  %val(n)= Weighted memberships of the ballots
     Type(Set_list), Allocatable :: DTV_mrg(:) ! (n0+1=> ind2) Membership for the deleted & merged clusters
                                               !    with cluster sizes %smb decreasing. 
                                               !  %n  = # ballot members of the cluster
                                               !  %smb= Cluster weight = Sum(DTV_mrg%val)
                                               !  %set(n)= Ballot  members in increasing order
                                               !  %val(n)= Weighted memberships of the ballots
     Real,    Allocatable :: rank_pt(:,:)    ! (mr,2) Point values (:,1) and their variances (:,2) 
                                             !        for the top 'mr' candidates, assuming standard point values.
     Real,    Allocatable :: unrank_pt(:,:)  ! (mr,2) Point values and variances of unranked candidates,
                                             !        as a function of the # ranked
     Real,    Allocatable :: rt_var(:)       ! (0:mt) Variances for rating points, by rating level, 
                                             !        with (0) = unrated variance

     Real,    Allocatable :: Borda1(:,:)     ! (nc,np+1) Borda tally for STV candidate clusters 
     Real,    Allocatable :: Borda2(:,:)     ! (nc,np+1) Borda tally for STV candidate clusters 
                             
     Real,    Allocatable :: Corr1(:,:)      ! (np,np) Correlations for STV candidate clusters 
     Real,    Allocatable :: Corr2(:,:)      ! (np,np) Correlations for DTV candidate clusters 

     Type(Set_list), Allocatable :: Clust1(:)  ! (ind1) Data for deleted & merged STV clusters
                                               !    with cluster sizes %smb decreasing. 
                                               !  %n  = # ballot members of the cluster
                                               !  %smb= Cluster weight = Sum(Clust1%val)
                                               !  %set(n)= Ballot  members in increasing order
                                               !  %val(n)= Weighted memberships of the ballots
     Type(Set_list), Allocatable :: Clust2(:)  ! (ind2) Data for deleted & merged DTV clusters
                                               !    with cluster sizes %smb decreasing. 
                                               !  %n  = # ballot members of the cluster
                                               !  %smb= Cluster weight = Sum(Clust2%val)
                                               !  %set(n)= Ballot  members in increasing order
                                               !  %val(n)= Weighted memberships of the ballots
     Integer, Allocatable :: key1(:), key2(:), cls1(:), cls2(:), electeds(:,:,:)
     Real,    Allocatable :: ptval(:)

     Character(80) :: msg
     Character(18) :: Dist     ! District of the election (inside the region)
     Character(16) :: Region   ! Region of the election
     Character(65) :: label
     Integer :: Year           ! Year of the election
     Integer :: UVP0           ! Under vote penalty if already used for the election file
     Real    :: x1, x2, fac
     Integer :: ind, ind1, ind2, ios, nb0, nb1, ncl, ndel, nmrg, nonC, nSTV, nDTV, tot_cnt
     Integer :: b, c, i, j, m, n, rt, mt, n1, n2, mr0, mrn, nbf, ntp
     
     If (Associated(orig0))     DeAllocate(orig0)
     If (Associated(global))    DeAllocate(global)

     If (Associated(wtb))       DeAllocate(wtb)
     If (Associated(ballot))    DeAllocate(ballot)
     If (Associated(ballot2))   DeAllocate(ballot2)

     If (Associated(Non_clust%L0)) DeAllocate(Non_clust%L0)
     If (Associated(Non_clust%L2)) DeAllocate(Non_clust%L2)
     If (Associated(Non_clust%Q0)) DeAllocate(Non_clust%Q0)
     If (Associated(Non_clust%T0)) DeAllocate(Non_clust%T0)
     If (Associated(Non_clust%T1)) DeAllocate(Non_clust%T1)

     If (Associated(Req_cand))  DeAllocate(Req_cand)
     If (Associated(Pos_cand))  DeAllocate(Pos_cand)
     
     Call Out ("Read ballots and simple consolidation for district #",idist,ln=1)
     
!    First read the ballot parameters

     If (Rating >= 1) then
       Open(9, File=Indat//District//'3', IOmsg=msg, IOstat=ios, Status='Old', Action='Read')
     Else if (Rating == 0) then
       Open(9, File=Indat//District//'2', IOmsg=msg, IOstat=ios, Status='Old', Action='Read')
     Else
       Open(9, File=Indat//District//'1', IOmsg=msg, IOstat=ios, Status='Old', Action='Read')
     End if

!             3 typical labels to be read:
!   'Rating type, # personal ballots, # processed, final # combined: '
!   '# candidates to be elected, # reduced, original #: '
!   'Max # ranked or rated, max # rating levels, undervote penalty: '

     Read(9,*) Dist, Region, Year
     Read(9,*) Label, rt, tot_cnt, nbf, nb
     Read(9,*) Label, np, nc, nc0
     Read(9,*) Label, mr0, mt, UVP0

     nc= Min(nc, nc_lim);  mr= Min(mr_spec, mr0, nc)
     mt= Max(mt,1)  ! mt = 1 for Rating = -1

     If (rt /= Rating) then
       Call Out("Error in 'Read_ballots'. Rating mismatch. Assumed value", Rating, &
                "read from file",rt, ln=1)
       Stop
     Else if (np < 2 .or. nc > nc0 .or. nc <= np .or. nb < nc .or. mr < 4) then
       Call Out ("Input parameter error in 'Read_ballots': Reduced # candidates",nc, &
                 "vs orig # cand",nc0, ln=1)
       Call Out ("# candidates to be elected",np, "Max # candidates to be ranked or rated",mr)
       Call Out ("# combined ballots",nb)
       Stop
     End if
       
!    Typcial label: 'Original candidates, as reordered'

     Allocate(orig0(nc));  Read(9,*) Label, orig0

!    'pt_val' input  point conversion values for rating or sometimes for weak ranking

     If (Rating == 1) then    ! Rating data
       Read(9,*);  Read(9,*) label
       Read(9,*) label, mrp, mrn

       n= mrp + mrn 
       If (n /= mr) then
         Call Out ("Warning in 'Read_ballots': mrp+mrn /= prior mr")
         Call Out ("Prior max # rated",mr, "vs max # pos. rated",mrp)  
         Call Out ("& max # neg. rated",mrn)
         If (mr > n) then
           mr= n;  Call Out ("Prior max # rated reduced to",mr)
         Else
           fac= mr / Real(n)
           mrp= Nint(fac * mrp);  mrn= mr - mrp
           Call Out ("Revised vs max # pos. rated",mrp, "and neg.",mrn)
         End if
       End if
       
       Allocate(pt_val(mt), rt_var(0:mt), ptval(mt), rank_pt(1,1), unrank_pt(1,1))
       Read(9,*) label, pt_val

       ptval= pt_val;  rank_pt= -1; unrank_pt= -1
       Call Rating0_pts (Parm_bal1,mt,ptval, rt_var)

     Else if (Rating == 0) then  ! Weak rankig
       mrp= mr;  mrn= 0

       If (mt > 3) then   ! Read in weak ranking points
         Allocate(pt_val(mt), rt_var(0:mt), ptval(mt), rank_pt(1,1), unrank_pt(1,1))
         Read(9,*) label, pt_val

         ptval= pt_val;  rank_pt= -1; unrank_pt= -1  
         Call Rating0_pts (Parm_bal1,mt,ptval, rt_var)
       Else               ! Use standard rank to point conversion 
         Allocate(rank_pt(mr,2), unrank_pt(mr,2), pt_val(1),  ptval(mr), rt_var(1))
         Call Ranking_pts (Parm_bal1,Max_pt,nc,mr, rank_pt,unrank_pt)
         ptval= rank_pt(:,1);  pt_val= -1;  rt_var= -1
       End if

     Else    ! Strong ranking
       mrp= mr;  mrn= 0;  mt= 1
       Allocate(rank_pt(mr,2), unrank_pt(mr,2), ptval(mr), pt_val(1), rt_var(1))
       Call Ranking_pts (Parm_bal1,Max_pt,nc,mr, rank_pt,unrank_pt)
       ptval= rank_pt(:,1);  pt_val= -1;  rt_var= -1
     End if
       
!    Memory allocation

     Allocate(wtb(nb), ballot(0:mr,nb), global(nc,5))

     If (Rating >= 0) then 
       Allocate(ballot2(0:mr,nb))
     Else
       Allocate(ballot2(0:1,1));  ballot2= -1
     End if

!    Now read the pre-processed ballots
       
     nb0= nb
     Call Read_ratings (Rating,np,mrp,mr,nc,nb0, rank_pt,unrank_pt, &
                        pt_val,rt_var,nb, wtb,ballot,ballot2, global)
     Close (9)

!    Noise levels for subtracting from mean rating vectors
!    for correlations (Noise_cor) cluster portions (Noise_por)

     If (Rating < 1) then
       x1= 0.85 * Dot_product(global(:,1), global(:,5))
       x2= Dot_product(global(:,3) + Noise_fac * global(:,4), global(:,5))

       Noise_cor= Min(x1,x2)
       Noise_cor= Max(Min(Noise_cor, Noise_max), Noise_min)
       Noise_por= Por_fac * Noise_cor

       Call Out ("1st estimate of noise to subtract for correlations", x1, &
                 "2nd estimate",x2)
     Else
       x2= Dot_product(global(:,4), global(:,5))
       Noise_cor= 0.0;  Noise_por= Por_fac * x2
     End if

     Call Out ("For district #", idist, ln=1)
     Call Out ("Noise to subtract for correlations", Noise_cor, &
                "noise cutoff for portions",Noise_por)

     Allocate(Pos_cand(nc));  Pos_cand(:np)= .true.;  Pos_cand(np+1:)= .false.

!    Next apply some non-clustering algorithms to these ballots.

     Allocate(Non_clust%Q0(np,0:2,nMt), STVmb(0:nb), DTVmb(0:nb), electeds(nc,nMt,2))

     Method= .true.;  Non_clust%Q0= 0
  
     Call Non_clustering (idist,Method,Rating, np,mr,nc,nb, pt_val, &
                          wtb(:nb),ballot(:,:nb),ballot2(:,:nb),    &
                          electeds,STVmb,DTVmb)
     
     Non_clust%Q0(:,0,:)= electeds(:np,:,2)
     Non_clust%Q0(:,1,:)= electeds(:np,:,1)

     Call Out (-1,"Non-clustering electeds", Non_clust%Q0(:np,:,2))

     Allocate(Req_cand(nc))
     Call Non_clustering_data (nc,np,nMt,Method, Non_clust%Q0(:,:1,:), &
                               nonC,Non_clust%L0,Non_clust%vl,Req_cand)

     Non_clust%n= nonC;  If (Rating == 1) Pos_cand= Pos_cand .or. Req_cand
     
     ind= np + 1
     Allocate (STV_mrg(ind), DTV_mrg(ind), Clust1(ind), Clust2(ind),       &
               Corr1(np,np), Corr2(np,np), Borda1(nc,ind), Borda2(nc,ind), &
               key1(np), key2(np), cls1(np), cls2(np))
        
     Call STV_clusters (np,nb,STVmb, STV_mrg,key1)

     If (Rating == 1) then
       Call Merge_STV2 (np,nc, ptval,STVmd_parm, ballot(:,:nb),ballot2(:,:nb),  &
                        STV_mrg, nSTV,cls1, Clust1,Borda1,Corr1)      
     Else if (Rating == 0) then
       Call Merge_STV2 (np,nc, ptval,STVmd_parm, ballot(:,:nb),ballot2(:,:nb), &
                        STV_mrg, nSTV,cls1, Clust1,Borda1,Corr1, Noise_cor)      
     Else
       Call Merge_STV (np,nc, Noise_cor,ptval, STVmd_parm,ballot(:,:nb), &
                       STV_mrg, nSTV,cls1, Clust1,Borda1,Corr1)
     End if
     Non_clust%Q0(:,2,1)= cls1(key1)

     Call STV_clusters (np,nb,DTVmb, DTV_mrg,key2)

     If (Rating == 1) then
       Call Merge_STV2 (np,nc, ptval,STVmd_parm, ballot(:,:nb),ballot2(:,:nb), &
                        DTV_mrg, nDTV,cls2, Clust2,Borda2,Corr2)      
     Else if (Rating == 0) then
       Call Merge_STV2 (np,nc, ptval,STVmd_parm, ballot(:,:nb),ballot2(:,:nb), &
                        DTV_mrg, nDTV,cls2, Clust2,Borda2,Corr2, Noise_cor)      
     Else
       Call Merge_STV (np,nc, Noise_cor,ptval, STVmd_parm,ballot(:,:nb), &
                       DTV_mrg, nDTV,cls2, Clust2,Borda2,Corr2)
     End if
     Non_clust%Q0(:,2,2)= cls2(key2)

     Non_clust%k= nSTV;  Non_clust%l= nDTV
     ind1= nSTV + 1;  ind2= nDTV + 1;
     ncl= Max(nSTV,nDTV);  ind= ncl + 1

     Allocate(Non_clust%T0(0:nc,0:ind,2), Non_clust%T1(ncl,ncl,2))
     Non_clust%T0= 0;  Non_clust%T1= 0

     Non_clust%T0(0,1:ind1,1)= Clust1(:ind1)%smb;  Non_clust%T0(1:,1:ind1,1)= Borda1(:,:ind1)
       Non_clust%T1(:nSTV,:nSTV,1)= Corr1(:nSTV,:nSTV)
     Non_clust%T0(0,1:ind2,2)= Clust2(:ind2)%smb;  Non_clust%T0(1:,1:ind2,2)= Borda2(:,:ind2)
       Non_clust%T1(:nDTV,:nDTV,2)= Corr2(:nDTV,:nDTV)

     Forall(c=1:nc) 
       Non_clust%T0(c,0,1)= Sum(Non_clust%T0(0,1:ind1,1) * Non_clust%T0(c,1:ind1,1)) / np
       Non_clust%T0(c,0,2)= Sum(Non_clust%T0(0,1:ind2,2) * Non_clust%T0(c,1:ind2,2)) / np
     End forall

     If (pr_out >= 1) then
       Call Out ("STV cluster data for # elected candidates",np, ln=1)
       Call Out ("Final # regular clusters",nSTV)
       Call Out ("Original candidate cluster order",STVmb(0)%set)
       Call Out ("Matching sizes",STVmb(0)%val)
       Call Out ("Original candidate cluster assignment to final",Non_clust%Q0(:,2,1))
       Call Out ("Final cluster sizes",Non_clust%T0(0,1:ind1,1))
       Call Out ("Final cluster-averaged mean vector",Non_clust%T0(1:,0,1))
       Call Out (-1,"Final STV cluster mean vectors",Borda1(:,:ind1))
       Call Out (-1,"Final STV cluster correlatrions",Corr1(:nSTV,:nSTV))

       Call Out ("DTV cluster data for # elected candidates",np, ln=1)
       Call Out ("Final # regular clusters",nDTV)
       Call Out ("Original candidate cluster order",DTVmb(0)%set)
       Call Out ("Matching sizes",DTVmb(0)%val)
       Call Out ("Original candidate cluster assignment to final",Non_clust%Q0(:,2,2))
       Call Out ("Final cluster sizes",Non_clust%T0(0,1:ind2,2))
       Call Out ("Final cluster-averaged mean vector",Non_clust%T0(1:,0,2))
       Call Out (-1,"Final DTV cluster mean vectors",Borda2(:,:ind2))
       Call Out (-1,"Final DTV cluster correlatrions",Corr2(:nDTV,:nDTV))
     End if

     Call Out ("Parameters for district #", idist)
     Call Out ("# original candidates",nc0, "# reduced",nc)
     Call Out ("# candidates to be elected",np, "max # ranked or rated",mr)
     Call Out ("# initial ballots",nb)

     If (Rating == 1) Call Out ("# pos rated candidates",mrp, &
                                "# neg rated candidates",mrn)
         
     Call Out ("Original candidates, as read in",orig0)
     Call Out (-1,"1:2(global mean vec & sig) 3:4(noise vec & sig) 5(mean - noise,normed)",global)

     If (Rating >= 0 .and. mt >= 4) then 
       Call Out ("Rating point conversion factors",pt_val)
     End if  
   End Subroutine Read_ballots

            
  Subroutine Read_ratings (Rating,np,mrp,mr, nc,nb0, rank_pt,unrank_pt, &
                           pt_val,rt_var, nb,wtb, ballot,ballot2, global)
  
!   Read the ranked or rated candidates, plus their rankings or ratings, 
!   from a ballot file pre-processed to be in the single line format 
!   used by the subroutine 'Read_vote'
  
    Integer,  Intent(in) :: Rating          ! = 1 for ratings 
                                            ! = 0 for weak ranking 
                                            ! = -1 for strong ranking
    Integer,  Intent(in) :: np              ! # candidates to be elected = total weight
    Integer,  Intent(in) :: mrp             ! Max # positively rated candidates per ballot
    Integer,  Intent(in) :: mr              ! Max # candidates to be ranked or rated
    Integer,  Intent(in) :: nc              ! Last candidate to be used
    Integer,  Intent(in) :: nb0             ! # ballots to be read
    Real,     Intent(in) :: rank_pt(:,:)    ! (mr,2) Point values (:,1) 
                                            !   and their variances (:,2) 
                                            !   for the top 'mr' candidates,
                                            !   assuming standard points
    Real,     Intent(in) :: unrank_pt(:,:)  ! (mr,2) Point values and
                                            !   variances of unranked cand
    Real,     Intent(in) :: pt_val(:)       ! (mt) Decreasing points for rating levels
    Real,     Intent(in) :: rt_var(0:)      ! (0:mt) Variances for rating 
                                            !   points, by rating level,
                                            !   with (0) = unrated variance
    
    Integer, Intent(out) :: nb              ! Final # ballots after any
                                            !   consolidation, if mr < mr0
    Real,    Intent(out) :: wtb(:)          ! (nb) Ballot weights summing to 'np'
    Integer, Intent(out) :: ballot(0:,:)    ! (0:mr,nb)  (1:nr,b) = candidates in preferential order
                                            !   (0,b) = 'nr' = # ranked or rated  
    Integer, Intent(out) :: ballot2(0:,:)   ! (0:mr,nb) (1:nr,b) =  increasing ranking or rating levels
                                            !   (0,b) = # positive ratings 
    Real,    Intent(out) :: global(:,:)     ! (nc,5)  1 = global mean rating vector, 2 = its sigma vector
                                            ! For rankings which are not top ranked, or ratings which are not top
                                            !   or bottom rated, representing mostly 'noise':
                                            ! 3 = noise mean rating vector, 4 = corresponding sigma vector
                                            ! 5 = mean rating above the noise, normalized
! Local:
    Real, Parameter :: eps= 0.000001

    Real    :: mean(nc), tmp(nc), smg(nc), smz(nc)
    Real    :: var_pt(nc), var_pt2(nc), wt_var(nc), wt_varz(nc)
    Integer :: bal(0:mr), bal2(0:mr), key(nc), Nrk(mr)

    Real    :: wt
    Integer :: b, i, l, n, ib, lp, mt, n1, nr, nz, nw, mrn, nzp, nzm, orig_nb

    Call Out ("Enter 'Read_ratings'")

!   Read in and consolidate any duplicate ballots, from a new 'nc' limit for example
 
    mt= Size(pt_val);  mrn= mr - mrp
    ballot= 0;  ballot2= 0;  nb= 0;  wtb= 0

    Ballot_loop1 : Do b= 1,nb0

      Call Read_vote (Rating,b,mrp,nc, nw, bal,bal2)

      nr= bal(0);  lp= bal2(0)
      If (nr < 1) Cycle Ballot_loop1
    
!     Combine initial ballots which are identical
      
      If (Rating >= 0) then
        Do ib= 1,nb
          n= ballot(0,ib);  l= ballot2(0,ib)
          If (n /= nr .or. l /= lp) Cycle

          If (All(bal(1:n)  == ballot(1:n,ib)) .and. &
              All(bal2(1:l) == ballot2(1:l,ib))) then
            wtb(ib)= wtb(ib) + nw;  Cycle Ballot_loop1
          End if
        End do

!       New ballot

        nb= nb + 1;  wtb(nb)= nw
        ballot(:nr,nb) = bal(:nr)
        ballot2(:nr,nb)= bal2(:nr)
      Else
        Do ib= 1,nb
          n= ballot(0,ib);  If (n /= nr) Cycle

          If (All(bal(1:n) == ballot(1:n,ib))) then
            wtb(ib)= wtb(ib) + nw;  Cycle Ballot_loop1
          End if
        End do

!       New ballot

        nb= nb + 1;  wtb(nb)= nw
        ballot(:nr,nb)= bal(:nr)
      End if

    End do Ballot_loop1

    orig_nb= Sum(wtb(:nb));  wt= Real(np);  wtb(:nb)= (wt/orig_nb) * wtb(:nb)
    
    Call Out ("Consolidated # ballots",nb, "vs prior #",nb0, ln=1)
    Call Out ("Original total # ballots", orig_nb, "normalized to weight",wt)

!   Summations for the 'global' data

    smz= 0;  smg= 0;  global= 0

    If (Rating < 1) then
      nz= Nint(0.5*mr)
    Else
      nzp= Nint(0.5*mrp);  nzm= Nint(0.5*mrn) 
    End if

    Ballot_loop2 : Do b= 1,nb
      nr= ballot(0,b)

      If (Rating == -1) then
        mean= unrank_pt(nr,1);  var_pt= unrank_pt(nr,2)
        mean(ballot(1:nr,b))  = rank_pt(:nr,1)
        var_pt(ballot(1:nr,b))= unrank_pt(:nr,2)

      Else if (Rating == 0 .and. mt < 4) then
        mean= unrank_pt(nr,1);  var_pt= unrank_pt(nr,2)
        mean(ballot(1:nr,b))  = rank_pt(ballot2(1:nr,b),1)
        var_pt(ballot(1:nr,b))= rank_pt(ballot2(1:nr,b),2)

      Else
        mean= 0;  var_pt= rt_var(0)
        mean(ballot(1:nr,b))  = pt_val(ballot2(1:nr,b))
        var_pt(ballot(1:nr,b))= rt_var(ballot2(1:nr,b))
      End if

      wt_var= wtb(b) / var_pt;  var_pt2= var_pt + mean**2
      smg= smg + wt_var

      global(:,1)= global(:,1) + mean * wt_var
      global(:,2)= global(:,2) + var_pt2 * wt_var
      
!     Focus on the noise by zeroing top rankings, also bottom for Rating = 1

      tmp= mean;  Call Sort (.false.,tmp,key)
      wt_varz= wt_var

      If (Rating <= 0) then
        If (nr > nz) then
          wt_varz(key(:nz))= 0;  smz= smz + wt_varz

          global(:,3)= global(:,3) + mean * wt_varz
          global(:,4)= global(:,4) + var_pt2 * wt_varz
        End if
      Else
        lp= ballot2(0,b);  n= Last_true(tmp > 0);  n1= n + 1

        If (n > nzp .and. n1 < nr - nzm) then
          wt_varz(key(:nzp))= 0;  wt_varz(key(nr-nzm:))= 0
          smz= smz + wt_varz

          global(:,3)= global(:,3) + Abs(mean) * wt_varz
          global(:,4)= global(:,4) + var_pt2 * wt_varz
        End if
      End if
    End do Ballot_loop2

!   Final 'global' computations

    smg= Max(smg,eps);  global(:,1)= global(:,1) / smg
    global(:,2)= Sqrt(global(:,2) / smg - global(:,1)**2)
    
    smz= Max(smz,eps);  global(:,3)= global(:,3) / smz
    global(:,4)= Sqrt(global(:,4) / smz - global(:,3)**2)

    tmp= global(:,1) - global(:,3);  global(:,5)= tmp / Sum(tmp)
    
    Do n= 1,mr
      Nrk(n)= Count(ballot(0,:) == n)
    End do

    Call Out ("# initial ballots for each # ranked",Nrk(:mr))
    Call Out (-1,"Global: mean & sig, noise & sig, mean - noise (normed)",global)
    
  End Subroutine Read_ratings
   
   
  Subroutine Read_vote (Rating,b, mrp,nc, nw, ballot,ballot2)
  
!   Read a single line 'bal' of the ballot file, assuming the following 
!   format:  bal(-1) = 'nw' = ballot weight = # original ballots represented, 
!   bal(0) = 'nr' = # candidates ranked or rated.

!   Then if Rating = -1 (strong ranking), bal(i) is the ith ranked candidate 'cn'
!   starting with the top ranked candidate, proceeding through the remaining positions
!   with '0': nw, nr, c1,...,ci,0...0 with i = nr, followed by mr - nr 0's.

!   If Rating >= -1, this is followed on the same line by bal(mr+1)= lp and 
!   bal(mr+1+i) = the rankng or rating level of candidate bal(i), proceeding
!   through the remaining mr - nr positions with 0. 
  
    Integer,  Intent(in) :: Rating       ! Rating data(1), weak ranking(0), strong ranking(-1)
    Integer,  Intent(in) :: b            ! Index of current ballot
    Integer,  Intent(in) :: mrp          ! Max # positive ratings permitted
    Integer,  Intent(in) :: nc           ! Last candidate to be used (for Rating <= 0)
    
    Integer, Intent(out) :: nw           ! ballot weight = # original ballots represented
    Integer, Intent(out) :: ballot(0:)   ! (0:mr)  (1:nr) = candidates in preferential order
                                         !   (0) = # ranked or rated  
    Integer, Intent(out) :: ballot2(0:)  ! (0:mr)  increasing ranking or rating levels
                                         !   with (0) = # positive levels, if Rating >= 0
! Local:
    Integer :: n, ln, lp, mr, nr, ios, mrn, lcn(nc)
    
    mr= Ubound(ballot,1);  ballot= 0;  ballot2= 0

    If (Rating >= 0) then
      Read(9,'(I7,2X,18I3)', iostat= ios) nw, ballot, ballot2
      nr= ballot(0);  lp= ballot2(0)
    Else
      Read(9,'(I7,2X,10I3)', iostat= ios) nw, ballot
      nr= ballot(0);  lp= nr
    End if

    If (ios /= 0 .or. nr <= 0 .or. lp <= 0 .or. nr > mr) then
      Call Out ("Error in 'Read_vote': at ballot",b, ln=1)
      Call Out ("ballot data",ballot)
      If (Rating == 1) Call Out ("ballot2 data",ballot2)
      Stop
    End if
    
    If (Rating == 1) then
      ln= nr - lp;  mrn= mr - mrp
      If (lp > mrp .or. ln > mrn .or. nr < 2) then
        Call Out ("Error in 'Read_vote': # pos ratings",lp, "# neg ratings",ln, ln=1)
        Call Out ("vs max # pos",mrp, "max # neg",mrn)
        Call Out ("ballot data",ballot)
        Call Out ("ballot2 data",ballot2)
        Stop
      End if
    Else
      Call List_of_true (ballot(1:nr) <= nc, n,lcn)

      If (n < nr) then
        ballot(1:n)= ballot(lcn(:n))
        ballot(0)= n; ballot(n+1:)= 0
        ballot2(1:n)= ballot2(lcn(:n))
        ballot2(0)= n; ballot2(n+1:)= 0
      End if
    End if
  End Subroutine Read_vote

  Subroutine Read_vote0 (b,nc, nw,ballot)
  
!   Read a single line 'bal' of the ballot file, assuming the following 
!   format:  bal(-1) = 'nw' = ballot weight = # original ballots represented, 
!   bal(0) = 'nr' = # candidates ranked or rated.
  
    Integer,  Intent(in) :: b          ! Ballot index being read
    Integer,  Intent(in) :: nc         ! Last candidate to be used

    Integer, Intent(out) :: nw         ! ballot weight = # original ballots represented
    Integer, Intent(out) :: ballot(0:) ! (0:mr)  (1:nr) = candidates in rank order
                                       !   (0) = 'nr' = # ranked 
! Local:
    Integer, Parameter :: mxc= 8
    Integer :: lcn(mxc)
    Integer :: n, mr, nr, ios
    
    Read(9,'(I7,2X,9I3)', iostat= ios) nw, ballot
    nr= ballot(0);  mr= Ubound(ballot,1)

    If (ios /= 0 .or. nw < 1 .or. nr < 1 .or. nr > mr .or. ballot(1) < 1) then
      Call Out ("Error in 'Read_vote0': at ballot",b, "with weight",nw, ln=1)
      Call Out ("and ballot data",ballot)
      Stop
    End if

    Call List_of_true (ballot(1:nr) <= nc, n,lcn)

    If (n < nr) then
      ballot(1:n)= ballot(lcn(:n)) 
      ballot(0)= n; ballot(n+1:)= 0
    End if

  End Subroutine Read_vote0


   Subroutine Non_clustering_data (nc,np,nMt,Method,Non_elc, nonC,Distinct,Non_map, Req_cand)
   
!    Determine the 'Distinct' sets of electeds candidates as computed by the 
!    non-clustering methods. Also determine any required candidates Req_cand - 
!    from the top ranked non-clustering methods. 
   
     Integer,  Intent(in) :: nc               ! # candidates
     Integer,  Intent(in) :: np               ! # candidates to be elected
     Integer,  Intent(in) :: nMt              ! Total # possible nonclustering methods

     Logical,  Intent(in) :: Method(:)        ! (nMt)    True for all methods being tested
     Integer,  Intent(in) :: Non_elc(:,0:,:)  ! (np,0:1,nMt) Elected candidates in increasing order (0)
                                              !   or preferential order (1), for the 'nMt' 
                                              !   non-clustering methods, with 1 = STV & 2 = DTV

     Integer, Intent(out) :: nonC             ! Total # possible nonclustering methods
     Integer,     Pointer :: Distinct(:,:)    ! (np,nonC) Distinct non-clustering sets of electeds

     Integer, Intent(out) :: Non_map(nMt)     ! For each method 'j', Non_map(j) = its index 'i' in Distinct

     Logical, Intent(out) :: Req_cand(:)      ! (nc)  True for the top candidate most often elected
                                              !       by the nonclustering methods
!  Local:
     Integer :: elc(np,nMt), count(nc)
     Real    :: xt
     Integer :: i, j, top
     
     Non_map= 0;  Req_cand= .false.
     count= 0;  elc= -1;  nonC= 0
     
     Do j= 1,nMt
       If (.not.Method(j)) Cycle
       
       If (nonC > 0) then
         i= Locate_set(Non_elc(:,0,j), elc(:,:nonC))
           
         If (i < 1) then
           nonC= nonC + 1;  Non_map(j)= nonC;  elc(:,nonC)= Non_elc(:,0,j)
         Else  
           Non_map(j)= i
         End if
       Else
         nonC= 1;  Non_map(j)= nonC;  elc(:,nonC)= Non_elc(:,0,j)
       End if
       
       top= Non_elc(1,1,j);  count(top)= count(top) + 1
     End do
     
     Allocate(Distinct(np,nonC));  Distinct= elc(:,:nonC)

     i= Maxloc(count,1);  Req_cand(i)= .true.

     If (pr_out > 1.5) then
       Call Out ("Counts for the top rated candidate of each Method",count)
       Call Out (-1,"The elected candidates for the non-clustering methods",Non_elc(:,0,:))
       Call Out (-1,"These candidates in preferential order",Non_elc(:,1,:))
     
       Call Out (-1,"The corresponding distinct sets of electeds in increasing order",Distinct)
       Call Out ("From the method to distinct set mapping",Non_map)
     End if
     
   End Subroutine Non_clustering_data


  Subroutine Domain_ballots (np,mr,nc,nco,domain, nb0,wtb0,ballot0, &
                             nb,tot_wt,wtb,ballot)

!   Read in all ballots but record only those which rank at least one candidate
!   in the domain 'domain', eliminating rankings of all candidates outside the domain.
  
    Integer,  Intent(in) :: np              ! # candidates to be elected = weight sum
    Integer,  Intent(in) :: mr              ! Max # candidates ranked per output ballot
                                            !   mr <= nc required
    Integer,  Intent(in) :: nc              ! Reduced # candidates
    Integer,  Intent(in) :: nco             ! Original # candidates
    Integer,  Intent(in) :: domain(:)       ! (nc) Subset of candidates to be used

    Integer,  Intent(in) :: nb0             ! original # ballots
    Real,     Intent(in) :: wtb0(:)         ! (nb0) Ballot weights summing to 'np'
    Integer,  Intent(in) :: ballot0(0:,:)   ! (0:,nb0)  (1:n,b) = candidates in preferential order
                                            !   (0,b) = 'n' = # ranked, possible n > mr 

    Integer, Intent(out) :: nb              ! Final # ballots after consolidation
    Real,    Intent(out) :: tot_wt          ! Total weight reduction fraction
    Real,    Intent(out) :: wtb(:)          ! (nb) Final ballot weights normalized to sum to 'np'
    Integer, Intent(out) :: ballot(0:,:)    ! (0:mr,nb)  (1:n,b) = candidates in preferential order
! Local:
    Real    :: rank_pt(mr), unrank_pt(mr)
    Integer :: bal0(mr), bal(0:mr), inv(nco)

    Logical :: ReOrd
    Real    :: wt
    Integer :: b, i, l, m, n, ib

    If (pr_out >= 1) then
      Call Out ("Enter 'Domain_ballots', reducing # candidates from",nco, "to",nc, ln=1)
      Call Out ("with max # ranked",mr, "for # open positions",np)
    End if

!   Combine prior ballots which are identical after removal of 
!   outside-the-domain candidates
 
    nb= 0;  wtb= 0;  ballot= 0
    Call Inverse_map(domain,inv)

    Call Ranking_pt0 (Parm_bal1,Max_pt,nc,mr, rank_pt,unrank_pt)

    Ballot_loop1 : Do b= 1,nb0
      n= Min(ballot0(0,b),mr);  bal0(:n)= inv(ballot0(1:n,b))
      Call List_of_true (bal0(:n) > 0, m,bal(1:))  ! bal(1:m)= cand inside the domain
      If (m < 1) Cycle Ballot_loop1

      bal(0)= m;  bal(1:m)= bal0(bal(1:m));  wt= wtb0(b)
      bal0(n+1:)= 0;  bal(m+1:)= 0

!     Match exactly to a prior ballot & add to its weight, if possible

      Do ib= 1,nb
        l= ballot(0,ib);  If (m /= l) Cycle

        If (All(bal(1:m) == ballot(1:m,ib))) then
          wtb(ib)= wtb(ib) + wt;  Cycle Ballot_loop1
        End if
      End do

!     Else make it a new ballot

      nb= nb + 1;  wtb(nb)= wt
      ballot(:m,nb)= bal(:m)
    End do Ballot_loop1
    
    tot_wt= Sum(wtb(:nb))

    If (pr_out >= 1) then 
      Call Out ("In 'Domain_ballots' total ballot weight reduced from",np, &
                "to",tot_wt,ln=1)
      Call Out ("with initial # ballots",nb0, "and final #",nb)
      Call Out ("for domain",domain)
    End if
    
  End Subroutine Domain_ballots


   Subroutine Read_ballots1 (District,idist,wt_dis, nc0,nc,np,mr,nb, &
                             orig0, wtb,ballot,ballot2, Avg)
   
!    Read ranking ballots only, depending on the format. Either preprocess (consolidate) 
!    and store the preprocessed ballots or read previously preprocessed ballots.

     Character(*), Intent(in) :: District     ! Voting district Name
     Integer,      Intent(in) :: idist        ! Voting district index
     Real,         Intent(in) :: wt_dis(0:,:) ! (mr_spec,mr_spec) Decreasing weight fractions for
                                              !  ballot rankings of candidates. (i,nr) = fraction
                                              !  for candidate ballot(i,b) when ballot(0,b) = nr,
                                              !  with wt_dis(0,nr) = Sum(wt_dis(1:nr,nr)) = 1 
                                              !  if nr >= UVP, else wt_dis(0,nr) = undervote penalty < 1.
                                              !  All weight on top ranked candidate if wt_dis(1,1) <= 0.

     Integer, Intent(out) :: nc0           ! original # candidates
     Integer, Intent(out) :: nc            ! # candidates used for the district
     Integer, Intent(out) :: np            ! # candidates to be elected
     Integer, Intent(out) :: mr            ! max # candidates ranked
     Integer, Intent(out) :: nb            ! # ballots read and consolidated

     Integer,     Pointer :: orig0(:)      ! (nc) Original 'nc0' candidates 
     
     Real,        Pointer :: wtb(:)        ! (nb) Ballot weights summing to 'np'
     Integer,     Pointer :: ballot(:,:)   ! (0:mr,nb)  (1:nr,b) = candidates in preferential order
                                           !   (0,b) = 'nr' = # ranked
     Real,        Pointer :: ballot2(:,:)  ! (mr,nb)  (:nr,b) = weight fractions for the candidate rankings, 
                                           !                    using wt_dis
     Real,        Pointer :: Avg(:)        ! (nc) Average candidate weights = global(:,1)
 
!  Local:
     Real,    Allocatable :: global(:,:)   ! (nc,5)  1 = global mean rating vector, 2 = its sigma vector
                                           ! For rankings which are not top ranked, or ratings which are not top
                                           !   or bottom rated, representing mostly 'noise':
                                           ! 3 = noise mean rating vector, 4 = corresponding sigma vector
                                           ! 5 = mean rating above the noise, normalized
     Character(80) :: msg
     Character(18) :: Dist     ! District of the election (inside the region)
     Character(16) :: Region   ! Region of the election
     Character(65) :: label
     Integer :: Year           ! Year of the election
     Integer :: tot_cnt        ! Total original # single ballots counted
     Integer :: nbf            ! Total # original ballots processed, possibly with duplicates combined

     Real    :: x1, x2
     Integer :: ios, nb0
     Integer :: rt, mt, n1, n2, mr0
     
     If (Associated(orig0))     DeAllocate(orig0)
     If (Associated(wtb))       DeAllocate(wtb)
     If (Associated(ballot))    DeAllocate(ballot)
     If (Associated(ballot2))   DeAllocate(ballot2)
     If (Associated(Avg))       DeAllocate(Avg)
     If (Associated(Req_cand))  DeAllocate(Req_cand)
     
     Call Out ("Read ballots and do simple consolidation for district #",idist,ln=1)
     
!    First read the ballot parameters

     If (Rating == 0) then
       Open(9, File=Indat//District//'2', IOmsg=msg, IOstat=ios, Status='Old', Action='Read')
     Else
       Open(9, File=Indat//District//'1', IOmsg=msg, IOstat=ios, Status='Old', Action='Read')
     End if

!             3 typical labels to be read:
!   'Rating type, # personal ballots, # processed, final # combined: '
!   '# candidates to be elected, # reduced, original #: '
!   'Max # ranked or rated, max # rating levels, undervote penalty: '

     Read(9,*) Dist, Region, Year
     Read(9,*) Label, rt, tot_cnt, nbf, nb0
     Read(9,*) Label, np, nc, nc0
     Read(9,*) Label, mr0, mt

     mr= Min(mr_spec, mr0, nc);  mt= Max(mt,1)  ! mt = 1 for Rating = -1

     If (rt /= Rating) then
       Call Out("Error in 'Read_ballots'. Rating mismatch. Assumed value", Rating, &
                "read from file",rt, ln=1)
       Stop
     Else if (np < 2 .or. nc > nc0 .or. nc <= np .or. nb0 < nc .or. mr < 3 .or. mr > nc) then
       Call Out ("Input parameter error in 'Read_ballots': Reduced # candidates",nc, &
                 "vs orig # cand",nc0, ln=1)
       Call Out ("# candidates to be elected",np, "Max # candidates to be ranked or rated",mr)
       Call Out ("# combined ballots",nb0)
       Stop
     End if
       
!    Typcial label: 'Original candidates, as reordered'

     Allocate(orig0(nc));  Read(9,*) Label, orig0

     If (Rating == 0) then  ! Weak rankig
       Read(9,*) label
     End if
       
!    Memory allocation

     Allocate(Req_cand(nc));  Req_cand= .false.

     Allocate(wtb(nb0), ballot(0:mr,nb0), ballot2(mr,nb0), global(nc,5))

!    Now read the pre-processed ballots

     Call Read_rankings2 (Rating,np,mr,nc,nb0, wt_dis, &
                          nb,wtb,ballot,ballot2, global)
     Close (9)

     Allocate(Avg(nc));  Avg= global(:,1)

!    Noise levels for subtracting from mean rating vectors
!    for correlations (Noise_cor) cluster portions (Noise_por)

     x1= Dot_product(global(:,1), global(:,5))
     x2= Dot_product(global(:,3) + Noise_fac * global(:,4), global(:,5))
     Noise_cor= x2;  Noise_por= Por_fac * Noise_cor

     Call Out ("For district "//Dist//" #", idist, ln=1)

     Call Out ("# original candidates",nc0, "# reduced",nc)
     Call Out ("# candidates to be elected",np, "max # ranked or rated",mr)
     Call Out ("# initial ballots",nb)
     Call Out ("Original candidates, as reordered",orig0)
     Call Out ("Their Avg values",Avg)

     Call Out ("1st estimate of noise to subtract for correlations", x1, &
               "2nd estimate (standard)",x2)
     Call Out ("Noise to subtract for correlations", Noise_cor, &
               "noise cutoff for portions",Noise_por)
   End Subroutine Read_ballots1


  Subroutine Domain_ballots1 (np,mr,nc,nco,domain, nb,wtb,ballot,ballot2, &
                              tot_wt, Avg,Ord)

!   Read in all ballots but record only those which rank at least one candidate
!   in the domain 'domain', eliminating rankings of all candidates outside the domain.
  
    Integer,  Intent(in) :: np            ! # candidates to be elected = weight sum
    Integer,  Intent(in) :: mr            ! Max # candidates ranked per output ballot
                                          !   mr <= nc required
    Integer,  Intent(in) :: nc            ! Reduced # candidates
    Integer,  Intent(in) :: nco           ! Original # candidates
    Integer,  Intent(in) :: domain(:)     ! (nc) Subset of candidates to be used

    Integer,  Intent(in) :: nb            ! original # ballots
    Real,     Intent(in) :: wtb(:)        ! (nb) Ballot weights summing to 'np'
    Integer,  Intent(in) :: ballot(0:,:)  ! (0:mr,nb)  (1:n,b) = candidates in preferential order
                                          !   (0,b) = 'n' = # ranked, possible n > mr 
    Real,     Intent(in) :: ballot2(:,:)  ! (mr,nb)  (:nr,b) = weight fractions for the candidate rankings, 
                                          !                    derived from wt_dis

    Real,    Intent(out) :: tot_wt        ! Total weight of reduced ballots = Sum(wtb(:nb))
    Real,    Intent(out) :: Avg(:)        ! (nc) Average candidate weights on the 'domain'
    Integer, Intent(out) :: Ord(:)        ! (nc) Ordering of these weights
! Local:
    Real, Parameter :: eps= 0.0001
    Real    :: mean(nc), balR(mr), w1(nb), ball2(mr,nb)
    Integer :: bal0(mr), bal(0:mr), inv(nco), ball1(0:mr,nb)

    Logical :: ReOrd
    Real    :: sm, wt, fac
    Integer :: b, i, l, m, n, ib, mb

    If (pr_out >= 1) then
      Call Out ("Enter 'Domain_ballots1', reducing # candidates from",nco, "to",nc, ln=1)
      Call Out ("with max # ranked",mr, "for # open positions",np)
      Call Out ("and domain",domain)
    End if

!   Combine prior ballots which are identical after removal of 
!   outside-the-domain candidates
 
    mb= 0;  w1= 0;  ball1= 0;  ball2= 0;  Avg= 0;  Ord= 0
    Call Inverse_map(domain,inv)

    Ballot_loop1 : Do b= 1,nb
      n= Min(ballot(0,b),mr);  bal0(:n)= inv(ballot(1:n,b))
      Call List_of_true (bal0(:n) > 0, m,bal(1:))  ! bal(1:m)= cand inside the domain
      If (m < 1) Cycle Ballot_loop1

      If (ballot2(2,b) > 0) then
        balR(:m)= ballot2(bal(1:m),b);  balR(m+1:)= 0;  fac= Sum(balR(:m))
      Else
        balR(1)= ballot2(1,b);  balR(2:)= 0;  fac= balR(1)
      End if

      bal(0)= m;  bal(1:m)= bal0(bal(1:m));  bal(m+1:)= 0

      wt= fac * wtb(b)

!     Match exactly to a prior ballot & add to its weight, if possible

      Do ib= 1,mb
        l= ball1(0,ib);  If (m /= l) Cycle
        If (Any(bal(1:m) /= ball1(1:m,ib))) Cycle

        If (balR(2) > 0) then
          sm= Sum(Abs(balR(:m) - ball2(:m,ib)))
        Else
          sm= Abs(balR(1) - ball2(1,ib))
        End if
        If (sm > eps) Cycle

        w1(ib)= w1(ib) + wt;  Cycle Ballot_loop1
      End do

!     Else make it a new ballot

      mb= mb + 1;  w1(mb)= wt
      ball1(:m,mb)= bal(:m)
      ball2(:m,mb)= balR(:m)
    End do Ballot_loop1
    
    tot_wt= Sum(w1(:mb))

    Ballot_loop2 : Do b= 1,mb
      n= ball1(0,b);  mean= 0
      mean(ball1(1:n,b))= ball2(1:n,b)
      Avg= Avg + mean * w1(b)
    End do Ballot_loop2

    Avg= Avg / tot_wt;  mean= Avg
    Call Sort (.false.,mean,Ord, ReOrd)

    If (pr_out >= 1) then 
      Call Out ("In 'Domain_ballots1' total ballot weight reduced from",np, &
                "to",tot_wt,ln=1)
      Call Out ("with # ballots",mb)
      Call Out ("Average candidate weights",Avg)
      If (ReOrd) Call Out ("Average weight reordering",Ord)
    End if
    
  End Subroutine Domain_ballots1

            
  Subroutine Read_rankings2 (Rating, np,mr,nc,nb0, wt_dis, &
                             nb,wtb, ballot,ballot2, global)
  
!   Read the ranked or rated candidates, plus their rankings or ratings, 
!   from a ballot file pre-processed to be in the single line format 
!   used by the subroutine 'Read_vote'
  
    Integer,  Intent(in) :: Rating        ! = 0 for weak ranking 
                                          ! = -1 for strong ranking
    Integer,  Intent(in) :: np            ! # candidates to be elected = total weight
    Integer,  Intent(in) :: mr            ! Max # candidates to be ranked or rated
    Integer,  Intent(in) :: nc            ! # candidates
    Integer,  Intent(in) :: nb0            ! # ballots to be read
    Real,     Intent(in) :: wt_dis(0:,:)  ! (0:mr_spec,mr_spec) Decreasing weight fractions for
                                          !  ballot rankings of candidates. (i,nr) = fraction
                                          !  for candidate ballot(i,b) when ballot(0,b) = nr,
                                          !  with wt_dis(0,nr) = Sum(wt_dis(1:nr,nr)) = 1 
                                          !  if nr >= UVP, else wt_dis(0,nr) = undervote penalty < 1.
                                          !  All weight on top ranked candidate if wt_dis(1,1) <= 0.

    Integer, Intent(out) :: nb            ! Final # ballots after any
                                          !   consolidation, if mr < mr0
    Real,    Intent(out) :: wtb(:)        ! (nb) Ballot weights summing to 'np'
    Integer, Intent(out) :: ballot(0:,:)  ! (0:mr,nb)  (1:nr,b) = candidates in preferential order
                                          !   (0,b) = 'nr' = # ranked or rated  
    Real,    Intent(out) :: ballot2(:,:)  ! (mr,nb) (1:nr,b) =  weight fractions for the restriccted 
                                          !                     candidate rankings
    Real,    Intent(out) :: global(:,:)   ! (nc,5)  1 = global mean rating vector, 2 = its sigma vector
                                          ! For rankings which are not top ranked, or ratings which are not top
                                          !   or bottom rated, representing mostly 'noise':
                                          ! 3 = noise mean rating vector, 4 = corresponding sigma vector
                                          ! 5 = mean rating above the noise, normalized
! Local:
    Real, Parameter :: eps= 0.0001

    Real    :: balR(mr), tmp1(nc), mean(nc)
    Integer :: bal(0:mr), bal2(0:mr), key(nc), Nrk(mr)

    Real    :: sm, wt
    Integer :: b, i, j, n, ib, lp, lv, ml, nr, nz, nw, orig_nb

    Call Out ("Enter 'Read_rankings2'")
 
    nb= 0;  wtb= 0;  ballot= 0;  ballot2= 0;  global= 0

!   Combine initial ballots that have identical candidate rankings
!   and compute their weight distribtions 'ballot2'

    Ballot_loop1 : Do b= 1,nb0

      Call Read_vote (Rating,b,mr,nc, nw, bal,bal2)
    
!     Ballot weights by rankings, as determined from wt_dis

      If (wt_dis(1,1) > 0) then
        If (Rating == 0) then
            lv= Maxval(bal2(1:nr));  ml= Max(lv,nr)
          balR(:nr)= wt_dis(bal2(1:nr),ml)
            sm= Sum(balR(:nr))
          balR(:nr)= (wt_dis(0,nr) / sm) * balR(:nr)
        Else
          balR(:nr)= wt_dis(1:nr,nr)
        End if
      Else
        balR(1)= wt_dis(0,nr);  balR(2:)= 0
      End if

      wt= nw * wt_dis(0,nr)

!     Combine initial ballots which are identical

      Do ib= 1,nb
        n= ballot(0,ib);  If (n /= nr) Cycle
        If (Any(bal(1:n) /= ballot(1:n,ib))) Cycle

        If (Rating == 0 .and. wt_dis(1,1) > 0 ) then
          sm= Sum(Abs(balR(:n) - ballot2(:n,ib)))
          If (sm > eps) Cycle
        End if

        wtb(ib)= wtb(ib) + wt;  Cycle Ballot_loop1
      End do

!     New ballot

      nb= nb + 1;  wtb(nb)= wt
      ballot(:nr,nb) = bal(:nr)
      ballot2(:nr,nb)= balR(:nr)

    End do Ballot_loop1

    orig_nb= Sum(wtb(:nb));  wt= Real(np);  wtb(:nb)= (wt/orig_nb) * wtb(:nb)
    
    Call Out ("Consolidated # ballots",nb, "vs prior #",nb0, ln=1)
    Call Out ("Original total # ballots", orig_nb, "normalized to weight",wt)

    nz= Nint(0.5*mr)

    Ballot_loop2 : Do b= 1,nb
      nr= ballot(0,b);  mean= 0;  mean(ballot(1:nr,b))= ballot2(:nr,b)

      global(:,1)= global(:,1) + mean * wtb(b)
      global(:,2)= global(:,2) + mean**2 * wtb(b)
      
!     Focus on the noise by zeroing top rankings

      tmp1= mean;  Call Sort (.false.,tmp1,key)
      mean(key(:nz))= 0

      global(:,3)= global(:,3) + mean * wtb(b)
      global(:,4)= global(:,4) + mean**2 * wtb(b)
    End do Ballot_loop2

    global(:,1)= global(:,1) / np
    global(:,2)= Sqrt(global(:,2) / np - global(:,1)**2)
    
    global(:,3)= global(:,3) / np
    global(:,4)= Sqrt(global(:,4) / np - global(:,3)**2)

    tmp1= global(:,1) - global(:,3);  global(:,5)= tmp1 / Sum(tmp1)

    Do n= 1,mr
      Nrk(n)= Count(ballot(0,:) == n)
    End do

    Call Out ("# initial ballots for each # ranked",Nrk(:mr))
    Call Out (-1,"Global: mean & sig, noise & sig, mean - noise (normed)",global)
    
  End Subroutine Read_rankings2
   

   Subroutine Read_ballots0 (District,idist,nc_lim, nc,np,mr,nb, wtb,ballot)
   
!    Read ranking ballots only, depending on the format. Either preprocess (consolidate) 
!    and store the preprocessed ballots or read previously preprocessed ballots.

     Character(*), Intent(in) :: District  ! Voting district Name
     Integer,      Intent(in) :: idist     ! Voting district index
     Integer,      Intent(in) :: nc_lim    ! Last candidate to be used

     Integer, Intent(out) :: nc            ! # candidates used for the district
     Integer, Intent(out) :: np            ! # candidates to be elected
     Integer, Intent(out) :: mr            ! max # candidates ranked
     Integer, Intent(out) :: nb            ! # ballots read and consolidated

     Real,        Pointer :: wtb(:)        ! (nb) Ballot weights summing to 'np'
     Integer,     Pointer :: ballot(:,:)   ! (0:mr,nb)  (1:nr,b) = candidates in preferential order
                                           !   (0,b) = 'nr' = # ranked
!  Local:
     Real,    Allocatable :: global(:,:)   ! (nc,5)  1 = global mean rating vector, 2 = its sigma vector
                                           ! For rankings which are not top ranked, or ratings which are not top
                                           !   or bottom rated, representing mostly 'noise':
                                           ! 3 = noise mean rating vector, 4 = corresponding sigma vector
                                           ! 5 = mean rating above the noise, normalized
     Real,    Allocatable :: rank_pt(:), unrank_pt(:) ! (mr)

     Character(80) :: msg
     Character(18) :: Dist     ! District of the election (inside the region)
     Character(16) :: Region   ! Region of the election
     Character(65) :: label
     Integer :: Year           ! Year of the election
     Integer :: tot_cnt        ! Total original # single ballots counted
     Integer :: nbf            ! Total # original ballots processed, possibly with duplicates combined
     Integer :: nb0            ! # ballots to be read

     Real    :: x1, x2
     Integer :: mt, rt, ios, mr0, nc0
     
     If (Associated(wtb))    DeAllocate(wtb)
     If (Associated(ballot)) DeAllocate(ballot)
     If (Associated(Req_cand)) DeAllocate(Req_cand)
     
     Call Out ("Read ballots and do simple consolidation for district #",idist,ln=1)
     
!    First read the ballot parameters

     Open(9, File=Indat//District//'1', IOmsg=msg, IOstat=ios, Status='Old', Action='Read')

!             3 typical labels to be read:
!   'Rating type, # personal ballots, # processed, final # combined: '
!   '# candidates to be elected, # reduced, original #: '
!   'Max # ranked or rated, max # rating levels, undervote penalty: '

     Read(9,*) Dist, Region, Year
     Read(9,*) Label, rt, tot_cnt, nbf, nb0
     Read(9,*) Label, np, nc, nc0
     Read(9,*) Label, mr0, mt
     Read(9,*)  ! Ignore list of original candidates

     nc= Min(nc,nc_lim);  mr= Min(mr_spec, mr0, nc)

     If (rt /= -1) then
       Call Out("Error in 'Read_ballots0'. Rating mismatch")
       Stop
     Else if (np < 2 .or. nc <= np .or. nb0 < nc .or. mr < 3) then
       Call Out ("Input parameter error in 'Read_ballots0': for reduced # candidates",nc, &
                 "or # to be elected",np, ln= 1) 
       Call Out ("Max # candidates to be ranked",mr, "# combined ballots",nb0)
       Stop
     End if

     Allocate(Req_cand(nc));  Req_cand= .false.

     Allocate(rank_pt(mr), unrank_pt(mr))
     Call Ranking_pt0 (Parm_bal1,Max_pt,nc,mr, rank_pt,unrank_pt)
       
     Allocate(wtb(nb0), ballot(0:mr,nb0), global(nc,5))

!    Read the pre-processed ballots

     Call Read_ranked0 (np,mr,nc,nb0, rank_pt,unrank_pt, nb,wtb, ballot, global)

     Close (9)

!    Noise levels for subtracting from mean rating vectors
!    for correlations (Noise_cor) cluster portions (Noise_por)

     x1= Dot_product(global(:,1), global(:,5))
     x2= Dot_product(global(:,3) + Noise_fac * global(:,4), global(:,5))
     Noise_cor= Min(x1,x2);  Noise_por= Por_fac * Noise_cor

     Call Out ("For district "//Dist//" #", idist, ln=1)
     Call Out ("reduced # canididates",nc)
     Call Out ("# candidates to be elected",np, "max # ranked",mr)
     Call Out ("# consolidated ballots",nb)

     Call Out ("1st estimate of noise to subtract for correlations", x1, &
               "2nd estimate (standard)",x2)
     Call Out ("Noise to subtract for correlations", Noise_cor, &
               "noise cutoff for portions",Noise_por)
   End Subroutine Read_ballots0
            
  Subroutine Read_ranked0 (np,mr,nc,nb0, rank_pt,unrank_pt, nb,wtb, ballot, global)
  
!   Read the strongly candidates from a preprocessed file using the subroutine 'Read_vote0'

!   Restrict all ballots to candidates <= 'nc', assuming mr <= nc
  
    Integer,  Intent(in) :: np              ! # candidates to be elected = total weight
    Integer,  Intent(in) :: mr              ! Max # candidates to be ranked
    Integer,  Intent(in) :: nc              ! Last candidate to be used
    Integer,  Intent(in) :: nb0             ! # ballots to be read

    Real,     Intent(in) :: rank_pt(:)      ! (mr) Point values for the top candidates
    Real,     Intent(in) :: unrank_pt(:)    ! (mr) Point values of unranked candidates

    Integer, Intent(out) :: nb              ! Final # ballots after any consolidation
    Real,    Intent(out) :: wtb(:)          ! (nb) Ballot weights summing to 'np'
    Integer, Intent(out) :: ballot(0:,:)    ! (0:mr,nb)  (1:nr,b) = candidates in preferential order
                                            !   (0,b) = 'nr' = # ranked or rated  
    Real,    Intent(out) :: global(:,:)     ! (nc,5)  1 = global mean rating vector, 2 = its sigma vector
                                            ! For rankings which are not top ranked, or ratings which are not top
                                            !   or bottom rated, representing mostly 'noise':
                                            ! 3 = noise mean rating vector, 4 = corresponding sigma vector
                                            ! 5 = mean rating above the noise, normalized
! Local:
    Real, Parameter :: eps= 0.000001

    Real    :: tmp1(nc), mean(nc)
    Integer :: bal(0:mr), key(nc), Nrk(mr)

    Real    :: smz, wt
    Integer :: b, n, ib, ml, nr, nz, nw, orig_nb

    Call Out ("Enter 'Read_ranked0'")
 
    nb= 0;  wtb= 0;  ballot= 0;  global= 0

!   Combine initial ballots that have identical candidate rankings

    Ballot_loop1 : Do b= 1,nb0

      Call Read_vote0 (b,nc, nw,bal)

      nr= bal(0);  If (nr < 1) Cycle Ballot_loop1

!     Combine initial ballots which are identical

      Do ib= 1,nb
        n= ballot(0,ib);  If (n /= nr) Cycle
        If (Any(bal(1:n) /= ballot(1:n,ib))) Cycle
        wtb(ib)= wtb(ib) + nw;  Cycle Ballot_loop1
      End do

!     New ballot

      nb= nb + 1;  wtb(nb)= nw
      ballot(:nr,nb)= bal(:nr)

    End do Ballot_loop1

    orig_nb= Sum(wtb(:nb));  wt= Real(np);  wtb(:nb)= (wt/orig_nb) * wtb(:nb)
    
    Call Out ("Consolidated # ballots",nb, "vs prior #",nb0, ln=1)
    Call Out ("Original total # ballots", orig_nb, "normalized to weight",wt)

    nz= Nint(0.50 * mr);  smz= 0

    Ballot_loop2 : Do b= 1,nb
      nr= ballot(0,b);  mean= unrank_pt(nr)
      mean(ballot(1:nr,b))= rank_pt(:nr)

      global(:,1)= global(:,1) + mean * wtb(b)
      global(:,2)= global(:,2) + mean**2 * wtb(b)
      
!     Focus on the noise by zeroing top rankings

      If (nr > nz) then
        tmp1= mean;  Call Sort (.false.,tmp1,key)
        mean(key(:nz))= 0;  smz= smz + wtb(b)

        global(:,3)= global(:,3) + mean * wtb(b)
        global(:,4)= global(:,4) + mean**2 * wtb(b)
      End if
    End do Ballot_loop2

    global(:,1)= global(:,1) / np
    global(:,2)= Sqrt(global(:,2) / np - global(:,1)**2)
    
    smz= Max(smz,eps)
    global(:,3)= global(:,3) / smz
    global(:,4)= Sqrt(global(:,4) / smz - global(:,3)**2)

    tmp1= global(:,1) - global(:,3);  global(:,5)= tmp1 / Sum(tmp1)

    Do n= 1,mr
      Nrk(n)= Count(ballot(0,:) == n)
    End do

    Call Out ("# initial ballots for each # ranked",Nrk(:mr))
    Call Out (-1,"Global: mean & sig, noise & sig, mean - noise (normed)",global)
    
  End Subroutine Read_ranked0

End Module Clusters1