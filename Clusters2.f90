
!            This module contains "Consolidate_ballots" and its subroutines. 
    
!  "Consolidate_ballots" is called from Clustering_PR to consolidate the pre-processed ballots
!   into slate consolidated ballots, finally into the slate clusters 'Memb' to be used for 
!   the clustering.
    
Module Clusters2
    
   Use Clusters0
   Use Clust_gen
   Use Clusters_support
   Use Graph_algorithms
   Use Factorials
   Use Util
   Use Output
   Use Types
   Use Precisn
   Use IEEE_Arithmetic
   Implicit None

   Logical,Parameter :: Test_var_1= .false.
    
 Contains  
    
   Subroutine Consolidate_ballots (idist, nc,np,mr,nb, pt_val, &
                                   wtb,ballot,ballot2, nr,nsl,Memb, Mean_rnd)
   
!    Consolidate ballots to slate ballots / slate clusters.
!    The membership data for the slate clusters are recorded in 'Memb' for subsequent use
!    in the clustering algorithm. In particular, the final clusters are consist of slate 
!    cluster members, though these slate clusters may have partial membership in more than 
!    one cluster of a cluster set (="fuzzy set" membership).

     Integer,   Intent(in) :: idist         ! Voting district index
     Integer,   Intent(in) :: nc, nb        ! # candidates and # initial ballots
     Integer,   Intent(in) :: np, mr        ! # candidates to be elected and max # ranked

     Real,      Intent(in) :: pt_val(:)     ! (mt) Decreasing point values for 
                                            !   rising rating or ranking levels
     
     Real,      Intent(in) :: wtb(:)        ! (nb) Input ballot weights
     Integer,   Intent(in) :: ballot(0:,:)  ! (0:mr,nb) For each consolidated ballot: (0,b) = # ranked 
                                            !   or rated, (1:mr,b) = candidates in preferential order
     Integer,   Intent(in) :: ballot2(0:,:) ! (0:mr,nb) (1:mr,b) = corresponding rating levels
                                            !   with (0,b) = # positively rated
     Integer,  Intent(out) :: nr            ! 
     Integer,  Intent(out) :: nsl           ! Final # slate clusters

     Type(Multi_listD), Pointer :: Memb(:)  ! (0:ns) Final slate cluster data
                                            !  %l = np = # candidates to be elected = 
                                            !            total weight in position units
                                            !  %n = nc = # candidates 
                                            !  %fsx = tot_wt = total weight of all slate clusters 
                                            !  %qx(2)  = correlation parameters
                                            !  %vl(mlv). (lv) = # slate clusters with
                                            !      # ranked or rated candidates <= lv

                                            ! 1:ns Case:
                                            !  %l   = # top slate candidates
                                            !  %m   = # bottom slate candidates (for Rating = 1)
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
                                            !  %tx(0:nc)= centered (subtract Noise_cor if Rating < 1) 
                                            !             and normalized mean vector
                                            !  %ux(nc) = variance normalized weight %fsx / %px

     Real,  Intent(out) :: Mean_rnd(0:,:,:) ! (0:nc,nr,nt) 'nr' mean rating vectors for
                                            !   'nt' random sets of clusters, 
                                            !   with cluster weights at 0
     
! Local:
     Real,    Allocatable :: rates(:,:)    ! (0:nc,nb) Ballot weights at 0, 
                                           !   mean rating vectors at 1:nc
     Real,    Allocatable :: cls_mean(:,:) ! (0:nc,lim) Mean rating vectors for the
                                           !   candidate clusters, cluster weight at 0
     Real,    Allocatable :: sl_wt(:)      ! (msl)     Weights of slate ballots
     Integer, Allocatable :: sl_lv(:,:)    ! (0:2,msl) Levels, etc, of slate ballots

     Integer, Allocatable :: Sl_nlv(:)     ! (mlv) Sum level counter: the last slate ballot 
                                           !       for each sum level
     Integer :: consol(nb)    ! Maps initial ballots to slate ballots
     Real    :: tot_wt
     Integer :: mlv, m1, m2, mt, nl, ns, nt, lim, msl
     
     Call Out ("")
     Call Out ("Enter 'Consolidate_ballots' for district #",idist,ln=1)
     Call Out ("Initial to slate ballots for rating type",Rating)
     
     mt= Size(pt_val);  nr= Size(Mean_rnd,2);  nt= Size(Mean_rnd,3)

     If (Rating == 1) then
       msl= Min(Max_sl1 + nc*(nc-1), nb)
       Min_slate= np / Real(msl)
     Else
       msl= Min(Max_sl, nb)
       Min_slate= np / Real(msl)
     End if
     
     Allocate(sl_wt(msl), sl_lv(0:2,msl))
     
       If (Rating == 1) then
         m1= Maxval(ballot2(0,:));  m2= Maxval(ballot(0,:) - ballot2(0,:))

         Call Consolidate_ratings (m1,m2, np,nb, wtb,ballot,ballot2, &
                                   sl_wt,sl_lv, tot_wt,nsl,consol)

         mlv= Maxval(sl_lv(0,:));  Allocate(Sl_nlv(mlv))
         Call Reord_ratings (nsl, sl_wt(:nsl),sl_lv(:,:nsl), consol, Sl_nlv)

       Else
         If (Rating == 0) then
           Call Consolidate_weak (mr, np,nb, wtb,ballot,ballot2, &
                                  sl_wt,sl_lv(0,:), tot_wt,nsl,consol)
         Else 
           Call Consolidate_strong (mr, np,nb, wtb,ballot, &
                                    sl_wt,sl_lv(0,:), tot_wt,nsl,consol)
         End if

         mlv= Maxval(sl_lv(0,:));  Allocate(Sl_nlv(mlv))
         Call Reord_rankings (nsl, sl_wt(:nsl),sl_lv(0,:nsl), consol, Sl_nlv)

         sl_lv(1,:nsl)= sl_lv(0,:nsl);  sl_lv(2,:nsl)= 0
       End if

       Call Out("For district #", idist, "# slate clusters",nsl, ln=1)
     
!      Compute the slate cluster membership data 'Memb'
       
       Allocate (Memb(0:nsl), Memb(0)%qx(2), Memb(0)%vl(mlv))
       Memb(0)%l= np;  Memb(0)%n= nc;  Memb(0)%fsx= tot_wt
       Memb(0)%qx= Parm_cor;  Memb(0)%vl= Sl_nlv
       
       Memb(1:)%fsx= sl_wt(:nsl)
       Memb(1:)%l= sl_lv(1,:nsl)
       Memb(1:)%m= sl_lv(2,:nsl)

       Call Slate_clust (nc,mr,mt, nb,nsl, pt_val,wtb, &
                          ballot,ballot2, consol, Memb)

!      Compute the clusters to be used to form the initial cluster sets 
    
       Allocate(cls_mean(0:nc,nsl))

       If (Rating == 1) then
         Call Candidate_clusters2 (nc,nsl, Memb, cls_mean)
       Else 
         Call Candidate_clusters (nc,nsl, Memb, cls_mean)
       End if
       
       nr= Min(nr,nsl-2)
       Call Select_init1 (np,nr, cls_mean, Mean_rnd(:,:nr,:))

   End Subroutine Consolidate_ballots

                          
   Subroutine Slate_clust (nc,mr,mt, nb,ns, pt_val, wtb,&
                           ballot,ballot2, consol, Memb)
                                     
!    Compute the slate ballot data structure 'Memb', computing 
!    the maximum liklihood estimates of the mean vectors of the 
!    slate clusters computed from the slate ballots.
  
     Integer,  Intent(in) :: nc             ! # candidates
     Integer,  Intent(in) :: mr             ! max # candidates ranked or rated
     Integer,  Intent(in) :: mt             ! # rating point values

     Integer,  Intent(in) :: nb             ! # ballots
     Integer,  Intent(in) :: ns             ! # slate ballots
     Real,     Intent(in) :: pt_val(:)      ! (mt) Decreasing point values of the rating levels

     Real,     Intent(in) :: wtb(:)         ! (nb) Ballot weights summing to 'np'
     Integer,  Intent(in) :: ballot(0:,:)   ! (0:mr,nb)  (1:nr,b) = candidates in preferential order
                                            !   (0,b) = 'nr' = # ranked or rated  
     Integer,  Intent(in) :: ballot2(0:,:)  ! (0:mr,nb) (1:nr,b) =  increasing ranking or rating levels
                                            !   (0,b) = # positive ratings 

     Integer,  Intent(in) :: consol(:)      ! (nb Maps initial ballots to slate ballots

     Type(Multi_listD), Intent(inout) :: Memb(0:) ! (0:ns) Final slate cluster data
                                                  !  %l = np = # candidates to be elected = 
                                                  !            total weight in position units
                                                  !  %n = nc = # candidates 
                                                  !  %qx(2)  = correlation parameters
                                                  !  %vl(mlv). (lv) = # slate clusters with
                                                  !      # ranked or rated candidates <= lv

                                                  ! 1:ns Case:
                                                  !  %l   = # top slate candidates
                                                  !  %m   = # bottom slate candidates (for Rating = 1)
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
                                                  !  %tx(0:nc)= centered (subtract Noise_cor if Rating < 1) 
                                                  !             and normalized mean vector
                                                  !  %ux(nc) = variance normalized weight %fsx / %px
! Local:
     Real, Parameter :: eps= 1.0E-13
     Real(Dblp) :: rt1, cor, dif, rate(0:nc,ns), sm(nc), sm1(nc), sm2(nc), wt_var(nc)
     Real    :: mean(nc), var_pt(nc), rt_var(0:mt)
     Real    :: rank_pt(mr,2), unrank_pt(mr,2)

     Integer :: lst(nb), key(ns)
     Integer :: b, j, n, lm, ln, lv, n1, n2, nr, s1, s2, sl, neg
     
     Call Out ("Enter 'Slate_clust' to compute the slate ballot membership data structure")

     If (Rating <= 0) then
       Call Ranking_pts (Parm_bal,Max_pt,nc,mr, rank_pt,unrank_pt)
       If (mt > 1)  Call Rating0_pts (Parm_bal1,mt,pt_val, rt_var)
     Else
       Call Rating_pts (Parm_bal,mt,pt_val, rt_var)
     End if
     
!    Loop over all slate ballots

     dif= Max_pt
     Slate_loop : Do sl= 1,ns

       Allocate(Memb(sl)%lt(nc), Memb(sl)%px(nc), Memb(sl)%rx(nc), &
                Memb(sl)%sx(nc), Memb(sl)%tx(0:nc), Memb(sl)%ux(nc))
       Memb(sl)%lt= 0; Memb(sl)%px= 0; Memb(sl)%rx= 0
       Memb(sl)%sx= 0; Memb(sl)%tx= 0; Memb(sl)%ux= 0
       
       sm= eps;  sm1= 0;  sm2= 0
       Call List_of_true (consol == sl, ln, lst)

       Combined_ballot_loop : Do j= 1,ln  ! Loop over the matching consolidated ballots
         b= lst(j);  nr= ballot(0,b)

         If (Rating == -1) then
           mean= unrank_pt(nr,1);  var_pt= unrank_pt(nr,2)
           mean(ballot(1:nr,b))  = rank_pt(:nr,1)
           var_pt(ballot(1:nr,b))= unrank_pt(:nr,2)
         Else if (Rating == 0 .and. mt == 1) then
           mean= unrank_pt(nr,1);  var_pt= unrank_pt(nr,2)
           mean(ballot(1:nr,b))  = rank_pt(ballot2(1:nr,b),1)
           var_pt(ballot(1:nr,b))= rank_pt(ballot2(1:nr,b),2)
         Else
           mean= 0;  var_pt= rt_var(0)
           mean(ballot(1:nr,b))  = pt_val(ballot2(1:nr,b))
           var_pt(ballot(1:nr,b))= rt_var(ballot2(1:nr,b))
         End if

         If (Test_var_1) var_pt= 1

         wt_var= wtb(b) / var_pt
         sm = sm + wt_var
         sm1= sm1 + mean * wt_var
         sm2= sm2 + (var_pt + mean**2) * wt_var
       End do Combined_ballot_loop
       
!      Compute the mean vector of the slate ballot and its variance
       sm1= sm1 / sm;  sm2= sm2 / sm
       Memb(sl)%sx= sm1
       Memb(sl)%px= sm2 - sm1**2

!      For ranking ballots, adjust the mean vectors to subtract
!      out the noise for correlation computations, using %tx

       If (Rating < 1) then
         rate(1:,sl)= sm1 - Noise_cor

         rt1= Maxval(rate(1:,sl))
         If (rt1 < Min_mx) then
           Call Out ("Warning in 'Slate_clust': Low noise adjusted rating",rt1, ln=1)
           Call Out ("at slate",sl)
           rate(1:,sl)= rate(1:,sl) + (Min_mx - rt1)
         End if
       Else
         rate(1:,sl)= sm1
       End if
     End do Slate_loop

!    Correct the noise significance level 'Noise_cor', if necessary
!    (unlikely) so that no slate cluster is left 'in the noise'

     Call Normalize_vec (rate)

     Do sl= 1,ns 
       Memb(sl)%tx= rate(:,sl)
       Memb(sl)%ux= Memb(sl)%fsx / Memb(sl)%px
       Memb(sl)%rx= Memb(sl)%sx * Memb(sl)%ux
     End do

!    List the slate ballots positively correlated with each slate ballot
!    and subsequent to it, plus the correlation, ordered by decreasing 
!    correlation

     lst= 0;  Memb(1:)%k= 0
     Do s1= 1,ns-1
       j= 0
       Do s2= s1+1,ns 
         Call Dot_product_M (Dot_fac, Memb(s1)%tx(1:),Memb(s2)%tx(1:), neg, cor)

         If (cor > Memb(0)%qx(1)) then
           j= j + 1;  lst(j)= s2;  rate(0,j)= cor
         End if
       End do

       If (j > 0) then
         Memb(s1)%k= j;  Allocate(Memb(s1)%vl(j), Memb(s1)%qx(j))
         Memb(s1)%qx= rate(0,:j)
         Call Sort (.false., Memb(s1)%qx, key(:j))
         Memb(s1)%vl= lst(key(:j))
       End if
     End do

     Call Form_Memb_sub (nc,ns, Memb)
     
   End Subroutine Slate_clust
                          
                          
   Subroutine Form_Memb_sub (nc,ns, Memb)
                                     
!    Compute the slate ballot data structure 'Memb', converting the ballot weights 
!    to position units and computing the maximum liklihood estimates
!    of the mean rating vectors, plus associated sigma vectors, for the 
!    final slate ballots as clusters.
   
!    Also estimate the noise level of the slate ballot mean rating vectors   
  
     Integer,              Intent(in) :: nc       ! # candidates
     Integer,              Intent(in) :: ns       ! # slate clusters

     Type(Multi_listD), Intent(inout) :: Memb(0:) ! (0:ns) Final slate cluster data
                                                  !  %l = np = # candidates to be elected = 
                                                  !            total weight in position units
                                                  !  %n = nc = # candidates 
                                                  !  %qx(2)  = correlation parameters
                                                  !  %vl(mlv). (lv) = # slate clusters with
                                                  !      # ranked or rated candidates <= lv

                                                  ! 1:ns Case:
                                                  !  %l   = # top slate candidates
                                                  !  %m   = # bottom slate candidates (for Rating = 1)
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
                                                  !  %tx(0:nc)= centered (subtract Noise_cor if Rating < 1) 
                                                  !             and normalized mean vector
                                                  !  %ux(nc) = variance normalized weight %fsx / %px
! Local:
     Real, Parameter :: del= 0.01
     Real(Dblp) :: wt(nc), var(nc), rate(0:nc,ns)
     Integer       :: nu, sl, key(nc)

     If (pr_out >= 1) Call Out ("Enter 'Form_Memb_sub'")

     nu= Memb(0)%l + 1
     
     Slate_loop : Do sl= 1,ns

!      Compute the preferential ordering %lt of the candiates, 
!      then the cluster width %fux from the top candidates if known
       
       wt= Abs(Memb(sl)%sx) / Max_pt + del
       Call Sort (.false.,wt, key)
       var(:nu)= Max(Memb(sl)%px(key(:nu)), del)

       Memb(sl)%fux= Sqrt( Sum(wt(:nu))/ Sum(wt(:nu)/ var(:nu)) )  ! Weighted harmonic mean

       wt= Memb(sl)%sx;  Call Sort (.false.,wt, Memb(sl)%lt)  ! Candidates in preferential order

       rate(1:,sl)= Memb(sl)%tx(0) * Memb(sl)%tx(1:);  rate(0,sl)= Memb(sl)%fsx
     End do Slate_loop

     If (pr_out > 1) Call Out (-1, "All slate ballots: weight, centered mean vector",Real(rate))
     
!    Compute the differenced slate ballot mean rating vectors Memb%tx     
     
     Do sl= 1,ns
       Memb(sl)%p= Count(rate(1:,sl) > 0)
       Memb(sl)%n= Count(rate(1:,sl) < 0)
     End do

!    Output slate membership data
     
     If (pr_out > 1.5) then
       Do sl= 1,ns
         Call Out ("Slate cluster",sl, "# positive candidates",Memb(sl)%l, ln=1)
         Call Out ("# negative candidates",Memb(sl)%m)
         Call Out ("# positive ratings",Memb(sl)%p, "# negative ratings",Memb(sl)%n)
         Call Out ("Membership weight",Real(Memb(sl)%fsx), "and width",Real(Memb(sl)%fux))
         Call Out ("Slate ballot mean rating vector", Real(Memb(sl)%sx))
         Call Out ("Corresponding candidates in preferential order", Memb(sl)%lt)                       
         Call Out ("Corresponding centered and normalized vector", Real(Memb(sl)%tx(1:)))
         Call Out ("Corresponding variance vector", Real(Memb(sl)%px))
       End do
     End if
   End Subroutine Form_Memb_sub
                          
    
  Subroutine Candidate_clusters (nc,ns, Memb, cls_mean)
    
!   Compute candidate set clusters (unordered) for ranking data, to initialize 
!   possible final clusters. These clusters top rank 1, 2, or 3 candidates, 
!   using slate cluster data:

!   The level 1 candidate clusters contain the level 2 and 3 slate clasters that
!   match the initial level 1 slate cluster. Likewise the level 2 candidate clusters
!   contain the matching level 3 slate clusters.
                                                       
    Integer,           Intent(in) :: nc       ! # candidates
    Integer,           Intent(in) :: ns       ! # slate clusters = # candidate clusters
    Type(Multi_listD), Intent(in) :: Memb(0:) ! (0:ns) Final slate cluster data
                                              ! 0 Case: 
                                              !   %n = nc = # candidates 
                                              !   %l = np = # candidates to be elected = 
                                              !             total weight in position units
                                              !   %qx(2)  = correlation parameters
                                              !   %vl(mlv). (lv) = # slate clusters with
                                              !      # ranked or rated candidates <= lv

                                              ! 1:ns Case:
                                              !   %fsx= slate membership
                                              !   %fux= slate width
                                              !   %l = # top slate candidates (= level)
                                              !   %m = # bottom slate candidates (for Rating = 1)
                                              !   %p = # positively rated candidates (for %tx)
                                              !   %n = # negatively rated candidates (for %tx)
                                              !   %lt(nc) Candidates in preferential order (from %sx)
                                              !   %px(nc) = variance vector for %sx
                                              !   %rx(nc) = mean vector with variance normalized weight= $sx*%ux
                                              !   %sx(nc) = mean rating vector
                                              !   %tx(0:nc)= centered (subtract Noise_cor if Rating < 1) 
                                              !              and normalized mean vector
                                              !   %ux(nc) = variance normalized weight = %fsx / %px

    Real,      Intent(out) :: cls_mean(0:,:)  ! (0:nc,ns)
                                              ! (1:,:) = candidate cluster mean vectors
                                              ! (0,:) = candidate cluster weights,
                                              !   with clusters sorted by decreasing weight
! Local:
    Type(Multi_listR) :: cand(ns)  ! (ns) Slate cluster candidate ranking
                                   ! %L0(lv,lv) where lv= level = Memb(sl)%l 
    Real    :: cls_mean2(nc,ns)
    Integer :: cnd(nc), lst(ns)
    Integer :: c, i, l1, l2, lv, s0, sl, mlv, sl0, sl1
    
    If (pr_out >= 1) Call Out ("Enter 'Candidate_clusters'")

    cls_mean= 0;  cls_mean2= 0;  cnd= 0

    Do sl= 1,ns
      lv= Memb(sl)%l;  Allocate(cand(sl)%L0(lv,lv));  cand(sl)%L0= 0
      
      cand(sl)%L0(1,1)= Memb(sl)%lt(1)
      Do i= 2,lv
        cnd(:i)= Memb(sl)%lt(:i);  Call Sort (.true.,cnd(:i))
        cand(sl)%L0(:i,i)= cnd(:i)
      End do
    End do

!   Search for all prior candidate clusters which are identical to the current
!   slate cluster at a lower level and add the current slate cluster to these.
!   Start by initializing candidate clusters by the level 1 slate clusters.
!   Then add the level 2 slate clusters, and, finally, the level 3 slate clusters

    Do sl= 1,ns
      cls_mean(0,sl) = Memb(sl)%fsx
      cls_mean(1:,sl)= Memb(sl)%rx
      cls_mean2(:,sl)= Memb(sl)%ux
    End do

    sl1= 0;  mlv= Size(Memb(0)%vl)
    Level_loop : Do lv= 1,mlv-1
      sl0= sl1;  sl1= Memb(0)%vl(lv)

      Do s0= sl0+1,sl1
        l1= Memb(s0)%l
        Do sl= sl1+1,ns
          l2= Memb(sl)%l
          If (All(cand(sl)%L0(:lv,l2) == cand(s0)%L0(:lv,l1))) then
            cls_mean(0,s0) = cls_mean(0,s0)  + Memb(sl)%fsx
            cls_mean(1:,s0)= cls_mean(1:,s0) + Memb(sl)%rx
            cls_mean2(:,s0)= cls_mean2(:,s0) + Memb(sl)%ux
          End if
        End do
      End do
    End do Level_loop
    
    Where (cls_mean2 > 0) cls_mean(1:,:)= cls_mean(1:,:) / cls_mean2

!   Sort the candidate clusters by weight

    Call Sort (.false., cls_mean(0,:), lst)
    Forall(c=1:nc) cls_mean(c,:)= cls_mean(c,lst)

    If (pr_out > 1) then
      Call Out ("Candidate cluster sizes", cls_mean(0,:))
      Call Out (-1,"Candidate cluster mean vectors", cls_mean(1:,:))
    End if
  End Subroutine Candidate_clusters

    
  Subroutine Candidate_clusters2 (nc,ns, Memb, cls_mean)
    
!   Compute candidate set clusters (unordered) for rating data, to initialize 
!   possible final clusters. These clusters top rank 1, 2, or 3 candidates (pos. ratings), 
!   also bottom rank 1, 2, or 3 candidates (neg. ratings), using slate cluster data:

!   The level 1 candidate clusters contain the level 2 and 3 slate clasters that
!   match the initial level 1 slate cluster. Likewise the level 2 candidate clusters
!   contain the matching level 3 slate clusters.
                                                       
    Integer,           Intent(in) :: nc       ! # candidates
    Integer,           Intent(in) :: ns       ! # slate clusters = # candidate clusters
    Type(Multi_listD), Intent(in) :: Memb(0:) ! (0:ns) Final slate cluster data
                                              ! 0 Case: 
                                              !   %n = nc = # candidates 
                                              !   %l = np = # candidates to be elected = 
                                              !             total weight in position units
                                              !   %qx(2)  = correlation parameters
                                              !   %vl(mlv). (lv) = # slate clusters with
                                              !      # ranked or rated candidates <= lv

                                              ! 1:ns Case:
                                              !   %fsx= slate membership
                                              !   %fux= slate width
                                              !   %l = # top slate candidates (= pos. level)
                                              !   %m = # bottom slate candidates (for Rating = 1)
                                              !   %p = # positively rated candidates (for %tx)
                                              !   %n = # negatively rated candidates (for %tx)
                                              !   %lt(nc) Candidates in preferential order (from %sx)
                                              !   %px(nc) = variance vector for %sx
                                              !   %rx(nc) = mean vector with variance normalized weight= $sx*%ux
                                              !   %sx(nc) = mean rating vector
                                              !   %tx(0:nc)= centered (subtract Noise_cor if Rating < 1) 
                                              !              and normalized mean vector
                                              !   %ux(nc) = variance normalized weight = %fsx / %px

    Real,      Intent(out) :: cls_mean(0:,:)  ! (0:nc,ns) Mean rating vectors for the candidate clusters
                                              !   with (0,sl) = weight of candidate cluster 'sl'
! Local:
    Type(Multi_listR) :: cand(ns)  ! (ns) Slate cluster candidate ranking
                                   ! %L1(lp,lp) where lp= pos level = Memb(sl)%l 
                                   ! %L2(ln,ln) where ln= neg level = Memb(sl)%m 
    Real    :: cls_mean2(nc,ns)
    Integer :: cnd(nc), key(ns)
    Integer :: l1p, l1n, l2p, l2n, mlv, sl0, sl1
    Integer :: c, i, j, lp, ln, s0, nm, sl, sm
    
    If (pr_out >= 1) Call Out ("Enter 'Candidate_clusters2")

    cls_mean= 0;  cls_mean2= 0;  cnd= 0

    Do sl= 1,ns
      lp= Memb(sl)%l;  ln= Memb(sl)%m
      Allocate(cand(sl)%L1(lp,lp), cand(sl)%L2(ln,ln))
      cand(sl)%L1= 0;  cand(sl)%L2= 0
      
      cand(sl)%L1(1,1)= Memb(sl)%lt(1)
      Do i= 2,lp
        cnd(:i)= Memb(sl)%lt(:i);  Call Sort (.true.,cnd(:i))
        cand(sl)%L1(:i,i)= cnd(:i)
      End do
      
      cand(sl)%L2(1,1)= Memb(sl)%lt(nc)
      Do i= 2,ln
        j= nc+1 - i
        cnd(:i)= Memb(sl)%lt(j:);  Call Sort (.true.,cnd(:i))
        cand(sl)%L2(:i,i)= cnd(:i)
      End do
    End do

!   Search for all prior candidate clusters which are identical to the current
!   slate cluster at a lower level and add the current slate cluster to these.
!   Start by initializing candidate clusters by the level 1 slate clusters.
!   Then add the level 2 slate clusters, and, finally, the level 3 slate clusters

    Do sl= 1,ns
      cls_mean(0,sl) = Memb(sl)%fsx
      cls_mean(1:,sl)= Memb(sl)%rx
      cls_mean2(:,sl)= Memb(sl)%ux
    End do

    sl1= 0;  mlv= Size(Memb(0)%vl)
    Level_loop : Do sm= 2,mlv-1
      sl0= sl1;  sl1= Memb(0)%vl(sm)

      Do s0= sl0+1,sl1
        l1p= Memb(s0)%l;  l1n= Memb(s0)%m  
        Do sl= sl1+1,ns
          l2p= Memb(sl)%l;  l2n= Memb(sl)%m
          If (l2p < l1p .or. l2n < l1n) Cycle

          If (All(cand(sl)%L1(:l1p,l2p) == cand(s0)%L1(:l1p,l1p)) .and. &
              All(cand(sl)%L2(:l1n,l2n) == cand(s0)%L2(:l1n,l1n))) then
            cls_mean(0,s0) = cls_mean(0,s0)  + Memb(sl)%fsx
            cls_mean(1:,s0)= cls_mean(1:,s0) + Memb(sl)%rx
            cls_mean2(:,s0)= cls_mean2(:,s0) + Memb(sl)%ux
          End if
        End do
      End do
    End do Level_loop
    
    Where (cls_mean2 > 0) cls_mean(1:,:)= cls_mean(1:,:) / cls_mean2

!   Sort the candidate clusters by weight

    Call Sort (.false., cls_mean(0,:), key)
    Forall(c=1:nc) cls_mean(c,:)= cls_mean(c,key)

    If (pr_out > 1) then
      Call Out ("Candidate cluster sizes", cls_mean(0,:))
      Call Out (-1,"Candidate cluster mean vectors", cls_mean(1:,:))
    End if
  End Subroutine Candidate_clusters2


   Subroutine Select_init1 (np,nr, cls_mean, Mean_rnd)
   
!    Select 'nt' random subsets 'Mean_rnd' of size 'nr' of the 
!    candidate clusters - those which top rate 1, or, 3 candidates. 
   
!    Except make the first subset of 'Mean_rnd' the top (= first)
!    'nr' candidate clusters from 'cls_mean'
   
!    Within each random subset, order the mean rating vectors 
!    by decreasing weight
   
     Integer, Intent(in) :: np               ! # candidates to be elected = total cluster set weight
     Integer, Intent(in) :: nr               ! # random mean vectors selected per set
     Real,    Intent(in) :: cls_mean(0:,:)   ! (0:nc,nl) Candidate clusters
                                             ! (1:,:) = mean vectors,  
                                             !   sorted by decreasing weight (0,:) 
     Real,   Intent(out) :: Mean_rnd(0:,:,:) ! (0:nc,nr,nt) Mean rating vectors for random initial clusters,
                                             !   ordered by decreasing cluster size, with sizes summing to 'np'
!  Local:
     Integer :: Selected(nr), key(nr)
     Real    :: wt, wt_fac, Urnd
     Integer :: i, j, cl, cn, nl, nt, Irnd

!     Integer :: iseed(2);  Call Random_seed (Get = iseed)

     If (pr_out >= 1) Call Out ("Enter 'Select_init1' to select random subsets of clusters")
     
     nt= Size(Mean_rnd,3); nl= Size(cls_mean,2);  Mean_rnd= 0
     
     Cluster_set : Do i= 1,nt
       Selected= 0  
       Clusters : Do cl= 1,nr
         Random_loop : Do
           Call Random_number (Urnd);  Irnd= Ceiling(nl*Urnd)
           
           j= First_true(Selected(:cl-1) == Irnd)
           If (j > 0) Cycle Random_loop  ! Repetition
           
           Selected(cl)= Irnd;  Cycle Clusters
         End do Random_loop
       End do Clusters
       
       Call Sort (.true.,Selected)
       Mean_rnd(:,:,i)= cls_mean(:,Selected)

       wt= Sum(Mean_rnd(0,:,i));  wt_fac= np / wt
       Mean_rnd(0,:,i)= wt_fac * Mean_rnd(0,:,i)

       Call Sort (.false.,Mean_rnd(0,:,i),key)
       Mean_rnd(1:,:,i)= Mean_rnd(1:,key,i)
       
       If (pr_out > 1) then
         Call Out ("Randomly selected candidate clusters",Selected)
         Call Out ("Cluster weights",Mean_rnd(0,:,i))
         Call Out (-1,"Cluster mean vectors",Mean_rnd(1:,:,i))
       End if
     End do Cluster_set
     
   End Subroutine Select_init1
 
   Subroutine Consolidate_strong (mr, np,nb, wtb,ballot, sl_wt,sl_lv, tot_wt,nsl,consol)
                                     
!    Consolidate the input ballots to slate ballots using the minimum 
!    weight criterion Min_slate. If a ballot 'b' ranking 'n' candidates 
!    does not already satisfy this weight, combine it with a ballot
!    ranking only n-1 candidates that match b's top n-1 candidates. 

!    Keep combining ballots this way until the combined ballot's weight 
!    exceeds Min_slate. When this gets down to a single candidate, 
!    n = 1, then any remaining combined ballots are designated 
!    as slate ballots regardless of weight.

!    However, start by identifying all initial ballots which already
!    qualify as slate ballots because their weight exceeds Min_slate
  
     Integer,  Intent(in) :: mr            ! Max # candidates to be ranked
     Integer,  Intent(in) :: np            ! # to be elected
     Integer,  Intent(in) :: nb            ! # input ballots

     Real,     Intent(in) :: wtb(:)        ! (nb) Initial ballot weights summing to 'np'
     Integer,  Intent(in) :: ballot(0:,:)  ! (0:mr,nb) Ranking data 
                                           !   (0,b)   = # candidates ranked
                                           !   (1:n,b) = candidates in preferential order
     Real,    Intent(out) :: sl_wt(:)      ! (msl=>nsl) ! Weights of the final slate ballots
     Integer, Intent(out) :: sl_lv(:)      ! (msl=>nsl) ! Level 'lv' (# ranked) 
                                           !              by the foundational ballot
       
     Real,    Intent(out) :: tot_wt        ! total weight of all slate clusters
     Integer, Intent(out) :: nsl           ! # slate ballots identified
     Integer, Intent(out) :: consol(:)     ! (nb) Maps initial ballots to final slate ballots
! Local:
     Real, Parameter :: eps= 0.001
     Real    :: wt(nb)       ! Weights for possible slate ballots 
     Integer :: consol0(nb)  ! Local mapping of initial ballots to provisional slate ballots
     Integer :: slt(mr,nb)   ! Candidates in increasing order for possible slate ballots

     Integer :: cn(mr)       ! Local candidates in increasing order
     Integer :: dt(0:mr,2)   ! Records 'nsl' (cumulative) and 'bk' (local)

     Real    :: tot_wt0
     Integer :: b, i, j, bk, lv, n0, mlv, ns0

     Call Out ("Enter 'Consolidate_strong': Combine strongly ranked initial ballots to slate ballots")

     nsl= 0;  consol= 0;  sl_wt= 0;  sl_lv= 0;  tot_wt0= Sum(wtb)

!    First add all individual initial ballots of sufficient weight

     Ballot_loop0 : Do b= 1,nb
       If (wtb(b) > Min_slate) then
         nsl= nsl + 1;  consol(b)= nsl
         sl_wt(nsl)= wtb(b);  sl_lv(nsl)= ballot(0,b)
       End if
     End do Ballot_loop0

     If (nsl == nb) Return

!    Next add combined initial ballots of sufficient weight,
!    Starting with those of which rank 'mlv' or fewer candidates

     j= 0;  dt= 0;  dt(0,1)= nsl

!    Combine ordered initial ballots

     Lev_loop : Do lv= mr-1,1,-1

       cn= 0;  wt= 0;  slt= 0;  bk= 0;  consol0= 0

       Ballot_loop1 : Do b= 1,nb
         If (consol(b) > 0 .or. ballot(0,b) < lv) Cycle Ballot_loop1

         cn(:lv)= ballot(1:lv,b)

         Do i= 1,bk  ! Try to match prior ballots at this level
           If (All(cn(:lv) == slt(:lv,i))) then
             wt(i)= wt(i) + wtb(b)   ! Add up weight
             consol0(b)= i           ! Record this membership
             Cycle Ballot_loop1
           End if
         End do

!        New slate ballot

         bk= bk + 1;  consol0(b)= bk
         wt(bk)= wtb(b);  slt(:lv,bk)= cn(:lv)
       End do Ballot_loop1

!      Record new the possible slate ballots as final slate ballots
!      if they are of sufficient weight

       Do i= 1,bk
         If (lv == 1 .or. wt(i) > Min_slate) then
           nsl= nsl + 1;  sl_wt(nsl)= wt(i);  sl_lv(nsl)= lv
           Where (consol0 == i) consol= nsl
         End if
       End do

       j= j + 1;  dt(j,1)= nsl;  dt(j,2)= bk
     End do Lev_loop

     If (pr_out >= 1) Call Out (-1, "Cumulative # slate ballots + local # provisional ballots",dt)

     n0= Count(consol <= 0)
     If (n0 > 0) then
       Call Out ("Error in 'Consolidate_strong': # initial ballots not mapped",n0, ln=1)
       Stop
     End if

     tot_wt= Sum(sl_wt(:nsl))
     If (Abs(tot_wt - tot_wt0) > eps) then
       Call Out ("Error in 'Consolidate_strong'")
       Call Out ("Weight sum of slate clusters",tot_wt, "versus prior",tot_wt0)
       Stop
     End if

   End Subroutine Consolidate_strong

   Subroutine Consolidate_weak (mr, np,nb, wtb,ballot,ballot2, &
                                sl_wt,sl_lv, tot_wt,nsl,consol)
                                     
!    Consolidate the input ballots to slate ballots using the minimum 
!    weight criterion Min_slate. If a ballot 'b' ranking 'n' candidates 
!    does not already satisfy this weight, combine it with a ballot
!    ranking only n-1 candidates that match b's top n-1 candidates. 

!    Keep combining ballots this way until the combined ballot's weight 
!    exceeds Min_slate. When this gets down to a single candidate, 
!    n = 1, then any remaining combined ballots are designated 
!    as slate ballots regardless of weight.

!    However, start by identifying all initial ballots which already
!    qualify as slate ballots because their weight exceeds Min_slate
  
     Integer,  Intent(in) :: mr            ! max # candidates ranked 
     Integer,  Intent(in) :: np            ! # candidates to be elected = total ballot weight
     Integer,  Intent(in) :: nb            ! # initial ballots

     Real,     Intent(in) :: wtb(:)        ! (nb) Initial ballot weights summing to 'np'
     Integer,  Intent(in) :: ballot(0:,:)  ! (0:mr,nb) For each consolidated ballot: (0,b) = # ranked 
                                           !   or rated, (1:mr,b) = candidates in preferential order
     Integer,  Intent(in) :: ballot2(:,:)  ! (mr,nb) (1:mr,b) = corresponding rating levels
     Real,    Intent(out) :: sl_wt(:)      ! (msl=>nsl) ! Weights of the slate ballots
     Integer, Intent(out) :: sl_lv(:)      ! (msl=>nsl) ! Their foundational levels (# ranked)
       
     Real,    Intent(out) :: tot_wt        ! total weight of all slate clusters
     Integer, Intent(out) :: nsl           ! # slate ballots identified
     Integer, Intent(out) :: consol(:)     ! (nb) Maps initial ballots to final slate ballots
! Local:
     Real, Parameter :: eps= 0.001
     Real    :: wt(nb)        ! Weights for possible slate ballots 
     Integer :: consol0(nb)   ! Local mapping of initial ballots to provisional slate ballots
     Integer :: slt(mr,2,nb)  ! Candidates in increasing order for possible slate ballots (:,1,b)
                              !   and corresponding ranking levels (:,2,b)

     Integer :: cn(mr,2)      ! Local candidates & corresponding ranking levels
     Integer :: key(mr)       ! Ordering key for cn(:,1)
     Integer :: dt(0:mr,2)    ! Records 'nsl' (cumulative) and 'bk' (local)

     Real    :: tot_wt0
     Integer :: b, i, j, n, bk, c1, lv, n0, mlv

     Call Out ("Enter 'Consolidate_weak': Combine weakly ranked initial ballots to slate ballots")

!    First add all individual initial ballots of sufficient weight

     nsl= 0;  consol= 0;  sl_wt= 0;  sl_lv= 0;  tot_wt0= Sum(wtb)

     Do b= 1,nb
       If (wtb(b) <= Min_slate) Cycle
       nsl= nsl + 1;  consol(b)= nsl
       sl_wt(nsl)= wtb(b);  sl_lv(nsl)= ballot(0,b)
     End do

     If (nsl == nb) Return

!    Next add combined initial ballots of sufficient weight,
!    Starting with those of which rank 'mlv' or more candidates

     j= 0;  dt= 0;  dt(0,1)= nsl

!    Combine unordered initial ballots

     Lev_loop : Do lv= mr-1,1,-1
       cn= 0;  key= 0;  wt= 0;  slt= 0;  bk= 0;  consol0= 0

       Ballot_loop2 : Do b= 1,nb
         If (consol(b) > 0 .or. ballot(0,b) < lv) Cycle Ballot_loop2

         If (lv == 1) then
           cn(1,1)= ballot(1,b);  cn(1,2)= ballot2(1,b)
         Else
           cn(:lv,1)= ballot(1:lv,b)
           Call Sort(.true.,cn(:lv,1), key(:lv));  cn(:lv,2)= ballot2(key(:lv),b)
         End if

         Do i= 1,bk  ! Try to match prior ballots at this level
           If (All(cn(:lv,:) == slt(:lv,:,i))) then
             wt(i)= wt(i) + wtb(b)   ! Add up weight
             consol0(b)= i           ! Record this membership
             Cycle Ballot_loop2
           End if
         End do

!        New slate ballot

         bk= bk + 1;  consol0(b)= bk
         wt(bk)= wtb(b);  slt(:lv,:,bk)= cn(:lv,:)
       End do Ballot_loop2

!      Record new the possible slate ballots as final slate ballots
!      if they are of sufficient weight or at the last level: lv = 1

       Do i= 1,bk
         If (lv == 1 .or. wt(i) > Min_slate) then
           nsl= nsl + 1;  sl_wt(nsl)= wt(i);  sl_lv(nsl)= lv
           Where (consol0 == i) consol= nsl
         End if
       End do

       j= j + 1;  dt(j,1)= nsl;  dt(j,2)= bk
     End do Lev_loop

     If (pr_out >= 1) Call Out (-1, "Cumulative # slate ballots + local # provisional ballots",dt)

     n0= Count(consol <= 0)
     If (n0 > 0) then
       Call Out ("Error in 'Consolidate_weak': # initial ballots not mapped",n0, ln=1)
       Stop
     End if

     tot_wt= Sum(sl_wt(:nsl))
     If (Abs(tot_wt - tot_wt0) > eps) then
       Call Out ("Error in 'Consolidate_weak'")
       Call Out ("Weight sum of slate clusters",tot_wt, "versus prior",tot_wt0)
       Stop
     End if
   End Subroutine Consolidate_weak

   Subroutine Reord_rankings (nsl, sl_wt,sl_lv, consol, Sl_nlv)

!    Reorder the slate ballots:  First by level (foundational # ranked)
!    then by ballot weights within levels.


     Integer,    Intent(in) :: nsl       ! # slate ballots
     Real,    Intent(inout) :: sl_wt(:)  ! (nsl) Weights of the reordered slate ballots
     Integer, Intent(inout) :: sl_lv(:)  ! (nsl) Foundational levels of the reordered slate ballots

     Integer, Intent(inout) :: consol(:) ! (nb) Maps initial ballots to reordered slate ballots

     Integer,   Intent(out) :: Sl_nlv(:) ! (mr) # slate ballots through each sum level
! Local:
     Integer :: key(nsl), key1(nsl), inv(nsl)
     Integer :: lv, sl

     If (pr_out >= 1) Call Out ("Enter 'Reord_rankings' to reorder slate ballots")

     Call Sort (.false., sl_wt, key)
     sl_lv= sl_lv(key)

     Call Sort (.true., sl_lv, key1)
     sl_wt= sl_wt(key1)

     key= key(key1);  Forall(sl=1:nsl) inv(key(sl))= sl
     consol= inv(consol)

!    The slate ballots are ordered by level. Identify the slate ballots
!    where the level increases.

     Do lv= 1,Size(Sl_nlv)
       Sl_nlv(lv)= Last_true(sl_lv <= lv)
     End do

     If (pr_out >= 1) then
       Call Out ("Cumulative # slate ballots <= each level",Sl_nlv)
       Call Out ("Slate ballot weights",sl_wt)
       Call Out ("Foundational # ranked candidates",sl_lv)
       If (pr_out > 1) Call Out ("The initial to slate ballot mapping",consol)
     End if
   End Subroutine Reord_rankings


   Subroutine Consolidate_ratings (m1,m2, np,nb, wtb,ballot,ballot2, &
                                   sl_wt,sl_lv, tot_wt,nsl,consol)
                                     
!    Consolidate the input ballots to slate ballots using the minimum 
!    weight criterion Min_slate. If a ballot 'b' rating 'l1' positive
!    candidates and 'l2' negative candidates does not already satisfy 
!    this weight, combine it with ballots with fewer positive or 
!    negative candidates which have the same positive and negative
!    ratings for these fewer numbers

!    Keep combining ballots this way until the combined ballot's weight 
!    exceeds Min_slate. When this gets down to a single positive 
!    and single negative candidate, then any remaining combined 
!    ballots are designated as slate ballots regardless of weight.

!    However, start by identifying all initial ballots which already
!    qualify as slate ballots because their weight exceeds Min_slate.
  
     Integer,  Intent(in) :: m1            ! computational limit on # positive ranked: l1 <= m1 
     Integer,  Intent(in) :: m2            ! computational limit on # negative ranked: l2 <= m2
     Integer,  Intent(in) :: np            ! # candidates to be elected = total weight
     Integer,  Intent(in) :: nb            ! # input ballots
                                           !   with l1 + l2 <= mr 
     Real,     Intent(in) :: wtb(:)        ! (nb) Input ballot weights
     Integer,  Intent(in) :: ballot(0:,:)  ! (0:mr,nb) For each consolidated ballot: (0,b) = # ranked 
                                           !   or rated, (1:mr,b) = candidates in preferential order
     Integer,  Intent(in) :: ballot2(0:,:) ! (0:mr,nb) (1:mr,b) = corresponding rating levels
                                           !   with (0,b) = # positively rated
       
     Real,    Intent(out) :: sl_wt(:)      ! (msl=>nsl) ! Weights of the final slate ballots
     Integer, Intent(out) :: sl_lv(0:,:)   ! (0:2,msl=>nsl) Levels of of these slates
                                           ! 0 = l1 + l2 = sum level
                                           ! 1 = l1 = # pos ratings   
                                           ! 2 = l1 = # negative ratings   
     Real,    Intent(out) :: tot_wt        ! total weight of all slate clusters
     Integer, Intent(out) :: nsl           ! # slate ballots identified
     Integer, Intent(out) :: consol(:)     ! (nb) Maps initial ballots to the slate ballots
! Local:
     Real, Parameter :: eps= 0.001
     Integer :: slt(m1+m2,2,nb)   ! (:,1,i) = candidates for possible slate ballot 'i'
                                  ! (:,2,i) = corresponding rating levels      
     Integer :: cn(m1+m2,2)       ! (:,1) = local candidates
                                  ! (:,2) = corresponding rating levels      
     Integer :: ky(m1+m2)         ! Ordering key for the candidates cn(:,1)
     Integer :: dt(0:m1*m2+1,2)   ! (0:ml,2) Records 'nsl' (cumulative) and 'bk' (local)

     Real    :: wt(nb)            ! Weights of possible slate ballots 
     Integer :: consol0(nb)       ! Intermediate mapping from initial to possible slate ballots
     Real    :: tot_wt0
     Integer :: b, i, j, n, bk, l1, l2, ln, lv, ml, mr, n0, n1, n2, nq, nr, nv, msl

     Call Out ("Enter 'Consolidate_ratings': Combine rating initial ballots to slate ballots")

     mr= Ubound(ballot,1);  msl= Size(sl_wt)
     nsl= 0;  consol= 0;  sl_wt= 0;  sl_lv= 0;  tot_wt0= Sum(wtb)

!    First add all individual initial ballots of sufficient weight

     Ballot_loop0 : Do b= 1,nb
       If (wtb(b) > Min_slate) then
         nsl= nsl + 1;  consol(b)= nsl;  sl_wt(nsl)= wtb(b)
         nr= ballot(0,b);  n1= ballot2(0,b);  n2= nr - n1  
         sl_lv(0,nsl)= nr;  sl_lv(1,nsl)= n1;  sl_lv(2,nsl)= n2
       End if
     End do Ballot_loop0

     If (nsl == nb) Return

     j= 0;  dt= 0;  dt(0,1)= nsl

     Neg_loop : Do l2= m2,1,-1
       Pos_loop : Do l1= m1,1,-1

         If (l1 == m1 .and. l2 == m2) Cycle Pos_loop

         wt= 0;  bk= 0;  slt= 0;  consol0= 0  
         ln= l1 + 1;  lv= l1 + l2;  cn= 0;  ky= 0 

         Ballot_loop : Do b= 1,nb
           nr= ballot(0,b);  n1= ballot2(0,b);  n2= nr - n1
           If (consol(b) > 0 .or. n1 < l1 .or. n2 < l2) Cycle Ballot_loop

           If (l1 > 1) then
             cn(:l1,1)= ballot(1:l1,b);  Call Sort(.true.,cn(:l1,1), ky(:l1))
             cn(:l1,2)= ballot2(ky(:l1),b)
           Else
             cn(1,1)= ballot(1,b);  cn(1,2)= ballot2(1,b)
           End if
          
           nq= n1 + 1;  nv= n1 + l2
           If (l2 > 1) then
             cn(ln:lv,1)= ballot(nq:nv,b);  Call Sort(.true.,cn(ln:lv,1), ky(:l2))
             cn(ln:lv,2)= ballot2(n1+ky(:l2),b)
           Else
             cn(ln,1)= ballot(nq,b);  cn(ln,2)= ballot2(nq,b)
           End if

           Do i= 1,bk  ! Try to match prior ballots at this level
             If (All(cn(:lv,:) == slt(:lv,:,i))) then
               wt(i)= wt(i) + wtb(b)   ! Add up weight
               consol0(b)= i           ! Record this membership
               Cycle Ballot_loop
             End if
           End do

!          New slate ballot

           bk= bk + 1;  wt(bk)= wtb(b);  consol0(b)= bk
           slt(:lv,:,bk)= cn(:lv,:)
         End do Ballot_loop

!        Record the new possible slate ballots as final slate ballots
!        if they are of sufficient weight or at the last level: l1 = l2 = 1

         Do i= 1,bk
           If (wt(i) > Min_slate) then
             nsl= nsl + 1;  sl_wt(nsl)= wt(i)
             sl_lv(0,nsl)= lv;  sl_lv(1,nsl)= l1;  sl_lv(2,nsl)= l2 
             Where (consol0 == i) consol= nsl
           End if
         End do

         j= j + 1;  dt(j,1)= nsl;  dt(j,2)= bk

!        Record the new possible slate ballots as final slate ballots, based on 
!        candidate matching only, if they are at the last level: l1 = l2 = 1

         If (l1 == 1 .and. l2 == 1) then
           wt= 0;  slt= 0;  bk= 0;  consol0= 0

           Ballot_loop1 : Do b= 1,nb
             If (consol(b) > 0) Cycle Ballot_loop1

             nr= ballot(0,b);  n1= ballot2(0,b);  nq= n1 + 1

             Do i= 1,bk  ! Try to match prior ballots at this level
               If (ballot(1,b) == slt(1,1,i) .and. ballot(nq,b) == slt(2,1,i)) then
                 wt(i)= wt(i) + wtb(b)   ! Add up weight
                 consol0(b)= i           ! Record this membership
                 Cycle Ballot_loop1
               End if
             End do

!            New slate ballot

             bk= bk + 1;  wt(bk)= wtb(b);  consol0(b)= bk
             slt(1,1,bk)= ballot(1,b);   slt(2,1,bk)= ballot(nq,b)
           End do Ballot_loop1

           n= nsl + bk
           If (n > msl) then
             Call Out ("Out of storage in 'Consolidate_ratings' beyond #",nsl,ln=1)
             Call Out ("since adding",bk, "exceeds the limit by",n - msl)
             Stop
           End if

           Do i= 1,bk
             nsl= nsl + 1;  sl_wt(nsl)= wt(i)
             sl_lv(0,nsl)= lv;  sl_lv(1,nsl)= l1;  sl_lv(2,nsl)= l2 
             Where (consol0 == i) consol= nsl
           End do

           j= j + 1;  dt(j,1)= nsl;  dt(j,2)= bk
         End if

       End do Pos_loop
     End do Neg_loop

     If (pr_out >= 1) Call Out (-1, "Cumulative # slate ballots + local # provisional ballots",dt)

     n0= Count(consol <= 0)
     If (n0 > 0) then
       Call Out ("Error in 'Consolidate_ratings': # initial ballots not mapped",n0, ln=1)
       Stop
     End if

     tot_wt= Sum(sl_wt(:nsl))
     If (Abs(tot_wt - tot_wt0) > eps) then
       Call Out ("Error in 'Consolidate_ratings'")
       Call Out ("Weight sum of slate clusters",tot_wt, "versus prior",tot_wt0)
       Stop
     End if

   End Subroutine Consolidate_ratings

   Subroutine Reord_ratings (nsl, sl_wt,sl_lv, consol, Sl_nlv)

!    Reorder the rating slate ballots, first sum level, 
!    then by ballot weight within levels.                                     
       
     Integer,    Intent(in) :: nsl          ! # slate ballots
     Real,    Intent(inout) :: sl_wt(:)     ! (nsl) ! Weights of the slate ballots
     Integer, Intent(inout) :: sl_lv(0:,:)  ! (0:2,nsl) Levels of of these slates
                                            ! 0 = l1 + l2 = sum level
                                            ! 1 = l1 = # pos ratings   
                                            ! 2 = l2 = # negative ratings   
     Integer, Intent(inout) :: consol(:)    ! (nb) Maps initial ballots to reordered slate ballots
     Integer,   Intent(out) :: Sl_nlv(:)    ! (mlv) # slate ballots through each sum level
! Local:
     Integer :: key(nsl), key1(nsl), inv(nsl)
     Integer :: sl, sm

     If (pr_out >= 1) Call Out ("Enter 'Reord_ratings' to consolidate rating ballots to slate ballots")

     Call Sort (.false., sl_wt, key)
     sl_lv(0,:)= sl_lv(0,key)

     Call Sort (.true., sl_lv(0,:), key1)
     sl_wt= sl_wt(key1)
  
     key= key(key1)
     sl_lv(1,:)= sl_lv(1,key);  sl_lv(2,:)= sl_lv(2,key)

     Forall(sl=1:nsl) inv(key(sl))= sl;  consol= inv(consol)

     Sl_nlv(1)= 0
     Do sm= 2,Size(Sl_nlv)
       Sl_nlv(sm)= Last_true(sl_lv(0,:) <= sm)
     End do

     If (pr_out >= 1) then
       Call Out ("Cumulative # slate ballots <= each level",Sl_nlv)
       Call Out ("Slate ballot weights",sl_wt)
       Call Out (-1,"Foundational # sum, pos, neg candidates",sl_lv)
     End if
   End Subroutine Reord_ratings


   Subroutine Consol1_ballots (nc,np,mr,nb, wtb,ballot, nr,nsl,Memb, Mean_rnd)
   
!    Consolidate ballots to slate ballots / slate clusters, assuming strong ranking ballots.
!    The membership data for the slate clusters are recorded in 'Memb' for subsequent use
!    in the clustering algorithm. In particular, the final clusters are consist of slate 
!    cluster members, though these slate clusters may have partial membership in more than 
!    one cluster of a cluster set (="fuzzy set" membership).

     Integer,   Intent(in) :: nc, nb        ! # candidates and # initial ballots
     Integer,   Intent(in) :: np, mr        ! # candidates to be elected and max # ranked

     
     Real,      Intent(in) :: wtb(:)        ! (nb) Input ballot weights
     Integer,   Intent(in) :: ballot(0:,:)  ! (0:mr,nb) For each consolidated ballot: (0,b) = # ranked 
                                            !   or rated, (1:mr,b) = candidates in preferential order
     Integer,  Intent(out) :: nr            ! 
     Integer,  Intent(out) :: nsl           ! Final # slate clusters

     Type(Multi_listD), Pointer :: Memb(:)  ! (0:ns) Final slate cluster data

     Real,  Intent(out) :: Mean_rnd(0:,:,:) ! (0:nc,nr,nt) 'nr' mean rating vectors for
                                            !   'nt' random sets of clusters, 
                                            !   with cluster weights at 0
     
! Local:
     Real,    Allocatable :: rates(:,:)    ! (0:nc,nb) Ballot weights at 0, 
                                           !   mean rating vectors at 1:nc
     Real,    Allocatable :: cls_mean(:,:) ! (0:nc,lim) Mean rating vectors for the
                                           !   candidate clusters, cluster weight at 0
     Real,    Allocatable :: sl_wt(:)      ! (msl)     Weights of slate ballots
     Integer, Allocatable :: sl_lv(:,:)    ! (0:2,msl) Levels, etc, of slate ballots

     Integer, Allocatable :: Sl_nlv(:)     ! (mlv) Sum level counter: the last slate ballot 
                                           !       for each sum level
     Integer :: consol(nb)    ! Maps initial ballots to slate ballots
     Real    :: tot_wt
     Integer :: mlv, m1, m2, mt, nl, ns, nt, lim, msl
     
     Call Out ("")
     Call Out ("Enter 'Consol1_ballots': initial strong ranking ballots to slate clusters")
     
     nr= Size(Mean_rnd,2);  nt= Size(Mean_rnd,3)

     msl= Min(Max_sl, nb)
     Min_slate= np / Real(msl)
     
     Allocate(sl_wt(msl), sl_lv(0:2,msl))
     
     Call Consolidate_strong (mr, np,nb, wtb,ballot, sl_wt,sl_lv(0,:), tot_wt,nsl,consol)

     mlv= Maxval(sl_lv(0,:));  Allocate(Sl_nlv(mlv))
     Call Reord_rankings (nsl, sl_wt(:nsl),sl_lv(0,:nsl), consol, Sl_nlv)

     sl_lv(1,:nsl)= sl_lv(0,:nsl);  sl_lv(2,:nsl)= 0

     Call Out(" # standard run slate clusters",nsl, ln=1)
     
!      Compute the slate cluster membership data 'Memb'
       
     Allocate (Memb(0:nsl), Memb(0)%qx(2), Memb(0)%vl(mlv))
     Memb(0)%l= np;  Memb(0)%n= nc;  Memb(0)%fsx= tot_wt
     Memb(0)%qx= Parm_cor;  Memb(0)%vl= Sl_nlv
       
     Memb(1:)%fsx= sl_wt(:nsl)
     Memb(1:)%l= sl_lv(1,:nsl)
     Memb(1:)%m= sl_lv(2,:nsl)

     Call Slate_clust1 (nc,mr, nb,nsl, wtb,ballot, consol, Memb)

!    Compute the clusters to be used to form the initial cluster sets 
    
     Allocate(cls_mean(0:nc,nsl))
     Call Candidate_clusters (nc,nsl, Memb, cls_mean)
       
     nr= Min(nr,nsl-1)
     Call Select_init1 (np,nr, cls_mean, Mean_rnd(:,:nr,:))

   End Subroutine Consol1_ballots

                          
   Subroutine Slate_clust1 (nc,mr,nb,ns, wtb,ballot, consol, Memb)
                                     
!    Compute the slate ballot data structure 'Memb', computing 
!    the maximum liklihood estimates of the mean vectors of the 
!    slate clusters computed from the slate ballots.
  
     Integer,  Intent(in) :: nc             ! # candidates
     Integer,  Intent(in) :: mr             ! max # candidates ranked or rated

     Integer,  Intent(in) :: nb             ! # ballots
     Integer,  Intent(in) :: ns             ! # slate ballots

     Real,     Intent(in) :: wtb(:)         ! (nb) Ballot weights summing to 'np'
     Integer,  Intent(in) :: ballot(0:,:)   ! (0:mr,nb)  (1:nr,b) = candidates in preferential order
                                            !   (0,b) = 'nr' = # ranked or rated  

     Integer,  Intent(in) :: consol(:)      ! (nb Maps initial ballots to slate ballots

     Type(Multi_listD), Intent(inout) :: Memb(0:) ! (0:ns) Final slate cluster data
! Local:
     Real, Parameter :: eps= 1.0E-12
     Real(Dblp) :: rt1, cor, dif, rate(0:nc,ns), sm(nc), sm1(nc), sm2(nc), wt_var(nc)
     Real    :: mean(nc), var_pt(nc)
     Real    :: rank_pt(mr,2), unrank_pt(mr,2)

     Integer :: lst(nb), key(ns)
     Integer :: b, j, n, lm, ln, lv, n1, n2, nr, s1, s2, sl, neg
     
     Call Out ("Enter 'Slate_clust1' to compute the slate ballot membership data structure")

     Call Ranking_pts (Parm_bal,Max_pt,nc,mr, rank_pt,unrank_pt)
     
!    Loop over all slate ballots

     dif= Max_pt
     Slate_loop : Do sl= 1,ns

       Allocate(Memb(sl)%lt(nc), Memb(sl)%px(nc), Memb(sl)%rx(nc), &
                Memb(sl)%sx(nc), Memb(sl)%tx(0:nc), Memb(sl)%ux(nc))
       Memb(sl)%lt= 0; Memb(sl)%px= 0; Memb(sl)%rx= 0
       Memb(sl)%sx= 0; Memb(sl)%tx= 0; Memb(sl)%ux= 0
       
       sm= eps;  sm1= 0;  sm2= 0
       Call List_of_true (consol == sl, ln, lst)

       Combined_ballot_loop : Do j= 1,ln  ! Loop over the matching consolidated ballots
         b= lst(j);  nr= ballot(0,b)

         mean= unrank_pt(nr,1);  var_pt= unrank_pt(nr,2)
         mean(ballot(1:nr,b))  = rank_pt(:nr,1)
         var_pt(ballot(1:nr,b))= unrank_pt(:nr,2)

         wt_var= wtb(b) / var_pt
         sm = sm + wt_var
         sm1= sm1 + mean * wt_var
         sm2= sm2 + (var_pt + mean**2) * wt_var
       End do Combined_ballot_loop
       
!      Compute the mean vector of the slate ballot and its variance
       sm1= sm1 / sm;  sm2= sm2 / sm
       Memb(sl)%sx= sm1
       Memb(sl)%px= sm2 - sm1**2

!      For ranking ballots, adjust the mean vectors to subtract
!      out the noise for correlation computations, using %tx

       rate(1:,sl)= sm1 - Noise_cor

       rt1= Maxval(rate(1:,sl))
       If (rt1 < Min_mx) then
         Call Out ("Warning in 'Slate_clust': Low noise adjusted rating",rt1, ln=1)
         Call Out ("at slate",sl)
         rate(1:,sl)= rate(1:,sl) + (Min_mx - rt1)
       End if
     End do Slate_loop

!    Correct the noise significance level 'Noise_cor', if necessary
!    (unlikely) so that no slate cluster is left 'in the noise'

     Call Normalize_vec (rate)

     Do sl= 1,ns 
       Memb(sl)%tx= rate(:,sl)
       Memb(sl)%ux= Memb(sl)%fsx / Memb(sl)%px
       Memb(sl)%rx= Memb(sl)%sx * Memb(sl)%ux
     End do

!    List the slate ballots positively correlated with each slate ballot
!    and subsequent to it, plus the correlation, ordered by decreasing 
!    correlation

     lst= 0;  Memb(1:)%k= 0
     Do s1= 1,ns-1
       j= 0
       Do s2= s1+1,ns 
         Call Dot_product_M (Dot_fac, Memb(s1)%tx(1:),Memb(s2)%tx(1:), neg, cor)

         If (cor > Memb(0)%qx(1)) then
           j= j + 1;  lst(j)= s2;  rate(0,j)= cor
         End if
       End do

       If (j > 0) then
         Memb(s1)%k= j;  Allocate(Memb(s1)%vl(j), Memb(s1)%qx(j))
         Memb(s1)%qx= rate(0,:j)
         Call Sort (.false., Memb(s1)%qx, key(:j))
         Memb(s1)%vl= lst(key(:j))
       End if
     End do

     Call Form_Memb_sub (nc,ns, Memb)
     
   End Subroutine Slate_clust1

End Module Clusters2