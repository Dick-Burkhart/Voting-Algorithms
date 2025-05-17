 Module Misc_methods

   Use Graph_algorithms
   Use Util
   Use Output
   Use Types
   Use Precisn
   Implicit None

 Contains

   Subroutine Non_clustering (idist,Method,Rating, np,mr,nc,nb, pt_val, &
                              wtb,ballot,ballot2,elected, STVmb,DTVmb)
   
!    Compute the elected sets by several non clustering methods:
!    STV, DTV, IRV, Top Vote, Borda Points, Borda, Condorcet, Approval,  
!    Borda Elimination, Majority Approval

!    The algorithms are constructed with options to minimize memory requirements.
   
     Integer,      Intent(in) :: idist          ! Voting district index
     Logical,      Intent(in) :: Method(:)      ! (nm) The methods to compute, nm <= 10
     Integer,      Intent(in) :: Rating         ! 1 for rating data
                                                ! 0 for ranking data with possible equal rankings
                                                ! -1 for candidates in ranking order
     Integer,      Intent(in) :: np             ! # to be elected
     Integer,      Intent(in) :: mr             ! max # ranked or rated candidates
     Integer,      Intent(in) :: nc             ! total # candidates
     Integer,      Intent(in) :: nb             ! total # ballots

     Real,         Intent(in) :: pt_val(:)      ! (mt) Decreasing point values of the rating levels
     
     Real,         Intent(in) :: wtb(:)         ! (nb) Ballot weights
     Integer,      Intent(in) :: ballot(0:,:)   ! (0:mr,nb). (0,:) = 'p' candidates in preferential order
     Integer,      Intent(in) :: ballot2(0:,:)  ! (0:mr,nb). corresponding ranking or rating levels
                                                !   with (0,:) = 'n' positively rated candidates

     Integer,     Intent(out) :: elected(:,:,:) ! (nc,nMt,2) Candidates in preferential order (1)
                                                !   or 1..np electeds in increasing order (2), for the 'nMt' 
                                                !   non-clustering methods, with 1 = STV & 2 = DTV

     Type(Set_list), Intent(out) :: STVmb(0:)   ! (0:nb) Ballot structure for clustering
                                                ! Ballot sets (0): 
                                                !   %n        = nc
                                                !   %p        = np
                                                !   %svl      = elected vote = Sum(%val(set(:np)))
                                                !   %smb      = independents vote = Sum(%val(set(np+1:)))
                                                !   %set(nc)  = candidates in elected ordering, with
                                                !               the first 'np' elected
                                                !   %lev(nc)  = candidate clusters by decreasing cluster vote
                                                !   %val(nc)  = candidate cluster votes, not reordered

                                                ! Ballots (1:nb):
                                                !   %n     = # elected candidates
                                                !   %p     = # positively rated candidates
                                                !   %r     = # ranked candidates
                                                !   %set(r)= 'n' elected candidates, then unelected,
                                                !            in original ranking order
                                                !   %val(r)= corresponding membership weights 
                                                !   %svl   = elected membership weight= Sum(%val(:n))
                                                !   %smb   = total membership weight  = Sum(%val)
     Type(Set_list), Intent(out) :: DTVmb(0:)   ! (0:nb) Ballot structure for clustering
                                                ! Ballot sets (0): 
                                                !   %n        = nc
                                                !   %p        = np
                                                !   %svl      = elected vote = Sum(%val(set(:np)))
                                                !   %smb      = independents vote = Sum(%val(set(np+1:)))
                                                !   %set(nc)  = candidates in elected ordering, with
                                                !               the first 'np' elected
                                                !   %lev(nc)  = candidate clusters by decreasing cluster vote
                                                !   %val(nc)  = candidate cluster votes, not reordered

                                                ! Ballots (1:nb):
                                                !   %n     = # elected candidates
                                                !   %p     = # positively rated candidates
                                                !   %r     = # ranked candidates
                                                !   %set(r)= 'n' elected candidates, then unelected,
                                                !            in original ranking order
                                                !   %val(r)= corresponding membership weights 
                                                !   %svl   = elected membership weight= Sum(%val(:n))
                                                !   %smb   = total membership weight  = Sum(%val)
!  Local:
     Logical,   Parameter :: Warren= .false.    ! Use the Warren method for STV, DTV, etc., or not
     Real,      Parameter :: Max_pt= 10.0      ! Maximum ranking point value
     Real,      Parameter :: S_parm(2)= 0.0    ! Standard Borda

     Real    :: S_rank(mr)         ! Point values for the top ranked candidates,
                                   !   assuming the standard parameters 'S_parm'
     Real    :: S_unrank(mr)       ! Point values of unranked candidates,
                                   !   as a function of the # ranked, 
                                   !   assuming the standard parameters 'S_parm'

     Real    :: mean(nc,nc+1-np)   ! Mean candidate rankings in decreasing order (for Borda_elimination)
     Integer :: ord(nc,nc+1-np)    ! Corresponding candidate ordering
     
     Real, Allocatable :: ptval(:), vec(:,:)

     Real    :: total, elc_vote, siz(nc)
     Integer :: Condorcet_vote(nc,4), elim(nc-np)
     
     Integer :: i, b, m, n, im, mb, mrn, mrp, mt, n1, n2, nm, iter
     
     nm= Size(Method);  mt= Size(pt_val)

!    Rating case: For rating (Rating = 1) and weak ranking (Rating = 0) 
!    use the rating or ranking levels in 'ballot'. Also use the level
!    to point conversions specified in 'pt_val' unless the size 'mt'
!    of 'pt_val' is 1 for weak ranking, in which case use the standard 
!    ranking level to point,conversion 'S_rank', which is also used 
!    for strong ranking (Rating = -1).

     If (mt >= 4) then
       Allocate(ptval(mt)); ptval= pt_val
     End if
     
     If (Rating == 0 .and. mt < 4) then
       mt= mr;  Allocate(ptval(mt));  ptval= S_rank
     End if

     If (Rating <= 0) then
       Call Ranking_pt0 (S_parm,Max_pt,nc,mr, S_rank,S_unrank)
     End if
     
     elected= -1;  siz= 0;  mean= 0;  ord= 0
     Condorcet_vote= 0;  elim= 0;  im= 1
     
     If (Method(im)) then
       Call Out ("STV algorithm")

       If (Warren) then
         If (Rating >= 0) then
           Call WSTV2 (nc,ptval, wtb,ballot,ballot2, elected(:,im,1), STVmb)
         Else
           Call WSTV (nc, wtb,ballot, elected(:,im,1), STVmb)
         End if
       Else
         If (Rating >= 0) then
           Call STV2 (nc,ptval, wtb,ballot,ballot2, elected(:,im,1), STVmb)
         Else
           Call STV (nc, wtb,ballot, elected(:,im,1), STVmb)
         End if
       End if

       elected(:np,im,2)= elected(:np,im,1)
       Call Sort (.true.,elected(:np,im,2))

       Call Out (-1,"STV", elected(:,im,:))
       Call Out ("Clusters sizes", STVmb(0)%val)
         total= Sum(STVmb(0)%val)
       Call Out ("Elected clusters weight",STVmb(0)%svl, &
                 "vs total weight",total)
     End if
     If (nm < 2) Return;  im= 2

     If (Method(im)) then
       Call Out ("DTV algorithm")

       If (Warren) then
         If (Rating >= 0) then
           Call WDTV2 (nc,ptval, wtb,ballot,ballot2,elected(:,im,1), DTVmb)
         Else
           Call WDTV (nc, wtb,ballot, elected(:,im,1), DTVmb)
         End if
       Else
         If (Rating >= 0) then
           Call DTV2 (nc,ptval, wtb,ballot,ballot2,elected(:,im,1), DTVmb)
         Else
           Call DTV (nc, wtb,ballot, elected(:,im,1), DTVmb)
         End if
       End if

       elected(:np,im,2)= elected(:np,im,1) 
       Call Sort (.true.,elected(:np,im,2))

       Call Out (-1,"DTV", elected(:,im,:))
       Call Out ("Clusters sizes", DTVmb(0)%val)
         total= Sum(DTVmb(0)%val)
       Call Out ("Elected clusters weight",DTVmb(0)%svl, &
                 "vs total weight",total)
     End if
     If (nm < 3) Return;  im= 3
     
     If (Method(im)) then
       Call Out ("Instant runoff voting algorithm")

       If (Rating >= 0) then
         Call IRV2 (nc, ptval, wtb,ballot,ballot2, elected(:,im,1),siz)
       Else
         Call IRV (nc, wtb,ballot, elected(:,im,1),siz)
       End if

       elected(:np,im,2)= elected(:np,im,1)
       Call Sort (.true.,elected(:np,im,2))
       Call Out (-1,"IRV", elected(:,im,:))

       Call Out ("Clusters sizes", siz)
         elc_vote= Sum(siz(elected(:np,im,1)));  total= Sum(siz)
       Call Out ("Elected clusters weight",elc_vote, &
                 "vs total weight",total)
     End if
     If (nm < 4) Return;  im= 4
     
     If (Method(im)) then 
       Call Out ("Top vote algorithm")
       Call Top_vote (np,wtb,ballot, mean(:,1),ord(:,1)) 

       elected(:,im,1)= ord(:,1)   
       elected(:np,im,2)= elected(:np,im,1)
       Call Sort (.true.,elected(:np,im,2))

       Call Out (-1,"Top Vote", elected(:,im,:))
       Call Out ("Decreasing mean rankings", mean(:,1))
     End if
     If (nm < 5) Return;  im= 5

     If (Method(im)) then 
       Call Out ("Standard Borda points algorithm")

       If (Rating >= 0) then
         Call Borda_points2 (np,ptval, wtb,ballot,ballot2(1:,:), &
                             mean(:,1),ord(:,1))
       Else
         Call Borda_points (np,S_rank,S_unrank, wtb,ballot,  mean(:,1),ord(:,1))
       End if

       elected(:,im,1)= ord(:,1)   
       elected(:np,im,2)= elected(:np,im,1)
       Call Sort (.true.,elected(:np,im,2))

       Call Out (-1,"Standard Borda candidates", elected(:,im,:))
       Call Out ("Decreasing mean rankings", mean(:,1))
     End if
     If (nm < 6) Return;  im= 6
     
     If (Method(im)) then 
       Call Out ("Condorcet algorithm")
       Allocate(vec(nc,nb))

       If (Rating >= 0) then
         Call Ballot_vec2 (ptval, ballot,ballot2(1:,:), vec)
       Else
         Call Ballot_vec (S_rank,S_unrank, ballot, vec) 
       End if

       Condorcet_vote= Condorcet0(np,nc, wtb,vec)

       elected(:,im,1)= Condorcet_vote(:,1)
       elected(:np,im,2)= elected(:np,im,1)

       Call Sort (.true.,elected(:np,im,2))
       Call Out (-1,"Condorcet candidates", elected(:,im,:))
       DeAllocate(vec)
     End if
        
     If (nm < 7) Return;  im= 7
     If (Method(im)) then
       Call Out ("Approval algorithm")

       Call Approval (np,wtb,ballot, mean(:,1),ord(:,1))

       elected(:,im,1)= ord(:,1)
       elected(:np,im,2)= elected(:np,im,1)
       Call Sort (.true.,elected(:np,im,2))

       Call Out (-1,"Approval candidate cluster", elected(:,im,:))
       Call Out ("Decreasing mean # approvals", mean(:,1))
     End if
     If (nm < 8) Return;  im= 8
     
     If (Method(im)) then
       Call Out ("Borda points elimination algorithm")

       If (Rating >= 0) then
         Call Borda_elimination2 (np,nc, ptval, wtb,ballot,ballot2(1:,:), &
                                  elim,mean,ord)
       Else
         Call Borda_elimination (np,nc, S_rank,S_unrank, wtb,ballot, &
                                  elim,mean,ord)
       End if
       
       elected(:,im,1)= ord(:,nc+1-np)    
       elected(:np,im,2)= elected(:np,im,1) 
       Call Sort (.true.,elected(:np,im,2))

       Call Out (-1,"Borda Eliminiation candidates", elected(:,im,:))
       Call Out (-1,"Decreasing mean rankings by elimination", mean)
       Call Out (-1,"Corresponding candidate orderings", ord)
     End if
     If (nm < 9) Return;  im= 9
     
     If (Method(im)) then
       Call Out ("Majority approval algorithm")

       If (Rating >= 0) then
         Call Majority_approval (np,wtb,ballot, mean(:,1),ord(:,1), ballot2(1:,:))
       Else
         Call Majority_approval (np,wtb,ballot, mean(:,1),ord(:,1))
       End if

       elected(:,im,1)= ord(:,1)
       elected(:np,im,2)= elected(:np,im,1)
       Call Sort (.true.,elected(:np,im,2))

       Call Out (-1,"Majority approval candidates", elected(:,im,:))
       Call Out ("The corresponding mean # approvals", mean(:,1))
     End if
     
    End Subroutine Non_clustering
    
    Subroutine Borda_points (np,rank_pt,unrank_pt, wtb,ballot, mean,elected)
      
!     Compute alternative point sums to the Borda Count, 
!     using one of several different point allocation methods.

      Integer,  Intent(in) :: np            ! # candidates to be elected
      Real,     Intent(in) :: rank_pt(:)    ! (mr) Point values for the 'mr' top ranked candidates
      Real,     Intent(in) :: unrank_pt(:)  ! (mr) Point values for the unranked candidates, depending on # ranked
     
      Real,     Intent(in) :: wtb(:)        ! (nb) Ballot weights, summing to 'np'
      Integer,  Intent(in) :: ballot(0:,:)  ! (0:nr,nb) Candidates (1:,1,:) ranked by each ballot
                                            ! in the order ranked, with ballot(0,1,b)= # ranked 
                                            ! by ballot 'b' up to mr. (1:,2,:) = the corresponding
                                            ! ranking or rating levels. 

      Real,    Intent(out) :: mean(:)       ! (nc) Mean rankings in decreasing order.
      Integer, Intent(out) :: elected(:)    ! (nc) Corresponding candidates. 
!  Local:
      Real, Parameter :: eps= 0.00001
      Real    :: mean_pt(Size(mean))
      Integer :: b, n, nb

      nb= Ubound(ballot,2);  mean= 0;  elected= 0
     
      Do b= 1,nb
        n= ballot(0,b);  mean_pt= unrank_pt(n)
        mean_pt(ballot(1:n,b))= rank_pt(:n)

        mean= mean + mean_pt * wtb(b)
      End do
       
      mean= mean / np
      Call Sort (.false., mean, elected, ez=eps)
    End Subroutine Borda_points
    
    Subroutine Borda_points2 (np,pt_val, wtb,ballot,ballot2, mean,elected)
      
!     Compute the Borda Count, given ranking or rating levels and point conversion.

      Integer,  Intent(in) :: np            ! # candidates to be elected
      Real,     Intent(in) :: pt_val(:)     ! (mt) Decreasing point values for ratings
      Real,     Intent(in) :: wtb(:)        ! (nb) Ballot weights, summing to 'np'

      Integer,  Intent(in) :: ballot(0:,:)  ! (0:mr,nb) Candidates (1:n,:) in preferential order
      Integer,  Intent(in) :: ballot2(:,:)  ! (mr,nb) Ranking or rating levels (1:n,:) that 
                                            !         correspond to the candidates ballot(1:n,:)

      Real,    Intent(out) :: mean(:)       ! (nc) Mean rankings in decreasing order.
      Integer, Intent(out) :: elected(:)    ! (nc) Corresponding candidates. 

!  Local:
      Real, Parameter :: eps= 0.00001
      Real    :: mean_pt(Size(mean))
      Integer :: b, n, nb

      nb= Ubound(ballot,2);  mean= 0;  elected= 0
     
      Do b= 1,nb
        n= ballot(0,b);  mean_pt= 0
        mean_pt(ballot(1:n,b))= pt_val(ballot2(:n,b))
        mean= mean + mean_pt * wtb(b)
      End do
       
      mean= mean / np
      Call Sort (.false., mean, elected, ez=eps)
    End Subroutine Borda_points2
    
    
   Subroutine Borda_elimination2 (np,nc, pt_val, wtb,ballot,ballot2, elim,mean,elected)
      
!    Baldwin voting method = Borda Count with elimination. 
!    This elects the Condorcet winner if one exists, and is
!    always in the Smith set. 
!    See (http://sites.google.com/site/thefell/statisticalsnipstprelections) 
!    for other properties.

!    If mr = nc, then this is the traditional Baldwin, otherwise it's a
!    rank restricted Baldwin.

     Integer,  Intent(in) :: np             ! # candidates to be elected
     Integer,  Intent(in) :: nc             ! # candidates
     Real,     Intent(in) :: pt_val(:)      ! (mt) Decreasing point values
     Real,     Intent(in) :: wtb(:)         ! (nb) Ballot weights, summing to 'np'

     Integer,  Intent(in) :: ballot(0:,:)   ! (0:mr,nb) Candidates (1:n,:) in preferential order
     Integer,  Intent(in) :: ballot2(:,:)   ! (mr,nb) Ranking or rating levels (1:n,:) that 
                                            !         correspond to the candidates ballot(1:n,:)

     Integer, Intent(out) :: elim(:)        ! (nc-np)      Eliminated candidates, in elimination order
     Real,    Intent(out) :: mean(:,:)      ! (nc,nc-np+1) Borda mean rankings in decreasing order at each elimination
     Integer, Intent(out) :: elected(:,:)   ! (nc,nc-np+1) Corresponding candidates. Initial np candidates
                                            !              are the Borda electeds. The final np candidates
                                            !              are the Baldwin electeds.
!  Local:
     Real,      Parameter :: eps= 0.00001
     Logical :: Keep(nc)
     Real    :: sm(nc), mean_pt(nc)
     Integer :: ls(nc), lev(nc), cand(nc) 
     Integer :: b, e,  k, n, r, nb, nr, maxe
     
     nb= Ubound(ballot,2);  elim= 0;  mean= 0;  elected= 0
     
     Keep= .true.;  maxe= nc - np

     Do r= 1,maxe+1  ! Elimination round 'r'

!      Do the Borda Count and sort it into decreasing order

       sm= 0;  nr= nc - r + 1
       
       Do b= 1,nb
         n= ballot(0,b);  Call List_of_true (Keep(ballot(1:n,b)), k,ls)
         If (k <= 0) Cycle

         cand(:k)= ballot(ls(:k),b);  lev(:k)= ballot2(ls(:k),b)
         mean_pt= 0;  mean_pt(cand(:k))= pt_val(lev(:k))

         Where (Keep) sm= sm + mean_pt * wtb(b)
       End do
       
       sm= sm / np;  mean(:,r)= sm
       Call Sort (.false., mean(:,r), elected(:,r), ez=eps)
       If (r > maxe) Exit
       
!      Eliminate the last place candidate
      
       e= elected(nr,r);  elim(r)= e;  Keep(e)= .false.
     End do

     If (pr_out > 1) then
       Call Out ("Borda points elimination data")
       Call Out (-1,"Decreasing mean rankings or ratings by elimination", mean)
       Call Out (-1,"Corresponding candidate orderings", elected)
     End if

   End Subroutine Borda_elimination2

   Subroutine Borda_elimination (np,nc, rank_pt,unrank_pt, wtb,ballot, elim,mean,elected)
      
!    Baldwin voting method = Borda Count with elimination. 
!    This elects the Condorcet winner if one exists, and is
!    always in the Smith set. 
!    See (http://sites.google.com/site/thefell/statisticalsnipstprelections) 
!    for other properties.

!    If mr = nc, then this is the traditional Baldwin, otherwise it's a
!    rank restricted Baldwin.

     Integer,  Intent(in) :: np             ! # candidates to be elected
     Integer,  Intent(in) :: nc             ! # candidates
     Real,     Intent(in) :: rank_pt(:)     ! (mr) Point values for the 'mr' top ranked candidates
     Real,     Intent(in) :: unrank_pt(:)   ! (mr) Point values for the unranked candidates, depending on # ranked
     Real,     Intent(in) :: wtb(:)         ! (nb) Ballot weights, summing to 'np'
     Integer,  Intent(in) :: ballot(0:,:)   ! (0:mr,nb) Ranking data
                                            !  (0,b) = # ranked.  (1:,b) The candidates in ranking order

     Integer, Intent(out) :: elim(:)        ! (nc-np)      Eliminated candidates, in elimination order
     Real,    Intent(out) :: mean(:,:)      ! (nc,nc-np+1) Borda mean rankings in decreasing order at each elimination
     Integer, Intent(out) :: elected(:,:)   ! (nc,nc-np+1) Corresponding candidates. Initial np candidates
                                            !              are the Borda electeds. The final np candidates
                                            !              are the Baldwin electeds.
!  Local:
     Real,      Parameter :: eps= 0.00001
     Logical :: Keep(nc)
     Real    :: sm(nc), mean_pt(nc)
     Integer :: ls(nc), cand(nc) 
     Integer :: b, e, k, n, r, nb, nr, maxe
     
     nb= Ubound(ballot,2);  elim= 0;  mean= 0;  elected= 0
     
     Keep= .true.;  maxe= nc - np

     Do r= 1,maxe+1  ! Elimination round 'r'

!      Do the Borda Count and sort it into decreasing order

       sm= 0;  nr= nc - r + 1
       
       Do b= 1,nb
         n= ballot(0,b);  Call List_of_true (Keep(ballot(1:n,b)), k,ls)
         If (k <= 0) Cycle

         cand(:k)= ballot(ls(:k),b);  mean_pt= unrank_pt(k)
         mean_pt(cand(:k))= rank_pt(:k)

         Where (Keep) sm= sm + mean_pt * wtb(b)
       End do
       
       sm= sm / np;  mean(:,r)= sm
       Call Sort (.false., mean(:,r), elected(:,r), ez=eps)
       If (r > maxe) Exit
       
!      Eliminate the last place candidate
      
       e= elected(nr,r);  elim(r)= e;  Keep(e)= .false.
     End do

     If (pr_out > 1) then
       Call Out ("Borda points elimination data")
       Call Out (-1,"Decreasing mean rankings or ratings by elimination", mean)
       Call Out (-1,"Corresponding candidate orderings", elected)
     End if

   End Subroutine Borda_elimination

   Subroutine Ballot_vec (rank_pt,unrank_pt, ballot, vec) 

!    Compute the rating point vector corresponding to each ballot,
!    given ballots data as ranked candidates.

     Real,    Intent(in) :: rank_pt(:)    ! (mr) Point values for the 'mr' top ranked candidates
     Real,    Intent(in) :: unrank_pt(:)  ! (mr) Point values and for the unranked candidates, depending on # ranked
     Integer, Intent(in) :: ballot(0:,:)  ! (0:mr,nb). (0,:) = 'n' candidates ranked in preferential order

     Real,   Intent(out) :: vec(:,:)      ! (nc,nb) The point vector corresponding to each ballot
!  Local:
     Integer :: b, n

     Do b= 1,Ubound(ballot,2)
       n= ballot(0,b);  vec(:,b)= unrank_pt(n)
       vec(ballot(1:n,b),b)= rank_pt(:n)
     End do

   End Subroutine Ballot_vec

   Subroutine Ballot_vec2 (pt_val, ballot,ballot2, vec) 

!    Compute the rating point vector corresponding to each ballot,
!    given ballots data as ranked candidates.

     Real,    Intent(in) :: pt_val(:)     ! (mt) Point values for the ranking or rating levels
     Integer, Intent(in) :: ballot(0:,:)  ! (0:mr,nb). (0,:) = 'n' candidates ranked in preferential order
     Integer, Intent(in) :: ballot2(:,:)  ! (mr,nb).   Corresponding ranking or rating levels

     Real,   Intent(out) :: vec(:,:)      ! (nc,nb) The point vector corresponding to each ballot
!  Local:
     Integer :: b, n

     vec= 0
     Do b= 1,Ubound(ballot,2)
       n= ballot(0,b);  vec(ballot(1:n,b),b)= pt_val(ballot2(:n,b))
     End do

   End Subroutine Ballot_vec2

    
   Function Condorcet0 (np,nc, wtb,points)
      
!    Condorcet ordering of the candidates (thus candidate A is Codorcet top ranked if A is ranked 
!    greater than each other candidates by the underlying system = the pairwise winner), with 
!    ambiguities in the Smith set resolved by the Borda ordering of the underlyng ratings.

     Integer, Intent(in) :: np           ! # candidates to be elected
     Integer, Intent(in) :: nc           ! # candidates
      
     Real,    Intent(in) :: wtb(:)       ! (nb) Ballot weights, summing to 'np'
     Real,    Intent(in) :: points(:,:)  ! (nc,nb) The point vector corresponding to each ballot

!    Result:
     Integer :: Condorcet0(nc,4)  ! 1 : Candidates in preferential order
                                  ! 2 : Candidates grouped as successive Smith sets:
                                  !     so that Condorcet(i,2)= ... Condorcet(j,2) = k
                                  !     if these are in the kth Smith set
                                  ! 3 : the Condorcet deficit, i.e., # candidates not exceeded
                                  !     for each member of the Smith set.
                                  ! 4 : Size of the Smith set for each successive Smith set
!  Local:
     Real, Parameter :: eps= 0.0001
     Real    :: Comp(nc,nc)                  ! Comp(i,j)= ratio of # ballots with ranking(i) >= ranking(j)
                                             !            to the # ballots with ranking(j) >= ranking(i)
     Logical :: Smith(nc)                    ! Smith set at the current iteration
     Logical :: Test(nc)                     ! Set of candidates which have been tested 
                                             !   for Smith set expansion at each iteration
     Logical :: UnOrd(nc)                    ! Candidates not yet ordered at each iteration
     Integer :: ngt(nc)                      ! ngt(i) = # candidates j for which Comp(i,j) >= 1
     Integer :: ls(nc), lst(nc)
     Integer :: i, j, k, l, m, n, p
     
!    First compute the number of ballots which rank candidate 'i' higher 
!    than candidate 'j' for each pair of candidates (i,i).

     Condorcet0= 0;  Comp= 0;  Forall(i=1:nc) Comp(i,i)= 1
       
     Do i= 1,nc-1
       Do j= i+1,nc
         Comp(i,j)= Sum(wtb, points(i,:) >= points(j,:)) / np
         Comp(j,i)= Sum(wtb, points(i,:) <= points(j,:)) / np
       End do
     End do

!   Now compute the ratio of the number of ballots that rank 'i' higher than 'j' 
!   to the ratio of rankings of 'j' over 'i' for each pair (i,j) of candiates.

     Do i= 1,nc-1
       Do j= i+1,nc
         Comp(i,j)= Comp(i,j) / (Comp(j,i) + eps)
         Comp(j,i)= 1 / (Comp(i,j) + eps)
       End do
     End do
     If (pr_out > 1) Call Out (1,"Condorcet comparison matrix", Comp)
     
!    Find Condorcet winner(s), omit these, then the Condorcet winner(s)
!    of the remaining candidates, and iterate until all candidates are
!    Condorcet ordered.

     UnOrd= .true.;  n= nc;  k= 0;  m= 0
     ls= 0;  lst= "ID"

     Do
       k= k + 1;  l= m + 1
       Do i= 1,n
         ngt(i)= Count(Comp(lst(i),lst(:n)) >= 1)
       End do
       If (pr_out > 1) Call Out ("Condorcet comparison count",ngt(:n))
       p= Maxloc(ngt(:n),1)
       
       If (ngt(p) < n) then
         Smith(:n)= .false.;  Test(:n)= .false.  
         Do 
           Test(p)= .true.
           Where (Comp(lst(1:n),lst(p)) >= 1) Smith(:n)= .true.
           p= First_true(Smith(:n) .and. .not.Test(:n))
           If (p <= 0) Exit
         End do
         
         Call List_of_true (Smith(:n), j,ls)
         m= l+j-1  
         Condorcet0(l:m,1)= lst(ls(:j))
         Condorcet0(l:m,2)= k
         Condorcet0(l:m,3)= n - ngt(ls(:j))
         Condorcet0(k,4)= j

         If (pr_out > 1) then
           Call Out ("At level",k, "Smith set has size",j)
           Call Out ("Smith set",Condorcet0(l:m,1))
           Call Out ("with deficits",Condorcet0(l:m,3))
         End if
       Else
         m= l  
         Condorcet0(m,1)= lst(p)
         Condorcet0(m,2)= k
         Condorcet0(m,3)= n - ngt(p)
         Condorcet0(k,4)= 1
       End if

       UnOrd(Condorcet0(l:m,1))= .false.
       Call List_of_true (UnOrd, n,lst)
       If (n <= 0) Exit
     End do

     If (pr_out >= 1) then
       Call Out ("Condorcet ordering",Condorcet0(:,1))
       Call Out ("Condorcet ambiguity groupings",Condorcet0(:,2))
       Call Out ("Condorcet deficits",Condorcet0(:,3))
       Call Out ("Condorcet Smith set sizes",Condorcet0(:k,4))
     End if
     
   End Function Condorcet0

   Subroutine Majority_approval (np,wtb,ballot,mean,elected, ballot2)
      
!    The Majority Approval method, is a modification of both Bucklin and  
!    approval voting. It may be used with rated or ranked ballots,
!    using the ranking or rating levels but not points.

!    The goal is for an elected representative to be approved by
!    a majority of voters, or to come as close as possible, with
!    approvals being used to rank to candidates.

!    If rating or ranking levels are present, they determine the #
!    top ranked, not candidate ranking.

!    See http://wiki.electorama.com/wiki/Majority_Approval_Voting.

     Integer,  Intent(in) :: np             ! # to be elected
     Real,     Intent(in) :: wtb(:)         ! (nb) Ballot weights, summing to 'np'
     Integer,  Intent(in) :: ballot(0:,:)   ! (0:mr,nb) Ranked or rated candidates in preferential order
                                            !  n= (0,b) = # ranked.  (1:n,b) the candidates

     Real,    Intent(out) :: mean(:)        ! (nc) Decreasing mean # approvals
     Integer, Intent(out) :: elected(:)     ! (nc) Candidates in order of election

     Integer, Optional, Intent(in) :: ballot2(:,:)  ! (mr,nb) Increasing ranking or rating levels
                                                    ! that match the candidates in ballot(1:n,:)  
!  Local:
     Real, Parameter :: eps= 0.00001
     Integer :: b, i, m, n, mr, nb

     mr= Ubound(ballot,1);  nb= Ubound(ballot,2) 
     mean= 0;  elected= 0
     
     Do i= 1,mr
       mean= 0
       Do b= 1,nb
         n= ballot(0,b)
         If (Present(ballot2)) then
           m= Last_true(ballot2(:n,b) <= i)
         Else
           m= Min(i,n)
         End if

         mean(ballot(1:m,b))= mean(ballot(1:m,b)) * wtb(b)
       End do

       mean= mean / np;  Call Sort (.false., mean, elected, ez=eps)

       n= Last_true (mean > 0.50)  ! # candidates with majority approval
       If (n >= np) Exit
     End do
     
   End Subroutine Majority_approval


   Subroutine Approval (np,wtb,ballot, mean,elected)
      
!    Compute the  'mean' rating vector for the candidates 
!    and the corresponding ordering 'elected', counting a ballot
!    as a vote for a candidate if the ballot ranks the candidate
!    That is, all ranked candidates are approved.

!    Here a maximum of 'mr' approved candidates is assumed.

     Integer,  Intent(in) :: np            ! # to be elected
     Real,     Intent(in) :: wtb(:)        ! (nb) Ballot weights, summing to 'np'
     Integer,  Intent(in) :: ballot(0:,:)  ! (0:mr,nb) ballot(1:n)= candidates in the order rated. 
                                           ! n= ballot(0)= # candidates rated.

     Real,    Intent(out) :: mean(:)       ! (nc) Decreasing mean # approvals
     Integer, Intent(out) :: elected(:)    ! (nc) Corresponding candidates. 
!  Local:
     Real, Parameter :: eps= 0.00001
     Integer :: b, n, nb
     
     nb= Size(ballot,2);  mean= 0;  elected= 0
     
     Do b= 1,nb
       n= ballot(0,b);  mean(ballot(1:n,b))= mean(ballot(1:n,b)) * wtb(b)
     End do
       
     mean= mean / np
     Call Sort (.false., mean, elected, ez=eps)
     
   End Subroutine Approval

    
   Subroutine Top_vote (np,wtb,ballot, mean,elected)
      
!    In the Top Vote method the each voter's top ranked candidate gets one point,
!    and the corresponding ordering 'elected'

     Integer,  Intent(in) :: np            ! # candidates to be elected
     Real,     Intent(in) :: wtb(:)        ! (nb) Ballot weights, summing to 'np'
     Integer,  Intent(in) :: ballot(0:,:)  ! (0:mr,nb) ballot(1:n)= candidates in the order ranked. 
                                           ! n= ballot(0)= # candidates ranked.

     Real,    Intent(out) :: mean(:)       ! (nc) Mean rankings in decreasing order.
     Integer, Intent(out) :: elected(:)    ! (nc) Corresponding candidates in preferential order
     
!  Local:
     Real,   Parameter :: eps= 0.00001
     Integer :: b, c, nb
     
     nb= Ubound(ballot,2);  mean= 0;  elected= 0
     
     Do b= 1,nb
       c= ballot(1,b);  mean(c)= mean(c) + wtb(b)
     End do
       
     mean= mean / np;  Call Sort (.false., mean, elected, ez=eps)
     
   End Subroutine Top_vote

   Subroutine Rem_frac (q,nc,ind, rateV, elect)

!    Compute the set of elected candidates by the traditional method 
!    of "Largest Remainder Method" (see Wikipedia) as adapted to
!    clusters characterized by candidate rating vectors.

     Integer,  Intent(in) :: q            ! cluster set index
     Integer,  Intent(in) :: nc           ! # candidates
     Integer,  Intent(in) :: ind          ! # clusters, including independents

     Real,     Intent(in) :: rateV(0:,:)  ! (0:nc,ind) Cluster data 
                                          !   (0,:) Cluster sizes
                                          !   (1:,:) = cluster mean rating vectors

     Integer, Intent(out) :: elect(:,0:)  ! (np,0:2) (:,0) Elected candidates - increasing order
                                          !    (:,1) Elected candidates, ordered by 
                                          !          decreasing cluster size and by rating
                                          !    (:,2) The corresponding clusters
!  Local:
     Integer :: ntb(ind)      ! # candidates to be assigned to each cluster = 
                              !   output of the traditional largest remainder method
     Integer :: ordZ(ind)     ! Stable sort of the clusters into increasing size
     Integer :: ordV(nc,ind)  !  Order of the candidates for each cluster 
                              !    by decreasing mean rating 

     Real, Parameter :: eps= 0.00001
     Logical :: Elected(nc)
     Real    :: tot_wt, siz(ind), rmd(ind), ratio(ind), tmp(nc)
     Integer :: ordR(ind)
     Integer :: i, j, k, l, n, p, cl, cn, mx,np, ncl
     
!    Elect candidates to fill whole number parts of the clusters

     ncl= ind-1;  np= Size(elect,1)
     elect= -1;  Elected= .false.
     
!    Stable sort the regular clusters by increasing size. For each cluster, 
!    sort its candidates by decreasing mean rating

     siz= rateV(0,:);  tot_wt= Sum(siz);  siz= (np / tot_wt) * siz

     ordZ= "ID";  If (ncl > 1) Call Sort (.true., siz(:ncl), ordZ(:ncl), ez=eps)
     
     Do k= 1,ind
       cl= ordZ(k);  tmp= rateV(1:,cl)
       Call Sort (.false., tmp, ordV(:,k), ez=eps)
     End do
     
!    Determine the number ntb(cl) to be elected per cluster 'cl', 
!    including independents.
     
     Do k= 1,ind
       ntb(k)= Floor(siz(k));  rmd(k)= siz(k) - ntb(k)
     End do
     
     If (pr_out > 1.0) then
       Call Out ("Rem_frac: For cluster set",q,"with # reg clusters",ncl, ln=1)
       Call Out ("Reg cluster reordering by increasing cluster size", ordZ)
       Call Out (-1,"Mean vector reordering (decreasing) for reordered clusters", ordV)
       Call Out ("Initial # elected from each reordered cluster",ntb)
     End if

     Call Sort (.false.,rmd, ordR, ez=eps);  p= Sum(ntb)
     
     Do k= 1,ind
       If (p < np) then
         ntb(ordR(k))= ntb(ordR(k)) + 1;  p= p + 1
       Else
         Exit
       End if
     End do
     
     If (pr_out > 1.0) then
       Call Out ("Decreasing remainders",rmd)
       Call Out ("corresponding cluster reordering",ordR)
       Call Out ("Final # elected from each reordered cluster",ntb)
     End if
     
     mx= Max(np,ind)
     Call Assign_to_clust (q,mx,ordV,ordZ, ntb, Elected,elect)
     
   End Subroutine Rem_frac
   
   Subroutine DHondt (q,nc,ind, rateV, elect)

!    Compute the set of elected candidates by the D'Hondt's Rule
!    (see Wikipedia) as adapted to candidate ratings.
!    See also www.ucl.ac.uk/~ucahhwi/dhondt.pdf.
   
!    This method may be described as the "largest remaining ratio method"
!    in contrast to the "largest reamining fraction" method. That is the
!    next candiate elected is the one with the largest ratio of actual
!    cluster size to its # elected representatives to date + 1. In fact DHondt
!    is a minimax method that minimizes the maximum value of this ratio 
!    over all voting blocks to prevent the worst disproportionality.
   
!    This may also be viewed as an attempt to maximize the average fractional 
!    # elected per weighted cluster as measured by a harmonic mean:  
!    (np / (Sum(siz(i)/(ntb(i)+1))) - 1 over clusters 'i' with sizes 'siz(i)' 
!    in position units and # to be elected per cluster 'ntb(i)'. 
!    The method is a greedy algorithm. 
   
!    DHondt sometimes favors larger clusters in comparison to the largest remainder 
!    method. The ratios focus on minimizing the disproportionality as a margin of 
!    error in the percentage of the vote (% of representatives assigned to a 
!    voting block versus its true %), instead of on the matching the exact size 
!    of the voting block, measured in units of elected positions.

     Integer,  Intent(in) :: q            ! cluster set index
     Integer,  Intent(in) :: nc           ! # candidates
     Integer,  Intent(in) :: ind          ! # clusters, including independents
     Real,     Intent(in) :: rateV(0:,:)  ! (0:nc,ind) Cluster data 
                                          !   (0,:) Cluster sizes
                                          !   (1:,:) = cluster mean rating vectors

     Integer, Intent(out) :: elect(:,0:)  ! (np,0:2) (:,0) Elected candidates - increasing order
                                          !    (:,1) Elected candidates, ordered by 
                                          !          decreasing cluster size and by rating
                                          !    (:,2) The corresponding clusters
!  Local:
     Integer :: ntb(ind)      ! # candidates to be assigned to each cluster = 
                              !   output of the traditional D'Hondt method
     Integer :: ordZ(ind)     ! Stable sort of the clusters into increasing size
     Integer :: ordV(nc,ind)  !  Order of the candidates for each cluster 
                              !    by decreasing mean rating 

     Real, Parameter :: eps= 0.00001
     Logical :: Elected(nc)
     Real    :: tot_wt, siz(ind), ratio(ind), tmp(nc)
     Integer :: i, j, k, l, n, p, cl, cn, mx, np, p0, ncl
     
     ncl= ind-1;  np= Size(elect,1)
     elect= -1;  Elected= .false.
     
!    Stable sort the regular clusters by increasing size. For each cluster, 
!    sort its candidates by decreasing mean rating

     siz= rateV(0,:);  tot_wt= Sum(siz);  siz= (np / tot_wt) * siz

     ordZ= "ID";  If (ncl > 1) Call Sort (.true., siz(:ncl), ordZ(:ncl), ez=eps)
     
     Do k= 1,ind
       cl= ordZ(k);  tmp= rateV(1:,cl)
       Call Sort (.false., tmp, ordV(1:,k), ez=eps)
     End do
     
     If (pr_out > 1.0) then
       Call Out ("DHondt: For cluster set",q,"with # reg clusters",ncl, ln=1)
       Call Out ("Reg cluster reordering by increasing cluster size", ordZ)
       Call Out (-1,"Mean vector reordering (decreasing) for reordered clusters", ordV)
     End if  
     
!    Determine the number of candidates to represent each cluster
     
     ratio= siz;  ntb= 0
     Do i= 1,np
       k= Maxloc(ratio,1);  ntb(k)= ntb(k) + 1
       ratio(k)= siz(k) / (ntb(k) + 1)
     End do
     
     If (pr_out > 1.0) then
       Call Out ("# elected from each reordered cluster",ntb)
       Call Out ("Average vote per # elected + 1, by reordered cluster",ratio)
     End if
     
     mx= Max(np,ind)
     Call Assign_to_clust (q,mx,ordV,ordZ, ntb, Elected,elect)
     
   End Subroutine DHondt
   
   Subroutine Assign_to_clust (q,mx,ordV,ordZ, ntb, Elected,elect)
   
!    Assign 'np' candidates to a prior allocation of 
!    positions 'ntb(k)' for each cluster 'k', 
!    assuming Sum(ntb)= np, using the candidate ordering 
!    ordV(:,cl) for each cluster 'cl'.
   
!    This is a "bottom-up, round-robin" algorithm.

!    Start with the smallest cluster getting its first choice, 
!    the next largest cluster its first choice, etc., in a sweep 
!    through the clusters, except that if a candidate was already 
!    chosen by a previous cluster then the current cluster switches 
!    to the next unchosen candidate in the mean rating vector. 
   
!    These sweeps through the clusters continue until all candidates 
!    have been assigned according to the cluster allotments ntb(k). 
!    Except that we also make provision for a group of voters 
!    that represent "none of the above". This group of independents 
!    comes last in the sweep but will still get its allotment of candidates.
   
     Integer,  Intent(in) :: q            ! cluster set index
     Integer,  Intent(in) :: mx           ! local array size = max(np,ind)
     Integer,  Intent(in) :: ordV(:,:)    ! (nc,ind)  Ordering of the candidates by decreasing rating
                                          !           for clusters ordered by increasing size
     Integer,  Intent(in) :: ordZ(:)      ! (ind)  Ordering of the reg clusters by increasing size
     Integer,  Intent(in) :: ntb(:)       ! (ind)  # candidates to be elected by increasing size cluster
     Logical, Intent(out) :: Elected(:)   ! (nc)   True for elected candidates
     Integer, Intent(out) :: elect(:,0:)  ! (np,0:2) (:,0) Elected candidates - increasing order
                                          !    (:,1) Elected candidates, ordered by 
                                          !          decreasing cluster size and by rating
                                          !    (:,2) The corresponding clusters
     
!  Local:
     Integer :: ls(mx), cls(mx), key(mx), nel(mx), nxt(mx), invZ(mx)
     Integer :: i, j, k, n, p, cl, cn, nc, np, ncl, ind
     
     nc= Size(ordV,1);  ind= Size(ordV,2);  np= Size(elect,1)
     ncl= ind - 1

     invZ= 0; Forall(i=1:ind) invZ(ordZ(i))= i
     
     Elected= .false.;  elect= 0
     p= 0;  nel= 0;  nxt= 1
     
     Top_loop : Do
       i= 0  
       Clust_loop : Do k= 1,ncl
         n= ntb(k) - nel(k)  ! 'n' = # remaining to be assigned to cluster 'k'
         If (n <= 0) Cycle Clust_loop
       
         i= i + 1;  cn= ordV(nxt(k),k)
         
         If (Elected(cn)) then
           j= First_true(.not.Elected(ordV(nxt(k)+1:,k)))
           nxt(k)= nxt(k) + j;  cn= ordV(nxt(k),k)
         End if
         
!        Elect 'cn'
            
         p= p + 1;  elect(p,1)= cn;  elect(p,2)= invZ(k)
         Elected(cn)= .true.;  nel(k)= nel(k) + 1;  nxt(k)= nxt(k) + 1
         If (p >= np) Exit Top_loop 
       End do Clust_loop

       If (i < 1) Exit Top_loop 
     End do Top_loop

!    Elect any remaining candidates from the set of independents

     If (p < np) then
       k= ind;  n= ntb(k) - nel(k)

       Ind_loop : Do i= 1,n
         cn= ordV(nxt(k),k)
         If (Elected(cn)) then
           j= First_true(.not.Elected(ordV(nxt(k)+1:,k)))
           nxt(k)= nxt(k) + j;  cn= ordV(nxt(k),k)
         End if
         
         p= p + 1;  elect(p,1)= cn;  elect(p,2)= ind
         Elected(cn)= .true.;  nel(k)= nel(k) + 1;  nxt(k)= nxt(k) + 1
         If (p >= np) Exit Ind_loop 
       End do Ind_loop
     End if

     elect(:,0)= elect(:,1);  Call Sort (.true.,elect(:,0))

!    Reverse the order of the regular clusters 
!    (from increasing cluster size to decreasing size)

     Call List_of_true (elect(:,2) < ind, n,ls)
     cls(:n)= elect(ls(:n),2);  Call Sort (.true., cls(:n), key(:n))

     elect(ls(:n),1)= elect(ls(key(:n)),1);  elect(ls(:n),2)= cls(:n)

     If (pr_out > 1.0) then
       Call Out ("For cluster set",q,ln=1)
       Call Out (-1,"Elected candidates for reg clusters of decreasing size",elect)
     End if
   End Subroutine Assign_to_clust


  Subroutine IRV (nc, wtb,ballot, elected,vote)

!   Standard "Instant Runoff Voting"
  
!   For each ballot, the total weight is always on the first "hopeful" candidate.
!   The elected candidates are the ones left after successive elimination of the
!   candidate cluster with the lowest size (= weight = vote)

! Input:
    Integer,    Intent(in) :: nc            ! # candidates
    Real,       Intent(in) :: wtb(:)        ! (nb) Ballot weights, summing to 'np'
    Integer,    Intent(in) :: ballot(0:,:)  ! (0:mr,nb) Candidates chosen by each ballot
                                            ! in the order ranked. ballot(0,nb)= # candidates
                                            ! remaining after elimination to date.
! Output:
    Integer,   Intent(out) :: elected(:)    ! (nc) List of elected candidates, in preferential order
    Real,      Intent(out) :: vote(:)       ! (nc) Candidate cluster vote sizes in position units
! Local:
    Real    :: vt_ord(nc)       ! Cluster vote sizes reordered to decreasing
    Integer :: Elect(nc)        ! = 0 if hopeful candidate, = -i if ith candidate eliminated
    Integer :: key(nc)          ! Re-ordering key
    Integer :: hope1(Size(wtb)) ! Index of the first hopeful candidate in the ballot, or 0 if none
    Real    :: swt
    Integer :: b, j, k, l, n, cn, nr, last, nb, nd, np

!   Initialize
    
    nb= Size(wtb);  np= Nint(Sum(wtb))
    elected= "ID";  vote= 0;  hope1= 1;  Elect= 0;  nd= 0

    Call Out ("Enter 'IRV' with # candidates",nc, "# to be elected",np,ln=1)  

    Do cn= 1,nc
      vote(cn)= Sum(wtb, mask= ballot(1,:) == cn)
    End do
    
!   Elimination loop
    
    Elim_cycle : Do nr= nc,np+1,-1
      vt_ord(:nr)= vote(elected(:nr))
      
      If (nr > 1) then
        Call Sort (.false., vt_ord(:nr), key(:nr))
        elected(:nr)= elected(key(:nr))
      End if

      If (pr_out > 1) then
        Call Out ("At elimination place",nr,ln=1)  
        Call Out ("Candidates sorted by decreasing cluster size", elected)
        Call Out ("from cluster sizes", vote)
      End if

!     Eliminate the last place candidate from all the higher place lists

      last= elected(nr);  l= 0;  swt= 0
      
      Elim_loop : Do b= 1,nb  
        k= hope1(b);      If (k <= 0)   Cycle Elim_loop
        cn= ballot(k,b);  If (cn /= last) Cycle Elim_loop
          
        n= ballot(0,b);  j= k + First_true(Elect(ballot(k+1:n,b)) == 0)

        If (j <= k) then  ! No more hopefuls. Cannot distribute the ballot weight
          hope1(b)= 0;  Cycle Elim_loop
        Else  ! Subsequent hopeful
          l= l + 1;  hope1(b)= j
        End if

        swt= swt + wtb(b);  cn= ballot(j,b)
        vote(last)= vote(last) - wtb(b)
        vote(cn)= vote(cn) + wtb(b)
      End do Elim_loop
      
      If (pr_out > 1) then
        Call Out ("Last place candidate - to be eliminated", last)
        Call Out ("From # ballot transfers",l, "with total weight",swt)
        Call Out ("Updated cluster sizes", vote)
      End if    
      
      nd= nd + 1;  Elect(last)= -nd
    End do Elim_cycle

    swt= Sum(vote)

    If (pr_out >= 1) then
      Call Out ("Final sum of clusters weight", swt)
      Call Out ("Final cluster sizes", vote)
      Call Out ("The electeds in preferential order", elected)
      Call Out ("Final status of hopefuls and deleted", Elect)
    End if    
  End Subroutine IRV

  Subroutine IRV2 (nc, pt_val, wtb,ballot,ballot2, elected,vote)

!   "Instant Runoff Voting" using rates to help determine the order of elimination:
!   candidate mean rating * vote size.
  
!   For each ballot, the total weight is always on the first "hopeful" candidate.
!   The elected candidates are the ones left after successive elimination of the
!   candidate cluster with the lowest size (= weight = vote)

! Input:
    Integer,  Intent(in) :: nc              ! # candidates
    Real,     Intent(in) :: pt_val(:)       ! (mt) Rating points
    Real,     Intent(in) :: wtb(:)          ! (nb) Ballot weights, summing to 'np'
    Integer,  Intent(in) :: ballot(0:,:)    ! (0:mr,nb). (0,:) = 'n' candidates in preferential order
    Integer,  Intent(in) :: ballot2(0:,:)   ! (0:mr,nb). corresponding ranking or rating levels
                                            !   (0,:) = 'p' = # positively rated canddidates        
! Output:
    Integer,   Intent(out) :: elected(:)    ! (nc) List of elected candidates, in preferential order
    Real,      Intent(out) :: vote(:)       ! (nc) Candidate cluster vote sizes in position units
! Local:
    Real    :: rate(nc)    ! (nc) Mean rating over the hopeful candidates
    Real    :: rateV(nc)   ! (nc) rate * vote over the hopeful candidates
    Real    :: rt_ord(nc)  ! (nc) rateV as reordered

    Integer :: Elect(nc)   ! (nc) = 0 if hopeful candidate, = -i if ith candidate eliminated
    Integer :: hope1(Size(wtb))  ! Index of the next positive hopeful candidate in the ballot, or 0 if none

    Integer :: lh(nc), key(nc), lst(nc)
    Real    :: swt, retain, trans_wt
    Integer :: b, j, k, l, m, n, p, r, cn, nb, nd, nh, np, nr, last

!   Initialize
    
    nb= Size(ballot,2);  np= Nint(Sum(wtb))
    elected= "ID";  hope1= 1;  vote= 0;  Elect= 0;  key= 0;  nd= 0
    
    Call Out ("Enter 'IRV2' with # candidates",nc, "# to be elected",np,ln=1)  

    Do cn= 1,nc
      vote(cn)= Sum(wtb, mask= ballot(1,:) == cn)
    End do
    
!   Elimination loop
    
    Elim_cycle : Do nr= nc,np,-1

!     Compute 'rateV' = 'rate' * 'vote' over the posivtively rated hopeful candidates

      rate= 0
      Rate_loop : Do b= 1,nb
        n= ballot2(0,b)
        Call List_of_true (Elect(ballot(1:n,b)) == 0, m,lh)
        If (m > 0) then
          lst(:m)= ballot(lh(:m),b)
          rate(lst(:m))= rate(lst(:m)) + pt_val(ballot2(lh(:m),b)) * wtb(b)
        End if
      End do Rate_loop 

      rate= rate / np;  rateV= rate * vote

      Call List_of_true (Elect == 0, nh,key)

      If (nh > 1) then  ! Sort the hopeful candidates by rateV
        rt_ord(:nh)= rateV(key(:nh))
        Call Sort (.false., rt_ord(:nh), lh(:nh))
        key(:nh)= key(lh(:nh))
      Else
        cn= key(1);  rt_ord(1)= rateV(cn)
      End if

      last= key(nh);  elected(:nr)= key(:nr)

      If (pr_out >= 1) then
        Call Out ("At elimination place",nr, "with last candidate", last, ln=1)  
        Call Out ("Candidates sorted by decreasing mean rate * cluster size", elected)
        Call Out ("from cluster sizes", vote)
      End if

      If (nr <= np) Exit Elim_cycle

!     Eliminate the last candidate, keeping all weight on the 1st hopeful of each ballot

      m= 0;  l= 0;  retain= 0;  swt= 0

      Elim_loop : Do b= 1,nb
        k= hope1(b);      If (k <= 0)     Cycle Elim_loop  ! No more positive hopefuls
        cn= ballot(k,b);  If (cn /= last) Cycle Elim_loop  ! 1st hopeful /= 'last'

        p= ballot2(0,b);  j= k + First_true(Elect(ballot(k+1:p,b)) == 0)
        trans_wt= wtb(b) 

        If (j <= k) then  ! No more positive hopefuls. Cannot distribute the ballot weight
          m= m + 1;  retain= retain + trans_wt
          hope1(b)= 0
        Else  ! Subsequent positive hopeful
          vote(last)= vote(last) - trans_wt
          cn= ballot(j,b);  vote(cn)= vote(cn) + trans_wt
          l= l + 1;  hope1(b)= j
        End if
      End do Elim_loop
      
      nd= nd + 1;  Elect(last)= -nd

      If (pr_out >= 1) then
        Call Out ("Deletion #",nd,"with candidate deleted", last, ln=1)
        Call Out ("# ballot transfers",l, "with total weight",swt)
        Call Out ("# ballots retained",m, "with total weight",retain)
        Call Out ("Updated cluster sizes", vote)
      End if    
    End do Elim_cycle

    swt= Sum(vote)

    If (pr_out >= 1) then
      Call Out ("Final sum of clusters weight", swt)
      Call Out ("Final cluster sizes", vote)
      Call Out ("The electeds in preferential order", elected)
      Call Out ("Final status of hopefuls and deleted", Elect)
    End if    
  End Subroutine IRV2

  Subroutine STV (nc, wtb,ballot, elected, STVmb)

!   This algorithm is an extension of STV, once the final 
!   candidate is elected, to produce better clustering. The
!   concept is to redistribute subsequent transfers
!   from the first 'hopeful' candidate remaining in a ballot 
!   to the first 'prior elected' or 'hopeful' candidate.

!   That is, it is no longer necessary to transfer the 
!   maximum possible in excess of the Droop quota 
!   to subsequent hopefuls in order to qualify one 
!   of them for election.
!
!   The effect is to redistribute some of the ballot weight 
!   to the candidate clusters of the prior electeds, so that
!   thier size may exceed the quota. Thus all subsequent 
!   hopeful candidate clusters, as they are deleted, will 
!   also redistribute some of their weight to the electeds
!   diminishing size of the pseuo-cluster of independents.

! Input:
    Integer,         Intent(in) :: nc           ! # candidates
    Real,            Intent(in) :: wtb(:)       ! (nb) Ballot weights, summing to 'np'
    Integer,         Intent(in) :: ballot(0:,:) ! (1:mr,nb) Candidates ranked by each ballot
                                                ! in the order ranked, with ballot(0,b)= # ranked 
                                                ! by ballot 'b' up to mr
! Output:
    Integer,        Intent(out) :: elected(:)   ! (nc) List of elected candidates in preferential order

    Type(Set_list), Intent(out) :: STVmb(0:)    ! (0:nb) Ballot ranking structure for clustering
                                                ! Ballot sets (0): 
                                                !   %n        = nc
                                                !   %p        = np
                                                !   %svl      = elected vote = Sum(%val(set(:np)))
                                                !   %smb      = independents vote = Sum(%val(set(np+1:)))
                                                !   %set(nc)  = candidates in elected ordering, with
                                                !               the first 'np' elected
                                                !   %lev(nc)  = candidate clusters by decreasing cluster vote
                                                !   %val(nc)  = candidate cluster votes, not reordered

                                                ! Ballots (1:nb):
                                                !   %n     = # elected candidates
                                                !   %p     = # positively rated candidates
                                                !   %r     = # ranked candidates
                                                !   %set(r)= 'n' elected candidates, then unelected,
                                                !            in original ranking order
                                                !   %val(r)= corresponding membership weights 
                                                !   %svl   = elected membership weight= Sum(%val(:n))
                                                !   %smb   = total membership weight  = Sum(%val)
! Local:
    Real    :: vote0(nc)     ! Size of initial candidate clusters
    Real    :: vote(nc)      ! Size of candidate clusters at current STV cycle
    Real    :: vt_ord(nc)    ! Hopeful candidate clusters, ordered by decreasing size
    Real    :: quota(nc)     ! Droop formula quotas by STV cycle

    Real    :: w_trans(nc,3) ! Records data by STV cycle
                             ! (:,1) = weight transfered to electeds, (:,2) = total weight transfered
                             ! (:,3) = weight of fully retained ballots
    Integer :: n_trans(nc,3) ! Records data by STV cycle
                             ! (:,1) = # ballot transfers, (:,2) = # ballots retained
                             ! (:,3) = deleted (neg) or elected (pos) candidate

    Integer :: Elect(nc)     ! (i) = i if ith elected, 0 if hopeful, -i if ith eliminated
    Integer :: Elect1(nc)    ! = 'Elect' except < 0 where 'Elect' > 0 but not all elected

    Real    :: vote1(nc)
    Integer :: lh(nc), key(nc)

    Real, Parameter :: min_trans= 0.001
    Logical :: Add

    Real    :: retain, surplus, trans_wt, transfer_fac, transferable, independents
    Real    :: vt, inc, nc1, swt, sm_pt, tot_hope, prior, total, vt_elc, vt_ind
    Integer :: c1, cn, h1, il, it, iv, mr, n1, nb, nh, ne, np, nr, top, nd, neh, np1, last
    Integer :: b, j, k, l, m, n, r

    nb= Size(wtb);  np= Nint(Sum(wtb));  mr= Ubound(ballot,1)
    np1= np + 1;  nc1= nc + 1
    
    Call Out ("STV algorithm for # electeds",np, "out of # candidates",nc,ln=1)

!   Initial ballot memberships in candidate clusters

    vote0= 0;  Elect= 0;  elected= 0;  independents= 0
    quota= 0;  n_trans= 0;  w_trans= 0

    Ballot_loop0 : Do b= 1,nb
      STVmb(b)%smb= wtb(b);  r= ballot(0,b)

      If (r < 1 .or. r > mr .or. &
        Any(ballot(1:r,b) < 1) .or. Any(ballot(1:r,b) > nc)) then
        Call Out ("Warning in STV: Bad ballot",b); Cycle Ballot_loop0
      End if

      If (Associated(STVmb(b)%set)) DeAllocate(STVmb(b)%set)
      If (Associated(STVmb(b)%val)) DeAllocate(STVmb(b)%val)
      Allocate(STVmb(b)%set(r), STVmb(b)%val(r))

      STVmb(b)%set= ballot(1:r,b);  
      STVmb(b)%val(1)= wtb(b);  STVmb(b)%val(2:)= 0

      cn= ballot(1,b);  vote0(cn)= vote0(cn) + wtb(b)
    End do Ballot_loop0

    If (pr_out >= 1) Call Out ("Initial STV candidate cluster sizes, unordered", vote0)
      
    vote= vote0
    key= 'ID'  ! Hopeful candidates, to be sorted by weight
    ne= 0      ! # currently elected
    nd= 0      ! # currenlty deleted
      
    STV_cycle : Do iv= 1,nc
        
      Call List_of_true (Elect == 0, nh,key)  ! nh = # hopeful candidates remaining
      If (nh <= 0) Exit STV_cycle             ! Done

      If (nh > 1) then   ! Sort the hopeful candidates from largest to smallest
        vt_ord(:nh)= vote(key(:nh))
        Call Sort (.false., vt_ord(:nh), lh(:nh))

        key(:nh)= key(lh(:nh));  tot_hope= Sum(vt_ord(:nh))
      Else
        vt_ord(1)= vote(key(1));  tot_hope= vt_ord(1)
      End if

      top= key(1);  last= key(nh)
      quota(iv)= tot_hope / (np1 - ne)  ! Droop formula for the quota

      prior= vt_ord(1);  surplus= prior - quota(iv)
      Add= ne < np .and. surplus >= 0    ! 'surplus' = weight to be transferred if > 0

      If (pr_out >= 1) then
        Call Out ("Hopeful candidates sorted by decreasing vote sizes", key(:nh))
        Call Out ("for vote sizes", vt_ord(:nh))
        Call Out ("With total hopeful weight",tot_hope, "quota",quota(iv), ln=1)
        Call Out ("& surplus",surplus)
      End if

      If (.not.Add) then ! Eliminate the bottom  hopeful 'last'
        Call Out ("Delete bottom hopeful candidate", last, ln=1)

        m= 0;  l= 0;  retain= 0;  swt= 0;  vt_elc= 0
        Elect1= Elect;  If (ne < np) Where(Elect > 0) Elect1= -10  ! Distribute only to next hopeful
        
        Elim_loop : Do b= 1,nb
          il= First_true(STVmb(b)%set == last);  If (il < 1) Cycle Elim_loop  ! 'last' not ranked
          trans_wt= STVmb(b)%val(il);  If (trans_wt <= 0) Cycle Elim_loop     ! No weight to transfer

          j= First_true(Elect1(STVmb(b)%set) >= 0 .and. STVmb(b)%set /= last) ! Next eligible candidate

          If (j > 0) then  ! Transfer weight & distribute to the 'j' candidate cluster
            STVmb(b)%val(il)= 0;  vote(last)= vote(last) - trans_wt

            STVmb(b)%val(j)= STVmb(b)%val(j) + trans_wt
            cn= STVmb(b)%set(j);  vote(cn)= vote(cn) + trans_wt
            If (Elect(STVmb(b)%set(j)) > 0) vt_elc= vt_elc + trans_wt

            l= l + 1;  swt= swt + trans_wt
          Else      ! No distribution: fully retain
            m= m + 1;  retain= retain + trans_wt
          End if
        End do Elim_loop

        nd= nd + 1;  Elect(last)= -nd;  elected(nc1 - nd)= last
        vt_ind= vote(last);  total= sum(vote)

        independents= independents + vt_ind

        n_trans(iv,1)= l;       n_trans(iv,2)= m;    n_trans(iv,3)= -last
        w_trans(iv,1)= vt_elc;  w_trans(iv,2)= swt;  w_trans(iv,3)= retain 
      
        If (pr_out > 1) then
          Call Out ("For elimination #",nd, "# ballot transfers",l, ln=1)
          Call Out ("Total weight transferred",swt, "weight to electeds",vt_elc)
          Call Out ("Remaining weight to independents",vt_ind, &
                    "Current independents weight",independents)
          Call Out ("Current vote total",total)
          Call Out ("Current cluster sizes",vote)
        End if

      Else  ! Elect the top candidate 'top'
        Call Out ("Add top hopeful candidate",top,ln=1)
        
!       Identify ballots with a transferable ballot weight for candidate 'top'

        ne= ne + 1;  m= 0;  l= 0;  retain= 0;  transferable= 0
        Elect1= Elect;  If (ne < np) Where(Elect > 0) Elect1= -10  ! Distribute only to next hopeful
          
        Add_loop0 : Do b= 1,nb
          STVmb(b)%n= 0;     it= First_true(STVmb(b)%set == top)
            STVmb(b)%r= it;  If (it < 1) Cycle Add_loop0  ! 'top' not ranked
          trans_wt= STVmb(b)%val(it);  If (trans_wt <= 0) Cycle Add_loop0     ! No weight to transfer

          j= First_true(Elect1(STVmb(b)%set) >= 0 .and. STVmb(b)%set /= top)
          STVmb(b)%n= j

          If (j > 0) then
            l= l + 1;  transferable= transferable + STVmb(b)%val(it)
          Else ! No distribution: fully retain
            m= m + 1;  retain= retain + STVmb(b)%val(it)
          End if
        End do Add_loop0

        If (pr_out > 1) then
          Call Out ("Transferable weight",transferable, "from # possible ballot transfers",l)
        End if

 !      Transfer excess ballot weight
         
        If (transferable > min_trans) then
          transfer_fac= Min(surplus / transferable, 1.0)
          l= 0;  swt= 0;  vt_elc= 0
          
          Add_loop : Do b= 1,nb
            If (STVmb(b)%n <= 0) Cycle Add_loop

!           Subtract from the 'top' candidate cluster

            it= STVmb(b)%r;  trans_wt= transfer_fac * STVmb(b)%val(it)
            STVmb(b)%val(it)= STVmb(b)%val(it) - trans_wt
            vote(top)= vote(top) - trans_wt

!           Distribute to the candidate cluster specified by 'j'

            j= STVmb(b)%n;  STVmb(b)%val(j)= STVmb(b)%val(j) + trans_wt
            cn= STVmb(b)%set(j);  vote(cn)= vote(cn) + trans_wt
            If (Elect(cn) > 0) vt_elc= vt_elc + trans_wt

            l= l + 1;  swt= swt + trans_wt
          End do Add_loop

          Elect(top)= ne;  elected(ne)= top
          n_trans(iv,1)= l;       n_trans(iv,2)= m;    n_trans(iv,3)= top
          w_trans(iv,1)= vt_elc;  w_trans(iv,2)= swt;  w_trans(iv,3)= retain

          vt= prior - vote(top);  total= Sum(vote)  

          If (pr_out >= 1) then
            Call Out ("# ballot transfers, # retained ballots, elected candidate",n_trans(iv,:))
            Call Out ("Weight transfered to electeds, total transfered, total fully retained",w_trans(iv,:))
            Call Out ("top candidate vote decrease",vt, "total vote",total)
            Call Out ("for current cluster sizes",vote)
          End if
        Else  ! Elect without transfering
          Elect(top)= ne;  elected(ne)= top
          n_trans(iv,2)= m;  n_trans(iv,3)= top;  w_trans(iv,3)= retain

          Call Out ("No transfers: STV elected candidate",top, &
                    "retains vote",vote(top), ln=1)
        End if
      End if
    End do STV_cycle

    If (pr_out >= 0.5) then
      Call Out (1,"STV # ballot transfers, # retained ballots, elected candidate",n_trans)
      Call Out (1,"Weight transfered to electeds, total transfered, total fully retained",w_trans)
      Call Out ("Final cluster sizes",vote)

      vt_elc= Sum(vote(elected(:np)));  vt= Sum(vote) 
      Call Out ("Final elected total",vt_elc, "vote total",vt)
    End if

!   Final STVmb output: Redo ballots so that the electeds 
!   come first in each ballo. Otherwise maintaining the original
!   ranking order.

    Call Final_STVmb (mr,nc, elected,Elect, STVmb, vote1)

  End Subroutine STV
   

  Subroutine STV2 (nc,pt_val, wtb,ballot,ballot2, elected, STVmb)

!   This algorithm is an extension of STV, once the final 
!   candidate is elected, to produce better clustering. The
!   concept is to redistribute subsequent transfers
!   from the first 'hopeful' candidate remaining in a ballot 
!   to the first 'prior elected' or 'hopeful' candidate.

!   That is, it is no longer necessary to transfer the 
!   maximum possible in excess of the Droop quota 
!   to subsequent hopefuls in order to qualify one 
!   of them for election.
!
!   The effect is to redistribute some of the ballot weight 
!   to the candidate clusters of the prior electeds, so that
!   thier size may exceed the quota. Thus all subsequent 
!   hopeful candidate clusters, as they are deleted, will 
!   also redistribute some of their weight to the electeds
!   diminishing size of the pseuo-cluster of independents.

! Input:
    Integer,  Intent(in) :: nc              ! # to be elected, # candidates, # ballots
    Real,     Intent(in) :: pt_val(:)       ! (mt) Rating or ranking points corresponding to 
                                            !      rating or ranking levels.
    Real,     Intent(in) :: wtb(:)          ! (nb) Ballot weights, summing to 'np'
    Integer,  Intent(in) :: ballot(0:,:)    ! (0:mr,nb). (0,:) = 'n' candidates in preferential order
    Integer,  Intent(in) :: ballot2(0:,:)   ! (mr,nb). corresponding ranking or rating levels
! Output:
    Integer,        Intent(out) :: elected(:)   ! (nc) List of elected candidates in preferential order

    Type(Set_list), Intent(out) :: STVmb(0:)    ! (0:nb) Ballot ranking structure for clustering
                                                ! Ballot sets (0): 
                                                !   %n        = nc
                                                !   %p        = np
                                                !   %svl      = elected vote = Sum(%val(set(:np)))
                                                !   %smb      = independents vote = Sum(%val(set(np+1:)))
                                                !   %set(nc)  = candidates in elected ordering, with
                                                !               the first 'np' elected
                                                !   %lev(nc)  = candidate clusters by decreasing cluster vote
                                                !   %val(nc)  = candidate cluster votes, not reordered

                                                ! Ballots (1:nb):
                                                !   %n     = # elected candidates
                                                !   %p     = # positively rated candidates
                                                !   %r     = # ranked candidates
                                                !   %set(r)= 'n' elected candidates, then unelected,
                                                !            in original ranking order
                                                !   %val(r)= corresponding membership weights 
                                                !   %svl   = elected membership weight= Sum(%val(:n))
                                                !   %smb   = total membership weight  = Sum(%val)
! Local:
    Real    :: vote0(nc)     ! Original size of each candidate cluster
    Real    :: vote(nc)      ! Currrent size of each candidate cluster
    Real    :: vt_ord(nc)    ! (nh) Hopeful candidate clusters, by decreasing size
    Real    :: quota(nc)     ! Droop formula quotas by STV cycle

    Real    :: rate(nc)      ! Mean rating vector of the hopefuls
    Real    :: rateV(nc)     ! Current 'rate' * 'vote'
    Real    :: rt_ord(nc)    ! (nh) Hopeful candidate clusters, by decreasing 'rateV'

    Real    :: w_trans(nc,3) ! Records data by STV cycle
                             ! (:,1) = weight transfered to electeds, (:,2) = total weight transfered
                             ! (:,3) = weight of fully retained ballots
    Integer :: n_trans(nc,3) ! Records data by STV cycle
                             ! (:,1) = # ballot transfers, (:,2) = # ballots retained
                             ! (:,3) = deleted (neg) or elected (pos) candidate

    Integer :: Elect(nc)   ! = i if ith elected, 0 if hopeful, -i if ith eliminated
    Integer :: Elect1(nc)  ! 'Elect' except < 0 where 'Elect' > 0 but not all elected

    Real    :: pt(nc), vote1(nc)
    Integer :: lh(nc), lst(nc), key(nc)

    Real, Parameter :: min_trans= 0.001, eps= 0.00001
    Logical :: Add
    Real    :: retain, surplus, trans_wt, transfer_fac, transferable, independents
    Real    :: vt, inc, swt, sm_pt, tot_hope, prior, total, vt_elc, vt_ind
    Integer :: mr, mt, nb, nd, n1, nh, ne, np, nr, nc1, ndf, np1, top, last
    Integer :: c1, cn, h1, il, it, iv 
    Integer :: b, j, l, m, n, p, r

    mr= Ubound(ballot,1);  mt= Size(pt_val);  np= Nint(Sum(wtb))
    np1= np + 1;  nc1= nc + 1;  nb= Size(wtb);  ndf= nc - np

    Call Out ("STV2 algorithm for # electeds",np, "out of # candidates",nc,ln=1)

!   Initial ballot memberships in candidate clusters

    elected= 0;  vote0= 0;  Elect= 0;  independents= 0
    quota= 0;  n_trans= 0;  w_trans= 0

    Ballot_loop0 : Do b= 1,nb
      STVmb(b)%smb= wtb(b);  r= ballot(0,b);  p= ballot2(0,b)
      STVmb(b)%p= p

      If (Associated(STVmb(b)%set)) DeAllocate(STVmb(b)%set)
      If (Associated(STVmb(b)%lev)) DeAllocate(STVmb(b)%lev)
      If (Associated(STVmb(b)%val)) DeAllocate(STVmb(b)%val)
      Allocate(STVmb(b)%set(r), STVmb(b)%lev(r), STVmb(b)%val(r))

      STVmb(b)%set= ballot(1:r,b) 
      STVmb(b)%lev= ballot2(1:r,b) 
      STVmb(b)%val(1)= wtb(b);  STVmb(b)%val(2:)= 0

      cn= ballot(1,b);  vote0(cn)= vote0(cn) + wtb(b)
    End do Ballot_loop0

    If (pr_out >= 1) Call Out ("Initial STV candidate cluster sizes, unordered", vote)
      
    vote= vote0
    ne= 0      ! # currently elected
    nd= 0      ! # currenlty deleted
      
    STV_cycle : Do iv= 1,nc
        
      Call List_of_true (Elect == 0, nh,key)  ! nh = # hopeful candidates remaining
      If (nh == 0) Exit STV_cycle  ! Done

!     Compute the mean ratings * vote over the positive hopeful candidates

      rate= 0
      Rate_loop : Do b= 1,nb
        p= STVmb(b)%p
        Call List_of_true (Elect(STVmb(b)%set(:p)) == 0, m,lh)

        If (m > 0) then
          lst(:m)= STVmb(b)%set(lh(:m))
          rate(lst(:m))= rate(lst(:m)) + pt_val(STVmb(b)%lev(lh(:m))) * wtb(b)
        End if
      End do Rate_loop 

      rate= rate / np;  rateV= rate * vote

      If (nh > 1) then  ! Sort the hopeful candidates from largest to smallest
                        ! based on the 'rateV' scaling
        rt_ord(:nh)= rateV(key(:nh))
        Call Sort (.false., rt_ord(:nh), lh(:nh))
        key(:nh)= key(lh(:nh))

        vt_ord(:nh)= vote(key(:nh))
        tot_hope= Sum(vt_ord(:nh))  ! Total vote for the hopefuls
      Else
        cn= key(1);  rt_ord(1)= rateV(cn)
        vt_ord(1)= vote(cn);  tot_hope= vote(cn)
      End if

      last= key(nh);  quota(iv)= tot_hope / (np1 - ne)  ! Droop formula for the quota

      pt(:nh)= vt_ord(:nh) - quota(iv)

      j= First_true(pt(:nh) >= 0);  Add= j > 0

      If (nd >= ndf .and. .not.Add) then
        j= Maxloc(pt(:nh), 1);  Add= .true.
      End if

      If (Add) then
        top= key(j);  prior= vote(top);  surplus= prior - quota(iv)
      End if

      If (pr_out >= 1) then
        Call Out ("Hopeful candidates sorted by decreasing rate-vote value", key(:nh))
        Call Out ("vs vote size", vt_ord(:nh))
        Call Out ("With total hopeful weight",tot_hope, "and quota",quota(iv))
      End if

      If (.not.Add) then ! Eliminate the bottom rate-vote candidate 'last' = key(nh)
        Call Out ("Delete bottom hopeful candidate", last, ln=1)

        m= 0;  l= 0;  retain= 0;  swt= 0;  vt_elc= 0

        Elect1= Elect;  If (ne < np) Where(Elect > 0) Elect1= -10  ! Distribute only to next hopeful
        
        Elim_loop : Do b= 1,nb
          il= First_true(STVmb(b)%set == last);  If (il < 1) Cycle Elim_loop  ! 'last' not ranked
          trans_wt= STVmb(b)%val(il);     If (trans_wt <= 0) Cycle Elim_loop  ! No weight to transfer or retain

          p= STVmb(b)%p;  j= First_true(Elect1(STVmb(b)%set(:p)) >= 0 .and. STVmb(b)%set(:p) /= last)

          If (j > 0) then  ! Transfer weight & distribute to the 'j' candidate cluster
            STVmb(b)%val(il)= 0;  vote(last)= vote(last) - trans_wt

            STVmb(b)%val(j)= STVmb(b)%val(j) + trans_wt
            cn= STVmb(b)%set(j);  vote(cn)= vote(cn) + trans_wt
            If (Elect(cn) > 0) vt_elc= vt_elc + trans_wt

            l= l + 1;  swt= swt + trans_wt
          Else      ! No distribution: fully retain
            m= m + 1;  retain= retain + trans_wt
          End if
        End do Elim_loop

        nd= nd + 1;  Elect(last)= -nd;  elected(nc1 - nd)= last

        Where(Abs(vote) <= eps) vote= 0
        vt_ind= vote(last);  total= sum(vote)

        independents= independents + vt_ind

        n_trans(iv,1)= l;  n_trans(iv,2)= m;  n_trans(iv,3)= -last
        w_trans(iv,1)= vt_elc;  w_trans(iv,2)= swt;  w_trans(iv,3)= retain 

        If (pr_out > 1) then
          Call Out ("# delete ballot transfers, # retained ballots, elected candidate",n_trans(iv,:))
          Call Out ("Weight transfered to electeds, total transfered, total fully retained",w_trans(iv,:))
          Call Out ("Remaining weight to independents",vt_ind, &
                    "current independents weight",independents)
          Call Out ("Current vote total",total)
          Call Out ("Current cluster sizes",vote)
        End if

      Else  ! ne < np .and. surplus >= 0, so elect the top candidate 'top'
        Call Out ("Add top hopeful candidate", top, ln=1)
        
!       Identify ballots with a transferable ballot weight for candidate 'top'

        ne= ne + 1;  m= 0;  l= 0;  retain= 0;  transferable= 0
        Elect1= Elect;  If (ne < np) Where(Elect > 0) Elect1= -10  ! Distribute only to next hopeful
          
        Add_loop0 : Do b= 1,nb
          STVmb(b)%n= 0;   it= First_true(STVmb(b)%set == top)
          STVmb(b)%r= it;  If (it < 1) Cycle Add_loop0  ! 'top' not ranked
          trans_wt= STVmb(b)%val(it);  If (trans_wt <= 0) Cycle Add_loop0  ! No weight to transfer or retain

          p= STVmb(b)%p;  j= First_true(Elect1(STVmb(b)%set(:p)) >= 0 .and. STVmb(b)%set(:p) /= top)
          STVmb(b)%n= j

          If (j > 0) then  ! Transferable to a positive rating level
            l= l + 1;  transferable= transferable + STVmb(b)%val(it)
          Else      ! No distribution: fully retain
            m= m + 1;  retain= retain + STVmb(b)%val(it)
          End if
        End do Add_loop0

        If (pr_out > 1) then
          Call Out ("Transferable weight",transferable, "from # possible ballot transfers",l)
        End if

 !      Transfer excess ballot weight
         
        If (transferable > min_trans) then
          transfer_fac= Min(surplus / transferable, 1.0)
          l= 0;  swt= 0;  vt_elc= 0
          
          Add_loop : Do b= 1,nb
            If (STVmb(b)%n <= 0) Cycle Add_loop

!           Subtract from the 'top' candidate cluster

            it= STVmb(b)%r;  trans_wt= transfer_fac * STVmb(b)%val(it)
            STVmb(b)%val(it)= STVmb(b)%val(it) - trans_wt
            vote(top)= vote(top) - trans_wt

!           Distribute to the candidate cluster specified by 'j'

            j= STVmb(b)%n;  STVmb(b)%val(j)= STVmb(b)%val(j) + trans_wt
            cn= STVmb(b)%set(j);  vote(cn)= vote(cn) + trans_wt
            If (Elect(cn) > 0) vt_elc= vt_elc + trans_wt

            l= l + 1;  swt= swt + trans_wt
          End do Add_loop

          Elect(top)= ne;  elected(ne)= top
          n_trans(iv,1)= l;       n_trans(iv,2)= m;    n_trans(iv,3)= top
          w_trans(iv,1)= vt_elc;  w_trans(iv,2)= swt;  w_trans(iv,3)= retain

          vt= prior - vote(top);  total= Sum(vote)  

          If (pr_out >= 1) then
            Call Out ("# ballot transfers, # retained ballots, elected candidate",n_trans(iv,:))
            Call Out ("Weight transfered to electeds, total transfered, total fully retained",w_trans(iv,:))
            Call Out ("top candidate vote decrease",vt, "total vote",total)
            Call Out ("for current cluster sizes",vote)
          End if
        Else  ! Elect without transfering
          Elect(top)= ne;  elected(ne)= top
          n_trans(iv,2)= m;  n_trans(iv,3)= top;  w_trans(iv,3)= retain

          If (pr_out >= 1) Call Out ("No transfers: STV elected candidate",top, &
                                     "retains vote",vote(top), ln=1)
        End if
      End if
    End do STV_cycle

    If (pr_out >= 0.5) then
      Call Out (1,"STV2 # ballot transfers, # retained ballots, elected candidate",n_trans)
      Call Out (1,"Weight transfered to electeds, total transfered, total fully retained",w_trans)
      Call Out ("Final cluster sizes",vote)

      vt_elc= Sum(vote(elected(:np)));  vt= Sum(vote)
      Call Out ("Final elected total",vt_elc, "vote total",vt)
    End if

!   Final STVmb output: Redo ballots so that the electeds 
!   come first in each ballo. Otherwise maintaining the original
!   ranking order.

    Call Final_STVmb (mr,nc, elected,Elect, STVmb, vote1)

  End Subroutine STV2


  Subroutine DTV (nc, wtb,ballot, elected, DTVmb)

!   DTV = "Distributed Transferable Vote" is a generalization of 
!   STV = "Single Transferable Vote". It's purpose to produce better
!   clustering results.
!   
!   The membership of a ballot 'b' with weight w(b) and candidate rankings 
!   (c1,...,cn) in a candidate cluster Ci is initially determined by the 
!   rank 'i' of ci in 'b', using a decreasing fractional series in the rank,
!   such as a geometric series:  Let g = a fraction such as 1/2 with 
!   dis(i,n) = g**i / (g**1 +...+ g**n) so that Sum(dis(:,n)) = 1. 
!   Then the membership w(Ci) of Ci is increased: w(Ci)= w(Ci) + w(ci,b)
!   where w(ci,b) = dis(i,n)*w(b). That is, the membership 
!   is distributed over the candidate clusters of the ranked candidates 
!   instead of being focused only on the top ranked cluster C1, as in STV.

!   When a candidate ci is deleted because w(Ci) is the smallest, then
!   its membership w(ci,b) in a ballot 'b' is redistributed to the 'k'
!   subsequent unelected, or hopeful, candidate clusters by the same 
!   fractional series:
!   Let j = 1,...,k be the indices of these subsequent hopefuls in 
!   ballot 'b', yielding the transfer w(Cj)= w(Cj) + dis(j,k)*w(ci,b) 
!   where w(ci,b) = current membership of 'b' in Ci.

!   When a candidate ci is elected because w(Ci) is the largest and 
!   exceeds the Droop formula, then, to reduce w(Ci) to the Droop value,
!   a fraction 'fr' of its membership w(ci,b) in a ballot 'b' is 
!   redistributed to the subsequent 'k' hopeful candidate clusters 
!   by the same fractional series:
!   Let j = 1,...,k be the indices of these hopefuls in this ballot, 
!   yielding the transfer w(Cj)= w(Cj) + dis(j,k)*fr*w(ci,b) where
!   w(ci,b) = current membership of 'b' in Ci. 

!   That is, fr = (surplus weight) / (transferable weight)

!   where the Droop formula is 
!        Quota = (total weight of the hopeful candidate clusters) / 
!                (1 + # candidates remaining to be elected)

!   Note that when a candidate cluster is eliminated, there may be 
!   ballots which rank that candidate with a positive weight but lack 
!   any subsequent hopefuls. In this case the weight is transferred 
!   to the cluster of "independents". Similarly, when a candidate is 
!   elected a ballot which weights that candidate may lack subsequent 
!   hopefuls, so its membership fraction is transferred to the independents.

! Input:
    Integer,  Intent(in) :: nc            ! # candidates
    Real,     Intent(in) :: wtb(:)        ! (nb) Ballot weights, summing to 'np'
    Integer,  Intent(in) :: ballot(0:,:)  ! (0:mr,nb) Candidates ranked by each ballot
                                          ! in the order ranked, with ballot(0,b)= # ranked 
                                          ! by ballot 'b' up to mr

! Output:
    Integer,        Intent(out) :: elected(:)   ! (nc) List of elected candidates in preferential order

    Type(Set_list), Intent(out) :: DTVmb(0:)    ! (0:nb) Ballot ranking structure for clustering
                                                ! Ballot sets (0): 
                                                !   %n        = nc
                                                !   %p        = np
                                                !   %svl      = elected vote = Sum(%val(set(:np)))
                                                !   %smb      = independents vote = Sum(%val(set(np+1:)))
                                                !   %set(nc)  = candidates in elected ordering, with
                                                !               the first 'np' elected
                                                !   %lev(nc)  = candidate clusters by decreasing cluster vote
                                                !   %val(nc)  = candidate cluster votes, not reordered

                                                ! Ballots (1:nb):
                                                !   %n     = # elected candidates
                                                !   %p     = # positively rated candidates
                                                !   %r     = # ranked candidates
                                                !   %set(r)= 'n' elected candidates, then unelected,
                                                !            in original ranking order
                                                !   %val(r)= corresponding membership weights 
                                                !   %svl   = elected membership weight= Sum(%val(:n))
                                                !   %smb   = total membership weight  = Sum(%val)
! Local:
    Real    :: vote0(nc)     ! Size of initial candidate clusters
    Real    :: vote(nc)      ! Size of candidate clusters at current STV cycle
    Real    :: vt_ord(nc)    ! Hopeful candidate clusters, ordered by decreasing size
    Real    :: quota(nc)     ! Droop formula quotas by STV cycle

    Real    :: w_trans(nc,3) ! Records data by STV cycle
                             ! (:,1) = weight transfered to electeds, (:,2) = total weight transfered
                             ! (:,3) = weight of fully retained ballots
    Integer :: n_trans(nc,3) ! Records data by STV cycle
                             ! (:,1) = # ballot transfers, (:,2) = # ballots retained
                             ! (:,3) = deleted (neg) or elected (pos) candidate

    Integer :: Elect(nc)     ! (i) = i if ith elected, 0 if hopeful, -i if ith eliminated
    Integer :: Elect1(nc)    ! = 'Elect' except < 0 where 'Elect' > 0 but not all elected

    Real    :: vote1(nc)
    Integer :: lh(nc), lst(nc), key(nc)

    Real, Allocatable :: dis(:,:)  ! (mr,mr) Decreasing point distribution, summing to 1 over (1:n,:) 
                                   !         for each # ranked candidates (:,n)

    Real, Parameter :: min_trans= 0.001, Dis_parm= 0.50
    Logical :: Add

    Real    :: retain, surplus, trans_wt, transfer_fac, transferable, independents
    Real    :: vt, gtr, inc, swt, sm_pt, tot_hope, prior, total, vt_elc, vt_ind
    Integer :: mr, n1, n2, nb, nh, ne, np, nr, nc1, np1, top, nd, neh, last
    Integer :: b, i, j, k, l, m, n, p, r, c1, cn, h1, ih, il, it, iv

    nb= Size(wtb);  np= Nint(Sum(wtb));  mr= Ubound(ballot,1)
    np1= np + 1;  nc1= nc + 1

    Call Out ("DTV algorithm for # electeds",np, "out of # candidates",nc,ln=1)

!   Geometric series for distribution of ballot membership to candidate clusters

    Allocate(dis(mr,mr));  dis= -1

    If (Dis_parm > 1) then  ! Geometric distribution
      dis(1,mr)= 1
      Do i= 2,mr
        dis(i,mr)= dis(i-1,mr) / Dis_parm
      End do
    Else                    ! Arithmetic distribution
      dis(1,mr)= 1;  inc= 1
      Do i= 2,mr
        inc= inc + Dis_parm
        dis(i,mr)= dis(i-1,mr) + inc
      End do
      dis(:,mr)= dis(mr:1:-1,mr)
    End if
  
    Do n= 1,mr
      sm_pt= Sum(dis(:n,mr))
      dis(:n,n)= dis(:n,mr) / sm_pt
    End do

    Call Out ("Enter 'DTV',with distribution parameter",Dis_parm, ln=1)
    If (pr_out >= 1) Call Out (-1,"DTV distribution fractions for each # ranked",dis)

!   Initial ballot memberships in candidate clusters

    elected= 0;  vote0= 0;  n_trans= 0;  w_trans= 0;  Elect= 0;  independents= 0;   quota= 0

    Ballot_loop0 : Do b= 1,nb
      DTVmb(b)%smb= wtb(b);  r= ballot(0,b)

      If (Associated(DTVmb(b)%set)) DeAllocate(DTVmb(b)%set, DTVmb(b)%val)
      Allocate(DTVmb(b)%set(r), DTVmb(b)%val(r))

      DTVmb(b)%set= ballot(1:r,b)
      DTVmb(b)%val= dis(:r,r) * wtb(b)

      vote0(DTVmb(b)%set)= vote0(DTVmb(b)%set) + DTVmb(b)%val
    End do Ballot_loop0

    If (pr_out >= 1) Call Out ("Initial DTV candidate cluster sizes, unordered", vote0)
      
    vote= vote0
    key= 'ID'  ! Hopeful candidates, to be sorted by weight
    ne= 0      ! # currently elected
    nd= 0      ! # currenlty deleted
      
    DTV_cycle : Do iv= 1,nc
        
      Call List_of_true (Elect == 0, nh,key)  ! nh = # hopeful candidates remaining
      If (nh == 0) Exit DTV_cycle  ! Done

      If (nh > 1) then  ! Sort the hopeful candidates from largest to smallest
        vt_ord(:nh)= vote(key(:nh))
        Call Sort (.false., vt_ord(:nh), lh(:nh))

        key(:nh)= key(lh(:nh));  tot_hope= Sum(vt_ord(:nh))
      Else
        vt_ord(1)= vote(key(1));  tot_hope= vt_ord(1)
      End if

      top= key(1);  last= key(nh)
      quota(iv)= tot_hope / (np1 - ne)  ! Droop formula for the quota
      
!     Elect candidates whose cluster size exceeds the quota, 
!     else remove the smallest hopeful candidate cluster

      prior= vote(top);  surplus= prior - quota(iv)
      Add= ne < np .and. surplus >= 0

      If (pr_out > 1) then
        Call Out ("Hopeful candidates sorted by decreasing cluster size", key(:nh))
        Call Out ("from updated decreasing cluster sizes", vt_ord(:nh))

        Call Out ("With total hopeful weight",tot_hope, "and quota",quota(iv))
        Call Out ("yields the surplus",surplus)
      End if

      If (.not.Add) then ! Eliminate the lowest value hopeful candidate 'last'
        Call Out ("Next eliminate the smallest candidate cluster", last, ln=1)

        m= 0;  l= 0;  retain= 0;  swt= 0;  vt_elc= 0
        Elect1= Elect;  If (ne < np) Where(Elect > 0) Elect1= -10  ! Distribute only to next hopeful
        
        Elim_loop : Do b= 1,nb

          il= First_true(DTVmb(b)%set == last);  If (il < 1) Cycle Elim_loop  ! 'last' not ranked
          trans_wt= DTVmb(b)%val(il);  If (trans_wt <= 0) Cycle Elim_loop     ! No weight to transfer

          Call List_of_true (Elect1(DTVmb(b)%set) >= 0 .and. DTVmb(b)%set /= last, ih,lh)

          If (ih > 0) then  ! Transfer weight & distribute to the 'ih' candidate clusters

            DTVmb(b)%val(il)= 0;  vote(last)= vote(last) - trans_wt

            Do i= 1,ih
              k= lh(i);  gtr= dis(i,ih) * trans_wt
              DTVmb(b)%val(k)= DTVmb(b)%val(k) + gtr
              cn= DTVmb(b)%set(k);  vote(cn)= vote(cn) + gtr
              If (Elect(cn) > 0) vt_elc= vt_elc + gtr
            End do

            l= l + 1;  swt= swt + trans_wt
          Else      ! No distribution: fully retain
            m= m + 1;  retain= retain + trans_wt
          End if
        End do Elim_loop

        nd= nd + 1;  Elect(last)= -nd;  elected(nc1 - nd)= last
        vt_ind= vote(last);  total= Sum(vote)

        independents= independents + vt_ind

        n_trans(iv,1)= l;       n_trans(iv,2)= m;     n_trans(iv,3)= -last
        w_trans(iv,1)= vt_elc;  w_trans(iv,2)= swt ;  w_trans(iv,3)= retain 
      
        If (pr_out > 1) then
          Call Out ("# delete ballot transfers, # retained ballots, elected candidate",n_trans(iv,:))
          Call Out ("Weight transfered to electeds, total transfered, total fully retained",w_trans(iv,:))
          Call Out ("Remaining weight to independents",vt_ind, &
                    "current independents weight",independents)
          Call Out ("Current vote total",total)
          Call Out ("Current cluster sizes",vote)
        End if

      Else  ! ne < np .and. surplus > 0, so elect the top candidate 'top'
        Call Out ("Next add the largest candidate cluster", top, ln=1)
        
!       Identify ballots with a transferable ballot weight for candidate 'top'

        ne= ne + 1;  m= 0;  l= 0;  retain= 0;  transferable= 0
        Elect1= Elect;  If (ne < np) Where(Elect > 0) Elect1= -10  ! Distribute only to next hopeful
          
        Add_loop0 : Do b= 1,nb
          DTVmb(b)%n= 0;     it= First_true(DTVmb(b)%set == top)
            DTVmb(b)%r= it;  If (it < 1) Cycle Add_loop0  ! 'top' not ranked
          trans_wt= DTVmb(b)%val(it);  If (trans_wt <= 0) Cycle Add_loop0     ! No weight to transfer

          j= First_true(Elect1(DTVmb(b)%set) >= 0 .and. DTVmb(b)%set /= top)
          DTVmb(b)%n= j

          If (j > 0) then
            l= l + 1;  transferable= transferable + DTVmb(b)%val(it)
          Else      ! No distribution: fully retain
            m= m + 1;  retain= retain + DTVmb(b)%val(it)
          End if
        End do Add_loop0

        If (pr_out >= 1) Call Out ("Transferable weight",transferable, &
                                  "from # possible ballot transfers",l, ln=1)
 
!       Transfer excess ballot weight
         
        If (transferable > min_trans) then
          transfer_fac= Min(surplus / transferable, 1.0)
          l= 0;  swt= 0;  vt_elc= 0
          
          Add_loop : Do b= 1,nb
            If (DTVmb(b)%n <= 0) Cycle Add_loop

            it= DTVmb(b)%r;  j= DTVmb(b)%n
            Call List_of_true (Elect1(DTVmb(b)%set(j:)) >= 0 .and. DTVmb(b)%set(j:) /= top, ih,lh)
            lh(:ih)= (j-1) + lh(:ih)

!           Subtract the 'trans_wt' weight from the 'top' candidate cluster

            trans_wt= transfer_fac * DTVmb(b)%val(it)
            DTVmb(b)%val(it)= DTVmb(b)%val(it) - trans_wt

            vote(top)= vote(top) - trans_wt
            l= l + 1;  swt= swt + trans_wt

!           Distrubte the 'trans_wt' weight to the 'lh' list

            Do i= 1,ih
              k= lh(i);  gtr= dis(i,ih) * trans_wt
              DTVmb(b)%val(k)= DTVmb(b)%val(k) + gtr
              cn= DTVmb(b)%set(k);  vote(cn)= vote(cn) + gtr
              If (Elect(cn) > 0) vt_elc= vt_elc + gtr
            End do
          End do Add_loop

          Elect(top)= ne;  elected(ne)= top
          n_trans(iv,1)= l;       n_trans(iv,2)= m;    n_trans(iv,3)= top
          w_trans(iv,1)= vt_elc;  w_trans(iv,2)= swt;  w_trans(iv,3)= retain

          vt= prior - vote(top);  total= Sum(vote)  

          If (pr_out >= 1) then
            Call Out ("# add ballot transfers, # retained ballots, elected candidate",n_trans(iv,:))
            Call Out ("Weight transfered to electeds, total transfered, total fully retained",w_trans(iv,:))
            Call Out ("top candidate vote decrease",vt, "total vote",total)
            Call Out ("for current cluster sizes",vote)
          End if
        Else  ! Elect without transfering
          Elect(top)= ne;  elected(ne)= top
          n_trans(iv,2)= m;  n_trans(iv,3)= top;  w_trans(iv,3)= retain

          If (pr_out >= 1)  Call Out ("No transfers: DTV elected candidate",top, &
                                      "retains vote",vote(top), ln=1)
        End if
      End if
    End do DTV_cycle

    If (pr_out >= 0.5) then
      Call Out (1,"DTV # ballot transfers, # retained ballots, elected candidate",n_trans)
      Call Out (1,"Weight transfered to electeds, total transfered, total fully retained",w_trans)
      Call Out ("Final cluster sizes",vote)

      vt_elc= Sum(vote(elected(:np)));  vt= Sum(vote) 
      Call Out ("Final elected total",vt_elc, "vote total",vt)
    End if

!   Final DTVmb output: Redo ballots so that the electeds 
!   come first in each ballo. Otherwise maintaining the original
!   ranking order.

    Call Final_STVmb (mr,nc, elected,Elect, DTVmb, vote1)

  End Subroutine DTV
   

  Subroutine DTV2 (nc,pt_val, wtb,ballot,ballot2,elected, DTVmb)

!   DTV = "Distributed Transferable Vote" is a generalization of 
!   STV = "Single Transferable Vote". It's purpose to produce better
!   clustering results.
!   
!   The membership of a ballot 'b' with weight w(b) and candidate rankings 
!   (c1,...,cn) in a candidate cluster Ci is initially determined by the 
!   rank 'i' of ci in 'b', using a decreasing fractional series in the rank.

!   Then the membership w(Ci) of Ci is increased: w(Ci)= w(Ci) + w(ci,b)
!   where w(ci,b) = wt(i)*w(b). That is, the membership 
!   is distributed over the candidate clusters of the ranked candidates 
!   instead of being focused only on the top ranked cluster C1, as in STV.

!   When a candidate ci is deleted because w(Ci) is the smallest, then
!   its membership w(ci,b) in a ballot 'b' is redistributed to the 'k'
!   subsequent unelected, or hopeful, candidate clusters by the same 
!   fractional series:
!   Let j = 1,...,k be the indices of these subsequent hopefuls in 
!   ballot 'b', yielding the transfer w(Cj)= w(Cj) + wt(j)*w(ci,b) 
!   where w(ci,b) = current membership of 'b' in Ci.

!   When a candidate ci is elected because w(Ci) is the largest and 
!   exceeds the Droop formula, then, to reduce w(Ci) to the Droop value,
!   a fraction 'fr' of its membership w(ci,b) in a ballot 'b' is 
!   redistributed to the subsequent 'k' hopeful candidate clusters 
!   by the same fractional series:
!   Let j = 1,...,k be the indices of these hopefuls in this ballot, 
!   yielding the transfer w(Cj)= w(Cj) + wt(j)*fr*w(ci,b) where
!   w(ci,b) = current membership of 'b' in Ci. 

!   That is, fr = (surplus weight) / (transferable weight)

!   where the Droop formula is 
!        Quota = (total weight of the hopeful candidate clusters) / 
!                (1 + # candidates remaining to be elected)

!   Note that when a candidate cluster is eliminated, there may be 
!   ballots which rank that candidate with a positive weight but lack 
!   any subsequent hopefuls. In this case the weight is transferred 
!   to the cluster of "independents". Similarly, when a candidate is 
!   elected a ballot which weights that candidate may lack subsequent 
!   hopefuls, so its membership fraction is transferred to the independents.

! Input:
    Integer, Intent(in) :: nc              ! # candidates
    Real,    Intent(in) :: pt_val(:)       ! (mt) Rating or ranking points corresponding to 
                                           !      rating or ranking levels.
    Real,    Intent(in) :: wtb(:)          ! (nb) Ballot weights, summing to 'np'
    Integer, Intent(in) :: ballot(0:,:)    ! (0:mr,nb). (0,:) = 'r' candidates in preferential order
    Integer, Intent(in) :: ballot2(0:,:)   ! (0:mr,nb). corresponding ranking or rating levels
                                           !   with (0,:) = 'n' positively rated candidates

! Output:
    Integer,        Intent(out) :: elected(:)   ! (nc) List of elected candidates in preferential order

    Type(Set_list), Intent(out) :: DTVmb(0:)    ! (0:nb) Ballot ranking structure for clustering (on output)
                                                ! Ballot sets (0): 
                                                !   %n        = nc
                                                !   %p        = np
                                                !   %svl      = elected vote = Sum(%val(set(:np)))
                                                !   %smb      = independents vote = Sum(%val(set(np+1:)))
                                                !   %set(nc)  = candidates in elected ordering, with
                                                !               the first 'np' elected
                                                !   %lev(nc)  = candidate clusters by decreasing cluster vote
                                                !   %val(nc)  = candidate cluster votes, not reordered

                                                ! Ballots (1:nb):
                                                !   %n     = # elected candidates
                                                !   %p     = # positively rated candidates
                                                !   %r     = # ranked candidates
                                                !   %set(r)= 'n' elected candidates, then unelected,
                                                !            in original ranking order
                                                !   %val(r)= corresponding membership weights 
                                                !   %svl   = elected membership weight= Sum(%val(:n))
                                                !   %smb   = total membership weight  = Sum(%val)
! Local:
    Real    :: vote0(nc)     ! Size of the initial candidate clusters
    Real    :: vote(nc)      ! Size of each candidate cluster at the 
                             !  current stage of the DTV cycle
    Real    :: vt_ord(nc)    ! (nh) Hopeful candidate clusters, by decreasing size
    Real    :: quota(nc)     ! Droop formula quotas by DTV cycle

    Real    :: rate(nc)      ! Mean rating vector of the hopefuls
    Real    :: rateV(nc)     ! Current 'rate' * 'vote'
    Real    :: rt_ord(nc)    ! (nh) Hopeful candidate clusters, by decreasing 'rateV'

    Real    :: w_trans(nc,3) ! Records data by STV cycle
                             ! (:,1) = weight transfered to electeds, (:,2) = total weight transfered
                             ! (:,3) = weight of fully retained ballots
    Integer :: n_trans(nc,3) ! Records data by STV cycle
                             ! (:,1) = # ballot transfers, (:,2) = # ballots retained
                             ! (:,3) = deleted (neg) or elected (pos) candidate

    Integer :: Elect(nc)     ! = i if ith elected, 0 if hopeful, -i if ith eliminated
    Integer :: Elect1(nc)    ! (nc) = 'Elect' except < 0 where 'Elect' > 0 but not all elected

    Real    :: pt(nc), wt(nc), vote1(nc)
    Integer :: lh(nc), lst(nc), key(nc)

    Real, Parameter :: min_trans= 0.001, eps= 0.0001
    Logical :: Add

    Real    :: retain, surplus, trans_wt, transfer_fac, transferable, independents
    Real    :: vt, tr, swt, sm_pt, tot_hope, prior, total, vt_elc, vt_ind
    Integer :: mr, mt, n1, n2, nb, nd, nh, ne, np, nr, nc1, ndf, np1, top, last
    Integer :: b, i, j, k, l, m, n, p, r, c1, cn, h1, ih, il, it, iv

    np= Nint(Sum(wtb));  np1= np + 1;  nc1= nc + 1;  mr= Ubound(ballot,1)
    nb= Size(wtb);  mt= Size(pt_val);  ndf= nc - np

    Call Out ("DTV2 algorithm for # electeds",np, "out of # candidates",nc,ln=1)

!   Initial ballot memberships in candidate clusters

    elected= 0;  vote0= 0;  Elect= 0;  independents= 0
    quota= 0;  n_trans= 0;  w_trans= 0
      
    Ballot_loop0 : Do b= 1,nb
      DTVmb(b)%smb= wtb(b);  r= ballot(0,b);  p= ballot2(0,b)
      DTVmb(b)%p= p
      
      If (Associated(DTVmb(b)%set)) DeAllocate(DTVmb(b)%set)
      If (Associated(DTVmb(b)%lev)) DeAllocate(DTVmb(b)%lev)
      If (Associated(DTVmb(b)%val)) DeAllocate(DTVmb(b)%val)
      Allocate(DTVmb(b)%set(r),DTVmb(b)%lev(r), DTVmb(b)%val(r))

      DTVmb(b)%set= ballot(1:r,b)  
      DTVmb(b)%lev= ballot2(1:r,b)

      wt(:p)= pt_val(DTVmb(b)%lev(:p)) 
      sm_pt= Sum(wt(:p));  wt(:p)= wt(:p) / sm_pt

      DTVmb(b)%val(:p)= wt(:p) * wtb(b);  DTVmb(b)%val(p+1:)= 0
      vote0(DTVmb(b)%set(:p))= vote0(DTVmb(b)%set(:p)) + DTVmb(b)%val(:p)

    End do Ballot_loop0

    If (pr_out >= 1) then
      total= Sum(vote0)
      Call Out ("Initial DTV vote total", total, ln=1)
      Call Out ("Corresponding candidate cluster votes, unordered", vote0)
    End if
      
    vote= vote0
    key= 'ID'  ! Hopeful candidates, to be sorted by weight
    ne= 0      ! # currently elected
    nd= 0      ! # currenlty deleted
      
    DTV_cycle : Do iv= 1,nc
        
      Call List_of_true (Elect == 0, nh,key)  ! nh = # hopeful candidates remaining
      If (nh == 0) Exit DTV_cycle  ! Done

!     Compute the mean ratings * vote over the positive hopeful candidates

      rate= 0
      Rate_loop : Do b= 1,nb
        p= DTVmb(b)%p;  Call List_of_true (Elect(DTVmb(b)%set(:p)) == 0, m,lh)

        If (m > 0) then
          lst(:m)= DTVmb(b)%set(lh(:m))
          rate(lst(:m))= rate(lst(:m)) + pt_val(DTVmb(b)%lev(lh(:m))) * wtb(b)
        End if
      End do Rate_loop 

      rate= rate / np;  rateV= rate * vote

      If (nh > 1) then  ! Sort the hopeful candidates from largest to smallest
                        ! based on the 'rateV' scaling
        rt_ord(:nh)= rateV(key(:nh))
        Call Sort (.false., rt_ord(:nh), lh(:nh))
        key(:nh)= key(lh(:nh))

        vt_ord(:nh)= vote(key(:nh))
        tot_hope= Sum(vt_ord(:nh))  ! Total vote for the hopefuls
      Else
        cn= key(1);  rt_ord(1)= rateV(cn)
        vt_ord(1)= vote(cn);  tot_hope= vote(cn)
      End if
      
      last= key(nh);  quota(iv)= tot_hope / (np1 - ne)  ! Droop formula for the quota

      pt(:nh)= vt_ord(:nh) - quota(iv)

      j= First_true(pt(:nh) >= 0);  Add= j > 0

      If (nd >= ndf .and. .not.Add) then
        j= Maxloc(pt(:nh), 1);  Add= .true.
      End if

      If (Add) then
        top= key(j);  prior= vote(top);  surplus= prior - quota(iv)
      End if

      If (pr_out >= 1) then
        Call Out ("Hopeful candidates sorted by decreasing rate-vote value", key(:nh))
        Call Out ("vs vote size", vt_ord(:nh))
        Call Out ("With total hopeful weight",tot_hope, "and quota",quota(iv))
      End if

      If (.not.Add) then ! Eliminate the bottom rate-vote candidate 'last' = key(nh)

        Call Out ("Next delete the bottom rate-vote candidate", last, ln=1)
        m= 0;  l= 0;  retain= 0;  swt= 0;  vt_elc= 0

        Elect1= Elect;  If (ne < np) Where (Elect > 0) Elect1= -10  ! Distribute only to next hopeful
        
        Elim_loop : Do b= 1,nb

          il= First_true(DTVmb(b)%set == last);  If (il < 1) Cycle Elim_loop  ! 'last' not ranked
          trans_wt= DTVmb(b)%val(il);  If (trans_wt <= 0) Cycle Elim_loop     ! No weight to transfer

          p= DTVmb(b)%p
          Call List_of_true (Elect1(DTVmb(b)%set(:p)) >= 0 .and. DTVmb(b)%set(:p) /= last, ih,lh)

          If (ih > 0) then  ! Transfer weight & distribute to the 'ih' candidate clusters

            DTVmb(b)%val(il)= 0;  vote(last)= vote(last) - trans_wt

            wt(:ih)= pt_val(DTVmb(b)%lev(lh(:ih))) 
            sm_pt= Sum(wt(:ih));  wt(:ih)= wt(:ih) / sm_pt

            Do i= 1,ih
              k= lh(i);  tr= wt(i) * trans_wt
              DTVmb(b)%val(k)= DTVmb(b)%val(k) + tr
              cn= DTVmb(b)%set(k);  vote(cn)= vote(cn) + tr
              If (Elect(cn) > 0) vt_elc= vt_elc + tr
            End do

            l= l + 1;  swt= swt + trans_wt
          Else      ! No distribution: fully retain
            m= m + 1;  retain= retain + trans_wt
          End if
        End do Elim_loop

        nd= nd + 1;  Elect(last)= -nd;  elected(nc1 - nd)= last
        vt_ind= vote(last);  total= Sum(vote)

        independents= independents + vt_ind

        n_trans(iv,1)= l;  n_trans(iv,2)= m;  n_trans(iv,3)= -last
        w_trans(iv,1)= vt_elc;  w_trans(iv,2)= swt;  w_trans(iv,3)= retain 
      
        If (pr_out > 1) then
          Call Out ("# delete ballot transfers, # retained ballots, elected candidate",n_trans(iv,:))
          Call Out ("Weight transfered to electeds, total transfered, total fully retained",w_trans(iv,:))
          Call Out ("Remaining weight to independents",vt_ind, &
                    "current independents weight",independents)
          Call Out ("Current vote total",total)
          Call Out ("Current cluster sizes",vote)
        End if

      Else  ! ne < np .and. surplus > 0, so elect the top candidate 'top'
        Call Out ("Next add the top rate-vote candidate", top, ln=1)
        ne= ne + 1;  m= 0;  l= 0;  retain= 0;  transferable= 0

!       Identify ballots with a transferable ballot weight for candidate 'top'

        Elect1= Elect;  If (ne < np) Where (Elect1 > 0) Elect1= -2
          
        Add_loop0 : Do b= 1,nb
          DTVmb(b)%n= 0;     it= First_true(DTVmb(b)%set == top)
            DTVmb(b)%r= it;  If (it < 1) Cycle Add_loop0  ! 'top' not ranked

          trans_wt= DTVmb(b)%val(it);  If (trans_wt <= 0) Cycle Add_loop0 ! No weight to transfer

          p= DTVmb(b)%p;  j= First_true(Elect1(DTVmb(b)%set(:p)) >= 0 .and. DTVmb(b)%set(:p) /= top)
          DTVmb(b)%n= j

          If (j > 0) then  ! Transferable to a positive rating level
            l= l + 1;  transferable= transferable + DTVmb(b)%val(it)
          Else      ! No distribution: fully retain
            m= m + 1;  retain= retain + DTVmb(b)%val(it)
          End if
        End do Add_loop0

        If (pr_out > 1) then
          Call Out ("Transferable weight",transferable, "from # possible ballot transfers",l)
        End if
 
!      Transfer excess ballot weight
         
        If (transferable > min_trans) then
          transfer_fac= Min(surplus / transferable, 1.0)
          l= 0;  swt= 0;  vt_elc= 0
          
          Add_loop : Do b= 1,nb
            If (DTVmb(b)%n == 0) Cycle Add_loop

!           Subtract the 'trans_wt' weight from the 'top' candidate cluster

            it= DTVmb(b)%r;  j= DTVmb(b)%n
            trans_wt= transfer_fac * DTVmb(b)%val(it)
            DTVmb(b)%val(it)= DTVmb(b)%val(it) - trans_wt

            vote(top)= vote(top) - trans_wt
            l= l + 1;  swt= swt + trans_wt

!           Distribute the 'trans_wt' weight to the list 'lh'
!           of eligible positively ranked candidates

            p= DTVmb(b)%p
            Call List_of_true (Elect1(DTVmb(b)%set(j:p)) >= 0 .and. &
                               DTVmb(b)%set(j:p) /= top, ih,lh)
            lh(:ih)= (j-1) + lh(:ih)

            wt(:ih)= pt_val(DTVmb(b)%lev(lh(:ih)))
            sm_pt= Sum(wt(:ih));  wt(:ih)= wt(:ih) / sm_pt

            If (Any(wt(:ih) <= 0.0)) then
              Call Out ("Error in 'DTV2': sum of weights not pos",wt(:ih))
            End if

            Do i= 1,ih
              k= lh(i);  tr= wt(i) * trans_wt
              DTVmb(b)%val(k)= DTVmb(b)%val(k) + tr
              cn= DTVmb(b)%set(k);  vote(cn)= vote(cn) + tr
              If (Elect(cn) > 0) vt_elc= vt_elc + tr
            End do
          End do Add_loop

          Elect(top)= ne;  elected(ne)= top
          n_trans(iv,1)= l;       n_trans(iv,2)= m;    n_trans(iv,3)= top
          w_trans(iv,1)= vt_elc;  w_trans(iv,2)= swt;  w_trans(iv,3)= retain

          vt= prior - vote(top);  total= Sum(vote)  

          If (pr_out >= 1) then
            Call Out ("# ballot transfers, # retained ballots, elected candidate",n_trans(iv,:))
            Call Out ("Weight transfered to electeds, total transfered, total fully retained",w_trans(iv,:))
            Call Out ("Current cluster sizes",vote)
            Call Out ("top candidate vote decrease",vt, "total vote", total)
          End if
        Else  ! Elect without transfering
          Elect(top)= ne;  elected(ne)= top
          n_trans(iv,2)= m;  n_trans(iv,3)= top;  w_trans(iv,3)= retain

          If (pr_out >= 1) then
            Call Out ("No transfers: DTV elected candidate",top, &
                      "retains vote",vote(top), ln=1)
          End if
        End if
      End if
    End do DTV_cycle

    If (pr_out >= 0.5) then
      Call Out (1,"DTV2 # ballot transfers, # retained ballots, elected candidate",n_trans)
      Call Out (1,"Weight transfered to electeds, total transfered, total fully retained",w_trans)
      Call Out ("Final cluster sizes",vote)

      vt_elc= Sum(vote(elected(:np)));  vt= Sum(vote) 
      Call Out ("Final elected total",vt_elc, "vote total",vt)
    End if

!   Final DTVmb output: Redo ballots so that the electeds 
!   come first in each ballo. Otherwise maintaining the original
!   ranking order.

    Call Final_STVmb (mr,nc, elected,Elect, DTVmb, vote1)

  End Subroutine DTV2
   
  Subroutine WSTV (nc, wtb,ballot, elected, STVmb)

!   This algorithm is the Warren version of STV.

!   The concept is to keep all features of STV except for
!   different quantities to be retained or transferred to cover 
!   the surplus for each ballot representing a newly elected 
!   candidate cluster 'cn'. This works by retaining the full 
!   ballot weight in the 'cn' cluster, with 0 transfer, for  
!   ballots of weight smaller than a computed value 'rtw',
!   then retaining 'rtw' from each of the ballots of weight
!   while transferring the rest: 'wt' - 'rtw'. 

!   The formula is
!       rtw = (Q - W(m)) / (N - m) where the 'cn'
!   candidate cluster contains N ballots with at least
!   partial membership, with m = # of these with weight
!   < rtw and W(m) = the sum of their weights. Thus the
!   sum of the retained ballots weights = the quota Q 
!   exactly, with each ballot weight - its retained weight
!   transferred according to the usual surplus rules.

!   The effect is that ballots of smaller ballot weight in a
!   candidate cluster 'cn' are will not be transferred at all, 
!   retaininng more diversity for that cluster, whereas 
!   those ballots of largest weight will spread more of 
!   the support for 'cn' to allied candidate clusters.
!   The net result will increase the dominance of the larger
!   voting blocks


! Input:
    Integer,         Intent(in) :: nc            ! # candidates
    Real,            Intent(in) :: wtb(:)        ! (nb) Ballot weights, summing to 'np'
    Integer,         Intent(in) :: ballot(0:,:)  ! (1:mr,nb) Candidates ranked by each ballot
                                                 ! in the order ranked, with ballot(0,b)= # ranked 
                                                 ! by ballot 'b' up to mr
! Output:
    Integer,        Intent(out) :: elected(:)   ! (nc) List of elected candidates in preferential order

    Type(Set_list), Intent(out) :: STVmb(0:)    ! (0:nb) Ballot ranking structure for clustering
                                                ! Ballot sets (0): 
                                                !   %n        = nc
                                                !   %p        = np
                                                !   %svl      = elected vote = Sum(%val(set(:np)))
                                                !   %smb      = independents vote = Sum(%val(set(np+1:)))
                                                !   %set(nc)  = candidates in elected ordering, with
                                                !               the first 'np' elected
                                                !   %lev(nc)  = candidate clusters by decreasing cluster vote
                                                !   %val(nc)  = candidate cluster votes, not reordered

                                                ! Ballots (1:nb):
                                                !   %n     = # elected candidates
                                                !   %p     = # positively rated candidates
                                                !   %r     = # ranked candidates
                                                !   %set(r)= 'n' elected candidates, then unelected,
                                                !            in original ranking order
                                                !   %val(r)= corresponding membership weights 
                                                !   %svl   = elected membership weight= Sum(%val(:n))
                                                !   %smb   = total membership weight  = Sum(%val)
! Local:
    Real    :: vote0(nc)     ! Size of initial candidate clusters
    Real    :: vote(nc)      ! Size of candidate clusters at current STV cycle
    Real    :: vt_ord(nc)    ! Hopeful candidate clusters, ordered by decreasing size
    Real    :: quota(nc)     ! Droop formula quotas by STV cycle

    Real    :: w_trans(nc,3) ! Records data by STV cycle
                             ! (:,1) = weight transfered to electeds, (:,2) = total weight transfered
                             ! (:,3) = weight of fully retained ballots
    Integer :: n_trans(nc,3) ! Records data by STV cycle
                             ! (:,1) = # ballot transfers, (:,2) = # ballots retained
                             ! (:,3) = deleted (neg) or elected (pos) candidate

    Integer :: Elect(nc)     ! (i) = i if ith elected, 0 if hopeful, -i if ith eliminated
    Integer :: Elect1(nc)    ! = 'Elect' except < 0 where 'Elect' > 0 but not all elected

    Real,    Allocatable :: wt(:,:)  ! (n,2)  (:,1) = increasing transferable ballots weights 
                                     !        (:,2) =
    Integer, Allocatable :: kw(:,:)  ! (n,2)  (:,1) = corresponding ballot indices
                                     !        (:,2) = reodrdering key for wt(:,1)
    Real     :: rt, ratio, trans, retain(3)
    Integer  :: nsp, ntr, nrt(2)

    Real     :: vote1(nc)
    Integer  :: lh(nc), key(nc)

    Real, Parameter :: min_trans= 0.001, eps= 0.0001
    Logical :: Add, Update
    Real    :: surplus, trans_wt, tot_hope, transferable, independents
    Real    :: vt, vt_elc, vt_ind, prior, total 

    Integer :: n1, nb, nh, ne, np, nd, nc1, np1, top, last
    Integer :: c1, cn, h1, il, it, iv, mr
    Integer :: b, h, i, j, k, l, m, n, p, r

    nb= Size(wtb);  np= Nint(Sum(wtb));  mr= Ubound(ballot,1)
    np1= np + 1;  nc1= nc + 1
    
    Call Out ("Warren STV method for # electeds",np, "out of # candidates",nc,ln=1)

!   Initial ballot memberships in candidate clusters

    elected= 0;  vote0= 0;  Elect= 0;  independents= 0
    quota= 0;   n_trans= 0;  w_trans= 0

    Ballot_loop0 : Do b= 1,nb
      STVmb(b)%smb= wtb(b);  r= ballot(0,b)

      If (Associated(STVmb(b)%set)) DeAllocate(STVmb(b)%set)
      If (Associated(STVmb(b)%val)) DeAllocate(STVmb(b)%val)
      Allocate(STVmb(b)%set(r), STVmb(b)%val(r))

      STVmb(b)%set= ballot(1:r,b);  
      STVmb(b)%val(1)= wtb(b);  STVmb(b)%val(2:)= 0

      cn= ballot(1,b);  vote0(cn)= vote0(cn) + wtb(b)
    End do Ballot_loop0

    If (pr_out >= 1) Call Out ("Initial STV candidate cluster sizes, unordered", vote0)
      
    vote= vote0
    key= 'ID'  ! Hopeful candidates, to be sorted by weight
    ne= 0      ! # currently elected
    nd= 0      ! # currenlty deleted
      
    STV_cycle : Do iv= 1,nc
        
      Call List_of_true (Elect == 0, nh,key)  ! nh = # hopeful candidates remaining
      If (nh == 0) Exit STV_cycle  ! Done

      If (nh > 1) then  ! Sort the hopeful candidates from largest to smallest
        vt_ord(:nh)= vote(key(:nh))
        Call Sort (.false., vt_ord(:nh), lh(:nh))

        key(:nh)= key(lh(:nh));  tot_hope= Sum(vt_ord(:nh))
      Else
        vt_ord(1)= vote(key(1));  tot_hope= vt_ord(1)
      End if

      top= key(1);  last= key(nh)

      quota(iv)= tot_hope / (np1 - ne)  ! Droop formula for the quota

      prior= vote(top);  surplus= prior - quota(iv)  ! Weight to be transferred if > 0

      Add= ne < np .and. surplus >= 0

      If (pr_out >= 1) then
        Call Out ("Hopeful candidates sorted by decreasing vote sizes", key(:nh))
        Call Out ("for vote sizes", vt_ord(:nh))
        Call Out ("With total hopeful weight",tot_hope, "quota",quota(iv), ln=1)
        Call Out ("& surplus",surplus)
      End if

      If (.not.Add) then ! Eliminate the bottom rate-vote candidate 'last' = key(nh)
        Call Out ("Delete bottom hopeful candidate", last, ln=1)

        m= 0;  l= 0;  rt= 0;  trans= 0;  vt_elc= 0
        Elect1= Elect;  If (ne < np) Where(Elect > 0) Elect1= -10  ! Distribute only to next hopeful
        
        Elim_loop : Do b= 1,nb
          il= First_true(STVmb(b)%set == last);  If (il < 1) Cycle Elim_loop  ! 'last' not ranked
          trans_wt= STVmb(b)%val(il);     If (trans_wt <= 0) Cycle Elim_loop  ! No weight to transfer or retain

          j= First_true(Elect1(STVmb(b)%set) >= 0 .and. STVmb(b)%set /= last) ! Next eligible candidate

          If (j > 0) then  ! Transfer weight & distribute to the 'j' candidate cluster
            STVmb(b)%val(il)= 0;  vote(last)= vote(last) - trans_wt

            STVmb(b)%val(j)= STVmb(b)%val(j) + trans_wt
            cn= STVmb(b)%set(j);  vote(cn)= vote(cn) + trans_wt
            If (Elect(STVmb(b)%set(j)) > 0) vt_elc= vt_elc + trans_wt

            l= l + 1;  trans= trans + trans_wt
          Else  ! No distribution: fully retain
            m= m + 1;  rt= rt + trans_wt
          End if
        End do Elim_loop

        nd= nd + 1;  Elect(last)= -nd;  elected(nc1 - nd)= last
        vt_ind= vote(last);  total= sum(vote)

        independents= independents + vt_ind

        n_trans(iv,1)= l;  n_trans(iv,2)= m;  n_trans(iv,3)= -last 
        w_trans(iv,1)= vt_elc;  w_trans(iv,2)= trans;  w_trans(iv,3)= rt
      
        If (pr_out > 1) then
          Call Out ("# delete ballot transfers, # retained ballots, elected candidate",n_trans(iv,:))
          Call Out ("Weight transfered to electeds, total transfered, total fully retained",w_trans(iv,:))
          Call Out ("Remaining weight to independents",vt_ind, &
                    "current independents weight",independents)
          Call Out ("Current vote total",total)
          Call Out ("Current cluster sizes",vote)
        End if

      Else  ! ne < np .and. surplus > 0, so elect the top candidate 'top'
        Call Out ("Add top hopeful candidate", top, ln=1)
        
!       Identify ballots with partially transferable vs fully retained 
!       ballot weight for candidate 'top'

        ne= ne + 1;  m= 0;  l= 0;  rt= 0;  transferable= 0

        Elect1= Elect;  If (ne < np) Where(Elect > 0) Elect1= -10  ! Distribute only to next hopeful
          
        Add_loop0 : Do b= 1,nb
          STVmb(b)%n= 0;   it= First_true(STVmb(b)%set == top)
          STVmb(b)%r= it;  If (it < 1) Cycle Add_loop0  ! 'top' not ranked

          trans_wt= STVmb(b)%val(it);  If (trans_wt <= 0) Cycle Add_loop0 ! No weight to transfer or retain

          j= First_true(Elect1(STVmb(b)%set) >= 0 .and. STVmb(b)%set /= top)
          STVmb(b)%n= j

          If (j > 0) then
            l= l + 1;  transferable= transferable + STVmb(b)%val(it)
          Else
            m= m + 1;  rt= rt + STVmb(b)%val(it)
          End if
        End do Add_loop0

        ntr= l;  nrt= 0;  retain= 0;  nrt(1)= m;  retain(1)= rt
        nsp= 0;  vt_elc= 0

        If (ntr > 0) then
          Call RetainW (top, ntr,quota(iv), Elect, STVmb,vote, nrt,retain, nsp,trans,vt_elc)
        
          Elect(top)= ne;  elected(ne)= top;  m= nrt(1) + nrt(2);  rt= Sum(retain)
          n_trans(iv,1)= nsp;       n_trans(iv,2)= m;    n_trans(iv,3)= top
          w_trans(iv,1)= vt_elc;  w_trans(iv,2)= trans;  w_trans(iv,3)= rt

          vt= prior - vote(top);  total= Sum(vote);  ratio= nsp / (nrt(2)+eps)  

          If (pr_out >= 1) then
            Call Out ("Quota",quota(iv), "vs actual retained weight",rt, ln=1)
            Call Out ("Surplus",surplus, "vs actually transfered",trans)
            Call Out ("vs actual top candidate decrease",vt, "total vote",total)
            Call Out ("Ratio of # transfers to # retentions in Warren",ratio)

            Call Out ("# transfers to surplus, # fully retained ballots, elected candidate",n_trans(iv,:))
            Call Out ("Weight transfered to electeds, total transfered, total fully retained",w_trans(iv,:))
            Call Out ("for current cluster sizes",vote)
          End if
        Else  ! Elect without transfering
          Elect(top)= ne;  elected(ne)= top
          n_trans(iv,2)= m;  n_trans(iv,3)= top;  w_trans(iv,3)= rt

          Call Out ("No transfers: WSTV elected candidate",top, &
                    "retains vote",vote(top), ln=1)
        End if
      End if
    End do STV_cycle

    If (pr_out >= 0.5) then
      Call Out (1,"WSTV # ballot transfers, # retained ballots, elected candidate",n_trans)
      Call Out (1,"Weight transfered to electeds, total transfered, total fully retained",w_trans)
      Call Out ("Final cluster sizes",vote)

      vt_elc= Sum(vote(elected(:np)));  vt= Sum(vote) 
      Call Out ("Final elected total",vt_elc, "vote total",vt)
    End if

!   Final STVmb output: Redo ballots so that the electeds 
!   come first in each ballo. Otherwise maintaining the original
!   ranking order.

    Call Final_STVmb (mr,nc, elected,Elect, STVmb, vote1)

  End Subroutine WSTV


  Subroutine WSTV2 (nc,pt_val, wtb,ballot,ballot2, elected, STVmb)

!   This algorithm is an extension of STV, once the final 
!   candidate is elected, to produce better clustering. The
!   concept is to redistribute subsequent transfers
!   from the first 'hopeful' candidate remaining in a ballot 
!   to the first 'prior elected' or 'hopeful' candidate.

!   That is, it is no longer necessary to transfer the 
!   maximum possible in excess of the Droop quota 
!   to subsequent hopefuls in order to qualify one 
!   of them for election.
!
!   The effect is to redistribute some of the ballot weight 
!   to the candidate clusters of the prior electeds, so that
!   thier size may exceed the quota. Thus all subsequent 
!   hopeful candidate clusters, as they are deleted, will 
!   also redistribute some of their weight to the electeds
!   diminishing size of the pseuo-cluster of independents.

! Input:
    Integer,  Intent(in) :: nc              ! # to be elected, # candidates, # ballots
    Real,     Intent(in) :: pt_val(:)       ! (mt) Rating or ranking points corresponding to 
                                            !      rating or ranking levels.
    Real,     Intent(in) :: wtb(:)          ! (nb) Ballot weights, summing to 'np'
    Integer,  Intent(in) :: ballot(0:,:)    ! (0:mr,nb). (0,:) = 'n' candidates in preferential order
    Integer,  Intent(in) :: ballot2(0:,:)   ! (mr,nb). corresponding ranking or rating levels
! Output:
    Integer,        Intent(out) :: elected(:)   ! (nc) List of elected candidates in preferential order

    Type(Set_list), Intent(out) :: STVmb(0:)    ! (0:nb) Ballot ranking structure for clustering
                                                ! Ballot sets (0): 
                                                !   %n        = nc
                                                !   %p        = np
                                                !   %svl      = elected vote = Sum(%val(set(:np)))
                                                !   %smb      = independents vote = Sum(%val(set(np+1:)))
                                                !   %set(nc)  = candidates in elected ordering, with
                                                !               the first 'np' elected
                                                !   %lev(nc)  = candidate clusters by decreasing cluster vote
                                                !   %val(nc)  = candidate cluster votes, not reordered

                                                ! Ballots (1:nb):
                                                !   %n     = # elected candidates
                                                !   %p     = # positively rated candidates
                                                !   %r     = # ranked candidates
                                                !   %set(r)= 'n' elected candidates, then unelected,
                                                !            in original ranking order
                                                !   %val(r)= corresponding membership weights 
                                                !   %svl   = elected membership weight= Sum(%val(:n))
                                                !   %smb   = total membership weight  = Sum(%val)
! Local:
    Real    :: vote0(nc)     ! Original size of each candidate cluster
    Real    :: vote(nc)      ! Currrent size of each candidate cluster
    Real    :: vt_ord(nc)    ! (nh) Hopeful candidate clusters, by decreasing size
    Real    :: quota(nc)     ! Droop formula quotas by STV cycle

    Real    :: rate(nc)      ! Mean rating vector of the hopefuls
    Real    :: rateV(nc)     ! Current 'rate' * 'vote'
    Real    :: rt_ord(nc)    ! (nh) Hopeful candidate clusters, by decreasing 'rateV'

    Real    :: w_trans(nc,3) ! Records data by STV cycle
                             ! (:,1) = weight transfered to electeds, (:,2) = total weight transfered
                             ! (:,3) = weight of fully retained ballots
    Integer :: n_trans(nc,3) ! Records data by STV cycle
                             ! (:,1) = # ballot transfers, (:,2) = # ballots retained
                             ! (:,3) = deleted (neg) or elected (pos) candidate

    Integer :: Elect(nc)   ! = i if ith elected, 0 if hopeful, -i if ith eliminated
    Integer :: Elect1(nc)  ! 'Elect' except < 0 where 'Elect' > 0 but not all elected

    Real,    Allocatable :: wt(:,:)  ! (n,2)  (:,1) = increasing transferable ballots weights 
                                     !        (:,2) =
    Integer, Allocatable :: kw(:,:)  ! (n,2)  (:,1) = corresponding ballot indices
                                     !        (:,2) = reodrdering key for wt(:,1)
    Real     :: rt, ratio, trans, retain(3)
    Integer  :: nsp, ntr, nrt(2)

    Real    :: pt(nc), vote1(nc)
    Integer :: lh(nc), lst(nc), key(nc)

    Real, Parameter :: min_trans= 0.001, eps= 0.0001
    Logical :: Add, Update

    Real    :: surplus, trans_wt, transfer_fac, transferable, independents
    Real    :: vt, sm_pt, prior, total, vt_elc, vt_ind, tot_hope

    Integer :: mr, mt, nb, nd, n1, nh, ne, np, nr, nc1, ndf, np1, top, last
    Integer :: c1, cn, h1, il, it, iv
    Integer :: b, j, l, m, n, p, r

    mr= Ubound(ballot,1);  mt= Size(pt_val);  np= Nint(Sum(wtb))
    np1= np + 1;  nc1= nc + 1;  nb= Size(wtb);  ndf= nc - np

    Call Out ("Warren STV2 algorithm for # electeds",np, "out of # candidates",nc,ln=1)

!   Initial ballot memberships in candidate clusters

    elected= 0;  vote0= 0;  Elect= 0;  independents= 0
    quota= 0;  n_trans= 0;  w_trans= 0

    Ballot_loop0 : Do b= 1,nb
      STVmb(b)%smb= wtb(b);  r= ballot(0,b);  p= ballot2(0,b)
      STVmb(b)%p= p

      If (Associated(STVmb(b)%set)) DeAllocate(STVmb(b)%set)
      If (Associated(STVmb(b)%lev)) DeAllocate(STVmb(b)%lev)
      If (Associated(STVmb(b)%val)) DeAllocate(STVmb(b)%val)
      Allocate(STVmb(b)%set(r), STVmb(b)%lev(r), STVmb(b)%val(r))

      STVmb(b)%set= ballot(1:r,b) 
      STVmb(b)%lev= ballot2(1:r,b) 
      STVmb(b)%val(1)= wtb(b);  STVmb(b)%val(2:)= 0

      cn= ballot(1,b);  vote0(cn)= vote0(cn) + wtb(b)
    End do Ballot_loop0

    If (pr_out >= 1) Call Out ("Initial STV candidate cluster sizes, unordered", vote)
      
    vote= vote0
    ne= 0      ! # currently elected
    nd= 0      ! # currenlty deleted
      
    STV_cycle : Do iv= 1,nc
        
      Call List_of_true (Elect == 0, nh,key)  ! nh = # hopeful candidates remaining
      If (nh == 0) Exit STV_cycle  ! Done

!     Compute the mean ratings * vote over the positive hopeful candidates

      rate= 0
      Rate_loop : Do b= 1,nb
        p= STVmb(b)%p
        Call List_of_true (Elect(STVmb(b)%set(:p)) == 0, m,lh)

        If (m > 0) then
          lst(:m)= STVmb(b)%set(lh(:m))
          rate(lst(:m))= rate(lst(:m)) + pt_val(STVmb(b)%lev(lh(:m))) * wtb(b)
        End if
      End do Rate_loop 

      rate= rate / np;  rateV= rate * vote

      If (nh > 1) then  ! Sort the hopeful candidates from largest to smallest
                        ! based on the 'rateV' scaling
        rt_ord(:nh)= rateV(key(:nh))
        Call Sort (.false., rt_ord(:nh), lh(:nh))
        key(:nh)= key(lh(:nh))

        vt_ord(:nh)= vote(key(:nh))
        tot_hope= Sum(vt_ord(:nh))  ! Total vote for the hopefuls
      Else
        cn= key(1);  rt_ord(1)= rateV(cn)
        vt_ord(1)= vote(cn);  tot_hope= vote(cn)
      End if

      last= key(nh);  quota(iv)= tot_hope / (np1 - ne)  ! Droop formula for the quota

      pt(:nh)= vt_ord(:nh) - quota(iv)

      j= First_true(pt(:nh) >= 0);  Add= j > 0

      If (nd >= ndf .and. .not.Add) then
        j= Maxloc(pt(:nh), 1);  Add= .true.
      End if

      If (Add) then
        top= key(j);  prior= vote(top);  surplus= prior - quota(iv)
      End if

      If (pr_out >= 1) then
        Call Out ("Hopeful candidates sorted by decreasing rate-vote value", key(:nh))
        Call Out ("vs vote size", vt_ord(:nh))
        Call Out ("With total hopeful weight",tot_hope, "and quota",quota(iv))
      End if

      If (.not.Add) then ! Eliminate the bottom rate-vote candidate 'last' = key(nh)
        Call Out ("Delete bottom hopeful candidate", last, ln=1)

        m= 0;  l= 0;  rt= 0;  trans= 0;  vt_elc= 0

        Elect1= Elect;  If (ne < np) Where(Elect > 0) Elect1= -10  ! Distribute only to next hopeful
        
        Elim_loop : Do b= 1,nb
          il= First_true(STVmb(b)%set == last);  If (il < 1) Cycle Elim_loop  ! 'last' not ranked
          trans_wt= STVmb(b)%val(il);     If (trans_wt <= 0) Cycle Elim_loop  ! No weight to transfer or retain

          p= STVmb(b)%p;  j= First_true(Elect1(STVmb(b)%set(:p)) >= 0 .and. STVmb(b)%set(:p) /= last)

          If (j > 0) then  ! Transfer weight & distribute to the 'j' candidate cluster
            STVmb(b)%val(il)= 0;  vote(last)= vote(last) - trans_wt

            STVmb(b)%val(j)= STVmb(b)%val(j) + trans_wt
            cn= STVmb(b)%set(j);  vote(cn)= vote(cn) + trans_wt
            If (Elect(cn) > 0) vt_elc= vt_elc + trans_wt

            l= l + 1;  trans= trans + trans_wt
          Else      ! No distribution: fully retain
            m= m + 1;  rt= rt + trans_wt
          End if
        End do Elim_loop

        nd= nd + 1;  Elect(last)= -nd;  elected(nc1 - nd)= last

        Where(Abs(vote) <= eps) vote= 0
        vt_ind= vote(last);  total= sum(vote)

        independents= independents + vt_ind

        n_trans(iv,1)= l;  n_trans(iv,2)= m;  n_trans(iv,3)= -last
        w_trans(iv,1)= vt_elc;  w_trans(iv,2)= trans;  w_trans(iv,3)= rt 

        If (pr_out > 1) then
          Call Out ("# delete ballot transfers, # retained ballots, elected candidate",n_trans(iv,:))
          Call Out ("Weight transfered to electeds, total transfered, total fully retained",w_trans(iv,:))
          Call Out ("Remaining weight to independents",vt_ind, &
                    "current independents weight",independents)
          Call Out ("Current vote total",total)
          Call Out ("Current cluster sizes",vote)
        End if

      Else  ! ne < np .and. surplus >= 0, so elect the top candidate 'top'
        Call Out ("Add top hopeful candidate", top, ln=1)
        
!       Identify ballots with partially transferable vs fully retained 
!       ballot weight for candidate 'top'

        ne= ne + 1;  m= 0;  l= 0;  rt= 0;  transferable= 0
        Elect1= Elect;  If (ne < np) Where(Elect > 0) Elect1= -10  ! Distribute only to next hopeful
          
        Add_loop0 : Do b= 1,nb
          STVmb(b)%n= 0;   it= First_true(STVmb(b)%set == top)
          STVmb(b)%r= it;  If (it < 1) Cycle Add_loop0  ! 'top' not ranked

          trans_wt= STVmb(b)%val(it);  If (trans_wt <= 0) Cycle Add_loop0  ! No weight to transfer or retain

          p= STVmb(b)%p;  j= First_true(Elect1(STVmb(b)%set(:p)) >= 0 .and. STVmb(b)%set(:p) /= top)
          STVmb(b)%n= j

          If (j > 0) then  ! Transferable to a positive rating level
            l= l + 1;  transferable= transferable + STVmb(b)%val(it)
          Else
            m= m + 1;  rt= rt + STVmb(b)%val(it)
          End if
        End do Add_loop0

        ntr= l;  nrt= 0;  retain= 0;  nrt(1)= m;  retain(1)= rt
        nsp= 0;  vt_elc= 0

        If (ntr > 0) then
          Call RetainW (top, ntr,quota(iv), Elect, STVmb,vote, nrt,retain, nsp,trans,vt_elc)
        
          Elect(top)= ne;  elected(ne)= top;  m= nrt(1) + nrt(2);  rt= Sum(retain)
          n_trans(iv,1)= nsp;       n_trans(iv,2)= m;    n_trans(iv,3)= top
          w_trans(iv,1)= vt_elc;  w_trans(iv,2)= trans;  w_trans(iv,3)= rt          

          vt= prior - vote(top);  total= Sum(vote);  ratio= nsp / (nrt(2)+eps)  

          If (pr_out >= 1) then
            Call Out ("Quota",quota(iv), "vs actual retained weight",rt, ln=1)
            Call Out ("Surplus",surplus, "vs actually transfered",trans)
            Call Out ("vs actual top candidate decrease",vt, "total vote",total)
            Call Out ("Ratio of # transfers to # retentions in Warren",ratio)
            Call Out ("# transfers to surplus, # fully retained ballots, elected candidate",n_trans(iv,:))
            Call Out ("Weight transfered to electeds, total transfered, total fully retained",w_trans(iv,:))
            Call Out ("for current cluster sizes",vote)
          End if
        Else  ! Elect without transfering
          Elect(top)= ne;  elected(ne)= top
          n_trans(iv,2)= m;  n_trans(iv,3)= top;  w_trans(iv,3)= rt

          If (pr_out >= 1) Call Out ("No transfers: STV elected candidate",top, &
                                     "retains vote",vote(top), ln=1)
        End if
      End if
    End do STV_cycle

    If (pr_out >= 0.5) then
      Call Out (1,"WSTV2 # ballot transfers, # retained ballots, elected candidate",n_trans)
      Call Out (1,"Weight transfered to electeds, total transfered, total fully retained",w_trans)
      Call Out ("Final cluster sizes",vote)

      vt_elc= Sum(vote(elected(:np)));  vt= Sum(vote)
      Call Out ("Final elected total",vt_elc, "vote total",vt)
    End if

!   Final STVmb output: Redo ballots so that the electeds 
!   come first in each ballo. Otherwise maintaining the original
!   ranking order.

    Call Final_STVmb (mr,nc, elected,Elect, STVmb, vote1)

  End Subroutine WSTV2

  Subroutine WDTV (nc, wtb,ballot, elected, DTVmb)

!   WDTV = "Distributed Transferable Vote" is a generalization of WSTV,
!   the Warren version of STV, = "Single Transferable Vote". It's purpose 
!   to produce better clustering results.
!   
!   The membership of a ballot 'b' with weight w(b) and candidate rankings 
!   (c1,...,cn) in a candidate cluster Ci is initially determined by the 
!   rank 'i' of ci in 'b', using a decreasing fractional series in the rank,
!   such as a geometric series:  Let g = a fraction such as 1/2 with 
!   dis(i,n) = g**i / (g**1 +...+ g**n) so that Sum(dis(:,n)) = 1. 
!   Then the membership w(Ci) of Ci is increased: w(Ci)= w(Ci) + w(ci,b)
!   where w(ci,b) = dis(i,n)*w(b). That is, the membership 
!   is distributed over the candidate clusters of the ranked candidates 
!   instead of being focused only on the top ranked cluster C1, as in STV.

!   When a candidate ci is deleted because w(Ci) is the smallest, then
!   its membership w(ci,b) in a ballot 'b' is redistributed to the 'k'
!   subsequent unelected, or hopeful, candidate clusters by the same 
!   fractional series:
!   Let j = 1,...,k be the indices of these subsequent hopefuls in 
!   ballot 'b', yielding the transfer w(Cj)= w(Cj) + dis(j,k)*w(ci,b) 
!   where w(ci,b) = current membership of 'b' in Ci.

!   When a candidate ci is elected because w(Ci) is the largest and 
!   exceeds the Droop formula, then, to reduce w(Ci) to the Droop value,
!   a fraction 'fr' of its membership w(ci,b) in a ballot 'b' is 
!   redistributed to the subsequent 'k' hopeful candidate clusters 
!   by the same fractional series:
!   Let j = 1,...,k be the indices of these hopefuls in this ballot, 
!   yielding the transfer w(Cj)= w(Cj) + dis(j,k)*fr*w(ci,b) where
!   w(ci,b) = current membership of 'b' in Ci. 

!   That is, fr = (surplus weight) / (transferable weight)

!   where the Droop formula is 
!        Quota = (total weight of the hopeful candidate clusters) / 
!                (1 + # candidates remaining to be elected)

!   Note that when a candidate cluster is eliminated, there may be 
!   ballots which rank that candidate with a positive weight but lack 
!   any subsequent hopefuls. In this case the weight is transferred 
!   to the cluster of "independents". Similarly, when a candidate is 
!   elected a ballot which weights that candidate may lack subsequent 
!   hopefuls, so its membership fraction is transferred to the independents.

! Input:
    Integer,  Intent(in) :: nc            ! # candidates
    Real,     Intent(in) :: wtb(:)        ! (nb) Ballot weights, summing to 'np'
    Integer,  Intent(in) :: ballot(0:,:)  ! (0:mr,nb) Candidates ranked by each ballot
                                          ! in the order ranked, with ballot(0,b)= # ranked 
                                          ! by ballot 'b' up to mr

! Output:
    Integer,        Intent(out) :: elected(:)   ! (nc) List of elected candidates in preferential order

    Type(Set_list), Intent(out) :: DTVmb(0:)    ! (0:nb) Ballot ranking structure for clustering
                                                ! Ballot sets (0): 
                                                !   %n        = nc
                                                !   %p        = np
                                                !   %svl      = elected vote = Sum(%val(set(:np)))
                                                !   %smb      = independents vote = Sum(%val(set(np+1:)))
                                                !   %set(nc)  = candidates in elected ordering, with
                                                !               the first 'np' elected
                                                !   %lev(nc)  = candidate clusters by decreasing cluster vote
                                                !   %val(nc)  = candidate cluster votes, not reordered

                                                ! Ballots (1:nb):
                                                !   %n     = # elected candidates
                                                !   %p     = # positively rated candidates
                                                !   %r     = # ranked candidates
                                                !   %set(r)= 'n' elected candidates, then unelected,
                                                !            in original ranking order
                                                !   %val(r)= corresponding membership weights 
                                                !   %svl   = elected membership weight= Sum(%val(:n))
                                                !   %smb   = total membership weight  = Sum(%val)
! Local:
    Real    :: vote0(nc)     ! Size of initial candidate clusters
    Real    :: vote(nc)      ! Size of candidate clusters at current STV cycle
    Real    :: vt_ord(nc)    ! Hopeful candidate clusters, ordered by decreasing size
    Real    :: quota(nc)     ! Droop formula quotas by STV cycle

    Real    :: w_trans(nc,3) ! Records data by STV cycle
                             ! (:,1) = weight transfered to electeds, (:,2) = total weight transfered
                             ! (:,3) = weight of fully retained ballots
    Integer :: n_trans(nc,3) ! Records data by STV cycle
                             ! (:,1) = # ballot transfers, (:,2) = # ballots retained
                             ! (:,3) = deleted (neg) or elected (pos) candidate

    Integer :: Elect(nc)     ! (i) = i if ith elected, 0 if hopeful, -i if ith eliminated
    Integer :: Elect1(nc)    ! = 'Elect' except < 0 where 'Elect' > 0 but not all elected

    Real,    Allocatable :: wt(:,:)  ! (n,2)  (:,1) = increasing transferable ballots weights 
                                     !        (:,2) =
    Integer, Allocatable :: kw(:,:)  ! (n,2)  (:,1) = corresponding ballot indices
                                     !        (:,2) = reodrdering key for wt(:,1)
    Real    :: rt, ratio, trans, retain(3)
    Integer :: nsp, ntr, nrt(2)

    Real    :: vote1(nc)
    Integer :: lh(nc), lst(nc), key(nc)

    Real, Allocatable :: dis(:,:)  ! (mr,mr) Decreasing point distribution, summing to 1 over (1:n,:) 
                                   !         for each # ranked candidates (:,n)

    Real, Parameter :: eps= 0.0001, min_trans= 0.001, Dis_parm= 0.50
    Logical :: Add

    Real    :: surplus, trans_wt, transfer_fac, transferable, independents
    Real    :: vt, gtr, inc, swt, sm_pt, tot_hope, prior, total, vt_elc, vt_ind
    Integer :: mr, n1, n2, nb, nh, ne, np, nr, nc1, np1, top, nd, neh, last
    Integer :: b, i, j, k, l, m, n, p, r, c1, cn, h1, ih, il, it, iv

    nb= Size(wtb);  np= Nint(Sum(wtb));  mr= Ubound(ballot,1)
    np1= np + 1;  nc1= nc + 1

    Call Out ("WDTV algorithm for # electeds",np, "out of # candidates",nc,ln=1)

!   Geometric series for distribution of ballot membership to candidate clusters

    Allocate(dis(mr,mr));  dis= -1

    If (Dis_parm > 1) then  ! Geometric distribution
      dis(1,mr)= 1
      Do i= 2,mr
        dis(i,mr)= dis(i-1,mr) / Dis_parm
      End do
    Else                    ! Arithmetic distribution
      dis(1,mr)= 1;  inc= 1
      Do i= 2,mr
        inc= inc + Dis_parm
        dis(i,mr)= dis(i-1,mr) + inc
      End do
      dis(:,mr)= dis(mr:1:-1,mr)
    End if
  
    Do n= 1,mr
      sm_pt= Sum(dis(:n,mr))
      dis(:n,n)= dis(:n,mr) / sm_pt
    End do

    Call Out ("Enter 'WDTV',with distribution parameter",Dis_parm, ln=1)
    If (pr_out >= 1) Call Out (-1,"WDTV distribution fractions for each # ranked",dis)

!   Initial ballot memberships in candidate clusters

    elected= 0;  vote0= 0;  n_trans= 0;  w_trans= 0;  Elect= 0;  independents= 0;   quota= 0

    Ballot_loop0 : Do b= 1,nb
      DTVmb(b)%smb= wtb(b);  r= ballot(0,b)

      If (Associated(DTVmb(b)%set)) DeAllocate(DTVmb(b)%set, DTVmb(b)%val)
      Allocate(DTVmb(b)%set(r), DTVmb(b)%val(r))

      DTVmb(b)%set= ballot(1:r,b)
      DTVmb(b)%val= dis(:r,r) * wtb(b)

      vote0(DTVmb(b)%set)= vote0(DTVmb(b)%set) + DTVmb(b)%val
    End do Ballot_loop0

    If (pr_out >= 1) Call Out ("Initial WDTV candidate cluster sizes, unordered", vote0)
      
    vote= vote0
    key= 'ID'  ! Hopeful candidates, to be sorted by weight
    ne= 0      ! # currently elected
    nd= 0      ! # currenlty deleted
      
    DTV_cycle : Do iv= 1,nc
        
      Call List_of_true (Elect == 0, nh,key)  ! nh = # hopeful candidates remaining
      If (nh == 0) Exit DTV_cycle  ! Done

      If (nh > 1) then  ! Sort the hopeful candidates from largest to smallest
        vt_ord(:nh)= vote(key(:nh))
        Call Sort (.false., vt_ord(:nh), lh(:nh))

        key(:nh)= key(lh(:nh));  tot_hope= Sum(vt_ord(:nh))
      Else
        vt_ord(1)= vote(key(1));  tot_hope= vt_ord(1)
      End if

      top= key(1);  last= key(nh)
      quota(iv)= tot_hope / (np1 - ne)  ! Droop formula for the quota
      
!     Elect candidates whose cluster size exceeds the quota, 
!     else remove the smallest hopeful candidate cluster

      prior= vote(top);  surplus= prior - quota(iv)
      Add= ne < np .and. surplus >= 0

      If (pr_out > 1) then
        Call Out ("Hopeful candidates sorted by decreasing cluster size", key(:nh))
        Call Out ("from updated decreasing cluster sizes", vt_ord(:nh))

        Call Out ("With total hopeful weight",tot_hope, "and quota",quota(iv))
        Call Out ("yields the surplus",surplus)
      End if

      If (.not.Add) then ! Eliminate the lowest value hopeful candidate 'last'
        Call Out ("Next eliminate the smallest candidate cluster", last, ln=1)

        m= 0;  l= 0;  rt= 0;  swt= 0;  vt_elc= 0
        Elect1= Elect;  If (ne < np) Where(Elect > 0) Elect1= -10  ! Distribute only to next hopeful
        
        Elim_loop : Do b= 1,nb

          il= First_true(DTVmb(b)%set == last);  If (il < 1) Cycle Elim_loop  ! 'last' not ranked
          trans_wt= DTVmb(b)%val(il);  If (trans_wt <= 0) Cycle Elim_loop     ! No weight to transfer

          Call List_of_true (Elect1(DTVmb(b)%set) >= 0 .and. DTVmb(b)%set /= last, ih,lh)

          If (ih > 0) then  ! Transfer weight & distribute to the 'ih' candidate clusters

            DTVmb(b)%val(il)= 0;  vote(last)= vote(last) - trans_wt

            Do i= 1,ih
              k= lh(i);  gtr= dis(i,ih) * trans_wt
              DTVmb(b)%val(k)= DTVmb(b)%val(k) + gtr
              cn= DTVmb(b)%set(k);  vote(cn)= vote(cn) + gtr
              If (Elect(cn) > 0) vt_elc= vt_elc + gtr
            End do

            l= l + 1;  swt= swt + trans_wt
          Else  ! No distribution: fully retain
            m= m + 1;  rt= rt + trans_wt
          End if
        End do Elim_loop

        nd= nd + 1;  Elect(last)= -nd;  elected(nc1 - nd)= last
        vt_ind= vote(last);  total= Sum(vote)

        independents= independents + vt_ind

        n_trans(iv,1)= l;       n_trans(iv,2)= m;     n_trans(iv,3)= -last
        w_trans(iv,1)= vt_elc;  w_trans(iv,2)= swt ;  w_trans(iv,3)= rt
      
        If (pr_out > 1) then
          Call Out ("# delete ballot transfers, # retained ballots, elected candidate",n_trans(iv,:))
          Call Out ("Weight transfered to electeds, total transfered, total fully retained",w_trans(iv,:))
          Call Out ("Remaining weight to independents",vt_ind, &
                    "current independents weight",independents)
          Call Out ("Current vote total",total)
          Call Out ("Current cluster sizes",vote)
        End if

      Else  ! ne < np .and. surplus > 0, so elect the top candidate 'top'
        Call Out ("Next add the largest candidate cluster", top, ln=1)
        
!       Identify ballots with partially transferable vs fully retained 
!       ballot weight for candidate 'top'

        ne= ne + 1;  m= 0;  l= 0;  rt= 0;  transferable= 0
        Elect1= Elect;  If (ne < np) Where(Elect > 0) Elect1= -10  ! Distribute only to next hopeful
          
        Add_loop0 : Do b= 1,nb
          DTVmb(b)%n= 0;     it= First_true(DTVmb(b)%set == top)
            DTVmb(b)%r= it;  If (it < 1) Cycle Add_loop0  ! 'top' not ranked
          trans_wt= DTVmb(b)%val(it);  If (trans_wt <= 0) Cycle Add_loop0     ! No weight to transfer

          j= First_true(Elect1(DTVmb(b)%set) >= 0 .and. DTVmb(b)%set /= top)
          DTVmb(b)%n= j

          If (j > 0) then  ! l = # ballots partially transferable
            l= l + 1;  transferable= transferable + DTVmb(b)%val(it)
          Else      ! No distribution: fully retain
            m= m + 1;  rt= rt + DTVmb(b)%val(it)
          End if
        End do Add_loop0

        ntr= l;  nrt= 0;  retain= 0;  nrt(1)= m;  retain(1)= rt
        nsp= 0;  vt_elc= 0

        If (ntr > 0) then
          Call RetainWD (top, ntr,quota(iv),dis, Elect,Elect1, DTVmb,vote, &
                         nrt,retain, nsp,trans,vt_elc)
        
          Elect(top)= ne;  elected(ne)= top;  m= nrt(1) + nrt(2);  rt= Sum(retain)
          n_trans(iv,1)= nsp;       n_trans(iv,2)= m;    n_trans(iv,3)= top
          w_trans(iv,1)= vt_elc;  w_trans(iv,2)= trans;  w_trans(iv,3)= rt

          vt= prior - vote(top);  total= Sum(vote);  ratio= nsp / (nrt(2)+eps)  

          If (pr_out >= 1) then
            Call Out ("Quota",quota(iv), "vs actual retained weight",rt, ln=1)
            Call Out ("Surplus",surplus, "vs actually transfered",trans)
            Call Out ("vs actual top candidate decrease",vt, "total vote",total)
            Call Out ("Ratio of # transfers to # retentions in Warren",ratio)

            Call Out ("# transfers to surplus, # fully retained ballots, elected candidate",n_trans(iv,:))
            Call Out ("Weight transfered to electeds, total transfered, total fully retained",w_trans(iv,:))
            Call Out ("for current cluster sizes",vote)
          End if
        Else  ! Elect without transfering
          Elect(top)= ne;  elected(ne)= top
          n_trans(iv,2)= m;  n_trans(iv,3)= top;  w_trans(iv,3)= rt

          Call Out ("No transfers: WDTV elected candidate",top, &
                    "retains vote",vote(top), ln=1)
        End if
      End if
    End do DTV_cycle

    If (pr_out >= 0.5) then
      Call Out (1,"WDTV # ballot transfers, # retained ballots, elected candidate",n_trans)
      Call Out (1,"Weight transfered to electeds, total transfered, total fully retained",w_trans)
      Call Out ("Final cluster sizes",vote)

      vt_elc= Sum(vote(elected(:np)));  vt= Sum(vote) 
      Call Out ("Final elected total",vt_elc, "vote total",vt)
    End if

!   Final DTVmb output: Redo ballots so that the electeds 
!   come first in each ballo. Otherwise maintaining the original
!   ranking order.

    Call Final_STVmb (mr,nc, elected,Elect, DTVmb, vote1)

  End Subroutine WDTV
   
  Subroutine WDTV2 (nc,pt_val, wtb,ballot,ballot2,elected, DTVmb)

!   WDTV2 = "Distributed Transferable Vote" is a generalization of WSTV2,
!   the Warren version of STV2, = "Single Transferable Vote". It's purpose 
!   to produce better clustering results.
!   
!   The membership of a ballot 'b' with weight w(b) and candidate rankings 
!   (c1,...,cn) in a candidate cluster Ci is initially determined by the 
!   rank 'i' of ci in 'b', using a decreasing fractional series in the rank.

!   Then the membership w(Ci) of Ci is increased: w(Ci)= w(Ci) + w(ci,b)
!   where w(ci,b) = wt(i)*w(b). That is, the membership 
!   is distributed over the candidate clusters of the ranked candidates 
!   instead of being focused only on the top ranked cluster C1, as in STV.

!   When a candidate ci is deleted because w(Ci) is the smallest, then
!   its membership w(ci,b) in a ballot 'b' is redistributed to the 'k'
!   subsequent unelected, or hopeful, candidate clusters by the same 
!   fractional series:
!   Let j = 1,...,k be the indices of these subsequent hopefuls in 
!   ballot 'b', yielding the transfer w(Cj)= w(Cj) + wt(j)*w(ci,b) 
!   where w(ci,b) = current membership of 'b' in Ci.

!   When a candidate ci is elected because w(Ci) is the largest and 
!   exceeds the Droop formula, then, to reduce w(Ci) to the Droop value,
!   a fraction 'fr' of its membership w(ci,b) in a ballot 'b' is 
!   redistributed to the subsequent 'k' hopeful candidate clusters 
!   by the same fractional series:
!   Let j = 1,...,k be the indices of these hopefuls in this ballot, 
!   yielding the transfer w(Cj)= w(Cj) + wt(j)*fr*w(ci,b) where
!   w(ci,b) = current membership of 'b' in Ci. 

!   That is, fr = (surplus weight) / (transferable weight)

!   where the Droop formula is 
!        Quota = (total weight of the hopeful candidate clusters) / 
!                (1 + # candidates remaining to be elected)

!   Note that when a candidate cluster is eliminated, there may be 
!   ballots which rank that candidate with a positive weight but lack 
!   any subsequent hopefuls. In this case the weight is transferred 
!   to the cluster of "independents". Similarly, when a candidate is 
!   elected a ballot which weights that candidate may lack subsequent 
!   hopefuls, so its membership fraction is transferred to the independents.

! Input:
    Integer, Intent(in) :: nc              ! # candidates
    Real,    Intent(in) :: pt_val(:)       ! (mt) Rating or ranking points corresponding to 
                                           !      rating or ranking levels.
    Real,    Intent(in) :: wtb(:)          ! (nb) Ballot weights, summing to 'np'
    Integer, Intent(in) :: ballot(0:,:)    ! (0:mr,nb). (0,:) = 'r' candidates in preferential order
    Integer, Intent(in) :: ballot2(0:,:)   ! (0:mr,nb). corresponding ranking or rating levels
                                           !   with (0,:) = 'n' positively rated candidates

! Output:
    Integer,        Intent(out) :: elected(:)   ! (nc) List of elected candidates in preferential order

    Type(Set_list), Intent(out) :: DTVmb(0:)    ! (0:nb) Ballot ranking structure for clustering (on output)
                                                ! Ballot sets (0): 
                                                !   %n        = nc
                                                !   %p        = np
                                                !   %svl      = elected vote = Sum(%val(set(:np)))
                                                !   %smb      = independents vote = Sum(%val(set(np+1:)))
                                                !   %set(nc)  = candidates in elected ordering, with
                                                !               the first 'np' elected
                                                !   %lev(nc)  = candidate clusters by decreasing cluster vote
                                                !   %val(nc)  = candidate cluster votes, not reordered

                                                ! Ballots (1:nb):
                                                !   %n     = # elected candidates
                                                !   %p     = # positively rated candidates
                                                !   %r     = # ranked candidates
                                                !   %set(r)= 'n' elected candidates, then unelected,
                                                !            in original ranking order
                                                !   %val(r)= corresponding membership weights 
                                                !   %svl   = elected membership weight= Sum(%val(:n))
                                                !   %smb   = total membership weight  = Sum(%val)
! Local:
    Real    :: vote0(nc)     ! Size of the initial candidate clusters
    Real    :: vote(nc)      ! Size of each candidate cluster at the 
                             !  current stage of the DTV cycle
    Real    :: vt_ord(nc)    ! (nh) Hopeful candidate clusters, by decreasing size
    Real    :: quota(nc)     ! Droop formula quotas by DTV cycle

    Real    :: rate(nc)      ! Mean rating vector of the hopefuls
    Real    :: rateV(nc)     ! Current 'rate' * 'vote'
    Real    :: rt_ord(nc)    ! (nh) Hopeful candidate clusters, by decreasing 'rateV'

    Real    :: w_trans(nc,3) ! Records data by STV cycle
                             ! (:,1) = weight transfered to electeds, (:,2) = total weight transfered
                             ! (:,3) = weight of fully retained ballots
    Integer :: n_trans(nc,3) ! Records data by STV cycle
                             ! (:,1) = # ballot transfers, (:,2) = # ballots retained
                             ! (:,3) = deleted (neg) or elected (pos) candidate

    Integer :: Elect(nc)     ! = i if ith elected, 0 if hopeful, -i if ith eliminated
    Integer :: Elect1(nc)    ! (nc) = 'Elect' except < 0 where 'Elect' > 0 but not all elected

    Real,    Allocatable :: wt(:,:)  ! (n,2)  (:,1) = increasing transferable ballots weights 
                                     !        (:,2) =
    Integer, Allocatable :: kw(:,:)  ! (n,2)  (:,1) = corresponding ballot indices
                                     !        (:,2) = reodrdering key for wt(:,1)
    Real    :: rt, ratio, trans, retain(3)
    Integer :: nsp, ntr, nrt(2)

    Real    :: pt(nc), vote1(nc)
    Integer :: lh(nc), lst(nc), key(nc)

    Real, Parameter :: min_trans= 0.001, eps= 0.0001
    Logical :: Add

    Real    :: surplus, trans_wt, transfer_fac, transferable, independents
    Real    :: vt, tr, swt, sm_pt, tot_hope, prior, total, vt_elc, vt_ind
    Integer :: mr, mt, n1, n2, nb, nd, nh, ne, np, nr, nc1, ndf, np1, top, last
    Integer :: b, i, j, k, l, m, n, p, r, c1, cn, h1, ih, il, it, iv

    np= Nint(Sum(wtb));  np1= np + 1;  nc1= nc + 1;  mr= Ubound(ballot,1)
    nb= Size(wtb);  mt= Size(pt_val);  ndf= nc - np

    Call Out ("WDTV2 algorithm for # electeds",np, "out of # candidates",nc,ln=1)

!   Initial ballot memberships in candidate clusters

    elected= 0;  vote0= 0;  Elect= 0;  independents= 0
    quota= 0;  n_trans= 0;  w_trans= 0
      
    Ballot_loop0 : Do b= 1,nb
      DTVmb(b)%smb= wtb(b);  r= ballot(0,b);  p= ballot2(0,b)
      DTVmb(b)%p= p
      
      If (Associated(DTVmb(b)%set)) DeAllocate(DTVmb(b)%set)
      If (Associated(DTVmb(b)%lev)) DeAllocate(DTVmb(b)%lev)
      If (Associated(DTVmb(b)%val)) DeAllocate(DTVmb(b)%val)
      Allocate(DTVmb(b)%set(r),DTVmb(b)%lev(r), DTVmb(b)%val(r))

      DTVmb(b)%set= ballot(1:r,b)  
      DTVmb(b)%lev= ballot2(1:r,b)

      pt(:p)= pt_val(DTVmb(b)%lev(:p)) 
      sm_pt= Sum(pt(:p));  pt(:p)= pt(:p) / sm_pt

      DTVmb(b)%val(:p)= pt(:p) * wtb(b);  DTVmb(b)%val(p+1:)= 0
      vote0(DTVmb(b)%set(:p))= vote0(DTVmb(b)%set(:p)) + DTVmb(b)%val(:p)

    End do Ballot_loop0

    If (pr_out >= 1) then
      total= Sum(vote0)
      Call Out ("Initial DTV vote total", total, ln=1)
      Call Out ("Corresponding candidate cluster votes, unordered", vote0)
    End if
      
    vote= vote0
    key= 'ID'  ! Hopeful candidates, to be sorted by weight
    ne= 0      ! # currently elected
    nd= 0      ! # currenlty deleted
      
    DTV_cycle : Do iv= 1,nc
        
      Call List_of_true (Elect == 0, nh,key)  ! nh = # hopeful candidates remaining
      If (nh == 0) Exit DTV_cycle  ! Done

!     Compute the mean ratings * vote over the positive hopeful candidates

      rate= 0
      Rate_loop : Do b= 1,nb
        p= DTVmb(b)%p;  Call List_of_true (Elect(DTVmb(b)%set(:p)) == 0, m,lh)

        If (m > 0) then
          lst(:m)= DTVmb(b)%set(lh(:m))
          rate(lst(:m))= rate(lst(:m)) + pt_val(DTVmb(b)%lev(lh(:m))) * wtb(b)
        End if
      End do Rate_loop 

      rate= rate / np;  rateV= rate * vote

      If (nh > 1) then  ! Sort the hopeful candidates from largest to smallest
                        ! based on the 'rateV' scaling
        rt_ord(:nh)= rateV(key(:nh))
        Call Sort (.false., rt_ord(:nh), lh(:nh))
        key(:nh)= key(lh(:nh))

        vt_ord(:nh)= vote(key(:nh))
        tot_hope= Sum(vt_ord(:nh))  ! Total vote for the hopefuls
      Else
        cn= key(1);  rt_ord(1)= rateV(cn)
        vt_ord(1)= vote(cn);  tot_hope= vote(cn)
      End if
      
      last= key(nh);  quota(iv)= tot_hope / (np1 - ne)  ! Droop formula for the quota

      pt(:nh)= vt_ord(:nh) - quota(iv)

      j= First_true(pt(:nh) >= 0);  Add= j > 0

      If (nd >= ndf .and. .not.Add) then
        j= Maxloc(pt(:nh), 1);  Add= .true.
      End if

      If (Add) then
        top= key(j);  prior= vote(top);  surplus= prior - quota(iv)
      End if

      If (pr_out >= 1) then
        Call Out ("Hopeful candidates sorted by decreasing rate-vote value", key(:nh))
        Call Out ("vs vote size", vt_ord(:nh))
        Call Out ("With total hopeful weight",tot_hope, "and quota",quota(iv))
      End if

      If (.not.Add) then ! Eliminate the bottom rate-vote candidate 'last' = key(nh)

        Call Out ("Next delete the bottom rate-vote candidate", last, ln=1)
        m= 0;  l= 0;  rt= 0;  swt= 0;  vt_elc= 0

        Elect1= Elect;  If (ne < np) Where (Elect > 0) Elect1= -10  ! Distribute only to next hopeful
        
        Elim_loop : Do b= 1,nb

          il= First_true(DTVmb(b)%set == last);  If (il < 1) Cycle Elim_loop  ! 'last' not ranked
          trans_wt= DTVmb(b)%val(il);  If (trans_wt <= 0) Cycle Elim_loop     ! No weight to transfer

          p= DTVmb(b)%p
          Call List_of_true (Elect1(DTVmb(b)%set(:p)) >= 0 .and. DTVmb(b)%set(:p) /= last, ih,lh)

          If (ih > 0) then  ! Transfer weight & distribute to positively rated candidates

            DTVmb(b)%val(il)= 0;  vote(last)= vote(last) - trans_wt

            pt(:ih)= pt_val(DTVmb(b)%lev(lh(:ih))) 
            sm_pt= Sum(pt(:ih));  pt(:ih)= pt(:ih) / sm_pt

            Do i= 1,ih
              k= lh(i);  tr= pt(i) * trans_wt
              DTVmb(b)%val(k)= DTVmb(b)%val(k) + tr
              cn= DTVmb(b)%set(k);  vote(cn)= vote(cn) + tr
              If (Elect(cn) > 0) vt_elc= vt_elc + tr
            End do

            l= l + 1;  swt= swt + trans_wt
          Else      ! No distribution: fully retain
            m= m + 1;  rt= rt + trans_wt
          End if
        End do Elim_loop

        nd= nd + 1;  Elect(last)= -nd;  elected(nc1 - nd)= last
        vt_ind= vote(last);  total= Sum(vote)

        independents= independents + vt_ind

        n_trans(iv,1)= l;  n_trans(iv,2)= m;  n_trans(iv,3)= -last
        w_trans(iv,1)= vt_elc;  w_trans(iv,2)= swt;  w_trans(iv,3)= rt 
      
        If (pr_out > 1) then
          Call Out ("# delete ballot transfers, # retained ballots, elected candidate",n_trans(iv,:))
          Call Out ("Weight transfered to electeds, total transfered, total fully retained",w_trans(iv,:))
          Call Out ("Remaining weight to independents",vt_ind, &
                    "current independents weight",independents)
          Call Out ("Current vote total",total)
          Call Out ("Current cluster sizes",vote)
        End if

      Else  ! elect the top candidate 'top'
        Call Out ("Next add the top rate-vote candidate", top, ln=1)

!       Identify ballots with partially transferable vs fully retained 
!       ballot weight for candidate 'top'

        ne= ne + 1;  m= 0;  l= 0;  rt= 0;  transferable= 0
        Elect1= Elect;  If (ne < np) Where(Elect > 0) Elect1= -10  ! Distribute only to next hopeful
          
        Add_loop0 : Do b= 1,nb
          DTVmb(b)%n= 0;     it= First_true(DTVmb(b)%set == top)
            DTVmb(b)%r= it;  If (it < 1) Cycle Add_loop0  ! 'top' not ranked
          trans_wt= DTVmb(b)%val(it);  If (trans_wt <= 0) Cycle Add_loop0     ! No weight to transfer

          p= DTVmb(b)%p
          j= First_true(Elect1(DTVmb(b)%set(:p)) >= 0 .and. DTVmb(b)%set(:p) /= top)
          DTVmb(b)%n= j

          If (j > 0) then  ! l = # ballots partially transferable to a positively rated candidate
            l= l + 1;  transferable= transferable + DTVmb(b)%val(it)
          Else             ! m = # ballots fully retained
            m= m + 1;  rt= rt + DTVmb(b)%val(it)
          End if
        End do Add_loop0

        ntr= l;  nrt= 0;  retain= 0;  nrt(1)= m;  retain(1)= rt
        nsp= 0;  vt_elc= 0

        If (ntr > 0) then
          Call RetainWD2 (top,mr, ntr,quota(iv),pt_val, Elect,Elect1, DTVmb,vote, &
                          nrt,retain, nsp,trans,vt_elc)
        
          Elect(top)= ne;  elected(ne)= top;  m= nrt(1) + nrt(2);  rt= Sum(retain)
          n_trans(iv,1)= nsp;       n_trans(iv,2)= m;    n_trans(iv,3)= top
          w_trans(iv,1)= vt_elc;  w_trans(iv,2)= trans;  w_trans(iv,3)= rt

          vt= prior - vote(top);  total= Sum(vote);  ratio= nsp / (nrt(2)+eps)  

          If (pr_out >= 1) then
            Call Out ("Quota",quota(iv), "vs actual retained weight",rt, ln=1)
            Call Out ("Surplus",surplus, "vs actually transfered",trans)
            Call Out ("vs actual top candidate decrease",vt, "total vote",total)
            Call Out ("Ratio of # transfers to # retentions in Warren",ratio)

            Call Out ("# transfers to surplus, # fully retained ballots, elected candidate",n_trans(iv,:))
            Call Out ("Weight transfered to electeds, total transfered, total fully retained",w_trans(iv,:))
            Call Out ("for current cluster sizes",vote)
          End if
        Else  ! Elect without transfering
          Elect(top)= ne;  elected(ne)= top
          n_trans(iv,2)= m;  n_trans(iv,3)= top;  w_trans(iv,3)= rt

          Call Out ("No transfers: WDTV elected candidate",top, &
                    "retains vote",vote(top), ln=1)
        End if
      End if
    End do DTV_cycle

    If (pr_out >= 0.5) then
      Call Out (1,"WDTV2 # ballot transfers, # retained ballots, elected candidate",n_trans)
      Call Out (1,"Weight transfered to electeds, total transfered, total fully retained",w_trans)
      Call Out ("Final cluster sizes",vote)

      vt_elc= Sum(vote(elected(:np)));  vt= Sum(vote) 
      Call Out ("Final elected total",vt_elc, "vote total",vt)
    End if

!   Final DTVmb output: Redo ballots so that the electeds 
!   come first in each ballo. Otherwise maintaining the original
!   ranking order.

    Call Final_STVmb (mr,nc, elected,Elect, DTVmb, vote1)

  End Subroutine WDTV2
   

  Subroutine Final_STVmb (mr,nc, elected,Elect, STVmb, vote)

!   Output STVmb so that each ballot ranks elected candidates first
!   in the original ranking order, then the independents in the 
!   original ranking order

    Integer,           Intent(in) :: mr         ! max # ranked candidates
    Integer,           Intent(in) :: nc         ! # candidates
    Integer,           Intent(in) :: elected(:) ! (nc) Candidates in final preferential order (see vote)
    Integer,           Intent(in) :: Elect(:)   ! (nc) Elected status (+) or deleted status (-) 
                                                !      of each candidate cluster. No 0 values.

    Type(Set_list), Intent(inout) :: STVmb(0:)  ! (0:nb) Ballot ranking structure for clustering
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
                                                !   %p     = # positively rated candidates
                                                !   %r     = # ranked or rated candidates
                                                !   %set(r)= 'n' elected candidates, then unelected,
                                                !            in original ranking order
                                                !   %val(r)= corresponding membership weights 
                                                !   %svl   = elected membership weight= Sum(%val(:n))
                                                !   %smb   = total membership weight  = Sum(%val)
    Real,             Intent(out) :: vote(:)    ! (nc) Final candidate cluster weights in decreasing order

! Local:
    Real    :: tmp1(mr), tmp2(mr), vot(nc)
    Integer :: key1(mr), key2(mr), elc(nc)
    Real    :: s1, s2, vt, vt_elc, vt_ind, total
    Integer :: b, k, m, n, p, r, cn, nb, np

    nb= Ubound(STVmb,1);  np= Count(Elect > 0);  vote= 0
    elc= 0;  tmp1= 0;  tmp2= 0;  key1= 0;  key2= 0

!   Reorder the ranked candidates of each ballot
!   to place the elected candidates first.
    
    Do b= 1,nb
      r= Size(STVmb(b)%set);  STVmb(b)%r= r
      n= 0;  m= 0;  STVmb(b)%svl= 0
        
      Do k= 1,r
        cn= STVmb(b)%set(k);  vt= STVmb(b)%val(k)
          
        If (Elect(cn) > 0) then
          n= n + 1;  key1(n)= cn;  tmp1(n)= vt  ! Elected vote
        Else
          m= m + 1;  key2(m)= cn;  tmp2(m)= vt  ! Independent vote
        End if
      End do

      STVmb(b)%n= n;  s1= 0;  s2= 0

      If (n > 0) then
        STVmb(b)%set(:n)= key1(:n);  STVmb(b)%val(:n)= tmp1(:n);  s1= Sum(tmp1(:n))
      End if

      If (m > 0) then
        STVmb(b)%set(n+1:)= key2(:m);  STVmb(b)%val(n+1:)= tmp2(:m);  s2= Sum(tmp2(:m))
      End if

      STVmb(b)%svl= s1;  STVmb(b)%smb= s1 + s2
    
      vote(STVmb(b)%set)= vote(STVmb(b)%set) + STVmb(b)%val
    End do

!   STVmb(0) output:

    vt_elc= Sum(vote(elected(:np)))   
    vt_ind= Sum(vote(elected(np+1:)))

    total= vt_elc + vt_ind
    vot= vote;  Call Sort (.false.,vot,elc)

    If (Associated(STVmb(0)%set)) DeAllocate(STVmb(0)%set, STVmb(0)%lev, STVmb(0)%val)
    Allocate(STVmb(0)%set(nc), STVmb(0)%lev(nc), STVmb(0)%val(nc))

    STVmb(0)%n  = nc;       STVmb(0)%p  = np
    STVmb(0)%svl= vt_elc;   STVmb(0)%smb= vt_ind 
    STVmb(0)%set= elected;  STVmb(0)%lev= elc;  STVmb(0)%val= vote

    If (pr_out >= 1) then
      Call Out ("For # elected",np, "with elected vote",vt_elc, ln=1)
      Call Out ("and independents vote",vt_ind, "with total vote",total)
      Call Out ("Candidates in elected order",elected)
      Call Out ("Candidate clusters ordered by decreasing vote",elc)
      Call Out ("Candidate cluster votes, not reordered", vote)
    End if

    vt= Sum(vote(elc(:np)))
    If (Abs(vt - vt_elc) > 0.0001) then
      Call Out ("Error in 'Final_STVmb': cluster elected vote total from elected order",vt_elc, &
                 "differs from the elected vote total from size order",vt)
    End if
  End Subroutine Final_STVmb
 
  Subroutine RetainW (top, ntr,quota, Elect, STVmb,vote, nrt,retain, nsp,trans,vt_elc)

!   Warren method of STV.

    Integer,           Intent(in) :: top        ! top candidate cluster: 'quota' to be retained, rest transfered
    Integer,           Intent(in) :: ntr        ! # transferable ballots
    Real,              Intent(in) :: quota      ! Min election weight
    Integer,           Intent(in) :: Elect(:)   ! (nc) Elected status (+), deleted status (-), or hopeful (0)

    Type(Set_list), Intent(inout) :: STVmb(0:)  ! (0:nb) 
    Real,           Intent(inout) :: vote(:)    ! (nc) Candidate cluster sizes

    Integer,        Intent(inout) :: nrt(:)     ! (2) Retained ballots:  1 = # non-transferable, 
                                                !   2 = # transferable but fully retained
    Real,           Intent(inout) :: retain(:)  ! (3) Retained ballot weights
                                                !   1 = non-transferable, 2 = transferable but fully retained
                                                !   3 = nsp * rtw = retained at 'rtw' per ballot

    Integer,         Intent (out) :: nsp        ! # ballots with partially transfered weight
    Real,            Intent (out) :: trans      ! Total transfered weight
    Real,            Intent (out) :: vt_elc     ! Weight transfered to elected candidates
! Local:
    Integer :: kw(ntr,2)
    Real    :: wt(ntr,2)

    Logical :: Update
    Real    :: sm, qrt, rtw, sm_rt, trans_wt
    Integer :: b, j, l, m, cn, it, l1

    l= 0
    Add_loop1 : Do b= 1,Ubound(STVmb,1)
      If (STVmb(b)%n <= 0) Cycle Add_loop1
      l= l + 1;  kw(l,1)= b;  it= STVmb(b)%r
      wt(l,1)= STVmb(b)%val(it)
    End do Add_loop1

    Call Sort (.true.,wt(:,1), kw(:,2));  kw(:,1)= kw(kw(:,2),1)

    qrt= quota - retain(1)

    If (qrt <= 0) then ! Retained weight may exceed the quota
      rtw= 0;  l1= ntr+1;   wt(:,2)= wt(:,1);  sm= Sum(wt(:,1))

    Else
      sm= 0;  wt(:,2)= 0;  l1= 1;  Update= .true.

      Do l= 1,ntr
        If (Update) then
          rtw= (qrt - sm) / (ntr+1 - l)
          rtw= Max(rtw,0.0)  ! If 0, then retained weight may exceed the quota

          If (wt(l,1) > rtw) then
            wt(l,2)= wt(l,1) - rtw;  Update= .false.
          Else
            l1= l+1;  sm= sm + wt(l,1)
          End if
        Else
          wt(l,2)= wt(l,1) - rtw
        End if
      End do
    End if

    m= l1 - 1;  nrt(2)= m;  nsp= ntr - m
    retain(2)= sm;  retain(3)= nsp * rtw
    trans= Sum(wt(l1:,2));  sm_rt= Sum(retain)

    vt_elc= 0 

    Do l= l1,ntr
      b= kw(l,1);  trans_wt= wt(l,2)

!     Subtract from the 'top' candidate cluster

      it= STVmb(b)%r;  STVmb(b)%val(it)= STVmb(b)%val(it) - trans_wt
      vote(top)= vote(top) - trans_wt

!     Distribute to the candidate cluster specified by 'j'

      j= STVmb(b)%n;  STVmb(b)%val(j)= STVmb(b)%val(j) + trans_wt
      cn= STVmb(b)%set(j);  vote(cn)= vote(cn) + trans_wt
      If (Elect(cn) > 0) vt_elc= vt_elc + trans_wt
    End do

    If (pr_out >= 1) then
      Call Out ("# transferable ballots",ntr, "Total weight transfered",trans, ln=1) 
      Call Out ("reduced quota",qrt, "& retained weight per ballot",rtw)

      Call Out ("# ballots without transfer",nrt(2), "their total weight",retain(2))
      Call Out ("total partially retained weight",retain(3))
      Call Out (-1,"Transferable ballot weights, actual transfers",wt)
    End if
  End Subroutine RetainW


  Subroutine RetainWD (top, ntr,quota,dis, Elect,Elect1, DTVmb,vote, nrt,retain, nsp,trans,vt_elc)

!   Warren method of DTV.

    Integer,           Intent(in) :: top        ! top candidate cluster: 'quota' to be retained, rest transfered
    Integer,           Intent(in) :: ntr        ! # transferable ballots
    Real,              Intent(in) :: quota      ! Min election weight
    Real,              Intent(in) :: dis(:,:)   ! (mr,mr) Point distribution for transfers
    Integer,           Intent(in) :: Elect(:)   ! (nc) Elected status (+), deleted status (-), or hopeful (0)
    Integer,           Intent(in) :: Elect1(:)  ! (nc) = 'Elect' except may include prior electeds

    Type(Set_list), Intent(inout) :: DTVmb(0:)  ! (0:nb) 
    Real,           Intent(inout) :: vote(:)    ! (nc) Candidate cluster sizes

    Integer,        Intent(inout) :: nrt(:)     ! (2) Retained ballots:  1 = # non-transferable, 
                                                !   2 = # transferable but fully retained
    Real,           Intent(inout) :: retain(:)  ! (3) Retained ballot weights
                                                !   1 = non-transferable, 2 = transferable but fully retained
                                                !   3 = nsp * rtw = retained at 'rtw' per ballot

    Integer,         Intent (out) :: nsp        ! # ballots with partially transfered weight
    Real,            Intent (out) :: trans      ! Total transfered weight
    Real,            Intent (out) :: vt_elc     ! Weight transfered to elected candidates
! Local:
    Integer :: kw(ntr,2), lh(Size(dis,1))
    Real    :: wt(ntr,2)

    Logical :: Update
    Real    :: sm, gtr, qrt, rtw, sm_rt, trans_wt
    Integer :: b, i, j, k, l, m, cn, ih, it, l1

    l= 0;  kw= 0;  wt= 0
    Add_loop1 : Do b= 1,Ubound(DTVmb,1)
      If (DTVmb(b)%n <= 0) Cycle Add_loop1
      l= l + 1;  kw(l,1)= b;  it= DTVmb(b)%r
      wt(l,1)= DTVmb(b)%val(it)
    End do Add_loop1

    Call Sort (.true.,wt(:,1), kw(:,2));  kw(:,1)= kw(kw(:,2),1)

    qrt= quota - retain(1)

    If (qrt <= 0) then ! Retained weight may exceed the quota
      rtw= 0;  l1= ntr+1;   wt(:,2)= wt(:,1);  sm= Sum(wt(:,1))

    Else
      sm= 0;  wt(:,2)= 0;  l1= 1;  Update= .true.

      Do l= 1,ntr
        If (Update) then
          rtw= (qrt - sm) / (ntr+1 - l)
          rtw= Max(rtw,0.0)  ! If 0, then retained weight may exceed the quota

          If (wt(l,1) > rtw) then
            wt(l,2)= wt(l,1) - rtw;  Update= .false.
          Else
            l1= l+1;  sm= sm + wt(l,1)
          End if
        Else
          wt(l,2)= wt(l,1) - rtw
        End if
      End do
    End if

    m= l1 - 1;  nrt(2)= m;  nsp= ntr - m
    retain(2)= sm;  retain(3)= nsp * rtw
    trans= Sum(wt(l1:,2));  sm_rt= Sum(retain)

    vt_elc= 0 

    Do l= l1,ntr
      b= kw(l,1);  trans_wt= wt(l,2)
      it= DTVmb(b)%r;  j= DTVmb(b)%n

!     Subtract from the 'top' candidate cluster

      DTVmb(b)%val(it)= DTVmb(b)%val(it) - trans_wt
      vote(top)= vote(top) - trans_wt

      Call List_of_true (Elect1(DTVmb(b)%set(j:)) >= 0 .and. DTVmb(b)%set(j:) /= top, ih,lh)
      lh(:ih)= (j-1) + lh(:ih)

!     Distrubte the 'trans_wt' weight to the 'lh' list

      Do i= 1,ih
        k= lh(i);  gtr= dis(i,ih) * trans_wt
        DTVmb(b)%val(k)= DTVmb(b)%val(k) + gtr
        cn= DTVmb(b)%set(k);  vote(cn)= vote(cn) + gtr
        If (Elect(cn) > 0) vt_elc= vt_elc + gtr
      End do
    End do

    If (pr_out >= 1) then
      Call Out ("# transferable ballots",ntr, "Total weight transfered",trans, ln=1) 
      Call Out ("reduced quota",qrt, "& retained weight per ballot",rtw)

      Call Out ("# ballots without transfer",nrt(2), "their total weight",retain(2))
      Call Out ("total partially retained weight",retain(3))
      Call Out (-1,"Transferable ballot weights, actual transfers",wt)
    End if
  End Subroutine RetainWD


  Subroutine RetainWD2 (top,mr, ntr,quota,pt_val, Elect,Elect1, DTVmb,vote, &
                        nrt,retain, nsp,trans,vt_elc)

!   Warren method of DTV2.

    Integer,           Intent(in) :: top        ! top candidate cluster
    Integer,           Intent(in) :: mr         ! max # rated ballots
    Integer,           Intent(in) :: ntr        ! # transferable ballots
    Real,              Intent(in) :: quota      ! Weight to be retained by the 'top' candidate
    Real,              Intent(in) :: pt_val(:)  ! (mt) Point values for rating / ranking levels

    Integer,           Intent(in) :: Elect(:)   ! (nc) Elected status (+), deleted status (-), or hopeful (0)
    Integer,           Intent(in) :: Elect1(:)  ! (nc) = 'Elect' except may include prior electeds

    Type(Set_list), Intent(inout) :: DTVmb(0:)  ! (0:nb) 
    Real,           Intent(inout) :: vote(:)    ! (nc) Candidate cluster sizes

    Integer,        Intent(inout) :: nrt(:)     ! (2) Retained ballots:  1 = # non-transferable, 
                                                !   2 = # transferable but fully retained
    Real,           Intent(inout) :: retain(:)  ! (3) Retained ballot weights
                                                !   1 = non-transferable, 2 = transferable but fully retained
                                                !   3 = nsp * rtw = retained at 'rtw' per ballot

    Integer,         Intent (out) :: nsp        ! # ballots with partially transfered weight
    Real,            Intent (out) :: trans      ! Total transfered weight
    Real,            Intent (out) :: vt_elc     ! Weight transfered to elected candidates
! Local:
    Integer :: kw(ntr,2), lh(mr)
    Real    :: wt(ntr,2), pt(mr)

    Logical :: Update
    Real    :: sm, tr, qrt, rtw, sm_rt, sm_pt, trans_wt
    Integer :: b, i, j, k, l, m, p, cn, ih, it, l1

    l= 0;  kw= 0;  wt= 0
    Add_loop1 : Do b= 1,Ubound(DTVmb,1)
      If (DTVmb(b)%n <= 0) Cycle Add_loop1
      l= l + 1;  kw(l,1)= b;  it= DTVmb(b)%r
      wt(l,1)= DTVmb(b)%val(it)
    End do Add_loop1

    Call Sort (.true.,wt(:,1), kw(:,2));  kw(:,1)= kw(kw(:,2),1)

    qrt= quota - retain(1)

    If (qrt <= 0) then ! Retained weight may exceed the quota
      rtw= 0;  l1= ntr+1;   wt(:,2)= wt(:,1);  sm= Sum(wt(:,1))

    Else
      sm= 0;  wt(:,2)= 0;  l1= 1;  Update= .true.

      Do l= 1,ntr
        If (Update) then
          rtw= (qrt - sm) / (ntr+1 - l)
          rtw= Max(rtw,0.0)  ! If 0, then retained weight may exceed the quota

          If (wt(l,1) > rtw) then
            wt(l,2)= wt(l,1) - rtw;  Update= .false.
          Else
            l1= l+1;  sm= sm + wt(l,1)
          End if
        Else
          wt(l,2)= wt(l,1) - rtw
        End if
      End do
    End if

    m= l1 - 1;  nrt(2)= m;  nsp= ntr - m
    retain(2)= sm;  retain(3)= nsp * rtw
    trans= Sum(wt(l1:,2));  sm_rt= Sum(retain)

    vt_elc= 0 

    Do l= l1,ntr
      b= kw(l,1);  trans_wt= wt(l,2)
      it= DTVmb(b)%r;  j= DTVmb(b)%n;  p= DTVmb(b)%p

!     Subtract from the 'top' candidate cluster

      DTVmb(b)%val(it)= DTVmb(b)%val(it) - trans_wt
      vote(top)= vote(top) - trans_wt

      Call List_of_true (Elect1(DTVmb(b)%set(j:p)) >= 0 .and. DTVmb(b)%set(j:p) /= top, ih,lh)
      lh(:ih)= (j-1) + lh(:ih)

      pt(:ih)= pt_val(DTVmb(b)%lev(lh(:ih))) 
      sm_pt= Sum(pt(:ih));  pt(:ih)= pt(:ih) / sm_pt

      Do i= 1,ih
        k= lh(i);  tr= pt(i) * trans_wt
        DTVmb(b)%val(k)= DTVmb(b)%val(k) + tr
        cn= DTVmb(b)%set(k);  vote(cn)= vote(cn) + tr
        If (Elect(cn) > 0) vt_elc= vt_elc + tr
      End do
    End do

    If (pr_out >= 1) then
      Call Out ("# transferable ballots",ntr, "Total weight transfered",trans, ln=1) 
      Call Out ("reduced quota",qrt, "& retained weight per ballot",rtw)

      Call Out ("# ballots without transfer",nrt(2), "their total weight",retain(2))
      Call Out ("total partially retained weight",retain(3))
      Call Out (-1,"Transferable ballot weights, actual transfers",wt)
    End if
  End Subroutine RetainWD2

End Module Misc_methods
    