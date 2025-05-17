
    
!   This module contains "Read_bal" to pre-process the district election ballot files.
!   This preprocessing includes re-ordering the candidates according to a restricted Borda Count.
!   It also includes the elimination of marginal candidates, where this fraction is so small 
!   that the candidate is not viable.
    
Module Clusters_pre

   Use Util
   Use Factorials
   Use Output
   Use Types
   Use Precisn
   Implicit None

   Real,    Parameter :: UVP_fac(5,5)= (/ (/    1.0,1.0,1.0,1.0,1.0 /), &
                                          (/0.67,   1.0,1.0,1.0,1.0 /), &
                                          (/0.50,0.80,  1.0,1.0,1.0 /), &
                                          (/0.40,0.65,0.85, 1.0,1.0 /), &
                                          (/0.33,0.55,0.75,0.90,1.0 /) /)       
                                            ! Undervote penalty factors
   Real, Parameter :: rand_sig= 1.0      ! Standard deviation for random number variations
   Real, Parameter :: Neg_frac1= 0.333   ! Fraction of 'mr' ratings to be used for negative ratings
                                         !   if Rating = 1 and if not otherwise specified
   
   Real, Parameter :: Max_pt= 10.0     ! Maximum point value for rating or ranking points
   Real, Parameter :: Marg_cut0= 0.10  ! Remove marginal candidates for ranking: 
                                       !    If Borda(c) < Marg_cut0 * Borda(np)
   Real, Parameter :: mem_fac= 1.5     ! Memory factor by which to increase the combinatorial limit
                                       !  for ratings due to inexactness of the point adjustment to neutrality
   Real, Parameter :: rank_parm= 0.50  ! Parameter for type ranking points
                                       !  > 1: Scale factor for a geometric decrease in points, 
                                       !        starting from Max_pt point for the top ranked candidate
                                       !  <= 1 & >= 0: Before rescaling to Max_pt, this is the point value
                                       !        for an arithmetic increase in the point increment above the 
                                       !        standard Borda Count. That is, x = 0 represents standard Borda
                                       !        with point values of (1,2,3,...,i) from the bottom ranked up, 
                                       !        while x > 0 represents point values of (1,2+x,3+3x,4+6x,...,i+(i-1)(i-2)x/2)
                                       !        corresponding to increments of (1,1+x,1+2x,1+3x,...,1+(i-1)x).
   Integer, Parameter :: Lim_nc= 10      ! Max reduced # candidate
   Integer, Parameter :: Lim_nb= 100000  ! Max # ballots

 Contains
    
   Subroutine Read_bal (Rating,mr_spec,mt_spec,UVP, id,District,Region,Year)
   
!    Read the ballots, depending on the format. Then pre-process and store in a standard format.

     Integer,       Intent(in) :: Rating     ! Specified rating type for the output disrict files: 
                                             !   -1 = strong ranking, 0 = weak ranking, 1 = rating. 
     Integer,       Intent(in) :: mr_spec    ! Max # ranked or rated candidates to be used for output
     Integer,       Intent(in) :: mt_spec    ! Max # rating levels to be used for output
     Integer,       Intent(in) :: UVP        ! "Under Vote Penalty", or scale factor < 1, on ballot weights.
                                             !   UVP = 1 means no penalty. In general the scale factor is nr/UVP
                                             !   when nr = # ranked or rated < UVP, and no penalty for nr >= UVP.
     Integer,       Intent(in) :: id         ! District #
     Character(*),  Intent(in) :: District   ! The electoral district
     Character(*),  Intent(in) :: Region     ! Region or jurisdiction containing the district
     Integer,       Intent(in) :: Year       ! Year of the election
!  Local:
    Integer, Allocatable :: wtb(:)      ! (nb) Raw ballot weights (# duplicates)
    Integer, Allocatable :: Bal(:,:,:)  ! (0:mr,2,nb)  Records all ballot data
                                        !  (0,1,b)  = 'nr' = # candidates ranked or rated
                                        !  (0,2,b)  = 'lp' = # pos candidates ranked or rated
                                        !  (i,1,b) = ith ranked or rated candidate: 1..nr
                                        !  (i,2,b) = corresponding ranking or rating level (nondecreasing)

     Integer :: bpm(0:9)                     ! Voting file input: 0: Rating0, 1: nbf, 2: nc0, 3: mr0, 4: mt0, 5: np, 
                                             !                    6: bal_fmt, 7: num, 8: dup, 9: cnt
                                             ! Rating0= 1 for rating levels
                                             !          0 for weak ranking  : ranking levels with equal rankings permitted
                                             !         -1 for strong ranking: candidates in ranking order
                                             ! nbf = # ballots recorded in the district file to be read, or determine this
                                             !       by reading to the end of the file, if nbf < 1.
                                             ! nc0 = original # candidates, before eliminating marginal candidates to get 'nc'
                                             ! mr0 = max # candidate rankings or ratings to be used (pos + neg for Rating0 = 1)
                                             ! mt0 = # rating levels with point conversion values 'pt_val' to be 
                                             !       read in if Rating0 >= 0 and mt0 >= 4. If mt0 < 4, 
                                             !       create default values 'rank_pt' of length 'mr' as needed.
                                             ! np  = # candidates to be elected
                                             ! bal_fmt = 1 for a single line of 'nr' candidates in ranking order,   
                                             !     beginning with the duplicate count, possibly terminated by a 0.
                                             !   = 2 for a single line, starting with the duplicate count, then 
                                             !     the ranking or rating level for each of the nc0 candidates
                                             !     in original order with 0 for unranked or rated candidates.
                                             !   = 3 same as bal_fmt = 1 except that the count of the # of ranked or
                                             !     rated candidates is inserted after the # duplicates.
                                             ! num = 1 for ballot numbering in each vote line (1st entry). Ignored.
                                             ! dup= first(1). 1 signifies # duplicate ballots (2nd entry)
                                             ! cnt= first(2). 1 signifies # candidates ranked or rated (3rd entry)
     Integer :: bal_fmt                      ! bpm(6)
     Integer :: Rating0                      ! bpm(0)

     Integer, Allocatable :: orig(:)         ! (nc0=>nc) orig(i)= the original candidate specified by the 
                                             !       reordered candidate i, according to the first Borda count (decreasing)
     Integer, Allocatable :: inv_orig(:)     ! (nc0) Inverse of 'orig'

     Real,        Pointer :: pt_val0(:)      ! (mt0) Decreasing point values for input rating or ranking levels
     Real,    Allocatable :: Borda0(:)       ! (nc0) Positive + negative Borda tally, if Rating = 1, else ranking tally
     Real,    Allocatable :: BordaP0(:)      ! (nc0) Borda tally of the positive ratings of each voter, if Rating = 1
     Real,    Allocatable :: BordaN0(:)      ! (nc0) Borda tally of the negative ratings of each voter, if Rating = 1

     Real,        Pointer :: pt_val(:)       ! (mt) Decreasing point values for output rating or ranking levels
     Real,    Allocatable :: Borda(:)        ! (nc) Positive + negative Borda tally, if Rating = 1, 
                                             !      else ranking tally, in decreasing order
     Integer, Allocatable :: ordBd(:)        ! (nc) Reordering of 'Borda' by 'Ballot_out'
     Integer, Allocatable :: key(:), cand(:), rate(:)
     Real,    Allocatable :: tmp(:)
     Character(80) :: label

     Integer :: nbf,nb, nc0,nc, np, mr0,mr, mrp,mrn
     Integer :: mt0, mt, mtp
     Real    :: dx, fac, inc, tot_wt

     Logical :: ReOrd
     Integer :: ios, tot_cnt, Alloc_err, first(2)
     Integer :: m0, m1, nBP, nBN, pl, nl, mb, np2, n10
     Integer :: b, i, j, k, l, m, n, cn, ib, n0
     
     Call Out ("Enter 'Read_bal':")
     Call Out ("Read and consolidate ballots for voting district "//District)
     Write (8, '(3A, I6)') "  region= ",Region, "  and year= ",Year
     
!    First read the ballot parameters

     label= "Elections0\"//District

     Open(9, File= Trim(Label), IOstat=ios, Status='Old', Action='Read')
         
     Read(9,*) bpm 
       
     Rating0= bpm(0);  nbf= bpm(1);  nc0= bpm(2);  mr0= bpm(3);  mt0= bpm(4)
     np= bpm(5);  bal_fmt= bpm(6);  first(1)= bpm(8);  first(2)= bpm(9)

     mr0= Min(mr0,mr_spec)
     np2= Min(np+1 + Nint(np/3.0),nc0)  ! Target min # candidates 'nc'
     
     Call Out ("From the input ballot file of Rating type", Rating0)
     Call Out ("# candidates to be elected", np, "out of # candidates tallied",nc0)
     Call Out ("Max # candidate rankings or ratings to be used from each ballot",mr0)
     Call Out ("# rankings or rating levels for specified point conversions",mt0)
     Call Out ("Target min # candidates to be considered for election",np2)
     
     If (np < 2 .or. nc0 < np2) then
       Call Out ("Voting file error for # to be elected",np, "vs # candidates",nc0)
       Stop
     End if

     If (Rating < Rating0) then
       Call Out ("Error: Target rating",Rating, " is less than file input rating", Rating0)
       Stop
     End if

!    Read in and adjust, or compute, the point conversion values 'pt_val0' 
!    in the case of ratings or weak rankings (Rating < Rating0)

     If (Rating0 >= 0) then
       If (mt0 < 4) mt0= mr0  ! Compute pt_val0
       Call Point_values0 (Rating0,Max_pt, mr0,mt0, pt_val0)  
     Else
       Call Out ("The file type specifies strong candidate ranking")
       mr0= Max(Min(mr0,nc0), 4);  
       mt0= 1;  Allocate(pt_val0(1))
       Read(9,*) pt_val0
     End if

!    First reading of the ballots. Use this to determine marginal candidates 
!    to be eliminated, also the candidate preferential reordering.
     
     Allocate(Borda0(nc0), BordaP0(nc0), BordaN0(nc0), orig(nc0), &
              inv_orig(nc0))

     n0= nbf
     Call Read0 (bpm(6:),Rating0,UVP, nc0,mr0,mt0, n0,pt_val0, &
                 nbf,tot_cnt, Borda0,BordaP0,BordaN0)

     Call Marginal_cand (Rating,Rating0,nc0, Borda0,BordaP0,BordaN0, &
                         nc,orig,inv_orig)

     If (Rating == 1) then
       Call Out ("The output data is for candidate rating")
     Else if (Rating == 0) then
       Call Out ("The output data is for weak candidate ranking")
     Else if (Rating == -1) then
       Call Out ("The output data is for strong candidate ranking")
     End if

!    Prepare to generate 'Rating' type data
     
     mr0= Min(mr0,nc)
     Call Memory_comp (Rating,nbf, nc,mr0, mr,mb)

     Allocate(wtb(mb), Bal(0:mr,2,mb))

     mrp= mr;  mrn= 0;  mt= 1;  mtp= 1  ! Default Rating < 1 values

     If (Rating >= 0) then
       If (Rating == Rating0) then
         mt= mt0;  Allocate(pt_val(mt));  pt_val= pt_val0
         mtp= Last_true(pt_val > 0)
       Else
         Call Point_values (Rating,Max_pt,mr, mt,mtp, pt_val)
       End if

       If (Rating == 1) then
         mrn= Ceiling(Neg_frac1 * mr);  mrp= mr - mrn
       End if
     End if

!    Now do the second reading of the ballots, storing the data in 'Bal'
     
     Close(9)
     Label= "Elections0\"//District
     Open(9, File= Trim(Label), IOstat=ios, Status='Old', Action='Read')
     Read(9,*);  Read(9,*)
     
     Call Read2 (nbf,bpm(6:),Rating0,Rating,UVP, nc0,nc,mr,mrp,&
                 mt, inv_orig,pt_val, tot_wt,nb, wtb,Bal)
     Close(9)

     Call Out ("Parameters for district "//District)
     Write (8, '(A,2I3,2I10)') "Rating, np, nb: ", Rating, np, nb
     Write (8, '(A,2I3, 2X,2I3, 2X,2I3)') "nc,nc0, mr,mrp, mt,mtp: ", &
                                           nc,nc0, mr,mrp, mt,mtp

!    Reorder ballots by ballot weight
     
     Allocate(key(nb));  Call Sort (.false.,wtb(:nb), key)

     Forall(i=0:mr) Bal(i,1,:nb)= Bal(i,1,key(:nb))
     Forall(i=0:mr) Bal(i,2,:nb)= Bal(i,2,key(:nb))
     
!    Write pre-processed ballot data

     If (Rating >= 1) then
       Label= "Elections123\"//District//'3'
       Open(9, File= Trim(Label), IOstat=ios, &
            Status='Replace', Action='Write')

     Else if (Rating == 0) then
       Label= "Elections123\"//District//'2'
       Open(9, File= Trim(Label), IOstat=ios, &
            Status='Replace', Action='Write')
     Else
       Label= "Elections123\"//District//'1'
       Open(9, File= Trim(Label), IOstat=ios, &
            Status='Replace', Action='Write')
     End if
       
     Write(9, '(A,3X,A,I6)')  District, Region, Year
     Write(9, '(A, I3, 3I10)') "'Rating type, # raw ballots, # processed, final # combined: '", &
                                 Rating, tot_cnt, nbf, nb
     Write(9, '(A, 3I3)')  "'# candidates to be elected, # reduced, original #: '", np, nc, nc0
     Write(9, '(A, 3I3)')  "'Max # ranked or rated, max # rating levels, undervote penalty: '", mr, mt, UVP

     Write(9, '(A,22I3)')  "'Original candidates, as reordered: '", orig(:nc)

     If (Rating == 0 .and. mt >= 4) then
       Write(9, '(A, 10F7.2)') "'Weak ranking point values: '", pt_val

     Else if (Rating == 1) then
       Write(9, '(/A)')        "'           Rating data: '"
       Write(9, '(A, 2I3)')    "'Max # pos & neg ratings: '", mrp, mrn 
       Write(9, '(A, 10F7.2)') "'Rating point values: '", pt_val
     End if
      
     If (Rating >= 0) then
       Do b= 1,nb
         Write(9, '(I7,2X, 20I3)') wtb(b), Bal(:,1,b), Bal(:,2,b)
       End do
     Else
       Do b= 1,nb
         Write(9, '(I7,2X, 10I3)') wtb(b), Bal(:,1,b)
       End do
     End if

     Write(9,'(I3)') 0;  Write(9,'(A)')  District
     Close(9)
     
   End Subroutine Read_bal

                           
  Subroutine Read0 (bpm,Rating0,UVP, nc0,mr0,mt0, n0,pt_val0, &
                    nbf,tot_cnt, Borda0,BordaP0,BordaN0)
  
!   Initial read of the ballot data, one ballot per line.

    Integer,  Intent(in) :: bpm(0:)    ! (0:3) 0 = bal_fmt, 1 = num, 2 = dup, 3 = cnt
                                       ! num = 1 for ballot numbering,     < 1 for no numbering
                                       ! dup = 1 for # duplicate ballots = 'wt', < 1 for assumed value of 1
                                       ! cnt = 1 for # ranked or rated candidates = 'nr', < 1 for no value
    Integer,  Intent(in) :: Rating0    ! Input file rating type:  -1 = strong ranking, 0 = weak ranking, 1 = rating
    Integer,  Intent(in) :: UVP        ! "Under Vote Penalty", or scale factor < 1, on ballot weights.
                                       !   UVP = 1 means no penalty. In general the scale factor is nr/UVP
                                       !   when nr = # ranked or rated < UVP, and no penalty for nr >= UVP.
    Integer,  Intent(in) :: nc0        ! # original candidates
    Integer,  Intent(in) :: mr0        ! max # candidates to be ranked or rated
    Integer,  Intent(in) :: mt0        ! # rating or ranking levels for Rating0 >= 0
    
    Integer,  Intent(in) :: n0         ! Prescribed # ballots in the input file, if > 0
    Real,     Intent(in) :: pt_val0(:) ! (mt0) Decreasing point values matching the increasing rating 
                                       !       levels (for Rating = 1), or ranking levels 
                                       !       (for Rating = 0 with mt0 = mr0)
    
    Integer, Intent(out) :: nbf        ! # valid ballots read
    Integer, Intent(out) :: tot_cnt    ! # original ballots represented
    Real,    Intent(out) :: Borda0(:)  ! (nc0) = BordaP0 + BordaN0 for Rating = 1, else = BordaP0
    Real,    Intent(out) :: BordaP0(:) ! (nc0) Borda tally of the rankings or positive ratings (Rating >= 0) 
    Real,    Intent(out) :: BordaN0(:) ! (nc0) Borda tally of the negative ratings (Rating = 1 only) 
    
! Local:
    Integer :: line(nc0), key(nc0), cnd(3+nc0), zl(3+nc0)
    Logical :: Last
    Integer :: bp(3), cand(nc0), rate(nc0)
    Real    :: wt, tot_wt, pts(nc0), rank_pt(mr0)
    Integer :: b, i, k, i2, n1, n2, nr, nw, bal_fmt
 
!   Read in each ballot, line by line for bal_fmt formatted ballots
    
    Call Rank_pts0 (0.50,Max_pt, rank_pt);  bal_fmt= bpm(0) 

    Borda0= 0;  BordaP0= 0;  BordaN0= 0
    b= 0;  tot_cnt= 0;  tot_wt= 0
    
    Ballot_loop : Do

      bp= bpm(1:)  
      Call Read_vote (b+1,nc0,mr0, bal_fmt,bp, Last, nr,cand,rate, &
                      line,key,cnd,zl)
      If (Last) Exit Ballot_loop  ! Reached end of data

      nw= bp(2)
      If (nr < 1 .or. nr > nc0 .or. nw < 1) then
        i= -1;  Cycle Ballot_loop  ! Invalid ballot
      End if

!     Accumulate the Borda tallies
      
      tot_cnt= tot_cnt + nw
      pts= -1;  wt= UVpen(UVP_fac,UVP,nr) * nw

      If (Rating0 == 1) then      ! Rating case
        pts(:nr)= pt_val0(rate(:nr))

        n1= Last_true(pts(:nr) > 0);  n2= n1 + 1
        i2= First_true(n2, pts(n2:nr) < 0)

        If (n1 >= 1)  BordaP0(cand(:n1))  = BordaP0(cand(:n1))   + wt * pts(:n1)
        If (i2 <= nr) BordaN0(cand(i2:nr))= BordaN0(cand(i2:nr)) + wt * pts(i2:nr)
        
      Else if (Rating0 == 0 .and. mt0 >= 4) then ! Weak ranking with pt_val
        pts(:nr)= pt_val0(rate(:nr))
        BordaP0(cand(:nr))= BordaP0(cand(:nr)) + wt * pts(:nr)
        
      Else  ! Strong ranking or weak without pt_val
        pts(:nr)= rank_pt(rate(:nr))
        Borda0(cand(:nr))= Borda0(cand(:nr)) + wt * pts(:nr)
      End if
      
      b= b + 1;  tot_wt= tot_wt + wt
    End do Ballot_loop

    nbf= b
    
    If (Rating0 == 1) then
      BordaP0= BordaP0 / tot_wt
      BordaN0= BordaN0 / tot_wt
      Borda0 = BordaP0 + BordaN0
    Else if (Rating0 == 0 .and. mr0 == mt0) then
      BordaP0= BordaP0 / tot_wt
      Borda0 = BordaP0
    Else
      Borda0= Borda0 / tot_wt
    End if
    
    Call Out ("Total # original ballots", tot_cnt, "with UVP weight",tot_wt, ln=1)
    Call Out ("# ballots read",nbf, "vs prescribed #",n0)

    Call Out ("Unordered overall Borda tally",Borda0)

    If (Rating0 == 1) then
      Call Out ("Borda tally of the positive ratings",BordaP0)
      Call Out ("Borda tally of the negative ratings",BordaN0)
    Else if (Rating0 == 0 .and. mr0 == mt0) then
      Call Out ("Borda tally of the weak rankings",BordaP0)
    End if
  End Subroutine Read0
  
  
  Subroutine Marginal_cand (Rating,Rating0,nc0, Borda0,BordaP0,BordaN0, &
                            nc,orig,inv_orig)
      
!   Eliminate marginal candidates

    Integer,    Intent(in) :: Rating       ! Output rating type: -1 / 0 = strong / weak ranking, 1 = rating
    Integer,    Intent(in) :: Rating0      ! Input rating type
    Integer,    Intent(in) :: nc0          ! original # candidates 
    Real,    Intent(inout) :: Borda0(:)    ! (nc0) = BordaP0 + BordaN0, if Rating0 = 1, else Borda for rankings
    Real,    Intent(inout) :: BordaP0(:)   ! (nc0) Borda positive ratings, if Rating0 = 1
    Real,    Intent(inout) :: BordaN0(:)   ! (nc0) Borda negative ratings, if Rating0 = 1

    Integer,   Intent(out) :: nc           ! # candidates to be retained
    Integer,   Intent(out) :: orig(:)      ! (nc0=>nc) Candidate ordering for 'Borda0' (decreasing)
    Integer,   Intent(out) :: inv_orig(:)  ! (nc0) Inverse of 'orig'. An original candidate mapped to 0 is deleted.
! Local:
    Logical :: Covered(nc0),  Cov0(nc0)
    Real    :: Bp1(nc0), Bn1(nc0)
    Real    :: cut, del, x1p, x1n
    Integer :: key(nc0), Removed(nc0)
    Integer :: i, n, nrm, lim

    lim= Lim_nc;  If (Rating >= 0) lim= Lim_nc - 1

    If (Rating0 == 1) then ! Rating data: remove marginal neutral candidates - 
                           ! when both BordaP0 and BordaN0 are near zero. 

      Bp1= BordaP0 / Maxval(BordaP0);  Bn1= BordaN0 / Minval(BordaN0)

      Covered= Bp1 > Marg_cut0 .or. Bn1 > Marg_cut0
      nc= Count(Covered)

      If (nc > lim) then
       cut= Marg_cut0;  del= .10 * cut
       n= Floor((1 - cut) / del)

        Do i= 1,n
            cut= cut + del
          Covered= Bp1 > cut .or. Bn1 > cut
          nc= Count(Covered);  If (nc <= lim) Exit
        End do
      End if

      Borda0= BordaP0 + BordaN0
      Call Sort (.false.,Borda0, key);  
      Cov0= Covered;  Covered= Covered(key)

      If (nc < nc0) then
        Removed= -1;    Call List_of_true (.not.Covered, nrm,Removed)
        orig= -1;  Call List_of_true (Covered, nc,orig)

        orig(:nc)= key(orig(:nc))
        Removed(:nrm) = key(Removed(:nrm))
      Else
        nrm= 0;  orig= key
      End if

    Else ! Ranking data: remove candidates of low Borda0

      Call Sort (.false.,Borda0,orig)
      Bp1= Borda0 / Borda0(1)

      cut= Marg_cut0;  nc= Last_true(Bp1 >= cut)

      If (nc > lim) then
       del= .10 * cut;  n= Floor((1 - cut) / del)

        Do i= 1,n
          cut= cut + del;  nc= Last_true(Bp1 >= cut)
          If (nc <= lim) Exit
        End do
      End if

      nrm= nc0 - nc;  Removed(:nrm)= orig(nc+1:)
    End if
    
    orig(nc+1:)= -1
    Call Inverse_map (orig(:nc), inv_orig)

    Call Out ("Final # candidates",nc, "out of",nc0, ln=1)
    Call Out ("List of candidates who were removed",Removed(:nrm))
    Call Out ("The original Borda candidate ordering",orig(:nc))
    Call Out ("The corresponding ordered Borda values", Borda0(:nc))
    Call Out ("The inverse of the candidate ordering",inv_orig)

    If (Rating0 == 1) then
      Call Out ("Normalized positive Borda values",  Bp1)
      Call Out ("Normalized negastive Borda values", Bn1)
      Call Out ("Retained clusters", Cov0)
    End if

  End Subroutine Marginal_cand


  Subroutine Read2 (nbf,bpm, Rating0,Rating,UVP, nc0,nc,mr,mrp,mt, &
                    inv_orig,pt_val, tot_wt,nb,wtb,Bal)
  
!   Read in the ballot data, one line per ballot. 
!   Adjust ratings as needed so that ratings average out to 0 while weak ranings 
!   sum to match the strong rankings sum

    Integer,  Intent(in) :: nbf         ! # ballots in the file
    Integer,  Intent(in) :: bpm(0:)     ! (0:3) 0 = bal_fmt, 1 = num, 2 = dup, 3 = cnt
                                        ! num = 1 for ballot numbering,  < 1 for no numbering
                                        ! dup = 1 for # duplicate ballots = 'wt', < 1 for assumed value of 1
                                        ! cnt = 1 for # ranked or rated candidates = 'nr', < 1 for no value
    Integer,  Intent(in) :: Rating0     ! Input voting file Rating type: -1 = strong ranking, 0 = weak ranking, 1 = rating
    Integer,  Intent(in) :: Rating      ! Output voting file rating type: -1 = strong ranking, 0 = weak ranking, 1 = rating,
                                        !   Rating >= Rating0 required
    Integer,  Intent(in) :: UVP         ! "Under Vote Penalty", or scale factor < 1, on ballot weights.
                                        !   UVP = 1 means no penalty. In general the scale factor is nr/UVP
                                        !   when nr = # ranked or rated < UVP, and no penalty for nr >= UVP.
    Integer,  Intent(in) :: nc0         ! # original candidates
    Integer,  Intent(in) :: nc          ! Revised # candidates
    Integer,  Intent(in) :: mr          ! Max # ranked or rated to be used upon output
    Integer,  Intent(in) :: mrp         ! Max # positively rated or ranked to be used upon output
    Integer,  Intent(in) :: mt          ! # rating levels for output file (for Rating == 1)

    Integer,  Intent(in) :: inv_orig(:) ! (nc0) Inverse of original Borda ordering 'orig'
    Real,     Intent(in) :: pt_val(:)   ! (mt) Decreasing point values of the 
                                        !   rating levels for Rating = 1, 
                                        !   or the ranking levels if Rating = 0 and mt = mr
    
    Real,    Intent(out) :: tot_wt      ! Total UVP weight of all ballots
    Integer, Intent(out) :: nb          ! # combined ballots recorded in 'Bal'
    Integer, Intent(out) :: wtb(:)      ! (nb) Raw ballot weights (# duplicates)
    Integer, Intent(out) :: Bal(0:,:,:) ! (0:mr,2,nb)  Records all ballot data
                                        !  (0,1,b)  = 'nr' = # candidates ranked or rated
                                        !  (0,2,b)  = 'lp' = # pos candidates ranked or rated
                                        !  (i,1,b) = ith ranked or rated candidate: 1..nr
                                        !  (i,2,b) = corresponding ranking or rating level (nondecreasing)
! Local:
    Integer :: ls(100), line(nc0), cnd(3+nc0), key(nc0), zl(3+nc0)
    Real    :: rnd(mr), frac_nr(mr,2)
    Real    :: wt, mt_fac, pt_sum, mean_num(2), dif
    Logical :: Last
    Integer :: bal0(0:mr,2), cand(nc0), rate(nc0), bp(3)
    Integer :: nw, mrn, mtp, bal_fmt
    Integer :: b, i, k, n, lp, ln, mb, mx, nr
 
    If (Rating > Rating0) mt_fac= Real(mt)/mr

    If (Rating == 1) then
      mtp= Last_true(pt_val > 0)
      mrn= mr - mrp;  mx= mt
    Else
      mx= mr
    End if

    bal_fmt= bpm(0);  mb= Ubound(Bal,3);  wtb= 0;  Bal= 0;  b= 0
    
    Ballot_loop :  Do k= 1,nbf
      bp= bpm(1:)  

      Call Read_vote (k,nc0,mr, bal_fmt,bp, Last, nr,cand,rate, &
                      line,key,cnd,zl)

      If (Last) Exit Ballot_loop  ! Reached end of data

 !    Convert to Borda ordered coordinates

      cand(:nr)= inv_orig(cand(:nr))

      If (Any(cand(:nr) < 1)) then
        n= nr;  Call List_of_true (cand(:n) > 0, nr,key)
        cand(:nr)= cand(key(:nr))
        rate(:nr)= rate(key(:nr))
      End if

!     Test for empty ballot

      nw= bp(2);  If (nw < 1 .or. nr < 1) Cycle Ballot_loop
      
!     If necessary, add a random candidate to get 2 ranked or 
!     rated candidates to approximate the sum restriction

      If (Rating >= 0 .and. nr == 1) then
        Call Single_rating (Rating, nc,mr,pt_val, cand(:2), rate(:2))
        lp= 1;  nr= 2
      End if
      
!     Add random perturbations

      If (Rating > Rating0) then
        rnd(:nr) = Ran_Gaus(nr)
        rate(:nr)= Nint(mt_fac*rate(:nr) + rand_sig*rnd(:nr))
        
        Where (rate(:nr) < 1)  rate(:nr)= 1
        Where (rate(:nr) > mx) rate(:nr)= mx
        Call Sort (.true.,rate(:nr),key(:nr))
        cand(:nr)= cand(key(:nr))

        If (Rating == 1) then  ! min lev must be pos, max must be neg
          rate(1)= Min(rate(1),mtp);  rate(nr)= Max(rate(nr),mtp+1)
          lp= Last_true(rate(:nr) <= mtp)
        Else
          lp= nr
        End if
      End if
      
!     For rating data, adjust the ratings to get the points
!     to average out to zero, at least approximately. For weak
!     ranking data, adjusts the ranking evels so that they average
!     out to the mean overall ranking level (mr+1)/2.
      
      If (Rating == 1) then  ! Rating
        Call Adjust_rating (k,b+1,mt,mtp,mrp,mrn,pt_val, nr, &
                            rate(:nr),cand(:nr), lp,pt_sum)

      Else if (Rating == 0) then  ! Weak ranking
        nr= Min(nr,mr);  lp= nr
        Call Adjust_ranking (k,nr,mt,pt_val, rate(:nr), dif)

      Else                        ! Strong ranking
        nr= Min(nr,mr);  lp= nr;  rate(:nr)= "ID"
      End if
      
      bal0= 0;  bal0(0,1)= nr;  bal0(1:nr,1)= cand(:nr)
      If (Rating >= 0) then
        bal0(0,2)= lp;  bal0(1:nr,2)= rate(:nr)
      End if 

!     Test for duplicate ballot. If found add the new ballot weight to the prior weight,
!     else record the new ballot.
      
      If (Rating >= 0) then
        Check_loop2 : Do i= 1,b
          If (Any(Bal(:,:,i) /= bal0)) Cycle Check_loop2

          wtb(i)= wtb(i) + nw;  Cycle Ballot_loop
        End do Check_loop2

        If (b >= mb) then
          Call Out ("Warning in 'Read1': Out of memory");  Exit Ballot_loop
        Else
          b= b + 1;  wtb(b)= nw;  Bal(:,:,b)= bal0
        End if
      Else
        Check_loop1 : Do i= 1,b
          If (Any(Bal(:,1,i) /= bal0(:,1))) Cycle Check_loop1

          wtb(i)= wtb(i) + nw;  Cycle Ballot_loop
        End do Check_loop1

        If (b >= mb) then
          Call Out ("Warning in 'Read1': Out of memory");  Exit Ballot_loop
        Else
          b= b + 1;  wtb(b)= nw;  Bal(:,1,b)= bal0(:,1)
        End if
      End if

!     New ballot

      If (b >= mb) then
        Call Out ("Warning in 'Read1': Out of memory");  Exit Ballot_loop
      End if
      
    End do Ballot_loop

    nb= b

!   Under voting penalties on the consolidated ballots

    If (UVP > 1) then
      Do b= 1,nb
        lp= Bal(0,2,b)   ! # pos ranked or rated
        nw= wtb(b)       ! Ballot weight = # duplicates

        wt= UVpen(UVP_fac,UVP,lp) * nw
        wtb(b)= Max(Nint(wt),1)
      End do
    End if

    tot_wt= Sum(wtb(:nb))

    Call Out ("Final # combined ballots",nb, "vs initial limit",nbf, ln=1)
    Call Out ("Total UVP ballot weight",tot_wt)

    mean_num(1)= Sum(Bal(0,1,:nb) * wtb(:nb)) / tot_wt
    Do n= 1,mr
     frac_nr(n,1)= Sum(wtb(:nb), Bal(0,1,:nb) == n) / tot_wt
    End do

    Call Out ("Mean # candidates ranked or rated per ballot",mean_num(1))
    Call Out ("Fraction of ballots for each # ranked or rated",frac_nr(:,1))

    If (Rating == 1) then
      mean_num(2)= Sum(Bal(0,2,:nb) * wtb(:nb)) / tot_wt
      Do n= 1,mr
        frac_nr(n,2)= Sum(wtb(:nb), Bal(0,2,:nb) == n) / tot_wt
      End do

      Call Out ("Mean # candidates positively ranked or rated per ballot",mean_num(2))
      Call Out ("Fraction of ballots for each # positively ranked or rated",frac_nr(:,2))
    End if

  End Subroutine Read2


  Subroutine Read_vote (b, nc,mr, bal_fmt, bpm, Last, nr,cand,rate, &
                        line,key,cnd,zl)

!   Read in a ballot from a file, using one of two formats,
!   yielding a list of 'nr' ranked or rated candidates 'cand'
!   with corresponding ranking or rating levels 'rate'

!   Note that the candidates will be ordered according to 
!   their increasing ranking or rating levels. If 2 or more
!   successive candidates have the same level, they will
!   have their natural order: by increasing candidate #.
 
  
    Integer,    Intent(in) :: b         ! ballot #
    Integer,    Intent(in) :: nc        ! total # candidates
    Integer,    Intent(in) :: mr        ! max # candidate ranking or ratings to be used
    Integer,    Intent(in) :: bal_fmt   ! Ballot format
                                        !   1 = variable # candidates in rank order
                                        !   2 = rate or rank level or 0, for each of 'nc' candidate

    Integer, Intent(inout) :: bpm(3)    ! 1 = num, 2 = dup, 3 = cnt on output, with
                                        ! < 1 on input if not to be read
    Logical,   Intent(out) :: Last      ! True for end of ballots in the file
    Integer,   Intent(out) :: nr        ! actual # candidate ranking or ratings to be used
    Integer, Intent(inout) :: cand(:)   ! (nc=>nr) List of the ranked or rated candidates in ranking
                                        !   or rating order
    Integer, Intent(inout) :: rate(:)   ! (nc=>nr) List of the matching ranking or rating levels, 
                                        !   by increasing levels, corresponding to decreasing points. 

    Integer, Intent(inout) :: line(:)   ! (nc)
    Integer, Intent(inout) :: key(:)    ! (nc)
    Integer, Intent(inout) :: cnd(:)    ! (3+nc)
    Integer, Intent(inout) :: zl(:)     ! (3+nc)
! Local:
!    Integer :: line(nc), key(nc)
    Integer :: i, j, n, n0, n1, ier
    
    Last= .false.;  nr= 0;  cand= 0;  rate= 0;  line= -1;  key= -1

    If (bal_fmt == 1) then
      Call Var_siz_record (b,bal_fmt,9,mr, bpm,n,line,ier, cnd,zl)
    Else
      Call Var_siz_record (b,bal_fmt,9,nc, bpm,n,line,ier, cnd,zl)
    End if

    Last= n < 1 .or. ier < 0;  If (Last) Return

    bpm(2)= Max(bpm(2),1)  !  Default # duplicates = 1

    If (bal_fmt == 2) then
      Call List_of_true (line(:nc) > 0, n0,cand)  ! Ordered by increasing candidate #
      Last= n0 < 1;  If (Last) Return
      
        rate(:n0)= line(cand(:n0))           ! Matching ranking or rating levels 
      Call Sort (.true.,rate(:n0), key(:n0)) ! This is a 'stable' sort: 'key' preserves
                                             !   the order for equal ranking or rating levels
      nr= Min(n0,mr)                         
      cand(:nr)= cand(key(:nr))              ! Natural candidate order is not changed for equal levels
    Else
      n0= n;  nr= Min(n0,mr);  cand(:nr)= line(:nr);  rate(:nr)= "ID"
    End if

    n1= bpm(3)
    If (n1 > 0 .and. n1 /= n0) then
      Call Out ("Warning in 'Read_vote': Mismatch in # ranked or rated")
      Call Out ("Specified #",n1, "Actual #",n0)
    End if

  End Subroutine Read_vote

   Subroutine Var_siz_record (b,form,unit,nc, v3,n,line, ier, cnd,zl)

!    Read a line from file 'unit', non-advancing, of up to 3 integers
!    of ballot data (num, dup, cnt), as specified by v3 > 0. 

!    Follow by reading in 'n' integers (ordered candidate indices 
!    when 'form' = 1) or 'nc' integers (ranking or rating level or 0
!    for each candidate index, when 'form' = 2).

!    When 'form' = 1, n = v3(3) if read in as > 0, else 'n' is 
!    determined by non-advancing input until either the 
!    terminating integer 0 is reach or an end-of-record.

!    This read is always advanced to the next line on output.

     Integer,     Intent(in) :: b        ! ballot #
     Integer,     Intent(in) :: form     ! Line format type. (1) Variable # candidates
                                         !                   (2) Fixed # candidates
     Integer,     Intent(in) :: unit     ! Input file unit #

     Integer,     Intent(in) :: nc       ! max # candidates in 'line'
     Integer,  Intent(inout) :: v3(:)    ! (3) initial data: v3(i) < 1 on input if 
                                         !   data 'i' is not present. The 'n3' data
                                         !   that are present read in to the first
                                         !   n3 positions of v3.
     Integer,    Intent(out) :: n        ! # integers read into 'line' (beyond 'v3')
     Integer,  Intent(inout) :: line(:)  ! (nc) Candidate data
     Integer,    Intent(out) :: ier      ! Error if < 0

     Integer,  Intent(inout) :: cnd(:)   ! (3+nc)
     Integer,  Intent(inout) :: zl(:)    ! (3+nc)
! Local:
!     Integer :: zl(3+nc), cnd(3+nc)
     Integer :: i, k, m, n3, ios, lst(3), vl(3)

!    Initial ballot data in 'v3'

     n= 0;  line= -1; ier= 1
     lst= -1;  vl= -1;  zl= -1;  cnd= -1

     Call List_of_true (v3 > 0, n3,lst)

     If (n3 > 0) then
       Read (unit,*, IOstat=ios) vl(:n3)

       If (ios /= 0 .or. Any(vl(:n3) < 1)) then
         Return
       End if

       v3(lst(:n3))= vl(:n3)  ! Output these values
       n= Max(v3(3), 0);  v3(2)= Max(v3(2),1)
     Else
       Return
     End if

!    The ballot rankings or ratings in 'line'

     If (form == 1) then ! Variable line length. Line may not include zeros before 'n'
       m= n3 + n

       If (n > 0) then   ! Read first 'n' candidate indices
         Backspace (unit, IOstat=ios);  Read(unit,*) cnd(:m)
         line(:n)= cnd(n3+1:m)

       Else  ! 'n' unknown, so do non-advancing input to determine it
          Backspace (unit, IOstat=ios);  m= n3+nc;  cnd= -1

         Read(unit,'(7I5)', Advance='No', Size=i, EOR=11, Err= 11, IOstat=ios) cnd(:m)

 11      k= Last_true_in(cnd > 0)
         n= k - n3;  line(:n)= cnd(n3+1:k)

         If (ios == 0) then
           Read(unit,*) ! Advance to the next line
         End if
       End if
     Else  ! Fixed line length. Line may include zeros
       n= nc;  m= n3 + n 
       Backspace (unit, IOstat=ios);  Read(unit,*) cnd(:m)
       line(:n)= cnd(n3+1:m)
     End if

   End Subroutine Var_siz_record
    
    Pure Subroutine Adjust_ranking (b,nr,mt,pt_val, rate, dif)
    
!     With ratings, adjust them so that their sum of point values approximates 0.
    
      Integer,    Intent(in) :: b          ! ballot #
      Integer,    Intent(in) :: nr         ! # candidates ranked
      Integer,    Intent(in) :: mt         ! # ranking levels
      Real,       Intent(in) :: pt_val(:)  ! (mt) Descending point values for the 
                                           !      rating levels 1...mt

      Integer, Intent(inout) :: rate(:)    ! (nr) Rating levels (increasing) for each candidate

      Real,      Intent(out) :: dif        ! Final point sum of differences
                                           ! the point average

! Local:
      Real, Parameter :: eps= 0.50
      Real    :: pts(nr), pt_gap(mt-1)
      Integer :: tmp(nr), rate0(nr)

      Logical :: Balanced, Improved
      Real    :: sm, mid, dif_up
      Integer :: i, k, n, l1, l2, ln, lp, r1, rv, mr1, mt1

      rate0= rate;  tmp= -1;  mid= Sum(pt_val) / mt
      
!     Initial difference point sum 'dif'

      pts= pt_val(rate) - mid;  dif= Sum(pts)

      Balanced= Abs(dif) <= eps;  If (Balanced) Return

!     Successive point decreases

      Forall(rv=1:mt-1) pt_gap(rv)= pt_val(rv) - pt_val(rv+1)

      If (dif > 0) then ! Lower point values by raising levels

        Loop1 : Do k= 1,mt-1

!         Look for the last level increase, from 'rv' to rv+1. 
!         Then raise this level.

          tmp(nr)= mt - rate(nr);  tmp(:nr-1)= rate(2:nr) - rate(1:nr-1)
          i= last_true(tmp(:nr) > 0);  rv= rate(i)

          If (i < nr) then
            dif_up= dif - pt_gap(rv);  r1= rv + 1
          Else
            dif_up= dif - (pt_val(rv) - pt_val(mt));  r1= mt
          End if

          Balanced= Abs(dif_up) < eps
          Improved= Abs(dif_up) <= Abs(dif)
            
          If (Balanced .or. Improved) then
            dif= dif_up;  rate(i)= r1
          End if
          If (Balanced .or. .not.Improved .or. dif_up < 0) Exit Loop1
        End do Loop1

      Else                 ! Raise point values by lowering levels

        Loop2 : Do k= 1,mt-1

!         Look for the first level increase 'i'. Then lower this level

          tmp(1)= rate(1) - 1;  tmp(2:nr)= rate(2:nr) - rate(1:nr-1)
          i= First_true(tmp(1:nr) > 0);  rv= rate(i)

          If (i > 1) then
            r1= rv - 1;  dif_up= dif + pt_gap(r1)
          Else
            dif_up= dif + (pt_val(1) - pt_val(rv));  r1= 1
          End if

          Balanced= Abs(dif_up) < eps
          Improved= Abs(dif_up) <= Abs(dif)
            
          If (Balanced .or. Improved) then
            dif= dif_up;  rate(i)= r1
          End if
          If (Balanced .or. .not.Improved .or. dif_up > 0) Exit Loop2
        End do Loop2
      End if

    End Subroutine Adjust_ranking

    
    Subroutine Adjust_rating (k,b,mt,mtp,mrp,mrn,pt_val, nr,rate,cand, lp,pt_sum)
    
!     With ratings, adjust them so that their sum of point values approximates 0.
     
      Integer,    Intent(in) :: k          ! initial ballot #
      Integer,    Intent(in) :: b          ! next consolidated ballot #
      Integer,    Intent(in) :: mt         ! # rating levels
      Integer,    Intent(in) :: mtp        ! # positive rating levels
      Integer,    Intent(in) :: mrp        ! Max # positive ratings permitted
      Integer,    Intent(in) :: mrn        ! Max # negative ratings permitted
      Real,       Intent(in) :: pt_val(:)  ! (mt) Descending point values for the 
                                           !      rating levels 1...mt
      Integer,    Intent(in) :: nr         ! # candidates rated

      Integer, Intent(inout) :: rate(:)    ! (nr) Rating levels (increasing)
      Integer, Intent(inout) :: cand(:)    ! (nr) Rating levels (increasing)

      Integer,   Intent(out) :: lp         ! # pos candidates rated
      Real,      Intent(out) :: pt_sum     ! Final point sum (target sum = 0)

! Local:
      Real, Parameter :: eps= 0.50
      Real    :: pts(nr), pt_gap(mt-1)
      Integer :: tmp(nr), rate0(nr)

      Logical :: Balanced, Improved
      Real    :: pt_up
      Integer :: i, j, l, l1, l2, lj, ln, n1, r1, rv, mr1, mt1

      rate0= rate;  tmp= -1

!     Initial point sum 'pt_sum'

      pts(:nr)= pt_val(rate(:nr));  pt_sum= Sum(pts(:nr))
      Balanced= Abs(pt_sum) <= eps

!     Successive point decreases

      Forall(rv=1:mt-1) pt_gap(rv)= pt_val(rv) - pt_val(rv+1)

      If (.not.Balanced) then
       If (pt_sum > 0) then ! Lower point values by raising levels

        Loop1 : Do j= 1,mt-1

!         Look for the last level increase, from 'rv' to rv+1. 
!         Then raise this level.

          tmp(nr)= mt - rate(nr);  tmp(:nr-1)= rate(2:nr) - rate(1:nr-1)
          i= last_true(tmp(:nr) > 0);  rv= rate(i)

          If (i < nr) then
            pt_up= pt_sum - pt_gap(rv);  r1= rv + 1
          Else
            pt_up= pt_sum - (pt_val(rv) - pt_val(mt));  r1= mt
          End if

          Balanced= Abs(pt_up) < eps
          Improved= Abs(pt_up) <= Abs(pt_sum)
            
          If (Balanced .or. Improved) then
            pt_sum= pt_up;  rate(i)= r1
          End if
          If (Balanced .or. .not.Improved .or. pt_up < 0) Exit Loop1
        End do Loop1

       Else                 ! Raise point values by lowering levels

        Loop2 : Do j= 1,mt-1

!         Look for the first level increase 'i'. Then lower this level

          tmp(1)= rate(1) - 1;  tmp(2:nr)= rate(2:nr) - rate(1:nr-1)
          i= First_true(tmp(1:nr) > 0);  rv= rate(i)

          If (i > 1) then
            r1= rv - 1;  pt_up= pt_sum + pt_gap(r1)
          Else
            pt_up= pt_sum + (pt_val(1) - pt_val(rv));  r1= 1
          End if

          Balanced= Abs(pt_up) < eps
          Improved= Abs(pt_up) <= Abs(pt_sum)
            
          If (Balanced .or. Improved) then
            pt_sum= pt_up;  rate(i)= r1
          End if
          If (Balanced .or. .not.Improved .or. pt_up > 0) Exit Loop2
        End do Loop2
       End if
      End if

!     Check for # pos ratings lp <= mrp and # negative ln <= mrn

      lp= Last_true(rate(:nr) <= mtp)  ! # pos rating levels
      ln= nr - lp                      ! # neg rating levels

      mr1= mrp + 1;  mt1= mtp + 1

      If (lp > mrp .and. ln < mrn) then ! Reduce pos rating levels by a swap
        j= Min(lp - mrp, mrn - ln);  rate(mr1:mrp+j)= mt1
        lp= lp - j;  ln= ln + j
        
        pt_sum= Sum(pt_val(rate(:nr)))

        If (pt_sum < -4) then
          l= First_true(rate(:lp) > 1)
          If (l > 0) then
            rv= rate(l);  rate(l)= rv - 1
            pt_sum= pt_sum +  pt_gap(rate(l))
          End if 
        End if

        If (pr_out > 1) then
          Write(8,'(A, 2I5, 2I3, F7.2, 15I3)') "Lower lp: k,b,lp,ln, pt_sum, rate0, rate", k,b,lp,ln, pt_sum, rate0, rate
        End if
      End if

      If (ln > mrn .and. lp < mrp) then  ! Reduce neg rating levels by a swap
        j= Min(ln - mrn, mrp - lp);  lj= lp + j;  rate(lp+1:lj)= mtp
        ln= ln - j;  lp= lp + j

        pt_sum= Sum(pt_val(rate(:nr)))

        If (pt_sum > 4) then
          l2= lj + 1;  l= lj + Last_true(rate(l2:nr) < mt)

          If (l > lj) then
            rv= rate(l);  rate(l)= rv + 1
            pt_sum= pt_sum -  pt_gap(rv)
          End if
        End if

        If (pr_out > 1) then
          Write(8,'(A, 2I5, 2I3, F7.2, 15I3)') "Raise lp: k,b,lp,ln, pt_sum, rate0, rate", k,b,lp,ln, pt_sum, rate0, rate
        End if
      End if

      If (lp > mrp .or. ln > mrn) then
        Call Out ("Warning in 'Adjust_rating': Bad lp",lp, "or ln",ln, ln=1)
      End if

      pts(:nr)= pt_val(rate(:nr));  pt_sum= Sum(pts(:nr))

    End Subroutine Adjust_rating

    Pure Subroutine Single_rating (Rating, nc,mr,pt_val, cand, rate)

!     If a Rating = 1 ballot has rated only one candidate, then add 
!     a different candidate, usually the last one, whose rating 
!     level is chosen so that the sum of the 2 ratings is close to 0.
!     Thus if one rating is positive, the other will be negative.

!     In the case of weak rankings (Rating = 0), the ranking level 
!     of the 2nd candidate is chosen so that the 2 ranking levels 
!     average to approximate the overall average ranking level (mr+1)/2.

      Integer,     Intent(in) :: Rating      ! Rating type = 1 or 0
      Integer,     Intent(in) :: nc          ! # candidates
      Integer,     Intent(in) :: mr          ! # ranking levels, for Rating = 0
      Real,        Intent(in) :: pt_val(:)   ! (mt) rating point values, for Rating = 1

      Integer,  Intent(inout) :: cand(:)     ! (2) The voted candidates, 1 input, 2 output
      Integer,  Intent(inout) :: rate(:)     ! (2) Their corresponding rating or ranking leves
!   Local:
      Integer :: cn1, cn2, rt1, key(2)

      cn1= cand(1);  rt1= rate(1)

      cn2= nc;  If (cn1 == nc) cn2= 1
      cand(2)= cn2

      If (Rating == 1) then
        rate(2)= Minloc(Abs(pt_val(rt1) + pt_val),1)
      Else if (Rating == 0) then
        rate(2)= mr+1 - rt1
      End if

      Call Sort (.true.,rate, key)
      cand= cand(key)

    End Subroutine Single_rating

    Subroutine Point_values0 (Rating0,Max_pt, mr0,mt0, pt_val0)
 
!     Adjust the rating point values 'pt_val' as read in for Rating0 = 1 
!     so that the neutral value is 0, the maximum rating is 10.0 and the 
!     minimum rating is >= -10.0, however excluding 0 as a point value.

!     Optionally, point values may be read into pt_val for weak rankings
!     in the Rating0 = 0 case (signaled by mt > 1), or if none are provided, 
!     then pt_val is assigned default values.

!     In the Rating0 = -1 case of strong rankings, pt_val is not used.
    
      Integer,    Intent(in) :: Rating0     ! File Rating type 
      Real,       Intent(in) :: Max_pt      ! Max possible point value
      Integer,    Intent(in) :: mr0         ! Max # ranked
      Integer, Intent(inout) :: mt0         ! # rating or weak ranking levels
      Real,          Pointer :: pt_val0(:)  ! (mt0) Point values read in for increasing rating or ranking levels.  
                                            !      Values must be decreasing, from positive to negative for ratings
                                            !      data, all positive for weak ranking data.  -1 if none.
!   Local:
      Real, Parameter :: eps= 0.10
      Logical :: Read_in, Error
      Real    :: fac
      Integer :: i, j, n, mp, n1, n2

      If (Rating0 >= 0 .and. mt0 >= 4) then
        Read_in= .true.
        Allocate (pt_val0(mt0));  Read(9,*) pt_val0
        
!       Point values must be strictly decreasing (> eps) from a positive value (> eps)

        Error= pt_val0(1) <= eps
        Do i= 2,mt0
          Error= Error .or. pt_val0(i) >= pt_val0(i-1) - eps
        End do          
        
!       For ratings the first 3 values mbust be positive and the last 2 negative.

        Error= Error .and. Rating0 == 1 .and. (pt_val0(3) < eps .or. pt_val0(mt0-1) > -eps)

!       For rankings all postive must be positive.

        Error= Error .and. Rating0 == 0 .and. pt_val0(mt0) <= eps

        If (Error) then
          Call Out ("Error in Point_values0, for voting data file with rating",Rating0, ln=1);  
          Call Out ("Invalid point values for rating or ranking levels",pt_val0)
          Stop
        End if

        If (Rating0 == 0) then
          fac= Max_pt / pt_val0(1);  If (fac /= 1) pt_val0= fac * pt_val0
        Else 
          mp= Last_true(pt_val0 > 0);  fac= Max_pt / pt_val0(1)
          If (fac /= 1) pt_val0(:mp)= fac * pt_val0(:mp)

          fac= -Max_pt / pt_val0(mt0)
          If (fac /= 1) pt_val0(mp+1:)= fac * pt_val0(mp+1:)
        End if
      Else
        Read_in= .false.
        Read(9,*)  ! Skip over any data in this line
      End if
      
!     Put in default rating or ranking points

      If (.not.Read_in) then
        mt0= mr0;  Allocate(pt_val0(mt0))

        If (Rating0 == 1) then
          n2= Ceiling(Neg_frac1*mt0);  n1= mt0 - n2
          Call Rank_pts0 (0.50,Max_pt, pt_val0(:n1))

          Call Rank_pts0 (0.50,Max_pt, pt_val0(n1+1:))
          pt_val0(n1+1:mt0)= -pt_val0(mt0:n1+1:-1)
        Else
          Call Rank_pts0 (0.50,Max_pt, pt_val0)
        End if
      End if
      
    End Subroutine Point_values0

    Subroutine Point_values (Rating,Max_pt,mr, mt,mtp, pt_val)
 
!     Assign default the rating point values to 'pt_val' when Rating = 1 
!     with 'mtn' of these positive and 'mtn' negative and mt= mtp + mtn.

!     Assign default the rating point values to 'pt_val' when Rating = 1 
!     with 'mtn' of these positive and 'mtn' negative and mt= mtp + mtn.

      Integer,    Intent(in) :: Rating     ! Output file Rating type, assuming Rating > Rating0
      Real,       Intent(in) :: Max_pt     ! Max possible point value
      Integer,    Intent(in) :: mr         ! Max # candidates to be ranked or rated

      Integer,   Intent(out) :: mt         ! # rating or weak ranking point levels
      Integer,   Intent(out) :: mtp        ! # positive ratings in pt_val if Rating = 1

      Real,          Pointer :: pt_val(:)  ! (mt) Point values read in for increasing rating or ranking levels.  
                                           !      Values must be decreasing, from positive to negative for ratings
                                           !      data, all positive for weak ranking data.  -1 if none.
!   Local:
      Integer :: mt0, mrn, mtn

!     Compute default rating point conversion values 'pt_val'
      
      If (Rating == 1) then
        mt= mr;  Allocate(pt_val(mt))
        mtn= Ceiling(Neg_frac1 * mt);  mtp= mt - mtn 

        Call Rank_pts0 (0.50,Max_pt, pt_val(:mtp))
        Call Rank_pts0 (0.50,Max_pt, pt_val(mtp+1:))
        pt_val(mtp+1:)= -pt_val(mt:mtp+1:-1)
       
        Call Out ("For total # rating levels",mt, ln=1)
        Call Out ("# positive rating levels",mtp, "# negative levels",mtn)
        Call Out ("Corresponding point values",pt_val)

      Else if (Rating == 0) then
        mt= mr;  mtp= mt;  Allocate(pt_val(mt))

        Call Rank_pts0 (0.50,Max_pt, pt_val)
       
        Call Out ("Total # ranking levels",mt, "is max # candidates to be ranked",mr, ln=1)
        Call Out ("Corresponding point values",pt_val)
      Else
        mt= 1;  mtp= 1;  Allocate(pt_val(1));  pt_val= -1
      End if

    End Subroutine Point_values

   Subroutine Memory_comp (Rating,nbf, nc,mr0, mr,mb)
  
!  Compute the maximum the number of ballots 'mb' that may be recorded one ballot 
!  structure: Bal(-1:nc,mb) for the given restritions on the number of candidates
!  and on the number of ranking or rating levels, adjusting the maximum number
!  of levels 'mr' as necessary, so that 'Bal' is allocatable. Take into account
!  the average restriction that will be imposed on the point values, which reduces
!  the dimension of the ranking / rating level space by 1 (and the # combinations
!  of rankings / ratings by more than this in general).
  
    Integer,  Intent(in) :: Rating  ! Rating vs ranking option
    Integer,  Intent(in) :: nbf     ! # counted consolidated ballots

    Integer,  Intent(in) :: nc      ! # candidates
    Integer,  Intent(in) :: mr0     ! Max # ranked or rated candidates

    Integer, Intent(out) :: mr      ! Max # ranked or rated candidates
    Integer, Intent(out) :: mb      ! Limit on # consolidated ballots
!  Local:
    Integer :: m, n, m1, m2, p1, p2, Alloc_err, mbl(mr0)

    Call Out ("Enter 'Memory_comp'")

    mbl= 0;  mb= 0;  p1= Floor(0.67*mr0);  p2= mr0 - p1

    Do m= 2,mr0
      If (Rating >= 0) then
        m1= Floor(0.667*m);  m2= m - m1;  m1= m1-1;  m2= m2-1
        n= N_subsets(nc,m) * N_subsets(m,m1) * (p1**m1) * (p2**m2)

        mb= mb + n
      Else 
        mb= mb + Factorial(nc,m)
      End if

      mbl(m)= mb
    End do

    Call Out ("Max # ballots for each max # ranked",mbl)

    mr= Max(Last_true(mbl < Lim_nb), 4)
    mb= Min(mbl(mr),nbf)

    If (mr < mr0) then
      Call Out ("Initial max # ranked or rated candidates",mr0, "reduced to",mr, ln=1)
    End if

    If (mb < nbf) then
      Call Out ("Initial # ballots",nbf, "reduced to # consolidated ballots",mb, ln=1)
    End if
    
  End Subroutine Memory_comp


  Pure Integer Function UVpen (UVP_fac, UVP, nr)
    Real,    Intent(in) :: UVP_fac(:,:) ! (5,5) Ballot weight 
                           ! scale factors that represent 
                           ! under vote penalties
    Integer, Intent(in) :: UVP ! Under vote penalties for 'nr' < UVP
    Integer, Intent(in) :: nr  ! # voted < UVP

    If (nr < UVP .and. nr > 0 .and. UVP <= Size(UVP_fac,2)) then
      UVpen= UVP_fac(nr,UVP)
    Else
      UVpen= 1.0
    End if
  End Function UVpen

End Module Clusters_pre