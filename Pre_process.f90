Program Pre_process 
    
!  This program does pre-processing on district ballot files to get them into a standard format 
!  for ranking or rating. It also may perform ballot consolidation, both to reduce the computational 
!  load by eliminating marginal candidates and to restrict the number of candidates ranked or rated 
!  by each ballot if that is specified by the user, or to reduce the computation.
 
!  For testing weak ranking in PR_clustering.f90, the strongly ranked ballots are perturbed 
!  so that sometimes more than 1 ballot may be assigned the same ranking.  Likewise, these 
!  strongly ranked ballots may be perturbed to produce ballots with ratings, given specified 
!  rating levels with point values from 10.0 to -10.0, which are adjusted to yield ballots 
!  whose rating point values average out close to 0.

   Use Clusters_pre
   Use Types
   Use Precisn
   Implicit None

   Integer :: Rating     ! Rating type: -1 = strong ranking, 0 = weak ranking, 1 = rating: For output files
   Integer :: mr_spec    ! Max # ranked or rated candidates to be used for the Rating type output
   Integer :: mt_spec    ! Max # rating levels to be used in pt_val for the Rating type (>= 0) output
   Integer :: UVP        ! "Under Vote Penalty", or scale factor < 1, on ballot weights.
                         !   UVP = 1 means no penalty. In general the scale factor is nr/UVP
                         !   when nr = # ranked or rated < UVP, and no penalty for nr >= UVP.
   Integer :: Dist1      ! First voting district to be processed
   Integer :: Dist2      ! Last  voting district to be processed
   Integer :: nDist      ! total # voting districts
   Character(18), Pointer :: District(:)=>Null() ! (nDist) The electoral district
   Character(16), Pointer :: Region(:)=>Null()   ! (nDist) The region of each district
   Integer, Pointer       :: Year(:)=>Null()     ! (nDist) The year of the election

   Character(80) :: msg
   Character(70) :: label
   Character(10) :: out_file
   Integer       :: i, id, ios, itmp(4)

!  Run input parameters and data. Original ballot data in 'Elections0'.
!  Other district data in 'Vote_files.txt'.

   label= "Pre_opt.txt"

!  Optionally process the ballots only from a range of districts (Dist1 to Dist2)
!  and limit the # ranked or rated to 'mr_spec' and # rating levels to 'mt_spec'

   Open(7, File= Trim(Label), IOmsg=msg, IOstat=ios, Status='Old', Action='Read')
     Read(7,*) label, out_file, pr_out, UVP, Dist1, Dist2
     Read(7,*) label, Rating, mr_spec, mt_spec
   Close(7)

   label= "Results\"//out_file

   Open(8, File= Trim(Label), IOmsg=msg, IOstat=ios, Status='Replace', Action='Write')
   
   label= "Elections0\"//"Vote_files.txt"

   Open(7, File= Trim(Label), IOmsg=msg, IOstat=ios, Status='Old', Action='Read')
     Read(7,'(I3)') nDist
     Allocate (District(nDist), Region(nDist), Year(nDist))
     Do id= 1,nDist
       Read(7,*) District(id), itmp, Region(id), Year(id)
       Read(7,*) 
     End do
   Close(7)

!  Read and preprocess the ballots in each district
   
   District_loop : Do id= Dist1,Dist2
       
     Call Read_bal (Rating,mr_spec,mt_spec,UVP, id,Trim(District(id)),Trim(Region(id)),Year(id))

   End do District_loop

   Close(8)

 End Program Pre_process   