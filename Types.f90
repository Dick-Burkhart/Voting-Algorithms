
Module Types
   Use Output
   Use Precisn
   Implicit None

!  A type for graph adjacency

   Type Adjacency
     Integer  :: nl = 0  ! Length of "ls"
     Integer  :: pos= 0  ! Next position in "ls" after "node" if ls lists the
                                         !   increasing nodes adjacent to "node".
     Integer, Pointer :: ls(:)=> Null()  ! List of indices
     Integer, Pointer :: vl(:)=> Null()  ! Integer edge values corresponding to "ls"
     Real,    Pointer :: wt(:)=> Null()  ! Real edge weights corresponding to "ls"
     
     Integer  :: sum_vl= 0  ! = Sum(vl) or other vertex value
     Real     :: sum_wt= 0  ! = Sum(wt) or other vertex value
     Real     :: fx= 0      ! Other real number

     Integer, Pointer :: sr(:)=> Null()  ! Other integer list
     Real,    Pointer :: wl(:)=> Null()  ! Other real list
   End Type Adjacency

   Type AdjacencyD
     Integer     :: nl = 0     ! Length of "ls"
     Integer     :: pos= 0     ! Next position in "ls" after "node" if ls lists the
                               !   increasing nodes adjacent to "node".
     Integer,    Pointer :: ls(:)=> Null()  ! List of indices
     Integer,    Pointer :: vl(:)=> Null()  ! Integer edge values corresponding to "ls"
     Real(Dblp), Pointer :: wt(:)=> Null()  ! Real edge weights corresponding to "ls"
     
     Integer     :: sum_vl= 0  ! = Sum(vl) or other vertex value
     Real(Dblp)  :: sum_wt= 0  ! = Sum(wt) or other vertex value

     Integer,    Pointer :: sr(:)=> Null()  ! Other integer list
     Real(Dblp), Pointer :: wl(:)=> Null()  ! Other real list
     Real(Dblp)  :: fx= 0      ! Other real number
   End Type AdjacencyD

!  A type for a list of sets which may have values and memberships

   Type Set_list
     Integer :: n= 0    ! # elements in 'set' or other value
     Integer :: p= 0    ! other value
     Integer :: r= 0    ! other value

     Integer, Pointer :: set(:)=> Null()  ! Elements of the set in increasing order
     Integer, Pointer :: lev(:)=> Null()  ! Associated levels
     Real,    Pointer :: val(:)=> Null()  ! Values corresponding to "set"
     Real,    Pointer :: mbr(:)=> Null()  ! Memberships corresponding to "set", or other values

     Real    :: svl= 0  ! = Sum(val) or other value
     Real    :: smb= 0  ! = Sum(mbr) or other value
   End Type Set_list

!  A multi-purpose list structure

   Type Multi_listR
     Integer  :: k= 0    
     Integer  :: l= 0    
     Integer  :: m= 0    
     Integer  :: n= 0    
     Integer  :: o= 0    
     Integer  :: p= 0    
     Integer  :: q= 0   
     Integer  :: sum_vl= 0
     
     Real     :: sum_wt= 0  
     Real     :: fsx= 0 
     Real     :: fux= 0 
     
     Integer,  Pointer :: ls(:)=>  Null() ! 
     Integer,  Pointer :: lt(:)=>  Null() ! 
     Integer,  Pointer :: vl(:)=>  Null() ! 
     
     Real,     Pointer :: wt(:)=>  Null() ! 
     Real,     Pointer :: px(:)=>  Null() ! 
     Real,     Pointer :: qx(:)=>  Null() ! 
     Real,     Pointer :: rx(:)=>  Null() ! 
     Real,     Pointer :: sx(:)=>  Null() ! 
     Real,     Pointer :: tx(:)=>  Null() ! 
     Real,     Pointer :: ux(:)=>  Null() ! 

     Integer,  Pointer :: L0(:,:)=>  Null() ! 
     Integer,  Pointer :: L1(:,:)=>  Null() ! 
     Integer,  Pointer :: L2(:,:)=>  Null() ! 
     Integer,  Pointer :: L3(:,:)=>  Null() ! 
     
     Real,     Pointer :: M0(:,:)=>  Null() ! 
     Real,     Pointer :: M1(:,:)=>  Null() ! 
     Real,     Pointer :: M2(:,:)=>  Null() ! 
     Real,     Pointer :: M3(:,:)=>  Null() ! 
     
     Integer,  Pointer :: Q0(:,:,:)=>  Null() ! 
     Integer,  Pointer :: Q1(:,:,:)=>  Null() ! 
     
     Real,     Pointer :: T0(:,:,:)=>  Null() ! 
     Real,     Pointer :: T1(:,:,:)=>  Null() ! 
   End Type Multi_listR

   Type List_of_MultiR
     Real     :: x= 0
     Real     :: y= 0
     Integer  :: m= 0
     Integer  :: n= 0             ! Size of 'ml'
     Type(Multi_listR), Pointer :: ml(:) =>  Null() ! A array of Multi_listR
   End Type List_of_MultiR

!  A list structure with double precision data

   Type Multi_listD
     Integer     :: k= 0            ! 
     Integer     :: l= 0            ! 
     Integer     :: m= 0            ! 
     Integer     :: n= 0            ! 
     Integer     :: o= 0            ! 
     Integer     :: p= 0            ! 
     Integer     :: q= 0            ! 
     Integer     :: sum_vl= 0       ! 
     
     Real(Dblp)  :: sum_wt= 0       !
     Real(Dblp)  :: fsx= 0          !
     Real(Dblp)  :: fux= 0          !
     
     Integer,    Pointer :: ls(:)=>  Null() !
     Integer,    Pointer :: lt(:)=>  Null() !
     Integer,    Pointer :: vl(:)=>  Null() !
     
     Real(Dblp), Pointer :: wt(:)=>  Null() !
     Real(Dblp), Pointer :: px(:)=>  Null() !
     Real(Dblp), Pointer :: qx(:)=>  Null() !
     Real(Dblp), Pointer :: rx(:)=>  Null() !
     Real(Dblp), Pointer :: sx(:)=>  Null() !
     Real(Dblp), Pointer :: tx(:)=>  Null() !
     Real(Dblp), Pointer :: ux(:)=>  Null() !
     
     Integer,    Pointer :: L0(:,:)=>  Null() ! 
     Integer,    Pointer :: L1(:,:)=>  Null() ! 
     Integer,    Pointer :: L2(:,:)=>  Null() ! 
     Integer,    Pointer :: L3(:,:)=>  Null() ! 
     
     Real(Dblp), Pointer :: M0(:,:)=>  Null() ! 
     Real(Dblp), Pointer :: M1(:,:)=>  Null() ! 
     Real(Dblp), Pointer :: M2(:,:)=>  Null() ! 
     Real(Dblp), Pointer :: M3(:,:)=>  Null() ! 
     
     Integer,    Pointer :: Q0(:,:,:)=>  Null() ! 
     Integer,    Pointer :: Q1(:,:,:)=>  Null() ! 
     
     Real(Dblp), Pointer :: T0(:,:,:)=>  Null() ! 
     Real(Dblp), Pointer :: T1(:,:,:)=>  Null() ! 
   End Type Multi_listD

   Type List_of_MultiD
     Real(Dblp)  :: x= 0d0
     Real(Dblp)  :: y= 0d0
     Integer     :: m= 0
     Integer     :: n= 0    ! Size of 'ml'
     Type(Multi_listD), Pointer :: ml(:) =>  Null() ! An array of Multi_listD
   End Type List_of_MultiD

!  Structure for linked lists, such an enumeration of the cliques of a graph.

   Type Linked
     Logical,      Pointer :: Comp(:)  => Null()  ! Complement of 'vertex'
     Integer,      Pointer :: vertex(:)=> Null()  ! List of vertices in the linked set
     Real,         Pointer :: mat(:,:) => Null()  ! Matrix for these vertices
     Integer               :: valI     => 0       ! Integer value associated to this set
     Real                  :: valR     => 0       ! Real value associated to this set
     Type(Linked), Pointer :: last     => Null()  ! Last linked set
     Type(Linked), Pointer :: next     => Null()  ! Next linked set
   End type Linked

   Interface Copy_pointer_ara
     Module procedure Copy_pointer_araI, Copy_pointer_araR, &
                      Copy_pointer_araD, Copy_pointer_araRD
   End Interface Copy_pointer_ara

   Interface Copy_pointer_mat
     Module procedure Copy_pointer_matI, Copy_pointer_matR, Copy_pointer_matD, &
                      Copy_pointer_matRD
   End Interface Copy_pointer_mat
   
   Interface Copy_pointer_ar3
     Module procedure Copy_pointer_ar3I, Copy_pointer_ar3R, Copy_pointer_ar3D,  &
                      Copy_pointer_ar3RD
   End Interface Copy_pointer_ar3
   
   Interface Copy_multi_List
     Module procedure Copy_multi_list_sngR, Copy_multi_list_sngD, Copy_multi_list_sngRD, &
                      Copy_multi_list_araR, Copy_multi_list_araD, Copy_multi_list_araRD
                      
   End Interface Copy_multi_list
 
   Interface DeAlloc_Adjacency
     Module procedure DeAlloc_adjacencyR, DeAlloc_adjacencyR_ar, DeAlloc_adjacencyD, DeAlloc_adjacencyD_ar
   End Interface DeAlloc_Adjacency
   
   Interface DeAlloc_multi_list
     Module procedure DeAlloc_multi_listR, DeAlloc_multi_listD
   End Interface DeAlloc_multi_list
   
   Interface DeAlloc_multi_list_ar
     Module procedure DeAlloc_multi_listR_ar, DeAlloc_multi_listD_ar
   End Interface DeAlloc_multi_list_ar
   
   Interface DeAlloc_List_of_multi
     Module procedure DeAlloc_List_of_multiR, DeAlloc_List_of_multiD
   End Interface DeAlloc_List_of_multi
   
Contains

   Subroutine DeAlloc_adjacencyR (A)
     Type(Adjacency), Intent(inout) :: A
     
     If (Associated(A%ls)) DeAllocate(A%ls)
     If (Associated(A%vl)) DeAllocate(A%vl)
     If (Associated(A%wt)) DeAllocate(A%wt)
     If (Associated(A%sr)) DeAllocate(A%sr)
     If (Associated(A%wl)) DeAllocate(A%wl)
   End Subroutine DeAlloc_adjacencyR   

   Subroutine DeAlloc_adjacencyR_ar (A)
     Type(Adjacency), Intent(inout) :: A(:)
     
     Integer :: i
     Do i= 1,Size(A)
       Call DeAlloc_adjacencyR (A(i))
     End do
   End Subroutine DeAlloc_adjacencyR_ar   

   Subroutine DeAlloc_adjacencyD (A)
     Type(AdjacencyD), Intent(inout) :: A
     
     If (Associated(A%ls)) DeAllocate(A%ls)
     If (Associated(A%vl)) DeAllocate(A%vl)
     If (Associated(A%wt)) DeAllocate(A%wt)
     If (Associated(A%sr)) DeAllocate(A%sr)
     If (Associated(A%wl)) DeAllocate(A%wl)
   End Subroutine DeAlloc_adjacencyD   

   Subroutine DeAlloc_adjacencyD_ar (A)
     Type(AdjacencyD), Intent(inout) :: A(:)
     Integer :: i

     Do i= 1,Size(A)
       Call DeAlloc_adjacencyD (A(i))
     End do
   End Subroutine DeAlloc_adjacencyD_ar   

   Subroutine DeAlloc_set_list (A)
     Type(Set_list), Intent(inout) :: A
     
     If (Associated(A%set)) DeAllocate(A%set)
     If (Associated(A%lev)) DeAllocate(A%lev)
     If (Associated(A%val)) DeAllocate(A%val)
     If (Associated(A%mbr)) DeAllocate(A%mbr)
   End Subroutine DeAlloc_set_list   

   Subroutine DeAlloc_set_list_ar (A)
     Type(Set_list), Intent(inout) :: A(:)
     Integer :: i

     Do i= 1,Size(A)
       Call DeAlloc_set_list (A(i))
     End do
   End Subroutine DeAlloc_set_list_ar   


   Subroutine DeAlloc_multi_listR (A)
     Type(Multi_listR), Intent(inout) :: A
     
     If (Associated(A%ls)) DeAllocate(A%ls)
     If (Associated(A%lt)) DeAllocate(A%lt)
     If (Associated(A%vl)) DeAllocate(A%vl)
     If (Associated(A%wt)) DeAllocate(A%wt)
     If (Associated(A%px)) DeAllocate(A%px)
     If (Associated(A%qx)) DeAllocate(A%qx)
     If (Associated(A%rx)) DeAllocate(A%rx)
     If (Associated(A%sx)) DeAllocate(A%sx)
     If (Associated(A%tx)) DeAllocate(A%tx)
     If (Associated(A%ux)) DeAllocate(A%ux)
     
     If (Associated(A%L0)) DeAllocate(A%L0)
     If (Associated(A%L1)) DeAllocate(A%L1)
     If (Associated(A%L2)) DeAllocate(A%L2)
     If (Associated(A%L3)) DeAllocate(A%L3)
     If (Associated(A%M0)) DeAllocate(A%M0)
     If (Associated(A%M1)) DeAllocate(A%M1)
     If (Associated(A%M2)) DeAllocate(A%M2)
     If (Associated(A%M3)) DeAllocate(A%M3)
     
     If (Associated(A%Q0)) DeAllocate(A%Q0)
     If (Associated(A%Q1)) DeAllocate(A%Q1)
     If (Associated(A%T0)) DeAllocate(A%T0)
     If (Associated(A%T1)) DeAllocate(A%T1)
   End Subroutine DeAlloc_multi_listR   

   Subroutine DeAlloc_multi_listR_ar (A)
     Type(Multi_listR), Intent(inout) :: A(:)
     Integer :: i

     Do i= 1,Size(A)
       Call DeAlloc_multi_listR (A(i))
     End do
   End Subroutine DeAlloc_multi_listR_ar   

   Subroutine DeAlloc_multi_listD (A)
     Type(Multi_listD), Intent(inout) :: A
     
     If (Associated(A%ls)) DeAllocate(A%ls)
     If (Associated(A%lt)) DeAllocate(A%lt)
     If (Associated(A%vl)) DeAllocate(A%vl)
     If (Associated(A%wt)) DeAllocate(A%wt)
     If (Associated(A%px)) DeAllocate(A%px)
     If (Associated(A%qx)) DeAllocate(A%qx)
     If (Associated(A%rx)) DeAllocate(A%rx)
     If (Associated(A%sx)) DeAllocate(A%sx)
     If (Associated(A%tx)) DeAllocate(A%tx)
     If (Associated(A%ux)) DeAllocate(A%ux)
     
     If (Associated(A%L0)) DeAllocate(A%L0)
     If (Associated(A%L1)) DeAllocate(A%L1)
     If (Associated(A%L2)) DeAllocate(A%L2)
     If (Associated(A%L3)) DeAllocate(A%L3)
     If (Associated(A%M0)) DeAllocate(A%M0)
     If (Associated(A%M1)) DeAllocate(A%M1)
     If (Associated(A%M2)) DeAllocate(A%M2)
     If (Associated(A%M3)) DeAllocate(A%M3)
     
     If (Associated(A%Q0)) DeAllocate(A%Q0)
     If (Associated(A%Q1)) DeAllocate(A%Q1)
     If (Associated(A%T0)) DeAllocate(A%T0)
     If (Associated(A%T1)) DeAllocate(A%T1)
   End Subroutine DeAlloc_multi_listD   

   Subroutine DeAlloc_multi_listD_ar (A)
     Type(Multi_listD), Intent(inout) :: A(:)
     Integer :: i
     Do i= 1,Size(A)
       Call DeAlloc_multi_listD (A(i))
     End do
   End Subroutine DeAlloc_multi_listD_ar   

   Subroutine DeAlloc_List_of_multiR (A)
     Type(List_of_MultiR), Intent(inout) :: A
     Integer :: i, i0, i1
     
     If (Associated(A%ml)) then
       i0= Lbound(A%ml,1);  i1= Ubound(A%ml,1)
       Do i= 1,i0,i1
         Call DeAlloc_multi_listR (A%ml(i))
       End do
       DeAllocate(A%ml)
     End if
   End Subroutine DeAlloc_List_of_multiR   

   Subroutine DeAlloc_List_of_multiD (A)
     Type(List_of_MultiD), Intent(inout) :: A
     Integer :: i, i0, i1
     
     If (Associated(A%ml)) then
       i0= Lbound(A%ml,1);  i1= Ubound(A%ml,1)
       Do i= 1,i0,i1
         Call DeAlloc_multi_listD (A%ml(i))
       End do
       DeAllocate(A%ml)
     End if
   End Subroutine DeAlloc_List_of_multiD   

   
   Subroutine Copy_adjacency (lst, A_in, A)
   
!    Copy one Adjacency structure to another
   
     Integer,            Intent(in) :: lst(:)   ! (n) List of indices of A_in to be copied from
     Type(Adjacency),    Intent(in) :: A_in(:)  ! (m) n <= m, Structure to be copied from
     Type(Adjacency), Intent(inout) :: A(:)     ! (n) Structure to be copied to
!  Local:
     Integer :: i, j, n, m, sv

     m= Size(A_in);  n= Size(lst)
     If (n > m .or. n > Size(A)) then
       Call Out ("Error in Copy_adjacency with Size(lst)",n, "but Size(A_in)",m)
       Call Out ("and Size(A)",Size(A));  Stop
     End if
     
     If (Any (lst < 1 .or. lst > m)) then
       Call Out ("Error in Copy_adjacency. lst",lst);  Stop
     End if
     
     Do j= 1,n
       i= lst(j)
       A(j)%nl    = A_in(i)%nl
       A(j)%pos   = A_in(i)%pos
       A(j)%sum_vl= A_in(i)%sum_vl
       A(j)%sum_wt= A_in(i)%sum_wt

       Call Copy_pointer_ara (A_in(i)%ls, A(j)%ls, sv)
         Call Test_copy (sv,i, 'ls', 'Copy_adjacency')
       Call Copy_pointer_ara (A_in(i)%vl, A(j)%vl, sv)
         Call Test_copy (sv,i, 'vl', 'Copy_adjacency')
       Call Copy_pointer_ara (A_in(i)%wt, A(j)%wt, sv)
         Call Test_copy (sv,i, 'wt', 'Copy_adjacency')
       Call Copy_pointer_ara (A_in(i)%sr, A(j)%sr, sv)
         Call Test_copy (sv,i, 'sr', 'Copy_adjacency')
       Call Copy_pointer_ara (A_in(i)%wl, A(j)%wl, sv)
         Call Test_copy (sv,i, 'sr', 'Copy_adjacency')
     End do
   End Subroutine Copy_adjacency

   Subroutine Copy1_set (A_in, A, n)
   
!    This copies one Set_list structure to another
   
     Type(Set_list),    Intent(in) :: A_in  ! Copy from
     Type(Set_list), Intent(inout) :: A     ! Copy to
     Integer,           Intent(in) :: n
!  Local:
     Integer :: sv
     
     A%n  = A_in%n
     A%p  = A_in%p
     A%r  = A_in%r
     A%svl= A_in%svl
     A%smb= A_in%smb

     Call Copy1_pointer_araI (A_in%set, A%set, sv, n)
       Call Test_copy (sv,0, 'set', 'Copy1_set')
     Call Copy1_pointer_araI (A_in%lev, A%lev, sv, n)
       Call Test_copy (sv,0, 'set', 'Copy1_set')
     Call Copy1_pointer_araR (A_in%val, A%val, sv, n)
       Call Test_copy (sv,0, 'val', 'Copy1_set')
     Call Copy1_pointer_araR (A_in%mbr, A%mbr, sv, n)
       Call Test_copy (sv,0, 'mbr', 'Copy1_set')
   End Subroutine Copy1_set
   
   Pure Subroutine Copy1_pointer_araI (ara_in, ara_out, sv, m)
     Integer,     Pointer :: ara_in(:)   ! Input array
     Integer,     Pointer :: ara_out(:)  ! Copy onto ara_out, starting from 
                             !  its lower bound if associated, else onto the 
                             !  array bounds of ara_in unless 'm' is present, 
                             !  in which case start at 1 
     Integer, Intent(out) :: sv ! = 0 if successful
     Integer,  Intent(in) :: m  ! Copy this # elements of ara_in if > 0 and 
                                ! up to Size(ara_in), or all elements 
                                ! if m <= 0 or m >= Size(ara_in)
!  Local:
     Integer :: n, l, u, lo, uo, szi, szo

     sv= 0
     If (.not.Associated(ara_in)) then
       If (Associated(ara_out)) DeAllocate(ara_out, stat=sv)
       Return
     End if
     
     l= Lbound(ara_in,1);  u= Ubound(ara_in,1);  szi= u - l + 1
     
     If (m > 0) then  
       n= Min(m,szi)
     Else
       n= szi
     End if
     u= l + n - 1
       
     If (Associated(ara_out)) then
       szo= Size(ara_out);  lo= Lbound(ara_out,1);  uo= lo + n - 1
         
       If (szo < n) then
         DeAllocate(ara_out, stat=sv);  Allocate(ara_out(lo:uo))
         ara_out= ara_in(l:u)
       Else
         ara_out(lo:uo)= ara_in(l:u);  ara_out(uo+1:)= 0
       End if
     Else
       Allocate(ara_out(n));  ara_out= ara_in(l:u)
     End if
   End Subroutine Copy1_pointer_araI

   Pure Subroutine Copy1_pointer_araR (ara_in, ara_out, sv, m)
     Real,        Pointer :: ara_in(:)   ! Input array
     Real,        Pointer :: ara_out(:)  ! Copy onto ara_out, starting from its lower bound if associated,
                                !   else onto the array bounds of ara_in unless 'm' is present,
                                !   in which case start at 1
     Integer, Intent(out) :: sv ! = 0 if successful
     Integer,  Intent(in) :: m  ! Copy this # elements of ara_in if > 0 and up to Size(ara_in),
                                !   or all elements if m <= 0 or m >= Size(ara_in)
!  Local:
     Integer :: n, l, u, lo, uo, szi, szo

     sv= 0
     If (.not.Associated(ara_in)) then
       If (Associated(ara_out)) DeAllocate(ara_out, stat=sv)
       Return
     End if
     
     l= Lbound(ara_in,1);  u= Ubound(ara_in,1);  szi= u - l + 1
     
     If (m > 0) then  
       n= Min(m,szi)
     Else
       n= szi
     End if
     u= l + n - 1
       
     If (Associated(ara_out)) then
       szo= Size(ara_out);  lo= Lbound(ara_out,1);  uo= lo + n - 1
         
       If (szo < n) then
         DeAllocate(ara_out, stat=sv);  Allocate(ara_out(lo:uo))
         ara_out= ara_in(l:u)
       Else
         ara_out(lo:uo)= ara_in(l:u);  ara_out(uo+1:)= 0
       End if
     Else
       Allocate(ara_out(n));  ara_out= ara_in(l:u)
     End if
   End Subroutine Copy1_pointer_araR


   Subroutine Copy_set (A_in, A)
   
!    This copy of one Set_list structure to another
   
     Type(Set_list),    Intent(in) :: A_in  ! Structure to be copied from
     Type(Set_list), Intent(inout) :: A     ! Structure to be copied to
!  Local:
     Integer :: sv
     
     A%n  = A_in%n
     A%p  = A_in%p
     A%r  = A_in%r
     A%svl= A_in%svl
     A%smb= A_in%smb

     Call Copy_pointer_ara (A_in%set, A%set, sv)
       Call Test_copy (sv,0, 'set', 'Copy_set')
     Call Copy_pointer_ara (A_in%lev, A%lev, sv)
       Call Test_copy (sv,0, 'set', 'Copy_set')
     Call Copy_pointer_ara (A_in%val, A%val, sv)
       Call Test_copy (sv,0, 'val', 'Copy_set')
     Call Copy_pointer_ara (A_in%mbr, A%mbr, sv)
       Call Test_copy (sv,0, 'mbr', 'Copy_set')
   End Subroutine Copy_set
   
   
   Subroutine Copy_Set_list (lst, A_in, A)
   
!    Copy one Set_list structure to another: A= A_in(lst)
   
     Integer,           Intent(in) :: lst(:)   ! (n) List of indices of A_in to be copied from
     Type(Set_list),    Intent(in) :: A_in(:)  ! (m) m >= n, Structure to be copied from
     Type(Set_list), Intent(inout) :: A(:)     ! (n) Structure to be copied to
!  Local:
     Integer :: i, j, n, m, sv

     m= Size(A_in);  n= Size(lst)
     If (n > m .or. n > Size(A)) then
       Call Out ("Error in Copy_Set_list with Size(lst)",n, "but Size(A_in)",m)
       Call Out ("and Size(A)",Size(A));  Stop
     End if
     
     If (Any (lst < 1 .or. lst > m)) then
       Call Out ("Error in Copy_Set_list. lst",lst);  Stop
     End if
     
     Do j= 1,n
       i= lst(j)
       A(j)%n  = A_in(i)%n
       A(j)%p  = A_in(i)%p
       A(j)%r  = A_in(i)%r
       A(j)%svl= A_in(i)%svl
       A(j)%smb= A_in(i)%smb

       Call Copy_pointer_ara (A_in(i)%set, A(j)%set, sv)
         Call Test_copy (sv,i, 'set', 'Copy_Set_list')
       Call Copy_pointer_ara (A_in(i)%lev, A(j)%lev, sv)
         Call Test_copy (sv,i, 'set', 'Copy_Set_list')
       Call Copy_pointer_ara (A_in(i)%val, A(j)%val, sv)
         Call Test_copy (sv,i, 'val', 'Copy_Set_list')
       Call Copy_pointer_ara (A_in(i)%mbr, A(j)%mbr, sv)
         Call Test_copy (sv,i, 'mbr', 'Copy_Set_list')
     End do
   End Subroutine Copy_Set_list

   Subroutine Copy_multi_list_sngR (Mls_in, Mls_out)
   
     Type(Multi_listR),  Intent(in) :: Mls_in  ! 
     Type(Multi_listR), Intent(out) :: Mls_out ! 
!  Local:
     Integer :: i, j, ni, nl, no, sv

     Mls_out%k     = Mls_in%k
     Mls_out%l     = Mls_in%l
     Mls_out%m     = Mls_in%m
     Mls_out%n     = Mls_in%n
     Mls_out%o     = Mls_in%o
     Mls_out%p     = Mls_in%p
     Mls_out%q     = Mls_in%q
     Mls_out%sum_vl= Mls_in%sum_vl
     Mls_out%sum_wt= Mls_in%sum_wt
     Mls_out%fsx   = Mls_in%fsx
     Mls_out%fux   = Mls_in%fux
     sv= 0

     Call Copy_pointer_ara (Mls_in%ls, Mls_out%ls, sv)
       Call Test_copy (sv,-1, 'ls', 'Copy_multi_list_sngR')
     Call Copy_pointer_ara (Mls_in%lt, Mls_out%lt, sv)
       Call Test_copy (sv,-1, 'lt', 'Copy_multi_list_sngR')
     Call Copy_pointer_ara (Mls_in%vl, Mls_out%vl, sv)
       Call Test_copy (sv,-1, 'vl', 'Copy_multi_list_sngR')
       
     Call Copy_pointer_ara (Mls_in%wt, Mls_out%wt, sv)
       Call Test_copy (sv,-1, 'wt', 'Copy_multi_list_sngR')
     Call Copy_pointer_ara (Mls_in%px, Mls_out%px, sv)
       Call Test_copy (sv,-1, 'px', 'Copy_multi_list_sngR')
     Call Copy_pointer_ara (Mls_in%qx, Mls_out%qx, sv)
       Call Test_copy (sv,-1, 'qx', 'Copy_multi_list_sngR')
     Call Copy_pointer_ara (Mls_in%rx, Mls_out%rx, sv)
       Call Test_copy (sv,-1, 'rx', 'Copy_multi_list_sngR')
     Call Copy_pointer_ara (Mls_in%sx, Mls_out%sx, sv)
       Call Test_copy (sv,-1, 'sx', 'Copy_multi_list_sngR')
     Call Copy_pointer_ara (Mls_in%tx, Mls_out%tx, sv)
       Call Test_copy (sv,-1, 'tx', 'Copy_multi_list_sngR')
     Call Copy_pointer_ara (Mls_in%ux, Mls_out%ux, sv)
       Call Test_copy (sv,-1, 'ux', 'Copy_multi_list_sngR')

    Call Copy_pointer_mat (Mls_in%L0, Mls_out%L0, sv)         
       Call Test_copy (sv,i, 'L0', 'Copy_multi_list_sngR')
     Call Copy_pointer_mat (Mls_in%L1, Mls_out%L1, sv)         
       Call Test_copy (sv,i, 'L1', 'Copy_multi_list_sngR')
     Call Copy_pointer_mat (Mls_in%L2, Mls_out%L2, sv)         
       Call Test_copy (sv,i, 'L2', 'Copy_multi_list_sngR')
     Call Copy_pointer_mat (Mls_in%L3, Mls_out%L3, sv)         
       Call Test_copy (sv,i, 'L3', 'Copy_multi_list_sngR')
     Call Copy_pointer_mat (Mls_in%M0, Mls_out%M0, sv)         
       Call Test_copy (sv,i, 'M0', 'Copy_multi_list_sngR')
     Call Copy_pointer_mat (Mls_in%M1, Mls_out%M1, sv)         
       Call Test_copy (sv,i, 'M1', 'Copy_multi_list_sngR')
     Call Copy_pointer_mat (Mls_in%M2, Mls_out%M2, sv)         
       Call Test_copy (sv,i, 'M2', 'Copy_multi_list_sngR')
     Call Copy_pointer_mat (Mls_in%M3, Mls_out%M3, sv)         
       Call Test_copy (sv,i, 'M3', 'Copy_multi_list_sngR')
       
     Call Copy_pointer_ar3 (Mls_in%Q0, Mls_out%Q0, sv)         
       Call Test_copy (sv,i, 'Q0', 'Copy_multi_list_sngR')
     Call Copy_pointer_ar3 (Mls_in%Q1, Mls_out%Q1, sv)         
       Call Test_copy (sv,i, 'Q1', 'Copy_multi_list_sngR')
     Call Copy_pointer_ar3 (Mls_in%T0, Mls_out%T0, sv)         
       Call Test_copy (sv,i, 'T0', 'Copy_multi_list_sngR')
     Call Copy_pointer_ar3 (Mls_in%T1, Mls_out%T1, sv)         
       Call Test_copy (sv,i, 'T1', 'Copy_multi_list_sngR')
   End Subroutine Copy_multi_list_sngR
     
   Subroutine Copy_multi_list_sngD (Mls_in, Mls_out)
   
     Type(Multi_listD),  Intent(in) :: Mls_in  ! 
     Type(Multi_listD), Intent(out) :: Mls_out ! 
!  Local:
     Integer :: i, j, ni, nl, no, sv

     Mls_out%k     = Mls_in%k
     Mls_out%l     = Mls_in%l
     Mls_out%m     = Mls_in%m
     Mls_out%n     = Mls_in%n
     Mls_out%o     = Mls_in%o
     Mls_out%p     = Mls_in%p
     Mls_out%q     = Mls_in%q
     Mls_out%sum_vl= Mls_in%sum_vl
     Mls_out%sum_wt= Mls_in%sum_wt
     Mls_out%fsx   = Mls_in%fsx
     Mls_out%fux   = Mls_in%fux
     sv= 0

     Call Copy_pointer_ara (Mls_in%ls, Mls_out%ls, sv)
       Call Test_copy (sv,-1, 'ls', 'Copy_multi_list_sngD')
     Call Copy_pointer_ara (Mls_in%lt, Mls_out%lt, sv)
       Call Test_copy (sv,-1, 'lt', 'Copy_multi_list_sngD')
     Call Copy_pointer_ara (Mls_in%vl, Mls_out%vl, sv)
       Call Test_copy (sv,-1, 'vl', 'Copy_multi_list_sngD')
       
     Call Copy_pointer_ara (Mls_in%wt, Mls_out%wt, sv)
       Call Test_copy (sv,-1, 'wt', 'Copy_multi_list_sngD')
     Call Copy_pointer_ara (Mls_in%px, Mls_out%px, sv)
       Call Test_copy (sv,-1, 'px', 'Copy_multi_list_sngD')
     Call Copy_pointer_ara (Mls_in%qx, Mls_out%qx, sv)
       Call Test_copy (sv,-1, 'qx', 'Copy_multi_list_sngD')
     Call Copy_pointer_ara (Mls_in%rx, Mls_out%rx, sv)
       Call Test_copy (sv,-1, 'rx', 'Copy_multi_list_sngD')
     Call Copy_pointer_ara (Mls_in%sx, Mls_out%sx, sv)
       Call Test_copy (sv,-1, 'sx', 'Copy_multi_list_sngD')
     Call Copy_pointer_ara (Mls_in%tx, Mls_out%tx, sv)
       Call Test_copy (sv,-1, 'tx', 'Copy_multi_list_sngD')
     Call Copy_pointer_ara (Mls_in%ux, Mls_out%ux, sv)
       Call Test_copy (sv,-1, 'ux', 'Copy_multi_list_sngD')

    Call Copy_pointer_mat (Mls_in%L0, Mls_out%L0, sv)         
       Call Test_copy (sv,i, 'L0', 'Copy_multi_list_sngD')
     Call Copy_pointer_mat (Mls_in%L1, Mls_out%L1, sv)         
       Call Test_copy (sv,i, 'L1', 'Copy_multi_list_sngD')
     Call Copy_pointer_mat (Mls_in%L2, Mls_out%L2, sv)         
       Call Test_copy (sv,i, 'L2', 'Copy_multi_list_sngD')
     Call Copy_pointer_mat (Mls_in%L3, Mls_out%L3, sv)         
       Call Test_copy (sv,i, 'L3', 'Copy_multi_list_sngD')
     Call Copy_pointer_mat (Mls_in%M0, Mls_out%M0, sv)         
       Call Test_copy (sv,i, 'M0', 'Copy_multi_list_sngD')
     Call Copy_pointer_mat (Mls_in%M1, Mls_out%M1, sv)         
       Call Test_copy (sv,i, 'M1', 'Copy_multi_list_sngD')
     Call Copy_pointer_mat (Mls_in%M2, Mls_out%M2, sv)         
       Call Test_copy (sv,i, 'M2', 'Copy_multi_list_sngD')
     Call Copy_pointer_mat (Mls_in%M3, Mls_out%M3, sv)         
       Call Test_copy (sv,i, 'M3', 'Copy_multi_list_sngD')
       
     Call Copy_pointer_ar3 (Mls_in%Q0, Mls_out%Q0, sv)         
       Call Test_copy (sv,i, 'Q0', 'Copy_multi_list_sngD')
     Call Copy_pointer_ar3 (Mls_in%Q1, Mls_out%Q1, sv)         
       Call Test_copy (sv,i, 'Q1', 'Copy_multi_list_sngD')
     Call Copy_pointer_ar3 (Mls_in%T0, Mls_out%T0, sv)         
       Call Test_copy (sv,i, 'T0', 'Copy_multi_list_sngD')
     Call Copy_pointer_ar3 (Mls_in%T1, Mls_out%T1, sv)         
       Call Test_copy (sv,i, 'T1', 'Copy_multi_list_sngD')
   End Subroutine Copy_multi_list_sngD
   
   Subroutine Copy_multi_list_sngRD (Mls_in, Mls_out)
   
     Type(Multi_listR),  Intent(in) :: Mls_in  ! 
     Type(Multi_listD), Intent(out) :: Mls_out ! 
!  Local:
     Integer :: i, j, ni, nl, no, sv

     Mls_out%k     = Mls_in%k
     Mls_out%l     = Mls_in%l
     Mls_out%m     = Mls_in%m
     Mls_out%n     = Mls_in%n
     Mls_out%o     = Mls_in%o
     Mls_out%p     = Mls_in%p
     Mls_out%q     = Mls_in%q
     Mls_out%sum_vl= Mls_in%sum_vl
     Mls_out%sum_wt= Mls_in%sum_wt
     Mls_out%fsx   = Mls_in%fsx
     Mls_out%fux   = Mls_in%fux
     sv= 0

     Call Copy_pointer_ara (Mls_in%ls, Mls_out%ls, sv)
       Call Test_copy (sv,-1, 'ls', 'Copy_multi_list_sngRD')
     Call Copy_pointer_ara (Mls_in%lt, Mls_out%lt, sv)
       Call Test_copy (sv,-1, 'lt', 'Copy_multi_list_sngRD')
     Call Copy_pointer_ara (Mls_in%vl, Mls_out%vl, sv)
       Call Test_copy (sv,-1, 'vl', 'Copy_multi_list_sngRD')
       
     Call Copy_pointer_ara (Mls_in%wt, Mls_out%wt, sv)
       Call Test_copy (sv,-1, 'wt', 'Copy_multi_list_sngRD')
     Call Copy_pointer_ara (Mls_in%px, Mls_out%px, sv)
       Call Test_copy (sv,-1, 'px', 'Copy_multi_list_sngRD')
     Call Copy_pointer_ara (Mls_in%qx, Mls_out%qx, sv)
       Call Test_copy (sv,-1, 'qx', 'Copy_multi_list_sngRD')
     Call Copy_pointer_ara (Mls_in%rx, Mls_out%rx, sv)
       Call Test_copy (sv,-1, 'rx', 'Copy_multi_list_sngRD')
     Call Copy_pointer_ara (Mls_in%sx, Mls_out%sx, sv)
       Call Test_copy (sv,-1, 'sx', 'Copy_multi_list_sngRD')
     Call Copy_pointer_ara (Mls_in%tx, Mls_out%tx, sv)
       Call Test_copy (sv,-1, 'tx', 'Copy_multi_list_sngRD')
     Call Copy_pointer_ara (Mls_in%ux, Mls_out%ux, sv)
       Call Test_copy (sv,-1, 'ux', 'Copy_multi_list_sngRD')

    Call Copy_pointer_mat (Mls_in%L0, Mls_out%L0, sv)         
       Call Test_copy (sv,i, 'L0', 'Copy_multi_list_sngRD')
     Call Copy_pointer_mat (Mls_in%L1, Mls_out%L1, sv)         
       Call Test_copy (sv,i, 'L1', 'Copy_multi_list_sngRD')
     Call Copy_pointer_mat (Mls_in%L2, Mls_out%L2, sv)         
       Call Test_copy (sv,i, 'L2', 'Copy_multi_list_sngRD')
     Call Copy_pointer_mat (Mls_in%L3, Mls_out%L3, sv)         
       Call Test_copy (sv,i, 'L3', 'Copy_multi_list_sngRD')
     Call Copy_pointer_mat (Mls_in%M0, Mls_out%M0, sv)         
       Call Test_copy (sv,i, 'M0', 'Copy_multi_list_sngRD')
     Call Copy_pointer_mat (Mls_in%M1, Mls_out%M1, sv)         
       Call Test_copy (sv,i, 'M1', 'Copy_multi_list_sngRD')
     Call Copy_pointer_mat (Mls_in%M2, Mls_out%M2, sv)         
       Call Test_copy (sv,i, 'M2', 'Copy_multi_list_sngRD')
     Call Copy_pointer_mat (Mls_in%M3, Mls_out%M3, sv)         
       Call Test_copy (sv,i, 'M3', 'Copy_multi_list_sngRD')
       
     Call Copy_pointer_ar3 (Mls_in%Q0, Mls_out%Q0, sv)         
       Call Test_copy (sv,i, 'Q0', 'Copy_multi_list_sngRD')
     Call Copy_pointer_ar3 (Mls_in%Q1, Mls_out%Q1, sv)         
       Call Test_copy (sv,i, 'Q1', 'Copy_multi_list_sngRD')
     Call Copy_pointer_ar3 (Mls_in%T0, Mls_out%T0, sv)         
       Call Test_copy (sv,i, 'T0', 'Copy_multi_list_sngRD')
     Call Copy_pointer_ar3 (Mls_in%T1, Mls_out%T1, sv)         
       Call Test_copy (sv,i, 'T1', 'Copy_multi_list_sngRD')
   End Subroutine Copy_multi_list_sngRD
   
   
 Subroutine Copy_multi_list_araR (lst, Mls_in, Mls_out)
   
!    This copy of one Multi_listR structure to another, 
!    may permute at the same time: Mls_out= Mls_in(lst)
   
     Integer,              Intent(in) :: lst(:)     ! (nl) List of indices of Mls_in to be copied from
     Type(Multi_listR),    Intent(in) :: Mls_in(:)  ! (ni) Structure to be copied from. 
                                                    !      The range of lst must be a subset of 1...ni
     Type(Multi_listR), Intent(inout) :: Mls_out(:) ! (no) Structure to be copied to (first 'nl' entries)
!  Local:
     Integer :: i, j, n, ni, nl, no, sv

     nl= Size(lst);  ni= Size(Mls_in);  no= Size(Mls_out)
     n= Min(nl,ni,no)
     
     If (nl < 1 .or. ni < nl .or. no < nl) then
       Call Out ("Warning in Copy_multi_list_araR with Size(lst)",nl, "but Size(Mls_in)",ni)
       Call Out ("and Size(Mls_out)",no);  Call Out ("with lst",lst)
     End if
     
     If (Any (lst(:n) < 1 .or. lst(:n) > ni)) then
       Call Out ("Error in Copy_multi_list_araR. lst",lst);  Call Out ("with lst",lst)
       Stop
     End if
     sv= 0

     Do j= 1,n
       i= lst(j)
       Mls_out(j)%k     = Mls_in(i)%k
       Mls_out(j)%l     = Mls_in(i)%l
       Mls_out(j)%m     = Mls_in(i)%m
       Mls_out(j)%n     = Mls_in(i)%n
       Mls_out(j)%o     = Mls_in(i)%o
       Mls_out(j)%p     = Mls_in(i)%p
       Mls_out(j)%q     = Mls_in(i)%q
       Mls_out(j)%sum_vl= Mls_in(i)%sum_vl
       Mls_out(j)%sum_wt= Mls_in(i)%sum_wt
       Mls_out(j)%fsx   = Mls_in(i)%fsx
       Mls_out(j)%fux   = Mls_in(i)%fux

       Call Copy_pointer_ara (Mls_in(i)%ls, Mls_out(j)%ls, sv)
         Call Test_copy (sv,i, 'ls', 'Copy_multi_list_araR')
       Call Copy_pointer_ara (Mls_in(i)%lt, Mls_out(j)%lt, sv)
         Call Test_copy (sv,i, 'lt', 'Copy_multi_list_araR')
       Call Copy_pointer_ara (Mls_in(i)%vl, Mls_out(j)%vl, sv)
         Call Test_copy (sv,i, 'vl', 'Copy_multi_list_araR')
       
       Call Copy_pointer_ara (Mls_in(i)%wt, Mls_out(j)%wt, sv)
         Call Test_copy (sv,i, 'wt', 'Copy_multi_list_araR')
       Call Copy_pointer_ara (Mls_in(i)%px, Mls_out(j)%px, sv)
         Call Test_copy (sv,i, 'px', 'Copy_multi_list_araR')
       Call Copy_pointer_ara (Mls_in(i)%qx, Mls_out(j)%qx, sv)
         Call Test_copy (sv,i, 'qx', 'Copy_multi_list_araR')
       Call Copy_pointer_ara (Mls_in(i)%rx, Mls_out(j)%rx, sv)
         Call Test_copy (sv,i, 'rx', 'Copy_multi_list_araR')
       Call Copy_pointer_ara (Mls_in(i)%sx, Mls_out(j)%sx, sv)
         Call Test_copy (sv,i, 'sx', 'Copy_multi_list_araR')
       Call Copy_pointer_ara (Mls_in(i)%tx, Mls_out(j)%tx, sv)
         Call Test_copy (sv,i, 'tx', 'Copy_multi_list_araR')
       Call Copy_pointer_ara (Mls_in(i)%ux, Mls_out(j)%ux, sv)
         Call Test_copy (sv,i, 'ux', 'Copy_multi_list_araR')

       Call Copy_pointer_mat (Mls_in(i)%L0, Mls_out(j)%L0, sv)         
         Call Test_copy (sv,i, 'L0', 'Copy_multi_list_araR')
       Call Copy_pointer_mat (Mls_in(i)%L1, Mls_out(j)%L1, sv)         
         Call Test_copy (sv,i, 'L1', 'Copy_multi_list_araR')
       Call Copy_pointer_mat (Mls_in(i)%L2, Mls_out(j)%L2, sv)         
         Call Test_copy (sv,i, 'L2', 'Copy_multi_list_araR')
       Call Copy_pointer_mat (Mls_in(i)%L3, Mls_out(j)%L3, sv)         
         Call Test_copy (sv,i, 'L3', 'Copy_multi_list_araR')
       Call Copy_pointer_mat (Mls_in(i)%M0, Mls_out(j)%M0, sv)         
         Call Test_copy (sv,i, 'M0', 'Copy_multi_list_araR')
       Call Copy_pointer_mat (Mls_in(i)%M1, Mls_out(j)%M1, sv)         
         Call Test_copy (sv,i, 'M1', 'Copy_multi_list_araR')
       Call Copy_pointer_mat (Mls_in(i)%M2, Mls_out(j)%M2, sv)         
         Call Test_copy (sv,i, 'M2', 'Copy_multi_list_araR')
       Call Copy_pointer_mat (Mls_in(i)%M3, Mls_out(j)%M3, sv)         
         Call Test_copy (sv,i, 'M3', 'Copy_multi_list_araR')
       
       Call Copy_pointer_ar3 (Mls_in(i)%Q0, Mls_out(j)%Q0, sv)         
         Call Test_copy (sv,i, 'Q0', 'Copy_multi_list_araR')
       Call Copy_pointer_ar3 (Mls_in(i)%Q1, Mls_out(j)%Q1, sv)         
         Call Test_copy (sv,i, 'Q1', 'Copy_multi_list_araR')
       Call Copy_pointer_ar3 (Mls_in(i)%T0, Mls_out(j)%T0, sv)         
         Call Test_copy (sv,i, 'T0', 'Copy_multi_list_araR')
       Call Copy_pointer_ar3 (Mls_in(i)%T1, Mls_out(j)%T1, sv)         
         Call Test_copy (sv,i, 'T1', 'Copy_multi_list_araR')
     End do
     
   End Subroutine Copy_multi_list_araR
   
   Subroutine Copy_multi_list_araD (lst, Mls_in, Mls_out)
   
!    This copy of one Multi_listD structure to another, 
!    may permute at the same time: Mls_out= Mls_in(lst)
   
     Integer,              Intent(in) :: lst(:)     ! (nl) List of indices of Mls_in to be copied from
     Type(Multi_listD),    Intent(in) :: Mls_in(:)  ! (ni) Structure to be copied from. 
                                                    !      The range of lst must be a subset of 1...ni
     Type(Multi_listD), Intent(inout) :: Mls_out(:) ! (no) Structure to be copied to (first 'nl' entries)
!  Local:
     Integer :: i, j, n, ni, nl, no, sv

     nl= Size(lst);  ni= Size(Mls_in);  no= Size(Mls_out)
     n= Min(nl,ni,no)
     
     If (nl < 1 .or. ni < nl .or. no < nl) then
       Call Out ("Warning in Copy_multi_list_araD with Size(lst)",nl, "but Size(Mls_in)",ni)
       Call Out ("and Size(Mls_out)",no);  Call Out ("with lst",lst)
     End if
     
     If (Any (lst(:n) < 1.or. lst(:n) > ni)) then
       Call Out ("Error in Copy_multi_list_araD. lst",lst);  Call Out ("with lst",lst)
       Stop
     End if
     sv= 0

     Do j= 1,n
       i= lst(j)
       Mls_out(j)%k     = Mls_in(i)%k
       Mls_out(j)%l     = Mls_in(i)%l
       Mls_out(j)%m     = Mls_in(i)%m
       Mls_out(j)%n     = Mls_in(i)%n
       Mls_out(j)%o     = Mls_in(i)%o
       Mls_out(j)%p     = Mls_in(i)%p
       Mls_out(j)%q     = Mls_in(i)%q
       Mls_out(j)%sum_vl= Mls_in(i)%sum_vl
       Mls_out(j)%sum_wt= Mls_in(i)%sum_wt
       Mls_out(j)%fsx   = Mls_in(i)%fsx
       Mls_out(j)%fux   = Mls_in(i)%fux

       Call Copy_pointer_ara (Mls_in(i)%ls, Mls_out(j)%ls, sv)
         Call Test_copy (sv,i, 'ls', 'Copy_multi_list_araD')
       Call Copy_pointer_ara (Mls_in(i)%lt, Mls_out(j)%lt, sv)
         Call Test_copy (sv,i, 'lt', 'Copy_multi_list_araD')
       Call Copy_pointer_ara (Mls_in(i)%vl, Mls_out(j)%vl, sv)
         Call Test_copy (sv,i, 'vl', 'Copy_multi_list_araD')
       
       Call Copy_pointer_ara (Mls_in(i)%wt, Mls_out(j)%wt, sv)
         Call Test_copy (sv,i, 'wt', 'Copy_multi_list_araD')
       Call Copy_pointer_ara (Mls_in(i)%px, Mls_out(j)%px, sv)
         Call Test_copy (sv,i, 'px', 'Copy_multi_list_araD')
       Call Copy_pointer_ara (Mls_in(i)%qx, Mls_out(j)%qx, sv)
         Call Test_copy (sv,i, 'qx', 'Copy_multi_list_araD')
       Call Copy_pointer_ara (Mls_in(i)%rx, Mls_out(j)%rx, sv)
         Call Test_copy (sv,i, 'rx', 'Copy_multi_list_araD')
       Call Copy_pointer_ara (Mls_in(i)%sx, Mls_out(j)%sx, sv)
         Call Test_copy (sv,i, 'sx', 'Copy_multi_list_araD')
       Call Copy_pointer_ara (Mls_in(i)%tx, Mls_out(j)%tx, sv)
         Call Test_copy (sv,i, 'tx', 'Copy_multi_list_araD')
       Call Copy_pointer_ara (Mls_in(i)%ux, Mls_out(j)%ux, sv)
         Call Test_copy (sv,i, 'ux', 'Copy_multi_list_araD')

       Call Copy_pointer_mat (Mls_in(i)%L0, Mls_out(j)%L0, sv)         
         Call Test_copy (sv,i, 'L0', 'Copy_multi_list_araD')
       Call Copy_pointer_mat (Mls_in(i)%L1, Mls_out(j)%L1, sv)         
         Call Test_copy (sv,i, 'L1', 'Copy_multi_list_araD')
       Call Copy_pointer_mat (Mls_in(i)%L2, Mls_out(j)%L2, sv)         
         Call Test_copy (sv,i, 'L2', 'Copy_multi_list_araD')
       Call Copy_pointer_mat (Mls_in(i)%L3, Mls_out(j)%L3, sv)         
         Call Test_copy (sv,i, 'L3', 'Copy_multi_list_araD')
       Call Copy_pointer_mat (Mls_in(i)%M0, Mls_out(j)%M0, sv)         
         Call Test_copy (sv,i, 'M0', 'Copy_multi_list_araD')
       Call Copy_pointer_mat (Mls_in(i)%M1, Mls_out(j)%M1, sv)         
         Call Test_copy (sv,i, 'M1', 'Copy_multi_list_araD')
       Call Copy_pointer_mat (Mls_in(i)%M2, Mls_out(j)%M2, sv)         
         Call Test_copy (sv,i, 'M2', 'Copy_multi_list_araD')
       Call Copy_pointer_mat (Mls_in(i)%M3, Mls_out(j)%M3, sv)         
         Call Test_copy (sv,i, 'M3', 'Copy_multi_list_araD')
       
       Call Copy_pointer_ar3 (Mls_in(i)%Q0, Mls_out(j)%Q0, sv)         
         Call Test_copy (sv,i, 'Q0', 'Copy_multi_list_araD')
       Call Copy_pointer_ar3 (Mls_in(i)%Q1, Mls_out(j)%Q1, sv)         
         Call Test_copy (sv,i, 'Q1', 'Copy_multi_list_araD')
       Call Copy_pointer_ar3 (Mls_in(i)%T0, Mls_out(j)%T0, sv)         
         Call Test_copy (sv,i, 'T0', 'Copy_multi_list_araD')
       Call Copy_pointer_ar3 (Mls_in(i)%T1, Mls_out(j)%T1, sv)         
         Call Test_copy (sv,i, 'T1', 'Copy_multi_list_araD')
     End do
     
   End Subroutine Copy_multi_list_araD
   
   Subroutine Copy_multi_list_araRD (lst, Mls_in, Mls_out)
   
!    This copy of a Multi_listR structure to a Multi_listD structure 
!    may permute at the same time: Mls_out= Mls_in(lst)
   
     Integer,              Intent(in) :: lst(:)     ! (nl) List of indices of Mls_in to be copied from
     Type(Multi_listR),    Intent(in) :: Mls_in(:)  ! (ni) Structure to be copied from. 
                                                    !      The range of lst must be a subset of 1...ni
     Type(Multi_listD), Intent(inout) :: Mls_out(:) ! (no) Structure to be copied to (first 'nl' entries)
!  Local:
     Integer :: i, j, n, ni, nl, no, sv

     nl= Size(lst);  ni= Size(Mls_in);  no= Size(Mls_out)
     n= Min(nl,ni,no)
     
     If (nl < 1 .or. ni < nl .or. no < nl) then
       Call Out ("Warning in Copy_multi_list_araRD with Size(lst)",nl, "but Size(Mls_in)",ni)
       Call Out ("and Size(Mls_out)",no);  Call Out ("with lst",lst)
     End if
     
     If (Any (lst(:n) < 1.or. lst(:n) > ni)) then
       Call Out ("Error in Copy_multi_list_araRD. lst",lst);  Call Out ("with lst",lst)
       Stop
     End if
     sv= 0

     Do j= 1,n
       i= lst(j)
       Mls_out(j)%k     = Mls_in(i)%k
       Mls_out(j)%l     = Mls_in(i)%l
       Mls_out(j)%m     = Mls_in(i)%m
       Mls_out(j)%n     = Mls_in(i)%n
       Mls_out(j)%o     = Mls_in(i)%o
       Mls_out(j)%p     = Mls_in(i)%p
       Mls_out(j)%q     = Mls_in(i)%q
       Mls_out(j)%sum_vl= Mls_in(i)%sum_vl
       Mls_out(j)%sum_wt= Mls_in(i)%sum_wt
       Mls_out(j)%fsx   = Mls_in(i)%fsx
       Mls_out(j)%fux   = Mls_in(i)%fux

       Call Copy_pointer_ara (Mls_in(i)%ls, Mls_out(j)%ls, sv)
         Call Test_copy (sv,i, 'ls', 'Copy_multi_list_araRD')
       Call Copy_pointer_ara (Mls_in(i)%lt, Mls_out(j)%lt, sv)
         Call Test_copy (sv,i, 'lt', 'Copy_multi_list_araRD')
       Call Copy_pointer_ara (Mls_in(i)%vl, Mls_out(j)%vl, sv)
         Call Test_copy (sv,i, 'vl', 'Copy_multi_list_araRD')
       
       Call Copy_pointer_ara (Mls_in(i)%wt, Mls_out(j)%wt, sv)
         Call Test_copy (sv,i, 'wt', 'Copy_multi_list_araRD')
       Call Copy_pointer_ara (Mls_in(i)%px, Mls_out(j)%px, sv)
         Call Test_copy (sv,i, 'px', 'Copy_multi_list_araRD')
       Call Copy_pointer_ara (Mls_in(i)%qx, Mls_out(j)%qx, sv)
         Call Test_copy (sv,i, 'qx', 'Copy_multi_list_araRD')
       Call Copy_pointer_ara (Mls_in(i)%rx, Mls_out(j)%rx, sv)
         Call Test_copy (sv,i, 'rx', 'Copy_multi_list_araRD')
       Call Copy_pointer_ara (Mls_in(i)%sx, Mls_out(j)%sx, sv)
         Call Test_copy (sv,i, 'sx', 'Copy_multi_list_araRD')
       Call Copy_pointer_ara (Mls_in(i)%tx, Mls_out(j)%tx, sv)
         Call Test_copy (sv,i, 'tx', 'Copy_multi_list_araRD')
       Call Copy_pointer_ara (Mls_in(i)%ux, Mls_out(j)%ux, sv)
         Call Test_copy (sv,i, 'ux', 'Copy_multi_list_araRD')

       Call Copy_pointer_mat (Mls_in(i)%L0, Mls_out(j)%L0, sv)         
         Call Test_copy (sv,i, 'L0', 'Copy_multi_list_araRD')
       Call Copy_pointer_mat (Mls_in(i)%L1, Mls_out(j)%L1, sv)         
         Call Test_copy (sv,i, 'L1', 'Copy_multi_list_araRD')
       Call Copy_pointer_mat (Mls_in(i)%L2, Mls_out(j)%L2, sv)         
         Call Test_copy (sv,i, 'L2', 'Copy_multi_list_araRD')
       Call Copy_pointer_mat (Mls_in(i)%L3, Mls_out(j)%L3, sv)         
         Call Test_copy (sv,i, 'L3', 'Copy_multi_list_araRD')
       Call Copy_pointer_mat (Mls_in(i)%M0, Mls_out(j)%M0, sv)         
         Call Test_copy (sv,i, 'M0', 'Copy_multi_list_araRD')
       Call Copy_pointer_mat (Mls_in(i)%M1, Mls_out(j)%M1, sv)         
         Call Test_copy (sv,i, 'M1', 'Copy_multi_list_araRD')
       Call Copy_pointer_mat (Mls_in(i)%M2, Mls_out(j)%M2, sv)         
         Call Test_copy (sv,i, 'M2', 'Copy_multi_list_araRD')
       Call Copy_pointer_mat (Mls_in(i)%M3, Mls_out(j)%M3, sv)         
         Call Test_copy (sv,i, 'M3', 'Copy_multi_list_araRD')
       
       Call Copy_pointer_ar3 (Mls_in(i)%Q0, Mls_out(j)%Q0, sv)         
         Call Test_copy (sv,i, 'Q0', 'Copy_multi_list_araRD')
       Call Copy_pointer_ar3 (Mls_in(i)%Q1, Mls_out(j)%Q1, sv)         
         Call Test_copy (sv,i, 'Q1', 'Copy_multi_list_araRD')
       Call Copy_pointer_ar3 (Mls_in(i)%T0, Mls_out(j)%T0, sv)         
         Call Test_copy (sv,i, 'T0', 'Copy_multi_list_araRD')
       Call Copy_pointer_ar3 (Mls_in(i)%T1, Mls_out(j)%T1, sv)         
         Call Test_copy (sv,i, 'T1', 'Copy_multi_list_araRD')
     End do
     
   End Subroutine Copy_multi_list_araRD
   

   Subroutine Reduce_Multi_list (lst,  Mls)
   
!    These 'lst' indices must be in increasing order for this routine to copy
!    from a Multi_listR structure to a proper subset of itself.
   
     Integer,              Intent(in) :: lst(:) ! (nl) Proper subset of 1...nio in increasing order
     Type(Multi_listR), Intent(inout) :: Mls(:) ! (nio) Move the structure at position j= lst(i) 
                                                !  to position i if j > i.
!  Local:
     Logical :: Err
     Integer :: i, j, nl, nio, sv

     nl= Size(lst);   nio= Size(Mls);  If (nl == 0) Return
     
     Err= lst(1) < 1 .or. lst(nl) > nio .or. nl > nio
     If (.not.Err) then
       Do j= 1,nl-1
         Err= lst(j+1) <= lst(j);  If (Err) Exit
       End do
     End if
     If (Err) then
       Call Out ("Error in Reduce_Multi_list. lst",lst);  Stop
     End if
     sv= 0

     Do j= 1,nl
       If (lst(j) <= j) Cycle
       
       i= lst(j)
       Mls(j)%k     = Mls(i)%k
       Mls(j)%l     = Mls(i)%l
       Mls(j)%m     = Mls(i)%m
       Mls(j)%n     = Mls(i)%n
       Mls(j)%o     = Mls(i)%o
       Mls(j)%p     = Mls(i)%p
       Mls(j)%q     = Mls(i)%q
       Mls(j)%sum_vl= Mls(i)%sum_vl
       Mls(j)%sum_wt= Mls(i)%sum_wt
       Mls(j)%fsx   = Mls(i)%fsx
       Mls(j)%fux   = Mls(i)%fux

       Call Copy_pointer_ara (Mls(i)%ls, Mls(j)%ls, sv)
         Call Test_copy (sv,i, 'ls', 'Reduce_Multi_list')
       Call Copy_pointer_ara (Mls(i)%lt, Mls(j)%lt, sv)
         Call Test_copy (sv,i, 'lt', 'Reduce_Multi_list')
       Call Copy_pointer_ara (Mls(i)%vl, Mls(j)%vl, sv)
         Call Test_copy (sv,i, 'vl', 'Reduce_Multi_list')
       
       Call Copy_pointer_ara (Mls(i)%wt, Mls(j)%wt, sv)
         Call Test_copy (sv,i, 'wt', 'Reduce_Multi_list')
       Call Copy_pointer_ara (Mls(i)%px, Mls(j)%px, sv)
         Call Test_copy (sv,i, 'px', 'Reduce_Multi_list')
       Call Copy_pointer_ara (Mls(i)%qx, Mls(j)%qx, sv)
         Call Test_copy (sv,i, 'qx', 'Reduce_Multi_list')
       Call Copy_pointer_ara (Mls(i)%rx, Mls(j)%rx, sv)
         Call Test_copy (sv,i, 'rx', 'Reduce_Multi_list')
       Call Copy_pointer_ara (Mls(i)%sx, Mls(j)%sx, sv)
         Call Test_copy (sv,i, 'sx', 'Reduce_Multi_list')
       Call Copy_pointer_ara (Mls(i)%tx, Mls(j)%tx, sv)
         Call Test_copy (sv,i, 'tx', 'Reduce_Multi_list')
       Call Copy_pointer_ara (Mls(i)%ux, Mls(j)%ux, sv)
         Call Test_copy (sv,i, 'ux', 'Reduce_Multi_list')
         
       Call Copy_pointer_mat (Mls(i)%L0, Mls(j)%L0, sv)         
         Call Test_copy (sv,i, 'L0', 'Reduce_Multi_list')
       Call Copy_pointer_mat (Mls(i)%L1, Mls(j)%L1, sv)         
         Call Test_copy (sv,i, 'L1', 'Reduce_Multi_list')
       Call Copy_pointer_mat (Mls(i)%L2, Mls(j)%L2, sv)         
         Call Test_copy (sv,i, 'L2', 'Reduce_Multi_list')
       Call Copy_pointer_mat (Mls(i)%L3, Mls(j)%L3, sv)         
         Call Test_copy (sv,i, 'L3', 'Reduce_Multi_list')
       Call Copy_pointer_mat (Mls(i)%M0, Mls(j)%M0, sv)         
         Call Test_copy (sv,i, 'M0', 'Reduce_Multi_list')
       Call Copy_pointer_mat (Mls(i)%M1, Mls(j)%M1, sv)         
         Call Test_copy (sv,i, 'M1', 'Reduce_Multi_list')
       Call Copy_pointer_mat (Mls(i)%M2, Mls(j)%M2, sv)         
         Call Test_copy (sv,i, 'M2', 'Reduce_Multi_list')
       Call Copy_pointer_mat (Mls(i)%M3, Mls(j)%M3, sv)         
         Call Test_copy (sv,i, 'M3', 'Reduce_Multi_list')
       
       Call Copy_pointer_ar3 (Mls(i)%Q0, Mls(j)%Q0, sv)         
         Call Test_copy (sv,i, 'Q0', 'Reduce_Multi_list')
       Call Copy_pointer_ar3 (Mls(i)%Q1, Mls(j)%Q1, sv)         
         Call Test_copy (sv,i, 'Q1', 'Reduce_Multi_list')
       Call Copy_pointer_ar3 (Mls(i)%T0, Mls(j)%T0, sv)         
         Call Test_copy (sv,i, 'T0', 'Reduce_Multi_list')
       Call Copy_pointer_ar3 (Mls(i)%T1, Mls(j)%T1, sv)         
         Call Test_copy (sv,i, 'T1', 'Reduce_Multi_list')
     End do
     
   End Subroutine Reduce_Multi_list

   
   Subroutine Test_copy (sv,index, array, routine)
     
     Integer,          Intent(in) :: sv
     Integer,          Intent(in) :: index
     Character(len=*), Intent(in) :: array
     Character(len=*), Intent(in) :: routine
     
     If (sv /= 0) then
       Call Out (array//" copy failed in "//routine)
       Call Out ("At input index",index);  Stop
     End if
   End Subroutine Test_copy
   
   Pure Subroutine Copy_pointer_araI (ara_in, ara_out, sv)
!      Make an exact copy
     Integer,     Pointer :: ara_in(:)   ! Input array
     Integer,     Pointer :: ara_out(:)  ! Copy onto ara_out
     Integer, Intent(out) :: sv          ! = 0 if successful
!  Local:
     Integer :: l, u

     sv= 0
     If (.not.Associated(ara_in)) then
       If (Associated(ara_out)) DeAllocate(ara_out,stat=sv)
       Return
     End if
     
     l= Lbound(ara_in,1);  u= Ubound(ara_in,1)
     
     If (Associated(ara_out)) then
       If (l /= Lbound(ara_out,1) .or. &
           u /= Ubound(ara_out,1)) then
         DeAllocate(ara_out, stat=sv) 
         Allocate(ara_out(l:u))
       End if
     Else
       Allocate(ara_out(l:u))
     End if

     ara_out= ara_in
   End Subroutine Copy_pointer_araI

   Pure Subroutine Copy_pointer_araR (ara_in, ara_out, sv)
!      Make an exact copy
     Real,        Pointer :: ara_in(:)   ! Input array
     Real,        Pointer :: ara_out(:)  ! Copy onto ara_out
     Integer, Intent(out) :: sv          ! = 0 if successful
!  Local:
     Integer :: l, u

     sv= 0
     If (.not.Associated(ara_in)) then
       If (Associated(ara_out)) DeAllocate(ara_out,stat=sv)
       Return
     End if
     
     l= Lbound(ara_in,1);  u= Ubound(ara_in,1)
     
     If (Associated(ara_out)) then
       If (l /= Lbound(ara_out,1) .or. &
           u /= Ubound(ara_out,1)) then
         DeAllocate(ara_out, stat=sv) 
         Allocate(ara_out(l:u))
       End if
     Else
       Allocate(ara_out(l:u))
     End if

     ara_out= ara_in
   End Subroutine Copy_pointer_araR


   Pure Subroutine Copy_pointer_araRD (ara_in, ara_out, sv)
!      Make an exact copy
     Real,          Pointer :: ara_in(:)  ! Input array
     Real(Dblp), Pointer :: ara_out(:) ! Copy onto ara_out
     Integer,   Intent(out) :: sv         ! = 0 if successful
!  Local:
     Integer :: l, u

     sv= 0
     If (.not.Associated(ara_in)) then
       If (Associated(ara_out)) DeAllocate(ara_out,stat=sv)
       Return
     End if
     
     l= Lbound(ara_in,1);  u= Ubound(ara_in,1)
     
     If (Associated(ara_out)) then
       If (l /= Lbound(ara_out,1) .or. &
           u /= Ubound(ara_out,1)) then
         DeAllocate(ara_out, stat=sv) 
         Allocate(ara_out(l:u))
       End if
     Else
       Allocate(ara_out(l:u))
     End if

     ara_out= ara_in
   End Subroutine Copy_pointer_araRD

   Pure Subroutine Copy_pointer_araD (ara_in, ara_out, sv)
!      Make an exact copy
     Real(Dblp), Pointer :: ara_in(:)  ! Input array
     Real(Dblp), Pointer :: ara_out(:) ! Copy onto ara_out
     Integer, Intent(out) :: sv           ! = 0 if successful
!  Local:
     Integer :: l, u

     sv= 0
     If (.not.Associated(ara_in)) then
       If (Associated(ara_out)) DeAllocate(ara_out,stat=sv)
       Return
     End if
     
     l= Lbound(ara_in,1);  u= Ubound(ara_in,1)
     
     If (Associated(ara_out)) then
       If (l /= Lbound(ara_out,1) .or. &
           u /= Ubound(ara_out,1)) then
         DeAllocate(ara_out, stat=sv) 
         Allocate(ara_out(l:u))
       End if
     Else
       Allocate(ara_out(l:u))
     End if

     ara_out= ara_in
   End Subroutine Copy_pointer_araD

   Pure Subroutine Copy_pointer_matI (ara_in, ara_out, sv)
     Integer,     Pointer :: ara_in(:,:)
     Integer,     Pointer :: ara_out(:,:)
     Integer, Intent(out) :: sv
!  Local:
     Integer :: l1, u1, l2, u2

     sv= 0
     If (.not.Associated(ara_in)) then
       If (Associated(ara_out)) DeAllocate(ara_out, stat=sv)
       Return
     End if
     
     l1= Lbound(ara_in,1);  u1= Ubound(ara_in,1)
     l2= Lbound(ara_in,2);  u2= Ubound(ara_in,2)
     
     If (Associated(ara_out)) then
       If (l1 /= Lbound(ara_out,1) .or. u1 /= Ubound(ara_out,1) .or. &
           l2 /= Lbound(ara_out,2) .or. u2 /= Ubound(ara_out,2)) then
         DeAllocate(ara_out, stat=sv);  Allocate(ara_out(l1:u1,l2:u2))
       End if
     Else
       Allocate(ara_out(l1:u1,l2:u2))
     End if

     ara_out= ara_in
   End Subroutine Copy_pointer_matI

   Pure Subroutine Copy_pointer_matR (ara_in, ara_out, sv)
     Real,        Pointer :: ara_in(:,:)
     Real,        Pointer :: ara_out(:,:)
     Integer, Intent(out) :: sv
!  Local:
     Integer :: l1, u1, l2, u2

     sv= 0
     If (.not.Associated(ara_in)) then
       If (Associated(ara_out)) DeAllocate(ara_out, stat=sv)
       Return
     End if
     
     l1= Lbound(ara_in,1);  u1= Ubound(ara_in,1)
     l2= Lbound(ara_in,2);  u2= Ubound(ara_in,2)
     
     If (Associated(ara_out)) then
       If (l1 /= Lbound(ara_out,1) .or. u1 /= Ubound(ara_out,1) .or. &
           l2 /= Lbound(ara_out,2) .or. u2 /= Ubound(ara_out,2)) then
         DeAllocate(ara_out, stat=sv);  Allocate(ara_out(l1:u1,l2:u2))
       End if
     Else
       Allocate(ara_out(l1:u1,l2:u2))
     End if

     ara_out= ara_in
   End Subroutine Copy_pointer_matR

   Pure Subroutine Copy_pointer_matD (ara_in, ara_out, sv)
     Real(Dblp),  Pointer :: ara_in(:,:)
     Real(Dblp),  Pointer :: ara_out(:,:)
     Integer, Intent(out) :: sv
!  Local:
     Integer :: l1, u1, l2, u2

     sv= 0
     If (.not.Associated(ara_in)) then
       If (Associated(ara_out)) DeAllocate(ara_out, stat=sv)
       Return
     End if
     
     l1= Lbound(ara_in,1);  u1= Ubound(ara_in,1)
     l2= Lbound(ara_in,2);  u2= Ubound(ara_in,2)
     
     If (Associated(ara_out)) then
       If (l1 /= Lbound(ara_out,1) .or. u1 /= Ubound(ara_out,1) .or. &
           l2 /= Lbound(ara_out,2) .or. u2 /= Ubound(ara_out,2)) then
         DeAllocate(ara_out, stat=sv);  Allocate(ara_out(l1:u1,l2:u2))
       End if
     Else
       Allocate(ara_out(l1:u1,l2:u2))
     End if

     ara_out= ara_in
   End Subroutine Copy_pointer_matD

   Pure Subroutine Copy_pointer_matRD (ara_in, ara_out, sv)
     Real,          Pointer :: ara_in(:,:)
     Real(Dblp), Pointer :: ara_out(:,:)
     Integer, Intent(out) :: sv
!  Local:
     Integer :: l1, u1, l2, u2

     sv= 0
     If (.not.Associated(ara_in)) then
       If (Associated(ara_out)) DeAllocate(ara_out, stat=sv)
       Return
     End if
     
     l1= Lbound(ara_in,1);  u1= Ubound(ara_in,1)
     l2= Lbound(ara_in,2);  u2= Ubound(ara_in,2)
     
     If (Associated(ara_out)) then
       If (l1 /= Lbound(ara_out,1) .or. u1 /= Ubound(ara_out,1) .or. &
           l2 /= Lbound(ara_out,2) .or. u2 /= Ubound(ara_out,2)) then
         DeAllocate(ara_out, stat=sv);  Allocate(ara_out(l1:u1,l2:u2))
       End if
     Else
       Allocate(ara_out(l1:u1,l2:u2))
     End if

     ara_out= ara_in
   End Subroutine Copy_pointer_matRD

   
   Pure Subroutine Copy_pointer_ar3I (ara_in, ara_out, sv)
     Integer,     Pointer :: ara_in(:,:,:)
     Integer,     Pointer :: ara_out(:,:,:)
     Integer, Intent(out) :: sv
!  Local:
     Integer :: l1, u1, l2, u2, l3, u3

     sv= 0
     If (.not.Associated(ara_in)) then
       If (Associated(ara_out)) DeAllocate(ara_out, stat=sv)
       Return
     End if
     
     l1= Lbound(ara_in,1);  u1= Ubound(ara_in,1)
     l2= Lbound(ara_in,2);  u2= Ubound(ara_in,2)
     l3= Lbound(ara_in,3);  u3= Ubound(ara_in,3)
     
     If (Associated(ara_out)) then
       If (l1 /= Lbound(ara_out,1) .or. u1 /= Ubound(ara_out,1) .or. &
           l2 /= Lbound(ara_out,2) .or. u2 /= Ubound(ara_out,2) .or. &
           l3 /= Lbound(ara_out,3) .or. u3 /= Ubound(ara_out,3)) then
         DeAllocate(ara_out, stat=sv);  Allocate(ara_out(l1:u1,l2:u2,l3:u3))
       End if
     Else
       Allocate(ara_out(l1:u1,l2:u2,l3:u3))
     End if

     ara_out= ara_in
   End Subroutine Copy_pointer_ar3I

   Pure Subroutine Copy_pointer_ar3R (ara_in, ara_out, sv)
     Real,        Pointer :: ara_in(:,:,:)
     Real,        Pointer :: ara_out(:,:,:)
     Integer, Intent(out) :: sv
!  Local:
     Integer :: l1, u1, l2, u2, l3, u3

     sv= 0
     If (.not.Associated(ara_in)) then
       If (Associated(ara_out)) DeAllocate(ara_out, stat=sv)
       Return
     End if
     
     l1= Lbound(ara_in,1);  u1= Ubound(ara_in,1)
     l2= Lbound(ara_in,2);  u2= Ubound(ara_in,2)
     l3= Lbound(ara_in,3);  u3= Ubound(ara_in,3)
     
     If (Associated(ara_out)) then
       If (l1 /= Lbound(ara_out,1) .or. u1 /= Ubound(ara_out,1) .or. &
           l2 /= Lbound(ara_out,2) .or. u2 /= Ubound(ara_out,2) .or. &
           l3 /= Lbound(ara_out,3) .or. u3 /= Ubound(ara_out,3)) then
         DeAllocate(ara_out, stat=sv);  Allocate(ara_out(l1:u1,l2:u2,l3:u3))
       End if
     Else
       Allocate(ara_out(l1:u1,l2:u2,l3:u3))
     End if

     ara_out= ara_in
   End Subroutine Copy_pointer_ar3R


   Pure Subroutine Copy_pointer_ar3D (ara_in, ara_out, sv)
     Real(Dblp), Pointer :: ara_in(:,:,:)
     Real(Dblp), Pointer :: ara_out(:,:,:)
     Integer,   Intent(out) :: sv
!  Local:
     Integer :: l1, u1, l2, u2, l3, u3

     sv= 0
     If (.not.Associated(ara_in)) then
       If (Associated(ara_out)) DeAllocate(ara_out, stat=sv)
       Return
     End if
     
     l1= Lbound(ara_in,1);  u1= Ubound(ara_in,1)
     l2= Lbound(ara_in,2);  u2= Ubound(ara_in,2)
     l3= Lbound(ara_in,3);  u3= Ubound(ara_in,3)
     
     If (Associated(ara_out)) then
       If (l1 /= Lbound(ara_out,1) .or. u1 /= Ubound(ara_out,1) .or. &
           l2 /= Lbound(ara_out,2) .or. u2 /= Ubound(ara_out,2) .or. &
           l3 /= Lbound(ara_out,3) .or. u3 /= Ubound(ara_out,3)) then
         DeAllocate(ara_out, stat=sv);  Allocate(ara_out(l1:u1,l2:u2,l3:u3))
       End if
     Else
       Allocate(ara_out(l1:u1,l2:u2,l3:u3))
     End if

     ara_out= ara_in
   End Subroutine Copy_pointer_ar3D

   Pure Subroutine Copy_pointer_ar3RD (ara_in, ara_out, sv)
     Real,          Pointer :: ara_in(:,:,:)
     Real(Dblp), Pointer :: ara_out(:,:,:)
     Integer, Intent(out) :: sv
!  Local:
     Integer :: l1, u1, l2, u2, l3, u3

     sv= 0
     If (.not.Associated(ara_in)) then
       If (Associated(ara_out)) DeAllocate(ara_out, stat=sv)
       Return
     End if
     
     l1= Lbound(ara_in,1);  u1= Ubound(ara_in,1)
     l2= Lbound(ara_in,2);  u2= Ubound(ara_in,2)
     l3= Lbound(ara_in,3);  u3= Ubound(ara_in,3)
     
     If (Associated(ara_out)) then
       If (l1 /= Lbound(ara_out,1) .or. u1 /= Ubound(ara_out,1) .or. &
           l2 /= Lbound(ara_out,2) .or. u2 /= Ubound(ara_out,2) .or. &
           l3 /= Lbound(ara_out,3) .or. u3 /= Ubound(ara_out,3)) then
         DeAllocate(ara_out, stat=sv);  Allocate(ara_out(l1:u1,l2:u2,l3:u3))
       End if
     Else
       Allocate(ara_out(l1:u1,l2:u2,l3:u3))
     End if

     ara_out= ara_in
   End Subroutine Copy_pointer_ar3RD

   
   Pure Subroutine Extend_Lnk (Lnk)
     Implicit None
     Type(Linked), Pointer :: Lnk
     Type(Linked), Pointer :: Tmp
     Integer :: sv

     If (Associated(Lnk)) then
       Allocate (Tmp, Stat=sv);  Lnk%next => Tmp
       Tmp%last => Lnk;  Nullify(Tmp%next);  Lnk => Tmp
     Else
       Allocate(Lnk, Stat=sv);  Nullify(Lnk%last, Lnk%next)
     End if
   End Subroutine Extend_Lnk

   Pure Subroutine Backward_Lnk (n, Lnk) ! Move backward 'n' times
     Implicit None
     Integer,   Intent(in) :: n
     Type(Linked), Pointer :: Lnk
     Integer :: i

     Do i= 1,n
       Lnk => Lnk%last
     End do
   End Subroutine Backward_Lnk

   Pure Subroutine Forward_Lnk (n, Lnk)  ! Move forward 'n' times
     Implicit None
     Integer,   Intent(in) :: n
     Type(Linked), Pointer :: Lnk
     Integer :: i

     Do i= 1,n
       Lnk => Lnk%next
     End do
   End Subroutine Forward_Lnk
   
   Pure Subroutine DeAllocate_LnkPtr (Lnk)
     Implicit None
     Type(Linked), Pointer :: Lnk
     
     If (.not.Associated(Lnk)) Return
     Do
       If (.not.Associated(Lnk%next)) Exit
       Lnk => Lnk%next
     End do

     Do
       If (Associated(Lnk%Comp))    DeAllocate(Lnk%Comp)
       If (Associated(Lnk%vertex)) DeAllocate(Lnk%vertex)
       If (Associated(Lnk%mat))    DeAllocate(Lnk%mat)
       If (.not.Associated(Lnk%last)) Exit
       Lnk => Lnk%last;  Nullify(Lnk%next)
     End do     
   End Subroutine DeAllocate_LnkPtr

End Module Types  

