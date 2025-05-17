 
 Module Graph_algorithms
 
   Use Util
   Use Factorials
   Use Precisn
   Use Types
   Use Output
   Implicit None

!  Diagnostic variables

   Integer :: n_greedy        ! # subroutine calls to Greedy_ind
   Integer :: n_not_so_greedy ! # recursive subroutine calls to Not_so_Greedy_ind

 Contains
   
   Subroutine Best_cliques (skip,cut,Value,Connect, n1,sing,svl, &
                            n2,pair,pvl, n3,Comp3, Best3)

!    List maximal cliques 'Best3' (highest sums of connectivity weighted values %smb)
!    for a connectivity matrix 'Connect', such as an overlap or correlation matrix,
!    using the vertex values 'Value'.

!    Do branch and bound, using %smb for bounding.

     Integer,      Intent(in) :: skip         ! Skip limit for the Not_so_Greedy search
     Real,         Intent(in) :: cut          ! Edge determined if connectivity > cut
     Real,         Intent(in) :: Value(:)     ! (nv) Decreasing vertex values
     Real,         Intent(in) :: Connect(:,:) ! (nv,nv) Connectivities, values <= 1.0,
                                              !         diagonal not used.
     
     Integer,     Intent(out) :: n1           ! # singleton components
     Integer,     Intent(out) :: sing(:)      ! (nv to n1) Singleton components, in order of decreasing value
     Real,        Intent(out) :: svl(:)       ! (nv to n1) Values of singleton components = vertex value
     
     Integer,     Intent(out) :: n2           ! # pair components
     Integer,     Intent(out) :: pair(:,:)    ! (2, nv/2 to n2) Pair components, in order 
                                              !   of decreasing connectivity weighted value
     Real,        Intent(out) :: pvl(:)       ! (nv/2 to n2) The connectivity values of pair components
     
     Integer,     Intent(out) :: n3           ! # graph components of size >= 3
     Type(Set_list),  Pointer :: Comp3(:)     ! (n3) Components of size >= 3
                                              !  %n     = component size
                                              !  %p     = index of the last clique in Best3 
                                              !           from this component
                                              !  %set(n)= list of component vertices
     Type(Set_list),  Pointer :: Best3(:)     ! (mb) Top cliques from all the components of size >= 3
                                              !      ordered by decreasing connectivitiy weighted 
                                              !      value %smb within each component
                                              !  %n      = clique size
                                              !  %p      = component index in Comp3 of each clique 'q'
                                              !  %svl    = clique value = sum of vertex values
                                              !  %smb    = clique connectivity weighted value = sum of 
                                              !            connectivities weighted by vertex values
                                              !  %set(n) = clique vertices
                                              !  %val(n) = clique vertex values
                                              !  %mbr(n) = clique vertex memberships = normalized 
                                              !            mean connectivities
!  Local:
     Type(Adjacency), Allocatable :: Graph(:)
     Type(Set_list),  Allocatable :: Clq(:), Clq_Cp(:)
     Logical, Allocatable :: Adj(:), Edge_Cp(:,:)
     Real,    Allocatable :: Conn0(:,:), val_Cp(:), Conn_Cp(:,:)
     Integer, Allocatable :: ls(:,:), cliq(:), mx_num(:)
     Integer, Allocatable :: ls_cp(:), vr_cp(:), compsize(:), Component(:)
     Logical :: ReOrd
     Real    :: sm, val_sum, max_corr
     Integer :: i, j, k, l, m, n, q, v, cp, nc, nl, nq, np, nv, q1, sv

     n1= 0;  sing= 0;  svl= 0;  n2= 0;  pair= 0;  pvl= 0;  n3= 0
     
     nv= Size(Value)  
     If (nv < 2) then
       If (nv == 1) then
         n1= 1;  sing(1)= 1
       End if  
       Return
     Else
       Allocate(Conn0(nv,nv));  Conn0= Connect
       Call Const_to_diag (-1.0,Conn0)
         max_corr= Maxval(Conn0)
       Call Const_to_diag (1.0,Conn0)

       If (max_corr < cut) then  ! All singletons
         n1= nv;  sing= "ID";  svl= Value
         Return
       End if
     End if
     
     Allocate(Graph(nv), ls_cp(nv), vr_cp(nv), compsize(nv), Component(nv))
     Graph(:)%nl= 0;  ls_cp= 0;  vr_cp= 0;  compsize= 0;  Component= 0
     
!    Compute the connected components

     Call Form_graph (cut, Conn0, Graph)
     
     If (pr_out > 1.5) then
       Call Out (0,"For connectivity matrix", Conn0)
       Call Out ("and cut value",cut)
       Do v= 1,nv
         n= Graph(v)%nl;  If (n < 1) Cycle
         Call Out ("For vertex",v,"# adjacent vertices",n, ln=1)
         Call Out ("Adjacent vertices",Graph(v)%ls)
       End do
     End if
     
!    Compute the # components 'nc', thier sizes compsize(:nc), and 
!    the vertex to component map 'Component':

     Call Connected_components (Graph, nc,compsize, Component, ls_cp)

!    Singletons
     
     Call List_of_true (compsize(:nc) == 1, n1,ls_cp) ! 'ls_cp' lists singletons
     
     If (n1 > 0) then
       Do cp= 1,n1
         sing(cp)= First_true(Component == ls_cp(cp))
         svl(cp) = Value(sing(cp))
       End do
       
       If (n1 > 1) then
         Call Sort (.false., svl(:n1), ls_cp(:n1),ReOrd)
         If (ReOrd) sing(:n1)= sing(ls_cp(:n1))
       End if
       
       If (pr_out > 1.5) then
         Call Out ("Singleton components", sing(:n1))
         Call Out ("Corresponding values", svl(:n1))
       End if
     End if

!    Pair components
     
     Call List_of_true (compsize(:nc) == 2, n2,ls_cp) ! 'ls_cp' lists pair components
     
     If (n2 > 0) then
       Do cp= 1,n2
         Call List_of_true (Component == ls_cp(cp), n,pair(:,cp))
         pvl(cp)= Sum(Value(pair(:,cp))) * Conn0(pair(1,cp),pair(2,cp))
       End do
       
       If (n2 > 1) then
         Call Sort (.false., pvl(:n2), ls_cp(:n2),ReOrd)
         If (ReOrd) pair(:,:n2)= pair(:,ls_cp(:n2))
       End if
       
       If (pr_out > 1.5) then
         Call Out (-1,"Pair components", pair(:,:n2))
         Call Out ("Connectivity value of each pair", pvl(:n2))
       End if
     End if

!    Components of size >= 3

     Call List_of_true (compsize(:nc) >= 3, n3,ls_cp)
     nq= 0

     If (n3 > 0) then
       If (Allocated(mx_num)) DeAllocate(mx_num);  Allocate (mx_num(n3))
       Do cp= 1,n3
         mx_num(cp)= compsize(ls_cp(cp)) ! Size of each component
       End do
       m= Sum(mx_num);  Allocate(Comp3(n3), Clq(m))
       
       If (pr_out >= 1.5) then
         Call Out ("For # components of size >= 3",n3, ln=1)
         Call Out ("Component sizes",compsize(ls_cp(:n3)))
       End if
     
       Component_loop : Do cp= 1,n3
         Call List_of_true (Component == ls_cp(cp), np,vr_cp) ! Vertices of the component

         Allocate(Comp3(cp)%set(np), val_Cp(np), Conn_Cp(np,np), &
                  cliq(np), Adj(np), Edge_Cp(np,np), Clq_Cp(mx_num(cp)))
         
         Comp3(cp)%n= np;  Comp3(cp)%set= vr_cp(:np)
         val_Cp = Value(Comp3(cp)%set)
         Conn_Cp= Conn0(Comp3(cp)%set,Comp3(cp)%set)
         Edge_Cp= Conn_Cp > cut

         If (pr_out >= 1.5) then
           Call Out ("For component",ls_cp(cp), "of size",np, ln=1)
           Call Out ("with vertices",Comp3(cp)%set)
           Call Out (-1,"connectivity matrix",Conn_Cp)
           Call Out (-1,"and edge matrix",Edge_Cp)
           Call Out ("Component vertex values",val_Cp)
         End if
       
         m= Min(skip,np)
         Call Not_so_Greedy_cliq (np,1,m, Edge_Cp,val_Cp,Conn_Cp, &
                                  Adj, -1,cliq, Clq_Cp)
         m= Last_true(Clq_Cp%n > 0) 
      
!        Add the maximal cliques from component 'cp' to 'Clq'

         Do q= 1,m
           nq= nq + 1;  n= Clq_Cp(q)%n;  Clq_Cp(q)%p= cp
           Clq_Cp(q)%set(:n)= vr_cp(Clq_Cp(q)%set(:n))

           Call Copy_set (Clq_Cp(q), Clq(nq))
!           Call Copy1_set (Clq_Cp(q), Clq(nq), n)
         End do
         Comp3(cp)%p= nq
         
         Call DeAlloc_set_list_ar (Clq_cp)
         DeAllocate(val_Cp, Conn_Cp, cliq, Adj, Edge_Cp, Clq_Cp)
       End do Component_loop
     
!      Copy 'Clq' to 'Best3' 
     
       If (nq > 0) then
         Allocate(Best3(nq))
         Do q= 1,nq
           Call Copy_set (Clq(q), Best3(q))
         End do
       End if
       Call DeAlloc_set_list_ar (Clq);  DeAllocate(Clq)
     End if
     
     Call DeAlloc_adjacency (Graph);  DeAllocate(Graph)
                        
   End Subroutine Best_cliques
                            
   Subroutine Form_graph (cut, Connect, Graph)

!    Form the 'Graph' on 'nv' vertices whose edges correspond to 
!    connectivity values 'Connect' >= cut, omitting diagonal values,
!    if cut > 0, or to  'Connect' < -cut, if cut < 0
   
     Real,             Intent(in) :: cut          ! Connectivity cutoff value for determining Graph edges
     Real,             Intent(in) :: Connect(:,:) ! (nv,nv) Connectivity matrix, diagonal not used
     Type(Adjacency), Intent(out) :: Graph(:)     ! (nv) The corresponding connectivity or disconnectivy Graph 
!  Local:
     Logical :: msk(Size(Graph))
     Integer :: i, n, v, nv, sv, lst(Size(Graph))
     
     nv= Size(Graph)
        
     Do v= 1,nv
       If (cut > 0) then
         Forall(i=1:nv) msk(i)= Connect(v,i) >= cut .and. i /= v
       Else
         Forall(i=1:nv) msk(i)= Connect(v,i) < -cut .and. i /= v
       End if
       Call List_of_true (msk, n,lst);  Graph(v)%nl= n
       
       If (n > 0) then
         Graph(v)%pos= First_true(1, lst(:n) > v)
         Allocate (Graph(v)%ls(n), Stat=sv);  Graph(v)%ls= lst(:n)
       Else
         Graph(v)%pos= 1
       End if
     End do
     
   End Subroutine Form_graph
    
   Subroutine Connected_components (Graph, nc, compsize, Component, Alloc_num)

!    "Graph" is decomposed into its connected components. When there is more than
!    one component, determining an optimal independent set for each component separately  
!    can give a big reduction in computational load, especially for branch-n-bound.
	
!    The algorithm to find a connected component uses several passes through an ordered list
!    of the vertices, each time finding vertices in a new "level" of the component. For example, 
!    to find the component of vertex 1, first set all the vertex labels to 0 (unlabeled) 
!    except set the first vertex's label to 1. Next label as level 1 all subsequent vertices 
!    adjacent to the first vertex. Then sweep through all the vertices beyond vertex 1,
!    checking for level 1 vertices. When found label all its subsequent unlabeled adjacent
!    vertices again with level 1 but all prior unlabeled adjacent vertices with level 2.
!    Continue this way through the end of the list of vertices. Then repeat this sweep
!    but now looking for level 2 vertices, so that this second sweep substitutes for 
!    backtracking during the level 1 sweep. 
   
!    During this second sweep label all the subsequent unlabeled adjacent vertices of a
!    level 2 vertex again with level 2 but all prior unlabeled adjacent vertices with level 3.
!    This is followed by similar sweep but for label 3 vertices, or until no new vertices of
!    higher level are found.  Then all the vertices of positive label form a connected component.
	
!    Now set all the labels of the vertices in the first component to -1 and find the first 
!    unlabeled vertex. This will be the start of the same kind of sweeps done for the first
!    component, automatically skipping over all negative labels. At the end change the 
!    component 2 labels to -2 and proceed to component 3. When all components have been found, 
!    negate all the labels so that each vertex is labeled by its (positive) component number.
	
!    After this, sort the component labels by increasing component size, so that singleton 
!    components are first, then pair components, etc., with the largest components last.

     Type(Adjacency),    Intent(in) :: Graph(:)     ! (nv) Graph to be analyzed

     Integer,           Intent(out) :: nc           ! # components
     Integer,           Intent(out) :: compsize(:)  ! (nv) compsize(c) = # vertices in the cth component
     Integer,           Intent(out) :: Component(:) ! (nv) Component(v) = c if vertex v is in the cth component
     Integer, Optional, Intent(out) :: Alloc_num(:) ! (nv) Alloc_num(c) a rough measure of the computational complexity 
                                                    !   of a maximal independent set for component c 
                                                    !   = (# edges in c) * (# edges in the complement graph of c).
                                                    !   Used for allocating CPU to the components
!  Local:
     Integer :: c, i, j, k, v, n1, n2, nv, v0, v1
     Integer :: lev, lev1
     Integer, Allocatable :: key(:), key_I(:)
     
     nv= Size(Graph);  If (nv < 1) Return
     
     v0= 1;  Component= 0;  compsize= 0

     Component_loop : Do nc= 1,nv

!     Determine the initial members of component nc.
!       (Assume that all vertices prior to v0 are already
!        assigned to components, plus vertices subsequent
!        to vertex v0 with negative Component value.)

      If (v0 == nv) then
        Component(v0)= -nc;  compsize(nc)= 1
        Exit Component_loop
      End if

      v1= v0 + 1;  lev= 1;  Component(v0)= lev

      Call ReAssign (0,lev, Graph(v0)%pos, Graph(v0)%ls, Component)

!     Search through the level 'lev' members and add their connections,
!     assigning level lev+1 to prior connections and level lev to subsequent
!     connections.  Halt the level loop when no lev+1 connections
!     are assigned.

      Level_loop : Do
        j= 0;  lev1= lev + 1

        Do v= v1,nv
          If (Component(v) == lev) then
            Call ReAssign (0,lev1, Graph(v)%pos, Graph(v)%ls, Component, j)
            Call ReAssign (0,lev, Graph(v)%pos, Graph(v)%ls, Component)
          End if
        End do

        If (j == 0)  Exit Level_loop

        lev= lev1
      End do Level_loop
      
      If (pr_out > 1.5) then
        Call Out ("For component",nc,"# levels",lev)
        Call Out ("Current component vector",Component)
      End if

!     Record Component "nc" as -nc, count it, and test for completion

      j= 0
      Do v= v0,nv
        If (Component(v) > 0) then
          j= j + 1;  Component(v)= -nc
        End if
      End do

      compsize(nc)= j;  k= v0 - 1
      v0= k + First_true(Component(v0:nv) == 0)
      If (v0 <= k)  Exit Component_loop
    End do Component_loop

    Component= -Component

!   Sort component labels by component size

    If (nc > 1) then
      Allocate (key(nc), key_I(nc))
      Call Sort (.true., compsize(:nc), key)
      Call Inverse_map (key, key_I)
      Component= key_I(Component)
      DeAllocate (key, key_I)
    End if

    If (pr_out > 1.5) then
      Call Out ("# connected components of the graph",nc,ln=1)
      If (nc > 1) then
        Call Out ("Component sizes",  compsize(:nc))
        Call Out ("Final component vector",Component)
      End if
    End if

!   Compute a rough measure of the computational complexity of a maximal independent set
!   for component c = (# edges in c) * (# edges in the complement graph of c)

    If (Present(Alloc_num)) then
      Alloc_num= 0
      If (nc <= 1) then
        Alloc_num(1)= Sum(Graph%nl) / 2.0
        Alloc_num(1)= Alloc_num(1) * (nv*(nv-1)/2 - Alloc_num(1))
        Alloc_num(1)= Max(Alloc_num(1), 1)
      Else
        Do v= 1,nv
          c= Component(v);  Alloc_num(c)= Alloc_num(c) + Graph(v)%nl
        End do
        Alloc_num(:nc)= Alloc_num(:nc) / 2.0
        Alloc_num(:nc)= Alloc_num(:nc) * (compsize(:nc)*(compsize(:nc)-1)/2 - Alloc_num(:nc))
      End if
      If (nc > 1 .and. pr_out > 1.5)  Call Out ("Allocation numbers",  Alloc_num(:nc))
    End if

   End Subroutine Connected_components

   Pure Subroutine ReAssign (val_in,val_out, i, ls, Comp, nre)

!    Re-assign "val_in" values of "Comp" to "val_out" 
!    where they are enumerated by the list "ls" for
!    indices < i if 'nre' is Present or indices >= i otherwise.

     Integer,    Intent(in) :: val_in, val_out
     Integer,    Intent(in) :: i
     Integer,    Intent(in) :: ls(:)
     Integer, Intent(inout) :: Comp(:)
     Integer, Optional, Intent(inout) :: nre
!  Local:
     Integer :: j, k

     If (Present(nre)) then
       Do j= 1,i-1
         k= ls(j)
         If (Comp(k) == val_in) then
           nre= nre + 1;  Comp(k)= val_out
         End if
       End do
     Else
       Do j= i,Size(ls)
         k= ls(j);  If (Comp(k) == val_in) Comp(k)= val_out
       End do
     End if
   End Subroutine ReAssign

    Subroutine Stable_matching (n,m, mxval, val, match)
    
!     Gale - Shapley algorithm for stable matching
!     given edge values (costs) that determine preferences.
    
!     The matching should minimize the costs = cut size.
    
      Integer,  Intent(in) :: n,m       ! Sizes of sets to match. n <= m
      Real,     Intent(in) :: mxval     ! Maximum value = infeasible match
      Real,     Intent(in) :: val(:,:)  ! (n,m)  Values that determine the preferences:
                                        !        lower value = more preferred
      Integer, Intent(out) :: match(:)  !  The matching (1st dimension to 2nd)
                                        !    0 if no matching 
!   Local:
      Logical :: Mat(n,m)
      Real    :: tmp1(n), tmp2(m)
      Integer :: l1(-1:m,n), l2(0:n,m)
      Integer :: i, j, k, r, p
      
      match= 0;  Mat= .false.;  l1= 0;  l2= 0
      
!     Compute l1(1:r,i)= feasible matches for man i in order of preference,
!     with r= l1(0,i), and l2(1:p,j)= feasible matches for woman j
!     in order of preference, with p= l2(0,j).

      Do i= 1,n
        tmp1= val(i,:);  Call Sort (.true., tmp1, l1(1:,i))
        l1(0,i)= Last_true (tmp1 < mxval)
        If (l1(0,i) == 0) Return
      End do
      
      Do j= 1,m
        tmp2= val(:,j);  Call Sort (.true., tmp2, l2(1:,j))
        l2(0,j)= Last_true (tmp2 < mxval)
        If (l2(0,j) == 0) Return
      End do
      
!     Test for a perfect matching already
      
      Do i= 1,n
        r= l1(0,i);  j= l1(1,i);  p= l2(0,j)
        If (r > 0 .and. p > 0 .and. l2(1,j) == i) match(i)= j
      End do
                                              
      If (All(match > 0)) Return  ! Stable match
      
!     Do the proposal and acceptance loop
      
      Do
        
!       Each unassigned man proposes to his next preference woman
        
        Mat= .false.
        Do i= 1,n
          If (match(i) > 0) then
            Mat(i,match(i))= .true.
          Else
            If (l1(-1,i) >= l1(0,i)) Return
            k= l1(-1,i) + 1;  Mat(i,l1(k,i))= .true.;  l1(-1,i)= k
          End if
        End do

!       Each woman accepts her best current proposal
        
        match= 0
        Do j= 1,m
          r= l2(0,j);  k= First_true(Mat(l2(1:r,j),j))
          If (k > 0)  match(l2(k,j))= j
        End do
      
        If (All(match > 0)) Return  ! Stable match
      End do
      
    End Subroutine Stable_matching
    
  Recursive Subroutine Not_so_Greedy_cliq (nv,v1,skip, Edge,Value,Connect, &
                                           Last_Adj, nl,Cur_cliq, Best_cliq)

!   This routine implements the not-so-greedy algorithm to find a nearly
!   optimal maximal clique, utilizing a skip parameter to search 
!   a neighborhood of the greedy solution.

!   Each partially constructed set is a maximal clique on the smallest 
!   initial segment of vertices 1,...,v that includes it. The best solutions 
!   are recorded in 'Best_cliq' according to decreasing sum of vertex values.

    Integer,            Intent(in) :: nv               ! Total # vertices in 'Edge'
    Integer,            Intent(in) :: v1               ! = Cur_cliq(nl) + 1 = next potential vertex 
                                                       !   in the clique. v1 <= nv required
    Integer,            Intent(in) :: skip             ! Skip parameter = "size" of the neighborhood 
                                                       ! of the greedy extension = max # branches to skip.  
                                                       ! 0 => ordinary greedy algorithm.
    Logical,            Intent(in) :: Edge(:,:)        ! (nv,nv) Edge matrix, diagonal not used
    Real,               Intent(in) :: Value(:)         ! (nv) Vertex values
    Real,               Intent(in) :: Connect(:,:)     ! (nv,nv) Connectivity matrix, diagonal not used
    Logical,            Intent(in) :: Last_Adj(v1:)    ! (v1:nv) Vertices adjacent to all members
                                                       !  of the current clique
    Integer,            Intent(in) :: nl               ! # vertices already in Cur_cliq
    
    Integer,         Intent(inout) :: Cur_cliq(:)      ! (nv) Current clique (through 'nl')
    Type(Set_list),  Intent(inout) :: Best_cliq(:)     ! (mb) Records up to 'mb' good independent sets
                                                       !   %n      = clique size
                                                       !   %svl    = clique value (= sum of vertex values)
                                                       !   %smb    = clique connectivity value 
                                                       !           = sum of mean connectivities, 
                                                       !             using vertex value weights
                                                       !   %set(n) = clique vertices
                                                       !   %val(n) = clique vertex values
                                                       !   %mbr(n) = clique vertex memberships 
                                                       !           = mean connectivities, using 
                                                       !             vertex value weights / %smb
! Local:
    Logical :: Cur_Adj(v1:nv)
    Real    :: max_val
    Integer :: i, j, v, n, mb, n1, nc, no, v2, vl, skp

    mb= Size(Best_cliq)
    
!   Initialize

    If (nl < 0) then
      n_not_so_greedy= 1;  n_greedy= 1;  
      Best_cliq%n= 0;  Best_cliq%svl= 0;  Best_cliq%smb= 0
      Do i= 1,mb
        If (Associated(Best_cliq(i)%set)) then
          DeAllocate(Best_cliq(i)%set, Best_cliq(i)%val, Best_cliq(i)%mbr)
        End if
      End do
      
!     Compute the Greedy solution
      
      Cur_Adj= .true.;  Cur_cliq= 0;  nc= 0

      Call Greedy_cliq (1,Edge, Cur_Adj, nc,Cur_cliq)

      Call Test_cliq ("Initial Greedy", Cur_cliq(:nc),Value, &
                      Connect, Best_cliq)

      If (nc >= nv .or. skip <= 0) Return;  Cur_Adj= .true.;  Cur_cliq= 0
      
      Call Not_so_Greedy_cliq (nv,1,skip, Edge,Value,Connect, Cur_Adj, &
                               0,Cur_cliq, Best_cliq)
      Return
    End if
    
!   Do the Bounds and Skip tests
    
    If (nl > 0) then
      max_val= Sum(Value(Cur_cliq(:nl))) + Sum(Value(v1:), mask=Last_Adj)

      If (max_val <= Best_cliq(mb)%svl) Return  ! Prune by bound test
      
      n_not_so_greedy= n_not_so_greedy + 1

!     Finish with a Greedy assignment for a skip of 0

      If (skip <= 0) then
        n_greedy= n_greedy + 1;  nc= nl
        Cur_Adj= Last_Adj

        Call Greedy_cliq (v1,Edge, Cur_Adj, nc,Cur_cliq)

        Call Test_cliq ("Not_so_Greedy", Cur_cliq(:nc),Value, &
                        Connect, Best_cliq)
        Return
      End if
    End if

    skp= skip;  n1= nl + 1
    
    Vertex_loop : Do v= v1,nv
      If (.not.Last_Adj(v)) Cycle
      
      Cur_cliq(n1)= v
        
      If (v < nv) then
        v2= v + 1;  Cur_Adj(v2:)= Last_Adj(v2:) .and. Edge(v,v2:)
        no= Count(Cur_Adj(v2:))

        If (no > 0) then
          Call Not_so_Greedy_cliq (nv,v2,skp, Edge,Value,Connect, &
                                   Cur_Adj(v2:), n1,Cur_cliq, Best_cliq)
        End if

        If (n1 > 1) skp= skp - 1
      Else
        no= 0
      End if
        
      If (no <= 0 .or. skp < 0) then
        Call Test_cliq ("Not_so_Greedy ", Cur_cliq(:n1),Value, &
                        Connect, Best_cliq)
      End if
    End do Vertex_loop
    
  End Subroutine Not_so_Greedy_cliq
  
  Pure Subroutine Greedy_cliq (v1,Edge, Adj, nl,cliq)

!   Greedy algorithm for computing a maximal clique with a high sum of values.  
!   The vertices are assumed to be ordered according to decreasing liklihood 
!   of being in an an optimal maximal clique, by value for example

    Integer,         Intent(in) :: v1         ! Add to the clique from this vertex forward
    Logical,         Intent(in) :: Edge(:,:)  ! (nv,nv) Connectivity graph, diagonal not used
    
    Logical,      Intent(inout) :: Adj(v1:)   ! (v1:nv) Adjacency set: these vertices, if true, are
                                              !   adjacent to all the vertices already in the clique
    Integer,      Intent(inout) :: nl         ! # vertices in the clique
    Integer,      Intent(inout) :: cliq(:)    ! (nv) List of vertices in the clique
! Local:
    Integer :: v, v2, nv
    
    nv= Size(Edge,1)

    Do v= v1,nv
      If (Adj(v)) then
        nl= nl + 1;  cliq(nl)= v
        If (v < nv) then
          v2= v + 1;  Adj(v2:)= Adj(v2:) .and. Edge(v,v2:)
        End if
      End if
    End do
  End Subroutine Greedy_cliq

  Subroutine Test_cliq (Label,cliq, Value,Connect, Best_cliq)

!   Test for a good or improved clique and insert it in Best_cliq so that 
!   the clique value Best_cliq%sum_vl is decreaseing.

    Character(*),      Intent(in) :: Label        ! Print out label
    Integer,           Intent(in) :: cliq(:)      ! (nl) The current clique, in increasing order
    Real,              Intent(in) :: Value(:)     ! (nv) Vertex values
    Real,              Intent(in) :: Connect(:,:) ! (nv,nv) Edge connectivity matrix,
                                                  !         diagonal not used
    Type(Set_list), Intent(inout) :: Best_cliq(:) ! (mb) Top cliques from all the components of size >= 3
                                                  !      ordered by decreasing connectivitiy weighted value %smb
                                                  !   %n      = clique size
                                                  !   %svl    = clique value = sum of vertex values
                                                  !   %smb    = clique connectivity weighted value = sum of 
                                                  !             connectivities weighted by vertex values
                                                  !   %set(n) = clique vertices
                                                  !   %val(n) = clique vertex values
                                                  !   %mbr(n) = clique vertex memberships = normalized 
                                                  !             mean connectivities
 ! Local:
    Logical, Allocatable :: Adj(:)
    Type(Set_list) :: Clq
    Integer :: i, j, n, mb, nb, i1, ni, nl, Intr(Size(cliq))

    nl= Size(cliq);  mb= Size(Best_cliq);  nb= Last_true(Best_cliq%n > 0)
    
!   First compute the data for the new clique 'Clq', including
!   connectivity weighted memberships %mbr and the total
!   connectivity weighted value %smb.

    Call Memb_cliq (nl,cliq, Value,Connect, Clq)                
    
!   Now test to see if the new clique is a subset of a prior clique
!   If so, skip it as repeated or not maximal.
    
    Do i= 1,nb
      If (nl <= Best_cliq(i)%n) then
        j= Subset(cliq, Best_cliq(i)%set, ni,Intr)
        If (j == 0 .or. j == 1) Return
      End if
    End do

!   Next: shift the current cliques so that the new clique 
!   may be inserted if feasible

    i= First_true(1, Best_cliq(:nb)%smb < Clq%smb)
    
    If (i <= mb) then
      nb= Min(nb + 1, mb);  i1= i + 1
      Do j= nb,i1,-1
        Call Copy_set (Best_cliq(j-1), Best_cliq(j))       
      End do

!     Insert the new set at position i
      Call Copy_set (Clq, Best_cliq(i))
!      Call Copy1_set (Clq, Best_cliq(i), -1)
    End if

    Call DeAlloc_set_list (Clq)

    If (pr_out > 1.5) then
      If (i <= mb) then
        Call Out (Label//" new best clique of size",nl, ln=1)
        Call Out ("and vertices",Best_cliq(i)%set)
        Call Out ("inserted at position",i)
        Call Out ("Current best clique sizes",Best_cliq(:nb)%n)
        Call Out ("Current best clique connectivity weighted values",Best_cliq(:nb)%smb)
      Else  
        Call Out (Label//" best clique list full")
        Call Out ("Current best clique sizes",Best_cliq(:nb)%n)
        Call Out ("Current best clique connectivity weighted values",Best_cliq(:nb)%smb)
      End if
    End if
  End Subroutine Test_cliq
    
  Subroutine Memb_cliq (n,cliq, Value,Connect, Clq)             
     
!   Compute the full membership data of a clique. The membership value of a vertex
!   of a clique is proportional to its average connectivity in the clique weighted 
!   by corresponding vertex values
 
    Integer,         Intent(in) :: n            ! # vertices in the clique
    Integer,         Intent(in) :: cliq(:)      ! (n)  Clique vertices
    Real,            Intent(in) :: Value(:)     ! (nv) Vertex values
    Real,            Intent(in) :: Connect(:,:) ! (nv,nv) Edge connectivity matrix,
                                                !         diagonal not used
    Type(Set_list), Intent(out) :: Clq          ! Full clique data
                                                !   %n      = clique size
                                                !   %svl    = clique value = sum of vertex values
                                                !   %smb    = clique connectivity weighted value = sum of 
                                                !             connectivities weighted by vertex values
                                                !   %set(n) = clique vertices
                                                !   %val(n) = clique vertex values
                                                !   %mbr(n) = clique vertex memberships = normalized 
                                                !             mean connectivities
!  Local:
     Real    :: conn(n)  ! Average weighted connectivity of each vertex 
                         ! to every other vertex, using 'Connect' for 
                         ! the connectivity and 'Value' for the weights
     Real    :: wt(n)    ! The mean edge connectivities 'conn' weighted by 
                         ! the vertex value, so %mbr= wt normalized by %smb= Sum(wt)
     Logical :: Msk(n)
     Integer :: i
     
     Allocate (Clq%set(n), Clq%val(n), Clq%mbr(n))
     
     Clq%n= n;  Clq%set= cliq(:n)
     Clq%val= Value(Clq%set);  Clq%svl= Sum(Clq%val)
     
     If (n <= 1) then
       Clq%smb= Clq%val(1);  Clq%mbr= 1
     Else if (n == 2) then
       wt= Connect(Clq%set(1),Clq%set(2)) * Clq%val
       Clq%smb= Sum(wt);  Clq%mbr= wt / Clq%smb
     Else
       Msk= .true.
       Do i= 1,n
         Msk(i)= .false.;  If (i > 1) Msk(i-1)= .true. 
         conn(i)= (Sum(Connect(Clq%set(i),Clq%set) * Clq%val, Msk)) / &
                  (Clq%svl-Clq%val(i))
       End do
       
       wt= conn * Clq%val;  Clq%smb= Sum(wt);  Clq%mbr= wt / Clq%smb
     End if
     
  End Subroutine Memb_cliq                      
    
End Module Graph_algorithms

    