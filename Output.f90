
!          R.H. Burkhart 12-3-2001

 Module Output

!  This module contains the generic interface for several output routines.
!  There is also an error handling routine.

   Use ISO_Fortran_Env
   Use Precisn
   Implicit None
   Integer :: ipu= 8         ! output file unit number
   Real    :: pr_out= 0      ! output level if positive, otherwise no output.
   Integer :: nerr= 0        ! error counter for 'hherr'
   Integer :: nerlim= 1000   ! limit on 'nerr' for 'hherr'
   Character(1) :: Nmr(9) = (/"1","2","3","4","5","6","7","8","9"/)
   Character*20 :: fmt

   Interface Out
     Module procedure sout, xout, xdout, iout, lout, cout, irout, riout
     Module procedure vxout, vxdout, viout, vlout, vcout
     Module procedure mxout, mxdout, miout, mlout, mcout 
   End Interface Out

 Contains

   Subroutine hherr (mode,subnam, ier)
  
!    Error handling routine.  Requires initialization.
  
     Implicit None
     Integer,          Intent(in) :: mode, ier
     Character(len=*), Intent(in) :: subnam

     If (pr_out <= 0)  Return

     nerr= nerr + 1
     If (nerr > nerlim) then
       Write (Unit=ipu, Fmt=99)
   99    Format (/' ***** Error handler (hherr) message limit', &
         & ' exceeded'/ T10,'One or more messages', ' suppressed')
  
     Else
      Select case (mode)
       Case (0)
         Write (Unit=ipu, Fmt=101) subnam, subnam, ier
   101     Format (/' ***** Warning reported by Subroutine ',A/ &
           & T10,'See ',A,' abstract (ier =',I6,')' ) 
       Case (1)
         Write (Unit=ipu, Fmt=103) subnam, subnam, ier
   103     Format (/' ***** Input argument error reported by', &
           & ' Subroutine ',A/ T10,'See ',A,' abstract (ier =',I6,')' ) 
       Case (3)
         Write (Unit=ipu, Fmt=107) subnam, subnam, ier
   107     Format (/' ***** Process error reported by', ' Subroutine ',A/ &
           & T10,'See ',A,' abstract (ier =',I6,')' ) 
       Case (4)
         Write (Unit=ipu, Fmt=109) subnam, subnam, ier
   109     Format (/' ***** Above was called by', ' Subroutine ',A/ &
           & T10,'See ',A,' abstract (ier =',I6,')' ) 
       Case default
         Write (Unit=ipu, Fmt=120) mode, subnam
   120     Format (/' ***** Error handler (hherr) used improperly'/ &
           & T10,'mode =',I7,' from ',A)
         Stop
       End select
     End if
   End subroutine hherr


   Subroutine sout (message)

!    Write out the string 'message'

     Character(len=*), Intent(in) :: message

     If (ipu <= 0 .or. pr_out <= 0)  Return
     Write(ipu,'(/A)') message
   End subroutine sout


   Subroutine xout (label1,x1, label2,x2, fm,ln)

!    Write out the real number 'x1', and optionally, 'x2'
!    fm(1) = full length, fm(2) = # decimal places

     Character(len=*),           Intent(in) :: label1
     Real,                       Intent(in) :: x1
     Character(len=*), Optional, Intent(in) :: label2
     Real,             Optional, Intent(in) :: x2

     Integer,       Optional, Intent(in) :: fm(:)
     Integer,       Optional, Intent(in) :: ln
!  Local:     
     Integer :: i
     If (ipu <= 0 .or. pr_out <= 0)  Return

     If (Present(ln)) then
       Do i=1,ln;  Write(ipu,*);  End do
       End if

     If (Present(fm)) then
       If (Present(x2) .and. Present(label2)) then
         fmt= "(2(2X,A,F"//Nmr(fm(1))//"."//Nmr(fm(2))//"))"      
         Write(ipu,fmt) label1//'=',x1, label2//'=',x2
       Else
         fmt= "(A,F"//Nmr(fm(1))//"."//Nmr(fm(2))//"))"
         Write(ipu,fmt) label1//'=',x1
       End if
     Else
       If (Present(x2) .and. Present(label2)) then
         Write(ipu,'(2(2X,A,F0.3))') &
               label1//'= ',x1, label2//'= ',x2
       Else
         Write(ipu,'(A,F0.3)') label1//'= ',x1
       End if
     End if
   End subroutine xout


   Subroutine xdout (label1,x1, label2,x2, fm,ln)

!    Write out the real number 'x1', and optionally, 'x2'
!    fm(1) = full length, fm(2) = # decimal places

     Character(len=*),           Intent(in) :: label1
     Real(Dblp),              Intent(in) :: x1
     Character(len=*), Optional, Intent(in) :: label2
     Real(Dblp),    Optional, Intent(in) :: x2

     Integer,       Optional, Intent(in) :: fm(:)
     Integer,       Optional, Intent(in) :: ln
!  Local:     
     Integer :: i
     If (ipu <= 0 .or. pr_out <= 0)  Return

     If (Present(ln)) then
       Do i=1,ln;  Write(ipu,*);  End do
     End if

     If (Present(x2) .and. Present(label2)) then
       Write(ipu,'(2(2X,A,F0.6))') &
                 label1//'= ',x1, label2//'= ',x2
     Else
       Write(ipu,'(A,F0.6)') label1//'= ',x1
     End if
   End subroutine xdout


   Subroutine iout (label1,i1, label2,i2, fm,ln)

!    Write out the integer 'i1', and optionally, 'i2'.
!    fm = integer length

     Character(len=*), Intent(in) :: label1
     Integer,          Intent(in) :: i1
     Character(len=*), Optional, Intent(in) :: label2
     Integer,          Optional, Intent(in) :: i2

     Integer, Optional, Intent(in) :: fm     ! Integer format
     Integer, Optional, Intent(in) :: ln
!  Local:     
     Integer :: i
     If (ipu <= 0 .or. pr_out <= 0)  Return

     If (Present(ln)) then
       Do i=1,ln;  Write(ipu,*);  End do
     End if

     If (Present(i2) .and. Present(label2)) then
       Write (ipu,'(2(2X,A,I0))') label1//'= ',i1, label2//'= ',i2
     Else
       Write (ipu,'(2X,A,I0)') label1//'= ',i1
     End if
   End subroutine iout


   Subroutine irout (label1,i1, label2,x2, fm,ln)

!    Write out the integer 'i1' followed by the real 'x2'.
!    fm(1) = full real length, fm(2) = # decimal places, 
!    fm(3) = integer length

     Character(len=*), Intent(in) :: label1
     Integer,          Intent(in) :: i1
     Character(len=*), Intent(in) :: label2
     Real,             Intent(in) :: x2

     Integer, Optional, Intent(in) :: fm(:)
     Integer, Optional, Intent(in) :: ln
!  Local:     
     Integer :: i, ni
     If (ipu <= 0 .or. pr_out <= 0)  Return

     If (Present(ln)) then
       Do i=1,ln;  Write(ipu,*);  End do
     End if

     If (Present(fm)) then
       ni= fm(1);  If (Size(fm) > 2) ni= fm(3)
       fmt= "(2X,A,I"//Nmr(ni)//", 2X,A,F"//Nmr(fm(1))//"."//Nmr(fm(2))//")"
       Write (ipu,fmt) label1//'=',i1, label2//'=',x2
     Else
       Write (ipu,'(2X,A,I0, 2X,A,F0.3)') label1//'= ',i1, &
                  label2//'= ',x2
     End if
   End subroutine irout


   Subroutine riout (label1,x1, label2,i2, fm,ln)

!    Write out the integer 'i1' followed by the real 'x2'.
!    fm(1) = full real length, fm(2) = # decimal places, 
!    fm(3) = integer length

     Character(len=*), Intent(in) :: label1
     Real,             Intent(in) :: x1
     Character(len=*), Intent(in) :: label2
     Integer,          Intent(in) :: i2

     Integer, Optional, Intent(in) :: fm(:)  ! (2) fm(1) = full length, fm(2) = # decimal places
     Integer, Optional, Intent(in) :: ln
!  Local:     
     Integer :: i, ni
     If (ipu <= 0 .or. pr_out <= 0)  Return

     If (Present(ln)) then
       Do i=1,ln;  Write(ipu,*);  End do
     End if

     If (Present(fm)) then
       ni= fm(1);  If (Size(fm) > 2) ni= fm(3)
       fmt= "(2X,A,F"//Nmr(fm(1))//"."//Nmr(fm(2))//", 2X,A,I"//Nmr(ni)//")"
       Write (ipu,fmt) label1//'=',x1, label2//'=',i2
     Else
       Write (ipu,'(2X,A,F0.3, 2X,A,I0)') label1//'= ',x1, &
                   label2//'= ',i2
     End if
   End subroutine riout


   Subroutine lout (label1,i1, label2,i2, ln)

!    Write out the integer 'i1', and optionally, 'i2'.

     Character(len=*), Intent(in) :: label1
     Logical,          Intent(in) :: i1
     Character(len=*), Optional, Intent(in) :: label2
     Logical,          Optional, Intent(in) :: i2

     Integer, Optional, Intent(in) :: ln
!  Local:     
     Integer :: i
     If (ipu <= 0 .or. pr_out <= 0)  Return

     If (Present(ln)) then
       Do i=1,ln;  Write(ipu,*);  End do
     End if

     If (Present(i2) .and. Present(label2)) then
       Write (ipu,'(2(4X,A,L3))') label1//'=',i1, label2//'=',i2
     Else
       Write (ipu,'(4X,A,L3)') label1//'=',i1
     End if
   End subroutine lout


   Subroutine cout (label1,c1, label2,c2, ln)

!    Write out the real number 'c1', and optionally, 'c2'

     Character(len=*),           Intent(in) :: label1
     Complex,                    Intent(in) :: c1
     Character(len=*), Optional, Intent(in) :: label2
     Complex,          Optional, Intent(in) :: c2

     Integer, Optional, Intent(in) :: ln
!  Local:     
     Integer :: i
     If (ipu <= 0 .or. pr_out <= 0)  Return

     If (Present(ln)) then
       Do i=1,ln;  Write(ipu,*);  End do
     End if

     If (Present(c2) .and. Present(label2)) then
       Write(ipu,79) label1//'=',c1, label2//'=',c2
     Else
       Write(ipu,81) label1//'=',c1
     End if

   79  Format (2(A,'(',F0.4,', ',F0.4,')'))
   81  Format (A,'(',F0.4,', ',F0.4,')')
   End subroutine cout


   Subroutine vxout (label, xx, m, fm)       ! Print out a real vector
!    fm(1) = full length, fm(2) = # decimal places

     Character(len=*),  Intent(in) :: label
     Real,              Intent(in) :: xx(:)
     Integer, Optional, Intent(in) :: m      ! Matrix dimension
     Integer, Optional, Intent(in) :: fm(:)  ! (2) fm(1) = full length, fm(2) = # decimal places
!  Local:
     Integer :: i, n, p, lm, nx
     
     If (ipu <= 0 .or. pr_out <= 0)  Return
     Write(ipu,*);  Write(ipu,*) '@'//label//':'

     If (Present(m)) then   ! 'xx' represents the upper triangle of a matrix
       lm= m*(m+1)/2;  nx= Size(xx);  p= 0
       If (lm < nx) then
         Call iout ("Error in vxout. Size of upper triangle",lm, &
                    "but Size(xx)",nx);  Stop
       End if
       
       Do i= 1,m
         n= p + 1;  p= n + m-i
         Write (ipu,'(I3)', Advance='NO') i
         If (Present(fm)) then
           fmt= "(10F"//Nmr(fm(1))//"."//Nmr(fm(2))//")"
           Write (ipu,fmt) xx(n:p)
         Else
           Write (ipu,'(8(2X,F0.3))') xx(n:p)
         End if
       End do
     Else                   ! Vector case
       If (Present(fm)) then
         fmt= "(10F"//Nmr(fm(1))//"."//Nmr(fm(2))//")"
         Write (ipu,fmt) xx
       Else 
         Write (ipu,'(8(2X,F0.3))') xx
       End if
     End if
   End subroutine vxout


   Subroutine vxdout (label, xd, m)          ! Print out a double precision vector
     Character(len=*),  Intent(in) :: label
     Real(Dblp),     Intent(in) :: xd(:)
     Integer, Optional, Intent(in) :: m      ! Matrix dimension
!  Local:
     Integer :: i, n, p, lm, nx
     
     If (ipu <= 0 .or. pr_out <= 0)  Return
     Write(ipu,*);  Write(ipu,*) '@'//label//':'

     If (Present(m)) then   ! 'xd' represents the upper triangle of a matrix
       lm= m*(m+1)/2;  nx= Size(xd);  p= 0
       If (lm < nx) then
         Call iout ("Error in vxdout. Size of upper triangle",lm, &
                    "but Size(xd)",nx);  Stop
       End if
       
       Do i= 1,m
         n= p + 1;  p= n + m-i
         Write (ipu,'(I3)', Advance='NO') i
         Write (ipu,'(6(2X,F0.6))') xd(n:p)
       End do
     Else                   ! Vector case
       Write (ipu,'(6(2X,F0.6))') xd
     End if
   End subroutine vxdout


   Subroutine viout (label, ix, fm)           ! Print out an integer vector
!    fm = integer length

     Character(len=*),  Intent(in) :: label
     Integer,           Intent(in) :: ix(:)
     Integer, Optional, Intent(in) :: fm  ! Matrix dimension
!  Local:
     Integer :: i, n, p, lm, nx
     
     If (ipu <= 0 .or. pr_out <= 0)  Return
     Write(ipu,*);  Write(ipu,*) '@'//label//':'

     If (Present(fm)) then   ! 'ix' represents the upper triangle of a matrix
       lm= fm*(fm+1)/2;  nx= Size(ix);  p= 0
       If (lm < nx) then
         Call iout ("Error in viout. Size of upper triangle",lm, &
                    "but Size(ix)",nx);  Stop
       End if
       
       Do i= 1,fm
         n= p + 1;  p= n + fm-i
         Write (ipu,'(I3)', Advance='NO') i
         Write (ipu,'(15(2X,I0))') ix(n:p)
       End do
     Else
       Write (ipu,'(15(2X,I0))') ix
     End if
   End subroutine viout


   Subroutine vlout (label, lx, fm)           ! Print out a logical vector
     Character(len=*),  Intent(in) :: label
     Logical,           Intent(in) :: lx(:)
     Integer, Optional, Intent(in) :: fm      ! Matrix dimension
!  Local:
     Integer :: i, n, p, lm, nx
     
     If (ipu <= 0 .or. pr_out <= 0)  Return
     Write(ipu,*);  Write(ipu,*) '@'//label//':'

     If (Present(fm)) then   ! 'lx' represents the upper triangle of a matrix
       lm= fm*(fm+1)/2;  nx= Size(lx);  p= 0
       If (lm < nx) then
         Call iout ("Error in vlout. Size of upper triangle",lm, &
                    "but Size(lx)",nx);  Stop
       End if
       
       Do i= 1,fm
         n= p + 1;  p= n + fm-i
         Write (ipu,'(I3)', Advance='NO') i
         Write (ipu,'(20L3)') lx(n:p)
       End do
     Else                   ! Vector case
       Write (ipu,'(20L3)') lx
     End if
   End subroutine vlout


   Subroutine vcout (label, cx, fm)           ! Print out a complex vector
     Character(len=*),  Intent(in) :: label
     Complex,           Intent(in) :: cx(:)
     Integer, Optional, Intent(in) :: fm      ! Matrix dimension
!  Local:
     Integer :: i, n, p, lm, nx
     
     If (ipu <= 0 .or. pr_out <= 0)  Return
     Write(ipu,*);  Write(ipu,*) '@'//label//':'

     If (Present(fm)) then  ! 'cx' represents the upper triangle of a matrix
       lm= fm*(fm+1)/2;  nx= Size(cx);  p= 0
       If (lm < nx) then
         Call iout ("Error in vcout. Size of upper triangle",lm, &
                    "but Size(cx)",nx);  Stop
       End if
       
       Do i= 1,fm
         n= p + 1;  p= n + fm-i
         Write (ipu,'(I3)', Advance='NO') i
         Write (ipu,79) cx(n:p)
       End do
     Else                  ! Vector case
       Write (ipu,79) cx
     End if

     79  Format (4('(',F0.4,', ',F0.4,')',2X))
   End subroutine vcout


   Subroutine mxout (iopt,label, xx, row1,fm)
  
!    Print out the real matrix 'xx' to unit ipu,
!    or the lower triangle of 'xx' as a symmetric matrix.
!    fm(1) = full length, fm(2) = # decimal places
  
     Integer,           Intent(in) :: iopt      ! Output by rows if 'iopt' >= 1, by columns if 'iopt' <= -1,
                                                !  starting with the first row if 'iopt' == 1, or
                                                !  with the first column if 'iopt' == -1.
                                                ! Output the lower triangle of 'xx' if iopt == 0.
     Character(len=*),  Intent(in) :: label
     Real,              Intent(in) :: xx(:,:)   ! (n,m)
     Integer, Optional, Intent(in) :: row1      ! Label of the first row (or column)
     Integer, Optional, Intent(in) :: fm(:)     ! (2) format specifier

!   Local:
     Integer :: i, j, n, m, r0
  
     If (ipu <= 0 .or. pr_out <= 0)  Return
     n= Ubound(xx,1);  m= Ubound(xx,2)
     
     Write(ipu,*);  Write(ipu,*) '@'//label//':'
     
     r0= 0;  If (Present(row1)) r0= row1 - 1
  
     If (iopt >= 1) then
       If (iopt == 1) then
         Do i= 1,n
           Write (ipu,'(I3,A)', Advance='NO') r0+i, ":"
           If (Present(fm)) then
             fmt= "(10F"//Nmr(fm(1))//"."//Nmr(fm(2))//")"
             Write (ipu,fmt) xx(i,:)
           Else 
             Write (ipu,'(8(2X,F0.3))') xx(i,:)
           End if  
         End do
       Else
         Do i= n,1,-1
           Write (ipu,'(I3,A)', Advance='NO') r0+i, ":"
           If (Present(fm)) then
             fmt= "(10F"//Nmr(fm(1))//"."//Nmr(fm(2))//")"
             Write (ipu,fmt) xx(i,:)
           Else 
             Write (ipu,'(8(2X,F0.3))') xx(i,:)
           End if  
         End do
       End if
     Else if (iopt == 0) then
       Do i= 1,n
         j= Min(i,m)
         Write (ipu,'(I3,A)', Advance='NO') r0+i, ":"
         If (Present(fm)) then
           fmt= "(10F"//Nmr(fm(1))//"."//Nmr(fm(2))//")"
           Write (ipu,fmt) xx(:i,j)
         Else 
           Write (ipu,'(8(2X,F0.3))') xx(:i,j)
         End if  
       End do
     Else if (iopt <= -1) then
       If (iopt == -1) then
         Do j= 1,m
           Write (ipu,'(I3,A)', Advance='NO') r0+j, ":"
           If (Present(fm)) then
             fmt= "(10F"//Nmr(fm(1))//"."//Nmr(fm(2))//")"
             Write (ipu,fmt) xx(:,j)
           Else 
             Write (ipu,'(8(2X,F0.3))') xx(:,j)
           End if  
         End do
       Else
         Do j= m,1,-1
           Write (ipu,'(I3,A)', Advance='NO') r0+j, ":"
           If (Present(fm)) then
             fmt= "(10F"//Nmr(fm(1))//"."//Nmr(fm(2))//")"
             Write (ipu,fmt) xx(:,j)
           Else 
             Write (ipu,'(8(2X,F0.3))') xx(:,j)
           End if  
         End do
       End if
     End if
   End subroutine mxout

   
   Subroutine mxdout (iopt,label, xd, row1)
  
!    Print out the double precision matrix 'xd' to unit ipu,
!    or the lower triangle of 'xd' as a symmetric matrix.
  
     Integer,           Intent(in) :: iopt      ! Output by rows if 'iopt' >= 1, by columns if 'iopt' <= -1,
                                                !  starting with the first row if 'iopt' == 1, or
                                                !  with the first column if 'iopt' == -1.
                                                ! Output the lower triangle of 'xd' if iopt == 0.
     Character(len=*),  Intent(in) :: label
     Real(Dblp),     Intent(in) :: xd(:,:)   ! (n,m)
     Integer, Optional, Intent(in) :: row1      ! First row (or column) to be output

!   Local:
     Integer :: i, j, n, m, r0
  
     If (ipu <= 0 .or. pr_out <= 0)  Return
     n= Ubound(xd,1);  m= Ubound(xd,2)
     
     Write(ipu,*);  Write(ipu,*) '@'//label//':'
     r0= 0;  If (Present(row1)) r0= row1 - 1
  
     If (iopt >= 1) then
       If (iopt == 1) then
         Do i= 1,n
           Write (ipu,'(I3,A)', Advance='NO') r0+i, ":"
           Write (ipu,'(6(2X,F0.6))') xd(i,:)
         End do
       Else
         Do i= n,1,-1
           Write (ipu,'(I3,A)', Advance='NO') r0+i, ":"
           Write (ipu,'(6(2X,F0.6))') xd(i,:)
         End do
       End if
     Else if (iopt == 0) then
       Do i= 1,n
         j= Min(i,m)
         Write (ipu,'(I3,A)', Advance='NO') r0+i, ":"
         Write (ipu,'(I6)') i
         Write (ipu,'(6(2X,F0.6))') xd(i,:j)
       End do
     Else if (iopt <= -1) then
       If (iopt == -1) then
         Do j= 1,m
           Write (ipu,'(I3,A)', Advance='NO') r0+j, ":"
           Write (ipu,'(6(2X,F0.6))') xd(:,j)
         End do
       Else
         Do j= m,1,-1
           Write (ipu,'(I3,A)', Advance='NO') r0+j, ":"
           Write (ipu,'(6(2X,F0.6))') xd(:,j)
         End do
       End if
     End if
   End subroutine mxdout

   
   Subroutine miout (iopt,label, ix, row1)
  
!    Print out the integer matrix 'ix' to unit ipu,
!    or the lower triangle of 'ix' as a symmetric matrix.
!    fm = integer length
  
     Integer,           Intent(in) :: iopt      ! Output by rows if 'iopt' >= 1, by columns if 'iopt' <= -1,
                                                !  starting with the first row if 'iopt' == 1, or
                                                !  with the first column if 'iopt' == -1.
                                                ! Output the lower triangle of 'ix' if iopt == 0.
     Character(len=*),  Intent(in) :: label
     Integer,           Intent(in) :: ix(:,:)   ! (n,m)
     Integer, Optional, Intent(in) :: row1      ! First row (or column) to be output
!   Local:
     Integer :: i, j, n, m, r0
  
     If (ipu <= 0 .or. pr_out <= 0)  Return
     n= Ubound(ix,1);  m= Ubound(ix,2)
     
     Write(ipu,*);  Write(ipu,*) '@'//label//':'
     r0= 0;  If (Present(row1)) r0= row1 - 1
  
      If (iopt >= 1) then
       If (iopt == 1) then
         Do i= 1,n
           Write (ipu,'(I3,A)', Advance='NO') r0+i, ":"
           Write (ipu,'(15(2X,I0))') ix(i,:)
         End do
       Else
         Do i= n,1,-1
           Write (ipu,'(I3,A)', Advance='NO') r0+i, ":"
           Write (ipu,'(15(2X,I0))') ix(i,:)
         End do
       End if
      Else if (iopt == 0) then
       Do i= 1,n
         j= Min(i,m)
         Write (ipu,'(I3,A)', Advance='NO') r0+i, ":"
         Write (ipu,'(15(2X,I0))') ix(i,:j)
       End do
      Else
       If (iopt == -1) then
         Do j= 1,m
           Write (ipu,'(I3,A)', Advance='NO') r0+j, ":"
           Write (ipu,'(15(2X,I0))') ix(:,j)
         End do
       Else
         Do j= m,1,-1
           Write (ipu,'(I3,A)', Advance='NO') r0+j, ":"
           Write (ipu,'(15(2X,I0))') ix(:,j)
         End do
       End if
      End if
   End subroutine miout

   
   Subroutine mlout (iopt,label, lx, row1)
  
!    Print out the logical matrix 'lx' to unit ipu,
!    or the lower triangle of 'lx' as a symmetric matrix.
  
     Integer,           Intent(in) :: iopt      ! Output by rows if 'iopt' >= 1, by columns if 'iopt' <= -1,
                                                !  starting with the first row if 'iopt' == 1, or
                                                !  with the first column if 'iopt' == -1.
                                                ! Output the lower triangle of 'lx' if iopt == 0.
     Character(len=*),  Intent(in) :: label
     Logical,           Intent(in) :: lx(:,:)   ! (n,m)
     Integer, Optional, Intent(in) :: row1      ! First row (or column) to be output

!   Local:
     Integer :: i, j, n, m, r0
  
     If (ipu <= 0 .or. pr_out <= 0)  Return
     n= Ubound(lx,1);  m= Ubound(lx,2)
     
     Write(ipu,*);  Write(ipu,*) '@'//label//':'
     r0= 0;  If (Present(row1)) r0= row1 - 1

     If (iopt >= 1) then
       If (iopt == 1) then
         Do i= 1,n
           Write (ipu,'(I3,A)', Advance='NO') r0+i, ":"
           Write (ipu,'(20L3)') lx(i,:)
         End do
       Else
         Do i= n,1,-1
           Write (ipu,'(I3,A)', Advance='NO') r0+i, ":"
           Write (ipu,'(20L3)') lx(i,:)
         End do
       End if
     Else if (iopt == 0) then
       Do i= 1,n
         j= Min(i,m)
         Write (ipu,'(I3,A)', Advance='NO') r0+i, ":"
         Write (ipu,'(20L3)') lx(i,:j)
       End do
     Else if (iopt <= -1) then
       If (iopt == -1) then
         Do j= 1,m
           Write (ipu,'(I3,A)', Advance='NO') r0+j, ":"
           Write (ipu,'(20L3)') lx(:,j)
         End do
       Else
         Do j= m,1,-1
           Write (ipu,'(I3,A)', Advance='NO') r0+j, ":"
           Write (ipu,'(20L3)') lx(:,j)
         End do
       End if
     End if
   End subroutine mlout

   
   Subroutine mcout (iopt,label, cx, row1)
  
!    Print out the complex matrix 'cx' to unit ipu,
!    or the lower triangle of 'cx' as a symmetric matrix.
  
     Integer,           Intent(in) :: iopt      ! Output by rows if 'iopt' >= 1, by columns if 'iopt' <= -1,
                                                !  starting with the first row if 'iopt' == 1, or
                                                !  with the first column if 'iopt' == -1.
                                                ! Output the lower triangle of 'cx' if iopt == 0.
     Character(len=*),  Intent(in) :: label
     Complex,           Intent(in) :: cx(:,:)   ! (n,m)
     Integer, Optional, Intent(in) :: row1      ! First row (or column) to be output

!   Local:
     Integer :: i, j, n, m, r0
  
     If (ipu <= 0 .or. pr_out <= 0)  Return
     n= Ubound(cx,1);  m= Ubound(cx,2)
     
     Write(ipu,*);  Write(ipu,*) '@'//label//':'
     r0= 0;  If (Present(row1)) r0= row1 - 1

     If (iopt >= 1) then
       If (iopt == 1) then
         Do i= 1,n
           Write (ipu,'(I3,A)', Advance='NO') r0+i, ":"
           Write (ipu,79) cx(i,:)
         End do
       Else
         Do i= n,1,-1
           Write (ipu,'(I3,A)', Advance='NO') r0+i, ":"
           Write (ipu,79) cx(i,:)
         End do
       End if
     Else if (iopt == 0) then
       Do i= 1,n
         j= Min(i,m)
         Write (ipu,'(I3,A)', Advance='NO') r0+i, ":"
         Write (ipu,79) cx(i,:j)
       End do
     Else if (iopt <= -1) then
       If (iopt == -1) then
         Do j= 1,m
           Write (ipu,'(I3,A)', Advance='NO') r0+j, ":"
           Write (ipu,79) cx(:,j)
         End do
       Else
         Do j= m,1,-1
           Write (ipu,'(I3,A)', Advance='NO') r0+j, ":"
           Write (ipu,79) cx(:,j)
         End do
       End if
     End if

     79  Format (4('(',F0.4,', ',F0.4,')',2X))
   End subroutine mcout

 End module Output
