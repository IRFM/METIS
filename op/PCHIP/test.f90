program toto
   implicit none
   integer,parameter    :: N=7,INCFD=1,NE=20
   real*8,dimension(N)  :: X, F, D
   real*8,dimension(NE) :: XE, FE
   integer :: IERR,i,j
   logical :: SKIP
   real*8  :: delta
   
   do i = 1, N
     X(i)=i-4
     if (i==4) F(i)=0
     if (i<4)  F(i)=-1
     if (i>4)  F(i)=1
   enddo   
   
   XE(1)=-3
   XE(NE)=3
   delta = 6.0/NE
   do i = 2, NE-1
     XE(i)=-3+delta*i
   enddo  
   call cos_pchip (N, X, F, NE, XE, FE)
   
   print*,"FE :",FE(1:NE)
   
end program

subroutine mexErrMsgTxt(MESSG)
        CHARACTER	:: MESSG
	STOP
end subroutine mexErrMsgTxt