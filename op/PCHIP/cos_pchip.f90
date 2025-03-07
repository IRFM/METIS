subroutine cos_pchip(N,x,y,Ni,xi,yi)
  implicit none
  integer,intent(in)   :: N,Ni
  real*8,dimension(n)  :: x,y,d
  real*8,dimension(ni) :: xi,yi
  
  integer,parameter    :: INCFD=1
  integer :: IERR
  logical :: SKIP=.true.
  
  call DPCHIM (N, X, y, d, INCFD, IERR)   
  !if (ierr > 0) print*,"PCHIP:DPCHIM warning:",ierr
  if (ierr < 0) print*,"PCHIP:DPCHIM error:",ierr
  
  call DPCHFE (N, X, y, d, INCFD, SKIP, NI, xi, yi, IERR)
  !if (ierr > 0) print*,"PCHIP:DPCHFE warning:",ierr
  if (ierr < 0) print*,"PCHIP:DPCHFE error:",ierr
  
end subroutine cos_pchip 
