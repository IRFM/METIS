program FACIT_interface
use constants
IMPLICIT NONE

integer           :: i, k, nx, nth, nvar, iref
real(rkind), allocatable, dimension(:)  :: RHON, TE, TI, NE, NI, Na, QPSI, Machi, TperpsTpar, FV, dpsidx
real(rkind), allocatable, dimension(:,:)  :: AsymPhi, AsymN
real(rkind), allocatable, dimension(:,:)  :: PhiV, NV, BV, RV, nn, jacob
real(rkind), allocatable, dimension(:)  :: gradNi, gradTi, gradNa, Za
real(rkind), allocatable, dimension(:)  :: theta
real(rkind), allocatable, dimension(:)  :: dmin, dmaj
real(rkind)              :: Ai, Aa, Zi, Za_axis
real(rkind)              :: tau, B0, R0, invaspect, Te0, Ni0, qa, tedge, nedge, cimp
real(rkind)              :: fH, bC, Zeff, sigH, TperpsTpar_axis, Machi_axis
real(rkind)              :: time1, time2, xref, amin, nsigH
real(rkind), dimension(5):: DUMMY
real(rkind), dimension(4):: regulopt
logical 	  	 :: pol_asym, rotation, full_geom
     ! OUTPUTS

      real(rkind), allocatable, dimension(:) :: Flux_imp, Da, Vconv ! main outputs: flux, transport coefficients and poloidal asymmetry
 
      real(rkind), allocatable, dimension(:) :: Da_PS, Da_BP, Da_CL, Ka_PS, Ka_BP, Ka_CL, Ha_PS, Ha_BP, Ha_CL, Va_PS,Va_BP, Va_CL    ! supplementary outputs: PS, BP, CL & centrifugal components

 open(10,file='facit_input.dat',status='OLD')
 read(10,101) DUMMY(1), DUMMY(2), Aa,Ai,Zi,B0,R0,invaspect,DUMMY(3),DUMMY(4),DUMMY(5),regulopt(1),regulopt(2),regulopt(3),regulopt(4)         ! npx   ng    small
 nx = INT(DUMMY(1))
 nth = INT(DUMMY(2))

 allocate(RHON(nx), TE(nx), TI(nx), NE(nx), NI(nx), Na(nx), QPSI(nx), Machi(nx), Flux_imp(nx), Vconv(nx), Da(nx), dmin(nx), dmaj(nx))
 allocate(theta(nth), BV(nx,nth), RV(nx,nth), jacob(nx,nth), nn(nx,nth), PhiV(nx,nth), NV(nx,nth), FV(nx), dpsidx(nx))
 allocate(AsymPhi(nx,2), AsymN(nx,2)) 
 allocate(gradNi(nx), gradTi(nx), gradNa(nx), Za(nx))
 allocate(Da_PS(nx), Da_BP(nx), Da_CL(nx), Ka_PS(nx), Ka_BP(nx), Ka_CL(nx), Ha_PS(nx), Ha_BP(nx), Ha_CL(nx),Va_PS(nx),Va_BP(nx), Va_CL(nx))

     ! model options
      pol_asym  = INT(DUMMY(3))
      full_geom = INT(DUMMY(4))
      rotation  = INT(DUMMY(5))

 do i=1,nx
   read(10,101) RHON(i),Za(i),TE(i),TI(i),NE(i),NI(i),Na(i),QPSI(i),Machi(i),gradNi(i),gradTi(i),gradNa(i),FV(i),dpsidx(i),AsymPhi(i,1),AsymPhi(i,2),AsymN(i,1),AsymN(i,2)
 enddo
 if (full_geom) then
	if (nth.ne.64) then
		write(6,*) 'nth should be equal to 64 when in.geom == 1'
		stop
	endif
 	read(10,102) theta
	do i = 1,nx
 		read(10,102) (BV(i,k), k=1,nth)
	enddo
	do i = 1,nx
 		read(10,102) (RV(i,k), k=1,nth)
	enddo
	do i = 1,nx
 		read(10,102) (jacob(i,k), k=1,nth)
	enddo
	do i = 1,nx
 		read(10,102) (PhiV(i,k), k=1,nth)
	enddo
	do i = 1,nx
 		read(10,102) (NV(i,k), k=1,nth)
	enddo
 endif
 close(10)
 101  format(18('  ',E14.7))
 102  format(64('  ',E14.7))
! write(*,*) 'close input file'

 
 call cpu_time(time1)
 !     print *, "Calling FACIT"

      call FACIT(nx, nth, RHON, theta,Za, Aa, Zi, Ai,TE, TI, NE, NI, Na, Machi,gradTi, gradNi, gradNa,invaspect, B0, R0, QPSI, dpsidx, FV, BV, RV, jacob,AsymPhi, AsymN,pol_asym, full_geom, rotation, regulopt,Flux_imp, Da, Vconv, nn, dmin, dmaj,Da_PS, Da_BP, Da_CL, Ka_PS, Ka_BP, Ka_CL,Ha_PS, Ha_BP,Ha_CL,Va_PS,Va_BP,Va_CL)               ! supplementary outputs: PS, BP, CL & centrifugal components
 !     print *, "Done with FACIT call (new)"
 call cpu_time(time2)
 time2 = time2-time1
 !write(6,*) time2

  !   print *, "Write output file"
 open(20,file='facit_output.dat',position='REWIND')
 if (full_geom) then
 do i=1,nx
   write(20,203) RHON(i), Za(i), Aa, Zi, Ai, TE(i), TI(i), NE(i), NI(i), Na(i), gradTi(i), gradNi(i), gradNa(i), QPSI(i), Machi(i), AsymPhi(i,1), AsymPhi(i,2), AsymN(i,1),AsymN(i,2), Vconv(i), Da(i), dmin(i), dmaj(i),Da_PS(i), Da_BP(i), Da_CL(i), Ka_PS(i), Ka_BP(i), Ka_CL(i),Ha_PS(i), Ha_BP(i), Ha_CL(i),Va_PS(i),Va_BP(i),Va_CL(i),time2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
 enddo
 do i=1,nx
   write(20,203) (nn(i,k),k=1,nth)
 enddo
 else
 do i=1,nx
   write(20,202) RHON(i), Za(i), Aa, Zi, Ai, TE(i), TI(i), NE(i), NI(i), Na(i), gradTi(i), gradNi(i), gradNa(i), QPSI(i), Machi(i), AsymPhi(i,1), AsymPhi(i,2), AsymN(i,1),AsymN(i,2), Vconv(i), Da(i), dmin(i), dmaj(i),Da_PS(i), Da_BP(i), Da_CL(i), Ka_PS(i), Ka_BP(i), Ka_CL(i),Ha_PS(i), Ha_BP(i), Ha_CL(i),Va_PS(i),Va_BP(i),Va_CL(i),time2
 enddo
 endif
 close(20)
 202  format(36('  ',E14.7))
 203  format(64('  ',E14.7))
  !   print *, "Write output done"

 deallocate(RHON, TE, TI, NE, NI, Na, QPSI, Machi, Flux_imp, Vconv, Da, dmin, dmaj)
 deallocate(theta, BV, RV, jacob, nn, PhiV, NV, FV, dpsidx) 
 deallocate(AsymPhi, AsymN) 
 deallocate(gradNi, gradTi, gradNa, Za)
 deallocate(Da_PS, Da_BP, Da_CL, Ka_PS, Ka_BP, Ka_CL, Ha_PS, Ha_BP, Ha_CL,Va_PS,Va_BP, Va_CL)

end
