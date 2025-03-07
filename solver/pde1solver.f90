module TransportData
  implicit none
  integer :: dimk,dimm,dimmat,dimbb,dimn
  real(kind=8),dimension(:,:),  pointer       :: f,fp,dfpdx,dfpdx2
  real(kind=8),dimension(:,:),allocatable     :: matd,matdp
  real(kind=8),dimension(:,:,:),pointer       :: T0,T1
  real(kind=8),dimension(:,:,:),allocatable   :: mata,matb,matc,    &
       matap,matbp,matcp
  real(kind=8),dimension(:,:),  pointer       :: V0,V1
  integer,dimension(:),allocatable            :: mode

  real(kind=8),dimension(:,:),allocatable     :: mats
  real(kind=8),dimension(:,:),allocatable     :: matu,matv,matw,matx,maty,matz
  real(kind=8),dimension(:),allocatable,target:: mat_coo,mat_csr
  real(kind=8),dimension(:),allocatable,target:: vectbb
  integer,dimension(:),allocatable,target     :: cooi,cooj,csri,csrj
  real(kind=8), pointer:: dx,dt,sca_f
  integer :: verbose = 0
end module TransportData

subroutine pde1solver(f_in,fp_in,          &
     mata_in,matb_in,matc_in,matd_in,      &
     matap_in,matbp_in,matcp_in,matdp_in,  &
     dimk_in,dimn_in,dimm_in,dx_in,dt_in,sca_f_in, &
     T0_in,T1_in,V0_in,V1_in,mode_in,dfpdx_in,dfpdx2_in)
  !*********************************************************************
  !solve the equations :
  !  Diff(f,t) = mata * diff(f,x,x) + matb * diff(f,x) + matc * f + matd 
  !
  !dt_in    : time step
  !dx_in    : space step
  !
  !f        : solution at the time t
  !fp       : solution at the time t+dt
  !dfpdx    : diff(fp,x)
  !dfpdx2   : diff(fp,x,x)
  !
  !dimk_in  : dimension of space index
  !  x=x0 + (k-1)*dx, k->[1..dimk]
  !dimm_in  : dimension of equation number 
  !dimn_in  : dimension of the solution at a one grid point
  !
  !sca_f_in : scalar 
  !           0   -> full implicit
  !           0.5 -> Crank Nicholson
  !           1   -> explicit (don't used it)
  !
  !
  !The big matrix
  !--------------
  !
  ! k=1               k=2  ----------------------------------------- k=K
  ! (1,1)--------(1,m)(1,1)--------(1,m)
  ! | b(1,m,n)   |     | c(1,m,n)   |
  ! |            |     |            |
  ! (n,1)--------(n,m)(n,1)--------(n,m)
  ! (1,1)--------(1,m)(1,1)--------(1,m)(1,1)--------(1,m)
  ! | a(2,m,n)   |     | b(2,m,n)   |    | c(2,m,n)   |
  ! |            |     |            |    |            |
  ! (n,m)--------(n,m)(n,m)--------(n,m)(n,m)--------(n,m)
  !       .                               
  !       .                                
  !       .                                 
  !       .                                  
  !       .                                   
  !       .                                    
  !       .                                    
  !       .                                      
  !       .                                       
  !
  !                                               (1,1)--------(1,m)(1,1)--------(1,m)
  !                                                |  a(K,m,n)  |     | b(K,m,n)   |
  !                                                |            |     |            |
  !                                               (n,m)--------(n,m)(n,m)--------(n,m)
  !
  ! where n=m if we are in predictive mode 
  !       n<m if we are in interpretative mode
  ! the first index is the space index
  ! the second index is the coupling coefficient between variables
  ! the third index is the rank of the equation
  ! 
  ! T0(l,m,j) : coefficients of the boundaries condition at x0
  ! V0(l,m)   : values of the boundaries condition at x0 
  !        T0(:,:,3) * diff(f,x,x) +  T0(:,:,2) * diff(f,x) + T0(:,:,1) * f = V0  
  !   where l=1 for t
  !         l=2 for t+dt
  ! T1(l,m,j) , V1(l,m) : idem as above but at x1
  !*************************************************************************************
  use TransportData
  implicit none
  integer*4, intent(in)            :: dimk_in,dimn_in,dimm_in
  real(kind=8),intent(in),target :: dt_in,dx_in,sca_f_in
  real(kind=8),dimension(dimk_in,dimm_in),target         :: f_in,fp_in
  real(kind=8),dimension(dimk_in,dimm_in),target         :: dfpdx_in,dfpdx2_in, &   
       matd_in,matdp_in
  real(kind=8),dimension(dimk_in,dimm_in,dimm_in),target :: mata_in,matb_in,matc_in, &
       matap_in,matbp_in,matcp_in
  real(kind=8),dimension(2,dimm_in,3),target   :: T0_in,T1_in
  real(kind=8),dimension(2,dimm_in),target     :: V0_in,V1_in

  real(kind=8),dimension(:,:,:), allocatable :: mat_aux

  real(kind=8),dimension(dimm_in)  :: mode_in
  integer,dimension(:),allocatable :: iwk,cooi_aux

  real(kind=8) :: dtx2,dt2x,aux
  integer :: i,j,k,m,n
  integer :: counteri,counterj,counter
  integer :: interpretatif=0,condlim=0

  interface
     integer function condlimord2_k1(interpretatif)
       integer,intent(in) :: interpretatif
     end function condlimord2_k1

     integer function condlimord2_kk(interpretatif)
       integer,intent(in) :: interpretatif
     end function condlimord2_kk
  end interface

  if (.false.) call mexPrintf("STARTING PDE1SOLVER")

  if (sca_f_in==1.0) then
     print*,"Error, the explicit Euler algorithm must not be used"
     stop
  endif

  !**********
  !allocation
  !**********
  dimk=dimk_in
  dimm=dimm_in
  dimn=dimn_in
  sca_f=>sca_f_in
  dx=>dx_in
  dt=>dt_in
  f=>f_in
  fp=>fp_in
  dfpdx=>dfpdx_in
  dfpdx2=>dfpdx2_in
  allocate(mata(dimk_in,dimm_in,dimm_in),matb(dimk_in,dimm_in,dimm_in),&
           matc(dimk_in,dimm_in,dimm_in),matd(dimk_in,dimm_in) )
  allocate(matap(dimk_in,dimm_in,dimm_in),matbp(dimk_in,dimm_in,dimm_in),&
           matcp(dimk_in,dimm_in,dimm_in),matdp(dimk_in,dimm_in) )
  allocate(mode(dimm_in))

  mata=mata_in
  matb=matb_in
  matc=matc_in
  matd=matd_in

  matap=matap_in
  matbp=matbp_in
  matcp=matcp_in
  matdp=matdp_in

  mode=int(mode_in)
  T0=>T0_in
  T1=>T1_in
  V0=>V0_in
  V1=>V1_in
  allocate(mats(dimk,dimm))
  allocate(matu(2,dimm),matv(2,dimm),matw(2,dimm))
  allocate(matx(2,dimm),maty(2,dimm),matz(2,dimm))
  allocate(mat_aux(dimk,dimm,dimm))

  dimmat=(dimk-2)*dimn*dimn*3 +4*dimn*dimn
  dimbb=dimn*dimk
  allocate(vectbb(dimbb),mat_coo(dimmat),cooi(dimmat),cooj(dimmat))
  if (.false.) call mexPrintf("ALLOCATED")

  !**************    
  !initialisation
  !**************
  interpretatif=sum(mode)
  if (interpretatif.ne.0) interpretatif=1
  dtx2  = dt/(dx*dx)
  dt2x  = dt/(2.0*dx)

  mat_aux = matb
  matb  = -2.0*dtx2*mata + dt*matc
  matc  = dtx2*mata      + dt2x*mat_aux
  mata  = dtx2*mata      - dt2x*mat_aux   
  matd  = dt*matd

  mat_aux = matbp
  matbp = -2.0*dtx2*matap + dt*matcp
  matcp = dtx2*matap      + dt2x*mat_aux
  matap = dtx2*matap      - dt2x*mat_aux
  matdp = dt*matdp

  deallocate(mat_aux)

  mats = sca_f*matd + (1.0-sca_f)*matdp + f; 

  !*********************
  !boundaries conditions
  !*********************
  condlim=sum(T0(:,:,3))
  if ( condlimord2_k1(interpretatif) > 0) then
     if (condlim>0) then
        print*," boundaries condition of order 2 is not implemented"
        stop
     else
        !print*,"boundaries condition of order 1 k=1"
        call condlimord1_k1
     endif
  endif

  condlim=sum(T1(:,:,3))
  if ( condlimord2_kk(interpretatif) > 0) then
     if (condlim>0) then
        print*," boundaries condition of order 2 is not implemented"
        stop
     else
        !print*,"boundaries condition of order 1 k=dimk"
        call condlimord1_kk
     endif
  endif
  if (.false.) call mexPrintf("Hello A1")


  !*******************
  !interpretative mode
  !*******************
  do m = 1,dimm
     if (mode(m)==1) cycle
     do k=2,dimk-1
        do n=1,dimm
           if (mode(n)==1) then
              mats(k,m)=mats(k,m) + sca_f*(mata(k,n,m)*f(k-1,n)+matb(k,n,m)*f(k,n)+matc(k,n,m)*f(k+1,n)) &
                   + (1.0-sca_f)*(matap(k,n,m)*fp(k-1,n)+matbp(k,n,m)*fp(k,n)+matcp(k,n,m)*fp(k+1,n))
           endif
        enddo
     enddo
  enddo

  !*************************
  !filling the second member
  !*************************
  !for k=1
  k=1
  counter=0
  do m=1,dimm
     if (mode(m)==0) then
        counter=counter+1
        vectbb(counter) = mats(k,m)
        do n=1,dimm
           if (mode(n)==0) then
              vectbb(counter) = vectbb(counter) + sca_f*(matb(k,n,m)*f(k,n) + matc(k,n,m)*f(k+1,n))
           endif
        enddo
     endif
  enddo

  !for k =[2,dimk-1]
  do k=2,dimk-1
     do m=1,dimm
        if (mode(m)==0) then 
           counter=counter+1
           vectbb(counter) = mats(k,m) 
           do n=1,dimm
              if (mode(n)==0) then
              vectbb(counter) = vectbb(counter) + sca_f*(mata(k,n,m)*f(k-1,n) + matb(k,n,m)*f(k,n) &
                   + matc(k,n,m)*f(k+1,n))
              endif
           enddo
        endif
     enddo
  enddo
  
  !for k=dimk
  k=dimk
  do m=1,dimm
     if (mode(m)==0) then
        counter=counter+1
        vectbb(counter) = mats(k,m) 
        do n=1,dimm
           if (mode(n)==0) then
           vectbb(counter) = vectbb(counter) + sca_f*(mata(k,n,m)*f(k-1,n) + matb(k,n,m)*f(k,n))
           endif
        enddo
     endif
  enddo

  !***********************
  !filling the big matrix 
  !***********************
  k=1
  counter =0
  counteri=0
  counterj=0
  call build_matcoo(matbp,1,k,counter,counteri,counterj)
  counteri=0
  counterj=dimn
  call build_matcoo(matcp,0,k,counter,counteri,counterj)
  do k=2,dimk-1
     counteri=(k-1)*dimn
     counterj=(k-2)*dimn
     call build_matcoo(matap,0,k,counter,counteri,counterj)
     counteri=(k-1)*dimn
     counterj=(k-1)*dimn
     call build_matcoo(matbp,1,k,counter,counteri,counterj)
     counteri=(k-1)*dimn
     counterj=k*dimn
     call build_matcoo(matcp,0,k,counter,counteri,counterj)
  enddo
  k=dimk
  counteri=(k-1)*dimn
  counterj=(k-2)*dimn
  call build_matcoo(matap,0,k,counter,counteri,counterj)
  counteri=(k-1)*dimn
  counterj=(k-1)*dimn
  call build_matcoo(matbp,1,k,counter,counteri,counterj)
  !*****************
  ! print the matrix
  !*****************
  if (verbose>0) then
     allocate(iwk(dimbb+1),cooi_aux(dimmat))
     allocate(mat_csr(dimmat),csri(dimmat),csrj(dimmat))
  endif
  dimmat=counter
  if (verbose>0) then
      print*,"SPARSE BIG MATRIX:"
     cooi_aux=cooi
     !  call coicsr (dimbb,dimmat,1,mat_coo,cooj,cooi,iwk)
     !call coocsr(dimbb,dimmat,mat_coo,cooi_aux,cooj,mat_csr,csrj,csri)
     !  call dump (1,dimbb,1,mat_coo,cooj,cooi,6)
     !call dump (1,dimbb,1,mat_csr,csrj,csri,6)
     print*,"RHS:"
     print*,vectbb
  endif
  !****************
  ! solve de matrix
  !****************
  if (.false.) call mexPrintf("Hello A2")
  call solvepde1d
  if (.false.) call mexPrintf("Hello A3")

  !*****************
  ! save de solution
  !*****************
  counter=0
  do k=1,dimk
     do m=1,dimm
        if (mode(m)==0) then
           counter=counter+1
           fp(k,m) = vectbb(counter)
        endif
     enddo
  enddo

  !fp(1,1)=mat_coo(1)
  !fp(2,1)=mat_coo(4)
  !fp(3,1)=aux
  !fp(4,1)=dt2x

  !do m=1,dimm
  !   if (mode(m) == 1) then
  !      do k=1,dimk
  !         fp(k,m) = fp_in(k,m)
  !      enddo
  !   endif
  !enddo
  if (verbose > 0) then
     print*,"SOLUTION:"
     print*,vectbb
  endif

  !************
  !finalisation
  !************
  deallocate(matu,matv,matw,matx,maty,matz,mats)
  deallocate(mata,matb,matc,matd,matap,matbp,matcp,matdp)
  deallocate(mode)
  if (verbose==0) then
     deallocate(mat_coo,cooi,cooj,vectbb)
  else
     deallocate(mat_coo,cooi,cooj,vectbb,iwk,cooi_aux,mat_csr,csri,csrj)
  endif

  if (dfpdx(1,1)==0.and.dfpdx2(1,1)==0) then
     print*,"dfpdx and dfpdx2 not available"
  endif

end subroutine pde1solver

subroutine build_matcoo(matp,option,k,counter,counteri,counterj)
  use TransportData
  implicit none
  integer,intent(in)            :: option,k
  integer,intent(inout)         :: counter,counteri,counterj
  integer                       :: counteri_aux,counterj_aux
  real(kind=8),dimension(dimk,dimm,dimm) :: matp

  integer :: m,n

  counteri_aux=counteri
  counterj_aux=counterj
  do m=1,dimm
     if (mode(m)==0) then     
        counteri=counteri+1
        counterj=counterj_aux
        do n=1,dimm
           if (mode(n)==0) then
              counterj=counterj+1
              counter=counter+1
              if (counteri==counterj) then
                 mat_coo(counter)=1.0-(1.0-sca_f)*matp(k,n,m)
              else
                 mat_coo(counter)=-(1.0-sca_f)*matp(k,n,m)
              endif
              if (mat_coo(counter)==0.0) then 
                 counter=counter-1
              else                
                 cooi(counter)=counteri
                 cooj(counter)=counterj
              endif
           endif
        enddo
     endif
  enddo
end subroutine build_matcoo

subroutine solvepde1d
  use TransportData
  implicit none
  include 'mpif.h'
  include 'dmumps_struc.h'
  type(DMUMPS_STRUC) :: mumps_par
  integer:: ierr, i

  call MPI_init(ierr)
  ! Define a communicator for the package.
  mumps_par%COMM = MPI_COMM_WORLD
  ! Initialize an instance of the package
  ! for L U factorization (sym = 0, with working host)
  mumps_par%JOB = -1
  mumps_par%SYM = 0
  mumps_par%PAR = 1
  call DMUMPS(mumps_par)
  if (.false.) call mexPrintf("Hello B1")

  !  Define problem on the host (processor 0)
  IF ( mumps_par%MYID .eq. 0 ) THEN
     !order of the matrix
     mumps_par%N = dimk*dimn
     !number of entries
     mumps_par%NZ = dimmat
     !allocate( mumps_par%IRN ( mumps_par%NZ ) )
     !allocate( mumps_par%JCN ( mumps_par%NZ ) )
     !allocate( mumps_par%A( mumps_par%NZ ) )
     !allocate( mumps_par%RHS ( mumps_par%N  ) )
     !DO I = 1, mumps_par%NZ
     !   mumps_par%IRN(I) = cooi(i)
     !   mumps_par%JCN(I) = cooj(i)
     !   mumps_par%A(I)   = mat_coo(i)
     !END DO
     !DO I = 1, mumps_par%N
     !   mumps_par%RHS(I) = vectbb(i)
     !END DO
     mumps_par%IRN=>cooi
     mumps_par%JCN=>cooj
     mumps_par%A=>mat_coo
     mumps_par%RHS=>vectbb
  END IF
  if (.false.) call mexPrintf("Hello B2")
  if (verbose <= 1) then
     mumps_par%ICNTL(4)=1
     mumps_par%ICNTL(3)=0
     !only because we have a problem with mexfile
     !mumps_par%ICNTL(1)=0
  endif

  if (.false.) call mexPrintf("Hello B2.5")
  !  do i = 1,40 
  !    print*,"icntl",i,mumps_par%ICNTL(i)
  !  enddo
  !  Call package for solution
  mumps_par%JOB = 6
  CALL DMUMPS(mumps_par)
  if (.false.) call mexPrintf("Hello B2.6")
  !  Solution has been assembled on the host
  IF ( mumps_par%MYID .eq. 0 .and. verbose > 1) THEN
     WRITE( 6, * ) ' MUMPS Solution is ',(mumps_par%RHS(I),I=1,mumps_par%N)
  END IF
  !vectbb=mumps_par%RHS
  !  Deallocate user data
  IF ( mumps_par%MYID .eq. 0 )THEN
     !DEALLOCATE( mumps_par%IRN )
     !DEALLOCATE( mumps_par%JCN )
     !DEALLOCATE( mumps_par%A   )
     !DEALLOCATE( mumps_par%RHS )
  END IF
  !  Destroy the instance (deallocate internal data structures)
  mumps_par%JOB = -2
  if (.false.) call mexPrintf("Hello B3")
  CALL DMUMPS(mumps_par)
  if (.false.) call mexPrintf("Hello B4")
  CALL MPI_FINALIZE(IERR)
end subroutine solvepde1d

integer function condlimord2_k1(interpretatif)

  !***********************
  !boundary conditions k=1
  !***********************

  use TransportData
  implicit none

  integer,intent(in) :: interpretatif
  integer :: dirichlet=0,neumann=0,m,i
  real(kind=8) :: aux

  condlimord2_k1=0

  matu = T0(:,:,3)/(dx*dx) + T0(:,:,2) / (2.0*dx)
  matv = T0(:,:,1)         - 2.0*T0(:,:,3) / (dx*dx)
  matw = T0(:,:,3)/(dx*dx) - T0(:,:,2) / (2.0*dx)

  do m=1,dimm
     if (mode(m)==1) cycle
     if (abs(matw(1,m))>1e-16) then
        if (verbose>1) print*,"cond1 k=1 m=",m
        neumann=1
        if (interpretatif==1.or.dirichlet==1) then
           if (verbose>1) print*,"boundary conditions not compatible"
           condlimord2_k1 = 1
        endif
     else         
        if (abs(matv(1,m))>1e-16) then
           if (verbose>1) print*,"cond2 k=1 m=",m
           dirichlet=1
           if (neumann==1) then
              if (verbose>1) print*,"boundary conditions not compatible"
              condlimord2_k1 = 1
           endif
        else
           if (verbose>1) print*,"boundary conditions not supported"
           condlimord2_k1 = 1
        endif
     endif

     if (abs(matw(2,m))>1e-16) then
        if (verbose>1) print*,"cond1p k=1 m=",m
        neumann=1
        if (interpretatif==1.or.dirichlet==1) then
           if (verbose>1) print*,"boundary conditions not compatible"
           condlimord2_k1 = 1
        endif
     else 
        if (abs(matv(2,m))>1e-16) then
           if (verbose>1) print*,"cond2p k=1 m=",m
           dirichlet=1
           if (neumann==1) then
              if (verbose>1) print*,"boundary conditions not compatible"
              condlimord2_k1 = 1
           endif
        else
           if (verbose>1) print*,"boundary conditions not supported"
           condlimord2_k1 =1
        endif
     endif
  enddo

  if (interpretatif==1) condlimord2_k1=1  

  if (condlimord2_k1 .ne. 0) return

  do m=1,dimm
     if (mode(m)==1) cycle
     if (abs(matw(1,m))>1e-16) then
        matu(1,m)   = - matu(1,m)/matw(1,m)
        matv(1,m)   = - matv(1,m)/matw(1,m)
        matw(1,m)   =   V0(1,m)  /matw(1,m)
        do i=1,dimm
           matb(1,m,i) =   matb(1,m,i) + matv(1,m)*mata(1,m,i)
           matc(1,m,i) =   matc(1,m,i) + matu(1,m)*mata(1,m,i)
           mats(1,m) = mats(1,m) + sca_f*matw(1,m)*mata(1,m,i)
        enddo
        !mata(1,m,:) = 0.0
     else         
        matx(1,m)   = -matu(1,m)/matv(1,m) 
        maty(1,m)   = V0(1,m) /matv(1,m)
        matb(1,:,m) = 0.0
        matc(1,:,m) = 0.0
        mats(1,m)   = 0.0
        do i=1,dimm
           matb(2,m,i) = matb(2,m,i) + matx(1,m)*mata(2,m,i)
           mats(2,i) =  mats(2,i) + sca_f*maty(1,m)*mata(2,m,i)
        enddo
        mata(2,m,:) = 0.0
     endif

     if (abs(matw(2,m))>1e-16) then
        matu(2,m)    = - matu(2,m)/matw(2,m)
        matv(2,m)    = - matv(2,m)/matw(2,m)
        matw(2,m)    =   V0(2,m)  /matw(2,m)
        do i=1,dimm
           matbp(1,m,i) =   matbp(1,m,i) + matv(2,m)*matap(1,m,i)
           matcp(1,m,i) =   matcp(1,m,i) + matu(2,m)*matap(1,m,i)
           mats(1,i) = mats(1,i) + (1.0-sca_f)*matw(2,m)*matap(1,m,i)
        enddo
        !matap(1,m,:) = 0.0
     else 
        matx(2,m)    = -matu(2,m)/matv(2,m) 
        maty(2,m)    = V0(2,m) /matv(2,m)
        matbp(1,:,m) = 0.0
        matcp(1,:,m) = 0.0 
        matcp(1,m,m) = matx(2,m)/(1.0 - sca_f)
        mats(1,m)    = maty(2,m)
        do i=1,dimm
           matbp(2,m,i) = matbp(2,m,i) + matx(2,m)*matap(2,m,i)
           mats(2,i) =  mats(2,i) + (1.0-sca_f)*maty(2,m)*matap(2,m,i)
        enddo
        matap(2,m,:) = 0.0
     endif
  enddo

end function condlimord2_k1

integer function condlimord2_kk(interpretatif)

  !***********************
  !boundary conditions k=dimk
  !***********************

  use TransportData
  implicit none

  integer,intent(in) :: interpretatif
  integer :: dirichlet=0,neumann=0,m,i
  real(kind=8) :: aux

  condlimord2_kk=0

  matu = T1(:,:,3)/(dx*dx) + T1(:,:,2) / (2.0*dx)
  matv = T1(:,:,1)         - 2.0*T1(:,:,3) / (dx*dx)
  matw = T1(:,:,3)/(dx*dx) - T1(:,:,2) / (2.0*dx)

  do m=1,dimm
     if (mode(m)==1) cycle
     if (abs(matu(1,m))>1e-16) then
        if (verbose>1) print*,"cond1 k=dimk m=",m
        neumann=1
        if (interpretatif==1.or.dirichlet==1) then
           if (verbose>1) print*,"boundary conditions not compatible"
           condlimord2_kk=1
        endif
     else 
        if (abs(matv(1,m))>1e-16) then
           if (verbose>1) print*,"cond2 k=dimk m=",m
           dirichlet=1
           if (neumann==1) then
              if (verbose>1) print*,"boundary conditions not compatible"
              condlimord2_kk=1
           endif
        else
           if (verbose>1) print*,"bondary conditions not supported"
           condlimord2_kk=1
        endif
     endif
     if (abs(matu(2,m))>1e-16) then
        if (verbose>1) print*,"cond1p k=dimk m=",m
        neumann=1
        if (interpretatif==1.or.dirichlet==1) then
           if (verbose>1) print*,"boundary conditions not compatible"
           condlimord2_kk=1
        endif
     else 
        if (abs(matv(2,m))>1e-16) then
           if (verbose>1) print*,"cond2p k=dimk m=",m
           if (neumann==1) then
              if (verbose>1) print*,"boundary conditions not compatible"
              condlimord2_kk=1
           endif
        else
           if (verbose>1) print*,"boundary conditions not supported"
           condlimord2_kk=1
        endif
     endif
  enddo
  if (interpretatif==1) condlimord2_kk=1

  if (condlimord2_kk.ne.0) return

  do m=1,dimm
     if (mode(m)==1) cycle
     if (abs(matu(1,m))>1e-16) then
        aux=matu(1,m)
        matu(1,m)   = - matw(1,m)/aux           
        matv(1,m)   = - matv(1,m)/aux
        matw(1,m)   =   V1(1,m)  /aux
        do i=1,dimm
           matb(dimk,m,i) =   matb(dimk,m,i) + matv(1,m)*matc(dimk,m,i)
           mata(dimk,m,i) =   mata(dimk,m,i) + matu(1,m)*matc(dimk,m,i)
           mats(dimk,i) = mats(dimk,i) + sca_f*matw(1,m)*matc(dimk,m,i)
        enddo
        !matc(dimk,m,:) = 0.0
     else 
        matx(1,m)   = -matw(1,m)/matv(1,m) 
        maty(1,m)   = V1(1,m) /matv(1,m)           
        mata(dimk,:,m)   = 0.0
        matb(dimk,:,m)   = 0.0
        mats(dimk,m)     = 0.0
        do i=1,dimm
           matb(dimk-1,m,i) = matb(dimk-1,m,i) + matx(1,m)*matc(dimk-1,m,i)
           mats(dimk-1,i) =  mats(dimk-1,i) +  sca_f*maty(1,m)*matc(dimk-1,m,i)
        enddo
        matc(dimk-1,m,:) = 0.0
     endif
     if (abs(matu(2,m))>1e-16) then
        aux=matu(2,m)
        matu(2,m)   = - matw(2,m)/aux           
        matv(2,m)   = - matv(2,m)/aux
        matw(2,m)   =   V1(2,m)  /aux
        do i=1,dimm
           matbp(dimk,m,i) =   matbp(dimk,m,i) + matv(2,m)*matcp(dimk,m,i)
           matap(dimk,m,i) =   matap(dimk,m,i) + matu(2,m)*matcp(dimk,m,i)
           mats(dimk,i) = mats(dimk,i) + (1.0-sca_f)*matw(2,m)*matcp(dimk,m,i)
        enddo
        !matcp(dimk,m,:) = 0.0
     else 
        matx(2,m)   = -matw(2,m)/matv(2,m) 
        maty(2,m)   = V1(2,m) /matv(2,m)
        matap(dimk,:,m)   = 0.0
        matap(dimk,m,m)   = matx(2,m)/(1.0 - sca_f)
        matbp(dimk,:,m)   = 0.0
        mats(dimk,m)      = maty(2,m)
        do i=1,dimm
           matbp(dimk-1,m,i) = matbp(dimk-1,m,i) + matx(2,m)*matcp(dimk-1,m,i)
           mats(dimk-1,i) =  mats(dimk-1,i) +  (1.0-sca_f)*maty(2,m)*matcp(dimk-1,m,i)
        enddo
        matcp(dimk-1,m,:) = 0.0
     endif
  enddo

end function condlimord2_kk

subroutine condlimord1_k1
  !***********************
  !boundary conditions k=1
  !***********************

  use TransportData
  implicit none

  integer :: m,i
  real(kind=8) :: aux

  matu = T0(:,:,2) 
  matv = T0(:,:,1)*dx - T0(:,:,2) 
  maty = V0(:,:)*dx
 
  do m=1,dimm
     if (mode(m)==1) cycle
     matx(1,m)   = -matu(1,m)/matv(1,m)
     maty(1,m)   =  maty(1,m)/matv(1,m)
     matb(1,:,m) = 0.0
     matc(1,:,m) = 0.0
     mats(1,m)   = 0.0    
     do i=1,dimm
        matb(2,m,i) = matb(2,m,i) + matx(1,m)*mata(2,m,i)
        mats(2,i) =  mats(2,i) + sca_f*maty(1,m)*mata(2,m,i)
     enddo
     mata(2,m,:) = 0.0    
     if (abs(matv(2,m))>1e-16) then
        matx(2,m)    = -matu(2,m)/matv(2,m)
        maty(2,m)   =  maty(2,m)/matv(2,m)
        matbp(1,:,m) = 0.0
        matcp(1,:,m) = 0.0
        matcp(1,m,m) = matx(2,m)/(1.0 - sca_f)
        mats(1,m)    = maty(2,m)
        do i=1,dimm
           matbp(2,m,i) = matbp(2,m,i) + matx(2,m)*matap(2,m,i)
           mats(2,i) =  mats(2,i) + (1.0-sca_f)*maty(2,m)*matap(2,m,i)
        enddo
        matap(2,m,:) = 0.0
     else
        print*,"Boundaries conditions error"
     endif        
  enddo
end subroutine condlimord1_k1

subroutine condlimord1_kk
  !**************************
  !boundary conditions k=dimk
  !**************************

  use TransportData
  implicit none

  integer :: m,i
  real(kind=8) :: aux

  matv = T1(:,:,2) + dx*T1(:,:,1)
  matw = - T1(:,:,2)
  maty = V1(:,:)*dx

  do m=1,dimm
     if (mode(m)==1) cycle
     matx(1,m)   = -matw(1,m)/matv(1,m)
     maty(1,m)   =  maty(1,m)/matv(1,m)           
     mata(dimk,:,m)   = 0.0
     matb(dimk,:,m)   = 0.0
     mats(dimk,m)     = 0.0
     do i=1,dimm
        matb(dimk-1,m,i) = matb(dimk-1,m,i) + matx(1,m)*matc(dimk-1,m,i)
        mats(dimk-1,i) =  mats(dimk-1,i) +  sca_f*maty(1,m)*matc(dimk-1,m,i)
     enddo
     matc(dimk-1,m,:) = 0.0
     if (abs(matv(2,m))>1e-16) then
        matx(2,m)   = -matw(2,m)/matv(2,m) 
        maty(2,m)   =  maty(2,m)/matv(2,m)
        matap(dimk,:,m)   = 0.0
        matap(dimk,m,m)   = matx(2,m)/(1.0 - sca_f)
        matbp(dimk,:,m)   = 0.0
        mats(dimk,m)      = maty(2,m)
        do i=1,dimm
           matbp(dimk-1,m,i) = matbp(dimk-1,m,i) + matx(2,m)*matcp(dimk-1,m,i)
           mats(dimk-1,i) =  mats(dimk-1,i) +  (1.0-sca_f)*maty(2,m)*matcp(dimk-1,m,i)
        enddo
        matcp(dimk-1,m,:) = 0.0
     else
        print*,"Boundaries conditions error"
     endif
  enddo
  
end subroutine condlimord1_kk
