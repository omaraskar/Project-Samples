!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!           SOLUTION OF THE 3D POISSON EQUATION BY NEW S.O.R.         !
!                             April 2013                              !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE sorNew3D(phi,phiv,                  &
                          rhs,Gamx,Gamy,Gamz,        &
                          xc,yc,sig,dsig,No_cp,nbe,  &
                          xv,yv,sigv,dsigv,No_vp,nbev)                

!---------------------------------------------------------------------!
!                                                                     !
!    This subroutine solves a linear system using the method based    !
!    in the steady-state of a parabolic diffusion equation and using  !
!    the S.O.R. technique for different relaxion factors to solve     !
!    the linear system.                                               !
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!    Variables:                                                       !
!   _______________________________________________________________   !
!  |     Name   |    Size     |  Description                       |  !  
!  |____________|_____________|____________________________________|  !
!  | <--> phi   |(N_CELL,NZ)  | Solution & initial guess           |  !
!  | <--> phiv  |(N_VERT,NZ-1)| Solution at the vertex values      |  !
!  |____________|_____________|____________________________________|  !
!  | ---> rhs   |(N_CELL,NZ)  | Right-hand side of the system      |  !
!  | ---> Gamx  |(N_CELL,NZ)  | Diffusive coefficient in x         |  ! 
!  | ---> Gamy  |(N_CELL,NZ)  | Diffusive coefficient in y         |  ! 
!  | ---> Gamz  |(N_CELL,NZ)  | Diffusive coefficient in z         |  ! 
!  |____________|_____________|____________________________________|  !
!  | ---> xc,yc | (N_CELL)    | Coordinates of the cell centers    |  !
!  | ---> sig   | (NZ)        | sigma value at the cell centers    |  !
!  | ---> dsig  | (NZ)        | = sig(k+1)-sig(k+1)                |  !
!  | ---> No_cp | (N_CELL,3)  | Numbering of surrounding 3 cell-cen|  !
!  | ---> nbe   | (N_CELL0)   | Tag type cell-center               |  !
!  |____________|_____________|____________________________________|  !
!  | ---> xv,yv | (N_VERT)    | Coordinates of the vertices        |  !
!  | ---> sigv  | (NZ-1)      | sigma value at the vertices        |  !
!  | ---> dsigv | (NZ-1)      | = sigv(k+1)-sigv(k)                |  !
!  | ---> No_vp | (N_CELL0,3) | Numbering of the 3 cell vertices   |  !
!  | ---> nbev  | (N_VERT)    | Tag type of cell vertex            |  !
!  |____________|_____________|____________________________________|  !
!                                                                     !
!---------------------------------------------------------------------!

!*********************************************************************!
!                                                                     !
!                           Definitions                               !
!                                                                     !
!*********************************************************************!
!      ____________________________________
!     |                                    |
!     |     Keys and common parameters     |
!     |____________________________________|

#     include "cppdefs.h"
!     ====================================
!     ==========  SEQUENTIAL =============
#     ifndef KeyParallel
            USE geometry
            implicit none
#     endif
!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
            USE parallel
            USE geometry
            implicit none
#     endif
!     =============== END ================    
!     ====================================  
!      ____________________________________
!     |                                    |
!     |      Declaration of variables      |
!     |____________________________________|

      real*8,dimension(:,:) :: phi(N_CELL,NZ)
      real*8,dimension(:,:) :: phiv(N_VERT,NZ-1)
!     -------------------------------------
      real*8,dimension(:,:) :: rhs(N_CELL,NZ)
      real*8,dimension(:,:) :: Gamx(N_CELL,NZ)
      real*8,dimension(:,:) :: Gamy(N_CELL,NZ)
      real*8,dimension(:,:) :: Gamz(N_CELL,NZ)
!     -------------------------------------
      real*8,dimension(:)   :: xc(N_CELL)
      real*8,dimension(:)   :: yc(N_CELL)
      real*8,dimension(:)   :: sig(NZ)
      real*8,dimension(:)   :: dsig(NZ)
      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:)  :: nbe(N_CELL0) 
!     -------------------------------------
      real*8,dimension(:)   :: xv(N_VERT)
      real*8,dimension(:)   :: yv(N_VERT)
      real*8,dimension(:)   :: sigv(NZ-1)
      real*8,dimension(:)   :: dsigv(NZ-1)
      integer,dimension(:,:):: No_vp(N_CELL0,3)
      integer,dimension(:)  :: nbev(N_VERT) 
!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|

      real*8,dimension(:,:) :: Am0(N_CELL0,NZ)
      real*8,dimension(:,:) :: Am1(N_CELL0,NZ)
      real*8,dimension(:,:) :: Am2(N_CELL0,NZ)
      real*8,dimension(:,:) :: Am3(N_CELL0,NZ)
      real*8,dimension(:,:) :: AmT(N_CELL0,NZ)
      real*8,dimension(:,:) :: AmB(N_CELL0,NZ)
      real*8,dimension(:,:) :: Bmv1T(N_CELL0,NZ)
      real*8,dimension(:,:) :: Bmv2T(N_CELL0,NZ)
      real*8,dimension(:,:) :: Bmv3T(N_CELL0,NZ)
      real*8,dimension(:,:) :: Bmv1B(N_CELL0,NZ)
      real*8,dimension(:,:) :: Bmv2B(N_CELL0,NZ)
      real*8,dimension(:,:) :: Bmv3B(N_CELL0,NZ)
      real*8,dimension(:,:) :: Newrhs(N_CELL0,NZ) 
      integer:: jc1,jc2,jc3
      integer:: jv1,jv2,jv3
!     --------------------------------------
      real*8,dimension(:,:) :: phin(N_CELL0,NZ)
      real*8,dimension(:,:) :: MAm0(N_CELL0,NZ)
      real*8,dimension(:,:) :: MAm1(N_CELL0,NZ)
      real*8,dimension(:,:) :: MAm2(N_CELL0,NZ)
      real*8,dimension(:,:) :: MAm3(N_CELL0,NZ)
      real*8,dimension(:,:) :: MAmT(N_CELL0,NZ)
      real*8,dimension(:,:) :: MAmB(N_CELL0,NZ)
      real*8,dimension(:,:) :: Vbm(N_CELL0,NZ)  
!     ----------------------------------------
      real*8  :: ErrorConv,som
      integer :: iterN
!     ----------------------------------------
      integer, parameter :: ChooseSolver = 1
!     ----------------------------------------

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '----> Begin subroutine: NewSOR 3D'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!*********************************************************************!
!                                                                     !
!                            Initialization                           !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |        Boundary Condition of the initial guess         |
!     |________________________________________________________|

!     ______________________________________________________
!     Cell-centers BC   
      call BCcellcenter3D(phi,xc,yc,sig,dsig,No_cp,nbe)

!     ________________________________________________________
!    |                                                        |
!    |            Matrix Am & Bm of the diffusion term        |
!    |________________________________________________________|

      do k=1,NZ
         do i=1,N_CELL0 	
            Am0(i,k)  = 0.0d0 
            Am1(i,k)  = 0.0d0
            Am2(i,k)  = 0.0d0
            Am3(i,k)  = 0.0d0
            AmT(i,k)  = 0.0d0
            AmB(i,k)  = 0.0d0
            Bmv1T(i,k)= 0.0d0 
            Bmv2T(i,k)= 0.0d0 
            Bmv3T(i,k)= 0.0d0 
            Bmv1B(i,k)= 0.0d0 
            Bmv2B(i,k)= 0.0d0 
            Bmv3B(i,k)= 0.0d0 
         enddo
      enddo

      call diffusion3D(Am0,Am1,Am2,Am3,AmT,AmB,             & 
                       Bmv1T,Bmv2T,Bmv3T,Bmv1B,Bmv2B,Bmv3B, &
                       Gamx,Gamy,Gamz,                      &
                       xc,yc,sig,dsig,No_cp,nbe,            &
                       xv,yv,sigv,dsigv,No_vp,nbev)

!*********************************************************************!
!                                                                     !
!                 Solution of the system (S.0.R.) 3D                  !
!                                                                     !
!*********************************************************************!


!      ________________________________________________________
!     |                                                        |
!     |                        New Matrix                      |
!     |________________________________________________________|

      do k=1,NZ
         do i=1,N_CELL0 	
            MAm0(i,k) = ttheta*Am0(i,k)-ttau  
            MAm1(i,k) = ttheta*Am1(i,k)
            MAm2(i,k) = ttheta*Am2(i,k)
            MAm3(i,k) = ttheta*Am3(i,k)
            MAmT(i,k) = ttheta*AmT(i,k)
            MAmB(i,k) = ttheta*AmB(i,k)
         enddo
      enddo
!      ________________________________________________________
!     |                                                        |
!     |             Initial calculations of the loop           |
!     |________________________________________________________|

      iterN=0
11    continue
      iterN=iterN+1

!     ________________________________________________________
!     Save current iteration  

      do k=1,NZ
         do i=1,N_CELL0	
            phin(i,k) = phi(i,k)
         enddo
      enddo
!      ________________________________________________________
!     |                                                        |
!     |                Update vertex values                    |
!     |________________________________________________________|
      
!     ------------------------------------------------
!     Interpolation  
      call interpolation3D(phiv,xv,yv,sigv,dsigv,No_vp,nbev,&
                           phi,xc,yc,sig,dsig,No_cp,nbe)
!     ------------------------------------------------
!     Vertex BC 
      call BCVertex3D(phiv,xv,yv,sigv,dsigv,No_vp,nbev, &
                      phi,xc,yc,sig,dsig,No_cp,nbe)
!      ________________________________________________________
!     |                                                        |
!     |                  New right-hand side                   |
!     |________________________________________________________|

      do k=2,NZ-1
         do i=1,N_CELL0 
!           ------------------------------
!           Cell-center terms
            jc1 = No_cp(i,1)
            jc2 = No_cp(i,2)
            jc3 = No_cp(i,3)	
  	    som =   Am0(i,k)*phi(i,k)   &
                  + Am1(i,k)*phi(jc1,k) &
                  + Am2(i,k)*phi(jc2,k) &
                  + Am3(i,k)*phi(jc3,k) &
                  + AmT(i,k)*phi(i,k+1) &       
                  + AmB(i,k)*phi(i,k-1) 
!           ------------------------------
!           Vertex terms
            jv1 = No_vp(i,1)
            jv2 = No_vp(i,2)
            jv3 = No_vp(i,3)
            Newrhs(i,k) = rhs(i,k)-( Bmv1T(i,k)*phiv(jv1,k)   &
                                    +Bmv2T(i,k)*phiv(jv2,k)   &
                                    +Bmv3T(i,k)*phiv(jv3,k)   &
                                    +Bmv1B(i,k)*phiv(jv1,k-1) &
                                    +Bmv2B(i,k)*phiv(jv2,k-1) &
                                    +Bmv3B(i,k)*phiv(jv3,k-1) )
!           ------------------------------
!           Final rhs at each iteration
            Vbm(i,k)  = -(1-ttheta)*som -ttau*phi(i,k) + Newrhs(i,k)
         enddo
      enddo
!      ________________________________________________________
!     |                                                        |
!     |             Solution of the system (S.0.R.)            |
!     |________________________________________________________|

!     ________________________________________________________
!     SOR
      if (ChooseSolver.eq.1) then 
          call solSOR3D(phi,                                &
                        MAm0,MAm1,MAm2,MAm3,MAmT,MAmB,Vbm,  & 
                        xc,yc,sig,dsig,No_cp,nbe)

!     ________________________________________________________
!     GMRES
      elseif (ChooseSolver.eq.2) then
         call  solGMRES3D(phi,                                &
                          MAm0,MAm1,MAm2,MAm3,MAmT,MAmB,Vbm,  & 
                          xc,yc,sig,dsig,No_cp,nbe)
      endif
!      ________________________________________________________
!     |                                                        |
!     |            Error of each convergence iteration         |
!     |________________________________________________________|

      ErrorConv=0.0d0
      do i=1,N_CELL0 
         do k=2,NZ-1
            ErrorConv = ErrorConv +dabs(phi(i,k)-phin(i,k))**2
         enddo
      enddo
      ErrorConv = dsqrt(ErrorConv)/float(N_CELL0*(NZ-2))   

!      ________________________________________________________
!     |                                                        |
!     |          Convergence criteria of the method            |
!     |________________________________________________________|

8     format(t10,a29,i5,a9,e10.3)
        
      if (ErrorConv.lt.epsConv) then
          write(*,*) ' '
          write(*,8) 'Sol Psedo-time Method: iters =',iterN,&
                     ', error =',ErrorConv
          write(*,*) ' '
      elseif (ErrorConv.gt.1.0d5) then
          write(*,*) ' '
          write(*,8) 'DIVERGENCE !!!!: iters =',iterN,&
                     ', error =',ErrorConv
          write(*,*) ' '
          stop
      elseif(iterN.gt.MaxConv) then
          write(*,*) ' '
          write(*,8) 'Non-convergence: iters =',iterN,&
                     ', error =',ErrorConv
          write(*,*) ' '
      else
     	  goto 11
      endif

19    continue

!*********************************************************************!
!                                                                     !
!                             Finalization                            !
!                                                                     !
!*********************************************************************!

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '<---- End   subroutine: NewSOR 3D'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	         
      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                      	  END OF New S.O.R. 3D                        !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
