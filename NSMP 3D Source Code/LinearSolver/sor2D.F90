!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!            1) SOLUTION OF THE 2D POISSON EQUATION BY S.O.R.         !
!                             April 2013                              !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE sor2D(phi2D,phiv2D,        &
                       rhs2D,Gamx2D,Gamy2D, & 
                       xc,yc,No_cp,nbe,     &
                       xv,yv,No_vp,nbev)                

!---------------------------------------------------------------------!
!                                                                     !
!    This subroutine solves a linear system using the S.O.R. techni-  !
!    que for different relaxion factors: relax. This system is good   !
!    for the Poisson problem because it consider the cell-center "Am" !
!    and the vertex "Bmv" matrix coefficients.                        !
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!    Variables:                                                       !
!   _______________________________________________________________   !
!  |     Name   |    Size    |  Description                        |  !  
!  |____________|____________|_____________________________________|  !
!  | <--> phi2D | (N_CELL)   | Solution & initial guess            |  !
!  | <--> phi2Dv| (N_VERT)   | Solution at the vertex values       |  !
!  |____________|____________|_____________________________________|  !
!  | ---> rhs2D | (N_CELL)   | Right-hand side of the system       |  ! 
!  | ---> Gamx2D| (N_CELL)   | Diffusive coefficients in x         |  !
!  | ---> Gamy2D| (N_CELL)   | Diffusive coefficients in y         |  !
!  |____________|____________|_____________________________________|  !
!  | ---> xc,yc | (N_CELL)   | Coordinates of the cell centers     |  !
!  | ---> No_cp | (N_CELL,3) | Numbering of surrounding 3 cell-cent|  !
!  | ---> nbe   | (N_CELL0)  | Tag type cell-center                |  !
!  |____________|____________|_____________________________________|  !
!  | ---> xv,yv | (N_VERT)   | Coordinates of the vertices         |  !
!  | ---> No_vp | (N_CELL0,3)| Numbering of the 3 cell vertices    |  !
!  | ---> nbev  | (N_VERT)   | Tag type of cell vertex             |  !
!  |____________|____________|_____________________________________|  !
!                                                                     !
!    Local parameters & variables:                                    !
!   _______________________________________________________________   !
!  |     Name   |    Size    |  Description                        |  !  
!  |____________|____________|_____________________________________|  !
!  |    Am0     | (N_CELL0)  | matrix coefficient of element i     |  !
!  |    Am1     | (N_CELL0)  | matrix coeff. horizontal neigborn 1 |  ! 
!  |    Am2     | (N_CELL0)  | matrix coeff. horizontal neigborn 2 |  ! 
!  |    Am3     | (N_CELL0)  | matrix coeff. horizontal neigborn 3 |  ! 
!  |____________|____________|_____________________________________|  !
!  |    Bmv1    | (N_CELL0)  | matrix coeff. vertex 1              |  ! 
!  |    Bmv2    | (N_CELL0)  | matrix coeff. vertex 2              |  ! 
!  |    Bmv3    | (N_CELL0)  | matrix coeff. vertex 3              |  ! 
!  |____________|____________|_____________________________________|  !  
!  |    Newrhs  | (N_CELL0)  | right hand side of the method       |  !  
!  |____________|____________|_____________________________________|  !
!                                                                     !
!    Subroutines used:                                                !
!   _______________________________________________________________   !
!  |                                                               |  !  
!  |   - diffusion2D                ( diffusion2D.F90 )            |  !
!  |   - interpolation2D            ( interpolation2D.F90 )        |  !
!  |   - BCcellcenter2D             ( BCcellcenter2D.F90 )         |  !
!  |   - BCvertex2D                 ( BCvertex2D.F90 )             |  !
!  |_______________________________________________________________|  !
!                                                                     !
!   --->  Input variables                                             !
!   <---  Output variables                                            !
!                                                                     !
!---------------------------------------------------------------------!
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

      real*8,dimension(:)   :: phi2D(N_CELL)
      real*8,dimension(:)   :: phiv2D(N_VERT)
!     ----------------------------------------
      real*8,dimension(:)   :: rhs2D(N_CELL)
      real*8,dimension(:)   :: Gamx2D(N_CELL)
      real*8,dimension(:)   :: Gamy2D(N_CELL)
!     ----------------------------------------
      real*8,dimension(:)   :: xc(N_CELL)
      real*8,dimension(:)   :: yc(N_CELL)
      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:)  :: nbe(N_CELL0) 
!     ----------------------------------------
      real*8,dimension(:)   :: xv(N_VERT)
      real*8,dimension(:)   :: yv(N_VERT)
      integer,dimension(:,:):: No_vp(N_CELL0,3)
      integer,dimension(:)  :: nbev(N_VERT) 
!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|

      real*8,dimension(:) :: Am0(N_CELL0)
      real*8,dimension(:) :: Am1(N_CELL0)
      real*8,dimension(:) :: Am2(N_CELL0)
      real*8,dimension(:) :: Am3(N_CELL0) 
      real*8,dimension(:) :: Bmv1(N_CELL0)
      real*8,dimension(:) :: Bmv2(N_CELL0)
      real*8,dimension(:) :: Bmv3(N_CELL0)
      real*8,dimension(:) :: Nm(N_CELL0) 
      real*8,dimension(:) :: Newrhs(N_CELL0) 
      integer:: jc1,jc2,jc3
      integer:: jv1,jv2,jv3
!     ----------------------------------------
      real*8  :: errorsys,residu,som
      real*8  :: errorNeum,errorsysOld
      integer :: iter
!     ---------------------------------------- 
      real*8, dimension(:) :: AA(3)
      real*8  :: x,y,dfBdn,Neumanndfdn2D
      real*8  :: nnx,nny
      integer :: jj
!     ----------------------------------------
      real*8  :: relax 
      relax = relaxSOR
!     ---------------------------------------- 

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '----> Begin subroutine: SOR 2D'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!*********************************************************************!
!                                                                     !
!                           Initialization                            !
!                                                                     !
!*********************************************************************!
!      ________________________________________________________
!     |                                                        |
!     |                     Initial guess                      |
!     |________________________________________________________|

!     ____________________________________
!     Cell-center BC
      call BCcellcenter2D(phi2D,xc,yc,No_cp,nbe)

!     _____________________________________
!     Vertex interpolation and vertex BC

      call interpolation2D(phiv2D,xv,yv,No_vp,nbev, &
                           phi2D,xc,yc,No_cp,nbe)

      call BCVertex2D(phiv2D,xv,yv,No_vp,nbev, &
                      phi2D,xc,yc,No_cp,nbe)

!      ________________________________________________________
!     |                                                        |
!     |            Matrix Am & Bm of the diffusion term        |
!     |________________________________________________________|
        
      do i=1,N_CELL0 	
         Am0(i)  = 0.0d0 
         Am1(i)  = 0.0d0
         Am2(i)  = 0.0d0
         Am3(i)  = 0.0d0
         Bmv1(i) = 0.0d0 
         Bmv2(i) = 0.0d0
         Bmv3(i) = 0.0d0
         Nm(i)   = 0.0d0
      enddo            

!     ________________________________________________________      
!     Dirichlet BC
      if (ChooseBoundary.eq.1) then
         call diffusion2D(Am0,Am1,Am2,Am3,Bmv1,Bmv2,Bmv3, &
                          Gamx2D,Gamy2D,                  &
                          xc,yc,No_cp,nbe,                &
                          xv,yv,No_vp,nbev)
!     ________________________________________________________
!     Neumann BC
      elseif (ChooseBoundary.eq.2) then
         call diffusionNeumann2D(Am0,Am1,Am2,Am3,Bmv1,Bmv2,Bmv3,Nm, &
                                 Gamx2D,Gamy2D,                     &
                                 xc,yc,No_cp,nbe,                   &
                                 xv,yv,No_vp,nbev)
         do i=1,N_CELL0
	    rhs2D(i) = rhs2D(i) - Nm(i)
         enddo
      endif

!*********************************************************************!
!                                                                     !
!                 Solution of the system (S.0.R.) 2D                  !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |             Initial calculations of the loop           |
!     |________________________________________________________|

      errorsys = 0.0d0

      iter=0
111   continue
      iter=iter+1 

      errorsysOld = errorsys
!      ________________________________________________________
!     |                                                        |
!     |         New right-hand side: rhs - Bm(vertex)          |
!     |________________________________________________________|

      do i=1,N_CELL0
         jv1 = No_vp(i,1)
         jv2 = No_vp(i,2)
         jv3 = No_vp(i,3)
         Newrhs(i) =  rhs2D(i)-( Bmv1(i)*phiv2D(jv1) &
                                +Bmv2(i)*phiv2D(jv2) &
                                +Bmv3(i)*phiv2D(jv3) ) 
      enddo
!      ________________________________________________________
!     |                                                        |
!     |             Solution of the system (S.0.R.)            |
!     |________________________________________________________|

      errorsys = 0.0d0
      do i=1,N_CELL0
         jc1 = No_cp(i,1)
         jc2 = No_cp(i,2)
         jc3 = No_cp(i,3)
	 som =  Am1(i)*phi2D(jc1)&
              + Am2(i)*phi2D(jc2)&
              + Am3(i)*phi2D(jc3)
	 residu = (Newrhs(i)-som)/Am0(i)-phi2D(i)
	 errorsys = errorsys + abs(residu)
	 phi2D(i) = phi2D(i) + relax*residu
      enddo

!     ________________________________________________________
!     Boundary conditions cell-centers 

      call BCcellcenter2D(phi2D,xc,yc,No_cp,nbe) 

!      ________________________________________________________
!     |                                                        |
!     |                Update vertex values                    |
!     |________________________________________________________|
        
!     ________________________________________________________
!     Interpolation of the inside vertex points         
      call interpolation2D(phiv2D,xv,yv,No_vp,nbev, &
                           phi2D,xc,yc,No_cp,nbe)

!     ________________________________________________________
!     Boundary Conditions of the vertex points 
      call BCVertex2D(phiv2D,xv,yv,No_vp,nbev, &
                      phi2D,xc,yc,No_cp,nbe)
!      ________________________________________________________
!     |                                                        |
!     |            Alternative Error (Neumann stop)            |
!     |________________________________________________________|

      errorNeum = abs(errorsys-errorsysOld)
!      ________________________________________________________
!     |                                                        |
!     |                  Convergence criteria                  |
!     |________________________________________________________|
        
      7 format(t10,a24,i5,a9,e10.3)
      6 format(t10,a14,i4,a9,e10.3)

      if ((errorsys.lt.eps).or.(errorNeum.lt.1d-5*eps)) then
         write(*,*) ' '
         write(*,7) 'Solution SOR 2D: iters =',iter,&
                    ', error =',errorsys
         write(*,*) ' '
      elseif (errorsys.gt.1.0d5) then
         write(*,*) ' '
         write(*,7) 'DIVERGENCE !!!!: iters =',iter,&
                    ', error =',errorsys
         write(*,*) ' '
         stop
      elseif(iter.gt.MaxIters) then
         write(*,*) ' '
         write(*,7) 'Non-convergence: iters =',iter,&
                    ', error =',errorsys
         write(*,*) ' '
      else
         goto 111
      endif

119   continue

!*********************************************************************!
!                                                                     !
!                           Finalization                              !
!                                                                     !
!*********************************************************************!

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '<---- End   subroutine: SOR 2D'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	         
      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!               2)  S.O.R. FOR CELL-CENTER VARIABLES                  !
!                             Oct 2013                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      subroutine solSOR2D(phi2D,                    &
                          MAm0,MAm1,MAm2,MAm3,Vbm,  & 
                          xc,yc,No_cp,nbe)

!---------------------------------------------------------------------!
!                                                                     !
!    This subroutine solves a linear system using the S.O.R. techni-  !
!    que for different relaxion factors: relax. This system is good   !
!    for the linear systems with coefficients "MAm" and righ-hand     !
!    side "Vbm" that only depends on cell-center phi values.          !
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!    Variables:                                                       !
!   _______________________________________________________________   !
!  |     Name   |    Size     |  Description                       |  !  
!  |____________|_____________|____________________________________|  !
!  | <--> phi2D | (N_CELL)    | Solution & initial guess           |  !
!  |____________|_____________|____________________________________|  !
!  | ---> MAm0  | (N_CELL0)   | matrix coefficient of element i    |  !
!  | ---> MAm1  | (N_CELL0)   | matrix coeff. horizontal neigborn 1|  ! 
!  | ---> MAm2  | (N_CELL0)   | matrix coeff. horizontal neigborn 2|  ! 
!  | ---> MAm3  | (N_CELL0)   | matrix coeff. horizontal neigborn 3|  ! 
!  | ---> Vbm   | (N_CELL0)   | right hand side of the method      |  !  
!  |____________|_____________|____________________________________|  !
!  | ---> xc,yc | (N_CELL)    | Coordinates of the cell centers    |  !
!  | ---> No_cp | (N_CELL,3)  | Numbering of surrounding 3 cell-cen|  !
!  | ---> nbe   | (N_CELL0)   | Tag type cell-center               |  !
!  |____________|_____________|____________________________________|  !
!                                                                     !
!    Subroutines used:                                                !
!   _______________________________________________________________   !
!  |                                                               |  !  
!  |   - BCcellcenter2D             ( BCcellcenter2D.F90 )         |  !
!  |_______________________________________________________________|  !
!                                                                     !
!   --->  Input variables                                             !
!   <---  Output variables                                            !
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

      real*8,dimension(:)   :: phi2D(N_CELL)
      real*8,dimension(:)   :: MAm0(N_CELL0)
      real*8,dimension(:)   :: MAm1(N_CELL0)
      real*8,dimension(:)   :: MAm2(N_CELL0)
      real*8,dimension(:)   :: MAm3(N_CELL0) 
      real*8,dimension(:)   :: Vbm(N_CELL0)
      real*8,dimension(:)   :: xc(N_CELL)
      real*8,dimension(:)   :: yc(N_CELL)
      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:)  :: nbe(N_CELL0)  
!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|

      real*8  :: errorsys,residu,som
      integer :: it
      integer:: jc1,jc2,jc3

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '----> Begin subroutine: SolSOR2D'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!*********************************************************************!
!                                                                     !
!                            Initialization                           !
!                                                                     !
!*********************************************************************!

      it=0
111   continue
      it=it+1
!     ________________________________________________________
!     SOR solution

      errorsys = 0.0d0
      do i=1,N_CELL0
         jc1 = No_cp(i,1)
         jc2 = No_cp(i,2)
         jc3 = No_cp(i,3)
	 som =  MAm1(i)*phi2D(jc1)&
              + MAm2(i)*phi2D(jc2)&
              + MAm3(i)*phi2D(jc3)
	 residu = (Vbm(i)-som)/MAm0(i)-phi2D(i)
	 errorsys = errorsys + abs(residu)
	 phi2D(i) = phi2D(i) + relaxSOR*residu
      enddo

!     ------------------------------------------------
!     Boundary conditions of the cell-centers 
      call BCcellcenter2D(phi2D,xc,yc,No_cp,nbe)

!     ________________________________________________________
!     Convergence criteria SOR   

6     format(t10,a26,i5,a9,e10.3)

      if (errorsys.lt.eps) then
         write(*,6) ' Solution S.0.R.: iters =',it,', error =',errorsys
      elseif (errorsys.gt.1.0d5) then
         write(*,6) ' DIVERGENCE !!!!: iters =',it,', error =',errorsys
         stop
      elseif(it.gt.MaxIters) then
         write(*,6) ' Non-convergence: iters =',it,', error =',errorsys
      else
         goto 111
      endif

119   continue

!*********************************************************************!
!                                                                     !
!                           Finalization                              !
!                                                                     !
!*********************************************************************!

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '<---- End   subroutine:SolSOR2D'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                        END OF S.O.R. SUBROUTINES                    !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
