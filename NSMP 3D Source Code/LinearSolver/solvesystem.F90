!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!               SOR Method (Successive-Over-Relaxation)               !
!                            Miguel Uh                                !
!                         September 2012                              !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE solvesystem(xg,Am0,Am1,Am2,Am3,AmT,AmB,bm,&
                             No_cp,nbe)

!---------------------------------------------------------------------!
!                                                                     !
!     This subroutine solves a linear system of the form:             !
!                                                                     !
!                            Am*xg = bm                               !
!                                                                     !
!     using the methods: SOR,  where:                                 !
!                                                                     !                
!           Am    is the matrix                                       !
!           bm    is the right-hand-side                              !
!           xg    is the solution and initial guess                   !
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!    Output variables:                                                !
!   _______________________________________________________________   !
!  |     Name    |    Size    | Description                        |  !  
!  |_____________|____________|____________________________________|  !
!  | <-- xg      |(N_CELL0,NZ)| Solution of the system             |  !  
!  |_____________|____________|____________________________________|  !
!                                                                     !
!    Input variables:                                                 !
!   _______________________________________________________________   !
!  |     Name    |    Size    | Description                        |  !  
!  |_____________|____________|____________________________________|  !  
!  | --> Am0     |(N_CELL0,NZ)| matrix coeff. horizontal point i   |  !
!  | --> Am1     |(N_CELL0,NZ)| matrix coeff. horizontal neigborn 1|  !
!  | --> Am2     |(N_CELL0,NZ)| matrix coeff. horizontal neigborn 2|  !
!  | --> Am3     |(N_CELL0,NZ)| matrix coeff. horizontal neigborn 3|  !
!  | --> AmT     |(N_CELL0,NZ)| matrix coeff. vertical top         |  ! 
!  | --> AmB     |(N_CELL0,NZ)| matrix coeff. vertical bottom      |  ! 
!  | --> bm      |(N_CELL0,NZ)| right-hand side                    |  !  
!  |_____________|____________|____________________________________|  !
!                                                                     !
!    Common parameters used:                                          !
!   _______________________________________________________________   !
!  |   Name       |                 Description                    |  !  
!  |______________|________________________________________________|  ! 
!  |--- N_CELL0   | Number of inside cell centers                  |  !
!  |--- N_CELL    | Total number of the cell centers               |  !
!  |    NZ        | Number of points in the vertical direction     |  !  
!  |______________|________________________________________________|  !
!                                                                     !
!    Common variables modified:                                       !
!   _______________________________________________________________   !
!  |   Name       |                 Description                    |  !  
!  |______________|________________________________________________|  !   
!  | * i,k        | Loop counters: vertices,cells, others          |  !
!  |______________|________________________________________________|  !
!                                                                     !
!    Local variables:                                                 !
!   _______________________________________________________________   !
!  |   Name      |                 Description                     |  !  
!  |_____________|_________________________________________________|  !   
!  | b0,b1,b2,b3 |  parallel auxiliar variables                    |  !
!  |_____________|_________________________________________________|  !
!  | jc1,jc2,jc3 |  cell-center index                              |  !
!  | iter        |  Current iteration                              |  !
!  | errorsys    |  current error of the iteration                 |  !
!  | errormax    |  shared error in the parallel option            |  !
!  | solF        |  auxiliar solution during iteration             |  !
!  |_____________|_________________________________________________|  !  
!  | residu      |  residual [b-(L+U+D)x]/D  in the SOR method     |  !
!  | som         |  = (L+U)x  in the SOR method                    |  !
!  | relax       |  relaxation parameter in the SOR method         |  !
!  |_____________|_________________________________________________|  !
!                                                                     !
!   -->  Input variables                                              !
!   <--  Output variables                                             !
!   ---  Parameters                                                   !
!        Common variables used                                        !
!    *   Common variables modified                                    !
!---------------------------------------------------------------------!

!*********************************************************************!
!                                                                     !
!                           Definitions                               !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |   Keys and common parameters                           |
!     |________________________________________________________|

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
!      ________________________________________________________
!     |                                                        |
!     |    Declaration of variables                            |
!     |________________________________________________________|

      real*8, dimension(:,:) ::  xg(N_CELL0,NZ)
      real*8, dimension(:,:) ::  Am0(N_CELL0,NZ)
      real*8, dimension(:,:) ::  Am1(N_CELL0,NZ)
      real*8, dimension(:,:) ::  Am2(N_CELL0,NZ)
      real*8, dimension(:,:) ::  Am3(N_CELL0,NZ)
      real*8, dimension(:,:) ::  AmT(N_CELL0,NZ)
      real*8, dimension(:,:) ::  AmB(N_CELL0,NZ)
      real*8, dimension(:,:) ::  bm(N_CELL0,NZ)
      integer,dimension(:,:) ::  No_cp(N_CELL,3)
      integer,dimension(:)   ::  nbe(N_CELL0)  
!      ________________________________________________________
!     |                                                        |
!     |    Declaration of local variables                      |
!     |________________________________________________________|

      real*8, dimension(:,:),allocatable :: solxg
      real*8, dimension(:,:),allocatable :: A0,A1,A2,A3,AT,AB
      real*8, parameter:: relax=0.9d0
      integer,parameter:: ChooseMethodSys = 1
      real*8  :: errorsys,errormax
      real*8  :: residu,som
      integer :: iter
      integer :: jc1,jc2,jc3

!     ~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '----> Begin subroutine: resolsystem'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~

!*********************************************************************!
!                                                                     !
!                          Initialization                             !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |                        Allocate                        |
!     |________________________________________________________|

      allocate(solxg(N_CELL0,NZ), &
               A0(N_CELL0,NZ),    &
               A1(N_CELL0,NZ),    &
               A2(N_CELL0,NZ),    &
               A3(N_CELL0,NZ),    &
               AT(N_CELL0,NZ),    &
               AB(N_CELL0,NZ))
!      ________________________________________________________
!     |                                                        |
!     |                   Initial guess  x0                    |
!     |________________________________________________________|

      do k=1,NZ
         do i=1,N_CELL0
            solxg(i,k) = xg(i,k)
         enddo
      enddo
!      ________________________________________________________
!     |                                                        |
!     |              Parallel sharing information              |
!     |________________________________________________________|

!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel3D
	  do k=1,NZ
             do i=1,N_CELL0
                A0(i,k)=Am0(i,k)
                A1(i,k)=Am1(i,k)
                A2(i,k)=Am2(i,k)
	        A3(i,k)=Am3(i,k)
	        AT(i,k)=AmT(i,k)
	        AB(i,k)=AmB(i,k)
             enddo
	  enddo

 	  call communication_variable(A0)
 	  call communication_variable(A1)
 	  call communication_variable(A2)
 	  call communication_variable(A3)
 	  call communication_variable(AT)
 	  call communication_variable(AB)

	  do k=1,NZ
             do i=1,N_CELL0
                Am0(i,k) = A0(i,k)
                Am1(i,k) = A1(i,k)
                Am2(i,k) = A2(i,k)
	        Am3(i,k) = A3(i,k)
	        AmT(i,k) = AT(i,k)
	        AmB(i,k) = AB(i,k)
             enddo
	  enddo

#     endif
!     =============== END ================    
!     ====================================

!      ________________________________________________________
!     |                                                        |
!     |                  Implicit solution                     |
!     |________________________________________________________|

      IF (ChooseMethodSys.eq.1) then

      iter=0
111   continue
      iter=iter+1
!      __________________________________
!     |                                  |
!     |         Solving by SOR           |
!     |__________________________________|  

!     ------------------------------------
!     Bottom boundary
      k=1
      errorsys = 0.0d0
      do i=1,N_CELL0
	 !if (nbe(i).ge.1) then	
            jc1 = No_cp(i,1)
            jc2 = No_cp(i,2)
            jc3 = No_cp(i,3)	       
	    som =  Am1(i,k)*solxg(jc1,k) &
	         + Am2(i,k)*solxg(jc2,k) &
		 + Am3(i,k)*solxg(jc3,k) &
                 + AmT(i,k)*solxg(i,k+1) 
	    residu = (bm(i,k)-som)/Am0(i,k)-solxg(i,k)
	    errorsys = errorsys + abs(residu)
	    solxg(i,k) = solxg(i,k) + relax*residu
         !endif
      enddo	


!     ------------------------------------
!     Inside points
      do k=2,NZ-1
	 do i=1,N_CELL0
	    !if (nbe(i).ge.1) then	
               jc1 = No_cp(i,1)
               jc2 = No_cp(i,2)
               jc3 = No_cp(i,3)	       
	       som =  Am1(i,k)*solxg(jc1,k) &
		    + Am2(i,k)*solxg(jc2,k) &
		    + Am3(i,k)*solxg(jc3,k) &
                    + AmT(i,k)*solxg(i,k+1) &
                    + AmB(i,k)*solxg(i,k-1) 
	       residu = (bm(i,k)-som)/Am0(i,k)-solxg(i,k)
	       errorsys = errorsys + abs(residu)
	       solxg(i,k) = solxg(i,k) + relax*residu
            !endif
        enddo	
      enddo

!     ------------------------------------
!     Top boundary
      k=NZ
      do i=1,N_CELL0
	 !if (nbe(i).ge.1) then	
            jc1 = No_cp(i,1)
            jc2 = No_cp(i,2)
            jc3 = No_cp(i,3)	       
	    som =  Am1(i,k)*solxg(jc1,k) &
	         + Am2(i,k)*solxg(jc2,k) &
		 + Am3(i,k)*solxg(jc3,k) &
                 + AmB(i,k)*solxg(i,k-1) 
	    residu = (bm(i,k)-som)/Am0(i,k)-solxg(i,k)
	    errorsys = errorsys + abs(residu)
	    solxg(i,k) = solxg(i,k) + relax*residu
         !endif
      enddo

!      __________________________________
!     |                                  |
!     |   communicate sol at each iter   |
!     |__________________________________| 

!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel3D
 	    call communication_pression(solxg)
	    call epsp_parallel(errorsys,errormax)
	    errorsys = errormax
#     endif
!     =============== END ================    
!     ====================================

!      __________________________________
!     |                                  |
!     |     Condition of convergence     |
!     |__________________________________| 


      if(iter.gt.MaxIters) then
          write(*,'(t22,a24,i5,a9,e10.3)') 'Non-convergence: iters =',iter,&
                                           ', error =',errorsys
	  stop
      else
     	  if (errorsys.gt.eps) goto 111
          write(*,'(t22,a24,i5,a9,e10.3)') 'Solution system: iters =',iter,&
                                           ', error =',errorsys
      endif

119   continue

      ENDIF 
!      ________________________________________________________
!     |                                                        |
!     |                    Final solution                      |
!     |________________________________________________________|

      do k=1,NZ
         do i = 1,N_CELL0
            xg(i,k) = solxg(i,k)
         enddo
      enddo


!*********************************************************************!
!                                                                     !
!                           Finalization                              !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |                         Deallocate                     |
!     |________________________________________________________|

      deallocate(solxg,A0,A1,A2,A3,AT,AB)

!     ~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '<---- End   subroutine: resolsystem'
         print*,' '
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~
	         
      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!                                                                     !
!                      	     END of resolsytem                        !
!                                                                     !
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

