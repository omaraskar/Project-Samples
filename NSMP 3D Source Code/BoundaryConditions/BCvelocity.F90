!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                  3D BOUNDARY CONDITION CELL-CENTERS                 !
!                             Nov 2013                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE BCvelcenter3D(phi,xc,yc,sig,dsig,No_cp,nbe,tagBC)
  
!---------------------------------------------------------------------!
!                                                                     !
!    This program calculates the correct cell-center boundary         !
!    condition of the velocity. The tag called "tagBC" is used to     !
!    choose between velocity components:                              !
!                          tagBC = 1 is u                             !
!                          tagBC = 2 is v                             !
!                          tagBC = 3 is w                             !
!---------------------------------------------------------------------!
!                                                                     !
!    Output variables:                                                !
!   _______________________________________________________________   !
!  |     Name    |    Size   | Description                         |  !  
!  |_____________|___________|_____________________________________|  !
!  | <--- phi    |(N_CELL,NZ)| Function at the cell-center         |  !   
!  |_____________|___________|_____________________________________|  !
!                                                                     !
!    Input variables:                                                 !
!   _______________________________________________________________   !
!  |     Name    |    Size   | Description                         |  !  
!  |_____________|___________|_____________________________________|  !  
!  | ---> xc,yc  |(N_CELL)   | Coordinates of the cell centers     |  !
!  | ---> sig    |(NZ)       | sigma value at the cell centers     |  !
!  | ---> dsig   |(NZ)       | = sig(k+1)-sig(k+1)                 |  !
!  | ---> No_cp  |(N_CELL,3) | Numbering of surrounding 3 cell-cent|  !
!  | ---> nbe    |(N_CELL0)  | Tag type cell-center                |  !
!  |_____________|___________|_____________________________________|  !
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

      real*8,dimension(:,:) :: phi(N_CELL,NZ)
      real*8, dimension(:)  :: xc(N_CELL)
      real*8, dimension(:)  :: yc(N_CELL)
      real*8, dimension(:)  :: sig(NZ)
      real*8, dimension(:)  :: dsig(NZ)
      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:)  :: nbe(N_CELL0)
      integer :: tagBC,elem
!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|

      real*8 :: x,y,z,fB
      real*8 :: z1,z2,z3,f1,f2,f3,a1,a2,a3
      real*8 :: funExamNSu,funExamNSv,funExamNSw

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
           print*, '      ----> Begin subroutine: BCvelcenter3D'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!*********************************************************************!
!                                                                     !
!           Dirichlet Boundary Condition for the velocity             !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |                  Velocity component u                  |
!     |________________________________________________________|

      IF (tagBC.eq.1) THEN
         DO i=1,N_CELL0
!           --------------------------------------------------
!           VERTICAL
            x = xc(i)
            y = yc(i)
!           ______                    
!           Bottom
            k = 1
            z = 0.5d0*(sig(1)+sig(2))
            fB = funExamNSu(x,y,z,time) 
            z1 = sig(1)
            z2 = sig(2)
            z3 = sig(3)
            f2 = phi(i,2)
            f3 = phi(i,3)
            a1 = (z-z2)*(z-z3)/((z1-z2)*(z1-z3))
            a2 = (z-z1)*(z-z3)/((z2-z1)*(z2-z3))
            a3 = (z-z1)*(z-z2)/((z3-z1)*(z3-z2))
            f1 = (fB-a2*f2-a3*f3)/a1
            phi(i,1) = f1
            !phi(i,1) = 2.0d0*fB-phi(i,2)
!           ______                   
!           Top
            k = NZ
            z = 0.5d0*(sig(NZ)+sig(NZ-1))
            fB = funExamNSu(x,y,z,time)
            z1 = sig(NZ-2)
            z2 = sig(NZ-1)
            z3 = sig(NZ)
            f1 = phi(i,NZ-2)
            f2 = phi(i,NZ-1)
            a1 = (z-z2)*(z-z3)/((z1-z2)*(z1-z3))
            a2 = (z-z1)*(z-z3)/((z2-z1)*(z2-z3))
            a3 = (z-z1)*(z-z2)/((z3-z1)*(z3-z2))
            f3 = (fB-a1*f1-a2*f2)/a3
            phi(i,NZ) = f3
            !phi(i,k)= 2.0d0*fB-phi(i,k-1)
!           --------------------------------------------------
!           HORIZONTAL
	    if (nbe(i).ne.0) then	
   	       do j=1,3
	          nc=No_cp(i,j)
!                 ====================================
!                 ==========  SEQUENTIAL =============
#                 ifndef KeyParallel
	          if (nc.lt.1.OR.nc.gt.N_CELL0) then
!                 ====================================
!                 =====  START PARALLEL OPTION =======
#                 else
                  elem = index_global(nc)
	          if (elem.lt.1.OR.elem.gt.N_CELL0global) then
#                 endif	
!                 =============== END ================    
!                 ====================================     
                     do k=1,NZ 
                        x = xe(i,j)
                        y = ye(i,j)
                        z = sig(k)
                        fB = funExamNSu(x,y,z,time)   
                        phi(nc,k) = 2.0d0*fB-phi(i,k)                    
                     enddo
                  endif 
	       enddo
            endif
         ENDDO
!      ________________________________________________________
!     |                                                        |
!     |                  Velocity component v                  |
!     |________________________________________________________|

      ELSEIF (tagBC.eq.2) THEN
         DO i=1,N_CELL0
!           --------------------------------------------------
!           VERTICAL
            x = xc(i)
            y = yc(i)
!           ______                    
!           Bottom
            k = 1
            z = 0.5d0*(sig(1)+sig(2))
            fB = funExamNSv(x,y,z,time) 
            z1 = sig(1)
            z2 = sig(2)
            z3 = sig(3)
            f2 = phi(i,2)
            f3 = phi(i,3)
            a1 = (z-z2)*(z-z3)/((z1-z2)*(z1-z3))
            a2 = (z-z1)*(z-z3)/((z2-z1)*(z2-z3))
            a3 = (z-z1)*(z-z2)/((z3-z1)*(z3-z2))
            f1 = (fB-a2*f2-a3*f3)/a1
            phi(i,1) = f1
            !phi(i,1) = 2.0d0*fB-phi(i,2)
!           ______                   
!           Top
            k = NZ
            z = 0.5d0*(sig(NZ)+sig(NZ-1))
            fB = funExamNSv(x,y,z,time)
            z1 = sig(NZ-2)
            z2 = sig(NZ-1)
            z3 = sig(NZ)
            f1 = phi(i,NZ-2)
            f2 = phi(i,NZ-1)
            a1 = (z-z2)*(z-z3)/((z1-z2)*(z1-z3))
            a2 = (z-z1)*(z-z3)/((z2-z1)*(z2-z3))
            a3 = (z-z1)*(z-z2)/((z3-z1)*(z3-z2))
            f3 = (fB-a1*f1-a2*f2)/a3
            phi(i,NZ) = f3
            !phi(i,k)= 2.0d0*fB-phi(i,k-1)
!           --------------------------------------------------
!           HORIZONTAL
	    if (nbe(i).ne.0) then	
   	       do j=1,3
	          nc=No_cp(i,j)
!                 ====================================
!                 ==========  SEQUENTIAL =============
#                 ifndef KeyParallel
	          if (nc.lt.1.OR.nc.gt.N_CELL0) then
!                 ====================================
!                 =====  START PARALLEL OPTION =======
#                 else
                  elem = index_global(nc)
	          if (elem.lt.1.OR.elem.gt.N_CELL0global) then
#                 endif	
!                 =============== END ================    
!                 ====================================             
                     do k=1,NZ 
                        x = xe(i,j)
                        y = ye(i,j)
                        z = sig(k)
                        fB = funExamNSv(x,y,z,time)   
                        phi(nc,k) = 2.0d0*fB-phi(i,k)                    
                     enddo
                  endif 
	       enddo
            endif
         ENDDO
!      ________________________________________________________
!     |                                                        |
!     |                  Velocity component w                  |
!     |________________________________________________________|

      ELSEIF (tagBC.eq.3) THEN
         DO i=1,N_CELL0
!           --------------------------------------------------
!           VERTICAL
            x = xc(i)
            y = yc(i)
!           ______                    
!           Bottom
            k = 1
            z = 0.5d0*(sig(1)+sig(2))
            fB = funExamNSw(x,y,z,time) 
            z1 = sig(1)
            z2 = sig(2)
            z3 = sig(3)
            f2 = phi(i,2)
            f3 = phi(i,3)
            a1 = (z-z2)*(z-z3)/((z1-z2)*(z1-z3))
            a2 = (z-z1)*(z-z3)/((z2-z1)*(z2-z3))
            a3 = (z-z1)*(z-z2)/((z3-z1)*(z3-z2))
            f1 = (fB-a2*f2-a3*f3)/a1
            phi(i,1) = f1
            phi(i,1) = 2.0d0*fB-phi(i,2)
!           ______                   
!           Top
            k = NZ
            z = 0.5d0*(sig(NZ)+sig(NZ-1))
            fB = funExamNSw(x,y,z,time)
            z1 = sig(NZ-2)
            z2 = sig(NZ-1)
            z3 = sig(NZ)
            f1 = phi(i,NZ-2)
            f2 = phi(i,NZ-1)
            a1 = (z-z2)*(z-z3)/((z1-z2)*(z1-z3))
            a2 = (z-z1)*(z-z3)/((z2-z1)*(z2-z3))
            a3 = (z-z1)*(z-z2)/((z3-z1)*(z3-z2))
            f3 = (fB-a1*f1-a2*f2)/a3
            phi(i,NZ) = f3
            phi(i,k)= 2.0d0*fB-phi(i,k-1)
!           --------------------------------------------------
!           HORIZONTAL
	    if (nbe(i).ne.0) then	
   	       do j=1,3
	          nc=No_cp(i,j)
!                 ====================================
!                 ==========  SEQUENTIAL =============
#                 ifndef KeyParallel
	          if (nc.lt.1.OR.nc.gt.N_CELL0) then
!                 ====================================
!                 =====  START PARALLEL OPTION =======
#                 else
                  elem = index_global(nc)
	          if (elem.lt.1.OR.elem.gt.N_CELL0global) then
#                 endif	
!                 =============== END ================    
!                 ====================================             
                     do k=1,NZ 
                        x = xe(i,j)
                        y = ye(i,j)
                        z = sig(k)
                        fB = funExamNSw(x,y,z,time)   
                        phi(nc,k) = 2.0d0*fB-phi(i,k)                    
                     enddo
                  endif 
	       enddo
            endif
         ENDDO
      ENDIF
 
!*********************************************************************!
!                                                                     !
!                           Finalization                              !
!                                                                     !
!*********************************************************************!

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
           print*, '      <----   End subroutine: BCvelcenter3D'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	         
      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                  VELOCITY VERTEX BOUNDARY CONDITION                 !
!                             Nov 2013                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE BCvelvertex3D(phiv,xv,yv,sigv,dsigv,No_vp,nbev,tagBC)

!---------------------------------------------------------------------!
!                                                                     !
!    This program calculates the correct vertex boundary condition    !
!    of the velocity components. The tag called "tagBC" is used to    !
!    choose between velocity components:                              !
!                          tagBC = 1 is u                             !
!                          tagBC = 2 is v                             !
!                          tagBC = 3 is w                             !
!---------------------------------------------------------------------!
!                                                                     !
!    Output variables:                                                !
!   _______________________________________________________________   !
!  |     Name    |     Size    | Description                       |  !  
!  |_____________|_____________|___________________________________|  !
!  | <--- phiv   |(N_VERT,NZ-1)| Function at the vertices          |  !
!  |_____________|_________________________________________________|  !
!                                                                     !
!    Input variables:                                                 !
!   _______________________________________________________________   !
!  |     Name    |    Size   | Description                         |  !  
!  |_____________|___________|_____________________________________|  !
!  | ---> xv,yv  |(N_VERT)   | Coordinates of the vertices         |  !
!  | ---> sigv   |(NZ-1)     | sigma value at the vertices         |  !
!  | ---> dsigv  |(NZ-1)     | = sigv(k+1)-sigv(k)                 |  !
!  | ---> No_vp  |(N_VERT,3) | Numbering of the 3 cell vertices    |  !
!  | ---> nbev   |(N_VERT)   | Tag type of cell vertex             |  !
!  |_____________|___________|_____________________________________|  !  
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

      real*8,dimension(:,:) :: phiv(N_VERT,NZ-1)
      real*8, dimension(:)  :: xv(N_VERT)
      real*8, dimension(:)  :: yv(N_VERT)
      real*8, dimension(:)  :: sigv(NZ-1)
      real*8, dimension(:)  :: dsigv(NZ-1)
      integer,dimension(:,:):: No_vp(N_CELL0,3)
      integer,dimension(:)  :: nbev(N_VERT)  
      integer :: tagBC
!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|

      real*8 :: x,y,z,fB
      real*8 :: funExamNSu,funExamNSv,funExamNSw

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
           print*, '      ----> Begin subroutine: BCvelvertex3D'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!*********************************************************************!
!                                                                     !
!              Dirichlet Boundary Condition for the velocity          !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |                  Velocity component u                  |
!     |________________________________________________________|

      IF (tagBC.eq.1) THEN
         DO nv=1,N_VERT
!           --------------------------------------------------
!           VERTICAL
            x = xv(nv)
            y = yv(nv)
!           _________                 
!           Bottom
            k = 1
            z  = sigv(k)
            phiv(nv,k) = funExamNSu(x,y,z,time) 
!           _________                
!           Top
            k = NZ-1
            z = sigv(k)  
            phiv(nv,k) = funExamNSu(x,y,z,time)
!           --------------------------------------------------
!           HORIZONTAL
            if (nbev(nv).ne.0) then
               do k=1,NZ-1 
                  x = xv(nv)
                  y = yv(nv)
                  z = sigv(k)
                  phiv(nv,k) = funExamNSu(x,y,z,time) 
               enddo
           endif
         ENDDO
!      ________________________________________________________
!     |                                                        |
!     |                  Velocity component v                  |
!     |________________________________________________________|

      ELSEIF (tagBC.eq.2) THEN
         DO nv=1,N_VERT
!           --------------------------------------------------
!           VERTICAL
            x = xv(nv)
            y = yv(nv)
!           _________                 
!           Bottom
            k = 1
            z  = sigv(k)
            phiv(nv,k) = funExamNSv(x,y,z,time) 
!           _________                
!           Top
            k = NZ-1
            z = sigv(k)  
            phiv(nv,k) = funExamNSv(x,y,z,time)
!           --------------------------------------------------
!           HORIZONTAL
            if (nbev(nv).ne.0) then
               do k=1,NZ-1 
                  x = xv(nv)
                  y = yv(nv)
                  z = sigv(k)
                  phiv(nv,k) = funExamNSv(x,y,z,time) 
               enddo
           endif
         ENDDO
!      ________________________________________________________
!     |                                                        |
!     |                  Velocity component w                  |
!     |________________________________________________________|

      ELSEIF (tagBC.eq.3) THEN
         DO nv=1,N_VERT
!           --------------------------------------------------
!           VERTICAL
            x = xv(nv)
            y = yv(nv)
!           _________                 
!           Bottom
            k = 1
            z  = sigv(k)
            phiv(nv,k) = funExamNSw(x,y,z,time) 
!           _________                
!           Top
            k = NZ-1
            z = sigv(k)  
            phiv(nv,k) = funExamNSw(x,y,z,time)
!           --------------------------------------------------
!           HORIZONTAL
            if (nbev(nv).ne.0) then
               do k=1,NZ-1 
                  x = xv(nv)
                  y = yv(nv)
                  z = sigv(k)
                  phiv(nv,k) = funExamNSw(x,y,z,time) 
               enddo
           endif
         ENDDO
      ENDIF

!*********************************************************************!
!                                                                     !
!                           Finalization                              !
!                                                                     !
!*********************************************************************!

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
           print*, '      <----   End subroutine: BCvelvertex3D'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	         
      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                END OF VELOCITY BOUNDARY CONDITIONS                  !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

!---------------------------------------------------------------------!
!                                                                     !
!                           OLD  BCvelocity                           !
!---------------------------------------------------------------------!

      SUBROUTINE BCvelocity(phi,phiv,                        &        
                            xc,yc,sig,dsig,No_cp,nbe,        &
                            xv,yv,sigv,dsigv,No_vp,nbev,     &
                            tagU)      

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

      real*8,dimension(:,:) :: phi(N_CELL,NZ) 
      real*8,dimension(:,:) :: phiv(N_VERT,NZ-1)
!     ----------------------------------------
      real*8,dimension(:)   :: xc(N_CELL)
      real*8,dimension(:)   :: yc(N_CELL)
      real*8,dimension(:)   :: sig(NZ)
      real*8,dimension(:)   :: dsig(NZ)
      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:)  :: nbe(N_CELL0) 
!     ----------------------------------------
      real*8,dimension(:)   :: xv(N_VERT)
      real*8,dimension(:)   :: yv(N_VERT)
      real*8,dimension(:)   :: sigv(NZ-1)
      real*8,dimension(:)   :: dsigv(NZ-1)
      integer,dimension(:,:):: No_vp(N_CELL0,3)
      integer,dimension(:)  :: nbev(N_VERT) 
!     ----------------------------------------
      integer:: tagU     

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t22,60a)'), '----> Begin subroutine: BC velocity'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!     old . . .

      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                      	      END OF BC                               !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
