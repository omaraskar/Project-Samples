!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                  3D BOUNDARY CONDITION CELL-CENTERS                 !
!                             Oct 2013                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE BCcellcenter3D(phi,xc,yc,sig,dsig,No_cp,nbe,Hpr,h,etan)
  
!---------------------------------------------------------------------!
!                                                                     !
!    This program calculates the correct cell-center boundary         !
!    condition. We have three choices:                                !
!                                                                     !      
!    ChooseBoundary = 0:  Combination of type of boundary. It is      !
!                         determintated by nbe & nbev.                !
!    ChooseBoundary = 1:  All boundaries of the problem are Dirichlet !
!    ChooseBoundary = 2:  All boundaries of the proble are Neumann    !
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!    Output variables:                                                !
!   _______________________________________________________________   !
!  |     Name    |    Size   | Description                         |  !  
!  |_____________|___________|_____________________________________|  !
!  | <-- phi     |(N_CELL,NZ)| Function at the cell-center         |  !   
!  |_____________|___________|_____________________________________|  !
!                                                                     !
!    Input variables:                                                 !
!   _______________________________________________________________   !
!  |     Name    |    Size   | Description                         |  !  
!  |_____________|___________|_____________________________________|  !  
!  | --> xc,yc   |(N_CELL)   | Coordinates of the cell centers     |  !
!  | --> sig     |(NZ)       | sigma value at the cell centers     |  !
!  | --> dsig    |(NZ)       | = sig(k+1)-sig(k+1)                 |  !
!  | --> No_cp   |(N_CELL,3) | Numbering of surrounding 3 cell-cent|  !
!  | --> nbe     |(N_CELL0)  | Tag type cell-center                |  !
!  | --> Hpr     |(N_CELL)   | Total depth = h + etan  cell-center |  !
!  | --> h       |(N_CELL)   | still depth             cell-center |  !
!  | --> etan    |(N_CELL)   | free surface            cell-center |  !
!  |_____________|___________|_____________________________________|  !
!                                                                     !
!    Local variables:                                                 !
!   _______________________________________________________________   !
!  |   Name      |                 Description                     |  !  
!  |_____________|_________________________________________________|  !   
!  | xee,yee     |  Perpendicular intersection of each edge        |  !
!  |_____________|_________________________________________________|  !
!                                                                     !
!   -->  Input variables                                              !
!   <--  Output variables                                             !
!   <->  Input and output variables                                   !
!   ---  Parameters                                                   !
!        Common variables used                                        !
!    *   Common variables modified                                    !
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
      real*8, dimension(:)  :: Hpr(N_CELL)
      real*8, dimension(:)  :: h(N_CELL)
      real*8, dimension(:)  :: etan(N_CELL)
!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|

      real*8 :: x,y,z,fB,dfBdn
      real*8 :: z1,z2,z3,f1,f2,f3,a1,a2,a3
      real*8 :: nnx,nny,nnz
      real*8 :: funSolExam3D,Neumanndfdn3D
      real*8 :: funExamNSp,funExamNSDp,NeumanndpdnNS
      integer :: elem,ii

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
           print*, '      ----> Begin subroutine: BCcellcenter3D'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!*********************************************************************!
!                                                                     !
!                        Navier-Stokes Equation                       !
!                                                                     !
!*********************************************************************!

      IF (RUNTestNSEqn.eq.1) THEN

!      ________________________________________________________
!     |                                                        |
!     |           N-S: Dirichlet Boundary Condition            |
!     |________________________________________________________|

      if (ChooseBoundary.eq.1) then
!        ________________________________________________________
!        Vertical
         do i=1,N_CELL0
            x = xc(i)
            y = yc(i)
!           ______                    
!           Bottom
            k = 1
            z1 = sig(1)
            z2 = sig(2)
            z3 = sig(3)
            z  = 0.5d0*(z1+z2)
            fB = funExamNSp(x,y,z,time) 
            f2 = phi(i,2)
            f3 = phi(i,3)
            a1 = (z-z2)*(z-z3)/((z1-z2)*(z1-z3))
            a2 = (z-z1)*(z-z3)/((z2-z1)*(z2-z3))
            a3 = (z-z1)*(z-z2)/((z3-z1)*(z3-z2))
            f1 = (fB-a2*f2-a3*f3)/a1
            !phi(i,1) = f1
            phi(i,1) = 2.0d0*fB-phi(i,2)
!           ______                   
!           Top
            k = NZ
            z1 = sig(NZ-2)
            z2 = sig(NZ-1)
            z3 = sig(NZ)
            z  = 0.5d0*(z2+z3)
            f1 = phi(i,NZ-2)
            f2 = phi(i,NZ-1)
            fB = funExamNSp(x,y,z,time)
            a1 = (z-z2)*(z-z3)/((z1-z2)*(z1-z3))
            a2 = (z-z1)*(z-z3)/((z2-z1)*(z2-z3))
            a3 = (z-z1)*(z-z2)/((z3-z1)*(z3-z2))
            f3 = (fB-a1*f1-a2*f2)/a3
            !phi(i,NZ) = f3
            phi(i,k)= 2.0d0*fB-phi(i,k-1)
         enddo
!        ________________________________________________________
!        Horizontal
         do i=1,N_CELL0
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
                        z = sig(k)*Hpr(i)-h(i)
                        fB = funExamNSp(x,y,z,time)          
                        phi(nc,k) = 2.0d0*fB-phi(i,k)                    
                     enddo
                  endif 
	       enddo
            endif
         enddo
      endif
!      ________________________________________________________
!     |                                                        |
!     |              N-S: Neumann BC approximation             |
!     |________________________________________________________|

      if (ChooseBoundary.eq.2) then
!        ________________________________________________________
!        Vertical
         do i=1,N_CELL0
            x = xc(i)
            y = yc(i)
!           ______                    
!           Bottom
            k = 1
            z = sig(k)+0.5d0*dsig(k)
            nnx =  0.0d0
            nny =  0.0d0
            nnz = -1.0d0
            dfBdn = NeumanndpdnNS(x,y,z,nnx,nny,nnz,time)
            phi(i,k) = dsig(k)*dfBdn + phi(i,k+1) 
!           ______                   
!           Top
            k = NZ
            z = sig(k)-0.5d0*dsig(k)
            nnx = 0.0d0
            nny = 0.0d0
            nnz = 1.0d0
	    dfBdn = NeumanndpdnNS(x,y,z,nnx,nny,nnz,time)
            phi(i,k) = dsig(k-1)*dfBdn + phi(i,k-1)
         enddo
!        ________________________________________________________
!        Horizontal
         do i=1,N_CELL0
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
                     nnx = normxc(i,j)
                     nny = normyc(i,j)
                     nnz = 0.0d0
                     do k=1,NZ 
                        x = xe(i,j)
                        y = ye(i,j)
                        z = sig(k)
                        dfBdn = NeumanndpdnNS(x,y,z,nnx,nny,nnz)
                        phi(nc,k) = 2.0d0*dlCE(i,j)*dfBdn + phi(i,k)
                     enddo
                  endif 
	       enddo
            endif
         enddo
      endif

!*********************************************************************!
!                                                                     !
!                            Other problems                           !
!                                                                     !
!*********************************************************************!
      ELSE
!      ________________________________________________________
!     |                                                        |
!     |     Dirichlet Boundary Condition (NEW FORMULATION)     |
!     |________________________________________________________|

      if (ChooseBoundary.eq.1) then
!        ________________________________________________________
!        Vertical
         do i=1,N_CELL0
            x = xc(i)
            y = yc(i)
!           ______                    
!           Bottom
            !k = 1
            !z1 = sig(1)
            !z2 = sig(2)
            !z3 = sig(3)
            z  = 0.5d0*(sig(1)+sig(2))*Hpr(i)-h(i)
            fB = funSolExam3D(x,y,z) 
            !f2 = phi(i,2)
            !f3 = phi(i,3)
            !a1 = (z-z2)*(z-z3)/((z1-z2)*(z1-z3))
            !a2 = (z-z1)*(z-z3)/((z2-z1)*(z2-z3))
            !a3 = (z-z1)*(z-z2)/((z3-z1)*(z3-z2))
            !f1 = (fB-a2*f2-a3*f3)/a1
            !phi(i,1) = f1
            phi(i,1) = 2.0d0*fB-phi(i,2)
!           ______                   
!           Top
            !k = NZ
            !z1 = sig(NZ-2)
            !z2 = sig(NZ-1)
            !z3 = sig(NZ)
            z  = 0.5d0*(sig(NZ-1)+sig(NZ))*Hpr(i)-h(i)
            !f1 = phi(i,NZ-2)
            !f2 = phi(i,NZ-1)
            fB = funSolExam3D(x,y,z)
            !a1 = (z-z2)*(z-z3)/((z1-z2)*(z1-z3))
            !a2 = (z-z1)*(z-z3)/((z2-z1)*(z2-z3))
            !a3 = (z-z1)*(z-z2)/((z3-z1)*(z3-z2))
            !f3 = (fB-a1*f1-a2*f2)/a3
            !phi(i,NZ) = f3
            phi(i,NZ)= 2.0d0*fB-phi(i,NZ-1)
         enddo
!        ________________________________________________________
!        Horizontal       
          do ii=N_CELL0+1,N_CELLexact
	     i = No_cp(ii,1)
	     j = No_cp(ii,2)
             do k=1,NZ 
                x = xe(i,j)
                y = ye(i,j)
                z = sig(k)*Hpr(i)-h(i)
                fB = funSolExam3D(x,y,z)   
                phi(ii,k) = 2.0d0*fB-phi(i,k)                    
             enddo
          enddo
      endif
          
!      ________________________________________________________
!     |                                                        |
!     |     Dirichlet Boundary Condition (OLD FORMULATION)     |
!     |________________________________________________________|

      if (ChooseBoundary.eq.1000) then
!        ________________________________________________________
!        Vertical

         do i=1,N_CELL0
            x = xc(i)
            y = yc(i)
!           ______                    
!           Bottom
            k = 1
            z1 = sig(1)
            z2 = sig(2)
            z3 = sig(3)
            z  = 0.5d0*(sig(1)+sig(2))*Hpr(i)-h(i)
            fB = funSolExam3D(x,y,z) 
            f2 = phi(i,2)
            f3 = phi(i,3)
            a1 = (z-z2)*(z-z3)/((z1-z2)*(z1-z3))
            a2 = (z-z1)*(z-z3)/((z2-z1)*(z2-z3))
            a3 = (z-z1)*(z-z2)/((z3-z1)*(z3-z2))
            f1 = (fB-a2*f2-a3*f3)/a1
            !phi(i,1) = f1
            phi(i,1) = 2.0d0*fB-phi(i,2)
            !phi(i,1)  = funSolExam3D(x,y,sig(1)*Hpr(i)-h(i))
!           ______                   
!           Top
            k = NZ
            z1 = sig(NZ-2)
            z2 = sig(NZ-1)
            z3 = sig(NZ)
            z  = 0.5d0*(sig(NZ-1)+sig(NZ))*Hpr(i)-h(i)
            f1 = phi(i,NZ-2)
            f2 = phi(i,NZ-1)
            fB = funSolExam3D(x,y,z)
            a1 = (z-z2)*(z-z3)/((z1-z2)*(z1-z3))
            a2 = (z-z1)*(z-z3)/((z2-z1)*(z2-z3))
            a3 = (z-z1)*(z-z2)/((z3-z1)*(z3-z2))
            f3 = (fB-a1*f1-a2*f2)/a3
            !phi(i,NZ) = f3
            phi(i,NZ)= 2.0d0*fB-phi(i,k-1)
            !phi(i,NZ)  = funSolExam3D(x,y,sig(NZ)*Hpr(i)-h(i))
         enddo
!        ________________________________________________________
!        Horizontal
         ii = N_CELL0
         do i=1,N_CELL0
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
                        z = sig(k)*Hpr(i)-h(i)
                        fB = funSolExam3D(x,y,z)   
                        phi(nc,k) = 2.0d0*fB-phi(i,k)                    
                     enddo
                     ii = ii + 1
                  endif 
	       enddo
            endif
         enddo
      endif          
!      ________________________________________________________
!     |                                                        |
!     |                 Neumann BC approximation               |
!     |________________________________________________________|

      if (ChooseBoundary.eq.2) then
!        ________________________________________________________
!        Vertical
         do i=1,N_CELL0
            x = xc(i)
            y = yc(i)
!           ______                    
!           Bottom
            k = 1
            z = sig(k)+0.5d0*dsig(k)
            nnx =  0.0d0
            nny =  0.0d0
            nnz = -1.0d0
            dfBdn = Neumanndfdn3D(x,y,z,nnx,nny,nnz)
            phi(i,k) = dsig(k)*dfBdn + phi(i,k+1) 
!           ______                   
!           Top
            k = NZ
            z = sig(k)-0.5d0*dsig(k)
            nnx = 0.0d0
            nny = 0.0d0
            nnz = 1.0d0
	    dfBdn = Neumanndfdn3D(x,y,z,nnx,nny,nnz)
            phi(i,k) = dsig(k-1)*dfBdn + phi(i,k-1)
         enddo
!        ________________________________________________________
!        Horizontal
         do i=1,N_CELL0
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
                     nnx = normxc(i,j)
                     nny = normyc(i,j)
                     nnz = 0.0d0
                     do k=1,NZ 
                        x = xe(i,j)
                        y = ye(i,j)
                        z = sig(k)
                        dfBdn = Neumanndfdn3D(x,y,z,nnx,nny,nnz)
                        phi(nc,k) = 2.0d0*dlCE(i,j)*dfBdn + phi(i,k)
                     enddo
                  endif 
	       enddo
            endif
         enddo
      endif
 
!*********************************************************************!
!                                                                     !
!                           Finalization                              !
!                                                                     !
!*********************************************************************!

      ENDIF

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
           print*, '      <----   End subroutine: BCcellcenter3D'
           print*,' '
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	         
      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                    END OF CELL-CENTER BOUNDARY CONDITION            !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
