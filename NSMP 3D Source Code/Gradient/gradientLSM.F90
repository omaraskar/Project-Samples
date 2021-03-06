!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                   INTERPOLATION OF MAIN VARIABLES                   !
!                             Jan 2013                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE grandientLSM(dfundx,dfundy,dfundsig,     &
                              fun,No_cp,nbe,sig)     
 
!---------------------------------------------------------------------!
!                                                                     !
!    This program approximates the gradient of a function fun using   !
!    the least square technique (Appendix A).                         !
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!    Output variables:                                                !
!   _______________________________________________________________   !
!  |     Name      |    Size   | Description                       |  !  
!  |_______________|___________|___________________________________|  !
!  | <-- dfundx    |(N_CELL,NZ)| d(fun)/dx   = gradient component x|  !
!  | <-- dfundy    |(N_CELL,NZ)| d(fun)/dy   = gradient component y|  ! 
!  | <-- dfundsig  |(N_CELL,NZ)| d(fun)/dsig = gradient component z|  !   
!  |_______________|___________|___________________________________|  !
!                                                                     !
!    Input variables:                                                 !
!   _______________________________________________________________   !
!  |     Name    |    Size   | Description                         |  !  
!  |_____________|___________|_____________________________________|  !  
!  | --> fun     |(N_CELL,NZ)| Function "fun" at the center element|  !
!  | --> nbe     | N_CELL    | Type of boundary cell (inside or bc)|  !
!  | --> No_cp   |(N_CELL,3) | Node No. of surrounding three cell  |  !
!  |_____________|___________|_____________________________________|  !
!                                                                     !
!    Common parameters used:                                          !
!   _______________________________________________________________   !
!  |   Name     |                   Description                    |  !  
!  |____________|__________________________________________________|  ! 
!  |--- N_CELL  |  Total number of the cells                       |  !
!  |--- N_CELL0 |  Number of the cell centers inside the domain    |  !
!  |    NZ      |  Number of vertical points                       |  ! 
!  |____________|__________________________________________________|  ! 
!  |  dxCC      |(N_CELL0,3)| = dx from cell center to cell neigb. |  !
!  |  dyCC      |(N_CELL0,3)| = dy from cell center to cell neigb. |  !
!  |  sum_xc2   | N_CELL0   | Sum(xc(i)-xc(neighborn))^2           |  !
!  |  sum_yc2   | N_CELL0   | Sum(yc(i)-yc(neighborn))^2           |  !
!  |  sum_xcyx  | N_CELL0   | Sum(xc(i)-xc(neig))*(yc(i)-yc(neig)) |  !
!  |____________|___________|______________________________________|  !
!                                                                     !
!    Common variables modified:                                       !
!   _______________________________________________________________   !
!  |   Name       |                 Description                    |  !  
!  |______________|________________________________________________|  !   
!  | * nv,nc,i,j,k| Loop counters: vertices,cells, other           |  !
!  |______________|________________________________________________|  !
!                                                                     !
!    Local parameters & variables:                                    !
!   _______________________________________________________________   !
!  |   Name      |                 Description                     |  !  
!  |_____________|_________________________________________________|  !   
!  | deter       |  determinant of the horizontal sums             |  !
!  | sumfx,sumfy |  real variables to add gradient components      |  !
!  | sumsigsig   |  sum of the sigma terms                         |  !
!  |_____________|_________________________________________________|  !   
!  | f0,f1,f2    |  Function at the last 3 points from the boundary|  !
!  | h1,h2       |  length dsig of last 2 boundary points          |  !
!  | deno        |  denominator of the finite difference formula   |  !
!  |_____________|_________________________________________________|  !   
!  | dfBdn       |  Value of the derivative at the boundary        |  !
!  |_____________|_________________________________________________|  !
!                                                                     !
!   -->  Input variables                                              !
!   <--  Output variables                                             !
!   ---  Parameters                                                   !
!    -   Common variables used                                        !
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

      real*8 ,dimension(:,:) :: dfundx(N_CELL,NZ)
      real*8 ,dimension(:,:) :: dfundy(N_CELL,NZ)
      real*8 ,dimension(:,:) :: dfundsig(N_CELL,NZ)
      real*8 ,dimension(:,:) :: fun(N_CELL,NZ)
      integer,dimension(:,:) :: No_cp(N_CELL,3)
      integer,dimension(:)   :: nbe(N_CELL0) 
      real*8, dimension(:)   :: sig(NZ)
!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|

      real*8 :: sumfx,sumfy,deter
      real*8 :: sumsigsig,dzp,dzm
      real*8 :: f0,f1,f2,h1,h2,deno
      real*8 :: dfBdn
      integer:: jc,jc2,elem

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t22,60a)'), '----> Begin subroutine: GradientLSM '
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!*********************************************************************!
!                                                                     !
!                         Horizontal gradient                         !
!                                                                     !
!*********************************************************************!

      do k=1,NZ
         do i=1,N_CELL	
            dfundx(i,k)  = 0.0d0
            dfundx(i,k)  = 0.0d0
            dfundsig(i,k)= 0.0d0
         enddo
      enddo

      DO k=1,NZ
         do i=1,N_CELL0	
!            __________________________________________________
!           |                                                  |
!           |          Inside & boundary cell-centers          |
!           |__________________________________________________|

	     sumfx = 0.0d0
	     sumfy = 0.0d0
	     do j=1,3
	        jc = No_cp(i,j)                
	        sumfx = sumfx + dxCC(i,j)*(fun(jc,k)-fun(i,k))
	        sumfy = sumfy + dyCC(i,j)*(fun(jc,k)-fun(i,k))
               
	     enddo
	     deter = sum_xc2(i)*sum_yc2(i)-sum_xcyc(i)*sum_xcyc(i)

             dfundx(i,k)=(sum_yc2(i)*sumfx-sum_xcyc(i)*sumfy)/deter
	     dfundy(i,k)=(sum_xc2(i)*sumfy-sum_xcyc(i)*sumfx)/deter
!            __________________________________________________
!           |                                                  |
!           |         Outside cell-centers  (ficticious)       |
!           |__________________________________________________|

!           ----------------------------------
!           Wall 
	    if (nbe(i).eq.1) then  
	       do j=1,3
	          nc=No_cp(i,j)
!                  ====================================
!                  ==========  SEQUENTIAL =============
#                  ifndef KeyParallel
	           if (nc.lt.1.OR.nc.gt.N_CELL0) then
!                  ====================================
!                  =====  START PARALLEL OPTION =======
#                  else
                   elem = index_global(nc)
	           if (elem.lt.1.OR.elem.gt.N_CELL0global) then
#                  endif	
!                  =============== END ================    
!                  ====================================
                     dfBdn = 0.0d0
		     dfundx(nc,k) = 2.0d0*dfBdn+dfundx(i,k) 
		     dfundy(nc,k) = 2.0d0*dfBdn+dfundy(i,k)	
	          endif
	       enddo 
	    endif  

!           ----------------------------------
!           Discharge normal 
	    if (nbe(i).eq.2) then  
	       do j=1,3
	          nc=No_cp(i,j)
!                  ====================================
!                  ==========  SEQUENTIAL =============
#                  ifndef KeyParallel
	           if (nc.lt.1.OR.nc.gt.N_CELL0) then
!                  ====================================
!                  =====  START PARALLEL OPTION =======
#                  else
                   elem = index_global(nc)
	           if (elem.lt.1.OR.elem.gt.N_CELL0global) then
#                  endif	
!                  =============== END ================    
!                  ====================================
                     dfBdn = 0.0d0
		     dfundx(nc,k) = 2.0d0*dfBdn+dfundx(i,k)
		     dfundy(nc,k) = 2.0d0*dfBdn+dfundy(i,k)	
	          endif
	       enddo 
	    endif  

!           ----------------------------------
!           Water level 
	    if (nbe(i).eq.3) then  
	       do j=1,3
	          nc=No_cp(i,j)
!                  ====================================
!                  ==========  SEQUENTIAL =============
#                  ifndef KeyParallel
	           if (nc.lt.1.OR.nc.gt.N_CELL0) then
!                  ====================================
!                  =====  START PARALLEL OPTION =======
#                  else
                   elem = index_global(nc)
	           if (elem.lt.1.OR.elem.gt.N_CELL0global) then
#                  endif	
!                  =============== END ================    
!                  ====================================
                     dfBdn = 0.0d0
		     dfundx(nc,k) = 2.0d0*dfBdn+dfundx(i,k)
		     dfundy(nc,k) = 2.0d0*dfBdn+dfundy(i,k)	
	          endif
	       enddo 
	    endif  

!           ----------------------------------
!           Special case: ?
	    if (nbe(i).eq.5) then  
	       dfundx(i,k) = 0.0d0
	       dfundy(i,k) = 0.0d0	
	       do j=1,3
	          nc=No_cp(i,j)
!                  ====================================
!                  ==========  SEQUENTIAL =============
#                  ifndef KeyParallel
	           if (nc.lt.1.OR.nc.gt.N_CELL0) then
!                  ====================================
!                  =====  START PARALLEL OPTION =======
#                  else
                   elem = index_global(nc)
	           if (elem.lt.1.OR.elem.gt.N_CELL0global) then
#                  endif	
!                  =============== END ================    
!                  ====================================
		     dfundx(nc,k) = 0.0d0
		     dfundy(nc,k) = 0.0d0	
	          endif
	       enddo 
	    endif  
         enddo
      ENDDO

!*********************************************************************!
!                                                                     !
!                         Vertical gradient                           !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |                    Inside elements                     |
!     |________________________________________________________|

      DO k=2,NZ-1
         do i=1,N_CELL0	
            dzp = sig(k+1)-sig(k)
            dzm = sig(k-1)-sig(k)
            sumsigsig = dzp**2 + dzm**2
            dfundsig(i,k) = ((fun(i,k+1)-fun(i,k))*dzp  &
                            +(fun(i,k-1)-fun(i,k))*dzm) &
                            /sumsigsig
         enddo
      ENDDO

!      ________________________________________________________
!     |                                                        |
!     |                   Boundary elements                    |
!     |________________________________________________________|

!      -----------------------------------------------
!      Bottom (k=1)          
       k=1      
       do i=1,N_CELL0
          f0 = fun(i,k) 
          f1 = fun(i,k+1)
          f2 = fun(i,k+2)
          h1 = sig(k+1)-sig(k)
          h2 = sig(k+2)-sig(k)
          deno = (h1*h2*h2-h2*h1*h1)
          dfundsig(i,k) = (-(h1**2)*f2+(h2**2)*f1-(h2*h2-h1*h1)*f0)/deno
       enddo

!      -----------------------------------------------
!      Free surface (k=NZ)          
       k=NZ
       do i=1,N_CELL0
          f0 = fun(i,k) 
          f1 = fun(i,k-1)
          f2 = fun(i,k-2)
          h1 = sig(k-1)-sig(k)
          h2 = sig(k-2)-sig(k)
          deno = (h1*h2*h2-h2*h1*h1)
          dfundsig(i,k) = (-(h1**2)*f2+(h2**2)*f1-(h2*h2-h1*h1)*f0)/deno
       enddo

!      ________________________________________________________
!     |                                                        |
!     |               Outside elements (ficticious)            |
!     |________________________________________________________|

      DO k=1,NZ
         do i=1,N_CELL0	
	    if ((nbe(i).ne.0)) then  
	       do j=1,3
	          nc=No_cp(i,j)
	          if (nc.gt.N_CELL0) then
                     dfBdn = 0.0d0
		     dfundsig(nc,k) = 2.0d0*dfBdn-dfundsig(i,k)	
	          endif
	       enddo 
	    endif  
         enddo
      ENDDO
 

!*********************************************************************!
!                                                                     !
!                   Finalization of the subroutine                    !
!                                                                     !
!*********************************************************************!

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t22,60a)'), '<---- End   subroutine: GradientLSM'
         print*,' '
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	         
      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!                                                                     !
!                      	   END OF INTERPOLATION                       !
!                                                                     !
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
