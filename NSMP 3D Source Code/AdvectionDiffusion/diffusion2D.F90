!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                             DIFUSSION                               !
!                             March 2013                              !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE diffusion2D(Am0,Am1,Am2,Am3,Bmv1,Bmv2,Bmv3, &
                             Gamx2D,Gamy2D,                  &
                             xc,yc,No_cp,nbe,                &
                             xv,yv,No_vp,nbev)     
 
!---------------------------------------------------------------------!
!                                                                     !
!    This subroutine approximates the diffusion contribution to the   !
!    general linear system.                                           !
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!    Output & Input variables:                                        !
!   _______________________________________________________________   !
!  |     Name   |    Size    |  Description                        |  !  
!  |____________|____________|_____________________________________|  !
!  | <--> Am0   | (N_CELL0)  | matrix coefficient of element i     |  !
!  | <--> Am1   | (N_CELL0)  | matrix coeff. horizontal neigborn 1 |  ! 
!  | <--> Am2   | (N_CELL0)  | matrix coeff. horizontal neigborn 2 |  ! 
!  | <--> Am3   | (N_CELL0)  | matrix coeff. horizontal neigborn 3 |  ! 
!  | <--> Bmv1  | (N_VERT)   | Matrix coeff. vertex value 1        |  !
!  | <--> Bmv2  | (N_VERT)   | Matrix coeff. vertex value 2        |  !  
!  | <--> Bmv3  | (N_VERT)   | Matrix coeff. vertex value 3        |  !  
!  |____________|____________|_____________________________________|  !
!                                                                     !
!    Input variables:                                                 !
!   _______________________________________________________________   !
!  |     Name    |    Size   | Description                         |  !   
!  |_____________|___________|_____________________________________|  !
!  | --> Gamx2D  |(N_CELL)   | Function diffusive coeff. Gamma_x   |  !
!  | --> Gamy2D  |(N_CELL)   | Function diffusive coeff. Gamma_x   |  !
!  |_____________|___________|_____________________________________|  !
!  | --> phi2D   |(N_CELL)   | Function phi at the cell-center     |  !
!  | --> nbe     |(N_CELL)   | Type of boundary cell (inside or bc)|  !
!  | --> No_cp   |(N_CELL,3) | Numbering of surrounding three cells|  !
!  | --> xc,yc   |(N_CELL)   | Coordinates of the cell centers     |  !
!  |_____________|___________|_____________________________________|  !
!  | --> phi2Dv  |(N_VERT)   | Function phi at the vertex          |  !
!  | --> xv,yv   |(N_VERT)   | Coordinates of the cell vertices    |  !
!  | --> No_vp   |(N_CELL0,3)| Numbering of the cell vertices      |  !
!  | --> nbev    |(N_VERT)   | Tag: Type of vertex (inside or bc)  |  !
!  |_____________|___________|_____________________________________|  !
!                                                                     !
!    Common parameters and variables used:                            !
!   _______________________________________________________________   !
!  |   Name     |                   Description                    |  !  
!  |____________|__________________________________________________|  ! 
!  |--- N_CELL  |  Total number of the cells                       |  !
!  |--- N_CELL0 |  Number of the cell centers inside the domain    |  !
!  |   AreaCell |  Area of the cell                                |  ! 
!  |____________|__________________________________________________|  !
!                                                                     !
!    Local parameters & variables:                                    !
!   _______________________________________________________________   !
!  |   Name        |               Description                     |  !    
!  |_______________|_______________________________________________|  !
!  | aa0,aa(1:3)   | Auxiliar horizontal matrix coefficients       |  !
!  | aaT,aaB       | Auxiliar vertical matrix coefficients         |  !
!  | C0j(1:3)      | Auxiliar mass flux values of the neighbors    |  !
!  | C0pos,C0neg   | Positive and negative values of the mass flux |  !
!  |_______________|_______________________________________________|  !
!  | jc,jj         | neighborn loop index                          |  !
!  |_______________|_______________________________________________|  !
!                                                                     !
!    Subroutines used:                                                !
!   _______________________________________________________________   !
!  |                                                               |  !  
!  |                    -  grandientLSM                            |  !
!  |                    -  massflux                                |  !
!  |_______________________________________________________________|  !
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

      real*8, dimension(:)  :: Am0(N_CELL0)
      real*8, dimension(:)  :: Am1(N_CELL0)
      real*8, dimension(:)  :: Am2(N_CELL0)
      real*8, dimension(:)  :: Am3(N_CELL0)
      real*8, dimension(:)  :: Bmv1(N_CELL0)
      real*8, dimension(:)  :: Bmv2(N_CELL0)
      real*8, dimension(:)  :: Bmv3(N_CELL0)
!     ------------------------------------
      real*8, dimension(:)  :: Gamx2D(N_CELL)
      real*8, dimension(:)  :: Gamy2D(N_CELL)
!     ------------------------------------
      real*8, dimension(:)  :: xc(N_CELL)
      real*8, dimension(:)  :: yc(N_CELL)
      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:)  :: nbe(N_CELL0) 
!     ------------------------------------
      real*8, dimension(:)  :: xv(N_VERT)
      real*8, dimension(:)  :: yv(N_VERT)
      integer,dimension(:,:):: No_vp(N_CELL0,3)
      integer,dimension(:)  :: nbev(N_VERT) 
!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|

      real*8 ,dimension(:) :: AA(1:3)
      real*8 :: nxLeng,nyLeng
      real*8 :: AA0,AAB,AAT
      real*8 :: BB1,BB2,BB3
      real*8 :: rhs0,Gamx_ej,Gamy_ej
      integer:: jc,jj,jc1,jc2,jc3
!     ------------------------------------
      real*8 :: sumax,sumay,sumfx
      real*8 :: sumbx,sumby,sumfy
      real*8 :: AreaReg,fv1,fv2
      integer:: jv1,jv2,jv3

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '----> Begin subroutine: Diffusion2D'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!*********************************************************************!
!                                                                     !
!                         Diffusion contributions                     !
!                                                                     !
!*********************************************************************!

         do i=1,N_CELL0	
!            __________________________________________________
!           |                                                  |
!           |                  Matrix: Am & Bmv                |
!           |__________________________________________________|

            AA0   = 0.0d0
            AA(1) = 0.0d0
            AA(2) = 0.0d0
            AA(3) = 0.0d0
            BB1   = 0.0d0
            BB2   = 0.0d0
            BB3   = 0.0d0

            do j=1,3
               nxLeng =  dyVV(i,j)
               nyLeng = -dxVV(i,j) 
!              --------------------------------
!              Cell-center index
	       jc = No_cp(i,j)
!              --------------------------------
!              Vertex index               
	       jj=j+1
	       if (jj.gt.3) jj=jj-3
	       jv1 = No_vp(i,j)
	       jv2 = No_vp(i,jj)
!              --------------------------------
!              Volume of the region               
               AreaReg = (areaCell(i)+areaCell(jc))/3.0d0
!              --------------------------------
!              Constants x
               sumax = -0.5d0*(yv(jv2)-yv(jv1))/AreaReg
               sumbx =  0.5d0*(yv(jv2)-yv(jv1))/AreaReg
               sumfx = -0.5d0*(yc(jc)-yc(i))/AreaReg
!              --------------------------------
!              Constants y
               sumay =  0.5d0*(xv(jv2)-xv(jv1))/AreaReg
               sumby = -0.5d0*(xv(jv2)-xv(jv1))/AreaReg
               sumfy =  0.5d0*(xc(jc)-xc(i))/AreaReg
!              ---------------------------------
!              Diffusive constants
               Gamx_ej = 0.5d0*(Gamx2D(i)+Gamx2D(jc))
               Gamy_ej = 0.5d0*(Gamy2D(i)+Gamy2D(jc))
               sumax = Gamx_ej*sumax
               sumbx = Gamx_ej*sumbx
               sumay = Gamy_ej*sumay
               sumby = Gamy_ej*sumby
!              ______________________________________________
!              Cell-center contribution
!              ---------------------------------
!              Matrix coeffients i
               AA0  =  sumax*nxLeng &
                     + sumay*nyLeng + AA0  
!              ---------------------------------
!              Matrix coeffients j
               AA(j)=  sumbx*nxLeng &
                     + sumby*nyLeng   
!              ______________________________________________
!              Vertex contribution 
               fv2 =  (Gamx_ej*sumfx*nxLeng+Gamy_ej*sumfy*nyLeng)    
               fv1 = -(Gamx_ej*sumfx*nxLeng+Gamy_ej*sumfy*nyLeng)  
               if (j.eq.1) then
                  BB2 = BB2 + fv2
                  BB1 = BB1 + fv1
               elseif (j.eq.2) then
                  BB3 = BB3 + fv2
                  BB2 = BB2 + fv1
               elseif (j.eq.3) then
                  BB1 = BB1 + fv2
                  BB3 = BB3 + fv1
               endif 
	    enddo
!            __________________________________________________
!           |                                                  |
!           |                  Contributions                   |
!           |__________________________________________________|

!           --------------------------------------
!           Matrix: cell-center coefficients
	    Am0(i) = Am0(i) + AA0    
            Am1(i) = Am1(i) + AA(1)   
            Am2(i) = Am2(i) + AA(2)   
            Am3(i) = Am3(i) + AA(3)
!           --------------------------------------
!           Matrix: vertex coefficients
            Bmv1(i) = Bmv1(i) + BB1
            Bmv2(i) = Bmv2(i) + BB2
            Bmv3(i) = Bmv3(i) + BB3
        enddo

!*********************************************************************!
!                                                                     !
!                   Finalization of the subroutine                    !
!                                                                     !
!*********************************************************************!

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '<---- End   subroutine: Diffusion2D'
         print*,' '
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	         
      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                    DIFUSSION WITH NEUMANN TREATMENT                 !
!                              Nov 2013                               !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE diffusionNeumann2D(Am0,Am1,Am2,Am3,Bmv1,Bmv2,Bmv3,Nm,&
                                    Gamx2D,Gamy2D,                    &
                                    xc,yc,No_cp,nbe,                  &
                                    xv,yv,No_vp,nbev)     
 
!---------------------------------------------------------------------!
!                                                                     !
!    This subroutine approximates the diffusion contribution term.    !
!    We have the cell-center coefficients Am, the vertex coefficients !
!    Bm and the extra term corresponding to the known BC named as Cm. !
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!    Output & Input variables:                                        !
!   _______________________________________________________________   !
!  |     Name   |    Size    |  Description                        |  !  
!  |____________|____________|_____________________________________|  !
!  | <--> Am0   | (N_CELL0)  | matrix coefficient of element i     |  !
!  | <--> Am1   | (N_CELL0)  | matrix coeff. horizontal neigborn 1 |  ! 
!  | <--> Am2   | (N_CELL0)  | matrix coeff. horizontal neigborn 2 |  ! 
!  | <--> Am3   | (N_CELL0)  | matrix coeff. horizontal neigborn 3 |  ! 
!  | <--> Bmv1  | (N_CELL0)  | Matrix coeff. vertex value 1        |  !
!  | <--> Bmv2  | (N_CELL0)  | Matrix coeff. vertex value 2        |  !  
!  | <--> Bmv3  | (N_CELL0)  | Matrix coeff. vertex value 3        |  ! 
!  | <--> Nm    | (N_CELL0)  | Term with Neumann BC                |  !  
!  |____________|____________|_____________________________________|  !
!                                                                     !
!    Input variables:                                                 !
!   _______________________________________________________________   !
!  |     Name    |    Size   | Description                         |  !   
!  |_____________|___________|_____________________________________|  !
!  | --> Gamx2D  |(N_CELL)   | Function diffusive coeff. Gamma_x   |  !
!  | --> Gamy2D  |(N_CELL)   | Function diffusive coeff. Gamma_x   |  !
!  |_____________|___________|_____________________________________|  !
!  | --> phi2D   |(N_CELL)   | Function phi at the cell-center     |  !
!  | --> nbe     |(N_CELL)   | Type of boundary cell (inside or bc)|  !
!  | --> No_cp   |(N_CELL,3) | Numbering of surrounding three cells|  !
!  | --> xc,yc   |(N_CELL)   | Coordinates of the cell centers     |  !
!  |_____________|___________|_____________________________________|  !
!  | --> phi2Dv  |(N_VERT)   | Function phi at the vertex          |  !
!  | --> xv,yv   |(N_VERT)   | Coordinates of the cell vertices    |  !
!  | --> No_vp   |(N_CELL0,3)| Numbering of the cell vertices      |  !
!  | --> nbev    |(N_VERT)   | Tag: Type of vertex (inside or bc)  |  !
!  |_____________|___________|_____________________________________|  !
!                                                                     !
!    Common parameters and variables used:                            !
!   _______________________________________________________________   !
!  |   Name     |                   Description                    |  !  
!  |____________|__________________________________________________|  ! 
!  |--- N_CELL  |  Total number of the cells:                      |  !
!  |--- N_CELL0 |  Number of the cell centers inside the domain    |  !
!  |   AreaCell |  Area of the cell                                |  ! 
!  |____________|__________________________________________________|  !
!                                                                     !
!    Local parameters & variables:                                    !
!   _______________________________________________________________   !
!  |   Name        |               Description                     |  !    
!  |_______________|_______________________________________________|  !
!  | aa0,aa(1:3)   | Auxiliar horizontal matrix coefficients       |  !
!  | aaT,aaB       | Auxiliar vertical matrix coefficients         |  !
!  | C0j(1:3)      | Auxiliar mass flux values of the neighbors    |  !
!  | C0pos,C0neg   | Positive and negative values of the mass flux |  !
!  |_______________|_______________________________________________|  !
!  | jc,jj         | neighborn loop index                          |  !
!  |_______________|_______________________________________________|  !
!                                                                     !
!    Subroutines used:                                                !
!   _______________________________________________________________   !
!  |                                                               |  !  
!  |                    -  grandientLSM                            |  !
!  |                    -  massflux                                |  !
!  |_______________________________________________________________|  !
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

      real*8, dimension(:)  :: Am0(N_CELL0)
      real*8, dimension(:)  :: Am1(N_CELL0)
      real*8, dimension(:)  :: Am2(N_CELL0)
      real*8, dimension(:)  :: Am3(N_CELL0)
      real*8, dimension(:)  :: Bmv1(N_CELL0)
      real*8, dimension(:)  :: Bmv2(N_CELL0)
      real*8, dimension(:)  :: Bmv3(N_CELL0)
      real*8, dimension(:)  :: Nm(N_CELL0)
!     ------------------------------------
      real*8, dimension(:)  :: Gamx2D(N_CELL)
      real*8, dimension(:)  :: Gamy2D(N_CELL)
!     ------------------------------------
      real*8, dimension(:)  :: xc(N_CELL)
      real*8, dimension(:)  :: yc(N_CELL)
      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:)  :: nbe(N_CELL0) 
!     ------------------------------------
      real*8, dimension(:)  :: xv(N_VERT)
      real*8, dimension(:)  :: yv(N_VERT)
      integer,dimension(:,:):: No_vp(N_CELL0,3)
      integer,dimension(:)  :: nbev(N_VERT) 
!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|

      real*8 ,dimension(:) :: AA(1:3)
      real*8 :: nxLeng,nyLeng
      real*8 :: AA0,AAB,AAT
      real*8 :: BB1,BB2,BB3
      real*8 :: rhs0,Gamx_ej,Gamy_ej
      integer:: jc,jj,jc1,jc2,jc3
!     ------------------------------------
      real*8 :: sumax,sumay,sumfx
      real*8 :: sumbx,sumby,sumfy
      real*8 :: AreaReg,fv1,fv2
      integer:: jv1,jv2,jv3
!     ------------------------------------
      real*8  :: x,y,nnx,nny,dfBdn,Neumanndfdn2D

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '----> Begin subroutine: Diffusion2D'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!*********************************************************************!
!                                                                     !
!                         Diffusion contributions                     !
!                                                                     !
!*********************************************************************!

         do i=1,N_CELL0	
!            __________________________________________________
!           |                                                  |
!           |       Matrix coefficients: Am, Bmv and Nm        |
!           |__________________________________________________|

            AA0   = 0.0d0
            AA(1) = 0.0d0
            AA(2) = 0.0d0
            AA(3) = 0.0d0
            BB1   = 0.0d0
            BB2   = 0.0d0
            BB3   = 0.0d0

            do j=1,3
               nxLeng =  dyVV(i,j)
               nyLeng = -dxVV(i,j) 
!              --------------------------------
!              Cell-center index
	       jc = No_cp(i,j)
!              --------------------------------
!              Vertex index               
	       jj=j+1
	       if (jj.gt.3) jj=jj-3
	       jv1 = No_vp(i,j)
	       jv2 = No_vp(i,jj)
!              --------------------------------
!              Volume of the region               
               AreaReg = (areaCell(i)+areaCell(jc))/3.0d0
!              --------------------------------
!              Constants x
               sumax = -0.5d0*(yv(jv2)-yv(jv1))/AreaReg
               sumbx =  0.5d0*(yv(jv2)-yv(jv1))/AreaReg
               sumfx = -0.5d0*(yc(jc)-yc(i))/AreaReg
!              --------------------------------
!              Constants y
               sumay =  0.5d0*(xv(jv2)-xv(jv1))/AreaReg
               sumby = -0.5d0*(xv(jv2)-xv(jv1))/AreaReg
               sumfy =  0.5d0*(xc(jc)-xc(i))/AreaReg
!              ---------------------------------
!              Diffusive constants
               Gamx_ej = 0.5d0*(Gamx2D(i)+Gamx2D(jc))
               Gamy_ej = 0.5d0*(Gamy2D(i)+Gamy2D(jc))
               sumax = Gamx_ej*sumax
               sumbx = Gamx_ej*sumbx
               sumay = Gamy_ej*sumay
               sumby = Gamy_ej*sumby
!              ________________________________________________________
!              BOUNDARY EDGES: Neumann BC 
               IF ((jc.lt.1).or.(jc.gt.N_CELL0)) THEN  
                  nnx = normxc(i,j)
                  nny = normyc(i,j)
                  x = 0.5d0*(xv(jv1)+xv(jv2))
                  y = 0.5d0*(yv(jv1)+yv(jv2))
                  dfBdn = Neumanndfdn2D(x,y,nnx,nny)
!                 ---------------------
                  Nm(i) = Nm(i) + dfBdn*Gamx2D(i)*dlVV(i,j)
!              ________________________________________________________
!              INSIDE EDGES: all the faces
               ELSE
!                 _________________________
!                 Cell-center contribution
!                 ---------------------
!                 Matrix coeffients i
                  AA0  =  sumax*nxLeng &
                        + sumay*nyLeng + AA0  
!                 ---------------------
!                 Matrix coeffients j
                  AA(j)=  sumbx*nxLeng &
                        + sumby*nyLeng   
!                 _________________________
!                 Vertex contribution 
                  fv2 =  (Gamx_ej*sumfx*nxLeng & 
                         +Gamy_ej*sumfy*nyLeng)
                  fv1 = -(Gamx_ej*sumfx*nxLeng &
                         +Gamy_ej*sumfy*nyLeng)  
                  if (j.eq.1) then
                     BB2 = BB2 + fv2
                     BB1 = BB1 + fv1
                  elseif (j.eq.2) then
                     BB3 = BB3 + fv2
                     BB2 = BB2 + fv1
                  elseif (j.eq.3) then
                     BB1 = BB1 + fv2
                     BB3 = BB3 + fv1
                  endif 
               ENDIF 
	    enddo
!            __________________________________________________
!           |                                                  |
!           |                  Contributions                   |
!           |__________________________________________________|

!           --------------------------------------
!           Matrix: cell-center coefficients
	    Am0(i) = Am0(i) + AA0    
            Am1(i) = Am1(i) + AA(1)   
            Am2(i) = Am2(i) + AA(2)   
            Am3(i) = Am3(i) + AA(3)
!           --------------------------------------
!           Matrix: vertex coefficients
            Bmv1(i) = Bmv1(i) + BB1
            Bmv2(i) = Bmv2(i) + BB2
            Bmv3(i) = Bmv3(i) + BB3
        enddo

!*********************************************************************!
!                                                                     !
!                   Finalization of the subroutine                    !
!                                                                     !
!*********************************************************************!

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '<---- End   subroutine: Diffusion2D'
         print*,' '
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	         
      RETURN
      END


!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!                                                                     !
!                      	    END OF DIFFUSION 2D                       !
!                                                                     !
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
