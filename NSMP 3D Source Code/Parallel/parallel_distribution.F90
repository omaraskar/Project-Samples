!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                        INPUT PARALLEL FILE                          !
!                      Miguel Angel Uh Zapata                         !
!                             Jun 2014                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!


      SUBROUTINE parallel_distribution                                       

!---------------------------------------------------------------------!
!                                                                     !
!    This program calculates all the values needed to define the      !
!    parallel topology and comunication between domains. It aslo      !
!    calculates the corresponding cell-centers and vertices to        ! 
!    each subdomain.                                                  !
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!    Output variables:                                                !
!   _______________________________________________________________   !
!  |      Name     |      Size       |        Description          |  !  
!  |_______________|_________________|_____________________________|  !   
!  |               |                 |                             |  !
!  | TotalSubdom   | Integer         | Total number of subdomains  |  !
!  |_______________|_________________|_____________________________|  !
!  |               |                 |                             |  ! 
!  | CCdom         |(N_CELL0global,1)| Cell-centers in each subdo. |  !
!  |               |(N_CELL0global,2)| Tag of the subsomain        |  !
!  | DimCCdom      |(TotalSub,1)     | No. of cells in each subdo. |  !
!  |               |(TotalSub,2)     | Initial index in global vec.|  !
!  |               |(TotalSub,3)     | Final index in global vector|  !
!  | BdyCCdom      |(N_CELL0global,1)| subdomain tag: s            |  !
!  |               |(N_CELL0global,2)| neighbor subsomain tag: sm  |  !
!  |               |(N_CELL0global,3)| Boundary cells in each subdo|  !
!  | DimBdyCCdom   |(TotalSub,1)     | No. of bdy cells in each sub|  !
!  |               |(TotalSub,2)     | Initial index in global vec.|  !
!  |               |(TotalSub,3)     | Final index in global vec.  |  !
!  | ExtraCCdom    |(N_CELL0global,1)| subdomain tag: s            |  !
!  |               |(N_CELL0global,2)| neighbor subsomain tag: sm  |  !
!  |               |(N_CELL0global,3)| overlaping boundary cells   |  !
!  | DimExtraCCdom |(TotalSub,1)     | No. of overla. cells in sub.|  !
!  |               |(TotalSub,2)     | Initial index in global vec.|  !
!  |               |(TotalSub,3)     | Final index in global vector|  !
!  |_______________|_________________|_____________________________|  ! 
!  |               |                 |                             |  !
!  | VVdom         |(N_VERTglobal,1) | Vertices in each subdo.     |  !
!  |               |(N_VERTglobal,2) | Tag of the subsomain        |  !
!  | DimVVdom      |(TotalSub,1)     | No. of vertices in each sub.|  !
!  |               |(TotalSub,2)     | Initial index in global vec.|  !
!  |               |(TotalSub,3)     | Final index in global vector|  !
!  | BdyVVdom      |(N_VERTglobal,1) | subdomain tag: s            |  !
!  |               |(N_VERTglobal,2) | neighbor subsomain tag: sm  |  !
!  |               |(N_VERTglobal,3) | Boundary vertices in subdo.n|  !
!  | DimBdyVVdom   |(TotalSub,1)     | No. of bdy. vertices in sub.|  !
!  |               |(TotalSub,2)     | Initial index in global vec.|  !
!  |               |(TotalSub,3)     | Final index in global vector|  !
!  |_______________|_______________________________________________|  ! 
!                                                                     !
!---------------------------------------------------------------------!

!*********************************************************************!
!                                                                     !
!                           Definitions                               !
!                                                                     !
!*********************************************************************!
!     ____________________________________
!    |                                    |
!    |     Keys and common parameters     |
!    |____________________________________|

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
!     ____________________________________
!    |                                    |
!    |   Declaration of local variables   |
!    |____________________________________|

      integer,dimension(:,:), allocatable :: BdyCCdom
      integer,dimension(:,:), allocatable :: DimBdyCCdom
      integer,dimension(:,:), allocatable :: BdyVVdom
      integer,dimension(:,:), allocatable :: DimBdyVVdom
      integer,dimension(:,:), allocatable :: InterCCdom
      integer,dimension(:,:), allocatable :: DimInterCCdom
      integer,dimension(:,:), allocatable :: BoundCCdom
      integer,dimension(:,:), allocatable :: DimBoundCCdom
!     ------------------------------
      integer,dimension(:,:), allocatable :: Eedges
      integer,dimension(:),   allocatable :: EdgesVec
      integer,dimension(:),   allocatable :: DimEdges
      integer,dimension(:),   allocatable :: Iindex
!     ------------------------------
      integer :: ii,s,s0,i0,j0,k0,m,sm,tag
      integer :: Jini,Jfin,Iini,Ifin
      integer :: Ielem,Jelem,Jvertex,Ivertex,kelem
      integer :: jc,jc1,jc2,jc3
      integer :: jv,jv1,jv2,jv3
      integer :: elem,irec,NNbefore
      real*8  :: xxc,yyc
!     ------------------------------
      character*80 title
      character*36 filen
      character*39 filen2

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
           print*, '      ----> Begin subroutine: parallel_distribution'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!*********************************************************************!
!                                                                     !
!                 Calculate all the information needed                !
!                                                                     !
!*********************************************************************!
!      ________________________________________________________
!     |                                                        |
!     |                   Number of subdomains                 |
!     |________________________________________________________|

      j0 = 1
      do i=2,N_CELL0global
         tag = 0
         do j=1,i-1
             if (TagParallel(i).eq.TagParallel(j)) tag=1
         enddo
         if (tag.eq.0) j0 = j0+1
      enddo
      TotalSubdomains = j0

      if (TotalSubdomains.ne.Nprocs) then
         print *,'                                                     '
         print*,'  ===================================================='
         print*,'                    FATAL ERROR !!!!!!                '
         print*,'      We must have:                                   '
         print*,'           # of subdomains = # of processors          '
         print*,'                                                      '
         print*,'      We have = ',Nprocs,    ' processors'
         print*,'      We have = ',TotalSubdomains,' subdomains'
         print*,'                                                      '
         print*,'  ===================================================='
         print *,'                                                     '
         stop
      endif
!      ________________________________________________________
!     |                                                        |
!     |           Cell-center indexes in each subdomain        |
!     |________________________________________________________|

      allocate(CCdom(1:N_CELL0global,2)) 
      allocate(DimCCdom(1:TotalSubdomains,3)) 

      i0 = 0
      do s=1,TotalSubdomains
         s0 = 0
         DimCCdom(s,2) = i0+1     !<----Initial index (subdomain)
         do i=1,N_CELL0global
            if (TagParallel(i).eq.(s-1)) then
                i0 = i0 + 1
                s0 = s0 + 1
                CCdom(i0,1) = i  
                CCdom(i0,2) = s   
            endif
         enddo
         DimCCdom(s,1) = s0       !<---- Number of components
         DimCCdom(s,3) = i0       !<---- Final index (subdomain)
      enddo

!     _________________________________________________________
!     Free memeory TagParallel

      deallocate(TagParallel)
!      ________________________________________________________
!     |                                                        |
!     |                  Boundary cell-centers                 |
!     |________________________________________________________|

      allocate(BdyCCdom(TotalSubdomains*N_CELL0global,4))
      allocate(DimBdyCCdom(1:TotalSubdomains,3)) 

      j0 = 0
      do s=1,TotalSubdomains           
         s0 = 0
         DimBdyCCdom(s,2) = j0 + 1     !<----Initial index (subdomain)
         Iini = DimCCdom(s,2)
         Ifin = DimCCdom(s,3)
!        --------------------------------------
         do m=1,TotalSubdomains-1         
            sm = mod(s+m,TotalSubdomains)
            if (sm==0) sm = TotalSubdomains
            Jini = DimCCdom(sm,2)
            Jfin = DimCCdom(sm,3)
            do i=Iini,Ifin   
               Ielem = CCdom(i,1)
               do j=Jini,Jfin
                  Jelem = CCdom(j,1)
                  jc1 = No_cp_global(Ielem,1)
                  jc2 = No_cp_global(Ielem,2)
                  jc3 = No_cp_global(Ielem,3)
                  if ((Jelem.eq.jc1).or. & 
                      (Jelem.eq.jc2).or. & 
                      (Jelem.eq.jc3)) then
!                     -------------------------
!                     Avoid repetition
                      tag = 0
                      do ii=DimBdyCCdom(s,2),j0
                         if (Jelem.eq.BdyCCdom(ii,4)) tag = 1
                      enddo
                      if (tag.eq.0) then
                         j0 = j0 + 1
                         s0 = s0 + 1
                         BdyCCdom(j0,1) = s-1
                         BdyCCdom(j0,2) = Ielem
                         BdyCCdom(j0,3) = sm-1
                         BdyCCdom(j0,4) = Jelem
                      endif
!                     -------------------------
                  endif
               enddo
            enddo
         enddo
!        --------------------------------------
         DimBdyCCdom(s,1) = s0       !<---- Number of components   
         DimBdyCCdom(s,3) = j0       !<---- Final index (subdomain)  
      enddo
!      ________________________________________________________
!     |                                                        |
!     |               Vertex indexes in each subdomain         |
!     |________________________________________________________|

      allocate(VVdom(2*N_VERTglobal,2)) 
      allocate(DimVVdom(TotalSubdomains,3)) 

      j0 = 0
      do s=1,TotalSubdomains
         s0 = 0
         DimVVdom(s,2) = j0+1       !<----Initial index (subdomain)
         Iini = DimCCdom(s,2)
         Ifin = DimCCdom(s,3)
!        -----------------------------------
!        First Triangle (three vertices)   
         s0 = s0 + 1 
         j0 = j0 + 1
         elem = CCdom(Iini,1)
         VVdom(j0,1) = No_vp_global(elem,1)
         VVdom(j0,2) = s
         s0 = s0 + 1
         j0 = j0 + 1
         VVdom(j0,1) = No_vp_global(elem,2)
         VVdom(j0,2) = s
         s0 = s0 + 1
         j0 = j0 + 1
         VVdom(j0,1) = No_vp_global(elem,3)
         VVdom(j0,2) = s         
!        -----------------------------------
!        Rest of triangles
         do i=Iini+1,Ifin
            do k=1,3
               elem = CCdom(i,1)
               jv = No_vp_global(elem,k)
!              ----------------
!              Avoid repetition
               tag = 0
               do j=DimVVdom(s,2),j0
                  if (jv.eq.VVdom(j,1)) tag = 1
               enddo
               if (tag.eq.0) then
                  j0 = j0 + 1
                  s0 = s0 + 1
                  VVdom(j0,1) = jv
                  VVdom(j0,2) = s  
               endif
!              ----------------
            enddo
         enddo
!        -----------------------------------
         DimVVdom(s,1) = s0       !<---- Number of components
         DimVVdom(s,3) = j0       !<---- Final index (subdomain)
      enddo
!      ________________________________________________________
!     |                                                        |
!     |                    Boundary vertex                     |
!     |________________________________________________________|

      allocate(BdyVVdom(N_VERTglobal,3)) 
      allocate(DimBdyVVdom(TotalSubdomains,3)) 

      j0 = 0
      do s=1,TotalSubdomains           
         s0 = 0
         DimBdyVVdom(s,2) = j0 + 1!<----Initial index (subdomain)
!        -----------------------------
         do m=1,TotalSubdomains-1         
            sm = mod(s+m,TotalSubdomains)
            if (sm==0) sm = TotalSubdomains
            do i=DimVVdom(s,2),DimVVdom(s,3)   
               Ivertex = VVdom(i,1)
               do j=DimVVdom(sm,2),DimVVdom(sm,3)
                  Jvertex = VVdom(j,1)         
                  if (Jvertex.eq.Ivertex) then
                     j0 = j0 + 1
                     s0 = s0 + 1
                     BdyVVdom(j0,1) = s-1
                     BdyVVdom(j0,2) = sm-1
                     BdyVVdom(j0,3) = Ivertex
                   endif
               enddo
            enddo
         enddo
!        ----------------------------- 
         DimBdyVVdom(s,1) = s0    !<---- Number of components     
         DimBdyVVdom(s,3) = j0    !<---- Final index (subdomain)   
      enddo
!      ________________________________________________________
!     |                                                        |
!     |              Extra cell-centers (overlaping)           |
!     |________________________________________________________|

      allocate(ExtraCCdom(TotalSubdomains*N_CELL0global,3))  
      allocate(DimExtraCCdom(TotalSubdomains,3)) 
      allocate(IniExtraFrom(TotalSubdomains,TotalSubdomains)) 
      allocate(FinExtraFrom(TotalSubdomains,TotalSubdomains)) 
      allocate(DimExtraFrom(TotalSubdomains,TotalSubdomains)) 

      j0 = 0
      do s=1,TotalSubdomains           
         s0 = 0
         DimExtraCCdom(s,2) = j0 + 1 !<----Initial index (subdomain)
!        --------------------------------------
         IniExtraFrom(s,s) = 1      
         DimExtraFrom(s,s) = 0      
!        --------------------------------------
         do m=1,TotalSubdomains-1         
            k0 = 0
            sm = mod(s+m,TotalSubdomains)
            if (sm.eq.0) sm = TotalSubdomains
!           <<<< Initial position of extra vectors coming from sm
            IniExtraFrom(s,sm) = DimExtraCCdom(s,2) + s0   
!           >>>>
            do nv=DimBdyVVdom(s,2),DimBdyVVdom(s,3)
               Ivertex = BdyVVdom(nv,3)         
               do j=DimCCdom(sm,2),DimCCdom(sm,3)   
                  Jelem = CCdom(j,1)
                  jv1 =  No_vp_global(Jelem,1)
                  jv2 =  No_vp_global(Jelem,2)
                  jv3 =  No_vp_global(Jelem,3)
                  if ((Ivertex.eq.jv1).or. & 
                      (Ivertex.eq.jv2).or. &
                      (Ivertex.eq.jv3)) then
!                     -------------------------
!                     Avoid repetition
                      tag = 0
                      do i=DimExtraCCdom(s,2),j0
                         if (Jelem.eq.ExtraCCdom(i,3)) tag = 1
                      enddo
                      if (tag.eq.0) then
                         j0 = j0 + 1
                         s0 = s0 + 1
                         k0 = k0 + 1
                         ExtraCCdom(j0,1) = s-1
                         ExtraCCdom(j0,2) = sm-1
                         ExtraCCdom(j0,3) = Jelem
                      endif
!                     -------------------------
                  endif
               enddo
            enddo
!           <<<< Number of extra vectors coming from sm
            DimExtraFrom(s,sm) = k0 
            FinExtraFrom(s,sm) = j0 
!           >>>>   
         enddo
!        ---------------------------------------    
         DimExtraCCdom(s,1) = s0  !<---- Number of components   
         DimExtraCCdom(s,3) = j0  !<---- Final index (subdomain)
      enddo
!      ________________________________________________________
!     |                                                        |
!     |   Interior cell-centers (without internal overlaping)  |
!     |________________________________________________________|

      allocate(InterCCdom(N_CELL0global,2))  
      allocate(DimInterCCdom(TotalSubdomains,3)) 
      allocate(BoundCCdom(N_CELL0global,2))  
      allocate(DimBoundCCdom(TotalSubdomains,3))  

      i0 = 0
      j0 = 0
      do s=1,TotalSubdomains  
         s0 = 0
         k0 = 0
         DimInterCCdom(s,2) = j0 + 1     !<----Initial index (subdomain)
         DimBoundCCdom(s,2) = i0 + 1     !<----Initial index (subdomain)
!        --------------------------------------
         do i=DimCCdom(s,2),DimCCdom(s,3)   
            Ielem = CCdom(i,1)
            tag = s
            do m=1,TotalSubdomains-1         
               sm = mod(s+m,TotalSubdomains)
               if (sm.eq.0) sm = TotalSubdomains
               do j=DimExtraCCdom(sm,2),DimExtraCCdom(sm,3)
                  Jelem = ExtraCCdom(j,3)
                  if (Ielem.eq.Jelem) then
                      tag = sm
                  endif
               enddo
            enddo
            if (tag.eq.s) then
               j0 = j0 + 1
               s0 = s0 + 1
               InterCCdom(j0,1) = Ielem
               InterCCdom(j0,2) = s
            else
               i0 = i0 + 1
               k0 = k0 + 1
               BoundCCdom(i0,1) = Ielem
               BoundCCdom(i0,2) = tag
            endif
         enddo
!        ---------------------------------------    
         DimInterCCdom(s,1) = s0         !<---- Number of components   
         DimInterCCdom(s,3) = j0         !<---- Final index (subdomain)

         DimBoundCCdom(s,1) = k0         !<---- Number of components   
         DimBoundCCdom(s,3) = i0         !<---- Final index (subdomain)
      enddo
      
!*********************************************************************!
!                                                                     !
!        Information for topology: neighbors, edges & indexes         !
!                                                                     !
!*********************************************************************!

      allocate(Eedges(TotalSubdomains,TotalSubdomains))  
      allocate(EdgesVec(TotalSubdomains*TotalSubdomains))  
      allocate(DimEdges(TotalSubdomains)) 
      allocate(Iindex(TotalSubdomains)) 

      NNbefore = 0
      j0 = 0
      do s=1,TotalSubdomains 
         Iini = DimExtraCCdom(s,2)
         Ifin = DimExtraCCdom(s,3)
!        ----------------------
         s0 = 1
         Eedges(1,s) = ExtraCCdom(Iini,2)
         j0 = j0 + 1
         EdgesVec(j0)= Eedges(1,s)
!        ----------------------
         do i=Iini+1,Ifin  
            Ielem = ExtraCCdom(i,2)
            tag = 0
            do j=1,s0
               if (Ielem.eq.Eedges(j,s)) tag = 1
            enddo
            if (tag.eq.0) then
                s0 = s0 + 1
                Eedges(s0,s) = Ielem
                j0 = j0 + 1
                EdgesVec(j0)= Eedges(s0,s)
            endif
         enddo
!        ----------------------
         DimEdges(s) = s0
         Iindex(s) = NNbefore + DimEdges(s)
         NNbefore = Iindex(s)
      enddo
      
!     _________________________________________________________________
!     Final calculation

      Nindex = TotalSubdomains 
      Nedges = Iindex(TotalSubdomains) 
!     __________________________
!     Allocate index and edges 
      allocate(indexval(0:Nindex-1),edges(0:Nedges-1))

!     __________________________
!     Index values of the graph (read from a file)
      do s=1,Nindex
         indexval(s-1) = Iindex(s)
      enddo
!     __________________________
!     Edge values of the graph (read from a file)
      do s=1,Nedges
         edges(s-1)  = EdgesVec(s)
      enddo

!*********************************************************************!
!                                                                     !
!                             Save results                            !
!                                                                     !
!*********************************************************************!

!     ____________________________________________________________
!     Save indixes  of cell-centers 

      irec = 100
      do s=1,TotalSubdomains
         filen='../output/MPImatlab/CC_proc   .txt'
         write(filen(28:30),'(i3.3)') s-1
         open(irec,file=filen)
         Iini = DimCCdom(s,2)
         Ifin = DimCCdom(s,3)
         do i=Iini,Ifin
            write(irec,*),CCdom(i,1:2)
         enddo
         close(irec)
      enddo
!     ____________________________________________________________
!     Save indixes of cell-centers boundary for each subdomain 

      irec = 110
      do s=1,TotalSubdomains
         filen2='../output/MPImatlab/CCbdy_proc   .txt'
         write(filen2(31:33),'(i3.3)') s-1
         open(irec,file=filen2)
         Iini = DimBdyCCdom(s,2)
         Ifin = DimBdyCCdom(s,3)
         do i=Iini,Ifin
            write(irec,*),BdyCCdom(i,1:4)
         enddo
         close(irec)
      enddo
!     ____________________________________________________________
!     Save indices of vextex values 

      irec = 120
      do s=1,TotalSubdomains
         filen='../output/MPImatlab/VV_proc   .txt'
         write(filen(28:30),'(i3.3)') s-1
         open(irec,file=filen)
         Iini = DimVVdom(s,2)
         Ifin = DimVVdom(s,3)
         do i=Iini,Ifin
            write(irec,*),VVdom(i,1:2)
         enddo
         close(irec)
      enddo
!     ____________________________________________________________
!     Save indices of vextex values 

      irec = 130
      do s=1,TotalSubdomains
         filen2='../output/MPImatlab/VVbdy_proc   .txt'
         write(filen2(31:33),'(i3.3)') s-1
         open(irec,file=filen2)
         Iini = DimBdyVVdom(s,2)
         Ifin = DimBdyVVdom(s,3)
         do i=Iini,Ifin
            write(irec,*),BdyVVdom(i,1:3)
         enddo
         close(irec)
      enddo
!     ____________________________________________________________
!     Save indices of vextex values 

      irec = 140
      do s=1,TotalSubdomains
         filen2='../output/MPImatlab/CCext_proc   .txt'
         write(filen2(31:33),'(i3.3)') s-1
         open(irec,file=filen2)
         Iini = DimExtraCCdom(s,2)
         Ifin = DimExtraCCdom(s,3)
         do i=Iini,Ifin
            write(irec,*),ExtraCCdom(i,1:3)
         enddo
         close(irec)
      enddo

!     ____________________________________________________________
!     Cell-center and vertex coordinates  

      open(150,file='../output/MPImatlab/xv.txt')
      open(160,file='../output/MPImatlab/yv.txt')
      open(170,file='../output/MPImatlab/xc.txt')
      open(180,file='../output/MPImatlab/yc.txt')
      open(190,file='../output/MPImatlab/jv.txt')
      do nv=1,N_VERTglobal
         write(150,*), xv_global(nv)
         write(160,*), yv_global(nv)
      enddo
      do i=1,N_CELL0global
	 jv1=No_vp_global(i,1)
	 jv2=No_vp_global(i,2)
	 jv3=No_vp_global(i,3)
         xxc = (xv_global(jv1)+xv_global(jv2)+xv_global(jv3))/3.0d0
         yyc = (yv_global(jv1)+yv_global(jv2)+yv_global(jv3))/3.0d0
         write(170,*),xxc
         write(180,*),yyc
         write(190,*),jv1,jv2,jv3
      enddo
      close(150)
      close(160)
      close(170)
      close(180)
      close(190)

!     ____________________________________________________________
!     Save Topology
      open(710,file='../output/Parallel/Topology.dat')
      write(710,*) Nprocs,Nindex,Nedges
      write(710,*) indexval(0:Nindex-1)
      write(710,*) edges(0:Nedges-1)
      close(710)

!*********************************************************************!
!                                                                     !
!                             Finalization                            !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |                       Display values                   |
!     |________________________________________________________|

#     ifdef KeyDisplayParallel
      IF (rang_topo.eq.0) THEN  
      print*,'_________________________________________________________'
      print*,'                                                         ' 
      print*,'       ===========================================       '
      print*,'              MPI: Distribution of elements             '
      print*,'       ===========================================       '
      print*,'                                                         ' 
      do s=1,TotalSubdomains 
         print*,' '
         print*,'SUBDOMAIN:',s-1
         print*,'   Dim cell-centers:        =',DimCCdom(s,1)
         print*,'                   interior =',DimInterCCdom(s,1) 
         print*,'                   boundary =',DimBoundCCdom(s,1) 
         print*,'   --------------------------------------'
         print*,'   Dim extra cell-centers   =',DimExtraCCdom(s,1) 
         do sm=1,TotalSubdomains
            if (DimExtraFrom(s,sm).ne.0)then
               print*,'   From neigh.:',sm-1,'=',&
               DimExtraFrom(s,sm)!,IniExtraFrom(s,sm),FinExtraFrom(s,sm)
            endif
         enddo
         print*,'   --------------------------------------'
         print*,'   Dim vertex               =',DimVVdom(s,1)
         !print*,'   Dim boundary cell-center =',DimBdyCCdom(s,1:3) 
         !print*,'   Dim boundary vertex      =',DimBdyVVdom(s,1:3) 
         print*,' '                                                     
      enddo
      !print *,'_________________________________________________________'
      !print *,'    '
      do s=1,TotalSubdomains 
         print*,'Subdomain:',s-1
         print*,'   No.Neighbors =',DimEdges(s)
         print*,'   Index        =',Iindex(s)
         print*,'   Neighbors    =',Eedges(1:DimEdges(s),s)
         print*,' '
      enddo
      ENDIF
#     endif

!      ________________________________________________________
!     |                                                        |
!     |                       Deallocate                       |
!     |________________________________________________________|

      !deallocate(CCdom)
      !deallocate(DimCCdom)
      !deallocate(VVdom) 
      !deallocate(DimVVdom)
      !deallocate(ExtraCCdom)
      !deallocate(DimExtraCCdom)

      deallocate(BdyCCdom)
      deallocate(DimBdyCCdom)
      deallocate(BdyVVdom) 
      deallocate(DimBdyVVdom)
      deallocate(InterCCdom)
      deallocate(DimInterCCdom)
      deallocate(BoundCCdom)
      deallocate(DimBoundCCdom)
!     --------------
      deallocate(Eedges)  
      deallocate(EdgesVec)  
      deallocate(DimEdges) 
      deallocate(Iindex) 

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
           print*, '      ----> End   subroutine: parallel_distribution'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END


!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                      	   END OF INPUT                               !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

