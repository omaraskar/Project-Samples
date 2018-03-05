!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                Choose a pre-defined model application               !
!                             Jan 2013                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

!---------------------------------------------------------------------!
!                                                                     !
!    This part of the program is called using the command:            !
!                    #include "cppdefs.h"                             !
!                                                                     !
!    Choose the C-preprocessing options by using the command          !
!                    #define   to activate option or                  !
!                    #undef    to deactivate option.                  !
!                                                                     !
!---------------------------------------------------------------------!
!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!                                                                     !
!      Detailed description of all available CPP options.             !
!                                                                     !
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!
!      ________________________________________________________
!     |                                                        |
!     |   All the test cases                                   |
!     |________________________________________________________|

#      undef KeyErosion         /* Application to the erosion Test case */
#      undef KeyClap            /* Application to clapage test-case*/
!      ________________________________________________________
!     |                                                        |
!     |   Compiler                                             |
!     |________________________________________________________|   

#      undef KeyXLF             /* XLF compiler */
#      undef KeyG77             /* GNU compiler */
!      ________________________________________________________
!     |                                                        |
!     |   Parallelization option                               |
!     |________________________________________________________| 

#      undef KeyParallel          /* Parallelization */
#      undef KeyDisplayParallel   /* Display parallel commends */
!      ________________________________________________________
!     |                                                        |
!     |   Mesh and Domain                                      |
!     |________________________________________________________| 

#      undef KeyMeshType1       /* Mesh type 1: square divided in two triangles */
#      undef KeyMeshType2       /* Mesh type 2: square divided in four triangles */
#      undef KeyMeshType3       /* Mesh type 3: BlueKenue unstructured grid */
#      undef KeyMesh1Example3   /* General Mesh designed for the Example 3: NN=16 */
#      undef KeyMesh2Example3   /* General Mesh designed for the Example 3: NN=32 */
#      undef KeyMesh3Example3   /* General Mesh designed for the Example 3: NN=64 */
#      undef KeyMesh4Example3   /* General Mesh designed for the Example 3: NN=128*/
#      undef KeyDomainGeneral   /* Mesh in a general domain */
#      undef KeyMeshFreeSurface /* Mesh used for the free-surface example */
!      ________________________________________________________
!     |                                                        |
!     |   Advection-Diffusion Scheme                           |
!     |________________________________________________________|

#      undef KeyAdvImplicit       /* Solving advection Implicit */
#      undef KeyAdvSemiImplicit   /* Solving advection Semi-implicit */
#      undef KeyAdvExplicit       /* Solving advection Explicit */
#      undef KeyAdvFullExplicit   /* Solving advection Full Explicit */
#      undef KeyAdvCenter         /* 2nd order Central Convection */

#      undef KeyDiffExplicit      /* Solving diffusion explicit */
#      undef KeyDiffImplicit      /* Solving diffusion implicit */
#      undef KeyDiffFullExplicit  /* Solving advection full explicit */
!      ________________________________________________________
!     |                                                        |
!     |   Linear system solvers                                |
!     |________________________________________________________|

#      undef KeySOR             /* Succesive-Over-Relaxation */
#      undef KeySiOR            /* Simultaneous-Over-Relaxation */
#      undef KeyGMRES           /* GMRES solver */
#      undef KeyPseudoTime      /* Pseudo time solver */
!      ________________________________________________________
!     |                                                        |
!     |   Computing options and parameters                     |
!     |________________________________________________________|

#      undef KeyMask            /* Mask velocity on a predifined area of the mesh */ 
#      undef KeyFluxLimiter     /* Using the flux limiter in the edge approx */
#      undef KeyRigid           /* No free surface calculation */
!      ________________________________________________________
!     |                                                        |
!     |   Interpolation 2D                                     |
!     |________________________________________________________|

#      undef KeyInterpoNew      /* Vertex interpolation: distance weighting, new formulation Jul2014*/
#      undef KeyInterpoArea     /* Vertex interpolation: area weighting*/
#      undef KeyInterpoDist     /* Vertex interpolation: distance C-V weighting*/
#      undef KeyInterpoDist2    /* Vertex interpolation: distance C-V weighting more points*/
#      undef KeyInterpoLSM      /* Vertex interpolation: using the LSM*/
#      undef KeyInterpoMix      /* Vertex interpolation: using a mix of LSM & dist*/
#      undef KeyUseInterGhost   /* Using ghost points during vertex interpolation*/
#      undef KeyUseVertexNeumann/* Interpolatiopn for Neumann BC*/
!      ________________________________________________________
!     |                                                        |
!     |   Display options                                      |
!     |________________________________________________________|

#      undef KeyDbg             /* Print information during execution */
#      undef KeyDisplay         /* Print lines-information during execution */

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!                                                                     !
!                           PROBLEM CHOICE                            !
!                                                                     !
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

!      ________________________________________________________
!     |                                                        |
!     |                  Choose a compiler                     |
!     |________________________________________________________|

#     define KeyG77
!      ________________________________________________________
!     |                                                        |
!     |   ***********    Choose a test Case   **************   |
!     |________________________________________________________|

!#     define KeyErosion

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!                                                                     !
!                   Options for each problem case.                    !
!                                                                     !
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

!      ________________________________________________________
!     |                                                        |
!     |                   EROSION TEST CASE                    |
!     |________________________________________________________|

!#if defined KeyErosion
!     --------------------------------------------------
!     Run in parallel:
!#            define KeyParallel
!     --------------------------------------------------
!     Linear Solver:
!#             define KeySOR
#             define KeySiOR
!#             define KeyGMRES
!#             define KeyPseudoTime
!     --------------------------------------------------
!     Type of interpolation:
#             define KeyInterpoNew
!#             define KeyInterpoDist
!#             define KeyInterpoLSM 
!#             define KeyInterpoMix 
!     Use ghost points during interpolation:
!#             define KeyUseInterGhost
!     --------------------------------------------------
!     Free-surface (advection) option:
#             define KeyFreeFullExplicit
!#             define KeyFreeExplicit
!#             define KeyFreeSemiImplicit
!#             define KeyFreeImplicit  
!     --------------------------------------------------
!     Advection:
!#            define KeyFluxLimiter
#             define KeyAdvCenter
#             define KeyAdvFullExplicit
!#             define KeyAdvExplicit
!#             define KeyAdvSemiImplicit
!#             define KeyAdvImplicit
!     --------------------------------------------------
!     Diffusion:
!#            define KeyDiffFullExplicit
!     --------------------------------------------------
!     Domain (only use to save data): 
!#            define KeyDomainGeneral  
!#            define KeyMeshType2
#             define KeyMeshType3
!#            define KeyMesh3Example3 
!#            define KeyMeshFreeSurface
!     --------------------------------------------------
!     Print information during execution
!#            define KeyDbg
!#            define KeyDisplay
!#            define KeyDisplayParallel  
!#endif

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!                                                                     !
!                                  END                                !
!                                                                     !
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

