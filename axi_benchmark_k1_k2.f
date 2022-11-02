C     **************************                                              
C     * PDE2D 9.7 MAIN PROGRAM *                                              
C     **************************                                              
C     *** 2D PROBLEM SOLVED (GALERKIN METHOD) ***                             
C##############################################################################
C     Is double precision mode to be used?  Double precision is recommended.  #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + If double precision mode is used, variables and functions assigned   +#
C     + names beginning with a letter in the range A-H or O-Z will be DOUBLE +#
C     + PRECISION, and you should use double precision constants and FORTRAN +#
C     + expressions throughout; otherwise such variables and functions will  +#
C     + be of type REAL.  In either case, variables and functions assigned   +#
C     + names beginning with I,J,K,L,M or N will be of INTEGER type.         +#
C     +                                                                      +#
C     + It is possible to convert a single precision PDE2D program to double +#
C     + precision after it has been created, using an editor.  Just change   +#
C     + all occurrences of "real" to "double precision"                      +#
C     +                    " tdp" to "dtdp"  (note leading blank)            +#
C     + Any user-written code or routines must be converted "by hand", of    +#
C     + course.  To convert from double to single, reverse the changes.      +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C##############################################################################
      implicit double precision (a-h,o-z)
      parameter (neqnmx=  99)
      parameter (ndelmx=  20)
      parameter (nxgrid=0,nygrid=0,npgrid=0,nqgrid=0)
C##############################################################################
C     NV0 = number of vertices in initial triangulation                       #
C##############################################################################
      PARAMETER (NV0 = 8)            
C##############################################################################
C     NT0 = number of triangles in initial triangulation                      #
C##############################################################################
      PARAMETER (NT0 = 8)            
C##############################################################################
C     How many differential equations (NEQN) are there in your problem?       #
C##############################################################################
      PARAMETER (NEQN = 2)            
C        DIMENSIONS OF WORK ARRAYS                                            
C        SET TO 1 FOR AUTOMATIC ALLOCATION                                    
      PARAMETER (IRWK8Z=           1)
      PARAMETER (IIWK8Z=           1)
      PARAMETER (NXP8Z=101,NYP8Z=101,KDEG8Z=1,NBPT8Z=51)
C##############################################################################
C     The solution is normally saved on a NX+1 by NY+1 rectangular grid of    #
C     points                                                                  #
C                   (XA + I*(XB-XA)/NX , YA + J*(YB-YA)/NY)                   #
C     I=0,...,NX, J=0,...,NY.  Enter values for NX and NY.  Suggested values  #
C     are NX=NY=25.                                                           #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + If you want to save the solution at an arbitrary user-specified      +#
C     + set of points, set NY=0 and NX+1=number of points.  In this case you +#
C     + can request tabular output of the solution, but you cannot make any  +#
C     + solution plots.                                                      +#
C     +                                                                      +#
C     + If you set NEAR8Z=1 in the main program, the values saved at each    +#
C     + output point will actually be the solution as evaluated at a nearby  +#
C     + integration point.  For most problems this obviously will produce    +#
C     + less accurate output or plots, but for certain (rare) problems, a    +#
C     + solution component may be much less noisy when plotted only at       +#
C     + integration points.                                                  +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C##############################################################################
      PARAMETER (NX = 40)             
      PARAMETER (NY = 40)             
C##############################################################################
C     The solution will be saved (for possible postprocessing) at the NSAVE+1 #
C     time points                                                             #
C                          T0 + K*(TF-T0)/NSAVE                               #
C     K=0,...,NSAVE.  Enter a value for NSAVE.                                #
C                                                                             #
C     If a user-specified constant time step is used, NSTEPS must be an       #
C     integer multiple of NSAVE.                                              #
C##############################################################################
      PARAMETER (NSAVE = 5)          
      common/parm8z/ pi,h1    ,h2    ,h3    ,h4    ,w1    ,w2   ,F0    
      dimension vxy(2,nv0+1),iabc(3,nt0+1),iarc(nt0+1),xgrid(nxgrid+1),y       
     &grid(nygrid+1),ixarc(2),iyarc(2),pgrid(npgrid+1),qgrid(nqgrid+1),i       
     &parc(2),iqarc(2),xbd8z(nbpt8z,nt0+4),ybd8z(nbpt8z,nt0+4),xout8z(0:       
     &nx,0:ny),yout8z(0:nx,0:ny),inrg8z(0:nx,0:ny),xcross(100),ycross(10       
     &0),tout8z(0:nsave),xd0(ndelmx),yd0(ndelmx)                               
C      dimension xres8z(nxp8z),yres8z(nyp8z),ures8z(neqn,nxp8z,nyp8z)          
      allocatable iwrk8z(:),rwrk8z(:)                                          
C      dimension iwrk8z(iiwk8z),rwrk8z(irwk8z)                                 
      character*40 title                                                       
      logical plot,symm,fdiff,evcmpx,crankn,noupdt,adapt,nodist,fillin,e       
     &con8z,ncon8z,restrt,gridid                                               
      common/dtdp14/ sint8z(20),bint8z(20),slim8z(20),blim8z(20)               
      common/dtdp15/ evlr8z,ev0r,evli8z,ev0i,evcmpx                            
      common/dtdp16/ p8z,evr8z(50),evi8z(50)                                   
      common/dtdp19/ toler(neqnmx),adapt                                       
      common/dtdp22/ nxa8z,nya8z,ifgr8z,kd8z,nbp8z                             
      common/dtdp23/ work8z(nxp8z*nyp8z+6)                                     
      common/dtdp30/ econ8z,ncon8z                                             
      common/dtdp46/ eps8z,cgtl8z,npmx8z,itype,near8z                          
      common/dtdp63/ amin8z(3*neqnmx),amax8z(3*neqnmx)                         
      common/dtdp65/ intri,iotri                                               
      common/dtdp76/ mdim8z,nx18z,ny18z,xa,xb,ya,yb,uout(0:nx,0:ny,3,neq       
     &n,0:nsave)                                                               
      pi = 4.0*atan(1.d0)                                                      
      nxa8z = nxp8z                                                            
      nya8z = nyp8z                                                            
      nx18z = nx+1                                                             
      ny18z = ny+1                                                             
      mdim8z = 3                                                               
      kd8z = kdeg8z                                                            
      nbp8z = nbpt8z                                                           
C##############################################################################
C     If you don't want to read the FINE PRINT, default NPROB.                #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + If you want to solve several similar problems in the same run, set   +#
C     + NPROB equal to the number of problems you want to solve.  Then NPROB +#
C     + loops through the main program will be done, with IPROB=1,...,NPROB, +#
C     + and you can make the problem parameters vary with IPROB.  NPROB      +#
C     + defaults to 1.                                                       +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C##############################################################################
      NPROB = 1                                                                
      do 78755 iprob=1,nprob                                                   
C##############################################################################
C     PDE2D solves the time-dependent system (note: U,A,B,F,FB,GB,U0 may be   #
C     vectors, C,RHO may be matrices):                                        #
C                                                                             #
C       C(X,Y,T,U,Ux,Uy)*d(U)/dT = d/dX* A(X,Y,T,U,Ux,Uy)                     #
C                                + d/dY* B(X,Y,T,U,Ux,Uy)                     #
C                                -       F(X,Y,T,U,Ux,Uy)                     #
C                                                                             #
C     or the steady-state system:                                             #
C                                                                             #
C                            d/dX* A(X,Y,U,Ux,Uy)                             #
C                          + d/dY* B(X,Y,U,Ux,Uy)                             #
C                                = F(X,Y,U,Ux,Uy)                             #
C                                                                             #
C     or the linear and homogeneous eigenvalue system:                        #
C                                                                             #
C                            d/dX* A(X,Y,U,Ux,Uy)                             #
C                          + d/dY* B(X,Y,U,Ux,Uy)                             #
C                                = F(X,Y,U,Ux,Uy) + lambda*RHO(X,Y)*U         #
C                                                                             #
C     in an arbitrary two-dimensional region, R, with 'fixed' boundary        #
C     conditions on part of the boundary:                                     #
C                                                                             #
C           U = FB(X,Y,[T])                                                   #
C                                                                             #
C     and 'free' boundary conditions on the other part:                       #
C                                                                             #
C         A*nx + B*ny = GB(X,Y,[T],U,Ux,Uy)                                   #
C                                                                             #
C     For time-dependent problems there are also initial conditions:          #
C                                                                             #
C            U = U0(X,Y)   at T=T0                                            #
C                                                                             #
C     Here Ux,Uy represent the (vector) functions dU/dX,dU/dY, and (nx,ny)    #
C     represents the unit outward normal to the boundary.                     #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + If your PDEs involve the solution at points other than (X,Y), the    +#
C     + function                                                             +#
C     +             (D)OLDSOL2(IDER,IEQ,XX,YY,KDEG)                          +#
C     + will interpolate (using interpolation of degree KDEG=1,2 or 3) to    +#
C     + (XX,YY) the function saved in UOUT(*,*,IDER,IEQ,ISET) on the last    +#
C     + time step or iteration (ISET) for which it has been saved.  Thus,    +#
C     + for example, if IDER=1, this will return the latest value of         +#
C     + component IEQ of the solution at (XX,YY), assuming this has not been +#
C     + modified using UPRINT...  If your equations involve integrals of the +#
C     + solution, for example, you can use (D)OLDSOL2 to approximate these   +#
C     + using the solution from the last time step or iteration.             +#
C     +                                                                      +#
C     + CAUTION: For a steady-state or eigenvalue problem, you must reset    +#
C     + NOUT=1 if you want to save the solution each iteration.              +#
C     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
C     + A system of NEQN complex partial differential equations must be      +#
C     + written as a system of 2*NEQN real equations, by separating the      +#
C     + equations into their real and imaginary parts.  However, note that   +#
C     + the complex arithmetic abilities of FORTRAN can be used to simplify  +#
C     + this separation.  For example, the complex PDE:                      +#
C     +       I*(Uxx+Uyy) = 1/(1+U**10),  where U = UR + UI*I                +#
C     + would be difficult to split up analytically, but using FORTRAN       +#
C     + expressions it is easy:                                              +#
C     +   A1 = -UIx,  B1 = -UIy,  F1 =  REAL(1.0/(1.0+CMPLX(UR,UI)**10))     +#
C     +   A2 =  URx,  B2 =  URy,  F2 = AIMAG(1.0/(1.0+CMPLX(UR,UI)**10))     +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C     You may now define global parameters, which may be referenced in any    #
C     of the "FORTRAN expressions" you input throughout the rest of this      #
C     interactive session.  You will be prompted alternately for parameter    #
C     names and their values; enter a blank name when you are finished.       #
C                                                                             #
C     Parameter names are valid FORTRAN variable names, starting in           #
C     column 1.  Thus each name consists of 1 to 6 alphanumeric characters,   #
C     the first of which must be a letter.  If the first letter is in the     #
C     range I-N, the parameter must be an integer.                            #
C                                                                             #
C     Parameter values are either FORTRAN constants or FORTRAN expressions    #
C     involving only constants and global parameters defined on earlier       #
C     lines.  They may also be functions of the problem number IPROB, if      #
C     you are solving several similar problems in one run (NPROB > 1).  Note  #
C     that you are defining global CONSTANTS, not functions; e.g., parameter  #
C     values may not reference any of the independent or dependent variables  #
C     of your problem.                                                        #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + If you define other parameters here later, using an editor, you must +#
C     + add them to COMMON block /PARM8Z/ everywhere this block appears, if  +#
C     + they are to be "global" parameters.                                  +#
C     +                                                                      +#
C     + The variable PI is already included as a global parameter, with an   +#
C     + accurate value 3.14159...                                            +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C##############################################################################
C        *******INITIAL TRIANGULATION OPTION                                   
      INTRI =          3                                                       
C        *******SET IOTRI = 1 TO DUMP FINAL TRIANGULATION TO FILE pde2d.tri    
      IOTRI = 0                                                                
C##############################################################################
C     For a general region, an initial triangulation is constructed           #
C     which generally consists of only as many triangles as needed to define  #
C     the region and to satisfy the following rules (triangles adjacent       #
C     to a curved boundary may be considered to have one curved edge):        #
C                                                                             #
C     1. The end points of each arc are included as vertices in the           #
C        triangulation.                                                       #
C                                                                             #
C     2. No vertex of any triangle may touch another in a point which is      #
C        not a vertex of the other triangle.                                  #
C                                                                             #
C     3. No triangle may have all three vertices on the boundary.             #
C                                                                             #
C     Now enter the number of vertices in the initial triangulation (NV0)     #
C     and the vertices (VXY(1,i),VXY(2,i)), i=1,...,NV0, when prompted.       #
C                                                                             #
C     The vertices may be numbered in any order, but that order will define   #
C     the vertex numbers referred to in the next list.                        #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + Normally, none of the elements of VXY, IABC, IARC should be          +#
C     + defaulted.  If any are defaulted, these initial triangulation arrays +#
C     + will be read from the file 'pde2d.tri'; in this case the values of   +#
C     + NV0 and NT0 entered below must be the same as those in this file.    +#
C     +                                                                      +#
C     + A file 'pde2d.tri' is normally created when the final triangulation  +#
C     + from another program (cases INTRI=1,2 or 3) is dumped.  Reset IOTRI  +#
C     + to 1 in the other program, to cause the final triangulation to be    +#
C     + dumped into pde2d.tri.                                               +#
C     +                                                                      +#
C     + An example situation where this is useful: your region is a thin     +#
C     + annulus, which is very difficult to triangulate by hand (INTRI=3)    +#
C     + but you can't use the INTRI=2 option because you need one of the     +#
C     + r=constant curves to be an interface, or because one of the boundary +#
C     + arcs is part fixed and part free.  In this case you can generate     +#
C     + the initial triangulation with an INTRI=2 program and dump it to     +#
C     + 'pde2d.tri,' then modify the arc numbers in the file appropriately   +#
C     + before reading with another, INTRI=3, program.                       +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C##############################################################################
      call dtdpwx(vxy,2*nv0,0)                                                 
C                                                                              
      VXY(1,1) =                                                               
     & 0                                                                       
      VXY(2,1) =                                                               
     & 0                                                                       
C                                                                              
      VXY(1,2) =                                                               
     & 4500                                                                    
      VXY(2,2) =                                                               
     & 0                                                                       
C                                                                              
      VXY(1,3) =                                                               
     & 4500                                                                    
      VXY(2,3) =                                                               
     & 1900                                                                    
C                                                                              
      VXY(1,4) =                                                               
     & 4500                                                                    
      VXY(2,4) =                                                               
     & 3000                                                                    
C                                                                              
      VXY(1,5) =                                                               
     & 0                                                                       
      VXY(2,5) =                                                               
     & 3000                                                                    
C                                                                              
      VXY(1,6) =                                                               
     & 0                                                                       
      VXY(2,6) =                                                               
     & 1900                                                                    
C                                                                              
      VXY(1,7) =                                                               
     & 2000                                                                    
      VXY(2,7) =                                                               
     & 2500                                                                    
C                                                                              
      VXY(1,8) =                                                               
     & 2000                                                                    
      VXY(2,8) =                                                               
     & 1000                                                                    
C##############################################################################
C     Now enter the number of triangles in the initial triangulation (NT0),   #
C     and for each triangle k, give the vertex numbers of vertices a,b,c:     #
C                 IABC(1,k) , IABC(2,k) , IABC(3,k)                           #
C     where the third vertex (c) is not on the boundary of R (or on a curved  #
C     interface arc), and the number,                                         #
C                             IARC(k)                                         #
C     of the arc cut off by the base, ab, of the triangle.  Put IARC(k)=0     #
C     if the base of triangle k does not intersect any boundary (or curved    #
C     interface) arc.  Recall that negative arc numbers correspond to         #
C     'fixed' boundary conditions, and positive arc numbers ( < 1000)         #
C     correspond to 'free' boundary conditions.                               #
C                                                                             #
C     The order of this list defines the initial triangle numbers.            #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + Any of the partial differential equation coefficients or other       +#
C     + functions of X and Y which you are asked to supply later may be      +#
C     + defined as functions of a variable 'KTRI', which holds the number    +#
C     + of the INITIAL triangle in which (X,Y) lies.  This is useful when    +#
C     + the region is a composite of different materials, so that some       +#
C     + material parameters have different values in different subregions.   +#
C     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
C     + If any interfaces between subregions are curved, the interface arcs  +#
C     + may be assigned unique arc numbers of 1000 and above, and treated    +#
C     + like boundary arcs in the initial triangulation definition.  (You    +#
C     + will not be allowed to specify boundary conditions on them, however.)+#
C     + The triangulation refinement will follow the interface arcs, so that +#
C     + no final triangles will straddle an interface.                       +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C##############################################################################
C                                                                              
      IABC(1,1) =                                                              
     & 1                                                                       
      IABC(2,1) =                                                              
     & 2                                                                       
      IABC(3,1) =                                                              
     & 8                                                                       
      IARC(1) =                                                                
     & 1                                                                       
C                                                                              
      IABC(1,2) =                                                              
     & 2                                                                       
      IABC(2,2) =                                                              
     & 3                                                                       
      IABC(3,2) =                                                              
     & 8                                                                       
      IARC(2) =                                                                
     & 2                                                                       
C                                                                              
      IABC(1,3) =                                                              
     & 3                                                                       
      IABC(2,3) =                                                              
     & 6                                                                       
      IABC(3,3) =                                                              
     & 8                                                                       
      IARC(3) =                                                                
     & 0                                                                       
C                                                                              
      IABC(1,4) =                                                              
     & 6                                                                       
      IABC(2,4) =                                                              
     & 1                                                                       
      IABC(3,4) =                                                              
     & 8                                                                       
      IARC(4) =                                                                
     & 4                                                                       
C                                                                              
      IABC(1,5) =                                                              
     & 6                                                                       
      IABC(2,5) =                                                              
     & 3                                                                       
      IABC(3,5) =                                                              
     & 7                                                                       
      IARC(5) =                                                                
     & 0                                                                       
C                                                                              
      IABC(1,6) =                                                              
     & 3                                                                       
      IABC(2,6) =                                                              
     & 4                                                                       
      IABC(3,6) =                                                              
     & 7                                                                       
      IARC(6) =                                                                
     & 2                                                                       
C                                                                              
      IABC(1,7) =                                                              
     & 4                                                                       
      IABC(2,7) =                                                              
     & 5                                                                       
      IABC(3,7) =                                                              
     & 7                                                                       
      IARC(7) =                                                                
     & -3                                                                      
C                                                                              
      IABC(1,8) =                                                              
     & 5                                                                       
      IABC(2,8) =                                                              
     & 6                                                                       
      IABC(3,8) =                                                              
     & 7                                                                       
      IARC(8) =                                                                
     & 4                                                                       
      call dtdpu(vxy,nv0,iabc,iarc,nt0)                                        
C##############################################################################
C     How many triangles (NTF) are desired for the final triangulation?       #
C##############################################################################
      NTF =        1000                                                         
C##############################################################################
C     If you don't want to read the FINE PRINT, enter ISOLVE = 4.             #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + The following linear system solvers are available:                   +#
C     +                                                                      +#
C     + 1. Band method                                                       +#
C     +               The band solver uses a reverse Cuthill-McKee ordering. +#
C     + 2. Frontal method                                                    +#
C     +               This is an out-of-core version of the band solver.     +#
C     + 3. Jacobi bi-conjugate gradient method                               +#
C     +               This is a preconditioned bi-conjugate gradient, or     +#
C     +               Lanczos, iterative method.  (This solver is MPI-       +#
C     +               enhanced, if MPI is available.)  If you want to        +#
C     +               override the default convergence tolerance, set a      +#
C     +               new relative tolerance CGTL8Z in the main program.     +#
C     + 4. Sparse direct method                                              +#
C     +               This is based on Harwell Library routines MA27/MA37,   +#
C     +               developed by AEA Industrial Technology at Harwell      +#
C     +               Laboratory, Oxfordshire, OX11 0RA, United Kingdom      +#
C     +               (used by permission).                                  +#
C     + 5. Local solver                                                      +#
C     +               Choose this option ONLY if alternative linear system   +#
C     +               solvers have been installed locally.  See subroutines  +#
C     +               (D)TD3M, (D)TD3N in file (d)subs.f for instructions    +#
C     +               on how to add local solvers.                           +#
C     + 6. MPI-based parallel band solver                                    +#
C     +               This is a parallel solver which runs efficiently on    +#
C     +               multiple processor machines, under MPI.  It is a       +#
C     +               band solver, with the matrix distributed over the      +#
C     +               available processors.  Choose this option ONLY if the  +#
C     +               solver has been activated locally.  See subroutine     +#
C     +               (D)TD3O in file (d)subs.f for instructions on how to   +#
C     +               activate this solver and the MPI-enhancements to the   +#
C     +               conjugate gradient solver.                             +#
C     +                                                                      +#
C     + Enter ISOLVE = 1,2,3,4,5 or 6 to select a linear system solver.      +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C                                                                             #
C     If you don't want to read the FINE PRINT, enter ISOLVE = 4.             #
C##############################################################################
      ISOLVE =          1                                                      
C##############################################################################
C     Enter the element degree (1,2,3 or 4) desired.  A suggested value is    #
C     IDEG = 3.                                                               #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + A negative value for IDEG can be entered, and elements of degree     +#
C     + ABS(IDEG) will be used, with a lower order numerical integration     +#
C     + scheme.  This results in a slight increase in speed, but negative    +#
C     + values of IDEG are normally not recommended.                         +#
C     +                                                                      +#
C     + The spatial discretization error is O(h**2), O(h**3), O(h**4) or     +#
C     + O(h**5) when IDEG = 1,2,3 or 4, respectively, is used, where h is    +#
C     + the maximum triangle diameter, even if the region is curved.         +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C##############################################################################
      IDEG =          2                                                        
C        *******TIME-DEPENDENT PROBLEM                                         
      itype = 2                                                                
C##############################################################################
C     Enter the initial time value (T0) and the final time value (TF), for    #
C     this time-dependent problem.  T0 defaults to 0.                         #
C                                                                             #
C     TF is not required to be greater than T0.                               #
C##############################################################################
      T0 = 0.0                                                                 
      T0 =                                                                     
     & 0                                                                       
      TF =                                                                     
     & 50000.*365.*24.*3600.                                                   
C##############################################################################
C     Do you want the time step to be chosen adaptively?  If you answer       #
C     'yes', you will then be prompted to enter a value for TOLER(1), the     #
C     local relative time discretization error tolerance.  The default is     #
C     TOLER(1)=0.01.  If you answer 'no', a user-specified constant time step #
C     will be used.  We suggest that you answer 'yes' and default TOLER(1)    #
C     (although for certain linear problems, a constant time step may be much #
C     more efficient).                                                        #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + If a negative value is specified for TOLER(1), then ABS(TOLER(1)) is +#
C     + taken to be the "absolute" error tolerance.  If a system of PDEs is  +#
C     + solved, by default the error tolerance specified in TOLER(1) applies +#
C     + to all variables, but the error tolerance for the J-th variable can  +#
C     + be set individually by specifying a value for TOLER(J) using an      +#
C     + editor, after the end of the interactive session.                    +#
C     +                                                                      +#
C     + Each time step, two steps of size dt/2 are taken, and that solution  +#
C     + is compared with the result when one step of size dt is taken.  If   +#
C     + the maximum difference between the two answers is less than the      +#
C     + tolerance (for each variable), the time step dt is accepted (and the +#
C     + next step dt is doubled, if the agreement is "too" good); otherwise  +#
C     + dt is halved and the process is repeated.  Note that forcing the     +#
C     + local (one-step) error to be less than the tolerance does not        +#
C     + guarantee that the global (cumulative) error is less than that value.+#
C     + However, as the tolerance is decreased, the global error should      +#
C     + decrease correspondingly.                                            +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C##############################################################################
      ADAPT = .TRUE.                                                           
      TOLER(1) = 0.01                                                          
      TOLER(1) =                                                               
     & 0.001                                                                   
      NOUPDT = .FALSE.                                                         
C##############################################################################
C     The time stepsize will be chosen adaptively, between an upper limit     #
C     of DTMAX = (TF-T0)/NSTEPS and a lower limit of 0.0001*DTMAX.  Enter     #
C     a value for NSTEPS (the minimum number of steps).                       #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + If you later turn off adaptive time step control, the time stepsize  +#
C     + will be constant, DT = (TF-T0)/NSTEPS.                               +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C##############################################################################
      NSTEPS =                                                                 
     & 50                                                                     
      dt = (tf-t0)/max(nsteps,1)                                               
C##############################################################################
C     If you don't want to read the FINE PRINT, enter 'no'.                   #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + Is the Crank-Nicolson scheme to be used to discretize time?  If you  +#
C     + answer 'no', a backward Euler scheme will be used.                   +#
C     +                                                                      +#
C     + If a user-specified constant time step is chosen, the second order   +#
C     + Crank Nicolson method is recommended only for problems with very     +#
C     + well-behaved solutions, and the first order backward Euler scheme    +#
C     + should be used for more difficult problems.  In particular, do not   +#
C     + use the Crank Nicolson method if the left hand side of any PDE is    +#
C     + zero, for example, if a mixed elliptic/parabolic problem is solved.  +#
C     +                                                                      +#
C     + If adaptive time step control is chosen, however, an extrapolation   +#
C     + is done between the 1-step and 2-step answers which makes the Euler  +#
C     + method second order, and the Crank-Nicolson method strongly stable.  +#
C     + Thus in this case, both methods have second order accuracy, and both +#
C     + are strongly stable.                                                 +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C##############################################################################
      CRANKN = .FALSE.                                                         
C##############################################################################
C     PDE2D solves the system of equations:                                   #
C                                                                             #
C        C11(X,Y,T,P,Px,Py,TMP,TMPx,TMPy)*d(P)/dT                             #
C      + C12(X,Y,T,P,Px,Py,TMP,TMPx,TMPy)*d(TMP)/dT =                         #
C                               d/dX* A1(X,Y,T,P,Px,Py,TMP,TMPx,TMPy)         #
C                             + d/dY* B1(X,Y,T,P,Px,Py,TMP,TMPx,TMPy)         #
C                             -       F1(X,Y,T,P,Px,Py,TMP,TMPx,TMPy)         #
C        C21(X,Y,T,P,Px,Py,TMP,TMPx,TMPy)*d(P)/dT                             #
C      + C22(X,Y,T,P,Px,Py,TMP,TMPx,TMPy)*d(TMP)/dT =                         #
C                               d/dX* A2(X,Y,T,P,Px,Py,TMP,TMPx,TMPy)         #
C                             + d/dY* B2(X,Y,T,P,Px,Py,TMP,TMPx,TMPy)         #
C                             -       F2(X,Y,T,P,Px,Py,TMP,TMPx,TMPy)         #
C                                                                             #
C     with 'fixed' boundary conditions:                                       #
C                                                                             #
C           P = FB1(X,Y,T)                                                    #
C           TMP = FB2(X,Y,T)                                                  #
C                                                                             #
C     or 'free' boundary conditions:                                          #
C                                                                             #
C          A1*nx + B1*ny = GB1(X,Y,T,P,Px,Py,TMP,TMPx,TMPy)                   #
C          A2*nx + B2*ny = GB2(X,Y,T,P,Px,Py,TMP,TMPx,TMPy)                   #
C                                                                             #
C     and initial conditions:                                                 #
C                                                                             #
C           P = P0(X,Y)     at T=T0                                           #
C           TMP = TMP0(X,Y)                                                   #
C                                                                             #
C     where P(X,Y,T) and TMP(X,Y,T) are the unknowns and C11,C12,F1,A1,B1,    #
C     C21,C22,F2,A2,B2,FB1,FB2,GB1,GB2,P0,TMP0 are user-supplied functions.   #
C                                                                             #
C     note:                                                                   #
C           (nx,ny) = unit outward normal to the boundary                     #
C           Px = d(P)/dX     TMPx = d(TMP)/dX                                 #
C           Py = d(P)/dY     TMPy = d(TMP)/dY                                 #
C                                                                             #
C     Is this problem symmetric?  If you don't want to read the FINE PRINT,   #
C     it is safe to enter 'no'.                                               #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + This problem is called symmetric if each of the matrices             +#
C     +                                                                      +#
C     +    F1.P    F1.Px    F1.Py    F1.TMP    F1.TMPx    F1.TMPy            +#
C     +    A1.P    A1.Px    A1.Py    A1.TMP    A1.TMPx    A1.TMPy            +#
C     +    B1.P    B1.Px    B1.Py    B1.TMP    B1.TMPx    B1.TMPy            +#
C     +    F2.P    F2.Px    F2.Py    F2.TMP    F2.TMPx    F2.TMPy            +#
C     +    A2.P    A2.Px    A2.Py    A2.TMP    A2.TMPx    A2.TMPy            +#
C     +    B2.P    B2.Px    B2.Py    B2.TMP    B2.TMPx    B2.TMPy            +#
C     +                                                                      +#
C     +         C11        C12                                               +#
C     +         C21        C22                                               +#
C     + and                                                                  +#
C     +         GB1.P     GB1.TMP                                            +#
C     +         GB2.P     GB2.TMP                                            +#
C     +                                                                      +#
C     + is always symmetric, where F1.P means d(F1)/d(P), and similarly      +#
C     + for the other terms.  In addition, GB1,GB2 must not depend on        +#
C     + Px,Py,TMPx,TMPy.                                                     +#
C     +                                                                      +#
C     + The memory and execution time are halved if the problem is known to  +#
C     + be symmetric.                                                        +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C                                                                             #
C     If you don't want to read the FINE PRINT, enter 'no'.                   #
C##############################################################################
      SYMM = .FALSE.                                                           
C##############################################################################
C     If you don't want to read the FINE PRINT, enter 'yes' (strongly         #
C     recommended).                                                           #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + The partial derivatives of some of the PDE and boundary condition    +#
C     + coefficients are required by PDE2D.  These may be calculated         +#
C     + automatically using a finite difference approximation, or supplied   +#
C     + by the user.  Do you want them to be calculated automatically?       +#
C     +                                                                      +#
C     + If you answer 'yes', you will not be asked to supply the derivatives,+#
C     + but there is a small risk that the inaccuracies introduced by the    +#
C     + finite difference approximation may cause the Newton iteration       +#
C     + to converge more slowly or to diverge, especially if low precision   +#
C     + is used.  This risk is very low, however, and since answering 'no'   +#
C     + means you may have to compute many partial derivatives, it is        +#
C     + recommended you answer 'yes' unless you have some reason to believe  +#
C     + there is a problem with the finite difference approximations.        +#
C     +                                                                      +#
C     + If you supply analytic partial derivatives, PDE2D will do some spot  +#
C     + checking and can usually issue a warning if any are supplied         +#
C     + incorrectly.                                                         +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C##############################################################################
      FDIFF = .TRUE.                                                           
C##############################################################################
C     You may calculate one or more integrals (over the entire region) of     #
C     some functions of the solution and its derivatives.  How many integrals #
C     (NINT), if any, do you want to calculate?                               #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + In the FORTRAN program created by the preprocessor, the computed     +#
C     + values of the integrals will be returned in the vector SINT8Z.  If   +#
C     + several iterations or time steps are done, only the last computed    +#
C     + values are saved in SINT8Z (all values are printed).                 +#
C     +                                                                      +#
C     + A limiting value, SLIM8Z(I), for the I-th integral can be set        +#
C     + below in the main program.  The computations will then stop          +#
C     + gracefully whenever SINT8Z(I) > SLIM8Z(I), for any I=1...NINT.       +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C##############################################################################
      NINT =          1                                                        
C##############################################################################
C     You may calculate one or more boundary integrals (over the entire       #
C     boundary) of some functions of the solution and its derivatives.  How   #
C     many boundary integrals (NBINT), if any, do you want to calculate?      #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + In the FORTRAN program created by the preprocessor, the computed     +#
C     + values of the integrals will be returned in the vector BINT8Z.  If   +#
C     + several iterations or time steps are done, only the last computed    +#
C     + values are saved in BINT8Z (all values are printed).                 +#
C     +                                                                      +#
C     + A limiting value, BLIM8Z(I), for the I-th boundary integral can be   +#
C     + set below in the main program.  The computations will then stop      +#
C     + gracefully whenever BINT8Z(I) > BLIM8Z(I), for any I=1...NBINT.      +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C##############################################################################
      NBINT =          0                                                       
      ndel = 0                                                                 
C##############################################################################
C     If you don't want to read the FINE PRINT, enter 'no'.                   #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + Do you want to read the initial conditions from the restart file,    +#
C     + if it exists (and use the conditions supplied above if it does not   +#
C     + exist)?                                                              +#
C     +                                                                      +#
C     + If so, PDE2D will dump the final solution at the end of each run     +#
C     + into a restart file "pde2d.res".  Thus the usual procedure for       +#
C     + using this dump/restart option is to make sure there is no restart   +#
C     + file in your directory left over from a previous job, then the       +#
C     + first time you run this job, the initial conditions supplied above   +#
C     + will be used, but on the second and subsequent runs the restart file +#
C     + from the previous run will be used to define the initial conditions. +#
C     +                                                                      +#
C     + You can do all the "runs" in one program, by setting NPROB > 1.      +#
C     + Each pass through the DO loop, T0,TF,NSTEPS and possibly other       +#
C     + parameters may be varied, by making them functions of IPROB.         +#
C     +                                                                      +#
C     + If the 2D or 3D collocation method is used, the coordinate           +#
C     + transformation should not change between dump and restart.           +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C##############################################################################
      RESTRT = .FALSE.                                                         
C     GRIDID = .FALSE. IF FINITE ELEMENT GRID CHANGES BETWEEN DUMP, RESTART    
      GRIDID = .TRUE.                                                          
C##############################################################################
C     The solution is saved on an NX+1 by NY+1 rectangular grid covering the  #
C     rectangle (XA,XB) x (YA,YB).  Enter values for XA,XB,YA,YB.  These      #
C     variables are usually defaulted.                                        #
C                                                                             #
C     The default is a rectangle which just covers the entire region.         #
C##############################################################################
C        defaults for xa,xb,ya,yb                                              
      call dtdpv   (xgrid,nxgrid,ygrid,nygrid,vxy,nv0,iarc,nt0,xa,xb,ya,       
     &yb)                                                                      
C        DEFINE XA,XB,YA,YB IMMEDIATELY BELOW:
      XA = 1
      XB = 2000
      call dtdpx2(nx,ny,xa,xb,ya,yb,hx8z,hy8z,xout8z,yout8z,npts8z)            
      call dtdpzz(ntf,ideg,isolve,symm,neqn,ii8z,ir8z)                         
      if (iiwk8z.gt.1) ii8z = iiwk8z                                           
      if (irwk8z.gt.1) ir8z = irwk8z                                           
C        *******allocate workspace                                             
      allocate (iwrk8z(ii8z),rwrk8z(ir8z))                                     
C        *******DRAW TRIANGULATION PLOTS (OVER                                 
C        *******RECTANGLE (XA,XB) x (YA,YB))?                                  
      PLOT = .TRUE.                                                            
C        *******call pde solver                                                
      call dtdp2x(xgrid, nxgrid, ygrid, nygrid, ixarc, iyarc, vxy, nv0,        
     &iabc, nt0, iarc, restrt, gridid, neqn, ntf, ideg, isolve, nsteps,        
     &nout, t0, dt, plot, symm, fdiff, itype, nint, nbint, ndel, xd0, yd       
     &0, crankn, noupdt, xbd8z, ybd8z, nbd8z, xout8z, yout8z, uout, inrg       
     &8z, npts8z, ny, tout8z, nsave, iwrk8z, ii8z, rwrk8z, ir8z)               
      deallocate (iwrk8z,rwrk8z)                                               
C        *******read from restart file to array ures8z                         
C      call dtdpr2(1,xres8z,nxp8z,yres8z,nyp8z,ures8z,neqn)                    
C        *******write array ures8z back to restart file                        
C      call dtdpr2(2,xres8z,nxp8z,yres8z,nyp8z,ures8z,neqn)                    
C        *******call user-written postprocessor                                
      call postpr(tout8z,nsave,xout8z,yout8z,inrg8z,nx,ny,uout,neqn)           
C        *******CONTOUR PLOT                                                   
C##############################################################################
C     Enter a value for IVAR, to select the variable to be plotted or         #
C     printed:                                                                #
C         IVAR = 1 means P  (possibly as modified by UPRINT,..)               #
C                2       A1                                                   #
C                3       B1                                                   #
C                4       TMP                                                  #
C                5       A2                                                   #
C                6       B2                                                   #
C##############################################################################
      IVAR =          4                                                        
      ivara8z = mod(ivar-1,3)+1                                                
      ivarb8z = (ivar-1)/3+1                                                   
C##############################################################################
C     If you don't want to read the FINE PRINT, default ISET1,ISET2,ISINC.    #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + The tabular output or plots will be made at times:                   +#
C     +        T(K) = T0 + K*(TF-T0)/NSAVE                                   +#
C     + for    K = ISET1, ISET1+ISINC, ISET1+2*ISINC,..., ISET2              +#
C     + Enter values for ISET1, ISET2 and ISINC.                             +#
C     +                                                                      +#
C     + The default is ISET1=0, ISET2=NSAVE, ISINC=1, that is, the tabular   +#
C     + output or plots will be made at all time values for which the        +#
C     + solution has been saved.                                             +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C##############################################################################
      ISET1 = 0                                                                
      ISET2 = NSAVE                                                            
      ISINC = 1                                                                
C##############################################################################
C     If you don't want to read the FINE PRINT, enter 'no'.                   #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + Do you want to scale the axes on the plot so that the region is      +#
C     + undistorted?  Otherwise the axes will be scaled so that the figure   +#
C     + approximately fills the plot space.                                  +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C##############################################################################
      NODIST = .FALSE.                                                          
C                                                                              
      alow = amin8z(ivar)                                                      
      ahigh = amax8z(ivar)                                                     
C##############################################################################
C     Enter lower (UMIN) and upper (UMAX) bounds for the contour values. UMIN #
C     and UMAX are often defaulted.                                           #
C                                                                             #
C     Labeled contours will be drawn corresponding to the values              #
C                                                                             #
C                  UMIN + S*(UMAX-UMIN),    for S=0.05,0.15,...0.95.          #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + By default, UMIN and UMAX are set to the minimum and maximum values  +#
C     + of the variable to be plotted.  For a common scaling, you may want   +#
C     + to set UMIN=ALOW, UMAX=AHIGH.  ALOW and AHIGH are the minimum and    +#
C     + maximum values over all output points and over all saved time steps  +#
C     + or iterations.                                                       +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C##############################################################################
      UMIN =  361.888                                                         
      UMAX =  584.111                                                          
C##############################################################################
C     Do you want two additional unlabeled contours to be drawn between each  #
C     pair of labeled contours?                                               #
C##############################################################################
      FILLIN = .TRUE.                                                         
C##############################################################################
C     Enter a title, WITHOUT quotation marks.  A maximum of 40 characters     #
C     are allowed.  The default is no title.                                  #
C##############################################################################
      TITLE = ' '                                                              
      TITLE = 'Temperature                             '                       
      call dtdprx(tout8z,nsave,iset1,iset2,isinc)                              
      do 78756 is8z=iset1,iset2,isinc                                          
      call dtdplg(uout(0,0,ivara8z,ivarb8z,is8z),nx,ny,xa,ya,hx8z,hy8z,i       
     &nrg8z,xbd8z,ybd8z,nbd8z,title,umin,umax,nodist,fillin,tout8z(is8z)       
     &)                                                                        
78756 continue                                                                 
C        *******CONTOUR PLOT                                                   
C##############################################################################
C     Enter a value for IVAR, to select the variable to be plotted or         #
C     printed:                                                                #
C         IVAR = 1 means P  (possibly as modified by UPRINT,..)               #
C                2       A1                                                   #
C                3       B1                                                   #
C                4       TMP                                                  #
C                5       A2                                                   #
C                6       B2                                                   #
C##############################################################################
      IVAR =          1                                                        
      ivara8z = mod(ivar-1,3)+1                                                
      ivarb8z = (ivar-1)/3+1                                                   
C##############################################################################
C     If you don't want to read the FINE PRINT, default ISET1,ISET2,ISINC.    #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + The tabular output or plots will be made at times:                   +#
C     +        T(K) = T0 + K*(TF-T0)/NSAVE                                   +#
C     + for    K = ISET1, ISET1+ISINC, ISET1+2*ISINC,..., ISET2              +#
C     + Enter values for ISET1, ISET2 and ISINC.                             +#
C     +                                                                      +#
C     + The default is ISET1=0, ISET2=NSAVE, ISINC=1, that is, the tabular   +#
C     + output or plots will be made at all time values for which the        +#
C     + solution has been saved.                                             +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C##############################################################################
      ISET1 = 0                                                                
      ISET2 = NSAVE                                                            
      ISINC = 1                                                                
C##############################################################################
C     If you don't want to read the FINE PRINT, enter 'no'.                   #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + Do you want to scale the axes on the plot so that the region is      +#
C     + undistorted?  Otherwise the axes will be scaled so that the figure   +#
C     + approximately fills the plot space.                                  +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C##############################################################################
      NODIST = .FALSE.                                                          
C                                                                              
      alow = amin8z(ivar)                                                      
      ahigh = amax8z(ivar)                                                     
C##############################################################################
C     Enter lower (UMIN) and upper (UMAX) bounds for the contour values. UMIN #
C     and UMAX are often defaulted.                                           #
C                                                                             #
C     Labeled contours will be drawn corresponding to the values              #
C                                                                             #
C                  UMIN + S*(UMAX-UMIN),    for S=0.05,0.15,...0.95.          #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + By default, UMIN and UMAX are set to the minimum and maximum values  +#
C     + of the variable to be plotted.  For a common scaling, you may want   +#
C     + to set UMIN=ALOW, UMAX=AHIGH.  ALOW and AHIGH are the minimum and    +#
C     + maximum values over all output points and over all saved time steps  +#
C     + or iterations.                                                       +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C##############################################################################
      UMIN =  33.888 D6                                                        
      UMAX =  56.111 D6                                                        
C##############################################################################
C     Do you want two additional unlabeled contours to be drawn between each  #
C     pair of labeled contours?                                               #
C##############################################################################
      FILLIN = .TRUE.                                                         
C##############################################################################
C     Enter a title, WITHOUT quotation marks.  A maximum of 40 characters     #
C     are allowed.  The default is no title.                                  #
C##############################################################################
      TITLE = ' '                                                              
      TITLE = 'Pressure                             '                       
      call dtdprx(tout8z,nsave,iset1,iset2,isinc)                              
      do 78796 is8z=iset1,iset2,isinc                                          
      call dtdplg(uout(0,0,ivara8z,ivarb8z,is8z),nx,ny,xa,ya,hx8z,hy8z,i       
     &nrg8z,xbd8z,ybd8z,nbd8z,title,umin,umax,nodist,fillin,tout8z(is8z)       
     &)                                                                        
78796 continue                                                                 
C        *******SURFACE PLOT                                                   
C##############################################################################
C     Enter a value for IVAR, to select the variable to be plotted or         #
C     printed:                                                                #
C         IVAR = 1 means P  (possibly as modified by UPRINT,..)               #
C                2       A1                                                   #
C                3       B1                                                   #
C                4       TMP                                                  #
C                5       A2                                                   #
C                6       B2                                                   #
C##############################################################################
      IVAR =          1                                                        
      ivara8z = mod(ivar-1,3)+1                                                
      ivarb8z = (ivar-1)/3+1                                                   
      ISET1 = 1                                                                
      ISET2 = NSAVE                                                            
      ISINC = 1                                                                
C##############################################################################
C     Enter the view latitude, VLAT, and the view longitude, VLON, desired    #
C     for this plot, in degrees.  VLAT and VLON must be between 10 and 80     #
C     degrees; each defaults to 45 degrees.  VLAT and VLON are usually        #
C     defaulted.                                                              #
C##############################################################################
      VLON = 45.0                                                              
      VLAT = 45.0                                                              
C                                                                              
      alow = amin8z(ivar)                                                      
      ahigh = amax8z(ivar)                                                     
C##############################################################################
C     Specify the range (UMIN,UMAX) for the dependent variable axis.  UMIN    #
C     and UMAX are often defaulted.                                           #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + By default, each plot will be scaled to just fit in the plot area.   +#
C     + For a common scaling, you may want to set UMIN=ALOW, UMAX=AHIGH.     +#
C     + ALOW and AHIGH are the minimum and maximum values over all output    +#
C     + points and over all saved time steps or iterations.                  +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C##############################################################################
      UMIN = alow                                                               
      UMAX = ahigh                                                             
C##############################################################################
C     Enter a title, WITHOUT quotation marks.  A maximum of 40 characters     #
C     are allowed.  The default is no title.                                  #
C##############################################################################
      TITLE = ' '                                                              
      TITLE = 'Pressure   '                       
      call dtdprx(tout8z,nsave,iset1,iset2,isinc)                              
      do 78766 is8z=iset1,iset2,isinc                                          
      call dtdpld(xout8z,yout8z,uout(0,0,ivara8z,ivarb8z,is8z),nx,ny,inr       
     &g8z,title,vlon,vlat,umin,umax,tout8z(is8z))                              
78766 continue                                                                 
C        *******SURFACE PLOT                                                   
C##############################################################################
C     Enter a value for IVAR, to select the variable to be plotted or         #
C     printed:                                                                #
C         IVAR = 1 means P  (possibly as modified by UPRINT,..)               #
C                2       A1                                                   #
C                3       B1                                                   #
C                4       TMP                                                  #
C                5       A2                                                   #
C                6       B2                                                   #
C##############################################################################
      IVAR =          4                                                        
      ivara8z = mod(ivar-1,3)+1                                                
      ivarb8z = (ivar-1)/3+1                                                   
      ISET1 = 0                                                                
      ISET2 = NSAVE                                                            
      ISINC = 1                                                                
C##############################################################################
C     Enter the view latitude, VLAT, and the view longitude, VLON, desired    #
C     for this plot, in degrees.  VLAT and VLON must be between 10 and 80     #
C     degrees; each defaults to 45 degrees.  VLAT and VLON are usually        #
C     defaulted.                                                              #
C##############################################################################
      VLON = 45.0                                                              
      VLAT = 45.0                                                              
C                                                                              
      alow = amin8z(ivar)                                                      
      ahigh = amax8z(ivar)                                                     
C##############################################################################
C     Specify the range (UMIN,UMAX) for the dependent variable axis.  UMIN    #
C     and UMAX are often defaulted.                                           #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + By default, each plot will be scaled to just fit in the plot area.   +#
C     + For a common scaling, you may want to set UMIN=ALOW, UMAX=AHIGH.     +#
C     + ALOW and AHIGH are the minimum and maximum values over all output    +#
C     + points and over all saved time steps or iterations.                  +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C##############################################################################
      UMIN = alow                                                               
      UMAX = ahigh                                                             
C##############################################################################
C     Enter a title, WITHOUT quotation marks.  A maximum of 40 characters     #
C     are allowed.  The default is no title.                                  #
C##############################################################################
      TITLE = ' '                                                              
      TITLE = 'Temperature   '                       
      call dtdprx(tout8z,nsave,iset1,iset2,isinc)                              
      do 78767 is8z=iset1,iset2,isinc                                          
      call dtdpld(xout8z,yout8z,uout(0,0,ivara8z,ivarb8z,is8z),nx,ny,inr       
     &g8z,title,vlon,vlat,umin,umax,tout8z(is8z))                              
78767 continue                                                                 
C        *******VECTOR FIELD PLOT                                              
C##############################################################################
C     Enter values for IVARX, IVARY to select the X and Y components of       #
C     the vector to be plotted.                                               #
C         IVARX or IVARY = 1 means P  (possibly as modified by UPRINT,..)     #
C                          2       A1                                         #
C                          3       B1                                         #
C                          4       TMP                                        #
C                          5       A2                                         #
C                          6       B2                                         #
C##############################################################################
      IVARX =          2                                                       
      IVARY =          3                                                       
      ivarxa8z = mod(ivarx-1,3)+1                                              
      ivarxb8z = (ivarx-1)/3+1                                                 
      ivarya8z = mod(ivary-1,3)+1                                              
      ivaryb8z = (ivary-1)/3+1                                                 
C##############################################################################
C     If you don't want to read the FINE PRINT, default ISET1,ISET2,ISINC.    #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + The tabular output or plots will be made at times:                   +#
C     +        T(K) = T0 + K*(TF-T0)/NSAVE                                   +#
C     + for    K = ISET1, ISET1+ISINC, ISET1+2*ISINC,..., ISET2              +#
C     + Enter values for ISET1, ISET2 and ISINC.                             +#
C     +                                                                      +#
C     + The default is ISET1=0, ISET2=NSAVE, ISINC=1, that is, the tabular   +#
C     + output or plots will be made at all time values for which the        +#
C     + solution has been saved.                                             +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C##############################################################################
      ISET1 = 1                                                                
      ISET2 = NSAVE                                                            
      ISINC = 1                                                                
C##############################################################################
C     If you don't want to read the FINE PRINT, enter 'no'.                   #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + Do you want to scale the axes on the plot so that the region is      +#
C     + undistorted?  Otherwise the axes will be scaled so that the figure   +#
C     + approximately fills the plot space.                                  +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C##############################################################################
      NODIST = .FALSE.                                                          
C                                                                              
      a1mag = max(abs(amin8z(ivarx)),abs(amax8z(ivarx)))                       
      a2mag = max(abs(amin8z(ivary)),abs(amax8z(ivary)))                       
C##############################################################################
C     For the purpose of scaling the arrows, the ranges of the two components #
C     of the vector are assumed to be (-VR1MAG,VR1MAG) and (-VR2MAG,VR2MAG).  #
C     Enter values for VR1MAG and VR2MAG.  VR1MAG and VR2MAG are often        #
C     defaulted.                                                              #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + By default, VR1MAG and VR2MAG are the maxima of the absolute values  +#
C     + of the first and second components.  For a common scaling, you may   +#
C     + want to set VR1MAG=A1MAG, VR2MAG=A2MAG.  A1MAG, A2MAG are the        +#
C     + maxima of the absolute values over all output points and over all    +#
C     + saved time steps or iterations.                                      +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C##############################################################################
      VR1MAG = 0.                                                            
      VR2MAG = 0.                                                        
C##############################################################################
C     Enter a title, WITHOUT quotation marks.  A maximum of 40 characters     #
C     are allowed.  The default is no title.                                  #
C##############################################################################
      TITLE = ' '                                                              
      TITLE = '(QR,QZ)                                 '                       
      call dtdprx(tout8z,nsave,iset1,iset2,isinc)                              
      do 78757 is8z=iset1,iset2,isinc                                          
      call dtdpla(uout(0,0,ivarxa8z,ivarxb8z,is8z),uout(0,0,ivarya8z,iva       
     &ryb8z,is8z),nx,ny,xa,ya,hx8z,hy8z,inrg8z,xbd8z,ybd8z,nbd8z,title,v       
     &r1mag,vr2mag,nodist,tout8z(is8z))                                        
78757 continue                                                                 
78755 continue                                                                 
      call endgks                                                              
      stop                                                                     
      end                                                                      
                                                                               
                                                                               
      subroutine tran8z(p,q,x,y)                                               
      implicit double precision (a-h,o-z)                                      
      x = p                                                                    
      y = q                                                                    
      return                                                                   
      end                                                                      
                                                                               
                                                                               
      subroutine xy8z(i8z,iarc8z,s,x,y,s0,sf)                                  
      implicit double precision (a-h,o-z)                                      
      dimension pxy(2,1000)                                                    
      common/parm8z/ pi,h1    ,h2    ,h3    ,h4    ,w1    ,w2   ,F0    
      x = 0.0                                                                  
      y = 0.0                                                                  
      return                                                                   
      end                                                                      
                                                                               
                                                                               
      subroutine dis8z(x,y,ktri,triden,shape)                                  
      implicit double precision (a-h,o-z)                                      
      logical adapt                                                            
      common/parm8z/ pi,h1    ,h2    ,h3    ,h4    ,w1    ,w2   ,F0    
C##############################################################################
C     Enter a FORTRAN expression for TRIDEN(X,Y), which controls the          #
C     grading of the triangulation.  TRIDEN should be largest where the       #
C     triangulation is to be most dense.  The default is TRIDEN(X,Y)=1.0      #
C     (a uniform triangulation).                                              #
C                                                                             #
C     TRIDEN may also be a function of the initial triangle number KTRI.      #
C##############################################################################
      TRIDEN = 1.0 
      IF (X.le.1000.) TRIDEN = 20.
C##############################################################################
C     If you don't want to read the FINE PRINT, enter 'no'.                   #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + Do you want the triangulation to be graded adaptively?               +#
C     +                                                                      +#
C     + If you answer "yes", make sure there is no "pde2d.adp" file in the   +#
C     + working directory the first time you run the program, then run the   +#
C     + program two or more times, possibly increasing NTF each time.  On    +#
C     + the first run, the triangulation will be graded as guided by         +#
C     + TRIDEN(X,Y), but on each subsequent run, information output by the   +#
C     + previous run to "pde2d.adp" will be used to guide the grading of the +#
C     + new triangulation.                                                   +#
C                                                                            +#
C     + If ADAPT=.TRUE., after each run a file "pde2d.adp" is written        +#
C     + which tabulates the values of the magnitude of the gradient of the   +#
C     + solution (at the last time step or iteration) at an output NXP8Z by  +#
C     + NYP8Z grid of points (NXP8Z and NYP8Z are set to 101 in a PARAMETER  +#
C     + statement in the main program, so they can be changed if desired).   +#
C     + If NEQN > 1, a normalized average of the gradients of the NEQN       +#
C     + solution components is used.                                         +#
C     +                                                                      +#
C     + You can do all the "runs" in one program, by setting NPROB > 1.      +#
C     + Each pass through the DO loop, PDE2D will read the gradient values   +#
C     + output the previous pass.  If RESTRT=.TRUE., GRIDID=.FALSE., and     +#
C     + T0,TF are incremented each pass through the DO loop, it is possible  +#
C     + in this way to solve a time-dependent problem with an adaptive,      +#
C     + moving, grid.                                                        +#
C     +                                                                      +#
C     + Increase the variable EXAG from its default value of 1.5 if you want +#
C     + to exaggerate the grading of an adaptive triangulation (make it less +#
C     + uniform).  EXAG should normally not be larger than about 2.0.        +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C                                                                             #
C     If you don't want to read the FINE PRINT, enter 'no'.                   #
C##############################################################################
      ADAPT = .FALSE.                                                          
C                                                                              
      EXAG = 1.5                                                               
      if (adapt) triden = dtdpgr2(x,y,triden**(1.0/exag))**exag                
C##############################################################################
C     If you don't want to read the FINE PRINT, default SHAPE.                #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + Enter a FORTRAN expression for SHAPE(X,Y), which controls the        +#
C     + approximate shape of the triangles.  The triangulation refinement    +#
C     + will proceed with the goal of generating triangles with an average   +#
C     + height to width ratio of approximately SHAPE(X,Y) near the point     +#
C     + (X,Y).  SHAPE must be positive.  The default is SHAPE(X,Y)=1.0.      +#
C     +                                                                      +#
C     + SHAPE may also be a function of the initial triangle number KTRI.    +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C##############################################################################
      SHAPE = 1.0                                                              
      return                                                                   
      end                                                                      
                                                                               
                                                                               
      subroutine pdes8z(yd8z,i8z,j8z,kint8z,idel8z,jdel8z,x,y,ktri,t)          
      implicit double precision (a-h,o-z)                                      
      parameter (neqnmx=  99)                                                  
      parameter (ndelmx=  20)                                                  
C        un8z(1,I),un8z(2,I),un8z(3,I) hold the (rarely used) values           
C        of UI,UIx,UIy from the previous iteration or time step                
      common/dtdp4/un8z(3,neqnmx),uu8z(3,neqnmx)                               
      common/dtdp17/normx,normy,iarc                                           
      double precision normx,normy,delamp(ndelmx,neqnmx)                       
      common/parm8z/ pi,h1    ,h2    ,h3    ,h4    ,w1    ,w2   ,F0    
      P  = uu8z(1, 1)                                                          
      Px = uu8z(2, 1)                                                          
      Py = uu8z(3, 1)                                                          
      Pnorm = Px*normx + Py*normy                                              
      TMP  = uu8z(1, 2)                                                        
      TMPx = uu8z(2, 2)                                                        
      TMPy = uu8z(3, 2)                                                        
      TMPnorm = TMPx*normx + TMPy*normy                                        
                          if (i8z.eq.0) then                                   
      yd8z = 0.0                                                               
C##############################################################################
C     Enter FORTRAN expressions for the functions whose integrals are to be   #
C     calculated and printed.  They may be functions of                       #
C                                                                             #
C                X,Y,P,Px,Py,TMP,TMPx,TMPy and (if applicable) T              #
C                                                                             #
C     The integrals may also contain references to the initial triangle       #
C     number KTRI.                                                            #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + If you only want to integrate a function over part of the region,    +#
C     + define that function to be zero in the rest of the region.           +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C##############################################################################
C                                                  INTEGRAL DEFINED            
      if (kint8z.eq.    1) yd8z =                                              
     & X*TMP                                                                   
C##############################################################################
C     Enter FORTRAN expressions for the functions whose integrals are to be   #
C     calculated and printed.  They may be functions of                       #
C                                                                             #
C                X,Y,P,Px,Py,TMP,TMPx,TMPy and (if applicable) T              #
C                                                                             #
C     The components (NORMx,NORMy) of the unit outward normal vector, and the #
C     initial triangle number KTRI, and the boundary arc number IARC may also #
C     be referenced.  You can also reference the normal derivatives Pnorm,    #
C     TMPnorm.                                                                #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + CURVED interior interface arcs are considered part of the boundary,  +#
C     + for the boundary integral computations, ONLY IF they have arc numbers+#
C     + in the range 8000-8999.  In this case, since an interface arc is     +#
C     + considered to be a boundary for both of the subregions it separates, +#
C     + the boundary integral will be computed twice on each curved          +#
C     + interface arc, once with (NORMx,NORMy) defined in each direction.    +#
C     +                                                                      +#
C     + If you only want to integrate a function over part of the boundary,  +#
C     + define that function to be zero on the rest of the boundary.  You    +#
C     + can examine the point (X,Y) to determine if it is on the desired     +#
C     + boundary segment, or the boundary arc number IARC, or the initial    +#
C     + triangle number KTRI (if INTRI=3),                                   +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C##############################################################################
C                                                  BND. INTEGRAL1 DEFINED      
C      if (kint8z.eq.-1) yd8z =                                                 
C     & 0.0                           
                          else                                                 
C##############################################################################
C     Now enter FORTRAN expressions to define the PDE coefficients, which     #
C     may be functions of                                                     #
C                                                                             #
C                       X,Y,T,P,Px,Py,TMP,TMPx,TMPy                           #
C                                                                             #
C     They may also be functions of the initial triangle number KTRI.         #
C                                                                             #
C     Recall that the PDEs have the form                                      #
C                                                                             #
C      C11*d(P)/dT + C12*d(TMP)/dT = d/dX*A1 + d/dY*B1 - F1                   #
C      C21*d(P)/dT + C22*d(TMP)/dT = d/dX*A2 + d/dY*B2 - F2                   #
C                                                                             #
C##############################################################################
      PERMK = 1.D-15
      IF (Y.GT.1900.) PERMK = 1.D-14
      COND = 2.0 
      G = 9.81
      CPR = 880
      RHOR = 2700.
C ***********************************************
      CALL IAPWS(TMP,P,RHO,CPW,VMU)
      EPS = 0.1
      HC1 = EPS*RHO*CPW + (1.-EPS)*RHOR*CPR
      HC2 = RHO*CPW
C   Equation (2)
      QR = -PERMK/VMU*Px
      QZ = -PERMK/VMU*(Py + RHO*G) 
                if (j8z.eq.0) then                                             
      yd8z = 0.0                                                               
C                                                  C(1,1) DEFINED              
      if (i8z.eq. -101) yd8z =                                                 
     & 0
C                                                  C(1,2) DEFINED              
      if (i8z.eq. -102) yd8z =                                                 
     & 0
C                                                  F1 DEFINED                  
      if (i8z.eq.    1) yd8z =                                                 
     & 0                                                                       
C                                                  A1 DEFINED                  
      if (i8z.eq.    2) yd8z =                                                 
     & X*RHO*QR                                                                 
C                                                  B1 DEFINED                  
      if (i8z.eq.    3) yd8z =                                                 
     & X*RHO*QZ                                                                 
C                                                  C(2,1) DEFINED              
      if (i8z.eq. -201) yd8z =                                                 
     & 0                                                                       
C                                                  C(2,2) DEFINED              
      if (i8z.eq. -202) yd8z =                                                 
     & X*HC1                                                                  
C                                                  F2 DEFINED                  
      if (i8z.eq.    4) yd8z =                                                 
     & X*HC2*(QR*TMPx+QZ*TMPy)                                                 
C                                                  A2 DEFINED                  
      if (i8z.eq.    5) yd8z =                                                 
     & X*COND*TMPx                                                 
C                                                  B2 DEFINED                  
      if (i8z.eq.    6) yd8z =                                                 
     & X*COND*TMPy                                                 
                else                                                           
                endif                                                          
                          endif                                                
      return                                                                   
      end                                                                      
                                                                               
                                                                               
      function u8z(i8z,x,y,ktri,t0)                                            
      implicit double precision (a-h,o-z)                                      
      common/parm8z/ pi,h1    ,h2    ,h3    ,h4    ,w1    ,w2   ,F0    
      u8z = 0.0                                                                
C##############################################################################
C     Now the initial values must be defined using FORTRAN expressions.       #
C     They may be functions of X and Y, and of the initial triangle number    #
C     KTRI.  They may also reference the initial time T0.                     #
C##############################################################################
C                                                  P0 DEFINED                  
      if (i8z.eq.    1) u8z =                                                  
     & 30.D6                                                                  
C                                                  TMP0 DEFINED                
      if (i8z.eq.    2) u8z =                                                  
     & 278.                                                                    
      return                                                                   
      end                                                                      
                                                                               
                                                                               
      function fb8z(i8z,iarc8z,ktri,s,x,y,t)                                   
      implicit double precision (a-h,o-z)                                      
      common/parm8z/ pi,h1    ,h2    ,h3    ,h4    ,w1    ,w2   ,F0    
      fb8z = 0.0                                                               
C        NO BOUNDARY CONDITIONS DEFINED ON NEGATIVE ARCS.                      
C        TO ADD BCs FOR NEGATIVE ARCS, USE BLOCK BELOW AS MODEL                
      IARC = -3                                                                 
                if (iarc8z.eq.iarc) then                                       
C##############################################################################
C     Enter FORTRAN expressions to define FB1,FB2 on this arc.  They may      #
C     be functions of X,Y and (if applicable) T.                              #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + These functions may also reference the initial triangle number KTRI, +#
C     + and the arc parameter S (S = P or Q when INTRI=2).                   +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C                                                                             #
C     Recall that fixed boundary conditions have the form                     #
C                                                                             #
C                    P = FB1                                                  #
C                    TMP = FB2                                                #
C                                                                             #
C##############################################################################
C                                                  FB1 DEFINED                 
      if (i8z.eq.    1) fb8z =                                                 
     & 30.D6
C                                                  FB2 DEFINED                 
      if (i8z.eq.    2) fb8z =                                                 
     & 278.
                endif                                                          
      return                                                                   
      end                                                                      
                                                                               
                                                                               
      subroutine gb8z(gd8z,i8z,j8z,iarc8z,ktri,s,x,y,t)                        
      implicit double precision (a-h,o-z)                                      
      parameter (neqnmx=  99)                                                  
C        un8z(1,I),un8z(2,I),un8z(3,I) hold the (rarely used) values           
C        of UI,UIx,UIy from the previous iteration or time step.               
C        (normx,normy) is the (rarely used) unit outward normal vector         
      common/dtdp4/un8z(3,neqnmx),uu8z(neqnmx,3)                               
      common/dtdp17/normx,normy,ibarc8z                                        
      common/dtdp49/bign8z                                                     
      double precision normx,normy                                             
      common/parm8z/ pi,h1    ,h2    ,h3    ,h4    ,w1    ,w2   ,F0    
      zero(f8z) = bign8z*f8z                                                   
      P  = uu8z( 1,1)                                                          
      Px = uu8z( 1,2)                                                          
      Py = uu8z( 1,3)                                                          
      TMP  = uu8z( 2,1)                                                        
      TMPx = uu8z( 2,2)                                                        
      TMPy = uu8z( 2,3)                                                        
      if (j8z.eq.0) gd8z = 0.0
C##############################################################################
C     Enter the arc number (IARC) of a boundary arc with positive arc number. #
C##############################################################################
      IARC =          1                                                        
                if (iarc8z.eq.iarc) then                                       
C##############################################################################
C     Enter FORTRAN expressions to define the following free boundary         #
C     condition functions on this arc.  They may be functions of              #
C                                                                             #
C                X,Y,P,Px,Py,TMP,TMPx,TMPy and (if applicable) T              #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + These functions may also reference the components (NORMx,NORMy) of   +#
C     + the unit outward normal vector, the initial triangle number KTRI,    +#
C     + and the arc parameter S (S = P or Q when INTRI=2).                   +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C                                                                             #
C     Recall that free boundary conditions have the form                      #
C                                                                             #
C                     A1*nx+B1*ny = GB1                                       #
C                     A2*nx+B2*ny = GB2                                       #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + Note that 'fixed' boundary conditions:                               +#
C     +                   Ui = FBi(X,Y,[T])                                  +#
C     + can be expressed as 'free' boundary conditions in the form:          +#
C     +                  GBi = zero(Ui-FBi(X,Y,[T]))                         +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C##############################################################################
                          if (j8z.eq.0) then                                   
C                                                  GB1 DEFINED                 
      if (i8z.eq.    1) gd8z =                                                 
     & 0                                                 
C                                                  GB2 DEFINED                 
      qmin = 0.05
      qmax = 5.0
      c = 500.
      F0 = qmin + (qmax-qmin)*exp(-x**2/(2*c**2))
      if (i8z.eq.    2) gd8z =                                            
     & X*F0                                                           
                          else                                                 
                          endif                                                
                endif                                                          
C##############################################################################
C     Enter the arc number (IARC) of a boundary arc with positive arc number. #
C##############################################################################
      IARC =          2                                                        
                if (iarc8z.eq.iarc) then                                       
C##############################################################################
C     Enter FORTRAN expressions to define the following free boundary         #
C     condition functions on this arc.  They may be functions of              #
C                                                                             #
C                X,Y,P,Px,Py,TMP,TMPx,TMPy and (if applicable) T              #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + These functions may also reference the components (NORMx,NORMy) of   +#
C     + the unit outward normal vector, the initial triangle number KTRI,    +#
C     + and the arc parameter S (S = P or Q when INTRI=2).                   +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C                                                                             #
C     Recall that free boundary conditions have the form                      #
C                                                                             #
C                     A1*nx+B1*ny = GB1                                       #
C                     A2*nx+B2*ny = GB2                                       #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + Note that 'fixed' boundary conditions:                               +#
C     +                   Ui = FBi(X,Y,[T])                                  +#
C     + can be expressed as 'free' boundary conditions in the form:          +#
C     +                  GBi = zero(Ui-FBi(X,Y,[T]))                         +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C##############################################################################
                          if (j8z.eq.0) then                                   
C                                                  GB1 DEFINED                 
      if (i8z.eq.    1) gd8z =                                                 
     & 0                                                               
C                                                  GB2 DEFINED                 
      if (i8z.eq.    2) gd8z =                                            
     & 0                                                      
                          else                                                 
                          endif                                                
                endif                                                            
C##############################################################################
C     Enter the arc number (IARC) of a boundary arc with positive arc number. #
C##############################################################################
      IARC =          4                                                        
                if (iarc8z.eq.iarc) then                                       
C##############################################################################
C     Enter FORTRAN expressions to define the following free boundary         #
C     condition functions on this arc.  They may be functions of              #
C                                                                             #
C                X,Y,P,Px,Py,TMP,TMPx,TMPy and (if applicable) T              #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + These functions may also reference the components (NORMx,NORMy) of   +#
C     + the unit outward normal vector, the initial triangle number KTRI,    +#
C     + and the arc parameter S (S = P or Q when INTRI=2).                   +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C                                                                             #
C     Recall that free boundary conditions have the form                      #
C                                                                             #
C                     A1*nx+B1*ny = GB1                                       #
C                     A2*nx+B2*ny = GB2                                       #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + Note that 'fixed' boundary conditions:                               +#
C     +                   Ui = FBi(X,Y,[T])                                  +#
C     + can be expressed as 'free' boundary conditions in the form:          +#
C     +                  GBi = zero(Ui-FBi(X,Y,[T]))                         +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C##############################################################################
                          if (j8z.eq.0) then                                   
C                                                  GB1 DEFINED                 
      if (i8z.eq.    1) gd8z =                                                 
     & 0                                                   
C                                                  GB2 DEFINED                 
      if (i8z.eq.    2) gd8z =                                            
     & 0                                                       
                          else                                                 
                          endif                                                
                endif                                                            
      return                                                                   
      end                                                                      
                                                                               
                                                                               
      subroutine pmod8z(x,y,ktri,t,a,b)                                        
      implicit double precision (a-h,o-z)                                      
      parameter (neqnmx=  99)                                                  
      common/dtdp4/un8z(3,neqnmx),uu8z(3,neqnmx)                               
      common/dtdp6/upr8z(3,neqnmx),uab8z(3,neqnmx)                             
      common/dtdp9/uprint(neqnmx),aprint(neqnmx),bprint(neqnmx)                
      common/dtdp14/sint(20),bint(20),slim8z(20),blim8z(20)                    
      common/parm8z/ pi,h1    ,h2    ,h3    ,h4    ,w1    ,w2   ,F0    
      P  = uu8z(1, 1)                                                          
      Px = uu8z(2, 1)                                                          
      Py = uu8z(3, 1)                                                          
      a1 = upr8z(2, 1)                                                         
      b1 = upr8z(3, 1)                                                         
      TMP  = uu8z(1, 2)                                                        
      TMPx = uu8z(2, 2)                                                        
      TMPy = uu8z(3, 2)                                                        
      a2 = upr8z(2, 2)                                                         
      b2 = upr8z(3, 2)                                                         
C##############################################################################
C     If you don't want to read the FINE PRINT, default all of the following  #
C     variables.                                                              #
C                                                                             #
C     +++++++++++++++ THE "FINE PRINT" (CAN USUALLY BE IGNORED) ++++++++++++++#
C     + Normally, PDE2D saves the values of P,A1,B1,TMP,A2,B2 at the         +#
C     + output points.  If different variables are to be saved (for later    +#
C     + printing or plotting) the following functions can be used to         +#
C     + re-define the output variables:                                      +#
C     +    define UPRINT(1) to replace P                                     +#
C     +           APRINT(1)            A1                                    +#
C     +           BPRINT(1)            B1                                    +#
C     +           UPRINT(2)            TMP                                   +#
C     +           APRINT(2)            A2                                    +#
C     +           BPRINT(2)            B2                                    +#
C     + Each function may be a function of                                   +#
C     +                                                                      +#
C     +    X,Y,P,Px,Py,A1,B1,TMP,TMPx,TMPy,A2,B2 and (if applicable) T       +#
C     +                                                                      +#
C     + Each may also be a function of the initial triangle number KTRI and  +#
C     + the integral estimates SINT(1),...,BINT(1),...                       +#
C     +                                                                      +#
C     + The default for each variable is no change, for example, UPRINT(1)   +#
C     + defaults to P.  Enter FORTRAN expressions for each of the            +#
C     + following functions (or default).                                    +#
C     ++++++++++++++++++++++++++ END OF "FINE PRINT" +++++++++++++++++++++++++#
C##############################################################################
C        DEFINE UPRINT(*),APRINT(*),BPRINT(*) HERE:
      CALL IAPWS(TMP,P,RHO,CPW,VMU)
      APRINT(1) = A1/RHO/X
      BPRINT(1) = B1/RHO/X
      return                                                                   
      end                                                                      
C        dummy routines                                                        
      function axis8z(i8z,x,y,z,ical8z)                                        
      implicit double precision (a-h,o-z)                                      
      axis8z = 0                                                               
      return                                                                   
      end                                                                      
                                                                               
                                                                               
      subroutine postpr(tout,nsave,xout,yout,inrg,nx,ny,uout,neqn)             
      implicit double precision (a-h,o-z)                                      
      dimension xout(0:nx,0:ny),yout(0:nx,0:ny),tout(0:nsave)                  
      dimension inrg(0:nx,0:ny),uout(0:nx,0:ny,3,neqn,0:nsave)                 
      common/parm8z/ pi,h1    ,h2    ,h3    ,h4    ,w1    ,w2   ,F0    
      common /dtdp27/ itask,npes,icomm                                         
      common /dtdp46/ eps8z,cgtl8z,npmx8z,itype,near8z                         
      data lun,lud/0,47/                                                       
      if (itask.gt.0) return                                                   
C     UOUT(I,J,IDER,IEQ,L) = U_IEQ, if IDER=1                                  
C                            A_IEQ, if IDER=2                                  
C                            B_IEQ, if IDER=3                                  
C       (possibly as modified by UPRINT,..)                                    
C       at the point (XOUT(I,J) , YOUT(I,J))                                   
C       at time/iteration TOUT(L).                                             
      open (19,file='axi_benchmark_k1_k2.out')
      write (19,9) tout(nsave)
    9 format (' T = ',E15.5,//,'          R    ','          Z    ',
     &  '          P    ','          TMP  ')
      do i=0,nx
      do j=0,ny
         write (19,10) xout(i,j),yout(i,j),uout(i,j,1,1,nsave),
     &                                     uout(i,j,1,2,nsave)
   10 format (4E15.5)
      enddo
      enddo 
      close (19)
C     INRG(I,J) = 1 if this point is in R                                      
C               = 0 otherwise                                                  
C       ******* ADD POSTPROCESSING CODE HERE:                                  
C                                                                              
C     SET MATLAB_PLOTS = 1 TO CREATE MATLAB PLOTFILES                          
C       pde2d.m, pde2d.rdm                                                     
      MATLAB_PLOTS = 1                                                         
      if (MATLAB_PLOTS .eq. 0) return                                          
      if (lun.eq.0) then                                                       
         lun = 46                                                              
         open (lun,file='pde2d.m')                                             
         open (lud,file='pde2d.rdm')                                           
         write (lun,*) 'fid = fopen(''pde2d.rdm'');'                           
      endif                                                                    
      do 78753 l=0,nsave                                                       
         if (tout(l).ne.dtdplx(2)) nsave0 = l                                  
78753 continue                                                                 
      write (lud,78754) nsave0                                                 
      write (lud,78754) neqn                                                   
      write (lud,78754) nx                                                     
      write (lud,78754) ny                                                     
78754 format (i8)                                                              
      do 78756 i=0,nx                                                          
      do 78755 j=0,ny                                                          
         write (lud,78762) xout(i,j),yout(i,j)                                 
78755 continue                                                                 
78756 continue                                                                 
      do 78761 l=0,nsave0                                                      
         write (lud,78762) tout(l)                                             
         do 78760 ieq=1,neqn                                                   
         do 78759 ider=1,3                                                     
         do 78758 i=0,nx                                                       
         do 78757 j=0,ny                                                       
            if (inrg(i,j).eq.1) then                                           
               write (lud,78762) uout(i,j,ider,ieq,l)                          
            else                                                               
               write (lud,78763)                                               
            endif                                                              
78757    continue                                                              
78758    continue                                                              
78759    continue                                                              
78760    continue                                                              
78761 continue                                                                 
78762 format (e16.8)                                                           
78763 format ('NaN')                                                           
C       ******* WRITE pde2d.m                                                  
      call mtdp2dg(itype,lun)                                                  
      return                                                                   
      end


      SUBROUTINE IAPWS(TMP,P,RHO,CPW,VMU)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPIX=800,NPIY=80)
      DIMENSION PIX(1:NPIX,0:NPIY,3)
      SAVE PIX
      LOGICAL READ 
      DATA READ /.FALSE./
      IF (.NOT.READ) THEN
C                        READ PIXEL FILE THE FIRST CALL ONLY
         OPEN (11,FILE='EOS_pure_water.dat')

         DO I=1,NPIX
         DO J=0,NPIY
            READ (11,*) XX,YY,PIX(I,J,1),PIX(I,J,2),PIX(I,J,3)
         END DO
         END DO 
         CLOSE (11)
         READ = .TRUE. 
      ENDIF
C                        ON ALL CALLS, USE BILINEAR INTERPOLATION 
C                          TO GET PIXEL VALUE AT (X,Y) 
C                        FOUR NEAREST NEIGHBORS TO (X,Y) ARE 
C                          (II,JJ)/NPIX,    (II+1,JJ)/NPIX,
C                          (II,JJ+1)/NPIX,  (II+1,JJ+1)/NPIX.
      HX = 800.D0/NPIX
      HY = 80.D0/NPIY
      X = min(TMP-273.d0,800.d0)
      X = max(X,5.d0)
      Y = min(P/1.d6,80.d0) 
      Y = max(Y,5.d0)
      II = MIN(INT(X/HX),NPIX-1)
      JJ = MIN(INT(Y/HY),NPIY-1)
      XX = X/HX - II
      YY = Y/HY - JJ
      DO K=1,3
         U1 = PIX(II,JJ,K)   + XX*(PIX(II+1,JJ,K)  -PIX(II,JJ,K))
         U2 = PIX(II,JJ+1,K) + XX*(PIX(II+1,JJ+1,K)-PIX(II,JJ+1,K))
         U = U1 + YY*(U2-U1)
         IF (K.EQ.1) RHO = U 
         IF (K.EQ.2) CPW = U
         IF (K.EQ.3) VMU = U
      END DO
      RETURN
      END
