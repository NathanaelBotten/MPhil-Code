! The following code was modified to be updated to the Fortran 90 standard
! February 2007
!
! CHANGES MADE
! ============
!
! Continuation characters: + at the beginning of lines was changed
! to & at the end of the line to be continued and the beginning of
! the continuation line.
! 
! Fixed form formatting: Removed the 'code starts in the 7th column' aspect
! that F77 requires, code now starts in 5th column, with each subsequent layer
! of code starting 3 columns further in.
!
! Commenting: Replaced all 1st column C with !, as C is not recognised as
! comment character in F90. Also introduced in line commenting with
! code sections that were updated
!
! Interface blocks: Commented out 'EXTERNAL' statement and replaced it with an 
! interface block that provides the complier the ability to check the function
! calls are correct.
!
! DOUBLE PRECISION changes: all declarations of DOUBLE PRECISION variables was
! changed to real(dp), and the introduction of 'use mykinds' statement.
!
! Modules: wrapped all of the subroutines DCUHRE uses in to a single 
! module, DCUHRE90, placed the test function into its own module (commented out
! below as to make DCUHRE useable with the myIntegrate interface).
!
! NUMFUN changes: Originally, DCUHRE accepted a NUMFUN arguement, which
! defined how many functions were to be integrated.  However, the intent
! was to modify DCUHRE to integrate only one function, so NUMFUN was
! removed from the input variables and defined as a parameter with a value
! of 1 in DCUHRE. The other subroutines are unchanged, accepting NUMFUN as 
! an input arguement, without allowing it to be changed.
!
! FUNSUB changes: In the original DCUHRE, FUNSUB was a subroutine designed to
! accept NUMFUN functions in NDIM dimensions, which were to be integrated over 
! the same region, with NUMFUN, X, NDIM and F passed as arguments. This
! was changed to make FUNSUB a function, taking one argument, the array X. This was 
! done to ensure compatibility with the myIntegrate interface. All code that
! relied on calling FUNSUB needed to be changed. All original code
! remains, but is commented out.

module DCUHRE90
  use nrtype

contains
  
  SUBROUTINE DCUHRE(NDIM,A,B,MINPTS,MAXPTS,FUNSUB,EPSABS,  &
       &            EPSREL,KEY,NW,RESTAR,RESULT,ABSERR,NEVAL,IFAIL,  &
       &            WORK)
!
!***BEGIN PROLOGUE DCUHRE
!***DATE WRITTEN   900116   (YYMMDD)
!***REVISION DATE  900116   (YYMMDD)
!***CATEGORY NO. H2B1A1
!***AUTHOR
!            Jarle Berntsen, The Computing Centre,
!            University of Bergen, Thormohlens gt. 55,
!            N-5008 Bergen, Norway
!            Phone..  47-5-544055
!            Email..  jarle@eik.ii.uib.no
!            Terje O. Espelid, Department of Informatics,
!            University of Bergen, Thormohlens gt. 55,
!            N-5008 Bergen, Norway
!            Phone..  47-5-544180
!            Email..  terje@eik.ii.uib.no
!            Alan Genz, Computer Science Department, Washington State
!            University, Pullman, WA 99163-2752, USA
!            Email..  acg@eecs.wsu.edu
!***KEYWORDS automatic multidimensional integrator,
!            n-dimensional hyper-rectangles,
!            general purpose, global adaptive
!***PURPOSE  The routine calculates an approximation to a given
!            vector of definite integrals
!
!      B(1) B(2)     B(NDIM)
!     I    I    ... I       (F ,F ,...,F      ) DX(NDIM)...DX(2)DX(1),
!      A(1) A(2)     A(NDIM)  1  2      NUMFUN
!
!       where F = F (X ,X ,...,X    ), I = 1,2,...,NUMFUN.
!              I   I  1  2      NDIM
!
!            hopefully satisfying for each component of I the following
!            claim for accuracy:
!            ABS(I(K)-RESULT(K)).LE.MAX(EPSABS,EPSREL*ABS(I(K)))
!***DESCRIPTION Computation of integrals over hyper-rectangular
!            regions.
!            DCUHRE is a driver for the integration routine
!            DADHRE, which repeatedly subdivides the region
!            of integration and estimates the integrals and the
!            errors over the subregions with greatest
!            estimated errors until the error request
!            is met or MAXPTS function evaluations have been used.
!
!            For NDIM = 2 the default integration rule is of 
!            degree 13 and uses 65 evaluation points.
!            For NDIM = 3 the default integration rule is of 
!            degree 11 and uses 127 evaluation points.
!            For NDIM greater then 3 the default integration rule
!            is of degree 9 and uses NUM evaluation points where
!              NUM = 1 + 4*2*NDIM + 2*NDIM*(NDIM-1) + 4*NDIM*(NDIM-1) +
!                    4*NDIM*(NDIM-1)*(NDIM-2)/3 + 2**NDIM
!            The degree 9 rule may also be applied for NDIM = 2
!            and NDIM = 3.
!            A rule of degree 7 is available in all dimensions.
!            The number of evaluation
!            points used by the degree 7 rule is
!              NUM = 1 + 3*2*NDIM + 2*NDIM*(NDIM-1) + 2**NDIM
!
!            When DCUHRE computes estimates to a vector of
!            integrals, all components of the vector are given
!            the same treatment. That is, I(F ) and I(F ) for
!                                            J         K 
!            J not equal to K, are estimated with the same
!            subdivision of the region of integration.
!            For integrals with enough similarity, we may save
!            time by applying DCUHRE to all integrands in one call.
!            For integrals that varies continuously as functions of
!            some parameter, the estimates produced by DCUHRE will
!            also vary continuously when the same subdivision is
!            applied to all components. This will generally not be
!            the case when the different components are given
!            separate treatment.
!
!            On the other hand this feature should be used with
!            caution when the different components of the integrals
!            require clearly different subdivisions.
!
!   ON ENTRY
!
!     NDIM   Integer.
!            Number of variables. 1 < NDIM <=  15.
!     NUMFUN Integer.
!            Number of components of the integral.
!     A      Real array of dimension NDIM.
!            Lower limits of integration.
!     B      Real array of dimension NDIM.
!            Upper limits of integration.
!     MINPTS Integer.
!            Minimum number of function evaluations.
!     MAXPTS Integer.
!            Maximum number of function evaluations.
!            The number of function evaluations over each subregion
!            is NUM.
!            If (KEY = 0 or KEY = 1) and (NDIM = 2) Then
!              NUM = 65
!            Elseif (KEY = 0 or KEY = 2) and (NDIM = 3) Then
!              NUM = 127
!            Elseif (KEY = 0 and NDIM > 3) or (KEY = 3) Then
!              NUM = 1 + 4*2*NDIM + 2*NDIM*(NDIM-1) + 4*NDIM*(NDIM-1) +
!                    4*NDIM*(NDIM-1)*(NDIM-2)/3 + 2**NDIM
!            Elseif (KEY = 4) Then
!              NUM = 1 + 3*2*NDIM + 2*NDIM*(NDIM-1) + 2**NDIM
!            MAXPTS >= 3*NUM and MAXPTS >= MINPTS
!            For 3 < NDIM < 13 the minimum values for MAXPTS are:
!             NDIM =    4   5   6    7    8    9    10   11    12
!            KEY = 3:  459 819 1359 2151 3315 5067 7815 12351 20235
!            KEY = 4:  195 309  483  765 1251 2133 3795  7005 13299
!     FUNSUB Externally declared subroutine for computing
!            all components of the integrand at the given
!            evaluation point.
!            It must have parameters (NDIM,X,NUMFUN,FUNVLS)
!            Input parameters:
!              NDIM   Integer that defines the dimension of the
!                     integral.
!              X      Real array of dimension NDIM
!                     that defines the evaluation point.
!              NUMFUN Integer that defines the number of
!                     components of I.
!            Output parameter:
!              FUNVLS Real array of dimension NUMFUN
!                     that defines NUMFUN components of the integrand.
!
!     EPSABS Real.
!            Requested absolute error.
!     EPSREL Real.
!            Requested relative error.
!     KEY    Integer.
!            Key to selected local integration rule.
!            KEY = 0 is the default value.
!                  For NDIM = 2 the degree 13 rule is selected.
!                  For NDIM = 3 the degree 11 rule is selected.
!                  For NDIM > 3 the degree  9 rule is selected.
!            KEY = 1 gives the user the 2 dimensional degree 13
!                  integration rule that uses 65 evaluation points.
!            KEY = 2 gives the user the 3 dimensional degree 11
!                  integration rule that uses 127 evaluation points.
!            KEY = 3 gives the user the degree 9 integration rule.
!            KEY = 4 gives the user the degree 7 integration rule.
!                  This is the recommended rule for problems that
!                  require great adaptivity.
!     NW     Integer.
!            Defines the length of the working array WORK.
!            Let MAXSUB denote the maximum allowed number of subregions
!            for the given values of MAXPTS, KEY and NDIM.
!            MAXSUB = (MAXPTS-NUM)/(2*NUM) + 1
!            NW should be greater or equal to
!            MAXSUB*(2*NDIM+2*NUMFUN+2) + 17*NUMFUN + 1
!            For efficient execution on parallel computers
!            NW should be chosen greater or equal to
!            MAXSUB*(2*NDIM+2*NUMFUN+2) + 17*NUMFUN*MDIV + 1
!            where MDIV is the number of subregions that are divided in
!            each subdivision step.
!            MDIV is default set internally in DCUHRE equal to 1.
!            For efficient execution on parallel computers
!            with NPROC processors MDIV should be set equal to
!            the smallest integer such that MOD(2*MDIV,NPROC) = 0.
!
!     RESTAR Integer.
!            If RESTAR = 0, this is the first attempt to compute
!            the integral.
!            If RESTAR = 1, then we restart a previous attempt.
!            In this case the only parameters for DCUHRE that may
!            be changed (with respect to the previous call of DCUHRE)
!            are MINPTS, MAXPTS, EPSABS, EPSREL and RESTAR.
!
!   ON RETURN
!
!     RESULT Real array of dimension NUMFUN.
!            Approximations to all components of the integral.
!     ABSERR Real array of dimension NUMFUN.
!            Estimates of absolute errors.
!     NEVAL  Integer.
!            Number of function evaluations used by DCUHRE.
!     IFAIL  Integer.
!            IFAIL = 0 for normal exit, when ABSERR(K) <=  EPSABS or
!              ABSERR(K) <=  ABS(RESULT(K))*EPSREL with MAXPTS or less
!              function evaluations for all values of K,
!              1 <= K <= NUMFUN .
!            IFAIL = 1 if MAXPTS was too small for DCUHRE
!              to obtain the required accuracy. In this case DCUHRE
!              returns values of RESULT with estimated absolute
!              errors ABSERR.
!            IFAIL = 2 if KEY is less than 0 or KEY greater than 4.
!            IFAIL = 3 if NDIM is less than 2 or NDIM greater than 15.
!            IFAIL = 4 if KEY = 1 and NDIM not equal to 2.
!            IFAIL = 5 if KEY = 2 and NDIM not equal to 3.
!            IFAIL = 6 if NUMFUN is less than 1.
!            IFAIL = 7 if volume of region of integration is zero.
!            IFAIL = 8 if MAXPTS is less than 3*NUM.
!            IFAIL = 9 if MAXPTS is less than MINPTS.
!            IFAIL = 10 if EPSABS < 0 and EPSREL < 0.
!            IFAIL = 11 if NW is too small.
!            IFAIL = 12 if unlegal RESTAR.
!     WORK   Real array of dimension NW.
!            Used as working storage.
!            WORK(NW) = NSUB, the number of subregions in the data
!            structure.
!            Let WRKSUB=(NW-1-17*NUMFUN*MDIV)/(2*NDIM+2*NUMFUN+2)
!            WORK(1),...,WORK(NUMFUN*WRKSUB) contain
!              the estimated components of the integrals over the
!              subregions.
!            WORK(NUMFUN*WRKSUB+1),...,WORK(2*NUMFUN*WRKSUB) contain
!              the estimated errors over the subregions.
!            WORK(2*NUMFUN*WRKSUB+1),...,WORK(2*NUMFUN*WRKSUB+NDIM*
!              WRKSUB) contain the centers of the subregions.
!            WORK(2*NUMFUN*WRKSUB+NDIM*WRKSUB+1),...,WORK((2*NUMFUN+
!              NDIM)*WRKSUB+NDIM*WRKSUB) contain subregion half widths.
!            WORK(2*NUMFUN*WRKSUB+2*NDIM*WRKSUB+1),...,WORK(2*NUMFUN*
!              WRKSUB+2*NDIM*WRKSUB+WRKSUB) contain the greatest errors
!              in each subregion.
!            WORK((2*NUMFUN+2*NDIM+1)*WRKSUB+1),...,WORK((2*NUMFUN+
!              2*NDIM+1)*WRKSUB+WRKSUB) contain the direction of
!              subdivision in each subregion.
!            WORK(2*(NDIM+NUMFUN+1)*WRKSUB),...,WORK(2*(NDIM+NUMFUN+1)*
!              WRKSUB+ 17*MDIV*NUMFUN) is used as temporary
!              storage in DADHRE.
!
!
!        DCUHRE Example Test Program
!
! 
!   DTEST1 is a simple test driver for DCUHRE. 
! 
!   Output produced on a SUN 3/50. 
! 
!       DCUHRE TEST RESULTS 
! 
!    FTEST CALLS = 3549, IFAIL =  0 
!   N   ESTIMATED ERROR    INTEGRAL 
!   1       0.000000       0.138508 
!   2       0.000000       0.063695 
!   3       0.000009       0.058618 
!   4       0.000000       0.054070 
!   5       0.000000       0.050056 
!   6       0.000000       0.046546 
!
!     PROGRAM DTEST1
!     EXTERNAL FTEST
!     INTEGER KEY, N, NF, NDIM, MINCLS, MAXCLS, IFAIL, NEVAL, NW
!     PARAMETER (NDIM = 5, NW = 5000, NF = NDIM+1)
!     real(dp) A(NDIM), B(NDIM), WRKSTR(NW)
!     real(dp) ABSEST(NF), FINEST(NF), ABSREQ, RELREQ
!     DO 10 N = 1,NDIM
!        A(N) = 0
!        B(N) = 1
!  10 CONTINUE
!     MINCLS = 0
!     MAXCLS = 10000
!     KEY = 0
!     ABSREQ = 0
!     RELREQ = 1E-3
!     CALL DCUHRE(NDIM, NF, A, B, MINCLS, MAXCLS, FTEST, ABSREQ, RELREQ,
!    * KEY, NW, 0, FINEST, ABSEST, NEVAL, IFAIL, WRKSTR)
!     PRINT 9999, NEVAL, IFAIL
!9999 FORMAT (8X, 'DCUHRE TEST RESULTS', //'     FTEST CALLS = ', I4,
!    * ', IFAIL = ', I2, /'    N   ESTIMATED ERROR    INTEGRAL')
!     DO 20 N = 1,NF
!        PRINT 9998, N, ABSEST(N), FINEST(N)
!9998    FORMAT (3X, I2, 2F15.6)
!  20 CONTINUE
!     END
!     SUBROUTINE FTEST(NDIM, Z, NFUN, F)
!     INTEGER N, NDIM, NFUN
!     real(dp) Z(NDIM), F(NFUN), SUM
!     SUM = 0
!     DO 10 N = 1,NDIM
!        SUM = SUM + N*Z(N)**2
!  10 CONTINUE
!     F(1) = EXP(-SUM/2)
!     DO 20 N = 1,NDIM
!        F(N+1) = Z(N)*F(1)
!  20 CONTINUE
!     END
!
!***LONG DESCRIPTION
!
!   The information for each subregion is contained in the
!   data structure WORK.
!   When passed on to DADHRE, WORK is split into eight
!   arrays VALUES, ERRORS, CENTRS, HWIDTS, GREATE, DIR,
!   OLDRES and WORK.
!   VALUES contains the estimated values of the integrals.
!   ERRORS contains the estimated errors.
!   CENTRS contains the centers of the subregions.
!   HWIDTS contains the half widths of the subregions.
!   GREATE contains the greatest estimated error for each subregion.
!   DIR    contains the directions for further subdivision.
!   OLDRES and WORK are used as work arrays in DADHRE.
!
!   The data structures for the subregions are in DADHRE organized
!   as a heap, and the size of GREATE(I) defines the position of
!   region I in the heap. The heap is maintained by the program
!   DTRHRE.
!
!   The subroutine DADHRE is written for efficient execution on shared
!   memory parallel computer. On a computer with NPROC processors we will
!   in each subdivision step divide MDIV regions, where MDIV is
!   chosen such that MOD(2*MDIV,NPROC) = 0, in totally 2*MDIV new regions.
!   Each processor will then compute estimates of the integrals and errors
!   over 2*MDIV/NPROC subregions in each subdivision step.
!   The subroutine for estimating the integral and the error over
!   each subregion, DRLHRE, uses WORK2 as a work array.
!   We must make sure that each processor writes its results to
!   separate parts of the memory, and therefore the sizes of WORK and
!   WORK2 are functions of MDIV.
!   In order to achieve parallel processing of subregions, compiler
!   directives should be placed in front of the DO 200
!   loop in DADHRE on machines like Alliant and CRAY.
!
!***REFERENCES
!   J.Berntsen, T.O.Espelid and A.Genz, An Adaptive Algorithm
!   for the Approximate Calculation of Multiple Integrals,
!   To be published.
!
!   J.Berntsen, T.O.Espelid and A.Genz, DCUHRE: An Adaptive
!   Multidimensional Integration Routine for a Vector of
!   Integrals, To be published.
!
!***ROUTINES CALLED DCHHRE,DADHRE
!***END PROLOGUE DCUHRE
!
!     The original file dcuhr!.txt follows this comment.
!
!     This is a summary of the routines in this file. 
!
!     DTEST1   Simple test driver for DCUHRE.
!              Sample output from a SUN 3/50 is included.
!     DCUHRE   Main Integrator. DCUHRE calls DCHHRE and DADHRE.
!     DCHHRE   Checks the input to DCUHRE.
!     DADHRE   The adaptive integration routine. 
!              DADHRE calls DTRHRE, DINHRE and DRLHRE.
!     DTRHRE   Maintaines the heap of subregions.
!     DINHRE   Computes weights and abscissas of the integration
!              rule. DINHRE calls D132RE, D112RE, D09HRE and D07HRE.
!     D132RE   Computes weights and abscissas for a 2-dimensional
!              rule of degree 13.
!     D113RE   Computes weights and abscissas for a 3-dimensional 
!              rule of degree 11.
!     D09HRE   Computes weights and abscissas for a degree 9 rule.
!     D07HRE   Computes weights and abscissas for a degree 7 rule.
!     DRLHRE   Computes estimates of integral and error over
!              subregions.
!     DFSHRE   Computes fully symmetric sums of function evaluations.


!   Global variables.
!
!      EXTERNAL FUNSUB
    interface
       function FUNSUB(Z)
         use nrtype
         real(dp), intent(in) :: Z(:)
         real(dp) :: FUNSUB
       end function FUNSUB
    end interface
    INTEGER NDIM,MINPTS,MAXPTS,KEY,NW,RESTAR
    INTEGER NEVAL,IFAIL
    real(dp) :: A(NDIM),B(NDIM),EPSABS,EPSREL
    integer, parameter :: NUMFUN = 1
    real(dp) :: RESULT(NUMFUN),ABSERR(NUMFUN),WORK(NW)
!
!   Local variables.
!
!   MDIV   Integer.
!          MDIV is the number of subregions that are divided in
!          each subdivision step in DADHRE.
!          MDIV is chosen default to 1.
!          For efficient execution on parallel computers
!          with NPROC processors MDIV should be set equal to
!          the smallest integer such that MOD(2*MDIV,NPROC) = 0.
!   MAXDIM Integer.
!          The maximum allowed value of NDIM.
!   MAXWT  Integer. The maximum number of weights used by the
!          integration rule.
!   WTLENG Integer.
!          The number of generators used by the selected rule.
!   WORK2  Real work space. The length
!          depends on the parameters MDIV,MAXDIM and MAXWT.
!   MAXSUB Integer.
!          The maximum allowed number of subdivisions
!          for the given values of KEY, NDIM and MAXPTS.
!   MINSUB Integer.
!          The minimum allowed number of subregions for the given
!          values of MINPTS, KEY and NDIM.
!   WRKSUB Integer.
!          The maximum allowed number of subregions as a function
!          of NW, NUMFUN, NDIM and MDIV. This determines the length
!          of the main work arrays.
!   NUM    Integer. The number of integrand evaluations needed
!          over each subregion.
!
    INTEGER MDIV,MAXWT,WTLENG,MAXDIM,LENW2,MAXSUB,MINSUB
    INTEGER NUM,NSUB,LENW,KEYF
    PARAMETER (MDIV=1)
    PARAMETER (MAXDIM=15)
    PARAMETER (MAXWT=14)
    PARAMETER (LENW2=2*MDIV*MAXDIM* (MAXWT+1)+12*MAXWT+2*MAXDIM)
    INTEGER WRKSUB,I1,I2,I3,I4,I5,I6,I7,I8,K1,K2,K3,K4,K5,K6,K7,K8
    real(dp) :: WORK2(LENW2)
!
!***FIRST EXECUTABLE STATEMENT DCUHRE
!
!   Compute NUM, WTLENG, MAXSUB and MINSUB,
!   and check the input parameters.
!
    CALL DCHHRE(MAXDIM,NDIM,NUMFUN,MDIV,A,B,MINPTS,MAXPTS,EPSABS,  &
         &            EPSREL,KEY,NW,RESTAR,NUM,MAXSUB,MINSUB,KEYF,  &
         &            IFAIL,WTLENG)
    WRKSUB = (NW - 1 - 17*MDIV*NUMFUN)/(2*NDIM + 2*NUMFUN + 2)
    IF (IFAIL.NE.0) THEN
       GO TO 999
    END IF
!
!   Split up the work space.
!
    I1 = 1
    I2 = I1 + WRKSUB*NUMFUN
    I3 = I2 + WRKSUB*NUMFUN
    I4 = I3 + WRKSUB*NDIM
    I5 = I4 + WRKSUB*NDIM
    I6 = I5 + WRKSUB
    I7 = I6 + WRKSUB
    I8 = I7 + NUMFUN*MDIV
    K1 = 1
    K2 = K1 + 2*MDIV*WTLENG*NDIM
    K3 = K2 + WTLENG*5
    K4 = K3 + WTLENG
    K5 = K4 + NDIM
    K6 = K5 + NDIM
    K7 = K6 + 2*MDIV*NDIM
    K8 = K7 + 3*WTLENG
!
!   On restart runs the number of subregions from the
!   previous call is assigned to NSUB.
!
    IF (RESTAR.EQ.1) THEN
       NSUB = WORK(NW)
    END IF
!
!   Compute the size of the temporary work space needed in DADHRE.
!
    LENW = 16*MDIV*NUMFUN
    CALL DADHRE(NDIM,NUMFUN,MDIV,A,B,MINSUB,MAXSUB,FUNSUB,EPSABS, &
         &            EPSREL,KEYF,RESTAR,NUM,LENW,WTLENG,  &
         &            RESULT,ABSERR,NEVAL,NSUB,IFAIL,WORK(I1),WORK(I2),  &
         &            WORK(I3),WORK(I4),WORK(I5),WORK(I6),WORK(I7),WORK(I8),  &
         &            WORK2(K1),WORK2(K2),WORK2(K3),WORK2(K4),WORK2(K5),  &
         &            WORK2(K6),WORK2(K7),WORK2(K8))
    WORK(NW) = NSUB
999 RETURN
!
!***END DCUHRE
!
  END SUBROUTINE DCUHRE
  SUBROUTINE DCHHRE(MAXDIM,NDIM,NUMFUN,MDIV,A,B,MINPTS,MAXPTS,  &
       &                  EPSABS,EPSREL,KEY,NW,RESTAR,NUM,MAXSUB,MINSUB,  &
       &                  KEYF,IFAIL,WTLENG)
!***BEGIN PROLOGUE DCHHRE
!***PURPOSE  DCHHRE checks the validity of the
!            input parameters to DCUHRE.
!***DESCRIPTION
!            DCHHRE computes NUM, MAXSUB, MINSUB, KEYF, WTLENG and 
!            IFAIL as functions of the input parameters to DCUHRE,
!            and checks the validity of the input parameters to DCUHRE.
!
!   ON ENTRY
!
!     MAXDIM Integer.
!            The maximum allowed number of dimensions.
!     NDIM   Integer.
!            Number of variables. 1 < NDIM <= MAXDIM.
!     NUMFUN Integer.
!            Number of components of the integral.
!     MDIV   Integer.
!            MDIV is the number of subregions that are divided in
!            each subdivision step in DADHRE.
!            MDIV is chosen default to 1.
!            For efficient execution on parallel computers
!            with NPROC processors MDIV should be set equal to
!            the smallest integer such that MOD(2*MDIV,NPROC) = 0.
!     A      Real array of dimension NDIM.
!            Lower limits of integration.
!     B      Real array of dimension NDIM.
!            Upper limits of integration.
!     MINPTS Integer.
!            Minimum number of function evaluations.
!     MAXPTS Integer.
!            Maximum number of function evaluations.
!            The number of function evaluations over each subregion
!            is NUM.
!            If (KEY = 0 or KEY = 1) and (NDIM = 2) Then
!              NUM = 65
!            Elseif (KEY = 0 or KEY = 2) and (NDIM = 3) Then
!              NUM = 127
!            Elseif (KEY = 0 and NDIM > 3) or (KEY = 3) Then
!              NUM = 1 + 4*2*NDIM + 2*NDIM*(NDIM-1) + 4*NDIM*(NDIM-1) +
!                    4*NDIM*(NDIM-1)*(NDIM-2)/3 + 2**NDIM
!            Elseif (KEY = 4) Then
!              NUM = 1 + 3*2*NDIM + 2*NDIM*(NDIM-1) + 2**NDIM
!            MAXPTS >= 3*NUM and MAXPTS >= MINPTS
!     EPSABS Real.
!            Requested absolute error.
!     EPSREL Real.
!            Requested relative error.
!     KEY    Integer.
!            Key to selected local integration rule.
!            KEY = 0 is the default value.
!                  For NDIM = 2 the degree 13 rule is selected.
!                  For NDIM = 3 the degree 11 rule is selected.
!                  For NDIM > 3 the degree  9 rule is selected.
!            KEY = 1 gives the user the 2 dimensional degree 13
!                  integration rule that uses 65 evaluation points.
!            KEY = 2 gives the user the 3 dimensional degree 11
!                  integration rule that uses 127 evaluation points.
!            KEY = 3 gives the user the degree 9 integration rule.
!            KEY = 4 gives the user the degree 7 integration rule.
!                  This is the recommended rule for problems that
!                  require great adaptivity.
!     NW     Integer.
!            Defines the length of the working array WORK.
!            Let MAXSUB denote the maximum allowed number of subregions
!            for the given values of MAXPTS, KEY and NDIM.
!            MAXSUB = (MAXPTS-NUM)/(2*NUM) + 1
!            NW should be greater or equal to
!            MAXSUB*(2*NDIM+2*NUMFUN+2) + 17*NUMFUN + 1
!            For efficient execution on parallel computers
!            NW should be chosen greater or equal to
!            MAXSUB*(2*NDIM+2*NUMFUN+2) + 17*NUMFUN*MDIV + 1
!            where MDIV is the number of subregions that are divided in
!            each subdivision step.
!            MDIV is default set internally in DCUHRE equal to 1.
!            For efficient execution on parallel computers
!            with NPROC processors MDIV should be set equal to
!            the smallest integer such that MOD(2*MDIV,NPROC) = 0.
!     RESTAR Integer.
!            If RESTAR = 0, this is the first attempt to compute
!            the integral.
!            If RESTAR = 1, then we restart a previous attempt.
!
!   ON RETURN
!
!     NUM    Integer.
!            The number of function evaluations over each subregion.
!     MAXSUB Integer.
!            The maximum allowed number of subregions for the
!            given values of MAXPTS, KEY and NDIM.
!     MINSUB Integer.
!            The minimum allowed number of subregions for the given
!            values of MINPTS, KEY and NDIM.
!     IFAIL  Integer.
!            IFAIL = 0 for normal exit.
!            IFAIL = 2 if KEY is less than 0 or KEY greater than 4.
!            IFAIL = 3 if NDIM is less than 2 or NDIM greater than
!                      MAXDIM.
!            IFAIL = 4 if KEY = 1 and NDIM not equal to 2.
!            IFAIL = 5 if KEY = 2 and NDIM not equal to 3.
!            IFAIL = 6 if NUMFUN less than 1.
!            IFAIL = 7 if volume of region of integration is zero.
!            IFAIL = 8 if MAXPTS is less than 3*NUM.
!            IFAIL = 9 if MAXPTS is less than MINPTS.
!            IFAIL = 10 if EPSABS < 0 and EPSREL < 0.
!            IFAIL = 11 if NW is too small.
!            IFAIL = 12 if unlegal RESTAR.
!     KEYF   Integer.
!            Key to selected integration rule.
!     WTLENG Integer.
!            The number of generators of the chosen integration rule.
!
!***ROUTINES CALLED-NONE
!***END PROLOGUE DCHHRE
!
!   Global variables.
!
    INTEGER NDIM,NUMFUN,MDIV,MINPTS,MAXPTS,KEY,NW,MINSUB,MAXSUB
    INTEGER RESTAR,NUM,KEYF,IFAIL,MAXDIM,WTLENG
    real(dp) :: A(NDIM),B(NDIM),EPSABS,EPSREL
!
!   Local variables.
!
    INTEGER LIMIT,J
!
!***FIRST EXECUTABLE STATEMENT DCHHRE
!
    IFAIL = 0
!
!   Check on legal KEY.
!
    IF (KEY.LT.0 .OR. KEY.GT.4) THEN
       IFAIL = 2
       GO TO 999
    END IF
!
!   Check on legal NDIM.
!
    IF (NDIM.LT.2 .OR. NDIM.GT.MAXDIM) THEN
       IFAIL = 3
       GO TO 999
    END IF
!
!   For KEY = 1, NDIM must be equal to 2.
!
    IF (KEY.EQ.1 .AND. NDIM.NE.2) THEN
       IFAIL = 4
       GO TO 999
    END IF
!
!   For KEY = 2, NDIM must be equal to 3.
!
    IF (KEY.EQ.2 .AND. NDIM.NE.3) THEN
       IFAIL = 5
       GO TO 999
    END IF
!
!   For KEY = 0, we point at the selected integration rule.
!
    IF (KEY.EQ.0) THEN
       IF (NDIM.EQ.2) THEN
          KEYF = 1
       ELSE IF (NDIM.EQ.3) THEN
          KEYF = 2
       ELSE
          KEYF = 3
       ENDIF
    ELSE
       KEYF = KEY
    ENDIF
!
!   Compute NUM and WTLENG as a function of KEYF and NDIM.
!
    IF (KEYF.EQ.1) THEN
       NUM = 65
       WTLENG = 14
    ELSE IF (KEYF.EQ.2) THEN
       NUM = 127
       WTLENG = 13
    ELSE IF (KEYF.EQ.3) THEN
       NUM = 1 + 4*2*NDIM + 2*NDIM* (NDIM-1) + 4*NDIM* (NDIM-1) +  &
            &          4*NDIM* (NDIM-1)* (NDIM-2)/3 + 2**NDIM
       WTLENG = 9
       IF (NDIM.EQ.2) WTLENG = 8
    ELSE IF (KEYF.EQ.4) THEN
       NUM = 1 + 3*2*NDIM + 2*NDIM* (NDIM-1) + 2**NDIM
       WTLENG = 6
    END IF
!
!   Compute MAXSUB.
!
    MAXSUB = (MAXPTS-NUM)/ (2*NUM) + 1
!
!   Compute MINSUB.
!
    MINSUB = (MINPTS-NUM)/ (2*NUM) + 1
    IF (MOD(MINPTS-NUM,2*NUM).NE.0) THEN
       MINSUB = MINSUB + 1
    END IF
    MINSUB = MAX(2,MINSUB)
!
!   Check on positive NUMFUN.
!
    IF (NUMFUN.LT.1) THEN
       IFAIL = 6
       GO TO 999
    END IF
!
!   Check on legal upper and lower limits of integration.
!
    DO 10 J = 1,NDIM
       IF (A(J)-B(J).EQ.0) THEN
          IFAIL = 7
          GO TO 999
       END IF
10  CONTINUE
!
!   Check on MAXPTS < 3*NUM.
!
    IF (MAXPTS.LT.3*NUM) THEN
       IFAIL = 8
       GO TO 999
    END IF
!
!   Check on MAXPTS >= MINPTS.
!
    IF (MAXPTS.LT.MINPTS) THEN
       IFAIL = 9
       GO TO 999
    END IF
!
!   Check on legal accuracy requests.
!
    IF (EPSABS.LT.0 .AND. EPSREL.LT.0) THEN
       IFAIL = 10
       GO TO 999
    END IF
!
!   Check on big enough double precision workspace.
!
    LIMIT = MAXSUB* (2*NDIM+2*NUMFUN+2) + 17*MDIV*NUMFUN + 1
    IF (NW.LT.LIMIT) THEN
       IFAIL = 11
       GO TO 999
    END IF
!
!    Check on legal RESTAR.
!
    IF (RESTAR.NE.0 .AND. RESTAR.NE.1) THEN
       IFAIL = 12
       GO TO 999
    END IF
999 RETURN
!
!***END DCHHRE
!
  END subroutine DCHHRE
  SUBROUTINE DADHRE(NDIM,NUMFUN,MDIV,A,B,MINSUB,MAXSUB,FUNSUB,  &
       &                  EPSABS,EPSREL,KEY,RESTAR,NUM,LENW,WTLENG,  &
       &                  RESULT,ABSERR,NEVAL,NSUB,IFAIL,VALUES,  &
       &                  ERRORS,CENTRS,HWIDTS,GREATE,DIR,OLDRES,WORK,G,W,  &
       &                  RULPTS,CENTER,HWIDTH,X,SCALES,NORMS)
!***BEGIN PROLOGUE DADHRE
!***KEYWORDS automatic multidimensional integrator,
!            n-dimensional hyper-rectangles,
!            general purpose, global adaptive
!***PURPOSE  The routine calculates an approximation to a given
!            vector of definite integrals, I, over a hyper-rectangular
!            region hopefully satisfying for each component of I the
!            following claim for accuracy:
!            ABS(I(K)-RESULT(K)).LE.MAX(EPSABS,EPSREL*ABS(I(K)))
!***DESCRIPTION Computation of integrals over hyper-rectangular
!            regions.
!            DADHRE repeatedly subdivides the region
!            of integration and estimates the integrals and the
!            errors over the subregions with  greatest
!            estimated errors until the error request
!            is met or MAXSUB subregions are stored.
!            The regions are divided in two equally sized parts along
!            the direction with greatest absolute fourth divided
!            difference.
!
!   ON ENTRY
!
!     NDIM   Integer.
!            Number of variables. 1 < NDIM <= MAXDIM.
!     NUMFUN Integer.
!            Number of components of the integral.
!     MDIV   Integer.
!            Defines the number of new subregions that are divided
!            in each subdivision step.
!     A      Real array of dimension NDIM.
!            Lower limits of integration.
!     B      Real array of dimension NDIM.
!            Upper limits of integration.
!     MINSUB Integer.
!            The computations proceed until there are at least
!            MINSUB subregions in the data structure.
!     MAXSUB Integer.
!            The computations proceed until there are at most
!            MAXSUB subregions in the data structure.
!
!     FUNSUB Externally declared subroutine for computing
!            all components of the integrand in the given
!            evaluation point.
!            It must have parameters (NDIM,X,NUMFUN,FUNVLS)
!            Input parameters:
!              NDIM   Integer that defines the dimension of the
!                     integral.
!              X      Real array of dimension NDIM
!                     that defines the evaluation point.
!              NUMFUN Integer that defines the number of
!                     components of I.
!            Output parameter:
!              FUNVLS Real array of dimension NUMFUN
!                     that defines NUMFUN components of the integrand.
!
!     EPSABS Real.
!            Requested absolute error.
!     EPSREL Real.
!            Requested relative error.
!     KEY    Integer.
!            Key to selected local integration rule.
!            KEY = 0 is the default value.
!                  For NDIM = 2 the degree 13 rule is selected.
!                  For NDIM = 3 the degree 11 rule is selected.
!                  For NDIM > 3 the degree  9 rule is selected.
!            KEY = 1 gives the user the 2 dimensional degree 13
!                  integration rule that uses 65 evaluation points.
!            KEY = 2 gives the user the 3 dimensional degree 11
!                  integration rule that uses 127 evaluation points.
!            KEY = 3 gives the user the degree 9 integration rule.
!            KEY = 4 gives the user the degree 7 integration rule.
!                  This is the recommended rule for problems that
!                  require great adaptivity.
!     RESTAR Integer.
!            If RESTAR = 0, this is the first attempt to compute
!            the integral.
!            If RESTAR = 1, then we restart a previous attempt.
!            (In this case the output parameters from DADHRE
!            must not be changed since the last
!            exit from DADHRE.)
!     NUM    Integer.
!            The number of function evaluations over each subregion.
!     LENW   Integer.
!            Defines the length of the working array WORK.
!            LENW should be greater or equal to
!            16*MDIV*NUMFUN.
!     WTLENG Integer.
!            The number of weights in the basic integration rule.
!     NSUB   Integer.
!            If RESTAR = 1, then NSUB must specify the number
!            of subregions stored in the previous call to DADHRE.
!
!   ON RETURN
!
!     RESULT Real array of dimension NUMFUN.
!            Approximations to all components of the integral.
!     ABSERR Real array of dimension NUMFUN.
!            Estimates of absolute accuracies.
!     NEVAL  Integer.
!            Number of function evaluations used by DADHRE.
!     NSUB   Integer.
!            Number of stored subregions.
!     IFAIL  Integer.
!            IFAIL = 0 for normal exit, when ABSERR(K) <=  EPSABS or
!              ABSERR(K) <=  ABS(RESULT(K))*EPSREL with MAXSUB or less
!              subregions processed for all values of K,
!              1 <=  K <=  NUMFUN.
!            IFAIL = 1 if MAXSUB was too small for DADHRE
!              to obtain the required accuracy. In this case DADHRE
!              returns values of RESULT with estimated absolute
!              accuracies ABSERR.
!     VALUES Real array of dimension (NUMFUN,MAXSUB).
!            Used to store estimated values of the integrals
!            over the subregions.
!     ERRORS Real array of dimension (NUMFUN,MAXSUB).
!            Used to store the corresponding estimated errors.
!     CENTRS Real array of dimension (NDIM,MAXSUB).
!            Used to store the centers of the stored subregions.
!     HWIDTS Real array of dimension (NDIM,MAXSUB).
!            Used to store the half widths of the stored subregions.
!     GREATE Real array of dimension MAXSUB.
!            Used to store the greatest estimated errors in
!            all subregions.
!     DIR    Real array of dimension MAXSUB.
!            DIR is used to store the directions for
!            further subdivision.
!     OLDRES Real array of dimension (NUMFUN,MDIV).
!            Used to store old estimates of the integrals over the
!            subregions.
!     WORK   Real array of dimension LENW.
!            Used  in DRLHRE and DTRHRE.
!     G      Real array of dimension (NDIM,WTLENG,2*MDIV).
!            The fully symmetric sum generators for the rules.
!            G(1,J,1),...,G(NDIM,J,1) are the generators for the
!            points associated with the Jth weights.
!            When MDIV subregions are divided in 2*MDIV
!            subregions, the subregions may be processed on different
!            processors and we must make a copy of the generators
!            for each processor.
!     W      Real array of dimension (5,WTLENG).
!            The weights for the basic and null rules.
!            W(1,1), ..., W(1,WTLENG) are weights for the basic rule.
!            W(I,1), ..., W(I,WTLENG) , for I > 1 are null rule weights.
!     RULPTS Real array of dimension WTLENG.
!            Work array used in DINHRE.
!     CENTER Real array of dimension NDIM.
!            Work array used in DTRHRE.
!     HWIDTH Real array of dimension NDIM.
!            Work array used in DTRHRE.
!     X      Real array of dimension (NDIM,2*MDIV).
!            Work array used in DRLHRE.
!     SCALES Real array of dimension (3,WTLENG).
!            Work array used by DINHRE and DRLHRE.
!     NORMS  Real array of dimension (3,WTLENG).
!            Work array used by DINHRE and DRLHRE.
!
!***REFERENCES
!
!   P. van Dooren and L. de Ridder, Algorithm 6, An adaptive algorithm
!   for numerical integration over an n-dimensional cube, J.Comput.Appl.
!   Math. 2(1976)207-217.
!
!   A.C.Genz and A.A.Malik, Algorithm 019. Remarks on algorithm 006:
!   An adaptive algorithm for numerical integration over an
!   N-dimensional rectangular region,J.Comput.Appl.Math. 6(1980)295-302.
!
!***ROUTINES CALLED DTRHRE,DINHRE,DRLHRE
!***END PROLOGUE DADHRE
!
!   Global variables.
! 
!      EXTERNAL FUNSUB
    interface
       function FUNSUB(Z)
         use nrtype
         real(dp), intent(in):: Z(:)
         real(dp) :: FUNSUB
       end function FUNSUB
    end interface
    INTEGER NDIM,NUMFUN,MDIV,MINSUB,MAXSUB,KEY,LENW,RESTAR
    INTEGER NUM,NEVAL,NSUB,IFAIL,WTLENG
    real(dp) :: A(NDIM),B(NDIM),EPSABS,EPSREL
    real(dp) :: RESULT(NUMFUN),ABSERR(NUMFUN)
    real(dp) :: VALUES(NUMFUN,MAXSUB),ERRORS(NUMFUN,MAXSUB)
    real(dp) :: CENTRS(NDIM,MAXSUB)
    real(dp) :: HWIDTS(NDIM,MAXSUB)
    real(dp) :: GREATE(MAXSUB),DIR(MAXSUB)
    real(dp) :: OLDRES(NUMFUN,MDIV)
    real(dp) :: WORK(LENW),RULPTS(WTLENG)
    real(dp) :: G(NDIM,WTLENG,2*MDIV),W(5,WTLENG)
    real(dp) :: CENTER(NDIM),HWIDTH(NDIM),X(NDIM,2*MDIV)
    real(dp) :: SCALES(3,WTLENG),NORMS(3,WTLENG)
!
!   Local variables.
!
!   INTSGN is used to get correct sign on the integral.
!   SBRGNS is the number of stored subregions.
!   NDIV   The number of subregions to be divided in each main step.
!   POINTR Pointer to the position in the data structure where
!          the new subregions are to be stored.
!   DIRECT Direction of subdivision.
!   ERRCOF Heuristic error coeff. defined in DINHRE and used by DRLHRE
!          and DADHRE.
!
    INTEGER I,J,K
    INTEGER INTSGN,SBRGNS
    INTEGER L1
    INTEGER NDIV,POINTR,DIRECT,INDEX
    real(dp) :: OLDCEN,EST1,EST2,ERRCOF(6)
!
!***FIRST EXECUTABLE STATEMENT DADHRE
!
!   Get the correct sign on the integral.
!
    INTSGN = 1
    DO 10 J = 1,NDIM
       IF (B(J).LT.A(J)) THEN
          INTSGN = - INTSGN
       END IF
10  CONTINUE
!
!   Call DINHRE to compute the weights and abscissas of
!   the function evaluation points.
!
    CALL DINHRE(NDIM,KEY,WTLENG,W,G,ERRCOF,RULPTS,SCALES,NORMS)
!
!   If RESTAR = 1, then this is a restart run.
!
    IF (RESTAR.EQ.1) THEN
       SBRGNS = NSUB
       GO TO 110
    END IF
!
!   Initialize the SBRGNS, CENTRS and HWIDTS.
!
    SBRGNS = 1
    DO 15 J = 1,NDIM
       CENTRS(J,1) = (A(J)+B(J))/2
       HWIDTS(J,1) = ABS(B(J)-A(J))/2
15  CONTINUE
!
!   Initialize RESULT, ABSERR and NEVAL.
!
    DO 20 J = 1,NUMFUN
       RESULT(J) = 0
       ABSERR(J) = 0
20  CONTINUE
    NEVAL = 0
!
!   Apply DRLHRE over the whole region.
!
    CALL DRLHRE(NDIM,CENTRS(1,1),HWIDTS(1,1),WTLENG,G,W,ERRCOF,NUMFUN,  &
         &            FUNSUB,SCALES,NORMS,X,WORK,VALUES(1,1),ERRORS(1,1),  &
         &            DIR(1))
    NEVAL = NEVAL + NUM
!
!   Add the computed values to RESULT and ABSERR.
!
    DO 55 J = 1,NUMFUN
       RESULT(J) = RESULT(J) + VALUES(J,1)
55  CONTINUE
    DO 65 J = 1,NUMFUN
       ABSERR(J) = ABSERR(J) + ERRORS(J,1)
65  CONTINUE
!
!   Store results in heap.
!
    INDEX = 1
    CALL DTRHRE(2,NDIM,NUMFUN,INDEX,VALUES,ERRORS,CENTRS,HWIDTS,  &
         &            GREATE,WORK(1),WORK(NUMFUN+1),CENTER,HWIDTH,DIR)
!
!***End initialisation.
!
!***Begin loop while the error is too great,
!   and SBRGNS+1 is less than MAXSUB.
!
110 IF (SBRGNS+1.LE.MAXSUB) THEN
!
!   If we are allowed to divide further,
!   prepare to apply basic rule over each half of the
!   NDIV subregions with greatest errors.
!   If MAXSUB is great enough, NDIV = MDIV
!
       IF (MDIV.GT.1) THEN
          NDIV = MAXSUB - SBRGNS
          NDIV = MIN(NDIV,MDIV,SBRGNS)
       ELSE
          NDIV = 1
       END IF
!
!   Divide the NDIV subregions in two halves, and compute
!   integral and error over each half.
!
       DO 150 I = 1,NDIV
          POINTR = SBRGNS + NDIV + 1 - I
!
!   Adjust RESULT and ABSERR.
!
          DO 115 J = 1,NUMFUN
             RESULT(J) = RESULT(J) - VALUES(J,1)
             ABSERR(J) = ABSERR(J) - ERRORS(J,1)
115       CONTINUE
!
!   Compute first half region.
!
          DO 120 J = 1,NDIM
             CENTRS(J,POINTR) = CENTRS(J,1)
             HWIDTS(J,POINTR) = HWIDTS(J,1)
120       CONTINUE
          DIRECT = DIR(1)
          DIR(POINTR) = DIRECT
          HWIDTS(DIRECT,POINTR) = HWIDTS(DIRECT,1)/2
          OLDCEN = CENTRS(DIRECT,1)
          CENTRS(DIRECT,POINTR) = OLDCEN - HWIDTS(DIRECT,POINTR)
!
!   Save the computed values of the integrals.
!
          DO 125 J = 1,NUMFUN
             OLDRES(J,NDIV-I+1) = VALUES(J,1)
125       CONTINUE
!
!   Adjust the heap.
!
          CALL DTRHRE(1,NDIM,NUMFUN,SBRGNS,VALUES,ERRORS,CENTRS,  &
               &      HWIDTS,GREATE,WORK(1),WORK(NUMFUN+1),CENTER,  &
               &      HWIDTH,DIR)
!
!   Compute second half region.
!
          DO 130 J = 1,NDIM
             CENTRS(J,POINTR-1) = CENTRS(J,POINTR)
             HWIDTS(J,POINTR-1) = HWIDTS(J,POINTR)
130       CONTINUE
          CENTRS(DIRECT,POINTR-1) = OLDCEN + HWIDTS(DIRECT,POINTR)
          HWIDTS(DIRECT,POINTR-1) = HWIDTS(DIRECT,POINTR)
          DIR(POINTR-1) = DIRECT
150    CONTINUE
!
!   Make copies of the generators for each processor.
!
       DO  I = 2,2*NDIV
           DO  J = 1,NDIM
               DO  K = 1,WTLENG
                  G(J,K,I) = G(J,K,1)
               end DO
            end DO
         end DO
!
!   Apply basic rule.
!
!vd$l cncall
       DO 200 I = 1,2*NDIV
          INDEX = SBRGNS + I
          L1 = 1 + (I-1)*8*NUMFUN
          CALL DRLHRE(NDIM,CENTRS(1,INDEX),HWIDTS(1,INDEX),WTLENG,  &
               &       G(1,1,I),W,ERRCOF,NUMFUN,FUNSUB,SCALES,NORMS,  &
               &       X(1,I),WORK(L1),VALUES(1,INDEX),  &
               &       ERRORS(1,INDEX),DIR(INDEX))
200    CONTINUE
       NEVAL = NEVAL + 2*NDIV*NUM
!
!   Add new contributions to RESULT.
!
       DO 220 I = 1,2*NDIV
          DO 210 J = 1,NUMFUN
             RESULT(J) = RESULT(J) + VALUES(J,SBRGNS+I)
210       CONTINUE
220    CONTINUE
!
!   Check consistency of results and if necessary adjust
!   the estimated errors.
!
       DO 240 I = 1,NDIV
          GREATE(SBRGNS+2*I-1) = 0
          GREATE(SBRGNS+2*I) = 0
          DO 230 J = 1,NUMFUN
             EST1 = ABS(OLDRES(J,I)- (VALUES(J,  &
                  &            SBRGNS+2*I-1)+VALUES(J,SBRGNS+2*I)))
             EST2 = ERRORS(J,SBRGNS+2*I-1) + ERRORS(J,SBRGNS+2*I)
             IF (EST2.GT.0) THEN
                ERRORS(J,SBRGNS+2*I-1) = ERRORS(J,SBRGNS+2*I-1)*  &
                     &                  (1+ERRCOF(5)*EST1/EST2)
                ERRORS(J,SBRGNS+2*I) = ERRORS(J,SBRGNS+2*I)*  &
                     &                       (1+ERRCOF(5)*EST1/EST2)
             END IF
             ERRORS(J,SBRGNS+2*I-1) = ERRORS(J,SBRGNS+2*I-1) +  &
                  &                           ERRCOF(6)*EST1
             ERRORS(J,SBRGNS+2*I) = ERRORS(J,SBRGNS+2*I) +  &
                  &                            ERRCOF(6)*EST1
             IF (ERRORS(J,SBRGNS+2*I-1).GT.  &
                  &      GREATE(SBRGNS+2*I-1)) THEN
                GREATE(SBRGNS+2*I-1) = ERRORS(J,SBRGNS+2*I-1)
             END IF
             IF (ERRORS(J,SBRGNS+2*I).GT.GREATE(SBRGNS+2*I)) THEN
                GREATE(SBRGNS+2*I) = ERRORS(J,SBRGNS+2*I)
             END IF
             ABSERR(J) = ABSERR(J) + ERRORS(J,SBRGNS+2*I-1) +  &
                  &       ERRORS(J,SBRGNS+2*I)
230       CONTINUE
240    CONTINUE
!
!   Store results in heap.
!
       DO 250 I = 1,2*NDIV
          INDEX = SBRGNS + I
          CALL DTRHRE(2,NDIM,NUMFUN,INDEX,VALUES,ERRORS,CENTRS,  &
               &      HWIDTS,GREATE,WORK(1),WORK(NUMFUN+1),CENTER,  &
               &      HWIDTH,DIR)
250    CONTINUE
       SBRGNS = SBRGNS + 2*NDIV
!
!   Check for termination.
!
       IF (SBRGNS.LT.MINSUB) THEN
          GO TO 110
       END IF
       DO 255 J = 1,NUMFUN
          IF (ABSERR(J).GT.EPSREL*ABS(RESULT(J)) .AND.  &
               &       ABSERR(J).GT.EPSABS) THEN
             GO TO 110
          END IF
255    CONTINUE
       IFAIL = 0
       GO TO 499
!
!   Else we did not succeed with the
!   given value of MAXSUB.
!
    ELSE
       IFAIL = 1
    END IF
!
!   Compute more accurate values of RESULT and ABSERR.
!
499 CONTINUE
    DO 500 J = 1,NUMFUN
       RESULT(J) = 0
       ABSERR(J) = 0
500 CONTINUE
    DO 510 I = 1,SBRGNS
       DO 505 J = 1,NUMFUN
          RESULT(J) = RESULT(J) + VALUES(J,I)
          ABSERR(J) = ABSERR(J) + ERRORS(J,I)
505    CONTINUE
510 CONTINUE
!
!   Compute correct sign on the integral.
!
    DO 600 J = 1,NUMFUN
       RESULT(J) = RESULT(J)*INTSGN
600 CONTINUE
    NSUB = SBRGNS
    RETURN
!
!***END DADHRE
!
  END subroutine DADHRE
  SUBROUTINE DINHRE(NDIM,KEY,WTLENG,W,G,ERRCOF,RULPTS,SCALES,NORMS)
!***BEGIN PROLOGUE DINHRE
!***PURPOSE DINHRE computes abscissas and weights of the integration
!            rule and the null rules to be used in error estimation.
!            These are computed as functions of NDIM and KEY.
!***DESCRIPTION DINHRE will for given values of NDIM and KEY compute or
!            select the correct values of the abscissas and
!            corresponding weights for different
!            integration rules and null rules and assign them to
!            G and W.
!            The heuristic error coefficients ERRCOF
!            will be computed as a function of KEY.
!            Scaling factors SCALES and normalization factors NORMS
!            used in the error estimation are computed.
!
!
!   ON ENTRY
!
!     NDIM   Integer.
!            Number of variables.
!     KEY    Integer.
!            Key to selected local integration rule.
!     WTLENG Integer.
!            The number of weights in each of the rules.
!
!   ON RETURN
!
!     W      Real array of dimension (5,WTLENG).
!            The weights for the basic and null rules.
!            W(1,1), ...,W(1,WTLENG) are weights for the basic rule.
!            W(I,1), ...,W(I,WTLENG), for I > 1 are null rule weights.
!     G      Real array of dimension (NDIM,WTLENG).
!            The fully symmetric sum generators for the rules.
!            G(1,J),...,G(NDIM,J) are the generators for the points
!            associated with the the Jth weights.
!     ERRCOF Real array of dimension 6.
!            Heuristic error coefficients that are used in the
!            error estimation in BASRUL.
!            It is assumed that the error is computed using:
!             IF (N1*ERRCOF(1) < N2 and N2*ERRCOF(2) < N3)
!               THEN ERROR = ERRCOF(3)*N1
!               ELSE ERROR = ERRCOF(4)*MAX(N1, N2, N3)
!             ERROR = ERROR + EP*(ERRCOF(5)*ERROR/(ES+ERROR)+ERRCOF(6))
!            where N1-N3 are the null rules, EP is the error for
!            the parent
!            subregion and ES is the error for the sibling subregion.
!     RULPTS Real array of dimension WTLENG.
!            A work array containing the number of points produced by
!            each generator of the selected rule.
!     SCALES Real array of dimension (3,WTLENG).
!            Scaling factors used to construct new null rules,
!            N1, N2 and N3,
!            based on a linear combination of two successive null rules
!            in the sequence of null rules.
!     NORMS  Real array of dimension (3,WTLENG).
!            2**NDIM/(1-norm of the null rule constructed by each of
!            the scaling factors.)
!
!***ROUTINES CALLED  D132RE,D113RE,D07HRE,D09HRE
!***END PROLOGUE DINHRE
!
!   Global variables.
!
    INTEGER NDIM,KEY,WTLENG
    real(dp) :: G(NDIM,WTLENG),W(5,WTLENG),ERRCOF(6)
    real(dp) :: RULPTS(WTLENG),SCALES(3,WTLENG)
    real(dp) :: NORMS(3,WTLENG)
!
!   Local variables.
!
    INTEGER I,J,K
    real(dp) :: WE(14)
!
!***FIRST EXECUTABLE STATEMENT DINHRE
!
!   Compute W, G and ERRCOF.
!
    IF (KEY.EQ.1) THEN
       CALL D132RE(WTLENG,W,G,ERRCOF,RULPTS)
    ELSE IF (KEY.EQ.2) THEN
       CALL D113RE(WTLENG,W,G,ERRCOF,RULPTS)
    ELSE IF (KEY.EQ.3) THEN
       CALL D09HRE(NDIM,WTLENG,W,G,ERRCOF,RULPTS)
    ELSE IF (KEY.EQ.4) THEN
       CALL D07HRE(NDIM,WTLENG,W,G,ERRCOF,RULPTS)
    END IF
!
!   Compute SCALES and NORMS.
!
    DO 100 K = 1,3
       DO 50 I = 1,WTLENG
          IF (W(K+1,I).NE.0) THEN
             SCALES(K,I) = - W(K+2,I)/W(K+1,I)
          ELSE
             SCALES(K,I) = 100
          END IF
          DO 30 J = 1,WTLENG
             WE(J) = W(K+2,J) + SCALES(K,I)*W(K+1,J)
30        CONTINUE
          NORMS(K,I) = 0
          DO 40 J = 1,WTLENG
             NORMS(K,I) = NORMS(K,I) + RULPTS(J)*ABS(WE(J))
40        CONTINUE
             NORMS(K,I) = 2**NDIM/NORMS(K,I)
50     CONTINUE
100 CONTINUE
    RETURN
!
!***END DINHRE
!
  END subroutine DINHRE
  SUBROUTINE D132RE(WTLENG,W,G,ERRCOF,RULPTS)
!***BEGIN PROLOGUE D132RE
!***AUTHOR   Jarle Berntsen, EDB-senteret,
!            University of Bergen, Thormohlens gt. 55,
!            N-5008 Bergen, NORWAY
!***PURPOSE D132RE computes abscissas and weights of a 2 dimensional
!            integration rule of degree 13.
!            Two null rules of degree 11, one null rule of degree 9
!            and one null rule of degree 7 to be used in error
!            estimation are also computed.
! ***DESCRIPTION D132RE will select the correct values of the abscissas
!            and corresponding weights for different
!            integration rules and null rules and assign them to
!            G and W. The heuristic error coefficients ERRCOF
!            will also be assigned.
!
!
!   ON ENTRY
!
!     WTLENG Integer.
!            The number of weights in each of the rules.
!
!   ON RETURN
!
!     W      Real array of dimension (5,WTLENG).
!            The weights for the basic and null rules.
!            W(1,1),...,W(1,WTLENG) are weights for the basic rule.
!            W(I,1),...,W(I,WTLENG), for I > 1 are null rule weights.
!     G      Real array of dimension (NDIM,WTLENG).
!            The fully symmetric sum generators for the rules.
!            G(1,J),...,G(NDIM,J) are the generators for the points
!            associated with the the Jth weights.
!     ERRCOF Real array of dimension 6.
!            Heuristic error coefficients that are used in the
!            error estimation in BASRUL.
!     RULPTS Real array of dimension WTLENG.
!            The number of points produced by each generator.
!***REFERENCES S.Eriksen,
!              Thesis of the degree cand.scient, Dept. of Informatics,
!              Univ. of Bergen,Norway, 1984.
!
!***ROUTINES CALLED-NONE
!***END PROLOGUE D132RE
!
!   Global variables
!
    INTEGER WTLENG
    real(dp) :: W(5,WTLENG),G(2,WTLENG),ERRCOF(6)
    real(dp) :: RULPTS(WTLENG)
!
!   Local variables.
!
    INTEGER I,J
    real(dp) :: DIM2G(16)
    real(dp) :: DIM2W(14,5)
!
    DATA (DIM2G(I),I=1,16)/0.2517129343453109D+00,  &
         &     0.7013933644534266D+00,0.9590960631619962D+00,  &
         &     0.9956010478552127D+00,0.5000000000000000D+00,  &
         &     0.1594544658297559D+00,0.3808991135940188D+00,  &
         &     0.6582769255267192D+00,0.8761473165029315D+00,  &
         &     0.9982431840531980D+00,0.9790222658168462D+00,  &
         &     0.6492284325645389D+00,0.8727421201131239D+00,  &
         &     0.3582614645881228D+00,0.5666666666666666D+00,  &
         &     0.2077777777777778D+00/
!
    DATA (DIM2W(I,1),I=1,14)/0.3379692360134460D-01,  &
         &     0.9508589607597761D-01,0.1176006468056962D+00,  &
         &     0.2657774586326950D-01,0.1701441770200640D-01,  &
         &     0.0000000000000000D+00,0.1626593098637410D-01,  &
         &     0.1344892658526199D+00,0.1328032165460149D+00,  &
         &     0.5637474769991870D-01,0.3908279081310500D-02,  &
         &     0.3012798777432150D-01,0.1030873234689166D+00,  &
         &     0.6250000000000000D-01/
!
    DATA (DIM2W(I,2),I=1,14)/0.3213775489050763D+00,  &
         &     - .1767341636743844D+00,0.7347600537466072D-01,  &
         &     - .3638022004364754D-01,0.2125297922098712D-01,  &
         &     0.1460984204026913D+00,0.1747613286152099D-01,  &
         &     0.1444954045641582D+00,0.1307687976001325D-03,  &
         &     0.5380992313941161D-03,0.1042259576889814D-03,  &
         &     - .1401152865045733D-02,0.8041788181514763D-02,  &
         &     - .1420416552759383D+00/
!
    DATA (DIM2W(I,3),I=1,14)/0.3372900883288987D+00,  &
         &     - .1644903060344491D+00,0.7707849911634622D-01,  &
         &     - .3804478358506310D-01,0.2223559940380806D-01,  &
         &     0.1480693879765931D+00,0.4467143702185814D-05,  &
         &     0.1508944767074130D+00,0.3647200107516215D-04,  &
         &     0.5777198999013880D-03,0.1041757313688177D-03,  &
         &     - .1452822267047819D-02,0.8338339968783705D-02,  &
         &     - .1472796329231960D+00/
!
    DATA (DIM2W(I,4),I=1,14)/ - .8264123822525677D+00,  &
         &     0.3065838614094360D+00,0.2389292538329435D-02,  &
         &     - .1343024157997222D+00,0.8833366840533900D-01,  &
         &     0.0000000000000000D+00,0.9786283074168292D-03,  &
         &     - .1319227889147519D+00,0.7990012200150630D-02,  &
         &     0.3391747079760626D-02,0.2294915718283264D-02,  &
         &     - .1358584986119197D-01,0.4025866859057809D-01,  &
         &     0.3760268580063992D-02/
!
    DATA (DIM2W(I,5),I=1,14)/0.6539094339575232D+00,  &
         &     - .2041614154424632D+00, - .1746981515794990D+00,  &
         &     0.3937939671417803D-01,0.6974520545933992D-02,  &
         &     0.0000000000000000D+00,0.6667702171778258D-02,  &
         &     0.5512960621544304D-01,0.5443846381278607D-01,  &
         &     0.2310903863953934D-01,0.1506937747477189D-01,  &
         &     - .6057021648901890D-01,0.4225737654686337D-01,  &
         &     0.2561989142123099D-01/
!
!***FIRST EXECUTABLE STATEMENT D132RE
!
!   Assign values to W.
!
    DO  I = 1,14
       DO  J = 1,5
          W(J,I) = DIM2W(I,J)
       end DO
    end DO
!
!   Assign values to G.
!
    DO I = 1,2
       DO  J = 1,14
          G(I,J) = 0
       end DO
    end DO

    G(1,2) = DIM2G(1)
    G(1,3) = DIM2G(2)
    G(1,4) = DIM2G(3)
    G(1,5) = DIM2G(4)
    G(1,6) = DIM2G(5)
    G(1,7) = DIM2G(6)
    G(2,7) = G(1,7)
    G(1,8) = DIM2G(7)
    G(2,8) = G(1,8)
    G(1,9) = DIM2G(8)
    G(2,9) = G(1,9)
    G(1,10) = DIM2G(9)
    G(2,10) = G(1,10)
    G(1,11) = DIM2G(10)
    G(2,11) = G(1,11)
    G(1,12) = DIM2G(11)
    G(2,12) = DIM2G(12)
    G(1,13) = DIM2G(13)
    G(2,13) = DIM2G(14)
    G(1,14) = DIM2G(15)
    G(2,14) = DIM2G(16)
!
!   Assign values to RULPTS.
!
    RULPTS(1) = 1
    DO 30 I = 2,11
       RULPTS(I) = 4
30  CONTINUE
    RULPTS(12) = 8
    RULPTS(13) = 8
    RULPTS(14) = 8
!
!   Assign values to ERRCOF.
!
    ERRCOF(1) = 10
    ERRCOF(2) = 10
    ERRCOF(3) = 1.
    ERRCOF(4) = 5.
    ERRCOF(5) = 0.5
    ERRCOF(6) = 0.25
!
!***END D132RE
!
    RETURN
  END subroutine D132RE
  SUBROUTINE D113RE(WTLENG,W,G,ERRCOF,RULPTS)
!***BEGIN PROLOGUE D113RE
!***AUTHOR   Jarle Berntsen, EDB-senteret,
!            University of Bergen, Thormohlens gt. 55,
!            N-5008 Bergen, NORWAY
!***PURPOSE D113RE computes abscissas and weights of a 3 dimensional
!            integration rule of degree 11.
!            Two null rules of degree 9, one null rule of degree 7
!            and one null rule of degree 5 to be used in error
!            estimation are also computed.
!***DESCRIPTION D113RE will select the correct values of the abscissas
!            and corresponding weights for different
!            integration rules and null rules and assign them to G
!            and W.
!            The heuristic error coefficients ERRCOF
!            will also be computed.
!
!
!   ON ENTRY
!
!     WTLENG Integer.
!            The number of weights in each of the rules.
!
!   ON RETURN
!
!     W      Real array of dimension (5,WTLENG).
!            The weights for the basic and null rules.
!            W(1,1),...,W(1,WTLENG) are weights for the basic rule.
!            W(I,1),...,W(I,WTLENG), for I > 1 are null rule weights.
!     G      Real array of dimension (NDIM,WTLENG).
!            The fully symmetric sum generators for the rules.
!            G(1,J),...,G(NDIM,J) are the generators for the points
!            associated with the the Jth weights.
!     ERRCOF Real array of dimension 6.
!            Heuristic error coefficients that are used in the
!            error estimation in BASRUL.
!     RULPTS Real array of dimension WTLENG.
!            The number of points used by each generator.
!
!***REFERENCES  J.Berntsen, Cautious adaptive numerical integration
!               over the 3-cube, Reports in Informatics 17, Dept. of
!               Inf.,Univ. of Bergen, Norway, 1985.
!               J.Berntsen and T.O.Espelid, On the construction of
!               higher degree three-dimensional embedded integration
!               rules, SIAM J. Numer. Anal.,Vol. 25,No. 1, pp.222-234,
!               1988.
!***ROUTINES CALLED-NONE
!***END PROLOGUE D113RE
!
!   Global variables.
!
    INTEGER WTLENG
    real(dp) :: W(5,WTLENG),G(3,WTLENG),ERRCOF(6)
    real(dp) :: RULPTS(WTLENG)
!
!   Local variables.
!
    INTEGER I,J
    real(dp) :: DIM3G(14)
    real(dp) :: DIM3W(13,5)
!
    DATA (DIM3G(I),I=1,14)/0.1900000000000000D+00,  &
         &     0.5000000000000000D+00,0.7500000000000000D+00,  &
         &     0.8000000000000000D+00,0.9949999999999999D+00,  &
         &     0.9987344998351400D+00,0.7793703685672423D+00,  &
         &     0.9999698993088767D+00,0.7902637224771788D+00,  &
         &     0.4403396687650737D+00,0.4378478459006862D+00,  &
         &     0.9549373822794593D+00,0.9661093133630748D+00,  &
         &     0.4577105877763134D+00/
!
    DATA (DIM3W(I,1),I=1,13)/0.7923078151105734D-02,  &
         &     0.6797177392788080D-01,0.1086986538805825D-02,  &
         &     0.1838633662212829D+00,0.3362119777829031D-01,  &
         &     0.1013751123334062D-01,0.1687648683985235D-02,  &
         &     0.1346468564512807D+00,0.1750145884600386D-02,  &
         &     0.7752336383837454D-01,0.2461864902770251D+00,  &
         &     0.6797944868483039D-01,0.1419962823300713D-01/  
!
    DATA (DIM3W(I,2),I=1,13)/0.1715006248224684D+01,  &
         &     - .3755893815889209D+00,0.1488632145140549D+00,   &
         &     - .2497046640620823D+00,0.1792501419135204D+00,  &
         &     0.3446126758973890D-02, - .5140483185555825D-02,  &
         &     0.6536017839876425D-02, - .6513454939229700D-03,  &
         &     - .6304672433547204D-02,0.1266959399788263D-01,  &
         &     - .5454241018647931D-02,0.4826995274768427D-02/
!
    DATA (DIM3W(I,3),I=1,13)/0.1936014978949526D+01,  &
         &     - .3673449403754268D+00,0.2929778657898176D-01,  &
         &     - .1151883520260315D+00,0.5086658220872218D-01,  &
         &     0.4453911087786469D-01, - .2287828257125900D-01,  &
         &     0.2908926216345833D-01, - .2898884350669207D-02,  &
         &     - .2805963413307495D-01,0.5638741361145884D-01,  &
         &     - .2427469611942451D-01,0.2148307034182882D-01/
!
    DATA (DIM3W(I,4),I=1,13)/0.5170828195605760D+00,  &
         &     0.1445269144914044D-01, - .3601489663995932D+00,  &
         &     0.3628307003418485D+00,0.7148802650872729D-02,  &
         &     - .9222852896022966D-01,0.1719339732471725D-01,  &
         &     - .1021416537460350D+00, - .7504397861080493D-02,  &
         &     0.1648362537726711D-01,0.5234610158469334D-01,  &
         &     0.1445432331613066D-01,0.3019236275367777D-02/
!
    DATA (DIM3W(I,5),I=1,13)/0.2054404503818520D+01,  &
         &     0.1377759988490120D-01, - .5768062917904410D+00,  &
         &     0.3726835047700328D-01,0.6814878939777219D-02,  &
         &     0.5723169733851849D-01, - .4493018743811285D-01,  &
         &     0.2729236573866348D-01,0.3547473950556990D-03,  &
         &     0.1571366799739551D-01,0.4990099219278567D-01,  &
         &     0.1377915552666770D-01,0.2878206423099872D-02/
!
!***FIRST EXECUTABLE STATEMENT D113RE
!
!   Assign values to W.
!
    DO I = 1,13
       DO  J = 1,5
          W(J,I) = DIM3W(I,J)
       end DO
    end DO

!
!   Assign values to G.
!
    DO  I = 1,3
       DO J = 1,13
          G(I,J) = 0
       end DO
    end DO

    G(1,2) = DIM3G(1)
    G(1,3) = DIM3G(2)
    G(1,4) = DIM3G(3)
    G(1,5) = DIM3G(4)
    G(1,6) = DIM3G(5)
    G(1,7) = DIM3G(6)
    G(2,7) = G(1,7)
    G(1,8) = DIM3G(7)
    G(2,8) = G(1,8)
    G(1,9) = DIM3G(8)
    G(2,9) = G(1,9)
    G(3,9) = G(1,9)
    G(1,10) = DIM3G(9)
    G(2,10) = G(1,10)
    G(3,10) = G(1,10)
    G(1,11) = DIM3G(10)
    G(2,11) = G(1,11)
    G(3,11) = G(1,11)
    G(1,12) = DIM3G(12)
    G(2,12) = DIM3G(11)
    G(3,12) = G(2,12)
    G(1,13) = DIM3G(13)
    G(2,13) = G(1,13)
    G(3,13) = DIM3G(14)
!
!   Assign values to RULPTS.
!
    RULPTS(1) = 1
    RULPTS(2) = 6
    RULPTS(3) = 6
    RULPTS(4) = 6
    RULPTS(5) = 6
    RULPTS(6) = 6
    RULPTS(7) = 12
    RULPTS(8) = 12
    RULPTS(9) = 8
    RULPTS(10) = 8
    RULPTS(11) = 8
    RULPTS(12) = 24
    RULPTS(13) = 24
!
!   Assign values to ERRCOF.
!
    ERRCOF(1) = 4
    ERRCOF(2) = 4.
    ERRCOF(3) = 0.5
    ERRCOF(4) = 3.
    ERRCOF(5) = 0.5
    ERRCOF(6) = 0.25
!
!***END D113RE
!
    RETURN
  END subroutine D113RE
  SUBROUTINE D09HRE(NDIM,WTLENG,W,G,ERRCOF,RULPTS)
!***BEGIN PROLOGUE D09HRE
!***KEYWORDS basic integration rule, degree 9
!***PURPOSE  To initialize a degree 9 basic rule and null rules.
!***AUTHOR   Alan Genz, Computer Science Department, Washington
!            State University, Pullman, WA 99163-1210 USA
!***LAST MODIFICATION 88-05-20
!***DESCRIPTION  D09HRE initializes a degree 9 integration rule,
!            two degree 7 null rules, one degree 5 null rule and one
!            degree 3 null rule for the hypercube [-1,1]**NDIM.
!
!   ON ENTRY
!
!   NDIM   Integer.
!          Number of variables.
!   WTLENG Integer.
!          The number of weights in each of the rules.
!
!   ON RETURN
!   W      Real array of dimension (5,WTLENG).
!          The weights for the basic and null rules.
!          W(1,1),...,W(1,WTLENG) are weights for the basic rule.
!          W(I,1),...,W(I,WTLENG), for I > 1 are null rule weights.
!   G      Real array of dimension (NDIM, WTLENG).
!          The fully symmetric sum generators for the rules.
!          G(1, J), ..., G(NDIM, J) are the are the generators for the
!          points associated with the Jth weights.
!   ERRCOF Real array of dimension 6.
!          Heuristic error coefficients that are used in the
!          error estimation in BASRUL.
!   RULPTS Real array of dimension WTLENG.
!          A work array.
!
!***REFERENCES A. Genz and A. Malik,
!             "An Imbedded Family of Fully Symmetric Numerical
!              Integration Rules",
!              SIAM J Numer. Anal. 20 (1983), pp. 580-588.
!***ROUTINES CALLED-NONE
!***END PROLOGUE D09HRE
!
!   Global variables
!
    INTEGER NDIM,WTLENG
    real(dp) :: W(5,WTLENG),G(NDIM,WTLENG),ERRCOF(6)
    real(dp) :: RULPTS(WTLENG)
!
!   Local Variables
!
    real(dp) :: RATIO,LAM0,LAM1,LAM2,LAM3,LAMP,TWONDM
    INTEGER I,J
!
!***FIRST EXECUTABLE STATEMENT D09HRE
!
!
!     Initialize generators, weights and RULPTS
!
    DO 30 J = 1,WTLENG
       DO 10 I = 1,NDIM
          G(I,J) = 0
10     CONTINUE
       DO 20 I = 1,5
          W(I,J) = 0
20     CONTINUE
       RULPTS(J) = 2*NDIM
30  CONTINUE
    TWONDM = 2**NDIM
    RULPTS(WTLENG) = TWONDM
    IF (NDIM.GT.2) RULPTS(8) = (4*NDIM* (NDIM-1)* (NDIM-2))/3
    RULPTS(7) = 4*NDIM* (NDIM-1)
    RULPTS(6) = 2*NDIM* (NDIM-1)
    RULPTS(1) = 1
!
!     Compute squared generator parameters
!
    LAM0 = 0.4707
    LAM1 = 4/ (15-5/LAM0)
    RATIO = (1-LAM1/LAM0)/27
    LAM2 = (5-7*LAM1-35*RATIO)/ (7-35*LAM1/3-35*RATIO/LAM0)
    RATIO = RATIO* (1-LAM2/LAM0)/3
    LAM3 = (7-9* (LAM2+LAM1)+63*LAM2*LAM1/5-63*RATIO)/  &
         &       (9-63* (LAM2+LAM1)/5+21*LAM2*LAM1-63*RATIO/LAM0)
    LAMP = 0.0625
!
!     Compute degree 9 rule weights
!
    W(1,WTLENG) = 1/ (3*LAM0)**4/TWONDM
    IF (NDIM.GT.2) W(1,8) = (1-1/ (3*LAM0))/ (6*LAM1)**3
    W(1,7) = (1-7* (LAM0+LAM1)/5+7*LAM0*LAM1/3)/  &
         &         (84*LAM1*LAM2* (LAM2-LAM0)* (LAM2-LAM1))
    W(1,6) = (1-7* (LAM0+LAM2)/5+7*LAM0*LAM2/3)/  &
         &         (84*LAM1*LAM1* (LAM1-LAM0)* (LAM1-LAM2)) -  &
         &         W(1,7)*LAM2/LAM1 - 2* (NDIM-2)*W(1,8)
    W(1,4) = (1-9* ((LAM0+LAM1+LAM2)/7- (LAM0*LAM1+LAM0*LAM2+  &
         &         LAM1*LAM2)/5)-3*LAM0*LAM1*LAM2)/  &
         &         (18*LAM3* (LAM3-LAM0)* (LAM3-LAM1)* (LAM3-LAM2))
    W(1,3) = (1-9* ((LAM0+LAM1+LAM3)/7- (LAM0*LAM1+LAM0*LAM3+  &
         &         LAM1*LAM3)/5)-3*LAM0*LAM1*LAM3)/  &
         &         (18*LAM2* (LAM2-LAM0)* (LAM2-LAM1)* (LAM2-LAM3)) -  &
         &         2* (NDIM-1)*W(1,7)
    W(1,2) = (1-9* ((LAM0+LAM2+LAM3)/7- (LAM0*LAM2+LAM0*LAM3+  &
         &         LAM2*LAM3)/5)-3*LAM0*LAM2*LAM3)/  &
         &         (18*LAM1* (LAM1-LAM0)* (LAM1-LAM2)* (LAM1-LAM3)) -  &
         &         2* (NDIM-1)* (W(1,7)+W(1,6)+ (NDIM-2)*W(1,8))
!
!     Compute weights for 2 degree 7, 1 degree 5 and 1 degree 3 rules
!
    W(2,WTLENG) = 1/ (108*LAM0**4)/TWONDM
    IF (NDIM.GT.2) W(2,8) = (1-27*TWONDM*W(2,9)*LAM0**3)/ (6*LAM1)**3
    W(2,7) = (1-5*LAM1/3-15*TWONDM*W(2,WTLENG)*LAM0**2* (LAM0-LAM1))/  &
         &          (60*LAM1*LAM2* (LAM2-LAM1))
    W(2,6) = (1-9* (8*LAM1*LAM2*W(2,7)+TWONDM*W(2,WTLENG)*LAM0**2))/  &
         &         (36*LAM1*LAM1) - 2*W(2,8)* (NDIM-2)
    W(2,4) = (1-7* ((LAM1+LAM2)/5-LAM1*LAM2/3+TWONDM*W(2,  &
         &         WTLENG)*LAM0* (LAM0-LAM1)* (LAM0-LAM2)))/  &
         &         (14*LAM3* (LAM3-LAM1)* (LAM3-LAM2))
    W(2,3) = (1-7* ((LAM1+LAM3)/5-LAM1*LAM3/3+TWONDM*W(2,  &
         &         WTLENG)*LAM0* (LAM0-LAM1)* (LAM0-LAM3)))/  &
         &         (14*LAM2* (LAM2-LAM1)* (LAM2-LAM3)) - 2* (NDIM-1)*W(2,7)
    W(2,2) = (1-7* ((LAM2+LAM3)/5-LAM2*LAM3/3+TWONDM*W(2,  &
         &         WTLENG)*LAM0* (LAM0-LAM2)* (LAM0-LAM3)))/  &
         &         (14*LAM1* (LAM1-LAM2)* (LAM1-LAM3)) -  &
         &         2* (NDIM-1)* (W(2,7)+W(2,6)+ (NDIM-2)*W(2,8))  
    W(3,WTLENG) = 5/ (324*LAM0**4)/TWONDM
    IF (NDIM.GT.2) W(3,8) = (1-27*TWONDM*W(3,9)*LAM0**3)/ (6*LAM1)**3
    W(3,7) = (1-5*LAM1/3-15*TWONDM*W(3,WTLENG)*LAM0**2* (LAM0-LAM1))/  &
         &          (60*LAM1*LAM2* (LAM2-LAM1))
    W(3,6) = (1-9* (8*LAM1*LAM2*W(3,7)+TWONDM*W(3,WTLENG)*LAM0**2))/  &
         &         (36*LAM1*LAM1) - 2*W(3,8)* (NDIM-2)
    W(3,5) = (1-7* ((LAM1+LAM2)/5-LAM1*LAM2/3+TWONDM*W(3,  &
         &         WTLENG)*LAM0* (LAM0-LAM1)* (LAM0-LAM2)))/  &
         &         (14*LAMP* (LAMP-LAM1)* (LAMP-LAM2))
    W(3,3) = (1-7* ((LAM1+LAMP)/5-LAM1*LAMP/3+TWONDM*W(3,  &
         &         WTLENG)*LAM0* (LAM0-LAM1)* (LAM0-LAMP)))/  &
         &         (14*LAM2* (LAM2-LAM1)* (LAM2-LAMP)) - 2* (NDIM-1)*W(3,7)
    W(3,2) = (1-7* ((LAM2+LAMP)/5-LAM2*LAMP/3+TWONDM*W(3,  &
         &         WTLENG)*LAM0* (LAM0-LAM2)* (LAM0-LAMP)))/  &
         &         (14*LAM1* (LAM1-LAM2)* (LAM1-LAMP)) -  &
         &         2* (NDIM-1)* (W(3,7)+W(3,6)+ (NDIM-2)*W(3,8))
    W(4,WTLENG) = 2/ (81*LAM0**4)/TWONDM
    IF (NDIM.GT.2) W(4,8) = (2-27*TWONDM*W(4,9)*LAM0**3)/ (6*LAM1)**3
    W(4,7) = (2-15*LAM1/9-15*TWONDM*W(4,WTLENG)*LAM0* (LAM0-LAM1))/  &
         &         (60*LAM1*LAM2* (LAM2-LAM1))
    W(4,6) = (1-9* (8*LAM1*LAM2*W(4,7)+TWONDM*W(4,WTLENG)*LAM0**2))/  &
         &         (36*LAM1*LAM1) - 2*W(4,8)* (NDIM-2)
    W(4,4) = (2-7* ((LAM1+LAM2)/5-LAM1*LAM2/3+TWONDM*W(4,  &
         &         WTLENG)*LAM0* (LAM0-LAM1)* (LAM0-LAM2)))/  &
         &         (14*LAM3* (LAM3-LAM1)* (LAM3-LAM2))
    W(4,3) = (2-7* ((LAM1+LAM3)/5-LAM1*LAM3/3+TWONDM*W(4,  &
         &         WTLENG)*LAM0* (LAM0-LAM1)* (LAM0-LAM3)))/  &
         &         (14*LAM2* (LAM2-LAM1)* (LAM2-LAM3)) - 2* (NDIM-1)*W(4,7)
    W(4,2) = (2-7* ((LAM2+LAM3)/5-LAM2*LAM3/3+TWONDM*W(4,  &
         &         WTLENG)*LAM0* (LAM0-LAM2)* (LAM0-LAM3)))/  &
         &         (14*LAM1* (LAM1-LAM2)* (LAM1-LAM3)) -  &
         &         2* (NDIM-1)* (W(4,7)+W(4,6)+ (NDIM-2)*W(4,8))
    W(5,2) = 1/ (6*LAM1)
!
!     Set generator values
!
    LAM0 = SQRT(LAM0)
    LAM1 = SQRT(LAM1)
    LAM2 = SQRT(LAM2)
    LAM3 = SQRT(LAM3)
    LAMP = SQRT(LAMP)
    DO 40 I = 1,NDIM
       G(I,WTLENG) = LAM0
40  CONTINUE
    IF (NDIM.GT.2) THEN
       G(1,8) = LAM1
       G(2,8) = LAM1
       G(3,8) = LAM1
    END IF
    G(1,7) = LAM1
    G(2,7) = LAM2
    G(1,6) = LAM1
    G(2,6) = LAM1
    G(1,5) = LAMP
    G(1,4) = LAM3
    G(1,3) = LAM2
    G(1,2) = LAM1
!
!     Compute final weight values.
!     The null rule weights are computed from differences between
!     the degree 9 rule weights and lower degree rule weights.
!
    W(1,1) = TWONDM
    DO 70 J = 2,5
       DO 50 I = 2,WTLENG
          W(J,I) = W(J,I) - W(1,I)
          W(J,1) = W(J,1) - RULPTS(I)*W(J,I)
50     CONTINUE
70  CONTINUE
    DO 80 I = 2,WTLENG
       W(1,I) = TWONDM*W(1,I)
       W(1,1) = W(1,1) - RULPTS(I)*W(1,I)
80  CONTINUE
!
!     Set error coefficients
!
    ERRCOF(1) = 5
    ERRCOF(2) = 5
    ERRCOF(3) = 1.
    ERRCOF(4) = 5
    ERRCOF(5) = 0.5
    ERRCOF(6) = 0.25
!
!***END D09HRE
!
  END subroutine D09HRE
  SUBROUTINE D07HRE(NDIM,WTLENG,W,G,ERRCOF,RULPTS)
!***BEGIN PROLOGUE D07HRE
!***KEYWORDS basic integration rule, degree 7
!***PURPOSE  To initialize a degree 7 basic rule, and null rules.
!***AUTHOR   Alan Genz, Computer Science Department, Washington
!            State University, Pullman, WA 99163-1210 USA
!***LAST MODIFICATION 88-05-31
!***DESCRIPTION  D07HRE initializes a degree 7 integration rule,
!            two degree 5 null rules, one degree 3 null rule and one
!            degree 1 null rule for the hypercube [-1,1]**NDIM.
!
!   ON ENTRY
!
!   NDIM   Integer.
!          Number of variables.
!   WTLENG Integer.
!          The number of weights in each of the rules.
!          WTLENG MUST be set equal to 6.
!
!   ON RETURN
!   W      Real array of dimension (5,WTLENG).
!          The weights for the basic and null rules.
!          W(1,1),...,W(1,WTLENG) are weights for the basic rule.
!          W(I,1),...,W(I,WTLENG), for I > 1 are null rule weights.
!   G      Real array of dimension (NDIM, WTLENG).
!          The fully symmetric sum generators for the rules.
!          G(1, J), ..., G(NDIM, J) are the are the generators for the
!          points associated with the Jth weights.
!   ERRCOF Real array of dimension 6.
!          Heuristic error coefficients that are used in the
!          error estimation in BASRUL.
!   RULPTS Real array of dimension WTLENG.
!          A work array.
!
!***REFERENCES A. Genz and A. Malik,
!             "An Imbedded Family of Fully Symmetric Numerical
!              Integration Rules",
!              SIAM J Numer. Anal. 20 (1983), pp. 580-588.
!***ROUTINES CALLED-NONE
!***END PROLOGUE D07HRE
!
!   Global variables
!
    INTEGER NDIM,WTLENG
    real(dp) :: W(5,WTLENG),G(NDIM,WTLENG),ERRCOF(6)
    real(dp) :: RULPTS(WTLENG)
!
!   Local Variables
!
    real(dp) :: RATIO,LAM0,LAM1,LAM2,LAMP,TWONDM
    INTEGER I,J
!
!***FIRST EXECUTABLE STATEMENT D07HRE
!
!
!     Initialize generators, weights and RULPTS
!
    DO 30 J = 1,WTLENG
       DO 10 I = 1,NDIM
          G(I,J) = 0
10     CONTINUE
       DO 20 I = 1,5
          W(I,J) = 0
20     CONTINUE
       RULPTS(J) = 2*NDIM
30  CONTINUE
    TWONDM = 2**NDIM
    RULPTS(WTLENG) = TWONDM
    RULPTS(WTLENG-1) = 2*NDIM* (NDIM-1)
    RULPTS(1) = 1
!
!     Compute squared generator parameters
!
    LAM0 = 0.4707
    LAMP = 0.5625
    LAM1 = 4/ (15-5/LAM0)
    RATIO = (1-LAM1/LAM0)/27
    LAM2 = (5-7*LAM1-35*RATIO)/ (7-35*LAM1/3-35*RATIO/LAM0)
!
!     Compute degree 7 rule weights
!
    W(1,6) = 1/ (3*LAM0)**3/TWONDM
    W(1,5) = (1-5*LAM0/3)/ (60* (LAM1-LAM0)*LAM1**2)
    W(1,3) = (1-5*LAM2/3-5*TWONDM*W(1,6)*LAM0* (LAM0-LAM2))/  &
         &         (10*LAM1* (LAM1-LAM2)) - 2* (NDIM-1)*W(1,5)
    W(1,2) = (1-5*LAM1/3-5*TWONDM*W(1,6)*LAM0* (LAM0-LAM1))/  &
         &         (10*LAM2* (LAM2-LAM1))
!
!     Compute weights for 2 degree 5, 1 degree 3 and 1 degree 1 rules
!
    W(2,6) = 1/ (36*LAM0**3)/TWONDM
    W(2,5) = (1-9*TWONDM*W(2,6)*LAM0**2)/ (36*LAM1**2)
    W(2,3) = (1-5*LAM2/3-5*TWONDM*W(2,6)*LAM0* (LAM0-LAM2))/  &
         &         (10*LAM1* (LAM1-LAM2)) - 2* (NDIM-1)*W(2,5)
    W(2,2) = (1-5*LAM1/3-5*TWONDM*W(2,6)*LAM0* (LAM0-LAM1))/  &
         &         (10*LAM2* (LAM2-LAM1))
    W(3,6) = 5/ (108*LAM0**3)/TWONDM
    W(3,5) = (1-9*TWONDM*W(3,6)*LAM0**2)/ (36*LAM1**2)
    W(3,3) = (1-5*LAMP/3-5*TWONDM*W(3,6)*LAM0* (LAM0-LAMP))/  &
         &         (10*LAM1* (LAM1-LAMP)) - 2* (NDIM-1)*W(3,5)
    W(3,4) = (1-5*LAM1/3-5*TWONDM*W(3,6)*LAM0* (LAM0-LAM1))/  &
         &         (10*LAMP* (LAMP-LAM1))
    W(4,6) = 1/ (54*LAM0**3)/TWONDM
    W(4,5) = (1-18*TWONDM*W(4,6)*LAM0**2)/ (72*LAM1**2)
    W(4,3) = (1-10*LAM2/3-10*TWONDM*W(4,6)*LAM0* (LAM0-LAM2))/  &
         &         (20*LAM1* (LAM1-LAM2)) - 2* (NDIM-1)*W(4,5)
    W(4,2) = (1-10*LAM1/3-10*TWONDM*W(4,6)*LAM0* (LAM0-LAM1))/  &
         &         (20*LAM2* (LAM2-LAM1))
!
!     Set generator values
!
    LAM0 = SQRT(LAM0)
    LAM1 = SQRT(LAM1)
    LAM2 = SQRT(LAM2)
    LAMP = SQRT(LAMP)
    DO 40 I = 1,NDIM
       G(I,WTLENG) = LAM0
40  CONTINUE
    G(1,WTLENG-1) = LAM1
    G(2,WTLENG-1) = LAM1
    G(1,WTLENG-4) = LAM2
    G(1,WTLENG-3) = LAM1
    G(1,WTLENG-2) = LAMP
!
!     Compute final weight values.
!     The null rule weights are computed from differences between
!     the degree 7 rule weights and lower degree rule weights.
!
    W(1,1) = TWONDM
    DO 70 J = 2,5
       DO 50 I = 2,WTLENG
          W(J,I) = W(J,I) - W(1,I)
          W(J,1) = W(J,1) - RULPTS(I)*W(J,I)
50     CONTINUE
70  CONTINUE
    DO 80 I = 2,WTLENG
       W(1,I) = TWONDM*W(1,I)
       W(1,1) = W(1,1) - RULPTS(I)*W(1,I)
80  CONTINUE
!
!     Set error coefficients
!
    ERRCOF(1) = 5
    ERRCOF(2) = 5
    ERRCOF(3) = 1
    ERRCOF(4) = 5
    ERRCOF(5) = 0.5
    ERRCOF(6) = 0.25
!
!***END D07HRE
!
  END subroutine D07HRE
  SUBROUTINE DRLHRE(NDIM,CENTER,HWIDTH,WTLENG,G,W,ERRCOF,NUMFUN,  &
       &                  FUNSUB,SCALES,NORMS,X,NULL,BASVAL,RGNERR,DIRECT)
!***BEGIN PROLOGUE DRLHRE
!***KEYWORDS basic numerical integration rule
!***PURPOSE  To compute basic integration rule values.
!***AUTHOR   Alan Genz, Computer Science Department, Washington
!            State University, Pullman, WA 99163-1210 USA
!***LAST MODIFICATION 90-02-06
!***DESCRIPTION DRLHRE computes basic integration rule values for a
!            vector of integrands over a hyper-rectangular region.
!            These are estimates for the integrals. DRLHRE also computes
!            estimates for the errors and determines the coordinate axis
!            where the fourth difference for the integrands is largest.
!
!   ON ENTRY
!
!   NDIM   Integer.
!          Number of variables.
!   CENTER Real array of dimension NDIM.
!          The coordinates for the center of the region.
!   HWIDTH Real Array of dimension NDIM.
!          HWIDTH(I) is half of the width of dimension I of the region.
!   WTLENG Integer.
!          The number of weights in the basic integration rule.
!   G      Real array of dimension (NDIM,WTLENG).
!          The fully symmetric sum generators for the rules.
!          G(1,J), ..., G(NDIM,J) are the are the generators for the
!          points associated with the Jth weights.
!   W      Real array of dimension (5,WTLENG).
!          The weights for the basic and null rules.
!          W(1,1),...,W(1,WTLENG) are weights for the basic rule.
!          W(I,1),...,W(I,WTLENG), for I > 1 are null rule weights.
!   ERRCOF Real array of dimension 6.
!          The error coefficients for the rules.
!          It is assumed that the error is computed using:
!           IF (N1*ERRCOF(1) < N2 and N2*ERRCOF(2) < N3)
!             THEN ERROR = ERRCOF(3)*N1
!             ELSE ERROR = ERRCOF(4)*MAX(N1, N2, N3)
!           ERROR = ERROR + EP*(ERRCOF(5)*ERROR/(ES+ERROR)+ERRCOF(6))
!          where N1-N4 are the null rules, EP is the error
!          for the parent
!          subregion and ES is the error for the sibling subregion.
!   NUMFUN Integer.
!          Number of components for the vector integrand.
!   FUNSUB Externally declared subroutine.
!          For computing the components of the integrand at a point X.
!          It must have parameters (NDIM,X,NUMFUN,FUNVLS).
!           Input Parameters:
!            X      Real array of dimension NDIM.
!                   Defines the evaluation point.
!            NDIM   Integer.
!                   Number of variables for the integrand.
!            NUMFUN Integer.
!                   Number of components for the vector integrand.
!           Output Parameters:
!            FUNVLS Real array of dimension NUMFUN.
!                   The components of the integrand at the point X.
!   SCALES Real array of dimension (3,WTLENG).
!          Scaling factors used to construct new null rules based
!          on a linear combination of two successive null rules
!          in the sequence of null rules.
!   NORMS  Real array of dimension (3,WTLENG).
!          2**NDIM/(1-norm of the null rule constructed by each of the
!          scaling factors.)
!   X      Real Array of dimension NDIM.
!          A work array.
!   NULL   Real array of dimension (NUMFUN, 8)
!          A work array.
!
!   ON RETURN
!
!   BASVAL Real array of dimension NUMFUN.
!          The values for the basic rule for each component
!          of the integrand.
!   RGNERR Real array of dimension NUMFUN.
!          The error estimates for each component of the integrand.
!   DIRECT Real.
!          The coordinate axis where the fourth difference of the
!          integrand values is largest.
!
!***REFERENCES
!   A.C.Genz and A.A.Malik, An adaptive algorithm for numerical
!   integration over an N-dimensional rectangular region,
!   J.Comp.Appl.Math., 6:295-302, 1980.
!
!   T.O.Espelid, Integration Rules, Null Rules and Error
!   Estimation, Reports in Informatics 33, Dept. of Informatics,
!   Univ. of Bergen, 1988.
!
!***ROUTINES CALLED: DFSHRE, FUNSUB
!
!***END PROLOGUE DRLHRE
!
!   Global variables.
!
!      EXTERNAL FUNSUB
    interface
       function FUNSUB(Z)
         use nrtype
         real(dp), intent(in) :: Z(:)
         real(dp) :: FUNSUB
       end function FUNSUB
    end interface
    INTEGER WTLENG,NUMFUN,NDIM
    real(dp) :: CENTER(NDIM),X(NDIM),HWIDTH(NDIM),BASVAL(NUMFUN),  &
         &                 RGNERR(NUMFUN),NULL(NUMFUN,8),W(5,WTLENG),  &
         &                 G(NDIM,WTLENG),ERRCOF(6),DIRECT,SCALES(3,WTLENG),  &
         &                 NORMS(3,WTLENG)
!
!   Local variables.
!
    real(dp) :: RGNVOL,DIFSUM,DIFMAX,FRTHDF
    INTEGER I,J,K,DIVAXN
    real(dp) :: SEARCH,RATIO
!
!***FIRST EXECUTABLE STATEMENT DRLHRE
!
!
!       Compute volume of subregion, initialize DIVAXN and rule sums;
!       compute fourth differences and new DIVAXN (RGNERR is used
!       for a work array here). The integrand values used for the
!       fourth divided differences are accumulated in rule arrays.
!
    RGNVOL = 1
    DIVAXN = 1
    DO 10 I = 1,NDIM
       RGNVOL = RGNVOL*HWIDTH(I)
       X(I) = CENTER(I)
       IF (HWIDTH(I).GT.HWIDTH(DIVAXN)) DIVAXN = I
10  CONTINUE
    RGNERR(1) = FUNSUB(X)
    DO 30 J = 1,NUMFUN
       BASVAL(J) = W(1,1)*RGNERR(J)
       DO 20 K = 1,4
          NULL(J,K) = W(K+1,1)*RGNERR(J)
20     CONTINUE
30  CONTINUE
    DIFMAX = 0
    RATIO = (G(1,3)/G(1,2))**2
    DO 60 I = 1,NDIM
       X(I) = CENTER(I) - HWIDTH(I)*G(1,2)
       NULL(1,5) = FUNSUB(X)
!          CALL FUNSUB(NDIM,X,NUMFUN,NULL(1,5))
       X(I) = CENTER(I) + HWIDTH(I)*G(1,2)
       NULL(1,6) = FUNSUB(X)
!          CALL FUNSUB(NDIM,X,NUMFUN,NULL(1,6))
       X(I) = CENTER(I) - HWIDTH(I)*G(1,3)
       NULL(1,7) = FUNSUB(X)
!          CALL FUNSUB(NDIM,X,NUMFUN,NULL(1,7))
       X(I) = CENTER(I) + HWIDTH(I)*G(1,3)
       NULL(1,8) = FUNSUB(X)
!          CALL FUNSUB(NDIM,X,NUMFUN,NULL(1,8))
       X(I) = CENTER(I)
       DIFSUM = 0
       DO 50 J = 1,NUMFUN
          FRTHDF = 2* (1-RATIO)*RGNERR(J) - (NULL(J,7)+NULL(J,8)) +  &
               &                 RATIO* (NULL(J,5)+NULL(J,6))
!
!           Ignore differences below roundoff
!
          IF (RGNERR(J)+FRTHDF/4.NE.RGNERR(J)) DIFSUM = DIFSUM +  &
               &            ABS(FRTHDF)
          DO 40 K = 1,4
             NULL(J,K) = NULL(J,K) + W(K+1,2)*  &
                  &      (NULL(J,5)+NULL(J,6)) +  &
                  &      W(K+1,3)* (NULL(J,7)+NULL(J,8))
40        CONTINUE
          BASVAL(J) = BASVAL(J) + W(1,2)* (NULL(J,5)+NULL(J,6)) +  &
               &     W(1,3)* (NULL(J,7)+NULL(J,8))
50     CONTINUE
       IF (DIFSUM.GT.DIFMAX) THEN
          DIFMAX = DIFSUM
          DIVAXN = I
       END IF
60  CONTINUE
    DIRECT = DIVAXN
!
!    Finish computing the rule values.
!
    DO 90 I = 4,WTLENG
       CALL DFSHRE(NDIM,CENTER,HWIDTH,X,G(1,I),NUMFUN,FUNSUB,RGNERR,  &
            &                NULL(1,5))
       DO 80 J = 1,NUMFUN
          BASVAL(J) = BASVAL(J) + W(1,I)*RGNERR(J)
          DO 70 K = 1,4
             NULL(J,K) = NULL(J,K) + W(K+1,I)*RGNERR(J)
70        CONTINUE
80     CONTINUE
90  CONTINUE
!
!    Compute errors.
!
    DO 130 J = 1,NUMFUN
!
!    We search for the null rule, in the linear space spanned by two
!    successive null rules in our sequence, which gives the greatest
!    error estimate among all normalized (1-norm) null rules in this
!    space.
!
       DO 110 I = 1,3
          SEARCH = 0
          DO 100 K = 1,WTLENG
             SEARCH = MAX(SEARCH,ABS(NULL(J,I+1)+SCALES(I,  &
                  &       K)*NULL(J,I))*NORMS(I,K))
100       CONTINUE
          NULL(J,I) = SEARCH
110    CONTINUE
       IF (ERRCOF(1)*NULL(J,1).LE.NULL(J,2) .AND.  &
            &  ERRCOF(2)*NULL(J,2).LE.NULL(J,3)) THEN
          RGNERR(J) = ERRCOF(3)*NULL(J,1)
       ELSE
          RGNERR(J) = ERRCOF(4)*MAX(NULL(J,1),NULL(J,2),NULL(J,3))
       END IF
       RGNERR(J) = RGNVOL*RGNERR(J)
       BASVAL(J) = RGNVOL*BASVAL(J)
130 CONTINUE
!
!***END DRLHRE
!
  END subroutine DRLHRE
  SUBROUTINE DFSHRE(NDIM,CENTER,HWIDTH,X,G,NUMFUN,FUNSUB,FULSMS,  &
          &            FUNVLS)
!***BEGIN PROLOGUE DFSHRE
!***KEYWORDS fully symmetric sum
!***PURPOSE  To compute fully symmetric basic rule sums
!***AUTHOR   Alan Genz, Computer Science Department, Washington
!            State University, Pullman, WA 99163-1210 USA
!***LAST MODIFICATION 88-04-08
!***DESCRIPTION DFSHRE computes a fully symmetric sum for a vector
!            of integrand values over a hyper-rectangular region.
!            The sum is fully symmetric with respect to the center of
!            the region and is taken over all sign changes and
!            permutations of the generators for the sum.
!
!   ON ENTRY
!
!   NDIM   Integer.
!          Number of variables.
!   CENTER Real array of dimension NDIM.
!          The coordinates for the center of the region.
!   HWIDTH Real Array of dimension NDIM.
!          HWIDTH(I) is half of the width of dimension I of the region.
!   X      Real Array of dimension NDIM.
!          A work array.
!   G      Real Array of dimension NDIM.
!          The generators for the fully symmetric sum. These MUST BE
!          non-negative and non-increasing.
!   NUMFUN Integer.
!          Number of components for the vector integrand.
!   FUNSUB Externally declared subroutine.
!          For computing the components of the integrand at a point X.
!          It must have parameters (NDIM, X, NUMFUN, FUNVLS).
!           Input Parameters:
!            X      Real array of dimension NDIM.
!                   Defines the evaluation point.
!            NDIM   Integer.
!                   Number of variables for the integrand.
!            NUMFUN Integer.
!                   Number of components for the vector integrand.
!           Output Parameters:
!            FUNVLS Real array of dimension NUMFUN.
!                   The components of the integrand at the point X.
!   ON RETURN
!
!   FULSMS Real array of dimension NUMFUN.
!          The values for the fully symmetric sums for each component
!          of the integrand.
!   FUNVLS Real array of dimension NUMFUN.
!          A work array.
!
!***ROUTINES CALLED: FUNSUB
!
!***END PROLOGUE DFSHRE
!
!   Global variables.
!
!      EXTERNAL FUNSUB
    interface
       function FUNSUB(Z)
         use nrtype
         real(dp), intent(in) :: Z(:)
         real(dp) ::FUNSUB
       end function FUNSUB
    end interface
    INTEGER NDIM,NUMFUN
    real(dp) :: CENTER(NDIM),HWIDTH(NDIM),X(NDIM),G(NDIM),  &
         &                 FULSMS(NUMFUN),FUNVLS(NUMFUN)
!
!   Local variables.
!
    INTEGER IXCHNG,LXCHNG,I,J,L
    real(dp) :: GL,GI
!
!***FIRST EXECUTABLE STATEMENT DFSHRE
!
    DO 10 J = 1,NUMFUN
       FULSMS(J) = 0
10  CONTINUE
!
!     Compute centrally symmetric sum for permutation of G
!
20  DO 30 I = 1,NDIM
       X(I) = CENTER(I) + G(I)*HWIDTH(I)
30  CONTINUE
!40    CALL FUNSUB(NDIM,X,NUMFUN,FUNVLS)
40  FUNVLS(1) = FUNSUB(X)
    DO 50 J = 1,NUMFUN
       FULSMS(J) = FULSMS(J) + FUNVLS(J)
50  CONTINUE
    DO 60 I = 1,NDIM
       G(I) = - G(I)
       X(I) = CENTER(I) + G(I)*HWIDTH(I)
       IF (G(I).LT.0) GO TO 40
60  CONTINUE
!
!       Find next distinct permutation of G and loop back for next sum.
!       Permutations are generated in reverse lexicographic order.
!
    DO 80 I = 2,NDIM
       IF (G(I-1).GT.G(I)) THEN
          GI = G(I)
          IXCHNG = I - 1
          DO 70 L = 1, (I-1)/2
             GL = G(L)
             G(L) = G(I-L)
             G(I-L) = GL
             IF (GL.LE.GI) IXCHNG = IXCHNG - 1
             IF (G(L).GT.GI) LXCHNG = L
70        CONTINUE
          IF (G(IXCHNG).LE.GI) IXCHNG = LXCHNG
          G(I) = G(IXCHNG)
          G(IXCHNG) = GI
          GO TO 20
       END IF
80  CONTINUE
!
!     Restore original order to generators
!
    DO 90 I = 1,NDIM/2
       GI = G(I)
       G(I) = G(NDIM-I+1)
       G(NDIM-I+1) = GI
90  CONTINUE
!
!***END DFSHRE
!
  END subroutine DFSHRE
  SUBROUTINE DTRHRE(DVFLAG,NDIM,NUMFUN,SBRGNS,VALUES,ERRORS,CENTRS,  &
     &                  HWIDTS,GREATE,ERROR,VALUE,CENTER,HWIDTH,DIR)
!***BEGIN PROLOGUE DTRHRE
!***PURPOSE DTRHRE maintains a heap of subregions.
!***DESCRIPTION DTRHRE maintains a heap of subregions.
!            The subregions are ordered according to the size
!            of the greatest error estimates of each subregion(GREATE).
!
!   PARAMETERS
!
!     DVFLAG Integer.
!            If DVFLAG = 1, we remove the subregion with
!            greatest error from the heap.
!            If DVFLAG = 2, we insert a new subregion in the heap.
!     NDIM   Integer.
!            Number of variables.
!     NUMFUN Integer.
!            Number of components of the integral.
!     SBRGNS Integer.
!            Number of subregions in the heap.
!     VALUES Real array of dimension (NUMFUN,SBRGNS).
!            Used to store estimated values of the integrals
!            over the subregions.
!     ERRORS Real array of dimension (NUMFUN,SBRGNS).
!            Used to store the corresponding estimated errors.
!     CENTRS Real array of dimension (NDIM,SBRGNS).
!            Used to store the center limits of the stored subregions.
!     HWIDTS Real array of dimension (NDIM,SBRGNS).
!            Used to store the hwidth limits of the stored subregions.
!     GREATE Real array of dimension SBRGNS.
!            Used to store the greatest estimated errors in
!            all subregions.
!     ERROR  Real array of dimension NUMFUN.
!            Used as intermediate storage for the error of a subregion.
!     VALUE  Real array of dimension NUMFUN.
!            Used as intermediate storage for the estimate
!            of the integral over a subregion.
!     CENTER Real array of dimension NDIM.
!            Used as intermediate storage for the center of
!            the subregion.
!     HWIDTH Real array of dimension NDIM.
!            Used as intermediate storage for the half width of
!            the subregion.
!     DIR    Integer array of dimension SBRGNS.
!            DIR is used to store the directions for
!            further subdivision.
!
!***ROUTINES CALLED-NONE
!***END PROLOGUE DTRHRE
!
!   Global variables.
!
    INTEGER :: DVFLAG,NDIM,NUMFUN,SBRGNS
    real(dp) :: VALUES(NUMFUN,*),ERRORS(NUMFUN,*)
    real(dp) :: CENTRS(NDIM,*)
    real(dp) :: HWIDTS(NDIM,*)
    real(dp) :: GREATE(*)
    real(dp) :: ERROR(NUMFUN),VALUE(NUMFUN)
    real(dp) :: CENTER(NDIM),HWIDTH(NDIM)
    real(dp) :: DIR(*)
!
!   Local variables.
!
!   GREAT  is used as intermediate storage for the greatest error of a
!          subregion.
!   DIRECT is used as intermediate storage for the direction of further
!          subdivision.
!   SUBRGN Position of child/parent subregion in the heap.
!   SUBTMP Position of parent/child subregion in the heap.
!
    INTEGER J,SUBRGN,SUBTMP
    real(dp) :: GREAT,DIRECT
!
!***FIRST EXECUTABLE STATEMENT DTRHRE
!
!   Save values to be stored in their correct place in the heap.
!
    GREAT = GREATE(SBRGNS)
    DIRECT = DIR(SBRGNS)
    DO 5 J = 1,NUMFUN
       ERROR(J) = ERRORS(J,SBRGNS)
       VALUE(J) = VALUES(J,SBRGNS)
5   CONTINUE
    DO 10 J = 1,NDIM
       CENTER(J) = CENTRS(J,SBRGNS)
       HWIDTH(J) = HWIDTS(J,SBRGNS)
10  CONTINUE
!
!    If DVFLAG = 1, we will remove the region
!    with greatest estimated error from the heap.
!
     IF (DVFLAG.EQ.1) THEN
        SBRGNS = SBRGNS - 1
        SUBRGN = 1
20      SUBTMP = 2*SUBRGN
        IF (SUBTMP.LE.SBRGNS) THEN
           IF (SUBTMP.NE.SBRGNS) THEN
!
!   Find max. of left and right child.
!
              IF (GREATE(SUBTMP).LT.GREATE(SUBTMP+1)) THEN
                 SUBTMP = SUBTMP + 1
              END IF
           END IF
!
!   Compare max.child with parent.
!   If parent is max., then done.
!
           IF (GREAT.LT.GREATE(SUBTMP)) THEN
!
!   Move the values at position subtmp up the heap.
!
              GREATE(SUBRGN) = GREATE(SUBTMP)
              DO 25 J = 1,NUMFUN
                 ERRORS(J,SUBRGN) = ERRORS(J,SUBTMP)
                 VALUES(J,SUBRGN) = VALUES(J,SUBTMP)
25            CONTINUE
              DIR(SUBRGN) = DIR(SUBTMP)
              DO 30 J = 1,NDIM
                 CENTRS(J,SUBRGN) = CENTRS(J,SUBTMP)
                 HWIDTS(J,SUBRGN) = HWIDTS(J,SUBTMP)
30            CONTINUE
                 SUBRGN = SUBTMP
                 GO TO 20
           END IF
        END IF
     ELSE IF (DVFLAG.EQ.2) THEN
!
!   If DVFLAG = 2, then insert new region in the heap.
!
        SUBRGN = SBRGNS
40      SUBTMP = SUBRGN/2
        IF (SUBTMP.GE.1) THEN
!
!   Compare max.child with parent.
!   If parent is max, then done.
!
           IF (GREAT.GT.GREATE(SUBTMP)) THEN
!
!   Move the values at position subtmp down the heap.
!
              GREATE(SUBRGN) = GREATE(SUBTMP)
              DO 45 J = 1,NUMFUN
                 ERRORS(J,SUBRGN) = ERRORS(J,SUBTMP)
                 VALUES(J,SUBRGN) = VALUES(J,SUBTMP)
45            CONTINUE
              DIR(SUBRGN) = DIR(SUBTMP)
              DO 50 J = 1,NDIM
                 CENTRS(J,SUBRGN) = CENTRS(J,SUBTMP)
                 HWIDTS(J,SUBRGN) = HWIDTS(J,SUBTMP)
50            CONTINUE
              SUBRGN = SUBTMP
              GO TO 40
           END IF
        END IF
     END IF
!
!    Insert the saved values in their correct places.
!
     IF (SBRGNS.GT.0) THEN
        GREATE(SUBRGN) = GREAT
        DO 55 J = 1,NUMFUN
           ERRORS(J,SUBRGN) = ERROR(J)
           VALUES(J,SUBRGN) = VALUE(J)
55      CONTINUE
        DIR(SUBRGN) = DIRECT
        DO 60 J = 1,NDIM
           CENTRS(J,SUBRGN) = CENTER(J)
           HWIDTS(J,SUBRGN) = HWIDTH(J)
60      CONTINUE
     END IF
!
!***END DTRHRE
!
    RETURN
  END subroutine DTRHRE
  
end module DCUHRE90

!!$module testfunc
!!$ 
!!$contains
!!$  
!!$  function FTEST(Z)
!!$    USE mykinds
!!$    INTEGER N
!!$    integer :: NDIM = 5
!!$    real(dp) :: Z(:), FTEST, SUM
!!$    SUM = 0
!!$    DO N = 1,NDIM
!!$       SUM = SUM + N*Z(N)**2
!!$    end DO
!!$    !   10 CONTINUE
!!$    FTEST = EXP(-SUM/2)
!!$!    DO N = 1,NDIM
!!$!       F(N+1) = Z(N)*F(1)
!!$    !   20 CONTINUE
!!$!    end DO
!!$  END function FTEST
!!$end module testfunc
!     Code for the paper "Calculating the proportion of triangles in 
!                         a Poisson-Voronoi tessellation of the plane"
!                         by A.Hayen and M.Quine
!     Andrew Hayen
!     School of Mathematics and Statistics 
!     University of Sydney NSW 2006
!     Australia
!     email: andrewh@maths.usyd.edu.au
!     Malcolm Quine
!     School of Mathematics and Statistics 
!     University of Sydney NSW 2006
!     Australia
!     email: malcolmq@maths.usyd.edu.au
!     Web: http://www.maths.usyd.edu.au/u/malcolmq/
!$$$
!$$$      PROGRAM DTEST1
!$$$      EXTERNAL FTEST
!$$$      INTEGER KEY, N, NF, NDIM, MINCLS, MAXCLS, IFAIL, NEVAL, NW
!$$$      PARAMETER (NDIM = 4, NW = 7500000, NF = 1)
!$$$      DOUBLE PRECISION A(NDIM), B(NDIM), WRKSTR(NW)
!$$$      DOUBLE PRECISION ABSEST(NF), FINEST(NF), ABSREQ, RELREQ
!$$$
!$$$
!$$$
!$$$c     maximum calls:
!$$$      MAXCLS=80000000
!$$$      KEY = 0
!$$$      ABSREQ = 0
!$$$      RELREQ = 1E-9
!$$$      ERROR=0.
!$$$      ANSWER=0.
!$$$      A1=1.
!$$$      A2=1.
!$$$      A3=1.
!$$$      A4=1.
!$$$      DO 10 n1 = 1,A1
!$$$      DO 10 N2 = 1,A2
!$$$      DO 10 N3 = 1,A3
!$$$      DO 10 N4 = 1,A4
!$$$
!$$$
!$$$
!$$$         
!$$$         a(1)=(n1-1)*(1/a1)
!$$$         b(1)=1/a1+A(1)
!$$$         a(2)=(n2-1)*(1/a2)
!$$$         b(2)=1/a2+A(2)
!$$$         a(3)=(n3-1)*(1/a3)
!$$$         b(3)=1/a3+A(3)
!$$$         a(4)=(n4-1)*(1/a4)
!$$$         b(4)=1/a4+A(4)
!$$$
!$$$
!$$$
!$$$
!$$$      PRINT 9995, A(1),A(2),A(3),A(4)
!$$$      PRINT 9994, B(1),B(2),B(3),B(4)
!$$$ 9995 FORMAT ('A: ',5F8.5)
!$$$ 9994 FORMAT ('B: ',5F8.5)
!$$$
!$$$      CALL DCUHRE(NDIM, NF, A, B, MINCLS, MAXCLS, FTEST, ABSREQ, RELREQ,
!$$$     * KEY, NW, 0, FINEST, ABSEST, NEVAL, IFAIL, WRKSTR)
!$$$
!$$$      PRINT 9999, NW,MAXCLS,NEVAL, IFAIL
!$$$ 9999 FORMAT ('  NW = ',I10,' ,  MAXCLS = ',I12,/
!$$$     * '     FTEST CALLS = ', I12,
!$$$     * ', IFAIL = ', I2, /'        ESTIMATED ERROR    INTEGRAL')
!$$$ 9998    FORMAT (5X,  2F15.9)
!$$$      ANSWER=ANSWER+FINEST(1)
!$$$      ERROR=ERROR+ABSEST(1)
!$$$      PRINT 9996, ERROR, ANSWER
!$$$ 9996 FORMAT(8X, 'ERROR = ',F15.12, ' ,ANSWER = ',F15.12)
!$$$   10 CONTINUE
!$$$      END
!$$$      SUBROUTINE FTEST(NDIM, X, NFUN, F)
!$$$      implicit real*8(a-h,o-z)
!$$$      INTEGER N, NDIM, NFUN
!$$$      DOUBLE PRECISION X(NDIM), F(NFUN), FINT
!$$$
!$$$      pi=4*atan(1.0)
!$$$c     change of variable
!$$$      a12=pi*X(1)/2
!$$$      a13=pi/2-a12+a12*X(2)
!$$$      t2=pi-2*a13+(2*a12-pi+2*a13)*X(3) 
!$$$      t3= pi -t2+ (2*a13+t2-pi)*X(4)
!$$$c     jacobian
!$$$      fjac=(pi/2)*a12*(2*a12-pi+2*a13)*(2*a13+t2-pi)
!$$$     
!$$$      f1=64*sin(t2+t3)**6*cos(a12)**3*cos(a13)**3*sin(t2)*sin(t3)
!$$$      f2=f1*cos(a12-t2)*cos(a13-t3)
!$$$      f3 = cos(a12)**2*sin(t2+t3)**2*(t3+cos(2*a13-t3)*sin(t3))
!$$$      f4 = cos(a13)**2*sin(t2+t3)**2*(t2+cos(2*a12-t2)*sin(t2))
!$$$      f5 = cos(t3-a13)**2*cos(a12)**2+cos(t2-a12)**2*cos(a13)**2
!$$$      f6 = 2*cos(t2-a12)*cos(t3-a13)*cos(a12)*cos(a13)*cos(t2+t3)
!$$$      f7 = (2*pi-(t2+t3))*(f5-f6)
!$$$      f8 = -2*sin(t2+t3)*cos(a13)*cos(t2-a12)
!$$$      f9 = cos(t3-a13)*cos(a12)-cos(t2+t3)*cos(t2-a12)*cos(a13)
!$$$      f10 = -2*sin(t2+t3)*cos(a12)*cos(t3-a13)
!$$$      f11 = cos(t2-a12)*cos(a13)-cos(t2+t3)*cos(t3-a13)*cos(a12)
!$$$      f12 = 0.5*(f8*f9+f10*f11)
!$$$      f13 = (f3+f4+f7+f12)
!$$$c     F(1) is the integrand
!$$$      F(1)=2*pi*f2*fjac/(f13**3)
!$$$       
!      end
!    SUBROUTINE DCUHRE(NDIM,NUMFUN,A,B,MINPTS,MAXPTS,FUNSUB,EPSABS,
!   +                  EPSREL,KEY,NW,RESTAR,RESULT,ABSERR,NEVAL,IFAIL,
!    +                  WORK)!
!


!!$!   DTEST1 is a simple test driver for DCUHRE.
!!$!
!!$!   Output produced on a SUN 3/50.
!!$!
!!$!       DCUHRE TEST RESULTS
!!$!
!!$!    FTEST CALLS = 3549, IFAIL =  0
!!$!   N   ESTIMATED ERROR    INTEGRAL
!!$!   1       0.000000       0.138508
!!$!   2       0.000000       0.063695
!!$!   3       0.000009       0.058618
!!$!   4       0.000000       0.054070
!!$!   5       0.000000       0.050056
!!$!   6       0.000000       0.046546
!!$!
!!$      PROGRAM DTEST1
!!$        use mykinds
!!$        use testfunc
!!$        use dcuhre90
!!$!        EXTERNAL FTEST
!!$        INTEGER KEY, N, NDIM, MINCLS, MAXCLS, IFAIL, NEVAL, NW
!!$        PARAMETER (NDIM = 5, NW = 5000, NF = 1)
!!$        real(dp) ::  A(NDIM), B(NDIM), WRKSTR(NW)
!!$        real(dp) ::  ABSEST(NF), FINEST(NF), ABSREQ, RELREQ
!!$        DO N = 1,NDIM
!!$           A(N) = 0
!!$           B(N) = 1
!!$        end DO
!!$!10      CONTINUE
!!$        MINCLS = 0
!!$        MAXCLS = 10000
!!$        KEY = 0
!!$        ABSREQ = 0
!!$        RELREQ = 1E-3
!!$        CALL DCUHRE(NDIM, A, B, MINCLS, MAXCLS, FTEST, ABSREQ, RELREQ,  &
!!$             & KEY, NW, 0, FINEST, ABSEST, NEVAL, IFAIL, WRKSTR)
!!$!      PRINT 9999, NEVAL, IFAIL
!!$! 9999 FORMAT (8X, 'DCUHRE TEST RESULTS', //'     FTEST CALLS = ', I4, &
!!$!     & ', IFAIL = ', I2, /'    N   ESTIMATED ERROR    INTEGRAL')
!!$        write(*,*), '       DCUHRE TEST RESULTS'
!!$        write(*,'(a,i4,a,i2)'), '     FTEST CALLS = ', NEVAL,', IFAIL = ', IFAIL
!!$        write(*,*), '   N   ESTIMATED ERROR    INTEGRAL'
!!$        
!!$        DO  N = 1,NF
!!$!         PRINT 9998, N, ABSEST(N), FINEST(N)
!!$! 9998    FORMAT (3X, I2, 2F15.6)
!!$         write(*,'(a,i2,2f15.6)'), '   ',N,ABSEST(N),FINEST(N)
!!$!20      CONTINUE
!!$      end DO
!!$    END PROGRAM DTEST1
