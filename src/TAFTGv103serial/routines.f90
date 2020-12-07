   MODULE routines
   !
   IMPLICIT NONE
   !
   !
	CONTAINS
	 
	 
	

	
	

!   SUBROUTINE FAST(b, N, Ne, k, p0, distobs, NoS, S, t, treal, splan, c, FSTobs, FSTT, CountFST)
!	   
!	   !This subprogram performs a single simulation with fixed initial allele frequencies (p0). 
!	   
!       USE irandom
!       
!       
!       IMPLICIT NONE
!     
!	   !Variables
!	   !Externals
!       INTEGER, INTENT(IN):: b				         !Size of the array N, adjusted for continuous or discrete generations
!	   INTEGER, INTENT(IN), DIMENSION(b+1):: N       !Total population size 
!	   INTEGER, INTENT(IN), DIMENSION(b):: Ne	     !Effective population size
!       INTEGER, INTENT(IN):: k                       !Nr of alleles
!       INTEGER, INTENT(IN), DIMENSION(k):: p0        !Initial allele frequencies
!       REAL(8), INTENT(IN):: distobs  	             !Distance between allele frequencies in consecutive samples
!       INTEGER, INTENT(IN)::  NoS                    !Nr of samples
!       INTEGER, INTENT(IN), DIMENSION(NoS):: S       !Sampling sizes
!       INTEGER, INTENT(IN), DIMENSION(NoS):: t	     !Number of generations between sample 0 and each sample
!       REAL(8), INTENT(IN), DIMENSION(NoS):: treal   !Same as t for semi-continuous generations
!       INTEGER, INTENT(IN):: splan			         !Switch para usar el plan de muestreo I o II
!       REAL(8), INTENT(INOUT) :: c      	         !Counter of positive cases
!	   REAL(8), INTENT(IN):: FSTobs					 !Onserved value of Fst
!	   REAL(8), INTENT(INOUT):: FSTT 			     !Temporal Fst 
!	   INTEGER, INTENT(INOUT):: countFST			 !Counter of positive temporal Fst's (significance)
!	
!       !Internals
!       INTEGER, ALLOCATABLE, DIMENSION (:):: Ninterna   ! Total population size for internal use (can change on the go)   
!	   INTEGER, DIMENSION (NoS,k):: xy                  ! Nr. of alleles in each sample 
!	   INTEGER, ALLOCATABLE, DIMENSION (:,:):: p, q     ! Nr. of alleles in the total population and in the effective population 
!	   REAL(8), ALLOCATABLE, DIMENSION (:):: qr         ! Relative allele frequencies at effective population   
!	   REAL(8):: Distsim, FSTTlocal                     ! Distance among simulated samples; FSTT of the current simulation
!	   INTEGER:: i, l, m, w, Total                      ! Counters
!	   REAL(8):: HT, HS, Hx, Hy 
!	   !-----------------------------------------------------------------------------------------------
!	   !-----------------------------------------------------------------------------------------------
!	   ALLOCATE(p(b+1, k), q(b, k), qr(k-1), Ninterna(b+1))
!	   
!	   !Initializing the array P(:,:)  
!	   Ninterna(1) = 0
!       DO i=1, k, 1	   
!	      p(1,i) = p0(i)
!	  	  Ninterna(1) = Ninterna(1) + p(1,i)
!	   ENDDO
!       DO i=2, b+1, 1
!	      Ninterna(i)=N(i)
!       ENDDO
!       !--------------------------------
!	   !S I M U L A T I O N   S T A R T 
!       !--------------------------------
!	   !Call the hypergeometric sampler for sample 0		   
!       CALL multhyper(k,p0,S(1),xy(1,:))
!	       
!	   DO w = 1,  NoS-1, 1
!	      !Adjustments for sampling plan 1 (before reproduction):
!	      !Organisms are not reeplaced for obtaining the effective population
!	      IF(splan==2) THEN
!	         DO l=1, k, 1
!		        p(t(w)+1,l)=p(t(w)+1,l)-xy(w,l)			   
!		     ENDDO
!		     Ninterna(t(w)+1)=Ninterna(t(w)+1)-S(w)
!	      ENDIF
!	      !This cycle generates all intermediate generations before generation t
!	      DO m=t(w)+1, t(w+1), 1		!m=t+1!
!	         IF (m==t(w)+1) THEN
!                !Call the hypergeometric sampler for the effective population
!				CALL multhyper(k,p(m,:),Ne(m),q(m,:)) 
!             ENDIF
!		     !Call the binomial sampler for next generation whole population 
!	         IF (m<t(w+1)) THEN
!		        DO i=1,k-1,1
!		           qr(i)=REAL(q(m,i))/REAL(Ne(m))
!		        ENDDO
!		        CALL genmul(Ne(m+1), qr, k, q(m+1,:))
!			 ELSE
!		        DO i=1, k-1, 1
!		           qr(i)=REAL(q(m,i))/REAL(Ne(m))
!		        ENDDO
!		        CALL genmul(N(m+1), qr, k, p(m+1,:))
!             ENDIF
!          ENDDO   
!          !Call the hypergeometric sampler for sample t			   
!          CALL multhyper(k,p(t(w+1)+1,:),S(w+1),xy(w+1,:)) 
!       ENDDO
!       !WE TAKE THE SUM OF DISTANCES BETWEEN FREQUENCIES OF CONSECUTIVE SAMPLES
!	   Distsim=0.0       
!	   DO l = 1, NoS-1, 1 
!	      DO i=1, k-1, 1
!	         Distsim = ABS( ((REAL(xy(l,i)))/(REAL(S(l))))-((REAL(xy(l+1,i)))/(REAL(S(l+1)))) ) + Distsim
!	      ENDDO
!	   ENDDO
!	   IF (Distsim>=Distobs) THEN
!	      c=c+1
!	   ENDIF
!	   !-----------------------
!       !-----------------------
!	   !FINALLY WE CALCULATE FSTT FOR THIS SIMULATION
!	   FSTTlocal = 0
!       DO l = 1, NoS-1, 1
!	      HT=1
!	      DO i=1, k, 1
!	         HT=HT-(((xy(l,i)/REAL(S(l)))+(xy(l+1,i)/REAL(S(l+1))))/2)**2
!	      ENDDO
!	      Hx=1
!	      Hy=1
!	      DO i=1, k, 1
!	         Hx=Hx-((xy(l,i)/REAL(S(l)))**2)
!	         Hy=Hy-((xy(l+1,i)/REAL(S(l+1)))**2)
!	      ENDDO
!	      HS = (Hx+Hy)/2
!		  IF ( HT /= 0 ) THEN	              			    
!	         FSTTlocal = FSTTlocal+((HT-HS)/HT)
!	      ENDIF
!       ENDDO
!	   FSTTlocal = FSTTlocal/(NoS-1)
!       IF (FSTobs>FSTTlocal) THEN
!	      countFST=countFST+1
!	   ENDIF
!       FSTT=FSTT+FSTTlocal	       
!	 
!	   DEALLOCATE(p, q, qr, Ninterna)
!	 
!500    FORMAT(' ',F10.4,'%')
!       
!       
!       RETURN
!	   
!   END SUBROUTINE FAST


   

        
   SUBROUTINE Dirichlet_param(k,NoS,xy,le,Ne,S,t,D)
     
         IMPLICIT NONE
     
         INTEGER, INTENT(IN):: k                         !Number of alleles (for allocating the arrays)
         INTEGER, INTENT(IN):: NoS                       !# of samples
         REAL(8), DIMENSION(NoS,k), INTENT(IN):: xy      !Observed allele frequencies                  
         INTEGER, INTENT(IN):: le                        !=t(NoS)+1; lenght of Ne array 
         INTEGER, DIMENSION(le), INTENT(IN):: Ne         !Effective population sizes
         INTEGER, DIMENSION(NoS), INTENT(IN):: S         !Sample sizes 
         INTEGER, DIMENSION(NoS), INTENT(IN):: t         !Number of generations until samples 
         REAL(8), DIMENSION(k), INTENT(OUT):: D          !Vector of parameters adjusted 
     
         INTEGER:: i, j, f
         REAL(16):: phi
         !-----------------------
         !Original Dirichlet    
         DO i=1,k,1   !First we collect freq at sample 0 (they don't have adjust by phi)) 
	        D(i)=xy(1,i)
         ENDDO
         DO i=1,k,1       !Then every posterior sample 
	        DO j=2,NoS,1
		       phi = 1-(1/DBLE(S(j)))
		       DO f = 1, t(j), 1
		          phi = phi * (1-(1/DBLE(Ne(f)))) 
		       ENDDO
		       D(i) = D(i) + xy(j,i)*phi
		    ENDDO
         ENDDO
		 
     
         RETURN
     
   ENDSUBROUTINE Dirichlet_param
         
         
         
         
         
   REAL FUNCTION Var_Cov(i,j,l,o,NoS,k,tt,S,FrecEst,t,N,Ne,splan)

       IMPLICIT NONE
       
       !IMPORTED VARIABLES
       INTEGER, INTENT(IN):: i, l  !INDICATE THE NUMBER OF ALLELE
       INTEGER, INTENT(IN):: j, o  !INDICATE THE NUMBER OF SAMPLE
       INTEGER, INTENT(IN):: NoS, k, tt  
       INTEGER, DIMENSION(NoS), INTENT(IN):: S, t 
       REAL(8), DIMENSION(k), INTENT(IN):: FrecEst
       INTEGER, DIMENSION(tt), INTENT(IN)::  Ne
       INTEGER, DIMENSION(tt+1), INTENT(IN)::	N
       INTEGER, INTENT (IN):: splan
       !INTERNAL VARIABLES
       INTEGER:: a
       REAL(8):: phiN		                     !COUNTER FOR LOOPS
       !++++++++++++++++++++++++
       !    S   T   A   R   T
       !++++++++++++++++++++++++
       !THE VAR_COV ESTIMATION COMES IN THE NEXT OPTIONS:
       
       !1. COVARIANCES BETWEEN AN ALLELE VERSUS ITSELF----------------------------------------
       IF (i==l) THEN
          !AT SAME GENERATION: ITSELF (THE ELEMENTS IN THE DIAGONAL -VARIANCES OF THEM-)
          IF (j==o) THEN
	         !IF THAT GENERATION IS ZERO
	         IF (j==1) THEN
	            Var_Cov = FrecEst(i)*(1-FrecEst(i))/REAL(S(1))	      	                              !ECUATION (4) OF WAPLES ARTICLE
	         !IF THAT GENERATION IS NOT ZERO
	         ELSE
	            phiN=1
		        DO a=1, t(j) 
	               phiN = ( 1-(1/REAL(Ne(a))) )*phiN
	            ENDDO
		        Var_Cov = FrecEst(i)*(1-FrecEst(i)) * ( 1 - ((1-(1/REAL(S(j))))*phiN) )                 !ECUATION (5) OF WAPLES ARTICLE
	         ENDIF
             !AT DIFFERENT GENERATIONS	
	      ELSE  
	         !COVARIANCES BETWEEN AN ALLELE IN GENERATION ZERO AND ONE IN OTHER GENERATION 
	         IF ((j==1).OR.(o==1)) THEN
	            IF (splan==1) THEN
			       Var_Cov = FrecEst(i)*(1-FrecEst(i))/N(1)                                             !ECUATION (A1a) APPENDIX WAPLES ARTICLE
		        ELSEIF (splan==2) THEN
	               Var_Cov = 0   						                                                  !ECUATION (A2a) APPENDIX WAPLES ARTICLE
		        ENDIF
	            !BETWEEN ALLELES IN GENERATIONS DIFFERENT THAN ZERO
	         ELSE
	            phiN=1	  !phiN CALCULUS IS COMMON TO BOTH OPTIONS
		        DO a=1,t(MIN(j,o)),1 
	               phiN = ( 1-(1/REAL(Ne(a))) )*phiN
	            ENDDO
	            IF (splan==1) THEN
	               Var_Cov = FrecEst(i)*(1-FrecEst(i)) * ( 1 - ( phiN * (1-(1/REAL(N(MIN(o,j))))) ) )   !ECUATION (A1b) APP. WAPLES ART.
		        ELSEIF (splan==2) THEN
		           Var_Cov = FrecEst(i)*(1-FrecEst(i)) * ( 1 - phiN )                                   !ECUATION (A2b) APP. WAPLES ART. 
		        ENDIF	      
	         ENDIF
	      ENDIF
       ENDIF
       !2. COVARIANCES BETWEEN DIFFERENT ALLELES----------------------------------------------
       IF (i/=l) THEN
          !AT SAME GENERATION
	      IF (j==o) THEN
             !IF THAT GENERATION IS ZERO
	         IF (j==1) THEN	 !(NO MATTER THE SAMPLING PLAN: ECUATION IS THE SAME)
	   	        Var_Cov = -FrecEst(i)*FrecEst(l)/S(1)                                                   !ECUATION (A1c & A2c) APP. WAPLES ART. 
	            !IF THAT GENERATION IS NOT ZERO
	         ELSE
	            phiN=1	  !phiN CALCULUS 
		        DO a=1,t(j),1 					   !(NO MATTER THE SAMPLING PLAN: ECUATION IS THE SAME)
	               phiN = ( 1-(1/REAL(Ne(a))) )*phiN
	            ENDDO
	            Var_Cov = -FrecEst(i)*FrecEst(l) * ( 1 - ( phiN * (1-(1/REAL(S(j)))) ) )                !ECUATION (A1e & A2e)
	         ENDIF 	  
	      !AT DIFFERENT GENERATIONS
          ELSE
	         !BETWEEN AN ALLELE IN GENERATION ZERO AND ONE IN OTHER GENERATION
	         IF ((j==1).OR.(o==1)) THEN
	            IF (splan==1) THEN
	               Var_Cov = -FrecEst(i)*FrecEst(l) / N(1)                                              !ECUATION (A1d)
		        ELSEIF (splan==2) THEN 
		           Var_Cov = 0                                                                          !ECUATION (A2d)
		        ENDIF
	         !BETWEEN ALLELES IN GENERATIONS DIFFERENT THAN ZERO
	         ELSE
	            phiN=1	  !phiN CALCULUS IS COMMON TO BOTH OPTIONS
		        DO a=1, t(MIN(j,o)) 
	               phiN = ( 1-(1/REAL(Ne(a))) )*phiN
	            ENDDO
	            IF (splan==1) THEN
	               Var_Cov = -FrecEst(i)*FrecEst(l) * ( 1 - ( phiN * (1-(1/REAL(N(MIN(o,j))))) ) )      !ECUATION (A1f)
		        ELSEIF (splan==2) THEN 
		           Var_Cov = -FrecEst(i)*FrecEst(l) * ( 1 - phiN  )                                     !ECUATION (A2f)
		        ENDIF
	         ENDIF
	      ENDIF
       ENDIF
       
       
       RETURN

   END FUNCTION Var_Cov



   

   SUBROUTINE Frec_Estim (k, NoS, b, xyobs, S, t, N, Ne, splan, FrecEst)

        IMPLICIT NONE

		!EXTERNAL VARIABLES
		
		INTEGER, INTENT(IN):: k, NoS, b
		REAL(8), INTENT(IN), DIMENSION(NoS,k):: xyobs
		INTEGER, INTENT(IN), DIMENSION(NoS):: S, t
		INTEGER, INTENT(IN), DIMENSION(b+1):: N
		INTEGER, INTENT(IN), DIMENSION(b):: Ne
		INTEGER, INTENT(IN):: splan
		REAL(8), INTENT(OUT), DIMENSION(k):: FrecEst
		!INTERNAL VARIABLES
		INTEGER:: i, j, l, a
		REAL(8):: phiN, M, SigFrec 				 ! Suma de N y de Ne en caso de N/Ne variable en el tiempo
        REAL(8), ALLOCATABLE :: Sigma(:,:), SigmaInv(:,:)  
	    REAL(8), ALLOCATABLE :: Frec(:,:), SigVector(:), Vector(:)  
	    REAL(8)::det(2)	
	    INTEGER:: info
		!----------------------------------------------------------------------------------------------------
		! I   N   I   C   I   O
		!----------------------------------------------------------------------------------------------------
	    ALLOCATE(  Sigma(NoS,NoS), SigmaInv(NoS,NoS), Frec(NoS,k)  ) 
	    ALLOCATE(  SigVector((NoS*(NoS+1))/2)  , Vector(NoS)  )
		!
		!PLACING THE ALLELE FREQUENCIES IN AN ARRAY
	    DO i=1, NoS
	       DO j=1, k
	          Frec(i,j)=xyobs(i,j)/S(i)
	       ENDDO
	    ENDDO
		!
		! MAKING THE SIGMA MATRIX
        Sigma(1,1)=1/REAL(S(1))
	    DO i=2, NoS, 1   
	       !
	       DO j=1, i-1, 1
	          ! CALCULATING phiN: A TERM THAT IS PRODUCT OF TERMS (1-1/Ne) WITH DIFFERENT Ne'S 
		      phiN=1
		      IF (j>1) THEN
		         DO l=1,t(j),1 
	                phiN = ( 1-(1/REAL(Ne(l))) )*phiN
	             ENDDO
		      ENDIF
		      !
		      ! ASIGNING VALUES TO OFF-DIAGONAL PLACES IN THE SIGMA MATRIX (COVARIANCES)
	  	      IF (splan==1) THEN
		         Sigma(i,j)=1-( phiN*(1-(1/REAL(N(j)))) )
	   	      ELSEIF (splan==2) THEN
	             Sigma(i,j)=1-phiN
	          ENDIF
		      !
		      Sigma(j,i)=Sigma(i,j)
	       ENDDO
	       ! ASIGNING DIAGONAL ELEMENTS OF THE SIGMA MATRIX (VARIANCES)
	       phiN=1
	       DO j=1, t(i) 	
	          phiN = ( 1-(1/REAL(Ne(j))) )*phiN
		      !WRITE(*,*)'phiN=',phiN,t(i),( 1-(1/REAL(Ne(j))) )
	       ENDDO	  
	       !
	       Sigma(i,i)=1-((1-1/REAL(S(i)))*(phiN))
	       !
	    ENDDO
	    !
	    !
	    !NOW WE OBTAIN THE INVERSE
	    a = 0
        DO j = 1, NoS  !FIRST WE OBTAIN A VECTOR WITH UPPER TRIANGLE ELEMENTS 
           DO i = 1, j !ARRANGED AS LINEAR VECTOR
              a = a + 1
              SigVector(a) = Sigma(i,j)
           ENDDO
        ENDDO
	    !--------------------------------------
	    CALL sppfa(SigVector,NoS,info)
	    CALL sppdi(SigVector,NoS,det,01)
	    !--------------------------------------
	    !NOW WE RETURN THE INVERSE VECTOR TO THE NORMAL FORM
	    a = 0
        DO j = 1, NoS  !FIRST WE OBTAIN A VECTOR WITH UPPER TRIANGLE ELEMENTS 
           DO i = 1, j !ARRANGED AS LINEAR VECTOR
              a = a + 1
              SigmaInv(i,j) = SigVector(a)
		      SigmaInv(j,i) = SigmaInv(i,j) 
           ENDDO
        ENDDO
	    !
	    !-----------------------------------------------------
	    !PRODUCT BETWEEN UNITARY VECTOR AND SIGMA INVERSE MATRIX
	    !
	    DO i=1, NoS
	       Vector(i)=0
	       DO j=1, NoS
	          Vector(i) = SigmaInv(j,i) + Vector(i)
	       ENDDO
	    !WRITE(*,*) 'Vector', Vector(i)
	    ENDDO
	    !------------------------------------------------------------------------------------
	    !WE OBTAIN THE SUM OF ALL THE ELEMENTS OF SIGMA INVERSE (M)
	    !
		M=0
	    DO i=1, NoS, 1
	       DO j=1, NoS, 1
	          M = SigmaInv(i,j) + M
	       ENDDO
	    ENDDO
	    !WRITE(*,*)'M=', M
	    !NOW WE OBTAIN THE DOT PRODUCT BETWEEN ALLELE FREQUENCIES VECTOR AND VECTOR (UNITARY VECTOR*SIGMA INVERSE)
	    !
		FrecEst(k)=1
	    DO i=1, k-1
	       SigFrec=0
	       CALL SDOT_1 (NoS,Vector,1,Frec(:,i),1,SigFrec)
	       !WRITE(*,*) SigFrec/M
	       !
	       FrecEst(i) = SigFrec/M
		   FrecEst(k) = FrecEst(k) - FrecEst(i)	
	    ENDDO
	    !-----------------------------------------------------------------------------------

		DEALLOCATE( Sigma, SigmaInv, Frec, SigVector, Vector )
	    
		!SSAY GODBAY
		!
   END SUBROUTINE Frec_Estim



   
   
   
   
   
   
  
   SUBROUTINE SPPFA (AP, N, INFO)
!***BEGIN PROLOGUE  SPPFA
!***PURPOSE  Factor a real symmetric positive definite matrix stored in
!            packed form.
!***LIBRARY   SLATEC (LINPACK)
!***CATEGORY  D2B1B
!***TYPE      SINGLE PRECISION (SPPFA-S, DPPFA-D, CPPFA-C)
!***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX FACTORIZATION, PACKED,
!             POSITIVE DEFINITE
!***AUTHOR  Moler, C. B., (U. of New Mexico)
!***DESCRIPTION
!
!     SPPFA factors a real symmetric positive definite matrix
!     stored in packed form.
!
!     SPPFA is usually called by SPPCO, but it can be called
!     directly with a saving in time if  RCOND  is not needed.
!     (Time for SPPCO) = (1 + 18/N)*(Time for SPPFA) .
!
!     On Entry
!
!        AP      REAL (N*(N+1)/2)
!                the packed form of a symmetric matrix  A .  The
!                columns of the upper triangle are stored sequentially
!                in a one-dimensional array of length  N*(N+1)/2 .
!                See comments below for details.
!
!        N       INTEGER
!                the order of the matrix  A .
!
!     On Return
!
!        AP      an upper triangular matrix  R , stored in packed
!                form, so that  A = TRANS(R)*R .
!
!        INFO    INTEGER
!                = 0  for normal return.
!                = K  if the leading minor of order  K  is not
!                     positive definite.
!
!
!     Packed Storage
!
!          The following program segment will pack the upper
!          triangle of a symmetric matrix.
!
!                K = 0
!                DO 20 J = 1, N
!                   DO 10 I = 1, J
!                      K = K + 1
!                      AP(K) = A(I,J)
!             10    CONTINUE
!             20 CONTINUE
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  SDOT
!***REVISION HISTORY  (YYMMDD)
!   780814  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  SPPFA
      INTEGER N,INFO
      REAL(8) AP(*)
!
      REAL(8) SDOT,T
      REAL(8) S
      INTEGER J,JJ,JM1,K,KJ,KK
!***FIRST EXECUTABLE STATEMENT  SPPFA
         JJ = 0
         DO 30 J = 1, N
            INFO = J
            S = 0.0E0
            JM1 = J - 1
            KJ = JJ
            KK = 0
            IF (JM1 .LT. 1) GO TO 20
            DO 10 K = 1, JM1
               KJ = KJ + 1
			   CALL SDOT_1(K-1,AP(KK+1),1,AP(JJ+1),1,SDOT)
               T = AP(KJ) - SDOT
               KK = KK + K
               T = T/AP(KK)
               AP(KJ) = T
               S = S + T*T
   10       CONTINUE
   20       CONTINUE
            JJ = JJ + J
            S = AP(JJ) - S
            IF (S .LE. 0.0E0) GO TO 40
            AP(JJ) = SQRT(S)
   30    CONTINUE
         INFO = 0
   40 CONTINUE
      RETURN
   END SUBROUTINE SPPFA


 


   SUBROUTINE SAXPY(N,SA,SX,INCX,SY,INCY)
!*     .. Scalar Arguments ..
      REAL(8) SA
      INTEGER INCX,INCY,N
!*     ..
!*     .. Array Arguments ..
      REAL(8) SX(*),SY(*)
!*     ..
!*
!*  Purpose
!*  =======
!*
!*     SAXPY constant times a vector plus a vector.
!*     uses unrolled loop for increments equal to one.
!*     jack dongarra, linpack, 3/11/78.
!*     modified 12/3/93, array(1) declarations changed to array(*)
!*
!*
!*     .. Local Scalars ..
      INTEGER I,IX,IY,M,MP1
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC MOD
!*     ..
      IF (N.LE.0) RETURN
      IF (SA.EQ.0.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
!*
!*        code for unequal increments or equal increments
!*          not equal to 1
!*
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX + 1
      IF (INCY.LT.0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
          SY(IY) = SY(IY) + SA*SX(IX)
          IX = IX + INCX
          IY = IY + INCY
   10 CONTINUE
      RETURN
!*
!*        code for both increments equal to 1
!*
!*
!*        clean-up loop
!*
   20 M = MOD(N,4)
      IF (M.EQ.0) GO TO 40
      DO 30 I = 1,M
          SY(I) = SY(I) + SA*SX(I)
   30 CONTINUE
      IF (N.LT.4) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
          SY(I) = SY(I) + SA*SX(I)
          SY(I+1) = SY(I+1) + SA*SX(I+1)
          SY(I+2) = SY(I+2) + SA*SX(I+2)
          SY(I+3) = SY(I+3) + SA*SX(I+3)
   50 CONTINUE
      RETURN
   END SUBROUTINE SAXPY



   

   SUBROUTINE SSCAL( N, SA, SX, INCX )
!****************************************************************************
!*                                                                          *
!*   DATA PARALLEL BLAS based on MPL                                        *
!*                                                                          *
!*   Version 1.0   1/9-92 ,                                                 *
!*   For MasPar MP-1 computers                                              *
!*                                                                          *
!*   para//ab, University of Bergen, NORWAY                                 *
!*                                                                          *
!*   These programs must be called using F90 style array syntax.            *
!*   Note that the F77 style calling sequence has been retained             *
!*   in this version for compatibility reasons, be aware that               *
!*   parameters related to the array dimensions and shape therefore may     *
!*   be redundant and without any influence.                                *
!*   The calling sequence may be changed in a future version.               *
!*   Please report any BUGs, ideas for improvement or other                 *
!*   comments to                                                            *
!*                    adm@parallab.uib.no                                   *
!*                                                                          *
!*   Future versions may then reflect your suggestions.                     *
!*   The most current version of this software is available                 *
!*   from netlib@nac.no , send the message `send index from maspar'         *
!*                                                                          *
!*   REVISIONS:                                                             *
!*                                                                          *
!****************************************************************************
      implicit none
!*
!*     scales a vector by a constant.
!*     uses unrolled loops for increment equal to 1.
!*     jack dongarra, linpack, 3/11/78.
!*
!*     .. Scalar Arguments ..
      INTEGER           INCX, N
      REAL(8)              SA
!*     ..
!*     .. Array Arguments ..
      REAL(8) :: SX(*)
!*     ..
!*     .. Local Scalars ..
      INTEGER           IX, jx
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC         ABS
!*     ..
!*     .. Executable Statements ..
!*
      IF( N.LE.0 ) return
!*
      ix = abs(incx)
      IF( iX.EQ.1 ) then
!*
!*        code for increment equal to 1
!*
        sx(1:n) = sa * sx(1:n)
      else
!*
!*        code for increment not equal to 1
!*
        jx = ix * ( n-1 ) + 1
        sx(1:jx:ix) = sa * sx(1:jx:ix)
      endif
!*
      RETURN
!*
!*     End of SSCAL .
!*
   END SUBROUTINE SSCAL





   SUBROUTINE SDOT_1 (N,SX,INCX,SY,INCY,SDOT)
!*     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
!*     ..
!*     .. Array Arguments ..
      REAL(8) SDOT
      REAL(8) SX(*),SY(*)
!*     ..
!*
!*  Purpose
!*  =======
!*
!*     forms the dot product of two vectors.
!*     uses unrolled loops for increments equal to one.
!*     jack dongarra, linpack, 3/11/78.
!*     modified 12/3/93, array(1) declarations changed to array(*)
!*
!*

!*     .. Local Scalars ..
      REAL(8) STEMP
      INTEGER I,IX,IY,M,MP1
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC MOD
!*     ..
      STEMP = 0.0e0
      SDOT = 0.0e0
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
!*
!*        code for unequal increments or equal increments
!*          not equal to 1
!*
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX + 1
      IF (INCY.LT.0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
          STEMP = STEMP + SX(IX)*SY(IY)
          IX = IX + INCX
          IY = IY + INCY
   10 CONTINUE
      SDOT = STEMP
      RETURN
!*
!*        code for both increments equal to 1
!*
!*
!*        clean-up loop
!*
   20 M = MOD(N,5)
      IF (M.EQ.0) GO TO 40
      DO 30 I = 1,M
          STEMP = STEMP + SX(I)*SY(I)
   30 CONTINUE
      IF (N.LT.5) GO TO 60
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
          STEMP = STEMP + SX(I)*SY(I) + SX(I+1)*SY(I+1) + &
                  SX(I+2)*SY(I+2) + SX(I+3)*SY(I+3) + SX(I+4)*SY(I+4)
   50 CONTINUE
   60 SDOT = STEMP
      RETURN
   END SUBROUTINE SDOT_1



   

   subroutine sppdi(ap,n,det,job)
      integer n,job
      real(8) ap(1)
      real(8) det(2)
!
!     sppdi computes the determinant and inverse
!     of a real symmetric positive definite matrix
!     using the factors computed by sppco or sppfa .
!
!     on entry
!
!        ap      real (n*(n+1)/2)
!                the output from sppco or sppfa.
!
!        n       integer
!                the order of the matrix  a .
!
!        job     integer
!                = 11   both determinant and inverse.
!                = 01   inverse only.
!                = 10   determinant only.
!
!     on return
!
!        ap      the upper triangular half of the inverse .
!                the strict lower triangle is unaltered.
!
!        det     real(2)
!                determinant of original matrix if requested.
!                otherwise not referenced.
!                determinant = det(1) * 10.0**det(2)
!                with  1.0 .le. det(1) .lt. 10.0
!                or  det(1) .eq. 0.0 .
!
!     error condition
!
!        a division by zero will occur if the input factor contains
!        a zero on the diagonal and the inverse is requested.
!        it will not occur if the subroutines are called correctly
!        and if spoco or spofa has set info .eq. 0 .
!
!     linpack.  this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.
!
!     subroutines and functions
!
!     blas saxpy,sscal
!     fortran mod
!
!     internal variables
!
      real(8) t
      real(8) s
      integer i,ii,j,jj,jm1,j1,k,kj,kk,kp1,k1
!
!     compute determinant
!
      if (job/10 .eq. 0) go to 70
         det(1) = 1.0e0
         det(2) = 0.0e0
         s = 10.0e0
         ii = 0
         do 50 i = 1, n
            ii = ii + i
            det(1) = ap(ii)**2*det(1)
!        ...exit
            if (det(1) .eq. 0.0e0) go to 60
   10       if (det(1) .ge. 1.0e0) go to 20
               det(1) = s*det(1)
               det(2) = det(2) - 1.0e0
            go to 10
   20       continue
   30       if (det(1) .lt. s) go to 40
               det(1) = det(1)/s
               det(2) = det(2) + 1.0e0
            go to 30
   40       continue
   50    continue
   60    continue
   70 continue
!
!     compute inverse(r)
!
      if (mod(job,10) .eq. 0) go to 140
         kk = 0
         do 100 k = 1, n
            k1 = kk + 1
            kk = kk + k
            ap(kk) = 1.0e0/ap(kk)
            t = -ap(kk)
            call sscal(k-1,t,ap(k1),1)
            kp1 = k + 1
            j1 = kk + 1
            kj = kk + k
            if (n .lt. kp1) go to 90
            do 80 j = kp1, n
               t = ap(kj)
               ap(kj) = 0.0e0
               call saxpy(k,t,ap(k1),1,ap(j1),1)
               j1 = j1 + j
               kj = kj + j
   80       continue
   90       continue
  100    continue
!
!        form  inverse(r) * trans(inverse(r))
!
         jj = 0
         do 130 j = 1, n
            j1 = jj + 1
            jj = jj + j
            jm1 = j - 1
            k1 = 1
            kj = j1
            if (jm1 .lt. 1) go to 120
            do 110 k = 1, jm1
               t = ap(kj)
               call saxpy(k,t,ap(j1),1,ap(k1),1)
               k1 = k1 + k
               kj = kj + 1
  110       continue
  120       continue
            t = ap(jj)
            call sscal(j,t,ap(j1),1)
  130    continue
  140 continue
      return
   end subroutine sppdi


   
   

   
   
   
   
   
   SUBROUTINE CHI2NC(X, F, THETA, Y, IFAULT)

	USE irandom
!
!       ALGORITHM AS 275 APPL.STATIST. (1992), VOL.41, NO.2
!
!       Computes the noncentral chi-square distribution function
!       with positive real degrees of freedom f and nonnegative
!       noncentrality parameter theta
!       
! NOTE OF ME (EDSON): ORIGINALLY WAS USED THE FUNCTION ALNGAM OR lngamma 
! PRESENT IN AS245 INSTEAD THE LNGAMMA OF RANDOM
!
      REAL(8), INTENT(IN):: X, F, THETA
      REAL(8), INTENT(OUT):: Y
      INTEGER, INTENT(OUT):: IFAULT
!
      INTEGER ITRMAX
      LOGICAL FLAG
      REAL(16):: ERRMAX, ZERO, ONE, TWO, LAM, N, U, V, X2, F2, T, TERM, BOUND !, ALNGAM 
!
      !EXTERNAL ALNGAM
!
      DATA ERRMAX, ITRMAX / 1.0E-6, 50 /
      DATA ZERO, ONE, TWO / 0.0, 1.0, 2.0 /
!
      Y = X
      IFAULT = 2
      IF (F .LE. ZERO .OR. THETA .LT. ZERO) RETURN
      IFAULT = 3
      IF (X .LT. ZERO) RETURN
      IFAULT = 0
      IF (X .EQ. ZERO) RETURN
      LAM = THETA / TWO
!
!       Evaluate the first term
!
      N = ONE
      U = EXP(-LAM)
      V = U
      X2 = X / TWO
      F2 = F / TWO
      T = X2 ** F2 * EXP(-X2) / EXP(lngamma(F2 + ONE))
!
!       There is no need to test IFAULT si
!       already been checked
!
      TERM = V * T
      Y = TERM
!
!       Check if (f+2n) is greater than x
!
      FLAG = .FALSE.
   10 IF ((F + TWO * N - X) .LE. ZERO) GO TO 30
!
!       Find the error bound and check for convergence
!
      FLAG = .TRUE.
   20 BOUND = T * X / (F + TWO * N - X)
      IF (BOUND .GT. ERRMAX .AND. INT(N) .LE. ITRMAX) GO TO 30
      IF (BOUND .GT. ERRMAX) IFAULT = 1
      RETURN
!
!       Evaluate the next term of the expansion and then the
!       partial sum
!
   30 U = U * LAM / N
      V = V + U
      T = T * X / (F + TWO * N)
      TERM = V * T
      Y = Y + TERM
      N = N + ONE
      IF (FLAG) GO TO 20
      GO TO 10
!
   END SUBROUTINE CHI2NC


   
   
   
   
   
   


   END MODULE routines