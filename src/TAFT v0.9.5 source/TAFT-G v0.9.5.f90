!****************************************************************************
!
!  PROGRAM: TAFTGv095
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program TAFTGv095


       use routines
	   use irandom
	   use luxury
    
       implicit none
     
       !Variables
       character(30):: text 
       logical:: SimTest, ChiTest, WapTest	  !For choosing tests
       character(20):: headline 
       integer:: algor                        !1=Empirical Bayes algorithm; 2=Full Bayesian 
       integer:: ploidia                      !Ploidia state (1=haploid,2=diploid,etc.)
       integer:: splan                        !Sampling plan (I or II following Waples))
       integer:: NoS                          !Nr of temporal samples 
       integer, allocatable:: t(:)            !Nr of generations between consecutive temporal samples
	   real, allocatable:: treal(:)           !Nr of generations between continuous temporal samples
       integer, allocatable:: nk(:)           !Vector with the nr of alleles 
       integer:: k, maxk                      !Number of alleles, maximum number of alleles 
	   logical:: Nevar                        !True if Ne is constant in time, if false a vector is readen    
       integer:: b                            !Total number of generations needed for array N
       integer, allocatable:: N(:), Ne(:)     !Total and effectivep opulation sizes
       integer, allocatable:: S(:)            !Sample sizes
       integer:: nloci                        !Nr of loci
       integer, allocatable:: nxyobs(:,:,:)   !Multilocus array with the allele frequencies
       real(8), allocatable:: xyobs(:,:)      !Nr of alleles in each sample
       real(8):: sim                          !Nr of simulations
       integer:: seed1                        !Seed for random number generation  
       
       integer:: i, j, l, m, o, a             !Counters
       integer:: y, z                         !Multipurpose  
       real:: x
       real(8):: x1
       integer:: ierror
       !For the tests
       integer, allocatable:: p0(:)           !Absolute allele freqs at population zero
       real, allocatable:: p0rel(:)           !Relative allele freqs at population zero
       real(8):: distobs                      !Observed distance between allele frequencies in samples
	   real(8):: phi, D_total                 !Used to estimate initial allele frequencies
       real(16), allocatable:: D(:)           !Parameter of a Dirichlet distribution
       real(8):: HT, HS, FSTobs, Hx, Hy       !Statistics used to calculate Fst
	   real:: FSTT, PositiveFSTT			  !Temporal Fst, count of positive cases of Fstt
	   integer:: countFST                     !Counter of Fstt's larger than observed Fst 
       real(8):: counter, c, cbis	          !Counters of simulations 
       integer:: adjust	                      !Adjustments   
	   logical:: ALLELE01, CONTIN01           !For choosing allele-by-allele, and semi-continuous generations 
       real:: Sump0rel                        !Sum of allele frequencies
       real(8):: proba                        !Final probability of the simulation test
       real:: Chi, ChiWap                     !Chi square, and Waples' adjusted Chi square
       real:: P_Chi, P_ChiWap				  !P-values of Chi sq. and Waples tests
       real, allocatable :: Totobs(:)     	  !Sum of frequencies for the contingency table
	   real:: TOT								
       real, allocatable:: SigVector(:), Vector(:), Vector2(:), Frec(:,:), FrecEst(:), P(:)
       real, allocatable:: estp(:), Sigma(:,:), SigmaInv(:,:)
       real, allocatable:: SIGMAc(:,:), SIGMAvector(:), SIGMAcInv(:,:)
       integer:: info, ifault
       real::det(2)	
       logical:: found
       integer, dimension(1):: xloc
	   !----------------------------------------------------------------------------------------------
	   ! I N I C I O
       !**********************************************************************************************
       !**********************************************************************************************
	   !(1) READING INPUT DATA
       !**********************************************************************************************
       !**********************************************************************************************
       Proba=0
       FSTT=0
       PositiveFSTT=0
       P_ChiWap=0
       P_Chi=0
       !++
       OPEN(UNIT=1, FILE='info.dat', STATUS='OLD', ACTION='READ', IOSTAT=ierror) ! Opening the file	    
       openif: if (ierror==0) then  
          !General data----------- 
          read(1,*,IOSTAT=ierror) text, SimTest
          if (SimTest) then
             read(1,*,IOSTAT=ierror) text, algor    
          endif    
          read(1,*,IOSTAT=ierror) text, WapTest
          read(1,*,IOSTAT=ierror) text, ChiTest 
          
          read(1,*,IOSTAT=ierror) text, ploidia
          read(1,*,IOSTAT=ierror) text, splan 
          !Samples----------------
          read(1,*,IOSTAT=ierror) text, NoS 
          !++
          allocate(t(NoS),treal(NoS),S(NoS)) 
          !++
          treal(1)=0
		  t(1)=0	   
		  do i=1,NoS-1,1
		     read(1,*,IOSTAT=ierror) text, t(i+1)
		     treal(i+1)=t(i+1)
          enddo
          !N/Ne-------------------
          b=t(NoS) 	         !b=t or b=t+1 if semi-continuous generations          
          allocate (N(b+1), Ne(b))
          read(1,*,IOSTAT=ierror) text, Nevar  
          if (Nevar) then
             read(1,*,IOSTAT=ierror) text, (N(i), i=1,b+1)  
             read(1,*,IOSTAT=ierror) text, (Ne(i), i=1,b) 
          else
             read(1,*,IOSTAT=ierror) text, z 
             N=z
             read(1,*,IOSTAT=ierror) text, z 
             Ne=z
          endif   
          !Adjusting for ploidia
	      do i=1,b+1,1
	         N(i)=ploidia*N(i) 
	      enddo
	      do i=1,b,1
	         Ne(i)=ploidia*Ne(i) 
          enddo 
          !Simulations and seed
          read(1,*,IOSTAT=ierror) text, sim 
          read(1,*,IOSTAT=ierror) text, seed1 
          !Number of loci--------------------
          read(1,*,IOSTAT=ierror) text, nloci
          !Maximum number of alleles---------
          read(1,*,IOSTAT=ierror) text, maxk
          
          allocate(nk(nloci),nxyobs(nloci,NoS,maxk))
          nxyobs=0
          do i=1,nloci,1
             !Number of alleles in this locus
             read(1,*,IOSTAT=ierror) text, nk(i)
             !Allele frequencies-----
             do j=1,NoS,1
                read(1,*,IOSTAT=ierror) text, (nxyobs(i,j,o), o=1,nk(i))
             enddo
          enddo   
          !++++++++++++++++++++++
          if (ierror/=0) then
             write(*,*) 'An error during input data file reading occurred'
             STOP
          endif
       else
          write(*,*) 'An error during input data file opening occurred'
          STOP
       endif openif         
       CLOSE(UNIT=1)
       
       !Checking for errors------------
       if ((.NOT.SimTest).and.(.NOT.ChiTest).and.(.NOT.WapTest)) then
          write(*,*)'Error: At least one analysis needs to be chosen' 
          STOP   
       endif	
       if (ploidia<1) then
          write(*,*)'Error: Ploidia state needs to be and integer > 0' 
          STOP   
       endif	
       do i=2,NoS,1
	      if ((t(i)-t(i-1))<1.0) then
             write(*,*)'ERROR!: Number of generations between samples cannot be lower than one' 
             STOP
	      endif
	   enddo
       if (ANY(nk<2)) then
           write(*,*)'ERROR!: Number of alleles must be 2 or more' 
           STOP    
       endif 
       if ((SimTest).or.(WapTest)) then
          do j=1, NoS, 1
		     if ((N(i)<=Ne(i)).OR.(N(i)<5).OR.(Ne(i)<5)) then 
	            write(*,*)'ERROR!: N, Ne cannot take those values'	
                STOP
             endif
          enddo
       endif
       do l=1,nloci,1
	      do i=1,NoS,1
	         do j=1,nk(l),1
	            if (nxyobs(l,i,j)<0) THEN
	               write(*,*)'ERROR!: Allele frequencies must be non-negative integers'
                   STOP
	            endif
             enddo
          enddo   
       enddo
       do l=1,nloci,1
	      do i=1,NoS,1
	         if (ALL(nxyobs(1,i,:)==0)) then
		        write(*,*)'ERROR!: Allele frequencies must be non-zero in at least one allele'
                STOP
             endif
          enddo   
       enddo
       do l=1,nloci,1
          do i=1,k,1
	         if (ALL(nxyobs(1,:,i)==0)) then
		        write(*,*)'ERROR!: Allele frequencies must be non-zero in at least one sample'
                STOP
             endif
          enddo
       enddo
       if ((algor==2).and.(sim<10)) then
	      write(*,*) 'ERROR!: Not enough number of simulations'   
          STOP
       endif
       !**********************************************************************************************
       !**********************************************************************************************
       !(2) ANALYSIS   
       !**********************************************************************************************
       !**********************************************************************************************
       CALL RLUXGO(4,seed1,0,0) 
       do l=1,nloci,1
           write (text,"(I15)") l 
           text=adjustl(text)
           write(*,*)'Analysing LOCUS ', text
           !----
           k=nk(l)
           allocate(xyobs(NoS,k),p0(k),p0rel(k))  
           do i=1,NoS,1  
              do j=1,k,1
                 xyobs(i,j)=nxyobs(l,i,j)
              enddo
           enddo
           S=0
           do i=1,NoS,1
              do j=1,k,1 
                 S(i)=S(i)+xyobs(i,j)
              enddo   
           enddo   
	       !(2.1) SIMULATION TEST
           !OBTAINING THE DISTANCE BETWEEN ALLELE FREQUENCIES OF SAMPLES
           Distobs=0
	       do i=1,NoS-1,1
	          do j=1,k-1,1
	             Distobs = ABS(((xyobs(i,j))/(REAL(S(i))))-((xyobs(i+1,j))/(REAL(S(i+1))))) + Distobs
	 	      enddo
           enddo
           !CALCULATING OBSERVED FST 
           if (SimTest) then 
	            FSTobs = 0 
	            do o=1,NoS-1,1
	               HT=1
	               do i=1,k,1
	                  HT=HT-(((xyobs(o,i)/REAL(S(o)))+(xyobs(o+1,i)/REAL(S(o+1))))/2)**2
	               enddo
	               Hx=1
	               Hy=1
	               do i=1,k,1
	                  Hx=Hx-((xyobs(o,i)/REAL(S(o)))**2)
	                  Hy=Hy-((xyobs(o+1,i)/REAL(S(o+1)))**2)
	               enddo
	               HS=(Hx+Hy)/2
		           if (HT/=0) then  	 
	                  FSTobs=((HT-HS)/HT)+FSTobs
		           endif 
                enddo
                FSTobs=FSTobs/(NoS-1)
	            FSTT=0
	            allocate(D(k))
	            !We calculate the parameter D 
                CALL Dirichlet_param(k,NoS,xyobs,b,Ne,S,t,D)
                D_total=SUM(D)
	            !SIMULATIONS ARE PERFORMED BY TWO WAYS: 
	            !1. EMPIRICAL BAYES: USING THE WAPLES ALLELE FREQUENCIES ESTIMATION AT GENERATION ZERO 
	            !2. FULL BAYESIAN: SAMPLE p0 FROM A DIRICHLET DISTRIBUTION CONDITIONAL TO THE SAMPLES 
                counter=0.0     ! COUNTER OF POSITIVE CASES (DISTANCE)
	            c=0.0		   ! COUNTER OF TOTAL NUMBER OF SIMULATIONS
	            countFST=0	   ! COUNTER OF POSITIVE CASES (FSTT)       
	            if (algor==1) then      
                    !EMPIRICAL BAYES ALGORITHM++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
	                !WE ESTIMATE ALLELE FREQUENCIES USING WAPLES METHOD
	                CALL Frec_Estim (k, NoS, b, xyobs, S, t, N, Ne, splan, p0rel)
                    !SINCE ALLELE FREQUENCIES CAN'T BE ZERO, UNOBSERVED ALLELES ARE ASSDIGNED A FREQUENCY = 1/2 OF THE MINIMUM OBSERVABLE    
                    adjust=0	
	                do i=1,k,1           !FIRST WE FIND HOW MANY ALLELES ARE ZERO
		               p0(i)=NINT(p0rel(i)*REAL(N(1)))
		               if (p0(i)<NINT((0.5*REAL(N(1)))/REAL(INT(D_total)))) then
		                  adjust=adjust + 1
		               endif			 
                    enddo 
	                if (adjust>0) then   ! THIS CONDITIONAL IS TRUE ONLY IF AT LEAST ONE FREQUENCY WAS ZERO
                       x1=MAX(NINT((0.5*REAL(N(1)))/REAL(INT(D_total))),1)  
	                   where (p0<NINT((0.5*REAL(N(1)))/REAL(INT(D_total))))
	                      p0=x1                          
                       elsewhere                          
		                  p0=MAX(1,p0-NINT(REAL(adjust*x1)/REAL(k-adjust)))
		               endwhere
	                endif
	                !FINALLY WE CALL THE SIMULATION PROGRAM--------------
	                CALL FAST(b,N,Ne,k,p0,distobs,NoS,S,t,treal,splan,sim,counter,c,CONTIN01,FSTobs,FSTT,countFST) 
                    !----------------------------------------------------
	                !EMPIRICAL BAYES ALGORITHM++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
                elseif (algor==2) then  
                    !FULL BAYESIAN ALGORITHM++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	                !FIRST WE GENERATE A VECTOR OF ALLELE FREQUENCIES AT GENERATION ZERO
                    do
                   
                       do   !This allows to keep sampling p0 until we all alleles get a non-null frequency 
                          CALL Dirichlet(k,real(D),p0rel)   
                          do i=1,k,1                   
		                     p0(i) = MAX( 1 , NINT(p0rel(i)*REAL(N(1))) )
                          enddo
                          adjust = SUM(p0)-N(1)
                          xloc=MAXLOC(p0)
                          p0(xloc(1))=p0(xloc(1))-adjust    
                          if (ALL(p0>0)) exit
                       enddo
                       Cbis=0
	                   !FINALLY WE CALL THE SIMULATION PROGRAM--------------
                       x1=1.0
		               CALL FAST(b,N,Ne,k,p0,distobs,NoS,S,t,treal,splan,x1,counter,cbis,CONTIN01,FSTobs,FSTT,countFST) 
	                   !----------------------------------------------------
                       if (MOD(INT(c),INT(sim/100))==0) then
                          WRITE(*,100) NINT((REAL(c)*100)/Sim)    
                          100 FORMAT(1x,I4,'%')                  
		               endif
		               c=c+cbis
		               if (c==sim) exit
	                enddo
	                !FULL BAYESIAN ALGORITHM++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	            endif
	            !FSTT------------------------------------------------------------
	            if (C>0) then
	               FSTT = FSTobs-(FSTT/REAL(C))
	               PositiveFSTT = (REAL(CountFST)/REAL(C))
	            else
	               FSTT = 0
		           PositiveFSTT = 0
	            endif
	            !----------------------------------------------------------------
	            !CALCULATING PROPORTIONS OBTAINED
	            if ((C>1).AND.(Counter>1)) then 
	               Proba = Counter/C
	            else
	               Proba = 0.0
                endif
	            !SIMULATION TEST 	 
	            !------------------------------------------------------------------------------------------
           endif
           if (ChiTest) then
	            !------------------------------------------------------------------------------------------
	            !CHI TEST 
	            ! 
	            Chi = 0.0
                ChiWap = 0.0
	            P_Chi = 0.0
	            P_ChiWap = 0.0
	            allocate(Totobs(k))
		        TOT=0
		        do i=1,k,1
		           Totobs(i)=0
		           do j=1,NoS,1
	                  Totobs(i)=xyobs(j,i)+Totobs(i)
			          TOT=TOT+xyobs(j,i)
	               enddo
                enddo
                do i=1,k,1
		           do j=1,NoS,1
	                  Chi=Chi+((((REAL(S(j))*Totobs(i))/TOT)- xyobs(j,i))**2)/((REAL(S(j))*Totobs(i))/TOT) 	         	   
	               enddo
                enddo
                CALL CHI2NC(Chi,REAL(k-1)*REAL(NoS-1),0.0,x,ifault)
	            P_Chi=1-x
                !------------------------------------------------------------------------------------------
           endif
           if (WapTest) then
	            !------------------------------------------------------------------------------------------
                !WAPLES TEST 
	            !
	            allocate(estp(k-1),Frec(NoS,k),FrecEst(k),Sigma(NoS,NoS),SigmaInv(NoS,NoS),Vector2(NoS*(k-1)))
	            allocate(SigVector((NoS*(NoS+1))/2),Vector(NoS),P(NoS*(k-1)),SIGMAc(NoS*(k-1),NoS*(k-1)))
	            allocate(SIGMAvector((NoS*(k-1))*((NoS*(k-1))+1)/2),SIGMAcInv(NoS*(k-1),NoS*(k-1)))
		        !PLACING THE ALLELE FREQUENCIES IN AN ARRAY
                do i=1,NoS,1
	               do j=1,k,1
	                  Frec(i,j)=xyobs(i,j)/S(i)
	               enddo
	            enddo
	            CALL Frec_Estim(k, NoS, b, xyobs, S, t, N, Ne, splan, FrecEst)
                !VECTOR OF DIFERENCES BETWEEN EXPECTED AND OBSERVED ALLELE FREQUENCIES 
	            do i=1,k-1,1
	               do j=1,NoS,1 
	                  P((NoS*(i-1))+j)=Frec(j,i)-FrecEst(i)
	               enddo
                enddo
                !CREATING THE SIGMA MATRIX (THE BIG ONE) BY PLACING ELEMENTS IN A '4D' ARRAY: 
	            !IT IS A 2D ARRAY OF MATRICES (EACH ONE ALSO a 2D)
                do m=1,k-1,1
	               do o=1,NoS,1
	                  do i=1,k-1,1  
		                 do j=1,NoS,1
			                SIGMAc((NoS*(i-1))+j,NoS*(m-1)+o)=Var_Cov(i,j,m,o,NoS,k,t(NoS),S,FrecEst,t,N,Ne,splan)
                         enddo						  
	 	              enddo
	               enddo											 
	            enddo
	            !NOW WE OBTAIN THE INVERSE
	            a=0			 !FIRST WE OBTAIN A VECTOR WITH UPPER TRIANGLE ELEMENTS
                do j=1,NoS*(k-1),1	 !ARRANGED AS A LINEAR VECTOR
                   do i=1,j,1
                      a=a+1
                      SIGMAvector(a)=SIGMAc(i,j)
                   enddo
                enddo     
                CALL sppfa(SIGMAvector,NoS*(k-1),info)
                CALL sppdi(SIGMAvector,NoS*(k-1),det,01)
	            a=0		 !WE RETURN THE VECTOR TO A MATRIX FORMAT
                do j=1,NoS*(k-1),1
                   do i=1,j,1
                      a=a+1
                      SIGMAcInv(i,j)=SIGMAvector(a) 
	 	              SIGMAcInv(j,i)=SIGMAcInv(i,j)
                   enddo
                enddo
	            !WE CALCULATE CHI BY MULTIPLYING THE P VECTOR BY SIGMAcINVERSE (RESULT IS Vector2) AND BY P AGAIN
	            do i=1,NoS*(k-1),1
	               CALL SDOT_1(NoS*(k-1),P,1,SIGMAcInv(:,i),1,Vector2(i)) 	     
	            enddo
	            CALL SDOT_1(NoS*(k-1),Vector2,1,P,1,ChiWap)
                CALL CHI2NC(ChiWap,REAL(k-1)*REAL(NoS-1),0.0,x,ifault)
	            P_ChiWap=1-x 
           endif
           !**********************************************************************************************
	       !  R E S U L T S 
           !**********************************************************************************************
           write (headline,"(I15)") l 
           headline=adjustl(headline)
           headline='Locus_'//TRIM(headline)
           headline=adjustl(headline)
           
           inquire(FILE='taft_out.csv', EXIST=found) 
           if (found) then
              if (l==1) then   !We replace (if found) only the first time; 
                 OPEN(UNIT=1, FILE='taft_out.csv', STATUS='REPLACE', ACTION='WRITE', IOSTAT=ierror) 
		      else
	             OPEN(UNIT=1, FILE='taft_out.csv', STATUS='OLD', ACTION='WRITE', ACCESS='APPEND', IOSTAT=ierror)  
              endif	   
	       else
              OPEN(UNIT=1, FILE='taft_out.csv', STATUS='NEW', ACTION='WRITE', IOSTAT=ierror) 			   
           endif 
           openif0: if (ierror==0) then
               if (l==1) then
                  write(1,*,iostat=ierror) 'Locus,Simulation_P,Fstt,Positive_Fstt,Waples_Test_P,Chi_Test_P'        
               endif
               write(1,1050,iostat=ierror) headline, Proba, FSTT, PositiveFSTT, P_ChiWap, P_Chi
               1050 format(1x,T1,A,',',5(F16.8,','))
               if (ierror/=0) then
                  write(*,*)'There was an error writting the results file'    
               endif    
           else
	          write(*,*)'There was an error in opening the results file'
           endif openif0
           CLOSE(UNIT=1)
	       !***********************************************************	
           !WE DEALLOCATE LOCUS-DEPENDENT ARRAYS
           if (SimTest) then
              deallocate(D)
           endif 
           if (ChiTest) then
              deallocate(Totobs)
           endif 
           if (WapTest) then
              deallocate(estp,Frec,FrecEst,Sigma,SigmaInv,Vector2,SigVector,Vector,P,SIGMAc,SIGMAvector,SIGMAcInv)
           endif 
           deallocate(xyobs,p0,p0rel)
       enddo
       !End the cycle of loci++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	   !WE DEALLOCATE THE REMAINING ARRAYS
       deallocate(t,treal,S,N,Ne,nk,nxyobs)
     
     


    end program TAFTGv095