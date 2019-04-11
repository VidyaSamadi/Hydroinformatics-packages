      function funcval(param_vl,nopt)

      include 'DimArrays.inc'
      include 'comVar.inc'
      include 'comRead.inc'
      include 'comRainRun.inc'
      include 'comAutCal.inc'
 
      integer nvent,NRIGHE,TotOss,totpl,chkPL,m,chkTE,KK,ztop
      integer manP(Nraing),manT(Ntermo),STpn(NRow,NCol),STtn(NRow,NCol)
      integer p,npf,PasPio,I,J,K,L,LL,Nstep,nums,nopt,totStep
      real tot,contg,totale,sumsq,sumsq1,contQtot
      real PioN(NRow,NCol),INVpn(NRow*Nraing,NCol)
      real*8 param_vl(16)
      character*1 dum


      call Initial_value ! inizialization of global arrays

      A0=param_vl(1)
      Wv=param_vl(2)
      Ksv(1)=param_vl(3)
	Ksv(2)=param_vl(4)
	Ksv(3)=param_vl(5)
      Ksat=param_vl(6)
      Hs=param_vl(7)
      H=param_vl(8)
      cint=param_vl(9)

      sumsq=0.  ! objective function
      sumsq1=0.
	totStep=0
      KK=1  
      Nstep=NOss*(StepPio/DT)
      K=1        
      contg=0.           ! the evapotranspiration formula requires the updating
      m=mm               ! of the current simulation month
      totale=totg 
	  
      do while (K.LE.Nstep)
         PasPio=(DT/StepPio)*(K-1.)+1
         chkPL=0
         totpl=0
         do p=1,Nplu  ! checking the availability of rainfall data for the current time step
            if (FGPIO(PasPio,p).EQ.1) then
C              rainfall data is missing at raingauge p
               chkPL=1
               totpl=totpl+1
               manP(totpl)=p
            end if
         end do
         if (chkPL.EQ.1) then
            if (risp2.EQ.1) then
C              application of the Thiessen method by excluding the raingauges with 
C              missing data whose number is stored in manP()
               call NewThiessen(Nplu,Coordp,totpl,manP,STpn,Slope)
            else
C              application of the inverse distance method by excluding the raingauges 
C              with missing data whose number is stored in manPp()
               call NewInvDist(Nplu,Coordp,totpl,manP,INVpn,Slope)
            end if
         end if

C        initialization of the array of the surface runoff intensity
         do I=1,idim
            do J=1,jdim
               PioN(I,J)=0.
            end do
         end do

C        continuous rainfall-runoff simulation
C        checking the availability of the temperature data
         call Check_temp(PasPio,K,Nstep,chkTE,STtn,manT)
         call Continuous(PasPio,K,Nstep,m,chkPL,chkTE,
     .                   STtn,PioN,STpn,INVpn)
         call Update_day(contg,totale,m)

         open(14,FILE='ztop.in',IOSTAT=ios,STATUS='OLD')
         read(14,*) ztop
         close(14)


         call Routing(PioN)   ! routing procedure

         if (Nmis.GT.0) then
            do LL=1,Nmis
               Qsup(KK,LL)=Q(outI(LL),outJ(LL))
               Qtot(KK,LL)=Q(outI(LL),outJ(LL))+QI(outI(LL),outJ(LL))                       
            end do
          KK=KK+1
         end if

C        updating the arrays of the surface and sub surface flow
         do I=1,idim
            do J=1,jdim
               Y(I,J)=Q(I,J); YY(I,J)=QQ(I,J)
               Q(I,J)=0.; QQ(I,J)=0.
               YI(I,J)=QI(I,J); YYI(I,J)=QQI(I,J)
               QI(I,J)=0.; QQI(I,J)=0.
            end do
         end do
         K=K+1                
      end do           

      totStep=K-1

	if (StepPio.ne.DT) then
	   do LL=1,Nmis
            K=1
            KK=1
            do while (K.LT.totStep)   
               I=K
               contQtot=0
               do while (I.LE.(K+(StepPio/DT-1)))
                  contQtot=contQtot+Qtot(I,LL)
                  I=I+1
               end do
               Qtot1(KK,LL)=contQtot/(StepPio/DT)
               KK=KK+1
               K=I
            end do
         end do
         totStep=totStep/(StepPio/DT)
	end if    

	if (Nmis.GT.0) then
	   do KK=1,totStep
            do LL=1,Nmis
               if (StepPio.ne.DT) then
 	         sumsq=sumsq+(Qtot1(KK,LL)-ObsDisc(KK,LL))**2 
 	       else
 	         sumsq=sumsq+(Qtot(KK,LL)-ObsDisc(KK,LL))**2 
 	       end if
            end do
         end do 
	end if    


      sumsq1=sumsq
      sumsq=0
      funcval=(sumsq1)

      open(unit=127,file='min-ris.txt',status='unknown')
120   read(127,'(a1)',end=128) 
	go to 120
128   backspace(127)

      write(*,*) 'Parameter values: ',(param_vl(i),i=1,nopt)
      write(*,*) 'Objective function = ',sumsq1
      write(127,*) 'Parameter values: ',(param_vl(i),i=1,nopt)  
      write(127,*) 'Objective function = ',sumsq1
      close(127)

      open(unit=31,file='outlet.out',status='replace')
C      write(31,*) '     Step   Total Flow' 

	if (Nmis.GT.0) then
	   do KK=1,totStep
            do LL=1,Nmis
               if (StepPio.ne.DT) then
               write(31,*) kk,Qtot1(kk,LL)
 	       else
               write(31,*) kk,Qtot(kk,LL)
 	       end if
            end do
         end do 
	end if    

	close (31)
	   
      return
      end