      subroutine Read_param(chkread,chkmet,npf)
		
      implicit none

	include 'DimArrays.inc'
      include 'comRead.inc'
	include 'comRainRun.inc'
      	
	integer I,J,ios,nvent,npf,NRighe,nums,TotOss
	integer chkmet,chkread,chka,NHRou
      real tot
      character*30 NameFl
	character*3 dum


      open(11,FILE='Affdef.in',IOSTAT=ios,STATUS='OLD')
	if (ios.EQ.0) then
         do I=1,9
            read(11,*)
         end do
         read(11,*) gg,mm,aa
         do I=1,3
            read(11,*)
         end do
	   read(11,*) NOss      ! number of observed rainfall data to be read (= length of the simulation)
	   read(11,*)
	   if (NOss.LE.NpOss) then
            read(11,*) StepPio   ! time step of the rainfall data [s]
	      read(11,*)
            read(11,*) idim,jdim ! dimension of the DEM
	      if ((idim.LE.NRow).AND.(jdim.LE.NCol)) then  
	         read(11,*)
               read(11,*) DX,DY     ! size of the cell [m]
	         read(11,*)
               read(11,*) DT        ! time step of the simulation [s]
	         read(11,*)
               read(11,*) A0        ! constant critical support area [km^2]
	         read(11,*)
               read(11,*) Wv        ! parameter for the hillslope (width/level of water)
	         read(11,*) 
	         read(11,*) NHRou
               if (NHRou.LE.NStrcl) then
	            read(11,*)
	            read(11,*) (Ksv(I),I=1,NHRou) ! Strickler roughness [m^1/3 s^-1] for each class
                  read(11,*)
	            read(11,*) NameFl ! file of the distribution of roughness classes for the hillslope 
	            open(UNIT=22,FILE=NameFl,IOSTAT=ios,STATUS='old')
	            if (ios.eq.0) then
	               do I=1,idim
                        read(22,*) (KsvDistr(I,J),J=1,jdim)
	               end do
	               close(22)
                     read(11,*)
                     read(11,*) Wr,Ksr0,Ksr1 ! parameter for the river network (width/level of water and 
	               read(11,*)	           ! maximum and minimum Strickler roughnesses [m^1/3 s^-1]
                     read(11,*) Ksat         ! saturated hydraulic conductivity [m/s]
	               read(11,*)
                     read(11,*) Bp        ! width of the rectangular cross section of the sub surface flow [m]
	               read(11,*)
                     read(11,*) Hs        ! parameter of the infiltration reservoir [s]
	               read(11,*)
                     read(11,*) H         ! parameter of the infiltration reservoir
	               read(11,*) 
                     read(11,*) cint      ! parameter for the interception reservoir
	               read(11,*)
                     read(11,*) Ichi,Jchi ! coordinates of the catchment outlet
C                    value of risp1
C                    0 - single event simulation with normal catchment
C                    1 - single event simulation with impervious catchment
                     read(11,*)
                     read(11,*) risp1
C                    value of risp2
C                    1 - spatial interpolation of the rainfall with the Thiessen method
C                    2 - spatial interpolation of the rainfall with the inverse distance method
                     read(11,*)
                     read(11,*) risp2
	               read(11,*)
                     read(11,'(a30)') NameFl ! file of the slope pointer

                     open(UNIT=22,FILE=NameFl,IOSTAT=ios,STATUS='OLD')
	               if (ios.eq.0) then 
                        do I=1,idim
                           read(22,*) (PUN(I,J),J=1,jdim)
                        end do
                        close(22)
1000	                  read(11,*) dum
                        if (dum.ne.'A51') then
	                     goto 1000
                        else
                           read(11,'(a30)') NameFl  ! file of the modified DEM
                           open(UNIT=22,FILE=NameFl,IOSTAT=ios,
     .      				      STATUS='OLD')
	                     if (ios.eq.0) then
                              do I=1,idim
                                 read(22,*) (QUO(I,J),J=1,jdim)
                              end do
                              close(22)
	                        read(11,*)
                              read(11,'(a30)') NameFl  ! file of the contributing area
                              open(UNIT=22,FILE=NameFl,IOSTAT=ios,
     .           				    		  STATUS='OLD')
	                        if (ios.eq.0) then
                                 do I=1,idim
                                    read(22,*) (AREE(I,J),J=1,jdim)
                                 end do
                                 close(22)
	                           read(11,*)
                                 read(11,'(a30)') NameFl ! file of the link property
                                 open(UNIT=22,FILE=NameFl,IOSTAT=ios,
     .						               STATUS='OLD')
                                 if (ios.eq.0) then
                                    nvent=22
                                    read(nvent,*)
                                    NRighe=1
300                                 NRighe=NRighe+1
                                    read(nvent,*,end=400)
                                    go to 300
400                                 rewind nvent
                                    NRighe=NRighe-1

C                                   computation of the slope of the links
                                    call Slope_comp(nvent,NRighe)
                                    close(22)

C                                   computation of the average slope of the catchment
                                    tot=0
                                    nums=0
                                    do I=1,idim
                                       do J=1,jdim
                                          if (Slope(I,J).GT.0) then
                                             nums=nums+1
                                             tot=tot+Slope(I,J)
                                          end if
                                       end do
                                    end do
                                    write(*,*)
                                    write(*,*) 'Average slope of the 
     .catchment',tot/nums
                                    write(*,*)
                                    read(11,*) 
                                    if ((risp.EQ.2).AND.(risp1.EQ.1)) 
     .                                                          then
C                                      single event simulation with impervious catchment
                                       do I=1,idim
                                          do J=1,jdim
                                             S(I,J)=0.00001
                                          end do
                                       end do
	                                 read(11,*) ! skipping the CN file              
                                    else
			                         read(11,'(a30)') NameFl
                                       open(UNIT=22,FILE=NameFl,
     .								     STATUS='OLD')
                                       call Soil_storativity(22)
                                       close(22)
                                    end if
C                                   TotPio reports the total rainfall [mm] for each cell
                                    do I=1,idim
                                       do J=1,jdim
                                          if (Slope(I,J).GT.0) then
                                             TotPio(I,J)=0
                                          else
                                             TotPio(I,J)=-1
                                          end if
                                       end do
                                    end do
                                    read(11,*)
                                    read(11,'(a30)') NameFl ! file of the precipitation data
                                    open(UNIT=22,FILE=NameFl,
     .                                           STATUS='OLD')
C                                   spatial interpolation of the rainfall according to the value of risp2
                                    call Input_rain(22,chkmet)
                                    close(22)
	                              if (chkmet.LT.98) then 
							         chka=0 
                                       if ((risp.EQ.3).or.(risp.EQ.4))  
     .                                                             then
	                                    read(11,*)
                                          read(11,'(a30)') NameFl ! file of the temperature data
                                          open(UNIT=22,FILE=NameFl,
     .									     STATUS='OLD')
                                          call Input_temp(nvent,chkmet)
                                          close(22)
	                                    if (chkmet.LT.98) then
	                                       totg=0       ! counter of the days
C                                            computation of the number of day from 1/1/aa to gg/mm/aa
                                             do I=1,mm-1  
                                                totg=totg+n(I)
                                             end do
                                             totg=totg+gg
								        else
								           close(11)
	                                       chka=1
                                             if (chkmet.EQ.99) then
										      call Warn_error(1004)
                                             end if							 
									    end if					  	
								     else
								        read(11,*)
	                                    read(11,*)								  
								     end if	 
								     if (chka.EQ.0) then	 																											   
	                                    read(11,*)
	                                    read(11,*) Nmis ! number of cross sections where the 
						                                ! simulated discharge is displayed								
								        if (Nmis.LE.Nsec) then
                                             read(11,*)
                                             do I=1,Nmis
                                                read(11,*) outI(I),
     .                                                     outJ(I)  ! coordinates of the output cross sections
										   end do
		                                   close(11) 
								        else
	                                       chkread=1
								           close(11)
	                                       close(22)
									       call Warn_error(1017) 				  
								        end if	
								     end if		 	 					                              
                                    else
	                                 chkread=1
                                       if (chkmet.EQ.99) then
								        call Warn_error(1005)
								     end if	   				  
					              end if
                                 else
	                              chkread=1
	                              call Warn_error(1006) 	
	                           end if
		                    else
	                           chkread=1
	                           call Warn_error(1007) 
                              end if                           
					     else
	                       chkread=1
	                       call Warn_error(1008) 
	                     end if
                        end if
                     else
	                 chkread=1
	                 call Warn_error(1009)
                     end if
	            else
	              chkread=1 
                   call Warn_error(1010)
	            end if
	         else
	            chkread=1
	            call Warn_error(1015)
	         end if
            else
               chkread=1
			 call Warn_error(1018)   
	      end if
	   else
	     chkread=1
	     close(11)
	     call Warn_error(1016)     
	   end if	  	   
	else 
	   chkread=1
	   call Warn_error(1011)
	end if     	      
			 
	close(11)
	close(22)		 
			 			 
      return
	end


	
C     ********************************************************************************
C
C     SUBROUTINE Slope_comp
C     Computation of the slope of the reaches
C
C     ********************************************************************************

      subroutine Slope_comp(nvent,NRighe)

      implicit none
      
	include 'DimArrays.inc'
	include 'comRead.inc'
	include 'comRainRun.inc'

      integer I0,J0,I1,J1,IS,JS,III,JJJ,chkok,contr,nvent,NRighe,IORD
      real XLEN,PEN


      NTRA=0
      contr=0
C     skipping the first two lines of the file
      read(nvent,*)
      contr=contr+1

      do III=1,idim
         do JJJ=1,jdim
            Slope(III,JJJ)=-1
         end do
      end do

      do while (contr.LT.NRighe)
         read(nvent,*) I0,J0,I1,J1,IORD,XLEN
         contr=contr+1
         NTRA=NTRA+1
         MI0(NTRA)=I0
         MJ0(NTRA)=J0
         MI1(NTRA)=I1
         MJ1(NTRA)=J1
         XLEN=XLEN*((DX+DY)/2)
         PEN=(QUO(I0,J0)-QUO(I1,J1))/XLEN  ! slope of the link
         if (PEN.LE.0) PEN=0.0001
         SLOPE(I0,J0)=PEN
         IS=(PUN(I0,J0)-1)/3-1
         JS=PUN(I0,J0)-5-3*IS
         IS=IS+I0
         JS=JS+J0
         chkok=0
         do while (chkok.EQ.0)
            if (JS.EQ.J1.AND.IS.EQ.I1) then
               chkok=1
            else
               J0=JS
               I0=IS
               SLOPE(I0,J0)=PEN
               IS=(PUN(I0,J0)-1)/3-1
               JS=PUN(I0,J0)-5-3*IS
               IS=IS+I0
               JS=JS+J0
            end if
         end do
      end do

      return
      end


C     *******************************************************************************
C
C     SUBROUTINE Input_rain
C     reading of the rainfall data file. Spatial interpolation of the rainfall data
C
C     *******************************************************************************

      subroutine Input_rain(nvent,chkmet) 

      implicit none

      include 'DimArrays.inc'
	include 'comRead.inc'
	include 'comRainRun.inc'

      integer TotOss,nvent,xpas,I,K,J,nskip,chkmet
      real THp(Nraing)

      read(nvent,*)
	read(nvent,*) Nplu ! number of raingauge stations
	read(nvent,*)
	if (Nplu.GT.Nraing) then
	   call Warn_error(1020)
	   chkmet=98
      else
         do I=1,Nplu
            read(nvent,*) Coordp(1,I),Coordp(2,I) ! coordinates of the raingauges
         end do
	   read(nvent,*)
         read(nvent,*) nskip  ! number of rainfall data to be skipped
	   read(nvent,*)
         read(nvent,*)
c        skipping the rainfall data
         if (nskip.NE.0) then
            do I=1,nskip
               read(nvent,*) 
            end do
         end if
C        reading of the rainfall data
         I=1
         do while (I.LE.NOss)
            read(nvent,*) xpas,(PIO(I,K),K=1,Nplu)
            I=I+1
         end do
         TotOss=I-1

         write(*,*)
         write(*,*) nskip,' steps of rainfall data have been skipped'
         write(*,*)
         write(*,*) NOss,' steps of rainfall data have been read'
         write(*,*)

C        initialization of FGPio that checks which raingauges are operative       
         I=1
	   do while (I.LE.TotOss)
	      chkmet=0
            do J=1,Nplu
               if (PIO(I,J).LT.0) then
                  FGPio(I,J)=1
	            chkmet=chkmet+1
               else
                  FGPio(I,J)=0
               end if
               if (chkmet.EQ.Nplu) then
		        I=TotOss+1   ! forcing the end of the loop   
		        chkmet=99                        
               end if
            end do
	      I=I+1
         end do

         if (chkmet.NE.99) then
            if (risp2.EQ.1) then
C              application of the  method considering the all raingauges
               call Thiessen(Nplu,Coordp,STp,Slope)
            else
C              application of the inverse distance method considering the all raingauges 
               call InvDist(Nplu,Coordp,INVp,Slope)
            end if
         end if
      end if

      return
      end


C     ************************************************************************************
C
C     SUBROUTINE Input_temp
C     reading of the temperature data file. Spatial interpolation of the temperature data
C
C     ************************************************************************************

      subroutine Input_temp(nvent,chkmet)

      implicit none

      include 'DimArrays.inc'
	include 'comRead.inc'
	include 'comRainRun.inc'

      integer I,J,K,nvent,xpas,TotOss,chkmet,nskip

C     reading of the parameters of the evapotranspiration formula
	read(nvent,*)
      read(nvent,*) aE,bE
	read(nvent,*)	
      read(nvent,*) (We(I),I=1,12)    ! monthly compensation factor
	read(nvent,*)	
      read(nvent,*) (Nsol(I),I=1,12)  ! monthly number of hours of sun
	read(nvent,*)
      read(nvent,*) (n(I),I=1,12)     ! number of days in each month
      read(nvent,*)
      read(nvent,*) Nter              ! number of temperature stations
	if (Nter.GT.Ntermo) then
	   call Warn_error(1021)
	   chkmet=98
	else
	   read(nvent,*)
         do I=1,Nter
            read(nvent,*) CoordT(1,I),CoordT(2,I),alt(I) ! coordinates and elevations of the stations
         end do
	   read(nvent,*)
         read(nvent,*) nskip     ! number of temperature data to be skipped
         read(nvent,*)
         read(nvent,*)
         if (nskip.NE.0) then
            do I=1,nskip
               read(nvent,*)
            end do
         end if

         I=1
         do while(I.LE.NOss)
            read(nvent,*) xpas,(gradi(I,K),K=1,Nter)
            I=I+1
         end do

         xpas=NOss

         write(*,*)
         write(*,*) nskip,' steps of temperature data have been skipped'
         write(*,*)
         write(*,*) NOss,' steps of temperature data have been read'
         write(*,*)

         TotOss=I-1
C        inizialization of FGTER that checks which temperature stations are operative
         I=1
         do while (I.LE.TotOss)
	      chkmet=0
            do J=1,Nter
               if (gradi(I,J).LT.-98.0) then
                  FGTer(I,J)=1
	            chkmet=chkmet+1
               else
                  FGTer(I,J)=0
               end if
	         if (chkmet.EQ.Nter) then
		        I=TotOss+1   ! forcing the end of the loop   
		        chkmet=99                        
               end if
            end do
	      I=I+1
         end do

         if (chkmet.NE.99) then
C           application of the Thiessen method the the whole temperature stations
            call Thiessen(Nter,CoordT,STt,Slope)
         end if
      end if
           
	return
      end


C     ********************************************************************************
C
C     SUBROUTINE Soil_Storativity
C     Computation of the S matrix (local soil storativity) starting from the
C     curve number of each cell given in the input file
C
C     ********************************************************************************

      subroutine Soil_storativity(nvent)

      implicit none
	
	include 'DimArrays.inc'
	include 'comRead.inc'
	include 'comRainRun.inc'
	
	integer I,J,nvent
      real TOTCN,SomCN,SomCN2,XMEAN,XVAR,XSTD

      TOTCN=0
      SomCN=0
      SomCN2=0

      do I=1,idim
         read(nvent,*) (CNum(J),J=1,jdim)
         do J=1,jdim
            if (CNum(J).GT.0) then
               if ((risp.EQ.3).OR.(risp.EQ.4)) then ! continuous simulation
C                 the CN numbers given for AMC=2 in the CN file is converted to AMC=1
                  CNum(J)=4.2*CNum(J)/(10-0.058*CNum(J))
               end if
               S(I,J)=254*(100./CNum(J)-1)         ! S=254*(100/CN-1) [mm]
               if (S(I,J).EQ.0) S(I,J)=0.0001
	         if (S(I,J).LT.0) then
			    write(*,*)
	            write(*,*) 'S < 0 for cell', I,J
	            write(*,*) 'Error in the CN file'
	            write(*,*)
			 end if 

               TOTCN=TOTCN+1
               SomCN=SomCN+CNum(J)
               SomCN2=SomCN2+CNum(J)*CNum(J)
            else
               S(I,J)=-1.
            end if
         end do
      end do

      XMEAN=SomCN/TOTCN
      XVAR=(SomCN2/TOTCN-XMEAN*XMEAN)
      XSTD=SQRT(abs(XVAR))
C	write(*,*)
C	write(*,*)
C	write(*,*) 'Curve Number statistics'
C	write(*,*)
C     write(*,*) ' Average', XMEAN
C      write(*,*)
C      write(*,*) ' Variance', XVAR
C      write(*,*)
C      write(*,*) ' Standard deviation', XSTD
C      write(*,*)

      return
      end



C     *******************************************************************************
C
C     SUBROUTINE Thiessen
C     Application of the Thiessen method considering all the raingauge stations
C
C     *******************************************************************************

      subroutine Thiessen(Nstr,Coord,Rif,Slope)

      implicit none
      
	include 'DimArrays.inc'
	include 'comRead.inc'

      real Coord(2,NumStr),TH(NumStr)
      integer Nstr,I,J,K,KK,Rif(Nrow,Ncol)
      real DMIN,Dist,totc,Slope(Nrow,Ncol)

      do I=1,Nstr
         TH(I)=0
      end do
      do I=1,idim
         do J=1,jdim
              Rif(I,J)=-1
              if (Slope(I,J).GT.0) then
                 DMIN=1.E20
                 do K=1,Nstr
                    Dist=(I-Coord(1,K))**2+(J-Coord(2,K))**2
                    if (Dist.LT.DMIN) then
                       DMIN=Dist
                       KK=K
                    end if
                 end do
                 Rif(I,J)=KK
                 TH(KK)=TH(KK)+1
              end if
         end do
      end do

      totc=0
      do I=1,Nstr
         totc=TH(I)+totc
      end do
      do I=1,Nstr
         TH(I)=TH(I)/totc
      end do

      return
      end

C     **********************************************************************************
C
C     SUBROUTINE InvDist
C     Application of the inverse distance method considering all the raingauges
C
C     **********************************************************************************

      subroutine InvDist(Nstr,Coord,INVp,Slope)

      implicit none
	
	include 'DimArrays.inc'
	include 'comRead.inc'
      
      real Coord(2,NumStr),TH(NumStr),sD2(NumStr)
      integer Nstr,I,J,K,indi,totc
      real totinvD2,D2,invD2,INVp(Nrow*Nraing,Ncol)
	real Slope(Nrow,Ncol)

      do I=1,Nstr
         TH(I)=0
      end do
      totc=0    

      do K=1,Nstr
         do J=1,jdim
            do I=1,idim
               indi=I+(K-1)*idim
               INVp(indi,J)=0
            end do
         end do
      end do
      do I=1,idim
         do J=1,jdim
            if (Slope(I,J).GT.0) then
               totc=totc+1
               totinvD2=0
               do K=1,Nstr
                  sD2(K)=0
               end do
               do K=1,Nstr
                  D2=(I-Coord(1,K))**2+(J-Coord(2,K))**2
                  if (D2.EQ.0) D2=0.0001
                  invD2=1/D2    
                  sD2(K)=invD2
                  totinvD2=totinvD2+invD2
               end do
               do K=1,Nstr
                  indi=I+(K-1)*idim
                  INVp(indi,J)=sD2(K)/totinvD2
               end do
            end if
         end do
      end do

      do K=1,Nstr
         do J=1,jdim
            do I=1,idim
               indi=I+(K-1)*idim
               TH(K)=TH(K)+INVp(indi,J)
            end do
         end do
      end do
      do I=1,Nstr
         TH(I)=TH(I)/totc
      end do

      return
      end