C     **************************************************************************
C
C     *************     S U B R O U T I N E     M O D 1            *************
C
C
C     **************************************************************************

      subroutine MOD1

      implicit none

      include 'DimArrays.inc'
      include 'VarMOD1.inc'
      include 'comread.inc'

      integer L,ios,chkgo,chk
      character*30 NameFl,NameFl2,NameFl3

      chkgo=0
      open(UNIT=30,FILE='Affdef.in',IOSTAT=ios,STATUS='OLD')
      if (ios.EQ.0) then
         do L=1,11
          read(30,*)
         end do  
         read(30,'(a30)') NameFl
C        opening of the DEM file
         open(UNIT=20,FILE=NameFl,IOSTAT=ios,STATUS='OLD')
         if (ios.EQ.0) then
            do L=1,5
               read(30,*)
            end do
            read(30,*) idim,jdim   ! reading of the dimension of the DEM
	      if ((idim.LE.NRow).AND.(jdim.LE.NCol)) then
               do L=1,27
                  read(30,*)
               end do
               read(30,*) Ichi,Jchi   ! reading of the coordinates of the outlet
               do L=1,5
                  read(30,*)
               end do
               read(30,'(a30)') NameFl  ! file of the slope pointer

               open(UNIT=19,FILE=NameFl,STATUS='UNKNOWN')
               call Flow_directions(20,19,chkgo,30) 
	         if (chkgo.EQ.0) then 
                  read(30,*)          
                  read(30,'(a30)') NameFl2  ! file of the contributing area

                  open(UNIT=13,FILE=NameFl2,STATUS='UNKNOWN')
                  call Contributing_Area(13,ORD)
                  close(13)

C                  write(*,*)
C                  write(*,*) 'Computation of the contributing area 
C     .completed'
C                  write(*,*)

                  open(UNIT=12,FILE='RiverOrder.out',STATUS='UNKNOWN') ! file of the link orders
C                 opening of the file that will contain information about the river network organised in
C                 reaches according to the Strahler classification. (LinkProperties_temp.out is a temporary 
C                 file; it will be stored later as NameFl3.out)
                  open(UNIT=13,FILE='LinkProperties_temp.out',STATUS=
     .                                          'UNKNOWN')
C                 Classification of the river network
                  call ORDINI(12,13,ORD,MORD,chk)
                  close(12)
	            if (chk.EQ.1) then
                     call Warn_error(1003)
	               close(13)
	            else
C                     write(*,*) 'Classification of the river network 
C     .completed'
C                     write(*,*)

C                    determination of the matrix which represents the mask of the catchment
                     call Mask(BAC,ORD)

                     read(30,*)
                     read(30,'(a30)') NameFl3  ! file of the link properties
                     open(UNIT=17,FILE=NameFl3,STATUS='UNKNOWN')
C                    determination of the properties of the river links
                     call Param(17,MORD,BAC,ORD)
                     close(17)
                     close(13,STATUS='DELETE')
C                     write(*,*) 'The river network determination was 
C     .successful'
C                     write(*,*)
                  end if
               end if
            else
               call Warn_error(1018)
               close(13)
	      end if
         else
           call Warn_error(1001)
         end if
      else
         call Warn_error(1011)
      end if

      close(19)
      close(20)
      close(30)

      return
      end


C     ********************************************************************************
C
C     subroutine Flow_directions
C     Determination of the flow direction (D8 method)
C
C     ********************************************************************************

      subroutine Flow_directions(flin,flout,chkgo,ntre)

      implicit none

      include 'DimArrays.inc'
      include 'comread.inc'

      integer flin,flout,ntre,ios,chkgo,I,J,K,II,JJ,I1,I2,J1,J2,npf
      integer npfi(npun),npfj(npun),numc,totq,MedQ
      real Z,ZI,PEN,P,MAXQUO
      character*30 NameFl

      numc=0
      totq=0

      do I=1,idim
         read(flin,*) (QUO(I,J),j=1,jdim)
         do J=1,jdim
C           if the elevation is less than 0, then the cell is
C           outside the catchment, therefore QUO=-1
            if (QUO(I,J).LT.0) then
               QUO(I,J)=-1.
            else
               totq=totq+QUO(I,J)
               numc=numc+1
            end if
         end do
      end do

C     computation of the average elevation of the catchment
      MAXQUO=-1000
      do I=1,idim
         do J=1,jdim
            if (QUO(I,J).GT.MAXQUO) then
               MAXQUO=QUO(I,J)
            end if
        end do
      end do

      MAXQUO=MAXQUO+100.0

C     average elevation of the catchment
      MedQ=totq/numc

C      write(*,*) 'Average elevation of the catchment - ',MedQ
C      write(*,*)

      read(ntre,*)
      read(ntre,*) npf  ! number of cells whose elevation has not to be increased
	if (npf.LE.npun) then       
         read(ntre,*)
         do I=1,npf
            read(ntre,*) npfi(I),npfj(I)  ! reading of the coordinates of the npf cells
         end do

C        subroutine Unpit to remove the pits
         call Unpit(npf,npfi,npfj,MAXQUO,chkgo)

         if (chkgo.EQ.0) then
            read(ntre,*)
            read(ntre,'(a30)') NameFl
C           File where the modified DEM is stored
            open(UNIT=130,FILE=NameFl,STATUS='UNKNOWN')
            do I=1,idim
               write(130,61) (QUO(I,J),J=1,jdim)
61             FORMAT(945F10.1)
            end do

C           8 possible outflow directions from a given cell
C
C                               1-----2-----3
C                               ! +   !   + !
C                               !   + ! +   !
C                          I    4-----+-----6
C                               !   + ! +   !
C                               ! +   !   + !
C                               7-----8-----9
C
C                                     J

            do I=1,idim
               I1=-1
               I2=1
               if (I.EQ.1) I1=0
               if (I.EQ.idim) I2=0
               do J=1,jdim
                  PUN(I,J)=0
                  if (QUO(I,J).GT.0) then
                     J1=-1
                     J2=1
                     if (J.EQ.1) J1=0
                     if (J.EQ.jdim) J2=0
                     Z=QUO(I,J)
                     PEN=-1.E20
                     do II=I1,I2
                        do JJ=J1,J2
                           if (QUO(I+II,J+JJ).GT.0) then
                              P=QUO(I+II,J+JJ)-Z
                              P=-P
                              if (P.GT.0) then
                                 if (II*JJ.EQ.0) then
                                              
C                                  | quo() |
C                            ------|-------|-------
C                             quo()|   z   | quo()
C                            ------|-------|-------
C                                  | quo() |

                                    if (P.GT.PEN) then
                                        PEN=P
                                        PUN(I,J)=3*II+JJ+5
                                    end if
                                 else
                                    if(QUO(I,J+JJ).LT.0.OR.
     .                                 QUO(I+II,J).LT.0) then

C                                        | -1 |
C                                   -----|----|-----
C                                     -1 |  z | -1
C                                   -----|----|-----
C                                        | -1 |

                                       P=P/sqrt(2.)
                                       GOTO 41
                                    end if

                                    if (QUO(II+I,J+JJ).GE.QUO(I,J+JJ))
     .                                 GO TO 7
                                    if (QUO(II+I,J+JJ).GE.QUO(II+I,J))
     .                                 GO TO 7

C                                         |  quo()| quo()
C                                    -----|-------|-------
C                                         |   z   | quo()
C                                    -----|-------|-------
C                                         |       |

                                    ZI=.25*(Z+QUO(I+II,J+JJ)+
     .                                  QUO(I,J+JJ)+QUO(I+II,J))
                                    P=(ZI-Z)/SQRT(2.)
                                    P=-P
41                                  if (P.GT.PEN) then
                                        PEN=P
                                        PUN(I,J)=3*II+JJ+5
                                    end if
                                 end if
                              end if
                           end if
7                       end do
                     end do
                  end if
              end do
            end do
            do I=1,idim
               write(flout,'(945I2)') (PUN(I,K),K=1,jdim)
            end do
            close(130)
          end if
      else
         chkgo=1
	   call Warn_error(1019)
	   close(13)
	end if
      return
      end


C     ********************************************************************************
C
C     subroutine Unpit
C     The subroutine removes the pits starting from the original DEM. It produces
C     a modified DEM stored in a file whose name is specified in the file
C     Affdef.in
C
C     *******************************************************************************

      subroutine Unpit(npf,npfi,npfj,MAXQUO,chkgo)

      implicit none

      include 'DimArrays.inc'
      include 'comread.inc'

      integer npf,vpf,chkpit,chkgreat,chkgo,npfi(npun),npfj(npun)
      integer I1,I2,J1,J2,KKK,I,L,II,III,J,JJ,JJJ
      real zpas,ZX,ZMIN,MAXQUO

      zpas=0.11 ! increasing rate of the elevation of the given pit
      I1=-1
      I2=1
      J1=-1
      J2=1
      I=2
      J=2
      KKK=0
      chkpit=1
      do while (chkpit.EQ.1)
         chkpit=0
         do while (I.LE.idim-1)
            do while (J.LE.jdim-1)
               vpf=0
               if (QUO(I,J).GE.0) then
                  if (npf.GT.0) then
                     do L=1,npf
                        if(I.EQ.npfi(l).AND.J.EQ.npfj(l)) vpf=1
                     end do
                  end if
                  if (vpf.EQ.0) then
                     ZMIN=1.E20
                     ZX=QUO(I,J)
                     chkgreat=0
                     do III=I1,I2
                        do JJJ=J1,J2
                           if (III.NE.0.OR.JJJ.NE.0) then
                              II=III+I
                              JJ=JJJ+J
                              if (QUO(II,JJ).GT.0) then
                                 if (ZX.GT.QUO(II,JJ)) then
                                    chkgreat=1
                                 else
                                    if(ZMIN.GT.QUO(II,JJ))
     .                                 ZMIN=QUO(II,JJ)
                                 end if
                              end if
                           end if
                       end do
                     end do
                     if (chkgreat.eq.0) then
                        if (ZX.LE.ZMIN) then
                           chkpit=1
                           QUO(I,J)=ZMIN+zpas
                           if (QUO(I,J).GE.MAXQUO) then
	                        write(*,*) 'Elevation of the cell ',I,J,
     .' greater than the maximum elevation of the catchment.'
                              call Warn_error(1002)
                              I=idim
                              J=jdim
                              chkgo=1
                           else
                              KKK=KKK+1
                              if (I.NE.2) I=I-1
                              if (J.NE.2) J=J-2
                              if (KKK/500*500.EQ.KKK) then
                                 write(*,'(I10,2I5,F30.2)') KKK,I,J,
     .                                                      QUO(I,J)
                              end if
                           end if
                        end if
                     end if
                  end if
               end if
               J=J+1
            end do
            J=2
            I=I+1
         end do
      end do

      return
      end
 


C     *******************************************************************************
C
C     subroutine Contributing_Area
C     Computation of the contributing areas on the basis of the flow paths
C
C     ********************************************************************************

      subroutine Contributing_Area(ntred,ORD)

      implicit none

      include 'DimArrays.inc'
      include 'comread.inc'

      integer I,J,I0,J0,KK,ntred,ORD(Nrow,Ncol)
      do I=1,idim
         do J=1,jdim
            AREE(I,J)=1
            ORD(I,J)=1
         end do
      end do

      KK=1
      do WHILE (KK.EQ.1)
         KK=0
         do I=1,idim
            do J=1,jdim
               if (ORD(I,J).NE.0) then
                  if (PUN(I,J).NE.0) then
                     I0=(PUN(I,J)-1)/3-1
                     J0=PUN(I,J)-5-3*I0
                     I0=I0+I
                     J0=J0+J
                     ORD(I0,J0)=ORD(I0,J0)+ORD(I,J)
                     AREE(I0,J0)=AREE(I0,J0)+ORD(I,J)
                     ORD(I,J)=0
                     KK=1
                  end if
               end if
            end do
         end do
      end do

      do I=1,idim
         do J=1,jdim
            if (ORD(I,J).EQ.1) AREE(I,J)=0
         end do
      end do

      do I=1,idim
         write(ntred,'(945(i7,1x))') (AREE(I,J),J=1,jdim)
      end do

      return
      end


C     *******************************************************************************
C
C     subroutine ORDINI
C     Classification of the river network according to the Strahler scheme
C
C     ********************************************************************************

      subroutine ORDINI(ndod,ntred,ORD,MORD,chk)

      implicit none

      include 'DimArrays.inc'
      include 'comread.inc'

      integer chk,ndod,ntred,punt0,chkend,IORD,IORD1,NMORD,ICONT,IERR
      integer aI2,aJ2,I,J,II,JJ,III,JJJ,IIN,JIN,IFIN,JFIN,I0,J0,I2,J2
      integer NN,IER,MORD(NRivi,4),ORD(Nrow,Ncol)
      real XLEN


      do I=1,idim
         do J=1,jdim
            if (AREE(I,J).EQ.1) then
               ORD(I,J)=1
               if (PUN(I,J).EQ.0) ORD(I,J)=-999
            else
               ORD(I,J)=0
            end if
         end do
      end do
      IORD=1
      IORD1=0
      NMORD=1
      ICONT=1
      do while (ICONT.EQ.1)
         IERR=0
         if (IORD1.NE.IORD) then
            write(ntred,'('' Order '',I2,/,''  IIN JIN IFIN JFIN  ORD'',
     .             ''    LEN       '')') IORD
            IERR=1
         end if
         IORD1=IORD
         ICONT=0
         I=1
         J=1
         do while (I.LE.idim)
            J=1
            do while (J.LE.jdim)
             chk=0
             if (NMORD.GT.NRivi) then
                I=idim+1  ! forcing the end of the loops
                J=jdim+1   
                chk=1            
             else
                if (ORD(I,J).EQ.IORD) then
                     ORD(I,J)=-IORD
                     ICONT=1
                     if (IERR.EQ.1) then
                        MORD(NMORD,1)=IORD
                        MORD(NMORD,2)=I
                        MORD(NMORD,3)=J
                        NMORD=NMORD+1
                     end if

                     III=I
                     JJJ=J
                     IIN=I
                     JIN=J
                     IFIN=I
                     JFIN=J
                     XLEN=0
                     if (PUN(III,JJJ).EQ.0) then
                        ORD(III,JJJ)=-1000
                     else
                        I0=(PUN(III,JJJ)-1)/3-1
                        J0=PUN(III,JJJ)-5-3*I0
                        if ((I0*J0).NE.0) then
                           XLEN=XLEN+SQRT(2.)
                        else
                           XLEN=XLEN+1
                        end if
                        IFIN=III+I0
                        JFIN=JJJ+J0
                        punt0=0
                        do while (AREE(IFIN,JFIN).EQ.AREE(III,JJJ)+1
     .                           .AND.punt0.EQ.0)
                           if (PUN(IFIN,JFIN).EQ.0) then
                              ORD(IFIN,JFIN)=-1000
                              punt0=1
                           else
                              ORD(IFIN,JFIN)=-IORD
                              III=IFIN
                              JJJ=JFIN
                              I0=(PUN(III,JJJ)-1)/3-1
                              J0=PUN(III,JJJ)-5-3*I0
                              if ((I0*J0).NE.0) then
                                 XLEN=XLEN+SQRT(2.)
                              else
                                 XLEN=XLEN+1
                              end if
                              IFIN=IFIN+I0
                              JFIN=JFIN+J0
                           end if
                        end do
                        if (punt0.EQ.0) then
                           if (PUN(IFIN,JFIN).EQ.0) then
                              ORD(IFIN,JFIN)=-1000
                           else
                              ORD(IFIN,JFIN)=IORD+1000
                           end if
                        end if
                     end if
                     write(ntred,'(2I4,1X,2I4,1X,I4,1X,F13.1)')IIN,JIN,
     .                                IFIN,JFIN,IORD,XLEN
                  end if
             end if
               J=J+1      
          end do
          I=I+1
         end do

         if (chk.EQ.0) then
          if (ICONT.NE.0) then
               if (IORD.EQ.1) then
                  IORD=IORD+1
                  do III=1,idim
                     do JJJ=1,jdim
                        if (ORD(III,JJJ).GT.1000) then
                           I=-1
                           if (III.EQ.1) I=0
                           I2=1
                           if (III.EQ.idim) I2=0
                           J2=1
                           if (JJJ.EQ.jdim) J2=0
                           chkend=0
                           aI2=I2
                           aJ2=J2
                           do while (I.LE.I2)
                              J=-1
                              if (JJJ.EQ.1) J=0
                              do while (J.LE.J2)
                                 NN=3*I+5+J
                                 II=III+I
                                 JJ=JJJ+J
                                 NN=10-NN
                                 if (PUN(II,JJ).EQ.NN) then
                                    if (ORD(II,JJ).LE.-999.OR.
     .                                 ORD(II,JJ).GE.0) then
                                       chkend=1
                                       I=aI2
                                       J=aJ2
                                    end if
                                 end if
                                 J=J+1
                              end do
                              I=I+1
                           end do
                           if (chkend.EQ.0) ORD(III,JJJ)=IORD
                        end if
                     end do
                  end do
               else
                  call TEST(IORD,IER,ORD)
                  if (IER.EQ.0)IORD=IORD+1
               end if
            end if
         else
          ICONT=99 ! forcing the end of the loop
         end if
      end do

      if (chk.EQ.0) then
         do I=1,idim
            do J=1,jdim
               if (ORD(I,J).GT.-999) then
                  ORD(I,J)=-ORD(I,J)
               else
                  ORD(I,J)=ORD(I,J)/10
                  if ((I.EQ.Ichi).AND.(J.EQ.Jchi)) then
                     ORD(I,J)=IORD-1
                  end if
               end if
            end do
            write(ndod,'(945(I2,1x))') (ORD(I,J),J=1,jdim)
         end do
         MORD(NMORD,1)=-1
      end if

      return
      end


C     *******************************************************************************
C
C     subroutine TEST
C     Classification of the river network according to the Strahler ordering system
C
C     ********************************************************************************

      subroutine TEST(IORD,IER,ORD)

      implicit none

      include 'DimArrays.inc'
      include 'comread.inc'

      integer IORD,IER,IPRI,ISEC,NN,chkend,ORD(Nrow,Ncol)
      integer I,J,I2,J2,II,JJ,III,JJJ,aI2,aJ2

      IER=0
      do III=1,idim
         do JJJ=1,jdim
            if (ORD(III,JJJ).GT.1000) then
               IPRI=0
               ISEC=0
               I=-1
               if (III.EQ.1) I=0
               I2=1
               if (III.EQ.idim) I2=0
               J2=1
               if (JJJ.EQ.jdim) J2=0
               aI2=I2
               aJ2=J2
               chkend=0
               do while (I.LE.I2)
                  J=-1
                  if (JJJ.EQ.1) J=0
                  do while (J.LE.J2)
                     NN=3*I+5+J
                     II=III+I
                     JJ=JJJ+J
                     NN=10-NN
                     if (PUN(II,JJ).EQ.NN) then
                        if (ORD(II,JJ).LE.-999.OR.ORD(II,JJ).GE.0) then
                           chkend=1
                           I=aI2
                           J=aJ2
                        else
                           if (-ORD(II,JJ).GT.IPRI) then
                              IPRI=-ORD(II,JJ)
                           else
                              if (-ORD(II,JJ).GT.ISEC) ISEC=-ORD(II,JJ)
                           end if
                        end if
                     end if
                     J=J+1
                  end do
                  I=I+1
               end do
               if (chkend.EQ.0) then
                  if (IPRI.EQ.ISEC) then
                     ORD(III,JJJ)=IPRI+1
                  else
                     ORD(III,JJJ)=IPRI
                     IER=1
                  end if
               end if
            end if
         end do
      end do

      return
      end


C     ********************************************************************************
C
C     subroutine Mask
C     Creation of the matrix BAC, which represents the mask of the catchment
C
C     ********************************************************************************

      subroutine Mask(BAC,ORD)

      implicit none

      include 'DimArrays.inc'
      include 'comread.inc'

      integer chkok,chkend,chkout,ORD(Nrow,Ncol),BAC(Nrow,Ncol)
      integer I,J,II,JJ,I0,J0

      do I=1,idim
         do J=1,jdim
            BAC(I,J)=0
         end do
      end do
      do I=1,idim
         do J=1,jdim
            chkok=0
            chkout=0
            if (AREE(I,J).EQ.1.AND.PUN(I,J).NE.0) then
               II=I
               JJ=J
               BAC(II,JJ)=-1
               do while (chkok.EQ.0.AND.chkout.EQ.0)
                  if (PUN(II,JJ).NE.0) then
                      I0=(PUN(II,JJ)-1)/3 -1
                      J0=PUN(II,JJ)-5-3*I0
                      JJ=JJ+J0
                      II=II+I0
                      if (BAC(II,JJ).EQ.-1) then
                         chkout=1
                      else
                         if (BAC(II,JJ).GT.0) then
                            chkok=1
                         else
                            BAC(II,JJ)=-1
                            if (II.EQ.Ichi.AND.JJ.EQ.Jchi) then
                               chkok=1
                            end if
                         end if
                      end if
                  else
                      chkout=1
                  end if
               end do
               if (chkout.EQ.0) then
                  II=I
                  JJ=J
                  BAC(II,JJ)=1
                  chkend=0
                  do while (chkend.EQ.0)
                     I0=(PUN(II,JJ)-1)/3 -1
                     J0=PUN(II,JJ)-5-3*I0
                     JJ=JJ+J0
                     II=II+I0
                     if (BAC(II,JJ).GT.0) then
                        chkend=1
                     else
                        BAC(II,JJ)=1
                        if (II.EQ.Ichi.AND.JJ.EQ.Jchi) then
                           chkend=1
                        end if
                     end if
                  end do
               end if
            else
               if (AREE(I,J).EQ.1) BAC(I,J)=-99
            end if
         end do
      end do
      do I=1,idim
         do J=1,jdim
            if (BAC(I,J).EQ.1) BAC(I,J)=ORD(I,J)
         end do
      end do

      return
      end


C     ********************************************************************************
C
C     subroutine Param
C     Computation of the characteristics of the links of the river network
C
C     ********************************************************************************

      subroutine Param(ndic,MORD,BAC,ORD)

      implicit none

      include 'DimArrays.inc'
      include 'comread.inc'

      integer NN,IIN,JIN,IFIN,JFIN,III,JJJ,I0,J0,aIORD,IORD,ndic
      integer ORD(Nrow,Ncol),BAC(Nrow,Ncol),Mord(NRivi,4)
      real XLEN

      NN=1
      write(ndic,'(''   IIN JIN   IFIN JFIN  ORD   LEN'')')
      do while (MORD(NN,1).NE.-1)
         IIN=MORD(NN,2)
         JIN=MORD(NN,3)
         if (BAC(IIN,JIN).GE.0) then
            III=IIN
            JJJ=JIN
            XLEN=0
            IORD=MORD(NN,1)
            aIORD=IORD
            do while (ORD(III,JJJ).EQ.aIORD)
               if (PUN(III,JJJ).NE.0) then
                  I0=(PUN(III,JJJ)-1)/3-1
                  J0=PUN(III,JJJ)-5-3*I0
                  if ((I0*J0).NE.0) then
                     XLEN=XLEN+SQRT(2.)
                  else
                     XLEN=XLEN+1.
                  end if
                  III=III+I0
                  JJJ=JJJ+J0
               else
                  aIORD=aIORD+1
               end if
            end do
            JFIN=JJJ
            IFIN=III
            NN=NN+1
            write(ndic,'(1X,I4,1X,I4,2X,I4,1X,I4,1X,I4,1X,F6.1)')
     .        IIN,JIN,IFIN,JFIN,IORD,XLEN
            if (IORD.NE.1) then
               call RIVER(IIN,JIN)
            end if
         else
            NN=NN+1
         end if
      end do
      return
      end


C     ********************************************************************************
C
C     subroutine RIVER
C
C     ********************************************************************************

      subroutine RIVER(III,JJJ)
      
      implicit none

      include 'DimArrays.inc'
      include 'comread.inc'

      integer chkout,I,I2,J,J2,NN,III,JJJ,IK,JK,II,JJ
      real ZX

      chkout=0
      do while (chkout.EQ.0)
         chkout=1
         I=-1
         I2=1
         if (III.EQ.1) I=0
         if (III.EQ.idim) I2=0
         ZX=-1.E36
         do while (I.LE.I2)
            J=-1
            J2=1
            if (JJJ.EQ.1) J=0
            if (JJJ.EQ.jdim) J2=0
            do while (J.LE.J2)
               NN=3*I+5+J
               II=III+I
               JJ=JJJ+J
               NN=10-NN
               if (PUN(II,JJ).EQ.NN) then
                  if (AREE(II,JJ).GE.ZX) then
                     ZX=AREE(II,JJ)
                     IK=II
                     JK=JJ
                     chkout=0
                     J=J+1
                  else
                     J=J+1
                  end if
               else
                  J=J+1
               end if
            end do
            I=I+1
         end do
         if (chkout.EQ.0) then
            III=IK
            JJJ=JK
          end if
      end do

      return
      end
