C     ********************************************************************************
C
C     subroutine Check_temp
C     Check of which temperature stations are not operative
C
C     ********************************************************************************

      subroutine Check_temp(PasPio,K,Nstep,chkTE,STtn,manT)

      implicit none

      include 'DimArrays.inc'
      include 'comRead.inc'
      include 'comRainRun.inc'

      integer chkTE,totTE,p,K,Nstep,PasPio
      integer STtn(Nrow,Ncol),manT(Ntermo)

      chkTE=0
      totTE=0
      do p=1,Nter
         if (FGTer(PasPio,p).EQ.1) then
            chkTE=1
            totTE=totTE+1
            manT(totTE)=p
         end if
      end do

      if (chkTE.EQ.1) then
C        application of the Thiessen method by escluding the temperature stations with missing data
         call NewThiessen(Nter,CoordT,totTE,manT,STtn,Slope)
      end if

      return
      end