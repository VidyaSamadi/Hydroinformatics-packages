C     ********************************************************************************
C
C     subroutine NewThiessen
C     Reapplication of the Thiessen method by excluding the raingauges with missing 
C     data, whose number is stored in manSt()
C
C     ********************************************************************************

      subroutine NewThiessen(Nstr,Coord,totstr,manSt,Rif,Slope)

      implicit none

      include 'DimArrays.inc'
      include 'comRead.inc'

      integer Nstr,totstr,I,J,L,K,KK,chkst,Rif(Nrow,Ncol),manSt(NumStr)
      real DMIN,Dist,Slope(Nrow,Ncol),Coord(2,NumStr)

      do I=1,idim
         do J=1,jdim
            Rif(I,J)=-1
            if (Slope(I,J).GT.0) then
               DMIN=1.E20
               do K=1,Nstr
                  chkst=0
                  do L=1,totstr
                     if (K.EQ.manSt(L)) chkst=1
                  end do
                  if (chkst.EQ.0) then
                     Dist=(I-Coord(1,K))**2+(J-Coord(2,K))**2
                     if (Dist.LT.DMIN) then
                        DMIN=Dist
                        KK=K
                    end if
                  end if
               end do
               Rif(I,J)=KK
            end if
         end do
      end do

      return
      end