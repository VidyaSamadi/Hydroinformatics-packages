C     ********************************************************************************
C
C     subroutine NewInvDist
C     Reapplication of the inverse distance by excluding the raingauges with missing 
C     data, whose number is stored in manSt()
C
C     ********************************************************************************

      subroutine NewInvDist(Nstr,Coord,totstr,manSt,INVpn,Slope)
      implicit none

      include 'DimArrays.inc'
      include 'comRead.inc'

      integer Nstr,totstr,I,J,K,L,chkst,indi,manSt(NumStr)
      real sD2(NumStr),totinvD2,D2,invD2,Coord(2,NumStr)
      real INVpn(Nrow*Nraing,Ncol),Slope(Nrow,Ncol)

   
      do K=1,Nstr
         do J=1,jdim
            do I=1,idim
               indi=I+(K-1)*idim
               INVpn(indi,J)=0
            end do
         end do
      end do

      do I=1,idim
         do J=1,jdim
            if (Slope(I,J).GT.0) then
               totinvD2=0
               do K=1,Nstr
                  sD2(K)=0
               end do
               do K=1,Nstr
                  chkst=0
                  do L=1,totstr
                     if (K.EQ.manSt(L)) chkst=1
                  end do
                  if (chkst.EQ.0) then
                     D2=(I-Coord(1,K))**2+(J-Coord(2,K))**2
                     if (D2.EQ.0) D2=0.0001
                     invD2=1/D2   
                     sD2(K)=invD2
                     totinvD2=totinvD2+invD2
                  end if
               end do

               do K=1,Nstr
                  chkst=0
                  do L=1,totstr
                     if (K.EQ.manSt(L)) chkst=1
                  end do
                  indi=I+(K-1)*idim
                  if (chkst.EQ.0) then
                     INVpn(indi,J)=sD2(K)/totinvD2
                  else
                     INVpn(indi,J)=0
                  end if
               end do
            end if
         end do
      end do

      return
      end