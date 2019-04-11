      subroutine Initial_value

      implicit none

      include 'DimArrays.inc'
      include 'comRead.inc'
      include 'comVar.inc'
      integer I,J


C     initialization of global array
      do I=1,idim
       do J=1,jdim
            PioL1(I,J)=0; F(I,J)=0.; ES(I,J)=0.; Inter(I,J)=0.
            ETp(I,J)=0; ETp1(I,J)=0; TE(I,J)=0.; W(I,J)=0.
            Q(I,J)=0.;QQ(I,J)=0.;QI(I,J)=0.;QQI(I,J)=0.
            Y(I,J)=0.;YY(I,J)=0.;YI(I,J)=0.;YYI(I,J)=0.
            C1(I,J)=0.;C2(I,J)=0.;C3(I,J)=0.;C4(I,J)=0.
            CI1(I,J)=0.;CI2(I,J)=0.;CI3(I,J)=0.;CI4(I,J)=0.
         end do
      end do

      return
      end