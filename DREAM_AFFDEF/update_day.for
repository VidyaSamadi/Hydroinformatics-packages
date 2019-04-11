      subroutine Update_day(contg,totale,m)

      implicit none

      include 'DimArrays.inc'
      include 'comRead.inc'

      integer m,n1,n2,chkm,L
      real contg,totale

      contg=contg+DT/StepPio  ! counter for the current day of simulation
      if ((contg-24.).EQ.0.) then
         contg=contg-24
         totale=totale+1
         if (totale.gt.365) totale=totale-365  ! new year
         if (totale.le.31) then
            m=1
         else
            if (totale.gt.334) then
               m=12
            else
               n1=n(1)
               n2=n(1)+n(2)
               chkm=0
               L=2
               do while ((chkm.EQ.0).AND.(L.LE.11))
                  if ((totale.GT.n1).AND.(totale.LE.n2)) then
                     m=L
                     chkm=1
                  else
                     n2=n2+n(L+1)
                     L=L+1
                  end if
               end do
            end if
         end if
      end if

      return 
      end