      subroutine write_output(npf,totKK)
   
      implicit none

      include 'DimArrays.inc'
      include 'comRead.inc'
      include 'comVar.inc'
      include 'comRainRun.inc'

      integer I,J,K,KK,LL,Iend,totKK,opzst1,opzst2,npf,totmass
      integer contQsup,contQtot,mass
      real Qsupst(Nout,Nsec),Qtotst(Nout,Nsec)
      character*30 NameFl 
	character*3 dum

      OPEN(UNIT=11,FILE='Affdef.in',STATUS='old')
1     read(11,*) dum  
      if (dum.ne.'A67') then
	   goto 1
      else
         if (Nmis.GT.0) then
            do I=1,Nmis
	         read(11,*) NameFl
               open(UNIT=30+I,FILE=NameFl,STATUS='UNKNOWN') 
            end do
            read(11,*)
            read(11,*) opzst1  ! print options
            read(11,*)
            read(11,*) opzst2
            if (opzst2.EQ.2) opzst1=1
               do I=1,Nmis
                  if (opzst2.EQ.1) then
                     write(30+I,*) 'Simulated river flow'
                     if (opzst1.EQ.1) then                
                        write(30+I,*) '   Step    Total Flow'
                     else
                        write(30+I,*) '  Step   Surface Flow    
     .Total Flow'
                     end if
                  else
                     write(30+I,*) 'Simulated annual maximum flow'
                     write(30+I,*) '   Year    Total Flow'
                  end if
              end do 
            close(11)

C           Printing of the results 
            if (StepPio-DT.NE.0.) then
C           the time step of the simulation is an entire submultiple of StepPio.
C           Therefore one computes the discharge at each StepPio step 
               do LL=1,Nmis
                  K=1
                  KK=1
                  do while (K.LT.totKK)   
                     I=K
                     contQsup=0
                     contQtot=0
                     do while (I.LE.(K+(StepPio/DT-1)))
                        contQsup=contQsup+Qsup(I,LL)
                        contQtot=contQtot+Qtot(I,LL)
                        I=I+1
                     end do
                     Qsupst(KK,LL)=contQsup/(StepPio/DT)
                     Qtotst(KK,LL)=contQtot/(StepPio/DT)
                     KK=KK+1
                     K=I
                  end do
               end do
               totKK=totKK/(StepPio/DT)  
            else
               do I=1,totKK
                  do LL=1,Nmis
                     Qsupst(I,LL)=Qsup(I,LL)
                     Qtotst(I,LL)=Qtot(I,LL)
                  end do
               end do
            end if
C           - printing of the simulated discharge
            if (opzst2.EQ.1) then
               do LL=1,Nmis
                  if (opzst1.EQ.1) then
                     do I=1,totKK            
                        write(30+LL,'(1x,I6,3x,F10.5)') 
     .                                I,Qtotst(I,LL)
                     end do
                  else
                     do I=1,totKK
                        write(30+LL,87) I,Qsupst(I,LL),
     .                                 Qtotst(I,LL)
87                      format(1(1x,I6.0),2(3x,F10.5))
                     end do
                  end if
               end do
            else
C           - printing of the simulated annual maximum discharge
            totmass=365*24  
            do LL=1,Nmis
               K=1
               KK=1
               do while (K.LT.totKK)
                  I=K
                  mass=-100000000
                  do while (I.LE.(K+(totmass-1)))
                     if (Qtotst(I,LL).GT.mass) then
                        mass=Qtotst(I,LL)
                     end if
                     I=I+1
                  end do
                  write(30+LL,'(1x,I6.0,3x,F10.5)') KK,mass
                  K=I
                  KK=KK+1
               end do
            end do 
           end if
         end if

C        file of the total local rainfall depth
         OPEN(UNIT=60,FILE='TotRain.out',STATUS='UNKNOWN')
         do I=1,idim
            write(60,'(1000(F12.2))') (TotPio(I,J),J=1,jdim)
         end do
         close(60)

         open(60,FILE='RaingThDistr.out',STATUS='UNKNOWN')
         do I=1,idim
            write(60,'(1000(I3))') (STp(I,J),J=1,jdim)
         end do
         close(60)

         if (risp.ne.2) then
	      open(60,FILE='TempStDistr.out',STATUS='UNKNOWN')
            do I=1,idim
               write(60,'(1000(I3))') (STt(I,J),J=1,jdim)
            end do
            close(60)
         end if
   
         open(UNIT=60,FILE='Slope.out',STATUS='UNKNOWN')
         do I=1,idim
             write(60,'(1000(F8.4,1x))') (Slope(I,J),J=1,jdim)
         end do
         close(60)
      end if

      return
      end
      