************************************************************************
*                                                                      *
*                            PROGRAMMA ARNO                            *
*     DAILY DATA VERSION                                               *
*----------------------------------------------------------------------*
*                                                                      *
*       IL PROGRAMMA ARNO COSTRUISCE I DEFLUSSI IN UNA SEZIONE         *
*        RICHIESTA PARTENDO DAGLI AFFLUSSI E DALLE TEMPERATURE         *
*                          MEDIE RAGGUAGLIATI                          *
*                                                                      *
*    DATI IN INGRESSO:                                                 *
*                                                                      *
*    IRAM         NUMERO DEI RAMI TROVATI                              *
*    NCOD         VETTORE CONTENENTE I CODICI ORDINATI                 *
*               (wetore containing the orderycods)                     *
*    NEL          NUMERO DI ELEMENTI PER OGNI RAMO                     *
*                 ((i number of elements for evry branch))             *
*    IMEM         NUMERO D'ORDINE DELL'ELEMENTO DI NCOD CUI            *
*                 SI ATTACCA IL PRIMO ELEMENTO DEL RAMO GENERICO       * 
*                  ((serial number of the element of ncod that         *
*                      it is attached the first element of             *
*                          the generic branch ))                       *
*    ISEZ         CODICI IN CUI SI VOGLIONO I DEFLUSSI                 *
*                                                                      *
*                 PIOGGE                                               *
*                 TEMPERATURE                                          *
*                 DATI DI BACINO                                       *
*                       ((codes in wich want themselves the deflussi   *
*                         rains field of data tempratore))             *
*                                                                      *
*                                                                      *
*    DATI IN USCITA: (data in supporting )                             *
*                                                                      *
*                 PORTATE NELLE SEZIONI RICHIESTE                      * 
*     ((exite in the section requested))                               *
*                                                                      *
************************************************************************
*                                                                      *
      REAL*8 DX,C,DC
      INTEGER*2 YY,MM,DD,HH,YYF,MMF,DDF,HHF
      CHARACTER*3 NCOD,NISEZ,NCODI
      CHARACTER NOME(30),TYP
*
      COMMON/ RUNF   /RUN(10000),XR(10000),XRB(10000),TYP
      COMMON/ NUMBER /N,NT,NS,ISERB
      COMMON/ STO    /STORAGE(10000,2),KSTOR,ICOL(30),IPOSOLD(30)
     .               ,NINOLD(30)
      COMMON/ INTERV /NIN,IRAM
      COMMON/ ALB    /NRAM,NEL(30),NCOD(30),IMEM(30)
      COMMON/ INP    /NISEZ(30),NSZ
      COMMON/ PRINT  /RUNOFF(10000)
      COMMON/ UP     /NCODI,W0,W0N,YINIZ
      COMMON/ CONVL  /DX(10,100),C(10,100),DC(10,100),NPAR,IUL,ICB
      COMMON/ DEFR   /NELINI,NELFIN,ISEZIONE
      COMMON/ STR    /ICOL1,IPOS
      COMMON/ DAT    /YY,MM,DD,HH,YYF,MMF,DDF,HHF
      COMMON/ EVAPO  /PET(10000),AET(10000)
      COMMON/ PIOGGE /PP(10000)
      COMMON/ NAME   /NOME
      COMMON/ NEVE   /TS
      COMMON/ PARNAN /BASE,WM,B,D,D1,AREA(100),NVET,ICOD,SOL,SOL1,PERC,
     $CESP,ASC,QC,PC,W0B
      COMMON/ DATAST /ISW(30),IORD(30)
*
******************************************************************
*
C DEBUG VARIABLES
      CALL READ
*
      CALL INPUT
*
      KSTOR=0
*
*---------> SI COMINCIA DALL'ULTIMO NODO DEL PRIMO RAMO CONTENENTE UNO
*           DEI CODICI RICHIESTI. QUESTO RAMO E' IL PRIMO DEI PERCORSI
*           PRINCIPALI  ((from the last knot of the first branch containing
*                           one of the codes richtesti .questo branch is the
*                             first one of the main route))
*
      JRAM=NIN+1
      IPOS=NEL(NIN)+1
10    CONTINUE
*
*---------> SI PERCORRE ALL'INDIETRO UN PERCORSO PRINCIPALE. SI PASSA
*           DA UN NODO AL PRECEDENTE SOLO QUANDO SONO STATI ESAURITI
*           TUTTI I PERCORSI A MONTE DEL NODO. SI ABBANDONA IL PERCORSO
*           PRINCIPALE SOLO SE SI TROVA UN SOTTOBACINO PARALLELO
* ((it is traveled to the back a rout principle 
*    it is gone from a knot alone prenedente guandi were
* exhausted all of the upstream rout of the knot 
*    it is left the rout alone principle if it is found
*             a sattobacino parallel.))
      IPOS=IPOS-1
      NELI=0
      IF(NIN.NE.1)NELI=NEL(NIN-1)
*
*---------> ((if it is exhausted a route principle is added the result of
*            such route to that of the route previous principle immagazzianto
*            in storage))
*
      IF(IPOS.LT.NELI+1) GO TO 80
      NCODI=NCOD(IPOS)
      ISWITCH=ISW(IORD(IPOS))
*
*---------> SI CONTROLLA SE IL BACINO DI CODICE "NCOD" E' DA SALTARE
*             (it is monitored if the field of code "ncod"e'da jump)

*
      IF(ISWITCH.EQ.0) GO TO 50
*
      CALL BACINO
*
      CALL UPDATE
*
      COMMON/ RUNF   /RUN(10000),XR(10000),XRB(10000),TYP
      COMMON/ NUMBER /N,NT,NS,ISERB
      COMMON/ ALB    /NRAM,NEL(30),NCOD(30),IMEM(30)
      COMMON/ UP     /NCODI,W0,W0N,YINIZ
      COMMON/ CONVL  /DX(10,100),C(10,100),DC(10,100),NPAR,IUL,ICB
      COMMON/ DAT    /YY,MM,DD,HH,YYF,MMF,DDF,HHF
      COMMON/ EVAPO  /PET(10000),AET(10000)
      COMMON/ PIOGGE /PP(10000)
      COMMON/ NAME   /NOME
      COMMON/ NEVE   /TS
      COMMON/ PARNAN /BASE,WM,B,D,D1,AREA(100),NVET,ICOD
     .                ,SOL,SOL1,PERC,CESP,ASC,QC,PC
     .                ,W0B
      COMMON/ DATAST /ISW(30),IORD(30)
      COMMON/ PAR    /DT,DDT,TASS,CLF,CLE,T0,MAXMEM,CSNOW,CSA,CSG
      COMMON/ STR    /ICOL1,IPOS
*
*
*
*---------> LEGGE SU FILE AFFDEF.txt   I VALORI DI CALIBRAZIONE
*               (( it reads on file affdef. dat the values of calibrazione  ))
*
      DO 5 I=1,N
        XR(I)=0.
        XRB(I) = 0
5     CONTINUE
      REWIND 1
      REWIND 2
      IPO=0
C
C  INITALISE CATCHMENT MB VARIABLES
C
      QC = 0.0
      PC = 0.0
      ASC=0.0
C
C  FILE CATCHDEF.txt
C
      REWIND(1)
      READ(1,*)
      READ(1,*)
10    READ(1,'(A3,2X,A,2X,I1,2X,I1,2X,F5.0)',END=30,ERR=8000)
     .IISEZ,TYP,NINGQ,ISERB,BASE
      READ(1,*)
      READ(1,*)
      IF (TYP.EQ.'q') TYP = 'Q'
      IF (TYP.EQ.'p') TYP = 'P'
      IF (TYP.EQ.'t') TYP = 'T'
      IF (TYP.EQ.'o') TYP = 'O'
c
      IPO=IPO+1
      IF(IISEZ.NE.NCODI) GO TO 10
C
C FILE STATOINI.txt
C
20    READ(2,'(A3,3F10.2)',ERR=8001) IISEZ1,W0,W0N,YINIZ
C     CALL MAIUSCOL(IISEZ1)
      IF(IISEZ1.NE.NCODI) GO TO 20
C
      IF(TYP.EQ.'Q'.OR.TYP.EQ.'T') GO TO 40
      IF(W0.GT.9000.) GO TO 35
      GO TO 40
C
30    WRITE(*,*) ' CANNOT FIND CATCHMENT ',NCODI
      STOP
35    WRITE(*,*) ' LO STATO DEL SISTEMA NON E'' AGGIORNATO PER TUTTI I
     .BACINI'
      STOP
40    CONTINUE
*
*---------> OPEN THE BASIN INPUT FILE 
*
      NOME1=NCODI//'.txt'
      OPEN(7,FILE=NOME1,STATUS='OLD',IOSTAT=IOCHECK)
      IF(IOCHECK.NE.0) GO TO 9999
      WRITE(*,*)  ' CATCHMENT :',NOME1
      WRITE(13,*)
      WRITE(13,*) ' CATCHMENT:',NOME1
      READ(7,*)
      READ(7,'(A3)') IISEZ1
      IF (IISEZ1.NE.NCODI) GOTO 8002
      READ(7,*)
      READ(7,*) NINQ
c      IF(NINGQ.GT.0) THEN     
c        READ(7,*)
c        READ(7,*) (C(I+2,1),DC(I+2,1),DX(I+2,1),I=1,NINGQ)
c      ENDIF
      READ(7,*)
      READ(7,*) NCB
      READ(7,*)
      READ(7,*) CC,DCC,DXX
C
C  COVER BLOCK DATA, FOR CATCHMENTS WITH RAINFALL INPUT
C
      AS = 0.0
      IF(TYP.EQ.'P'.OR.TYP.EQ.'T') THEN
        READ(7,*)    
        READ(7,*) C(1,1),DC(1,1),DX(1,1)
        DX(1,1)=DX(1,1)*1000.
        READ(7,*) C(2,1),DC(2,1),DX(2,1)
        DX(2,1) =DX(2,1)*1000.
C
C  CHINESE MODEL INPUT DATA
C
      READ(7,*)
      READ(7,*) WM,SOL,SOL1,B,D,D1,PERC,CESP,Ts

C
C   RAINFALL INPUT ON COVER BLOCK   ((rainfall input on cover block))
C
          READ(7,*) AREA(1)

C  READ IN BASIN PET DATA AND RAINFALL DATA HERE
       NOME1 = 'hisrain25.txt'
       OPEN (UNIT=4,FILE=NOME1,IOSTAT=IOCHECK)
c       READ(4,*)
c       IF (IOCHECK.NE.0) GOTO 9999
c       DO 64 J=1,NS
c         READ(4,*,ERR=8003) RDUM
c 64    CONTINUE
       DO 65 J=1,N
          READ(4,*,ERR=8003) PP(J)
 65    CONTINUE
       CLOSE(4)
       NOME1='hispan25.txt'
       OPEN (UNIT=4,FILE=NOME1,IOSTAT=IOCHECK)
c       READ(4,*)
c       IF (IOCHECK.NE.0) GOTO 9999
c       DO 66 J=1,NS
c         READ(4,*,ERR=8004) RDUM
c 66    CONTINUE
       DO 68 J=1,N
          READ(4,*,ERR=8004) PET(J)
 68    CONTINUE
       CLOSE(4)
C
C  COVER BLOCK LOOP
C       
      DO 75 IC=1,NCB
C
C   START OF MAINLOOP
C
c  W0B = COVER BLOCK WATER CONTENT ( W0 = GLOBAL INITIAL VALUE)
c
      W0B = W0
      IF(W0B.GT.WM) W0B=WM
      IF (SOL.GE.WM) SOL=WM-1.
      TOTQ=W0B 
C
C  OPEN AET FILE
C
        FILNAM=NCODI//'.aet'
        OPEN (UNIT=12,FILE=FILNAM)
c        DO 84 J=1,5
c 84        WRITE(12,*)
C
C   TIMESTEP LOOP  CALL SOIL/ET MODEL
C
        DO 85 IHR=1,N
C  
C
C  CALL SOIL / ET ACCOUNTING MODEL (CHINESE)
C
          CALL CHINESE(IC,IHR)
	open (19,file='wob.txt')
	write(19,*) W0B
C
C  WRITE AET TO FILE
C
          WRITE(12,83) AET(IHR)
 83       FORMAT (F12.5)
 85     CONTINUE
C
C  RUNOFF CONVOLUTION
C
        NPAR=1
        IUL=1
*
        ICB=IC
        CALL CONV (XRB)
        DO 120 I=1,N
          XR(I) = XR(I)+XRB(I)
 120    CONTINUE
 75   CONTINUE 
       DO 125 J=1,2+NINGQ
          DX(J,1)=DX(J,1)*1000.
125    CONTINUE
       ENDIF
C
C  EXTERNAL INFLOW, ONLY 1 ALLOWED IN THIS VERSION
C
       IF(NINGQ.GT.0) THEN
C
C  OPEN FLOW FILE
C
        OPEN (5,FILE='basin/'//NCODI//'.inf',ERR=9999)
        NPAR=2
        READ(5,*)
        READ(5,*)
        READ(5,*)
        READ(5,*)
        READ(5,*)
        DO 137 J=1,NS
137       READ(5,*,END=8007,ERR=8006) DUM,DUM         
        DO 141 J=1,N
          READ(5,*,END=8007,ERR=8006) DUM,PP(J)
 141    CONTINUE
        NPAR=NPAR+1
 142    FORMAT(f8.1,f9.3)
        ICB=1
        CLOSE(5)
        CALL CONV (PP)
        DO 135 J=1,N
          XR(J)=XR(J)+PP(J)
135     CONTINUE
150     CONTINUE
      ENDIF
      GO TO 9998
9999  WRITE(*,*) ' NON ESISTE IL FILE', NOME1
      STOP
9998  CONTINUE
      ICB=1
      DX(2,1)=DXX*1000.
      DC(2,1)=DCC
      C(2,1) =CC
*
      CLOSE(12)
      RETURN
8000  WRITE(*,*) ' FORMAT ERROR IN "affdef.txt"'
      STOP
8001  WRITE(*,*) ' FORMAT ERROR IN "STATOINI.txt"'
      STOP
8002  WRITE(*,*) ' FORMAT ERROR IN "',NOME1,'"'
      STOP
8003  WRITE(*,*) ' FORMAT ERROR IN RAINFALL FILE ',FILNAM
      STOP
8004  WRITE(*,*) ' FORMAT ERROR IN PET FILE ',NOME1
      STOP
8005  WRITE(*,*) ' ERRORE NEL FORMATO DEL FILE "TEMPMED.TPR"'
      STOP
8006  WRITE(*,*) ' FORMAT ERROR IN INFLOW FILE '
      STOP
8007  WRITE(*,*) ' END OF FILE REACHED IN INFLOW FILE'
      STOP
      END
*
************************************************************************
*
*     SUBROUTINE CONV (Y)
*
*     ESEGUE LA CONVOLUZIONE DELL'IDROGRAMMA UNITARIO CON LE PORTATE
*     CALCOLATE
*  ((it performs the convoluzione dell' unitary idrogramma with the courses))
*
*
************************************************************************
*
      SUBROUTINE CONV (Y)
      REAL*8 U,DX,C,DC,S,DT
	REAL*4 Z(10000),Y(10000)
*
      COMMON / IUHPAR / U(1000),NM
      COMMON / CONVL  / DX(10,100),C(10,100),DC(10,100),NPAR,IUL,ICB
      COMMON/  NUMBER / N,NT,NS,ISERB
      COMMON/ PAR    /DT,DDT,TASS,CLF,CLE,T0,MAXMEM,CSNOW,CSA,CSG
*
*---------> CALCOLA LE ORDINATE DELL'IDROGRAMMA UNITARIO
*              ((it calculates the ordinates of the unitary idrogrimma))
*
      IF(IUL.EQ.1) CALL UXLAT
      IF(IUL.NE.1) CALL UX
      WRITE(13,1000) (U(J),J=1,NM)
*
*---------> ESEGUE LA CONVOLUZIONE PER LA PROPAGAZIONE
*           IN ALVEO O LUNGO I VERSANTI

*  ((it performs the convoluzione for the propagation in riverbed or long the sides))
*
      DO 40 I=1,N
        S=0.
        DO 30 J=1,NM
          IND=I-J+1
          IF(IND) 10,10,20
10        IND=1
20        S=S+U(J)*Y(IND)
30      CONTINUE
        Z(I)=S
40    CONTINUE
      IF (NPAR.GE.2) GO TO 90
*
*---------> CALCOLA LE ORDINATE DELL'IDROGRAMMA UNITARIO LATERALE
*                ((it calculates the ordinates of the unitary idrogramma side))
*
      NPAR=2
5     CALL UXLAT
*
      WRITE(13,1000) (U(J),J=1,NM)
*
*---------> ESEGUE LA CONVOLUZIONE PER GLI AFFLUSSI LATERALI
*               ((it performs the convoluzione for the side influxes))
*
      DO 80 I=1,N
        S=0.
        DO 70 J=1,NM
          IND=I-J+1
          IF(IND) 50,50,60
50        IND=1
60        S=S+U(J)*Z(IND)
70      CONTINUE
        Y(I)=S
80    CONTINUE
      GO TO 110
*
*---------> RIASSEGNA IL RISULTATO AL VETTORE Y IN
*           INGRESSO
*                 ((raissegna the result to the carrier y in entrance))
*
90    DO 100 I=1,N
      Y(I)=Z(I)
100   CONTINUE
110   IUL=0
*
      RETURN
*
1000  FORMAT (' IDROGRAMMA UNITARIO :',8F6.4,/,10(22X,8F6.4,/))
*
      END
***** SUBROUTINE UX ****************************************************
*                                                                      *
*     COMPUTES THE ORDINATES OF THE DISCRETIZED PARABOLIC RESPONSE     *
*     FUNCTION FOR UPSTREAM CONCENTRATED INFLOWS 
*
*        ((computes the ordinates of the discretized parabolic
*              answers function for upstream concentrated inflows))
*
*
*
*                                                                      *
********************************************************** START *******
      SUBROUTINE UX
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*4 DDT,TASS,CLF,CLE,T0,CSNOW,CSA,CSG
      COMMON / CONVL  / XX(10,100),CONV(10,100),DIFF(10,100),ICOD,IUL
     $,IUX
      COMMON / IUHPAR / U(1000),NM
      COMMON/ PAR    /DT,DDT,TASS,CLF,CLE,T0,MAXMEM,CSNOW,CSA,CSG
      DATA UMAX/.999D0/,XMAX/88.D0/
      MEM=MAXMEM-1
      DELTAX=XX(ICOD,IUX)
      C=CONV(ICOD,IUX)
      D=DIFF(ICOD,IUX)
      COST=DELTAX*C/D
      IF (COST.GT.XMAX) COST=XMAX
      COST=DEXP(COST)
      DSC=D/C/C
      V=DELTAX/C
      U(1)=0.
      DO 10 I=1,MEM
        T=I*DT
        CXT=C*T
        DEN=DSQRT(2.*D*T)
        ARG=-(DELTAX-CXT)/DEN
        X1=CDNOR(ARG)
        ARG=-(DELTAX+CXT)/DEN
        X2=CDNOR(ARG)*COST
        U1=X1+X2
        U2=X1-X2
        U(I+1)=T*U1-V*U2
        U(I)=(U(I+1)-U(I))/DT
        IF(I.GT.1.AND.U(I).GT.UMAX) GO TO 20
10    CONTINUE
      WRITE(*,*)' ===> ERRORE NEI PARAMETRI DELL'' IDROGRAMMA UNITARIO'
      WRITE(*,*)' ===> MEMORIA =  ',MEM
      WRITE(*,*)' ===> C = ',C,' D = ',D,' DX = ',DELTAX
      WRITE(13,*)' ===> ','U'
      WRITE(13,*)' ===> ',(U(J),J=1,MAXMEM)
      STOP
20    S=0.0
      NM=I
      DO 30 I=2,NM
        J=NM-I+2
        U(J)=(U(J)-U(J-1))
        IF(U(J).LT.0.D0) U(J)=0.D0
        S=S+U(J)
30    CONTINUE
      IF(U(1).LT.0.D0) U(1)=0.D0
      S=S+U(1)
      DO 40 I=1,NM
        U(I)=U(I)/S
40    CONTINUE
      RETURN
      END
***** SUBROUTINE UXLAT *************************************************
*                                                                      *
*   computes the ordinates of the discretized parabolic response
*         function for lateral distributed inflows
*
*
*                        
*                                                                      *
********************************************************** START *******
      SUBROUTINE UXLAT
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*4 DDT,TASS,CLF,CLE,T0,CSNOW,CSA,CSG
      COMMON / CONVL  / XX(10,100),CONV(10,100),DIFF(10,100),ICOD,IUL
     $,ICB
      COMMON / IUHPAR / ULAT(1000),NM
      COMMON/ PAR    /DT,DDT,TASS,CLF,CLE,T0,MAXMEM,CSNOW,CSA,CSG
      EQUIVALENCE (C,E)
      DATA UMAX/.999D0/,XMAX/88.D0/,
     +     PIG/3.141592654D0/,SQPIG/1.772453851D0/,SQ2/1.414213562D0/
      MEM=MAXMEM-1
      DELTAX=XX(ICOD,ICB)
      C=CONV(ICOD,ICB)
      D=DIFF(ICOD,ICB)
      COST=DELTAX*C/D
      IF (COST.GT.XMAX) COST=XMAX
      COST=DEXP(COST)
      DSC=D/C/C
      V=DELTAX/C
      V2=V*V
      ULAT(1)=0.
      DO 10 I=1,MEM
        T=I*DT
        T2=T*T
        CXT=C*T
        DXT=D*T
        D1=DELTAX-CXT
        D2=DELTAX+CXT
        SDXT=DSQRT(DXT)
        DEN=SQ2*SDXT
        DEN1=4.*DXT
        A1=T2+V2
        A2=2.*V*(T-DSC)
        A3=D1*D1/DEN1
        IF (A3.GT.XMAX) A3=XMAX
        A3=2.*V2*SDXT/(DEXP(A3)*DELTAX*SQPIG)
        ARG=-D1/DEN
        X1=CDNOR(ARG)
        ARG=-D2/DEN
        X2=CDNOR(ARG)*COST
        U1=X1+X2
        U2=X1-X2
	
        ULAT(I+1)=E*(T2-A1*U1+A2*U2+A3)/(2.*DELTAX)
        ULAT(I)=(ULAT(I+1)-ULAT(I))/DT
        IF(I.GT.1.AND.ULAT(I).GT.UMAX) GO TO 20
10    CONTINUE
      WRITE(*,*)' ===> ERRORE NEI PARAMETRI DELL'' IDROGRAMMA UNITARIO'
      WRITE(*,*)' ===> MEMORIA =  ',MEM
      WRITE(*,*)' ===> C = ',C,' D = ',D,' DX = ',DELTAX
      WRITE(13,*)' ===> ','ULAT'
      WRITE(13,*)' ===> ',(ULAT(J),J=1,MAXMEM)
      STOP
20    SL=0.0
      NM=I
      DO 30 I=2,NM
        J=NM-I+2
        ULAT(J)=(ULAT(J)-ULAT(J-1))
        IF(ULAT(J).LT.0.D0) ULAT(J)=0.D0
        SL=SL+ULAT(J)
30    CONTINUE
      IF(ULAT(1).LT.0.D0) ULAT(1)=0.D0
      SL=SL+ULAT(1)
      DO 40 I=1,NM
        ULAT(I)=ULAT(I)/SL
40    CONTINUE
      RETURN
      END
************************************************************************
      SUBROUTINE DEFRAMO
*
      REAL*8 DX,C,DC
      CHARACTER*3 NCOD,NCODI,NISEZ
      CHARACTER*1 TYP
*
      COMMON/ RUNF   /RUN(10000),XR(10000),XRB(10000),TYP
      COMMON/ NUMBER /N,NT,NS,ISERB
      COMMON/ STO    /STORAGE(10000,2),KSTOR,ICOL(30),IPOSOLD(30)
     .               ,NINOLD(30)
      COMMON/ INTERV /NIN,IRAM
      COMMON/ ALB    /NRAM,NEL(30),NCOD(30),IMEM(30)
      COMMON/ INP    /NISEZ(30),NSZ
      COMMON/ UP     /NCODI,W0,W0N,YINIZ
      COMMON/ CONVL  /DX(10,100),C(10,100),DC(10,100),NPAR,IUL,ICB
      COMMON/ DEFR   /NELINI,NELFIN,ISEZIONE
      COMMON/ DATAST /ISW(30),IORD(30)
      COMMON/ STR    /ICOL1,IPOS
*
      IPOS2=IPOS
          DO 10 K=1,N
          RUN(K)=0.
10        CONTINUE
        DO 50 J=NELINI,NELFIN
        IPOS2=NELFIN-J+NELINI
        NCODI=NCOD(IPOS2)
        ISWITCH=ISW((IPOS2))
        IF(ISWITCH.EQ.0) GO TO 50
*
        CALL BACINO
*
        CALL UPDATE
*
          DO 20 K=1,N
          RUN(K)=RUN(K)+XR(K)
20        CONTINUE
C      IF(ISERB.NE.0) CALL SERB(RUN)
      DO 30 I=1,NSZ
      IF(NCODI.NE.NISEZ(I)) GO TO 30
*
      OPEN(11,FILE=NCODI//'flow2190.txt')
      CALL PRNTC(I,RUN)
      CLOSE(11)
*
      NISEZ(I)=NISEZ(NSZ)
      NSZ=NSZ-1
        IF(NSZ.EQ.0) STOP
      GO TO 40
30    CONTINUE
40    CONTINUE
      NPAR=2
      ICB=1
*
      CALL CONV (RUN)
*
50    CONTINUE
*
       RETURN
       END
*
************************************************************************
*
      SUBROUTINE CHINESE(ICB,I)
*
*    LA SUBROUTINE NAN CALCOLA IL DEFLUSSO USCENTE DA UN SOTTOBACINO
*    DATE LE PIOGGE TENENDO CONTO ANCHE DELLA EVENTUALE NEVE
*
*      (( the subroutine nan calculates the ending deflusso from a sottobacino
*        given the rains holding account also of the eventual snow))
*
*
************************************************************************
*
      REAL*8 DT
      INTEGER*2 YY,MM,DD,HH,YYF,MMF,DDF,HHF
      CHARACTER*3 NCODI
      CHARACTER*1 TYP
*
      COMMON/ RUNF   /RUN(10000),XR(10000),XRB(10000),TYP
      COMMON/ NUMBER /N,NT,NS,ISERB
      COMMON/ EVAPO  /PET(10000),AET(10000)
      COMMON/ UP     /NCODI,W0,W0N,YINIZ
      COMMON/ PIOGGE /PP(10000)
      COMMON/ NEVE   /TS0
      COMMON/ PARNAN /BASE,WM,B,D,D1,AREA(100),NVET,ICOD,SOL,SOL1,PERC,
     $CESP,ASC,QC,PC,W0B
      COMMON/ PAR    /DT,DDT,TASS,CLF,CLE,T00,MAXMEM,CSNOW,CSA,CSG
      COMMON/ DAT    /YY,MM,DD,HH,YYF,MMF,DDF,HHF
*
      B1=1.D0+B
      BM=1.D0/B1
C
C  SET INITIAL ET VALUE TO POTENTIAL VALUE 
C
      ETX = PET(I)
      PX=PP(I)
      RX=0.
      IF(PX.GT.0.) THEN
        C=(1.D0-W0B/WM)**BM-PX/(WM*B1)
         RX=PX-(WM-W0B)
          IF(C.GT.0.D0) THEN
            RX=RX+WM*C**B1
          ENDIF
         IF(RX.LT.0.) RX=0.
        PX=PX-RX
      ENDIF
      UNSA=(1.-W0B/WM)**B
      SA=1.-UNSA
*
*---------> CALCOLA L'ESAURIMENTO DI FALDA
*              ((it calculates the exhaustion of layer))
*
***************************************
      ETX=ETX*(SA+UNSA*W0B/WM)
***************************************
      AFFNET=PX-ETX
      W0B=W0B+AFFNET
      IF(W0B.LT.0.) THEN
          W0B=0.
          DRAIN=D1*W0B/WM
        ELSE
          DRAIN=D1*W0B/WM
          IF(W0B.GT.SOL) THEN
              DRAIN=DRAIN+(D-D1)*((W0B-SOL)/(WM-SOL))**CESP
          ENDIF

           IF(W0B.GE.DRAIN) THEN
              W0B=W0B-DRAIN
                 ELSE
              DRAIN=W0B
            W0B=0.
           ENDIF
*
*---------> CALCOLA IL PERCOLATO PROFONDO
*             ((calculates the deep percolato))
*
          IF(W0B.GT.SOL1)then
            W0B=W0B-PERC*(W0B-SOL1)/(WM-SOL1)
          endif
          IF(W0B.LT.0.) W0B=0.
      ENDIF
*
*---------> CALCOLA IL RUNOFF E LO TRASFERISCE ALLA MATRICE XRB
*              ((calculates the runoff and the transfer to the matrix xrb))
*
        IF(W0B.GT.WM) W0B=WM
        RX=RX+DRAIN
c
c  ALLOCATE AET VALUE TO ARRAY AET AND RUNOFF VALUE TO ARRAY XRB
c
        AET(I)=ETX
        XRB(I)=BASE+RX*AREA(ICB)*1000./DT
      RETURN
      END
****************************************************************************
                           SUBROUTINE READ
*                                                                      *
*     LA SUBROUTINE READ LEGGE I COEFFICIENTI PER IL CALCOLO DI ETP,   *
*     I PARAMETRI DEL MODELLO NEVE, I DATI DI PIOGGIA E I DATI DI      *
*     TEMPERATURA RELATIVI AL PERIODO SCELTO, L'OPZIONE SERBATOI. 
*      ((the subroutine read reads the coefficienti for the calculation
*            of etp, the parameters of the I model snow, the data of rain
*         and the relevant data of temperature to the choice period,))
*         l' option tanks.  
*    
*                                                                      *
************************************************************************
*
      REAL*8 DT
      INTEGER*2 YY,MM,DD,HH,YYF,MMF,DDF,HHF
      INTEGER*4 TIN,TFIN,TIME
      CHARACTER*3 NCOD,ISEZ,CODI(30),NISEZ
      CHARACTER TYP
*
      COMMON/ NUMBER /N,NT,NS,ISERB
      COMMON/ ALB    /NRAM,NEL(30),NCOD(30),IMEM(30)
      COMMON/ DAT    /YY,MM,DD,HH,YYF,MMF,DDF,HHF
      COMMON/ EVAPO  /PET(10000),AET(10000)
      COMMON/ DATAST /ISW(30),IORD(30)
      COMMON/ PAR    /DT,DDT,TASS,CLF,CLE,T0,MAXMEM,CSNOW,CSA,CSG
      COMMON/INP/NISEZ(30),NSZ
*
*
*---------> APRE IL FILE ALBERO.txt E LEGGE LA DESCRIZIONE DELL'ALBERO
*   ((opens the file albero.txt and law the description of the albero))
*
      OPEN(1,FILE='leggi.txt',STATUS='OLD')
      READ(1,*,ERR=8005) NSZ
      READ(1,'(30A3)',ERR=8005) (NISEZ(I),I=1,NSZ)
      CLOSE (1)
C
C  TREE STRUCTURE
C
      OPEN(1,FILE='albero.txt',STATUS='OLD',IOSTAT=IOCHECK)
      IF(IOCHECK.NE.0) GO TO 9990
      READ(1,*,ERR=8000) NRAM
      READ(1,*,ERR=8000) (NEL(I),I=1,NRAM)
      READ(1,'(16A3)',ERR=8000) (NCOD(I),I=1,NEL(NRAM))
      READ(1,*,ERR=8000) (IMEM(I),I=1,NRAM)
      CLOSE (1)
*
*---------> APRE IL FILE PARAM.txt E LEGGE IL GRADIENTE DI TEMPERATURA E I
*           COEFFICIENTI PER IL CALCOLO DELL'ETP
*           opens the file param.txt and read the gradient of temperature 
*                  and the coefficient for the calculation dell'etp
*
      OPEN(1,FILE='param.txt',STATUS='OLD',IOSTAT=IOCHECK)
      IF(IOCHECK.NE.0) GO TO 9995
      READ(1,*,ERR=8001) DT,DDT,TASS,CLF,CLE,T0,MAXMEM,CSNOW,CSA,CSG
      CLOSE (1)
C
      OPEN(13,FILE='idrogram.txt',STATUS='UNKNOWN',IOSTAT=IOCHECK)
C
      OPEN(2,FILE='statoini.txt',STATUS='OLD',IOSTAT=IOCHECK)
      IF(IOCHECK.NE.0) GO TO 9999
c      READ(2,*,ERR=8002) OY,OM,OD,OH
*
*---------> APRE I FILES PIOGGE.TPR E TEMPER.TPR CONTROLLA SE I FILES
*           CONTENGONO LE PIOGGE AI PLUVIOMETRI O GLI AFFLUSSI RAGG. E
*           LEGGE DATA INIZIO, DATA FINE, NUMERO INTERVALLI DI TEMPO, PIOGGE
*   opens the files rains.tpr and temper.tpr monitors if the files  contain the
*     rains to the pluviometri or the influxes ragg. and law given beginning,
*            given fine, I number intervals of time, rains))
*
*
      OPEN(1,FILE='sim.txt',STATUS='OLD')
      READ(1,*,ERR=8003) YY,MM,DD
      READ(1,*,ERR=8003) YYF,MMF,DDF
      CLOSE(1)
*
*
C      DSTA=TIME(YY,MM,DD,HH)-TIME(OY,OM,OD,OH)
C      IF(YY.NE.OY) DSTA=DSTA+TIME(OY,12,31,24)
C      IF(DSTA.NE.1) GO TO 9993
C
      OPEN(9,FILE='w0w0nyin.txt',STATUS='UNKNOWN')
      WRITE(9,'(4I3)') YYF,MMF,DDF,HHF
*
*---------> APRE IL FILE TEMPER.TPR E LEGGE DATA INIZIO, DATA FINE,
*           NUMERO INTERVALLI DI TEMPO, TEMPERATURE
*   ((opens the file temper. tpr and law given beginning, given fine, I number intervals
*         of time, temperatures))
*  SET NS TO DAYS ELAPSED SINCE START DATE
*
*
      HH=0
      HHF=0
      TIN = TIME(YY,MM,DD,HH)
      TFIN = TIME(YYF,MMF,DDF,HHF)
      N=4200
      NS=0
20    CONTINUE
*
*---------> APRE IL FILE STAZILIV.COD
*               ((opens the file staziltv. cod))
*
      OPEN(1,FILE='affdef.txt',STATUS='OLD',IOSTAT=IOCHECK)
      IF(IOCHECK.NE.0) GO TO 9998
      READ(1,*) NBAC
      READ(1,*)
*
*---> COSTRUISCE IL NOME E APRE IL FILE STATO INIZIALE DEL BACINO
*     ((builds the name and opens the file with initial of the catchment))
*
      DO 30 I=1,30
        READ(1,'(A3,1X,A)',END=40,ERR=8004) CODI(I),TYP
        ISW(I)=1
        IF(TYP.EQ.'O'.OR.TYP.EQ.'o') ISW(I)=0
        IORD(I)=I
30    CONTINUE
40    CONTINUE
      IF(NEL(NRAM).NE.NBAC) GO TO 9992
      DO 60 I=1,NBAC
        DO 50 J=I,NBAC
          IF(NCOD(I).NE.CODI(J)) GO TO 50
          ISEZ=CODI(J)
          CODI(J)=CODI(I)
          CODI(I)=ISEZ
          IO=(J)
          IORD(J)=IORD(I)
          IORD(I)=IO
          GO TO 60
50      CONTINUE
        WRITE(*,*) ' NON ESISTE IL CODICE ',NCOD(I)
	  STOP
60    CONTINUE
*
      GO TO 9994
*
8000  WRITE(*,*) ' ERRORE NEL FORMATO DEL FILE "ALBERO.txt"'
      STOP
8001  WRITE(*,*) ' ERRORE NEL FORMATO DEL FILE "PARAM.txt"'
      STOP
8002  WRITE(*,*) ' ERRORE NEL FORMATO DEL FILE "STATOINI.txt"'
      STOP
8003  WRITE(*,*) ' ERRORE NEL FORMATO DEL FILE "SIM.txt"'
      STOP
8004  WRITE(*,*) ' ERRORE NEL FORMATO DEL FILE "AFFDEF.txt"'
      STOP
8005  WRITE(*,*) ' ERRORE NEL FORMATO DEL FILE "LEGGI.txt"'
      STOP
9991  WRITE(*,*) ' NON ESISTE IL FILE PORTAT.TPR'
      STOP
9888  WRITE(*,*) ' NON ESISTE IL FILE TEMPMED.TPR'
      STOP
9990  WRITE(*,*) ' MANCA IL FILE "ALBERO.txt "'
      STOP
9992  WRITE(*,*) 'IL NUMERO DI BACINI CONTENUTO NEL FILE "AFFDEF.txt" N
     .ON '
      WRITE(*,*)' CORRISPONDE A QUELLO CONTENUTO NELLO SCHEMA AD ALBERO'
      STOP
9993  WRITE(*,*) ' LO STATO DEL SISTEMA NON E'' AGGIORNATO'
      STOP
9995  WRITE(*,*) ' NON ESISTE IL FILE PARAM.txt'
      STOP
9996  WRITE(*,*) ' NON ESISTE IL FILE PIOGGE.TPR'
      STOP
9997  WRITE(*,*) ' NON ESISTE IL FILE TEMPER.TPR'
      STOP
9998  WRITE(*,*) ' NON ESISTE IL FILE AFFDEF.txt'
      STOP
9999  WRITE(*,*) ' NON ESISTE IL FILE STATOINI.txt'
      STOP
9994  CONTINUE
      RETURN
      END
************************************************************************
      SUBROUTINE STOR
      CHARACTER*3 NCOD,NISEZ
*
      COMMON/ NUMBER /N,NT,NS,ISERB
      COMMON/ STO    /STORAGE(10000,2),KSTOR,ICOL(30),IPOSOLD(30)
     .               ,NINOLD(30)
      COMMON/ ALB    /NRAM,NEL(30),NCOD(30),IMEM(30)
      COMMON/ INP    /NISEZ(30),NSZ
      COMMON/ PRINT  /RUNOFF(10000)
      COMMON/ STR    /ICOL1,IPOS
      COMMON/ INTERV /NIN,IRAM
*
       KSTOR=KSTOR+1
       IF(IPOS.NE.ICOL(KSTOR)) GO TO 20
        DO 10 I=1,N
          STORAGE(I,KSTOR)= RUNOFF(I)+STORAGE(I,KSTOR)
          RUNOFF(I)=0.
10      CONTINUE
20      CONTINUE
        DO 30 I=1,N
          STORAGE(I,KSTOR)= RUNOFF(I)
          RUNOFF(I)=0.
30      CONTINUE
      ICOL(KSTOR)=ICOL1
      NINOLD(KSTOR)=NIN
      IPOSOLD(KSTOR)=IPOS
       RETURN
       END
*
************************************************************************
      SUBROUTINE INTER

      INTEGER*2 NNIN(30),NIRAM(30)
      CHARACTER*3 NCOD,ISEZ,NISEZ
*
      COMMON/ INTERV /NIN,IRAM
      COMMON/ ALB    /NRAM,NEL(30),NCOD(30),IMEM(30)
      COMMON/ INP    /NISEZ(30),NSZ
*
      DO 60 JK=1,NSZ
      ISEZ=NISEZ(JK)
      NELINI=1
        DO 20 I=1,NRAM
        NELFIN=NEL(I)
          DO 10 J=NELINI,NELFIN
          IF(NCOD(J).NE.ISEZ) GO TO 10
          NIN=I
          NPOS=IMEM(I)
          GO TO 30
10        CONTINUE
        NELINI=NELFIN+1
20      CONTINUE
        WRITE(*,*) ' NON ESISTE IL CODICE ',ISEZ
        STOP
30    CONTINUE
*
      DO 40 I=NIN+1,NRAM+1
      IF(IMEM(I).LE.NPOS) GO TO 40
      IRAM=I-1
      GO TO 50
40    CONTINUE
50    CONTINUE
*
      IF(IRAM.EQ.0) IRAM=NIN
*
      NIRAM(JK)=IRAM
      NNIN(JK)=NIN
60    CONTINUE
      DO 80 I=1,NSZ
      DO 70 J=I+1,NSZ
      IF(NNIN(I).LT.NNIN(J)) GO TO 70
      ISEZ=NISEZ(I)
      NISEZ(I)=NISEZ(J)
      NISEZ(J)=ISEZ
      NIN=NNIN(I)
      NNIN(I)=NNIN(J)
      NNIN(J)=NIN
70    CONTINUE
80    CONTINUE
      IRAM=NIRAM(1)
      DO 90 I=2,NSZ
      IF(NIRAM(I).GT.IRAM) IRAM=NIRAM(I)
90    CONTINUE
      NIN=NNIN(1)
      RETURN
      END
*
************************************************************************
      SUBROUTINE INPUT

      INTEGER*2 IORD(30)
      CHARACTER*3 NCOD,ISEZ,NISEZ

      COMMON/INTERV/NIN,IRAM
      COMMON/ ALB    /NRAM,NEL(30),NCOD(30),IMEM(30)
      COMMON/INP/NISEZ(30),NSZ
*
      DO 30 I=1,NSZ
        DO 20 J=1,NEL(NRAM)
          IF(NISEZ(I).EQ.NCOD(J)) IORD(I)=J
20      CONTINUE
30    CONTINUE
      DO 50 I=1,NSZ
        DO 40 J=I+1,NSZ
          IF(IORD(I).LT.IORD(J)) GO TO 40
          II=IORD(I)
          IORD(I)=IORD(J)
          IORD(J)=II
          ISEZ=NISEZ(I)
          NISEZ(I)=NISEZ(J)
          NISEZ(J)=ISEZ
40       CONTINUE
50     CONTINUE
      CALL INTER
      RETURN
      END
*
************************************************************************
      SUBROUTINE PRNTC(I,Y)

      REAL*4  Y(10000)
      CHARACTER*3 NISEZ

      COMMON/NUMBER/N,NT,NS,ISERB
      COMMON/INP/NISEZ(30),NSZ

        WRITE(11,*) ' Sim '
c        WRITE(11,*)
c        WRITE(11,*)
c        WRITE(11,*)
c        WRITE(11,*)
        DO 100 K=1,N
          WRITE(11,'(F12.4)') Y(K)
 100    CONTINUE
        IF(NSZ.EQ.0) STOP

        RETURN

        END
************************************************************************
*
      FUNCTION TIME(YY,MO,DD,HH)
*
*     LA FUNCTION TIME TRASFORMA UNA DATA ESPRESSA IN AA MM HH IN UN
*     INTERO CHE RAPPRESENTA IL NUMERO DI ORE DALL'INIZIO DELL'ANNO
*     ALLA DATA STESSA
*
*         ((the function time transforms an express date in aa mm hh in a
*          whole one that represents the number of hours from the beginning 
*              of the year to the same date))
*
************************************************************************
*
      INTEGER*4 TIME
      INTEGER*2 HH,DD,MO,YY,TDD,YEAREL
      INTEGER*2 MONTH(12)
*
      DATA MONTH/31,30,31,31,28,31,30,31,30,31,31,30/
*
      MONTH(5)=28
c      IF(MOD(YY,4).EQ.0) MONTH(5)=29
      TDD=DD-1
      MOM=MO-1
*
      DO 10 I=1,MOM
        TDD=TDD+MONTH(I)
10    CONTINUE
*
      YEAREL = YY-52
      TIME=HH+TDD+YEAREL*365
*
      RETURN
      END
************************************************************************
*
                      SUBROUTINE UPDATE
*
*     LA SUBROUTINE UPDATE SCRIVE SU FILE LO STATO DEL SISTEMA AL
*     DELLA SIMULAZIONE
*
*
*          ((the subroutine update writes on file the state of the system 
*           to the of the simulation))
*
************************************************************************
*
      CHARACTER*3 NCODI
      COMMON/UP/NCODI,W0,W0N,YINIZ
*
      WRITE(9,'(A3,3F10.4)') NCODI,W0,W0N,YINIZ
*
      RETURN
*
      END
*
***** FUNCTION CDNOR ***************************************************
*                                                                      *
*     COMPUTES THE CUMULATIVE STANDARD NORMAL DISTRIBUTION 
*
*            ((computes the cumulative standard normal distribution))
*                                                                      *
********************************************************** START *******
      FUNCTION CDNOR(X)
      IMPLICIT REAL*8 (A-H,O-Z)
      DATA XMIN/-690.D0/
      AX=DABS(X)
      ARG=AX*AX/2.D0
      R=1.0/(1.D0+.2316419D0*AX)
      CDNOR=DLOG(((((1.330274429D0*R-1.821255978D0)*R+1.781477937D0)*
     *R-.356563782D0)*R+.319381530D0)*R*.39894228)
      DIF=CDNOR-ARG
      IF(DIF.LT.XMIN) DIF=XMIN
      CDNOR=DEXP(DIF)
      IF(X.GT.0.D0) CDNOR=1.D0-CDNOR
      RETURN
      END
