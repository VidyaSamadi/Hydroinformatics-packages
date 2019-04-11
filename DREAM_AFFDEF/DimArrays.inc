      integer NRow,NCol,Nraing,Ntermo,Nsec,NpOss,NRivi
      integer NumStr,SPio,SDT,Nout,npun,NStrcl
      
c     Maximum dimension of the Digital Elevation Model
      PARAMETER(NRow=1000,NCol=1000) 

c     Maximum number of cells whose elevation has not to be increased
      PARAMETER(npun=10)
      
c     Maximum number of raingauges
      PARAMETER(Nraing=4) 
c     Interval of registration of the rainfall data [s]
      PARAMETER(SPio=900)
      
c     Maximum number of temperature stations      
      PARAMETER(Ntermo=1)
        
c     Maximum number of simulation steps
      PARAMETER(NpOss=40000) 
c     Time step of the simulation [s]      
      PARAMETER(SDT=900) 
      
c     Maximum number of links of the river network
      PARAMETER(NRivi=10000)

c     Maximum number of classes of the Strickler roughness for the hillslope
      PARAMETER(NStrcl=3)
      
c     Maximum number of cross-sections where to display the discharge (manual calib) or 
c     to be considered for the automatic calibration    
      PARAMETER(Nsec=1)    

c     Dimension of the discharge output file
      PARAMETER(Nout=NpOss*SPio/SDT)
c     Maximum dimension of arrays used for the raingauge distribution      
      PARAMETER(NumStr=MAX(Nraing,Ntermo))	                                    