      subroutine Warn_error(ncase)

      implicit none

      integer ncase

      write(*,*)
      select case (ncase)

             case(1001) 
               write(*,*) 'DEM not found'
               write(*,*)
               write(*,*) 'The river network determination failed'

             case(1002)                       
               write(*,*) 'Error in removing the pit'
               write(*,*)
               write(*,*) 'The river network determination failed'

             case(1003) 
               write(*,*) 'The maximum number of links of the river 
     .network has been exceeded'
	         write(*,*)
			 write(*,*) 'Only the temporary file of the link 
     .properties has been produced'

             case(1004)
               write(*,*) 'Error in temperature data'

             case(1005)
               write(*,*) 'Error in rainfall data'

             case(1006)
               write(*,*) 'File of link property not found'

             case(1007) 
               write(*,*) 'File of the contributing area not found' 

             case(1008)
               write(*,*) 'File of the modified DEM not found'

             case(1009)
               write(*,*) 'File of the slope pointer not found'

             case(1010)
               write(*,*) 'File of the classes of the hillslope 
     .roughness not found'

             case(1011)
               write(*,*) 'File Affdef.in not found'

             case(1012)
               write(*,*) 'Simulation time step > time step of the 
     .rainfall data'

             case(1013)
               write(*,*) 'The time step of the simulation is wrong'

             case(1014)
               write(*,*) 'The number of discharge data exceeds the 
     .dimension of the file'

	       case(1015)
               write(*,*) 'The number of classes of hillslope roughness  
     .exceeds the dimension of the file'

	       case(1016)
             write(*,*) 'The length of the simulation exceeds 
     .the dimension of the file'

	       case(1017)
             write(*,*) 'Number of cross-section where to display 
     .the river flow exceeds the maximum number allowed'

             case(1018)
             write(*,*) 'The dimension of the DEM exceeds the 
     .maximum allowed'

	       case(1019)
             write(*,*) 'The number of pits whose elevation has not to 
     .be raised exceeds the maximum allowed'

		   case(1020)
             write(*,*) 'The number of raingauges exceeds the maximum 
     .allowed'

		   case(1021)
             write(*,*) 'The number of temperature stations exceeds the 
     .maximum allowed'

	       case(1022)
	       write(*,*) 'Authomatic calibration cannot be performed
     . since more than 1 output cross-section has been specified'
      end select                  

      write(*,*)

      return
      end