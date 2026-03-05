#ifdef HAVE_GRIB
module grib_reader

	use grib_api
	
	implicit none

contains
	subroutine read_grib_file(filename, what_data, met_grid, met_time, met_month)
	
		character(*), intent(in) :: filename, what_data
		real, dimension(:,:,:), intent(inout) :: met_grid
		integer, intent(out) :: met_time, met_month
		
		integer :: grib_unit
		integer :: igrib, iret, total_size, grib_count
		real, dimension(:), allocatable :: single_variable, lats, lons
		character(len=20)  :: name_space, met_date, tmp_str
	  	integer            :: kiter
	  	character(len=256) :: key
	  	character(len=256) :: name_value, lev_value, value, lev_type
	  	character(len=512) :: all
		integer :: t_index, rh_index, o3_index, level, dummy
		character(len=4) :: year
	
		call grib_open_file(grib_unit, filename, 'r')
	
	    call grib_new_from_file(grib_unit,igrib, iret)

		grib_count = 0

     ! valid name_spaces are ls and mars
			name_space = 'mars'
			call grib_keys_iterator_new(igrib,kiter,name_space)
			do 
        		call grib_keys_iterator_next(kiter, iret) 
        
        		if (iret .ne. 1) exit
        
        		call grib_keys_iterator_get_name(kiter,key)
				call grib_get(igrib,trim(key), value)
				
				if (key == "time") then 
					read(value, *) met_time
					met_time = met_time / 100
				endif
                if (key == "date") then
                    read(value,'(a4,i2,i2)') year, met_month, dummy
                endif

			end do
			call grib_keys_iterator_delete(kiter)
				
		if (what_data == "NCEP_GDAS") then 
		
			total_size = 360 * 181
			allocate(single_variable(total_size), lats(total_size), lons(total_size))
		
			name_space = 'ls'
		
			t_index = 26
			rh_index = 47
			o3_index = 60
			
  			do while (iret /= GRIB_END_OF_FILE)

				grib_count=grib_count+1
     
			    call grib_keys_iterator_new(igrib,kiter,name_space)
     			level = -1
     
     			do
        			call grib_keys_iterator_next(kiter, iret) 
        
        			if (iret .ne. 1) exit
        
        			call grib_keys_iterator_get_name(kiter,key)
				    call grib_get(igrib,trim(key), value)

        			if (key == "shortName") then
						name_value = value
					else
						name_value = " "
					endif
					
					if (key == "level") then
						read(value, *) level
					endif

					if (key == "typeOfLevel") then 
						read(value, *) lev_type
					endif

					if (name_value == 't') then
	 	        		call grib_get_data(igrib, lats, lons, single_variable)	 	        
						print*, grib_count, level, lev_type
						if (level >= 10 .and. level <= 1000) then 
							met_grid(:,:,t_index) = reshape(single_variable, (/360,181/) )
							t_index = t_index - 1
						endif
						if (level == 0 .and. trim(lev_type) == "surface" ) then
							met_grid(:,:,51) = reshape(single_variable, (/360,181/) )
						endif
					endif	

					if (name_value == 'r') then
	 	        		call grib_get_data(igrib, lats, lons, single_variable)	 	        
						if (level >= 100 .and. level <= 1000) then 
							met_grid(:,:,rh_index) = reshape(single_variable, (/360,181/) )
							rh_index = rh_index - 1
						endif
						if (level == 2) then
							met_grid(:,:,52) = reshape(single_variable, (/360,181/) )					
						endif
					endif	
					
					if (name_value == "10u") then 
	 	        		call grib_get_data(igrib, lats, lons, single_variable)	 
						met_grid(:,:,48) = reshape(single_variable, (/360,181/) )
					endif							
					if (name_value == "10v") then 
	 	        		call grib_get_data(igrib, lats, lons, single_variable)	 	        
						met_grid(:,:,49) = reshape(single_variable, (/360,181/) )
					endif							
					if (name_value == "sp") then 
	 	        		call grib_get_data(igrib, lats, lons, single_variable)	 	        
						met_grid(:,:,50) = reshape(single_variable, (/360,181/) )
					endif							
					if (name_value == "msl") then 
	 	        		call grib_get_data(igrib, lats, lons, single_variable)	 	        
						met_grid(:,:,53) = reshape(single_variable, (/360,181/) )
					endif							
					if (name_value == "tco3") then 
	 	        		call grib_get_data(igrib, lats, lons, single_variable)	 	        
						met_grid(:,:,54) = reshape(single_variable, (/360,181/) )
					endif							
#ifdef CT_CODE
					if (grib_count >= 121 .and. grib_count <= 126) then

	 	        			call grib_get_data(igrib, lats, lons, single_variable)	 	        
							met_grid(:,:,o3_index) = reshape(single_variable, (/360,181/) )
							o3_index = o3_index - 1

					endif
#endif					
				end do
			    
			    call grib_keys_iterator_delete(kiter)
    			call grib_release(igrib)
    			call grib_new_from_file(grib_unit,igrib, iret)
				
			end do			
		
	
		else if (what_data == "ECMWF") then 
			total_size = 360 * 181
			allocate(single_variable(total_size), lats(total_size), lons(total_size))
		
! variables and levels to be read		
!			(/ "t", "r", "10u", "10v", "sp", "2d", "skt", "lsm", "o3"/)
!			for o3 = (/"18", "22", "25", "30", "34", "38"/)

			name_space = 'ls'
		
			t_index = 20
			rh_index = 40
			o3_index = 46
			
  			do while (iret /= GRIB_END_OF_FILE)

				grib_count=grib_count+1
     
			    call grib_keys_iterator_new(igrib,kiter,name_space)
     
     			do
        			call grib_keys_iterator_next(kiter, iret) 
        
        			if (iret .ne. 1) exit
        
        			call grib_keys_iterator_get_name(kiter,key)
				    call grib_get(igrib,trim(key), value)

        			if (key == "shortName") then
						name_value = value
					else
						name_value = " "
					endif

					if (key == "level") then
						read(value, *) level
					endif

					if (name_value == 't') then
	 	        		call grib_get_data(igrib, lats, lons, single_variable)	 	        
						if (level >= 10) then 
							met_grid(:,:,t_index) = reshape(single_variable, (/360,181/) )
							t_index = t_index - 1
						endif
					endif	

					if (name_value == 'r') then
	 	        		call grib_get_data(igrib, lats, lons, single_variable)	 	        
						if (level >= 10) then 
							met_grid(:,:,rh_index) = reshape(single_variable, (/360,181/) )
							rh_index = rh_index - 1
						endif
					endif	

					if (name_value == "10u") then 
	 	        		call grib_get_data(igrib, lats, lons, single_variable)	 
						met_grid(:,:,41) = reshape(single_variable, (/360,181/) )
					endif							
					if (name_value == "10v") then 
	 	        		call grib_get_data(igrib, lats, lons, single_variable)	 	        
						met_grid(:,:,42) = reshape(single_variable, (/360,181/) )
					endif							
					if (name_value == "sp") then 
	 	        		call grib_get_data(igrib, lats, lons, single_variable)	 	        
						met_grid(:,:,43) = reshape(single_variable, (/360,181/) )
					endif							
					if (name_value == "skt") then 
	 	        		call grib_get_data(igrib, lats, lons, single_variable)	 	        
						met_grid(:,:,44) = reshape(single_variable, (/360,181/) )
					endif							
					if (name_value == "2d") then 
	 	        		call grib_get_data(igrib, lats, lons, single_variable)	 	        
						met_grid(:,:,45) = reshape(single_variable, (/360,181/) )
					endif							

#ifdef CT_CODE
					if (name_value == 'o3') then
						if (level == 18 .or. level == 22 .or. level == 25 .or. &
	    					level == 30 .or. level == 34 .or. level == 38) then  

	 	        			call grib_get_data(igrib, lats, lons, single_variable)	 	        
							met_grid(:,:,o3_index) = reshape(single_variable, (/360,181/) )
							o3_index = o3_index + 1
						endif	
					endif
#endif					
				end do
			    
			    call grib_keys_iterator_delete(kiter)
    			call grib_release(igrib)
    			call grib_new_from_file(grib_unit,igrib, iret)
				
			end do			
		
! there is a single grid, nothing more in each of these products, so no worries, really		
		else if (what_data == "SEA_ICE") then 

			total_size = 720 * 360
			allocate(single_variable(total_size), lats(total_size), lons(total_size))

			call grib_get_data(igrib, lats, lons, single_variable)	 	        
			met_grid(:,:,1) = reshape(single_variable, (/ 720, 360 /) )
		
		else if (what_data == "OZONE") then 

			total_size = 360 * 181
			allocate(single_variable(total_size), lats(total_size), lons(total_size))

			call grib_get_data(igrib, lats, lons, single_variable)	 	        
			met_grid(:,:,1) = reshape(single_variable, (/ 360,181 /) )

		
		endif		
		
		deallocate(single_variable, lats, lons)
			
		call grib_close_file(grib_unit)


	end subroutine read_grib_file


end module grib_reader
#endif
