subroutine open_backup

  use parameters
  use scalar, only : t_init
  
  implicit none


! Remember when you use RESTART option, you must assignate a value to "sigma" in TMs test this value is ussually "sigma_0"
  
  
     if (t_init == 0 .and. RESTART == 1) then

        open (unit = 222, file = "all_variables.dat", status='old',            &
              access='sequential', action='read', IOSTAT=ierror)

        print*, ""
        print*, ""
        print*, "succesfull opened all_variables.dat with RESTART == 1"
        print*, ""
        print*, ""


      
     end if

  
   end subroutine open_backup
