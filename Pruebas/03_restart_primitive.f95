    subroutine init_restart

    use scalar
    use parameters
    use threevectors
    use funciones

    implicit none
    
! Remember when you use RESTART option, you must assignate a value to "sigma" in TMs test this value is ussually "sigma_0"

    if (DIM == 1) then
       
!:::::::::::OPEN BACK UP FILE:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!    open (unit = 222, action="read", file = "all_variables_2.dat", IOSTAT=ierror)

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


       do i=-6,imax+6

!::::::::: READ BACK UP FILE::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

          read (222,*) psi(i,1,1), phi(i,1,1), Bx (i,1,1), By (i,1,1), Bz (i,1,1), Vx (i,1,1), Vy (i,1,1), Vz (i,1,1), &
                       Ex (i,1,1), Ey (i,1,1), Ez (i,1,1), q  (i,1,1), rho(i,1,1), p  (i,1,1)

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      end do ! for i

    
    else if (DIM ==2) then
       
!:::::::::::OPEN BACK UP FILE:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!           call open_backup

!!$    open (unit = 222, file = "all_variables_2.dat", status='old',            &
!!$         access='sequential', action='read', IOSTAT=ierror)

! form='formatted',     
    
!!$    open (unit = 223, file = "file2a.dat")
    
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

! Remember when you use RESTART option, you must assignate a value to "sigma" in TMs test this value is ussually "sigma_0"

        do i=-6,imax+6
           do j=-6,jmax+6

!::::::::: READ BACK UP FILE::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::                                                                                          

             read (222,*) psi(i,j,1), phi(i,j,1), Bx (i,j,1), By (i,j,1), Bz (i,j,1), Vx (i,j,1), Vy (i,j,1), Vz (i,j,1), &
                          Ex (i,j,1), Ey (i,j,1), Ez (i,j,1), q  (i,j,1), rho(i,j,1), p  (i,j,1)

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::                                                                                          

            end do ! for i
         end do ! for i


!!$         do i=-4,imax+4
!!$            do j=-4,jmax+4
!!$
!!$!::::::::: READ BACK UP FILE::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::                                                                                          
!!$
!!$       write (223,*) i, j, By (i,j,1)
!!$
!!$!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::                                                                                          
!!$
!!$              end do ! for i
!!$          end do ! for i      

!         close(68)

         
  else

        write(*,*) "STOP: init_restart"
        write(*,*) "DIM is not valid"
        stop

  end if
     



!--------FORMAT----------
333  format (E17.5,E17.5,E17.5,E17.5,E17.5,E17.5,E17.5,E17.5,E17.5,E17.5,E17.5,E17.5,E17.5,E17.5)
!--------FORMAT----------
     

   end subroutine init_restart
