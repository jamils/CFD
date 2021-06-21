program readin
  REAL, DIMENSION(53,17) :: x
  REAL, DIMENSION(53,17) :: y
  open(unit=12,file='Project_Files/Grids/Inlet-Grids/Inlet.33x17.grd')
  read(12,*) nzones
  read(12,*) imax, jmax, kmax 

  read(12,*) (((x(i,j),i=1,imax),j=1,jmax),k=1,kmax),  &
             (((y(i,j),i=1,imax),j=1,jmax),k=1,kmax),  &
             (((zztemp,i=1,imax),j=1,jmax),k=1,kmax)

!  ! Read in x-coordinate
!  do k = 1, kmax
!  do j = 1, jmax
!  do i = 1, imax
!  read(12,*) x(i,j)
!  enddo
!  enddo
!  enddo
!  ! Read in y-coordinate
!  do k = 1, kmax
!  do j = 1, jmax
!  do i = 1, imax
!    read(12,*) y(i,j)
!  enddo
!  enddo
!  enddo
  
  do i=1,imax
     do j=1,jmax
        write (*,*) y(i,j)
     enddo
  enddo
  
end program readin
