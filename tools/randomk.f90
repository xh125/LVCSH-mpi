program main
  implicit none
  
  character(len=20) :: filename = 'Random_K_list.in'
  integer :: fileunit = 20
  
  integer :: nk
  integer :: ipol,ik
  real(kind=8) :: kpoint(3)
  real(kind=8) :: x
  
  write(*,*) "Please input the Number of K-points (<40000)"
  read(*,*) nk
  
  if(nk>40000) then
    write(*,*) "The max number of k-points in QE is 40000, need set a small Number of K-Points"
    write(*,*) "Please input the Number again:"
    read(*,*) nk
  endif
  
  open(unit=fileunit,file=filename)
  write(fileunit,"("K_POINTS crystal")") 
  write(fileunit,"(I5)") nk
  
  call random_seed()
  
  do ik=1,nk
    kpoint = 0.0
    do ipol=1,3
      call random_number(x)
      kpoint(ipol) = x
    enddo
    write(fileunit,"(3F13.6,2X,I5)") (kpoint(ipol),ipol=1,3),1
  enddo
  
  close(fileunit)
  
end program