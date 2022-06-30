module supercellfold
  use kinds,only : dp,dpc
  use constants,only : maxlen
  implicit none
  integer :: sc1,sc2,sc3
  ! Supercell range A,B,C
  real(kind=dp) :: fsthick
  ! In supercell the width of the Fermi surface window to take into accout states to output.
  character(len=maxlen) :: dirgmnvkqsc
  ! The director to save the gmnvkq file in supercell
  logical :: scread 
  ! read the gmnvkq in supercell
  
  contains
  
  

end module 