program symmetry
implicit none

! gfortran -ffree-line-length-none raman_tensor_phonopy.f90 -o raman_tensor_phonopy

integer :: i,j,k,m,n
integer :: n_type,total_atom,n_mode,n_nonequivalent
real*8 :: a_x,a_y,a_z
real*8 :: b_x,b_y,b_z
real*8 :: c_x,c_y,c_z
integer,allocatable :: n_atom(:),atomic_number(:)
real*8,allocatable :: x_frac(:),y_frac(:),z_frac(:)
real*8,allocatable :: x(:),y(:),z(:),x_tmp(:),y_tmp(:),z_tmp(:)
real*8,allocatable :: dx(:,:),dy(:,:),dz(:,:),displace(:,:),d_max(:)
real*8,allocatable :: frequency(:),species_mass(:),reduced_raman(:),reduced_raman_average(:),cross_section(:)
real*8,allocatable :: occupation(:)
real*8,allocatable :: eps_der(:,:,:,:), raman_tensor(:,:,:)
real*8 :: polarization_incident(3), polarization_scattered(3),intensity_tmp(3)
real*8 :: z_max,z_min,surface_depth,tmp,frequency_cutoff
character*100 :: line,line1,line2,line3
character*2 ::  at_type(100)=" "
character*30 :: char_1,char_2,char_3
character*30 :: filename,axis
character*2,allocatable :: species(:),species_atom(:)
character*30, allocatable :: representation(:)
integer :: option, flag
logical :: symbol_line,too_large,existence
! define the functions
integer :: get_z
real*8 :: atomic_mass
character*2 ::  symbol
real*8,parameter :: constant=0.004824125  ! constant=hcm-1/KbT

frequency_cutoff=0.03
! frequency_cutoff=2.0

! write(*,"(A,F5.1,A)") "Frequency under",frequency_cutoff*33.356407889," cm-1 is for reduced Raman while above is for cross section"

!!!!!! At room temperature 298 K, KbT = 25.7 meV; 1 cm-1 equals to 0.12398 meV 
!!!!!! thus hcm-1/KbT=0.12398/25.7=0.004824125 

inquire(file="input",exist=existence)
if(existence) then
  open(unit=11,file="input")
  read(11,*) polarization_incident(:)
  read(11,*) polarization_scattered(:)
  read(11,*) axis
  close(11) 
else
  write(*,*) "please input the polarization of the incident light:"
  read(*,*) polarization_incident(:)
  write(*,*) "please input the polarization of the scattered light:"
  read(*,*) polarization_scattered(:)
  write(*,*) "please indicate the direction perpendicular to the sample surface: x, y or z"
  read(*,*) axis
  open(unit=11,file="input")
  write(11,"(3F4.1,A)") polarization_incident(:), "    ! the polarization of the incident light"
  write(11,"(3F4.1,A)") polarization_scattered(:),"    ! the polarization of the scattered light"
  write(11,"(A2,A)")      axis, "             ! the surface normal (or out-of-plane) direction, generall it is z"
  close(11)
end if

open(unit=1,file="CONTCAR")
open(unit=2,file="band.yaml")


do i=1,7
  read(1,*)
end do

read(1,*) line

if(line .eq. "Direct" .or. line .eq. "cart") then
  symbol_line=.True.
else
  symbol_line=.False.
end if

rewind(1)

if(symbol_line) then
  do i=1,5
    read(1,*) 
  end do
  read(1,"(A)") line
  line=trim(line)
  line=adjustl(line)
  read(line,*,end=10) at_type(:)
  10 n_type=0
  do i=1,100
    if(at_type(i) .ne. " ") then
      n_type=n_type+1
    end if
  end do
else
  read(1,"(A)") line
  read(line,*,end=100) at_type(:)
  100 n_type=0
  do i=1,100
    if(at_type(i) .ne. " ") then
      n_type=n_type+1
    end if
  end do
end if


allocate(n_atom(n_type),atomic_number(n_type))

rewind(1)


read(1,"(A)") line1
read(1,*) 
read(1,*) a_x,a_y,a_z
read(1,*) b_x,b_y,b_z
read(1,*) c_x,c_y,c_z 
if(symbol_line) then
  read(1,"(A)") line3
else
  line3=line1
end if
read(1,*) n_atom(:)
read(1,*) 

total_atom=0
do i=1,n_type
  total_atom=total_atom+n_atom(i)
end do

allocate(x_frac(total_atom),y_frac(total_atom),z_frac(total_atom),species_mass(total_atom))
allocate(x(total_atom),y(total_atom),z(total_atom))
allocate(x_tmp(total_atom),y_tmp(total_atom),z_tmp(total_atom))


do i=1,total_atom
  read(1,*) x_frac(i),y_frac(i),z_frac(i)
end do

!!determine the mass for each atom

k=0
do i=1,n_type
  do j=1,n_atom(i)
    k=k+1
    atomic_number(i)=get_z(at_type(i))
    species_mass(k)=atomic_mass(atomic_number(i))
  end do
end do

    
n_mode=total_atom*3

allocate(dx(n_mode,total_atom),dy(n_mode,total_atom),dz(n_mode,total_atom))

allocate(frequency(n_mode))

do    
  read(2,"(A)") line  
  line=trim(adjustl(line))    
  read(line,*,end=201) char_1, char_2
  201 if(char_1 .eq. "-" .and. char_2 .eq. "#") then    
   exit    
  end if    
end do   

backspace(2)

do i=1,n_mode
  read(2,*) 
  read(2,*) char_1,frequency(i)
!   write(*,*) frequency(i)
  read(2,*) 
  do j=1,total_atom
    read(2,*) 
    read(2,*) char_1,char_2,dx(i,j)
    read(2,*) char_1,char_2,dy(i,j)
    read(2,*) char_1,char_2,dz(i,j)
  end do
end do

!!!!! normalization of the phonon eigenvectors by the root of the atomic mass to become phonon eigendisplacements 
do i=1,n_mode
  do j=1,total_atom
    dx(i,j)=dx(i,j)/sqrt(species_mass(j))
    dy(i,j)=dy(i,j)/sqrt(species_mass(j))
    dz(i,j)=dz(i,j)/sqrt(species_mass(j))
  end do
end do

close(1)
close(2)


allocate(eps_der(total_atom,3,3,3),raman_tensor(n_mode,3,3))

! open(unit=3,file="eps_der")
! 
! read(3,*)
! 
! do i=1,total_atom
!   read(3,*)
!   read(3,*) eps_der(i,1,1,1),eps_der(i,1,1,2),eps_der(i,1,1,3), eps_der(i,2,1,1),eps_der(i,2,1,2),eps_der(i,2,1,3),eps_der(i,3,1,1),eps_der(i,3,1,2),eps_der(i,3,1,3)
!   read(3,*) eps_der(i,1,2,1),eps_der(i,1,2,2),eps_der(i,1,2,3), eps_der(i,2,2,1),eps_der(i,2,2,2),eps_der(i,2,2,3),eps_der(i,3,2,1),eps_der(i,3,2,2),eps_der(i,3,2,3)
!   read(3,*) eps_der(i,1,3,1),eps_der(i,1,3,2),eps_der(i,1,3,3), eps_der(i,2,3,1),eps_der(i,2,3,2),eps_der(i,2,3,3),eps_der(i,3,3,1),eps_der(i,3,3,2),eps_der(i,3,3,3)
! end do
! 
! close(3)


open(unit=3,file="RAMFILE")

do    
  read(3,"(A)") line    
  line=trim(adjustl(line))    
  read(line,*,end=101) char_1,char_2    
  101 if(char_2 .eq. "Name") then    
   exit    
  end if    
end do    
    
do i=1,total_atom
  read(3,*) 
  read(3,*) eps_der(i,1,1,1),eps_der(i,1,1,2),eps_der(i,1,1,3), eps_der(i,2,1,1),eps_der(i,2,1,2),eps_der(i,2,1,3),eps_der(i,3,1,1),eps_der(i,3,1,2),eps_der(i,3,1,3)
  read(3,*) eps_der(i,1,2,1),eps_der(i,1,2,2),eps_der(i,1,2,3), eps_der(i,2,2,1),eps_der(i,2,2,2),eps_der(i,2,2,3),eps_der(i,3,2,1),eps_der(i,3,2,2),eps_der(i,3,2,3)
  read(3,*) eps_der(i,1,3,1),eps_der(i,1,3,2),eps_der(i,1,3,3), eps_der(i,2,3,1),eps_der(i,2,3,2),eps_der(i,2,3,3),eps_der(i,3,3,1),eps_der(i,3,3,2),eps_der(i,3,3,3)
end do

close(3)


do i=1,n_mode
  raman_tensor(i,:,:)=0.0D0
  do j=1,total_atom
    raman_tensor(i,:,:)=raman_tensor(i,:,:)+eps_der(j,1,:,:)*dx(i,j)+eps_der(j,2,:,:)*dy(i,j)+eps_der(j,3,:,:)*dz(i,j)
  end do
end do 


do i=1,n_mode
  do j=1,3
    do k=1,3
      if(abs(raman_tensor(i,j,k))>10.0**(6.0)) then
        too_large= .true.
      end if
    end do
  end do
end do


allocate(reduced_raman(n_mode),reduced_raman_average(n_mode),cross_section(n_mode))
allocate(representation(n_mode),occupation(n_mode))

do i=1,n_mode
  if(frequency(i) < frequency_cutoff ) then  !!! frequency is in the unit of THz
    occupation(i)=0.0  !! frequency is too small to be considered
  else
    occupation(i)=1.0/(exp(frequency(i)*33.356407889*constant)-1.0)
  end if 
end do  


do i=1,n_mode
  intensity_tmp(1)=polarization_incident(1)*raman_tensor(i,1,1)+polarization_incident(2)*raman_tensor(i,2,1)+polarization_incident(3)*raman_tensor(i,3,1)
  intensity_tmp(2)=polarization_incident(1)*raman_tensor(i,1,2)+polarization_incident(2)*raman_tensor(i,2,2)+polarization_incident(3)*raman_tensor(i,3,2)
  intensity_tmp(3)=polarization_incident(1)*raman_tensor(i,1,3)+polarization_incident(2)*raman_tensor(i,2,3)+polarization_incident(3)*raman_tensor(i,3,3)
  reduced_raman(i)=intensity_tmp(1)*polarization_scattered(1)+intensity_tmp(2)*polarization_scattered(2)+intensity_tmp(3)*polarization_scattered(3)
  reduced_raman(i)=reduced_raman(i)**2
  if(frequency(i) < frequency_cutoff ) then  !!! frequency is in the unit of THz
    reduced_raman(i)=reduced_raman(i)  !! frequency is too small to be divided
  else
    reduced_raman(i)=reduced_raman(i)/frequency(i)*(1.0+occupation(i)) ! reduced Raman divided by the frequency and multiplied by (1+n), so more like cross section now
  end if
end do

!!! average all polarizations of the incident and scattered light in the back-scattering set-up 
!!! (the default surface normal direction is the z direction)

if(axis .eq. "z") then
  do i=1,n_mode
    reduced_raman_average(i)=abs(raman_tensor(i,1,1))**2+abs(raman_tensor(i,1,2))**2+abs(raman_tensor(i,2,1))**2+abs(raman_tensor(i,2,2))**2
    if(frequency(i) < frequency_cutoff ) then  !!! frequency is in the unit of THz
      reduced_raman_average(i)=reduced_raman_average(i)  !! frequency is too small to be divided
    else
      reduced_raman_average(i)=reduced_raman_average(i)/frequency(i)*(1.0+occupation(i)) ! reduced Raman divided by the frequency and multiplied by (1+n), so more like cross section now
    end if
  end do
else if(axis .eq. "y") then
  do i=1,n_mode
    reduced_raman_average(i)=abs(raman_tensor(i,1,1))**2+abs(raman_tensor(i,1,3))**2+abs(raman_tensor(i,3,1))**2+abs(raman_tensor(i,3,3))**2
    if(frequency(i) < frequency_cutoff ) then  !!! frequency is in the unit of THz
      reduced_raman_average(i)=reduced_raman_average(i)  !! frequency is too small to be divided
    else
      reduced_raman_average(i)=reduced_raman_average(i)/frequency(i)*(1.0+occupation(i)) ! reduced Raman divided by the frequency and multiplied by (1+n), so more like cross section now
    end if
  end do
else if(axis .eq. "x") then
  do i=1,n_mode
    reduced_raman_average(i)=abs(raman_tensor(i,2,2))**2+abs(raman_tensor(i,2,3))**2+abs(raman_tensor(i,3,2))**2+abs(raman_tensor(i,3,3))**2
    if(frequency(i) < frequency_cutoff ) then  !!! frequency is in the unit of THz
      reduced_raman_average(i)=reduced_raman_average(i)  !! frequency is too small to be divided
    else
      reduced_raman_average(i)=reduced_raman_average(i)/frequency(i)*(1.0+occupation(i)) ! reduced Raman divided by the frequency and multiplied by (1+n), so more like cross section now
    end if
  end do
else 
  write(*,*) "axis input wrong!"
  stop
end if



inquire(file="irreps.yaml",exist=existence)
if(existence) then
  open(unit=2,file="irreps.yaml")
  do
    read(2,*) char_1
    if(char_1 .eq. "normal_modes:") then
      exit
    end if
  end do
  do i=1,n_mode
    read(2,*)
    read(2,*)
    read(2,*) char_2, representation(i)
    read(2,*)
  end do
  close(2)
  
  open(unit=1,file="Raman_tensor")

  do i=1,n_mode
    write(1,"(I5,A,F10.3,A,F10.3,A,A)") i,"   ",frequency(i)," THz ",frequency(i)*33.356407889," cm-1    ",representation(i)
    do j=1,3
      if(too_large) then
        write(1,"(3F16.3)") raman_tensor(i,j,:) 
      else
        write(1,"(3F12.3)") raman_tensor(i,j,:)
      end if
    end do
  end do
  close(1)
  
  open(unit=2,file="Raman_intensity")
  open(unit=3,file="Raman_intensity_polarization_averaged") 
  do i=1,n_mode
    write(2,"(F9.3,F16.3,A,A)") frequency(i)*33.356407889,reduced_raman(i),"       ",representation(i)
    write(3,"(F9.3,F16.3,A,A)") frequency(i)*33.356407889,reduced_raman_average(i),"       ",representation(i)
  end do 
  close(unit=2)
  close(unit=3)
else
  open(unit=1,file="Raman_tensor")

  do i=1,n_mode
    write(1,"(I5,A,F10.3,A,F10.3,A)") i,"   ",frequency(i)," THz ",frequency(i)*33.356407889," cm-1 "
    do j=1,3
      if(too_large) then
        write(1,"(3F16.3)") raman_tensor(i,j,:) 
      else
        write(1,"(3F12.3)") raman_tensor(i,j,:)
      end if
    end do
  end do
  close(1)
  open(unit=2,file="Raman_intensity") 
  open(unit=3,file="Raman_intensity_polarization_averaged") 
  do i=1,n_mode
    write(2,"(F9.3,F16.3)") frequency(i)*33.356407889,reduced_raman(i)
    write(3,"(F9.3,F16.3)") frequency(i)*33.356407889,reduced_raman_average(i)
  end do 
  close(unit=2)
  close(unit=3)
end if
  

end program symmetry

!--------------------------------------------------------------------
  function get_z(Sym)
!given the symbol, gets the atomic number
    integer get_z
    character(len=*), intent(in):: Sym
    integer, parameter  :: NZ=103
    integer ii
    logical :: found=.false.
    character(len=2)    :: SYMBOL  
    character(len=2), parameter :: NAME(NZ) =                       &
     &         (/'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne', &
     &           'Na','Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca', &
     &           'Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn', &
     &           'Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y ','Zr', &
     &           'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn', &
     &           'Sb','Te','I ','Xe','Cs','Ba','La','Ce','Pr','Nd', &
     &           'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb', &
     &           'Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg', &
     &           'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th', &
     &           'Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm', &
     &           'Md','No','Lr'/)
    
    get_z=0
    do ii=1,NZ
       if ( trim(Sym) .eq. trim( NAME(ii) ) ) then
          found=.true.
          get_z=ii
          return
       end if
    enddo
    if (.not.found) then
       write (6,*) "get_z : Could not find Atomic Number for element ",trim(Sym)
       stop
    endif
        
  end function get_z
!----------------------------------------------------------------------
  FUNCTION atomic_mass( Z )
! Given the atomic number, returns the atomic mass 
! written by E. Cruz
! Atomic mass
    real*8 :: atomic_mass
! Atomic number
    integer, intent(in) :: Z       

    integer, parameter  :: NZ=103
    real*8, parameter :: Mass(NZ) =                                  &
     &     (/ 1.00794,  4.002602,  6.941,  9.012182,  10.811,        &
     &       12.0107,  14.0067,  15.9994,  18.9984032,  20.1797,     &
     &       22.98976928,  24.305,  26.9815386,  28.0855, 30.973762, &
     &       32.065,  35.453,  39.948,  39.0983,  40.078,            &
     &       44.955912,  47.867,  50.9415,  51.9961,  54.938045,     &
     &       55.845,  58.933195,  58.6934,  63.546,  65.409,         &
     &       69.723,  72.64,  74.9216,  78.96,  79.904,              &
     &       83.798,  85.4678,  87.62,  88.90585,  91.224,           &
     &       92.90638,  95.94,  98.0,  101.07,  102.9055,            &
     &      106.42,  107.8682,  112.411,  114.818,  118.71,          &
     &      121.76,  127.6,  126.90447,  131.293,  132.9054519,      &
     &      137.327,  138.90547,  140.116,  140.90765,  144.242,     &
     &      145.0,  150.36,  151.964,  157.25,  158.92535,           &
     &      162.5,  164.93032,  167.259,  168.93421,  173.04,        &
     &      174.967,  178.49,  180.94788,  183.84,  186.207,         &
     &      190.23,  192.217,  195.084,  196.966569,  200.59,        &
     &      204.3833,  207.2,  208.9804,  210.0,  210.0,             &
     &      222.0,  223.0, 226.0,  227.0,  232.03806,                &
     &      231.03588,  238.02891,  237.0,  244.0, 243.0,            &
     &      247.0,  247.0,  251.0,  252.0,  257.0,                   &
     &      258.0,  259.0,  262.0 /)  

    IF ( Z.EQ.0 ) THEN
       atomic_mass = 0.0d0
    ELSE IF (ABS(Z).LE.NZ) THEN
       atomic_mass = Mass(ABS(Z))
    ELSE
       WRITE(6,*) 'Atomic_mass: ERROR: No data for Z =', Z
       atomic_mass=0.0d0
    ENDIF

  END function atomic_mass


