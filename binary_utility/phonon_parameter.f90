program phonon_parameter
implicit none

!!!!!!!!!!!!!!!!!!!!!!!!
!!! compilation command
! ifort discar.f90 lib/lapack.a lib/blas.a -o discar
!!!!!!!!!!!!!!!!!!!!!
integer :: i,j,k,tmp,m
integer :: n_trans_a,n_trans_b,n_trans_c
integer :: n_type,total_atom,total_atom_super,n_atom_ineq,n_spacegroup,n_crystal
integer,allocatable :: n_atom(:),n_atom_super(:),index_ineq(:),index_all(:),n_atom_ineq_type(:)
integer :: super(3)
real*8 :: lattice
real*8 :: a_x,a_y,a_z
real*8 :: b_x,b_y,b_z
real*8 :: c_x,c_y,c_z
real*8 :: lattice_x,lattice_y,lattice_z
real*8 :: a_x_super,a_y_super,a_z_super
real*8 :: b_x_super,b_y_super,b_z_super
real*8 :: c_x_super,c_y_super,c_z_super
real*8 ::  pi=DACOS(-1.0D0)
real*8 :: displace
real*8,allocatable :: x_frac(:),y_frac(:),z_frac(:)
real*8,allocatable :: x(:),y(:),z(:),at_mass(:),atom_ineq_mass(:)
real*8,allocatable :: x_super(:),y_super(:),z_super(:)
real*8,allocatable :: x_super_frac(:),y_super_frac(:),z_super_frac(:)
real*8 :: r(6,3),r_frac(6,3),angle(3)
! real*8 :: allocatable :: x_bottom(:),y_bottom(:),z_bottom(:)
integer :: n_bottom
character*200 :: line1,line2,line3,line,string, at_type(100)=" ", project
character*10, allocatable :: at_symbol(:),at_ineq_symbol(:),at_ineq_type(:)
character*100  :: char1, char2, char3 
real*8, allocatable ::  A(:,:), Ainv(:,:)
real*8, allocatable :: work(:)
integer, allocatable :: ipiv(:)
integer :: n, info  
character*10 :: style 
logical :: symbol
! define the functios
integer :: get_z
real*8 :: atomic_mass

 

open(unit=1,file="CONTCAR")


!!!! to determine whether we have the line showing atom symbols

do i=1,7
  read(1,*)
end do

read(1,*) line

if(line .eq. "Direct" .or. line .eq. "cart") then
  symbol=.True.
else
  symbol=.False.
end if

rewind(1)

if(symbol) then
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

allocate(n_atom(n_type))
 

rewind(1)  

read(1,"(A)") line1
read(1,*) lattice
read(1,*) a_x,a_y,a_z
read(1,*) b_x,b_y,b_z
read(1,*) c_x,c_y,c_z 
if(symbol) then
  read(1,"(A)") line3
else
  line3=line1
end if
read(1,*) n_atom(:)
backspace(1)
read(1,"(A)") line2
read(1,*) line


total_atom=0
do i=1,n_type
  total_atom=total_atom+n_atom(i)
end do

 a_x=a_x*lattice; a_y=a_y*lattice; a_z=a_z*lattice
 b_x=b_x*lattice; b_y=b_y*lattice; b_z=b_z*lattice
 c_x=c_x*lattice; c_y=c_y*lattice; c_z=c_z*lattice



allocate(x_frac(total_atom),y_frac(total_atom),z_frac(total_atom))
allocate(x(total_atom),y(total_atom),z(total_atom))
allocate(at_symbol(total_atom),at_mass(total_atom))

if(line .eq. "Direct") then
  do i=1,total_atom
    read(1,*) x_frac(i),y_frac(i),z_frac(i)
  end do

  do i=1,total_atom
    x(i)=a_x*x_frac(i) + b_x*y_frac(i) + c_x*z_frac(i)
    y(i)=a_y*x_frac(i) + b_y*y_frac(i) + c_y*z_frac(i)
    z(i)=a_z*x_frac(i) + b_z*y_frac(i) + c_z*z_frac(i)
  end do
else
  do i=1,total_atom
    read(1,*) x(i),y(i),z(i)
  end do
end if

lattice_x=dsqrt(a_x**2+a_y**2+a_z**2)
lattice_y=dsqrt(b_x**2+b_y**2+b_z**2)
lattice_z=dsqrt(c_x**2+c_y**2+c_z**2)

angle(1)=dacos((b_x*c_x+b_y*c_y+b_z*c_z)/(lattice_y*lattice_z))/pi*180.0
angle(2)=dacos((a_x*c_x+a_y*c_y+a_z*c_z)/(lattice_x*lattice_z))/pi*180.0
angle(3)=dacos((a_x*b_x+a_y*b_y+a_z*b_z)/(lattice_x*lattice_y))/pi*180.0


k=0
do i=1,n_type
  do j=1,n_atom(i)
    k=k+1
    write(at_symbol(k),"(A2)") at_type(i)
    at_mass(k)=atomic_mass(get_z(at_symbol(k)))
  end do
end do
 
 
open(unit=3,file="symmetry")
read(3,*)
read(3,"(A)") line
line=trim(line)
line=adjustl(line)
i=index(line, "(")    !!! index: determine the location of a specific character
j=index(line, ")")
read(line(i+1:j-1),*)  n_spacegroup

!!! !!! initial value 
if(n_spacegroup== 62) then
  n_crystal= 78              
else if(n_spacegroup== 143) then
  n_crystal= 174  
else if(n_spacegroup== 156) then
  n_crystal= 190 
else if(n_spacegroup== 164) then
  n_crystal= 200 
else if(n_spacegroup== 187) then
  n_crystal= 225 
else if(n_spacegroup== 194) then
  n_crystal= 232 
else 
  n_crystal= n_spacegroup 
end if
  
do
  read(3,*) char1
  if(char1 == "atom_mapping:") then
    exit
  end if
end do

allocate(index_all(total_atom))
do i=1,total_atom
  read(3,*) char1, index_all(i)
end do

n_atom_ineq=0
do i=1,total_atom
  do j=1,total_atom
    if(index_all(j) .eq. i) then
      n_atom_ineq=n_atom_ineq+1
      exit
    end if
  end do
end do

allocate(index_ineq(n_atom_ineq),at_ineq_symbol(n_atom_ineq))



index_ineq(1)=index_all(1)

do i=2,n_atom_ineq
  do j=2,total_atom
    tmp=0
    do k=1,i-1
      if(index_all(j) .eq. index_ineq(k)) then
        tmp=tmp+1
        exit
      end if
    end do
    if(tmp .eq. 0) then
      index_ineq(i)=index_all(j)
      exit
    end if
  end do
end do

allocate(n_atom_ineq_type(n_type))

n_atom_ineq_type=0
k=0
do i=1,n_type
  do j=1,n_atom(i)
    k=k+1
    do m=1,n_atom_ineq
      if(index_ineq(m)==k) then
        n_atom_ineq_type(i)=n_atom_ineq_type(i)+1
        exit
      end if
    end do
  end do
end do 




open(unit=2,file="atom_nonequivalent")

write(2,"(I0,A,I0)") n_spacegroup,"  ",n_crystal
write(2,"(I0)") n_atom_ineq

k=0
do i=1,n_type
  do j=1,n_atom_ineq_type(i)
    k=k+1
    write(2,"(A,I0,A,I0,A,F8.3)") trim(at_symbol(index_ineq(k))),j,"  ",index_ineq(k),"  ",at_mass(index_ineq(k))
  end do
end do



end program phonon_parameter


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
     &       32.066,  35.453,  39.948,  39.0983,  40.078,            &
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












