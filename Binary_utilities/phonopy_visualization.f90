program symmetry
implicit none

integer :: i,j,k,m,n
integer :: n_trans_a,n_trans_b,n_trans_c
integer :: n_type,total_atom,n_mode,total_atom_super
integer :: super(3)
real*8 :: a_x,a_y,a_z
real*8 :: b_x,b_y,b_z
real*8 :: c_x,c_y,c_z
real*8 :: a_x_super,a_y_super,a_z_super
real*8 :: b_x_super,b_y_super,b_z_super
real*8 :: c_x_super,c_y_super,c_z_super
integer,allocatable :: n_atom(:),n_atom_super(:),atomic_number(:)
real*8,allocatable :: x_frac(:),y_frac(:),z_frac(:)
real*8,allocatable :: x_super(:),y_super(:),z_super(:)
real*8,allocatable :: x(:),y(:),z(:),x_tmp(:),y_tmp(:),z_tmp(:)
real*8,allocatable :: dx(:,:),dy(:,:),dz(:,:),displace(:,:),d_max(:)
real*8,allocatable :: frequency(:),species_mass(:)
real*8 :: l_cylinder,l_cone,factor
real*8 :: z_max,z_min,surface_depth
character*100 :: line,line1,line2,line3,at_type(100)=" "
character*30 :: char_1,char_2,char_3
character*30 :: filename,supercell
integer :: option, flag
logical :: symbol 
logical :: existence
! define the functions
integer :: get_z
real*8 :: atomic_mass

open(unit=1,file="CONTCAR")
open(unit=2,file="band.yaml")


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


allocate(n_atom(n_type),atomic_number(n_type))

rewind(1)



read(1,"(A)") line1
read(1,*) 
read(1,*) a_x,a_y,a_z
read(1,*) b_x,b_y,b_z
read(1,*) c_x,c_y,c_z 
if(symbol) then
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

do i=1,total_atom
  x(i)=a_x*x_frac(i) + b_x*y_frac(i) + c_x*z_frac(i)
  y(i)=a_y*x_frac(i) + b_y*y_frac(i) + c_y*z_frac(i)
  z(i)=a_z*x_frac(i) + b_z*y_frac(i) + c_z*z_frac(i)
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

! z_min=minval(z(:))
! z_max=maxval(z(:))

n_mode=total_atom*3

allocate(dx(n_mode,total_atom),dy(n_mode,total_atom),dz(n_mode,total_atom))

allocate(frequency(n_mode))

do    
  read(2,"(A)") line  
  line=trim(adjustl(line))    
  read(line,*,end=101) char_1, char_2
  101 if(char_1 .eq. "-" .and. char_2 .eq. "#") then    
   exit    
  end if    
end do   

backspace(2)

do i=1,n_mode
  read(2,*) 
  read(2,*) char_1,frequency(i)
  read(2,*) 

  do j=1,total_atom
    read(2,*) 
    read(2,*) char_1,char_2,dx(i,j)
    read(2,*) char_1,char_2,dy(i,j)
    read(2,*) char_1,char_2,dz(i,j)
  end do
end do

!!!!! normalization of the phonon eigenvectors (eigenmodes) by the root of the atomic mass to become phonon eigendisplacements (atomic displacements) 
do i=1,n_mode
  do j=1,total_atom
    dx(i,j)=dx(i,j)/sqrt(species_mass(j))
    dy(i,j)=dy(i,j)/sqrt(species_mass(j))
    dz(i,j)=dz(i,j)/sqrt(species_mass(j))
  end do
end do


write(*,*) "Please input the normalized factor to get proper cylinder length: "
read(*,*)  factor

!! after normalized, dx, dy, dz are too small, so they are enlarged by the same factor

factor=factor*maxval(sqrt(species_mass(:)))

! write(*,*) "Please indicate whether you want to draw in white or black: 1 or 2"
! read(*,*) flag

! write(*,*) "Please input the unit you want for the frequency: 1 for cm-1 and 2 for meV:"
! read(*,*) option

l_cylinder=4.0*factor
l_cone=1.5*factor

do i=1,n_mode
  write(filename,"(A,I0)") "mode",i
  open(unit=3,file=filename)
!   if(flag==1) then
!     write(3,*) "draw color white"
!   else
!     write(3,*) "draw color black"
!   end if
  write(3,*) "draw color blue"
  
  do j=1,total_atom
    x_tmp(j)=x(j)+l_cylinder*dx(i,j)
    y_tmp(j)=y(j)+l_cylinder*dy(i,j)
    z_tmp(j)=z(j)+l_cylinder*dz(i,j)
    write(3,*) "draw cylinder {", x(j), y(j), z(j),"} {", x_tmp(j), y_tmp(j), z_tmp(j),"} radius 0.1 resolution 80 filled yes"
    write(3,*) "draw cone  {",x_tmp(j), y_tmp(j), z_tmp(j),"} {", &
    x_tmp(j)+l_cone*dx(i,j),y_tmp(j)+l_cone*dy(i,j),z_tmp(j)+l_cone*dz(i,j),"} &
radius 0.2 resolution 80"
  end do
  close(unit=3)
    
end do

close(unit=1)

option=1

inquire(file="all_mode.txt",exist=existence)
if(.not. existence) then
open(unit=4,file="all_mode.txt") 

if(option==1) then
!   write(4,*) "# unit cm-1"
  do i=1,n_mode
    write(4,"(I3,A,F10.2)") i,"   ",frequency(i)*33.356407889
  end do
else
!   write(4,*) "# unit meV"
  do i=1,n_mode
  write(4,"(I3,A,F10.2)") i,"   ",frequency(i)*4.135668892
  end do 
end if

close(unit=4)
end if

write(*,*) "Do you want to visualize the vibrations in a supercell instead of the default primitive unit cell: type yes or no"
read(*,*) supercell

if(supercell .eq. "yes") then
  write(*,*) "Please input the dimension of supercell:"
  read(*,*) super(:)
  
  a_x_super=a_x*super(1); a_y_super=a_y*super(1); a_z_super=a_z*super(1)
  b_x_super=b_x*super(2); b_y_super=b_y*super(2); b_z_super=b_z*super(2)
  c_x_super=c_x*super(3); c_y_super=c_y*super(3); c_z_super=c_z*super(3)
  
  allocate(n_atom_super(n_type))
  
  n_atom_super(:)=n_atom(:)*super(1)*super(2)*super(3)
  total_atom_super=total_atom*super(1)*super(2)*super(3)
  
  deallocate(x_tmp,y_tmp,z_tmp)
  allocate(x_super(total_atom_super),y_super(total_atom_super),z_super(total_atom_super))
  allocate(x_tmp(total_atom_super),y_tmp(total_atom_super),z_tmp(total_atom_super))

  
  open(unit=1,file="POSCAR_supercell")
  write(1,"(A)") trim(line1)//"  "
  write(1,"(A)") "  1.0000000000000"
  write(1,"(3F15.9)") a_x_super,a_y_super,a_z_super
  write(1,"(3F15.9)") b_x_super,b_y_super,b_z_super
  write(1,"(3F15.9)") c_x_super,c_y_super,c_z_super
  write(1,"(A)") trim(line3)//"  "
  write(1,*) n_atom_super(:)
  write(1,"(A)") "cart"
  
  j=0
  do i=1,total_atom  
    do n_trans_c=0,super(3)-1  
      do n_trans_b=0,super(2)-1  
        do n_trans_a=0,super(1)-1  
          j=j+1  
          x_super(j)=x(i)+a_x*n_trans_a+b_x*n_trans_b+c_x*n_trans_c  
          y_super(j)=y(i)+a_y*n_trans_a+b_y*n_trans_b+c_y*n_trans_c  
          z_super(j)=z(i)+a_z*n_trans_a+b_z*n_trans_b+c_z*n_trans_c  
        end do  
      end do  
    end do  
  end do  
  
  do i=1,total_atom_super
    write(1,"(3F15.9)") x_super(i),y_super(i),z_super(i)
  end do
  
  close(1)
  
!   rewind(2)
!   do    
!     read(2,"(A)") line  
!     line=trim(adjustl(line))    
!     read(line,*,end=102) char_1, char_2
!     102 if(char_1 .eq. "-" .and. char_2 .eq. "#") then    
!      exit    
!     end if    
!   end do   
! 
!   backspace(2)
  
  do i=1,n_mode
!     read(2,*) 
!     read(2,*) char_1,frequency(i)
!     read(2,*) 
  
    write(filename,"(A,I0)") "mode_super",i
    open(unit=3,file=filename)
    write(3,*) "draw color blue"
    
    k=0
    do j=1,total_atom
!       read(2,*) 
!       read(2,*) char_1,char_2,dx(i,j)
!       read(2,*) char_1,char_2,dy(i,j)
!       read(2,*) char_1,char_2,dz(i,j)
      
      do n_trans_c=0,super(3)-1      
        do n_trans_b=0,super(2)-1      
          do n_trans_a=0,super(1)-1      
            k=k+1      
            x_tmp(k)=x_super(k)+l_cylinder*dx(i,j)
            y_tmp(k)=y_super(k)+l_cylinder*dy(i,j)
            z_tmp(k)=z_super(k)+l_cylinder*dz(i,j)
            write(3,*) "draw cylinder {", x_super(k), y_super(k), z_super(k),"} {", &
            x_tmp(k), y_tmp(k), z_tmp(k),"} radius 0.1 resolution 80 filled yes"
            write(3,*) "draw cone  {",x_tmp(k), y_tmp(k), z_tmp(k),"} {", &
            x_tmp(k)+l_cone*dx(i,j),y_tmp(k)+l_cone*dy(i,j),z_tmp(k)+l_cone*dz(i,j),"} &
            radius 0.2 resolution 80"
          end do      
        end do      
      end do      
    end do
    close(unit=3)
      
  end do
  
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
