program phonon_initialization
implicit none

!!!!!!!!!!!!!!!!!!!!!!!!
!!! compilation command
! ifort discar.f90 lib/lapack.a lib/blas.a -o discar
!!!!!!!!!!!!!!!!!!!!!
integer :: i,j,k
integer :: n_trans_a,n_trans_b,n_trans_c
integer :: n_type,total_atom,total_atom_super,n_atom_ineq,n_spacegroup,n_crystal
integer,allocatable :: n_atom(:),n_atom_super(:),index_ineq(:)
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
real*8,allocatable :: x(:),y(:),z(:),atom_mass(:)
real*8,allocatable :: x_super(:),y_super(:),z_super(:)
real*8,allocatable :: x_super_frac(:),y_super_frac(:),z_super_frac(:)
real*8 :: r(6,3),r_frac(6,3),angle(3)
! real*8 :: allocatable :: x_bottom(:),y_bottom(:),z_bottom(:)
integer :: n_bottom
character*200 :: line1,line2,line3,line,string, at_type(100)=" ", project
character*10, allocatable :: at_symbol(:),at_ineq_symbol(:)
real*8, allocatable ::  A(:,:), Ainv(:,:)
real*8, allocatable :: work(:)
integer, allocatable :: ipiv(:)
integer :: n, info  
character*10 :: style 
logical :: symbol
 
external DGETRF
external DGETRI

write(*,*) "Please input the project name:"
read(*,*) project


open(unit=1,file="BPOSCAR")
open(unit=2,file="atom_nonequivalent")

read(2,*) n_spacegroup,n_crystal
! read(2,*) angle(:)

read(2,*) n_atom_ineq

allocate(index_ineq(n_atom_ineq),at_ineq_symbol(n_atom_ineq),atom_mass(n_atom_ineq))
do i=1,n_atom_ineq
  read(2,*) at_ineq_symbol(i),index_ineq(i),atom_mass(i)
end do




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
allocate(at_symbol(total_atom))

if(line .eq. "Direct") then
  do i=1,total_atom
    read(1,*) x_frac(i),y_frac(i),z_frac(i)
!     read(1,*) x_frac(i),y_frac(i),z_frac(i), at_symbol(i)
  end do

  do i=1,total_atom
    x(i)=a_x*x_frac(i) + b_x*y_frac(i) + c_x*z_frac(i)
    y(i)=a_y*x_frac(i) + b_y*y_frac(i) + c_y*z_frac(i)
    z(i)=a_z*x_frac(i) + b_z*y_frac(i) + c_z*z_frac(i)
  end do
else
  do i=1,total_atom
    read(1,*) x(i),y(i),z(i)
!      read(1,*) x(i),y(i),z(i), at_symbol(i)
  end do
end if

lattice_x=dsqrt(a_x**2+a_y**2+a_z**2)
lattice_y=dsqrt(b_x**2+b_y**2+b_z**2)
lattice_z=dsqrt(c_x**2+c_y**2+c_z**2)

angle(1)=dacos((b_x*c_x+b_y*c_y+b_z*c_z)/(lattice_y*lattice_z))/pi*180.0
angle(2)=dacos((a_x*c_x+a_y*c_y+a_z*c_z)/(lattice_x*lattice_z))/pi*180.0
angle(3)=dacos((a_x*b_x+a_y*b_y+a_z*b_z)/(lattice_x*lattice_y))/pi*180.0



open(unit=7,file=trim(project)//".d00")
open(unit=8,file=trim(project)//".d01")

write(7,"(2A)") trim(adjustl(project)),"                             "
write(7,"(I6,A,I6,A)") n_spacegroup," 2 ",n_crystal," 1 Number of space group, ab initio/modeling"
write(7,"(A)")  " Unit Cell  " 
write(7,"(F18.12,A)") lattice_x, " Lattice constant A"
write(7,"(F18.12,A)") lattice_y, " Lattice constant B"
write(7,"(F18.12,A)") lattice_z, " Lattice constant C"
write(7,"(F18.12,A)") angle(1), " Lattice angle alpha"
write(7,"(F18.12,A)") angle(2), " Lattice angle beta"
write(7,"(F18.12,A)") angle(3), " Lattice angle gamma"
write(7,"(A)") "      0 Number of non-equiv.elastic"
write(7,"(I7,A)") n_atom_ineq, " Number of non-equiv.displacive"

write(8,"(A)") "  2     Spheres/supercell"
do i=1,n_atom_ineq
  write(8,"(2A)") at_ineq_symbol(i),"    Name of displacive site"
  write(8,"(3F18.12,F11.3,A)") x_frac(index_ineq(i)),y_frac(index_ineq(i)),z_frac(index_ineq(i)), atom_mass(i),"   Positions of displacive site and mass"
  write(8,"(A)")    "   1  1  1  0  0  0    Degree of freedom X Y Z Rx Ry Rz"
end do
do i=1,n_atom_ineq
  write(8,"(2A)") at_ineq_symbol(i),"    0 Number of mom.inert.elements"
end do
do i=1,n_atom_ineq
  write(8,"(A)") "         1.0000        0.0000        0.0000        0.0000        0.0000        0.0000   Elements. of mom.inert"
end do
write(8,"(A)") "        1.0000000    10.0000000   0       1.000000 Diel.con rho neut power"
do i=1,n_atom_ineq
  write(8,"(A,I4,A)") at_ineq_symbol(i), n_atom_ineq, "   Number of eff.charge elements"
end do
do i=1,n_atom_ineq
  write(8,"(A)") "   0.0000000     0.0000000     0.0000000     0.0000000     0.0000000     0.0000000     0.0000000     0.0000000     0.0000000   Elements. of eff.charge"
end do
write(8,"(A)") "      1     Number of shells for displacive p."










end program phonon_initialization














