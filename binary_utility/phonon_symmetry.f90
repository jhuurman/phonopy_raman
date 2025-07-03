program symmetry
implicit none

integer :: i,j,k
integer :: n_type,total_atom,n_mode
real*8 :: a_x,a_y,a_z
real*8 :: b_x,b_y,b_z
real*8 :: c_x,c_y,c_z
integer,allocatable :: n_atom(:)
real*8,allocatable :: x_frac(:),y_frac(:),z_frac(:)
real*8,allocatable :: x(:),y(:),z(:),x_tmp(:),y_tmp(:),z_tmp(:)
real*8,allocatable :: dx(:),dy(:),dz(:)
real*8,allocatable :: frequency(:)
real*8 :: l_cylinder,l_cone,factor
character*100 :: line,line1,line2,line3,at_type(100)=" "
character*30 :: char_1,char_2,char_3
character*30 :: filename
character*30, allocatable :: representation(:)
logical :: symbol



open(unit=1,file="CONTCAR")


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


total_atom=0
do i=1,n_type
  total_atom=total_atom+n_atom(i)
end do

! allocate(x_frac(total_atom),y_frac(total_atom),z_frac(total_atom))
! allocate(x(total_atom),y(total_atom),z(total_atom))
! allocate(x_tmp(total_atom),y_tmp(total_atom),z_tmp(total_atom))
! allocate(dx(total_atom),dy(total_atom),dz(total_atom))
! 
! do i=1,total_atom
!   read(1,*) x_frac(i),y_frac(i),z_frac(i)
! end do
! 
! do i=1,total_atom
!   x(i)=a_x*x_frac(i) + b_x*y_frac(i) + c_x*z_frac(i)
!   y(i)=a_y*x_frac(i) + b_y*y_frac(i) + c_y*z_frac(i)
!   z(i)=a_z*x_frac(i) + b_z*y_frac(i) + c_z*z_frac(i)
! end do



n_mode=total_atom*3


allocate(representation(n_mode))
allocate(frequency(n_mode))

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

open(unit=2,file="all_mode.txt")
read(2,"(A)") line1
do i=1,n_mode
  read(2,*) j, frequency(i)
end do 
close(2)

open(unit=2,file="all_mode.txt") 
write(2,"(A)") line1
do i=1,n_mode
  write(2,"(I0,A,F6.2,A,A)") i,"   ",frequency(i),"      ",representation(i)
end do 
close(unit=2)

 

end program symmetry