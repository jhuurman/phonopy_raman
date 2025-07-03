program broadening
implicit none

integer :: i,j,n,n_point,n_mode,n_samplings,n_insertion,n_normal
integer :: reason
real*8  :: width,HWHM,x_min,x_max,y_max,y0_max,y_total_max,intensity !!! full width at half maximum 
real*8,allocatable :: x0(:),y0(:),x(:),y(:),x_total(:),y_total(:)
real*8,allocatable :: line_width(:)
real*8 :: integral0, integral
character*80 :: filename,filename_broadening,line

open(unit=1,file="broadening_input")
read(1,*) filename
read(1,*) n_mode
read(1,*) HWHM
read(1,*) n_insertion
read(1,*) n_normal
close(1)



open(unit=2,file=trim(filename))
n=0
do
  read(2,*,end=10) line
  if(line .ne. "") then
    n=n+1
  else 
    exit
  end if
end do

10 rewind(2)

n=n+1

n_samplings=n*n_insertion

allocate(x0(n),y0(n),x(n_samplings),y(n_samplings))
allocate(x_total(n+n_samplings),y_total(n+n_samplings))
allocate(line_width(n))

if(n_normal==3) then
  do i=1,n-1
    read(2,*) x0(i),y0(i),line_width(i)
  end do
else
  do i=1,n-1
    read(2,*) x0(i),y0(i)
  end do
end if

x0(n)=x0(n-1)+50.0  !! so we can have a zero baseline 
y0(n)=0.0           !! so we can have a zero baseline 

close(2)

x_min=x0(1);x_max=x0(n)

do i=1,n_samplings
  x(i)=x_min+(x_max-x_min)/(n_samplings-1)*(i-1)
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! x_total is to ensure original data points included
do i=1,n
  x_total(i)=x0(i)
end do

do i=n+1,n+n_samplings
  x_total(i)=x(i-n)
end do

call sort(n+n_samplings,x_total)

!!!!!!!!!!!!!!!!!!!

if(n_mode==1) then   !!! Gussian broadening
  if(n_normal==3) then
    line_width=line_width/sqrt(log(2.0))
    do i=1,n_samplings
      y(i)=0
      do j=1,n
        y(i)=y(i)+y0(j)*exp(-(x(i)-x0(j))**2/line_width(j)**2)
      end do
    end do
  
    do i=1,n+n_samplings
      y_total(i)=0
      do j=1,n
        y_total(i)=y_total(i)+y0(j)*exp(-(x_total(i)-x0(j))**2/line_width(j)**2)
      end do
    end do
    line_width=line_width*sqrt(log(2.0))
  else 
    width=HWHM/sqrt(log(2.0))  !!! for Gussian broadening
    do i=1,n_samplings
      y(i)=0
      do j=1,n
        y(i)=y(i)+y0(j)*exp(-(x(i)-x0(j))**2/width**2)
      end do
    end do
  
    do i=1,n+n_samplings
      y_total(i)=0
      do j=1,n
        y_total(i)=y_total(i)+y0(j)*exp(-(x_total(i)-x0(j))**2/width**2)
      end do
    end do
  end if
else if(n_mode==2) then  !!! Lorentzian broadening
  if(n_normal==3) then
    do i=1,n_samplings
      y(i)=0
      do j=1,n
        y(i)=y(i)+y0(j)*(line_width(j))**2/((line_width(j))**2+(x(i)-x0(j))**2)
      end do
    end do
  
    do i=1,n+n_samplings
      y_total(i)=0
      do j=1,n
        y_total(i)=y_total(i)+y0(j)*(line_width(j))**2/((line_width(j))**2+(x_total(i)-x0(j))**2)
      end do
    end do  
  else
    width=HWHM
    do i=1,n_samplings
      y(i)=0
      do j=1,n
        y(i)=y(i)+y0(j)*(width)**2/((width)**2+(x(i)-x0(j))**2)
      end do
    end do
  
    do i=1,n+n_samplings
      y_total(i)=0
      do j=1,n
        y_total(i)=y_total(i)+y0(j)*(width)**2/((width)**2+(x_total(i)-x0(j))**2)
      end do
    end do
  end if
else
  write(*,*) "Warning, please choose the right mode!!!!"
end if


open(unit=3,file=trim(filename)//"_broadening")


if(n_normal==1) then                           !!!for dos, keep the integral of dos the same as number of electrons
  call integration(x0(1),x0(n),n,y0,integral0)
  call integration(x(1),x(n_samplings),n_samplings,y,integral)
  do i=1,n_samplings
    write(3,*) x(i),y(i)/integral*integral0    
  end do
else if(n_normal==2) then                    !!!for raman, keep the peak intensity the same
  y0_max=maxval(y0)
  y_total_max=maxval(y_total)
  do i=1,n+n_samplings
    write(3,*) x_total(i),y_total(i)
  end do
 else if(n_normal==3) then
  do i=1,n+n_samplings
    write(3,*) x_total(i),y_total(i)
  end do 
  
  open(unit=9,file=trim(filename)//"_intensity")
  do j=1,n
    n_point=0
    intensity=0.0
    do i=1,n+n_samplings
      if(abs(x_total(i)-x0(j))<=2.0*line_width(j)) then
        intensity=intensity+y_total(i)
        n_point=n_point+1
      end if
    end do
    intensity=intensity/dfloat(n_point)*4.0*line_width(j)
    write(9,"(2F16.3)") x0(j),intensity
  end do
  close(9)  
 else 
   write(*,*) "Warning, please choose the right normalization mode!!!!"
end if

end program broadening


subroutine integration(a,b,n,f,y)
integer :: i,n
real*8 :: a,b,h,y
real*8 :: f(n)

h=(b-a)/dfloat(n-1)
y=0.0
do i=2,n-1 
  y=y+f(i)*h
end do
y=y+(f(1)+f(n))*h*0.5

end subroutine integration



subroutine sort(n,f)
implicit none

integer :: i,j,n
integer :: auxiliary
real*8 :: first
real*8 :: f(n)

do i=1,n
  auxiliary=i
  first=f(auxiliary)
  do j=i,n
    if(f(j).lt.first) then
      auxiliary=j
      first=f(j)
    end if
  end do
  if(auxiliary.ne.i) then
    f(auxiliary)=f(i) 
    f(i)=first
  end if
end do


end subroutine sort