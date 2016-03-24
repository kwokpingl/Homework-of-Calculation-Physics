real*8 function G(x)		!Define function x=G(x)
implicit none
real*8 x
G=(x**3)/3
!g=(3*x)**(1/3)
return
end function

real*8 function Gx(x)		!!Define function G'(x)
implicit none
real*8 x
Gx=x**2
!gx=(3*x)**(-2/3)
return
end function

real*8 function root(x,L)
implicit none
real*8 x,L
root=( (x**3)/3 - L*x )/(1-L)
!root=( (3*x)**(1/3) - L*x )/(1-L)
return
end function


program main

real*8 Xstar,Xk,L,Size,nowSize,count
real*8,external::G,Gx,root
print *,"Post Accelarate"
print *,"Please input the origin value for loop："

read *,Xk
nowSize=999		!nowSize is unknow in the first time，initialize a value which is bigger than supposed
Size=0.001		!Set the supposed precision
count=0			!set a loop-time counter

do while( nowSize > Size )		!Keep loop until get the root in precision

count=count+1
L=Gx(Xk)
Xstar=root(Xk,L)
nowSize=abs( G(Xstar)-Xstar )
Xk=Xstar

print *,"The root is： ",Xk,"Loop Times: ",count

enddo
pause
end program
