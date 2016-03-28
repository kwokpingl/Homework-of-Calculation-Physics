real*8 function G(x)	!Define function x=G(x)
implicit none
real*8 x
G=(x**3)/3
!G=(3*x)**(1/3)
return
end function

program main

real*8 Xstar,Xk,Size,nowSize,count		!Xk is the root in loop K，Xstar is the root in loop K+1
real*8,external::G						!Size is supposed Precision，nowSize is now precision
print *,"Jacobi Iteration"
print *,"Please input the origin value for loop："

read *,Xk
nowSize=999		!nowSize is unknow in the first time，initialize a value which is bigger than supposed
Size=0.001		!Set the supposed precision
count=0			!set a loop-time counter

do while( nowSize > Size )	!Keep loop until get the root in precision

count=count+1
Xstar=G(Xk)
nowSize=abs( Xk - Xstar )
Xk=Xstar		!One loop  finish,update the root

print *,"The root is： ",Xk,"Loop Times: ",count

enddo
end program
