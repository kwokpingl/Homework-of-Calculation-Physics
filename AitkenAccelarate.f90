real*8 function G(x)            !Define function x=G(x)
implicit none
real*8 x
G=(x**3)/3
!g=(3*x)**(1/3)
return
end function

real*8 function root(Xk,Xba,Xflow)
implicit none
real*8 Xk,Xba,Xflow
root = Xflow - (Xflow-Xba)**2 / (Xflow-2*Xba+Xk)
return
end function


program main

real*8 Xstar,Xk,Xba,Xflow,L,Size,nowSize,count
real*8,external::G,root
print *,"Aitken Accelarate"
print *,"Please input the origin value for loop："

read *,Xk
nowSize=999
Size=0.001
count=0

do while( nowSize > Size )

count=count+1
Xba=G(Xk)
Xflow=G(Xba)
Xstar=root(Xk,Xba,Xflow)
nowSize=abs( G(Xstar)-Xstar )
Xk=Xstar

print *,"The Root is： ",Xk,"Loop Times: ",count

enddo
end program
