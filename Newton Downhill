real*8 function Fun(x)		!Define the original function
implicit none
real*8 x
Fun=x**3/3-x
return
end

real*8 function root(x,A)	!Define the iteration function，set two index :Xk,λ
implicit none
real*8 x,A
root=x-A*(x**3/3-x)/(x*x-1.0)
return
end

program main
!Xk is the root in loop-K
real*8 Xk,A,delta,xSize,ySize,Asize,usefulA,count
real*8 xTemp,fTemp,f1Temp
real*8,external::Fun,root

print *,"Newton Downhill"
print *,"Please input the origin value for iteration : "

A=1			    	!A is the Downhill-Weight λ,initializ as 1
Asize=0.01		!Asize is the Downhill-Weight lower bound，-
		        	!which means when λ is lower ，you can judge that  it's the root which lead the function to diverge，-
		        	!and you should correct the root value x=x+delta
delta=0.1		  !delta is the small-value for correcting
xSize=0.0001	!Set the root precision
ySize=0.0001	!Set the function value precision for check
count=0		!Set a loop-time counter

read *,Xk

do
	count=count+1
	xTemp=root(Xk,A)
	fTemp=abs(Fun(Xk))
	f1Temp=abs(Fun(xTemp))

	if( f1Temp < fTemp ) then
	
		Xk=xTemp
		print *,"The root is： ",Xk,"Loop Times:",count
	!	print *,"Downhill-weight A = ",A
		if(fTemp<ySize .and. abs(xTemp-Xk)<xSize ) then
			stop
		endif
	
	else
	
		if( A > Asize ) then
			A=A/2
		else
			Xk=Xk+delta
			A=1
		endif
	endif

enddo

endprogram
