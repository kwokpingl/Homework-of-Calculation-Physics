real*8 function f(x,y)      ! Define f=yx
    implicit none
    real*8 :: x ,y
    f = -x*x*y*y
    return
end function f

real*8 function RungeKutta(y,K1,K2,K3,K4,h)    ! Define RungeKutta Method
    implicit none
    real*8 :: y ,K1,K2,K3,K4 ,h
    RungeKutta = y + h/6.0d0*(K1+2.0d0*K2+2.0d0*K3+K4)
    return
end function RungeKutta

program main
    implicit none
    integer i ,j ,N
    real*8,external::f ,RungeKutta
    real*8 a ,b ,h ,x ,y ,K1 ,K2 ,K3 ,K4
    ! Initialize the index
    a=0d0
    b=1.5d0
    h=0.1d0/8d0
    N=(b-a)/h
    y=3.0d0
    ! Calculate the y
    open(1,file='RungeKuttaout.txt')
    do i=0,N
        x=a+i*h
        K1=f(x,y)
        K2=f(x+h/2.0d0,y+h/2.0d0*K1)
        K3=f(x+h/2.0d0,y+h/2.0d0*K2)
        K4=f(x+h,y+h*K3)
        y=RungeKutta(y,K1,K2,K3,K4,h)
        write(1,*) x ,y ,y-3.0d0/(1.0d0+x*x*x)
    enddo
    close(1)
end program main