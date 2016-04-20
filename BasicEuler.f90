real*8 function f(x,y)      ! Define f=yx
    implicit none
    real*8 :: x ,y
    f = -x*x*y*y
    return
end function f

real*8 function Euler(y,f,h)    ! Define Euler Method
    implicit none
    real*8 :: y ,f ,h
    Euler = y + h*f
    return
end function Euler

program main
    implicit none
    integer i ,j ,N
    real*8,external::f ,Euler
    real*8 a ,b ,h ,x ,y
    ! Initialize the index
    a=0d0
    b=1.5d0
    h=0.1d0/8d0
    N=(b-a)/h
    y=3.0d0
    ! Calculate the y
    open(1,file='bEulerout.txt')
    do i=0,N
        x=a+i*h
        y=Euler(y,f(x,y),h)
        write(1,*) x ,y
    enddo
    close(1)
end program main