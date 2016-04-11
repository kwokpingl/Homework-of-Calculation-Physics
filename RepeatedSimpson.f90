real*8 function f(x)
    implicit none
    real*8 :: x
    f = sin(x)
    return
end function f

program main
    implicit none
    integer i ,m
    real*8 h ,a ,b ,result ,x1 ,x2
    real*8,external::f
    h=0.1
    a=1
    b=5
    m=(b-a)/(2.0*h)
    result = h/3.0*( f(a) + f(b) + f(a) )
    do i=1,m-1
        x1 = a+(2*i+1)*h
        x2 = a+(2*i)*h
        result = result + h/3.0* ( 2.0 * (2.0* f(x1) + f(x2) ) )
    enddo
    write(*,*) "The result is : ",result
end program main