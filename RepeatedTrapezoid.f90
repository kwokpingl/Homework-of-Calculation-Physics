real*8 function f(x)
    implicit none
    real*8 :: x
    f = sin(x)
    return
end function f

program main
    implicit none
    integer i ,n
    real*8 h ,a ,b ,result ,xi
    real*8,external::f
    h=0.1
    a=1.0
    b=5.0
    n=(b-a)/h
    print *,n
    result =h/2.0*( f(a) + f(b) )
    xi = a
    do i=1,n-1
        xi = a+i*h
        result = result + h/2.0*2.0*f(xi)
    enddo
    write(*,*) "The result is : ",result
end program main