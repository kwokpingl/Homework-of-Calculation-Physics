real*8 function f(x)            ! Define the origin function f(x)
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
        write(*,*) "Repeated Trapezoid integration :"
        ! Initialize the index
        h=0.1d0
        a=1.0d0
        b=5.0d0
        n=(b-a)/h
        result =h/2.0*( f(a) + f(b) )
        ! Calculate the integrate
        do i=1,n-1
        xi = a+i*h
        result = result + h*f(xi)
        enddo
        write(*,*) "The result is : ",result
end program main