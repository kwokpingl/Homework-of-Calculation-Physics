real*8 function f(x)        ! Define the origin function f(x)
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
        write(*,*) "Repeated Simpson integration :"
        ! Initialize the index
        h=0.1d0
        a=1.0d0
        b=5.0d0
        m=(b-a)/(2.0d0*h)
        result = h/3.0d0*( f(a) - f(b) )
        ! Calculate the integrate
        do i=1,m
                x1 = a+(2*i-1)*h
                x2 = a+(2*i)*h
                result = result + h/3.0d0* ( 2.0d0 * (2.0d0* f(x1) + f(x2) ) )
        enddo
        write(*,*) "The result is : ",result
end program main