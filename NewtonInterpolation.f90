real*8 function fun(x)          ! Define f(x)
    implicit none
    real*8 :: x
    fun = 1/(1+x*x)
    return
end function fun

program main
        implicit none
        integer i,j
        real*8 a(0:15) ,Y(0:15) ,X(0:15)
        real*8 xi ,result ,temp
        real*8,external::fun
        write(*,*) "Newton Interpolation :"
        write(*,*) "Please input the x :"
        !Read data knew
        open(55,file='dotKnew.txt')
        do i=0,15
                read(55,*) X(i) ,Y(i)
        enddo
        close(55)

        open(1,file='Newtonout.txt')
        !Calculate the a(i)
        xi=-5.0d0
        do while(xi<=5.0d0)
            do i=0,15
                    temp=fun(X(i))
                    do j=0,15
                            if (j/=i) then
                                    temp = temp/(X(i)-X(j))
                            endif
                    enddo
                    a(0:i) = a(0:i) + temp
            enddo
            !Calculate the N(x)
            do i=0,15
                    temp=a(i)
                    do j=0,15
                            temp = temp*(xi-X(j))
                    enddo
                    result = result +temp
            enddo
            write(1,*) xi,result
            xi=xi+0.0001
        enddo
        close(1)
end program main