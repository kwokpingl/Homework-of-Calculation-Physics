real*8 function fun(x)          ! Define f(x)
    implicit none
    real*8 :: x
    fun = 1/(1+x*x)
    return
end function fun

program main
        implicit none
        integer i,j
        real*8 a(15) ,Y(15) ,X(15)
        real*8 xi ,result ,temp
        real*8,external::fun
        !Read data knew
        open(55,file='dotKnew.txt')
        do i=1,15
                read(55,*) X(i) ,Y(i)
        enddo
        close(55)
        !Calculate the a(i)
        do i=1,15
                temp=fun(X(i))
                do j=1,15
                        if (j/=i) then
                                temp = temp/(X(i)-X(j))
                        endif
                enddo
                a(1:i) = a(1:i) + temp
        enddo
        !Calculate the N(x)
        result = 
end program main