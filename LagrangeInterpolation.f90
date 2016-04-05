program main
        implicit none
        integer i,j
        real*8 l(15) ,Y(15) ,X(15)
        real*8 xi ,result
        !Read data knew
        open(55,file='dotKnew.txt')
        do i=1,15
                read(55,*) X(i) ,Y(i)
        enddo
        close(55)
        !Calculate the L(x)
        write(*,*) "Lagrange Interpolation :"
        write(*,*) "Please input the x :"
        read *,xi
        do i=1,15
                l(i)=1
                do j=1,15
                        if (j /= i) then
                                l(i)=l(i) * (xi-X(j)) / (X(i)-X(j))
                        endif
                enddo
        enddo
        result = dot_product(Y(1:15),l(1:15))
        write(*,*) "n = ", 15 ,"L(x) = ",result
end program main