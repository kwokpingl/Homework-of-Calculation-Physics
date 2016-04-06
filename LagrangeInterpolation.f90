program main
        implicit none
        integer i,j
        real*8 l(0:15) ,Y(0:15) ,X(0:15)
        real*8 xi ,result
        write(*,*) "Lagrange Interpolation :"
        write(*,*) "Please input the x :"
        !Read data knew
        open(55,file='dotKnew.txt')
        do i=0,15
                read(55,*) X(i) ,Y(i)
        enddo
        close(55)
        
        open(1,file='Lagrangeout.txt')
        !Calculate the l(x)
        xi=-5.0
        do while(xi <=5.0d0)
                do i=0,15
                        l(i)=1
                        do j=0,15
                                if (j /= i) then
                                        l(i)=l(i) * (xi-X(j)) / (X(i)-X(j))
                                endif
                        enddo
                enddo
                !Calculate the L(x)
                result = dot_product(Y(0:15),l(0:15))
                write(1,*) xi ,result
                xi=xi+0.0001
        enddo
        close(1)
end program main