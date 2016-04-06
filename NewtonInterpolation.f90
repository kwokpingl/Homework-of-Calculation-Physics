program main
    implicit none
    integer i,j,k
    real*8 a(0:15) ,Y(0:15) ,X(0:15)
    real*8 xi ,result ,temp
    real*8,external::fun
    write(*,*) "Newton Interpolation :"
    write(*,*) "Reading the x and y from file..."
    !Read data knew
    open(55,file='dotKnew.txt')
    do i=0,15
        read(55,*) X(i) ,Y(i)
    enddo
    close(55)
    !Create a new file to write data
    open(1,file='Newtonout.txt')
    xi=-5.0d0
    write(*,*) "Calculating the N(x)..."
    do while(xi<=5.0d0)
        a(:)=0
        result=0
        !Calculate the a(i)
        do i=0,15
            do j=0,i
                temp=Y(j)
                do k=0,i
                    if (k/=j) then
                        temp = temp/( X(j)-X(k) )
                    endif
                enddo
                a(i)=a(i)+temp
            enddo
        enddo
        !Calculate the N(x)
        do i=0,15
            temp=a(i)
            do j=0,i-1
                temp = temp*( xi-X(j) )
            enddo
            result = result +temp
        enddo
        write(1,*) xi,result
        xi=xi+0.01
    enddo
    close(1)
    write(*,*) "Succed .The data is written in 'Newtonout.txt' ."
end program main