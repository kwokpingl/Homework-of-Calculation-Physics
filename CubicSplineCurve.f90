program main
    implicit none
    integer i,j,k
    real*8 Y(0:15) ,X(0:15) ,m(0:15) ,M1(0:15) ,a(14) ,b(14) ,h(0:14) ,a1(14) ,b1(14)
    real*8 xi ,result
    write(*,*) "Cubic Spline Curve :"
    write(*,*) "Please input the x :"
    read *,xi
    !Read data knew
    open(55,file='dotKnew.txt')
    do i=0,15
            read(55,*) X(i) ,Y(i)
    enddo
    close(55)
    !Calculate the h(i)
    do i=0,14
        h(i) = X(i+1)-X(i)
    enddo
    !Calculate the a(i) & b(i)
    do k=1,14       ! San Wan Ju Equation
        a(k) = h(k-1)/(h(k-1)+h(k))
        b(k) = 3*( (1-a(k))*(Y(k)-Y(k-1))/h(k-1) + a(k)*(Y(k+1)-Y(k))/h(k) )
    enddo
    do k=1,14       ! San Zhuan Jiao Equation
        a1(k) = h(k)/(h(k-1)+h(k))
        b1(k) = 6/(h(k)+h(k-1))*( (Y(k+1)-Y(k))/h(k) + (Y(k)-Y(k-1))/h(k-1) )
    enddo
    !Build the matrix A ,for AM=b
    
end program main