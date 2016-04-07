program main
    implicit none
    integer i,j,k
    real*8 Y(0:15) ,X(0:15) ,m(0:15) ,M1(0:15) ,h(0:14)
    real*8 a(0:14) ,b(14) ,a1(0:14) ,b1(14)     !Define a(0) to default ,which is of no use in calculation ,while calculating Q(i,j) for QM=b
    real*8 Q(15,0:16) ,Q1(15,0:16)                 !Define Q(0) & Q(16) to default ,in case of over boundary
    real*8 xi ,result ,result1
    write(*,*) "Cubic Spline Curve :"
    write(*,*) "Please input the x :"
    read *,xi
    !Read data knew
    open(55,file='dotKnew.txt')
    do i=0,15
        read(55,*) X(i) ,Y(i)       !Read X,Y
    enddo
    close(55)
    !Calculate the h(i)
    do i=0,14
        h(i) = X(i+1)-X(i)
    enddo
    !Calculate the a(i) & b(i)
    do k=1,14
        ! San Wan Ju Equation
        a(k) = h(k-1)/(h(k-1)+h(k))
        b(k) = 3*( (1-a(k))*(Y(k)-Y(k-1))/h(k-1) + a(k)*(Y(k+1)-Y(k))/h(k) )
        ! San Zhuan Jiao Equation
        a1(k) = h(k)/(h(k-1)+h(k))
        b1(k) = 6/(h(k)+h(k-1))*( (Y(k+1)-Y(k))/h(k) + (Y(k)-Y(k-1))/h(k-1) )
    enddo
    !Build the matrix Q ,for Qm=b & QM=b
    do i=1,14
        j=i-1
        Q(i,j)     = 1-a(i-1)
        Q(i,j+1) = 2
        Q(i,j+2) = a(i)
        Q1(i,j)     = 1-a1(i-1)
        Q1(i,j+1) = 2
        Q1(i,j+2) = a1(i)
    enddo
    !Build new b ,for Qm=b & QM=b
    b(1)  = b(1) - (1-a(1))*m(0)
    b(14)= b(14) - a(14)*m(15)
    b1(1)  = b1(1) - (1-a1(1))*M1(0)
    b1(14)= b1(14) - a1(14)*M1(15)

    !Slove equations by Chasing Elimination
    !Build Triangular matrix Q
    do i=1,13
            do j=i+1,14
                    Q(j,i:i+2) = Q(j,i:i+2) - Q(i,i:i+2)/Q(i,i)*Q(j,i)
                    b(j) = b(j) - b(i)/Q(i,i)*Q(j,i)
                    Q1(j,i:i+2) = Q1(j,i:i+2) - Q1(i,i:i+2)/Q1(i,i)*Q1(j,i)
                    b1(j) = b1(j) - b1(i)/Q1(i,i)*Q1(j,i)
            enddo
    enddo
    !Slove Equations
    do i = 14 ,1 ,-1
            m(i) = (Q(i,10)-dot_product(Q(i, i+1:9),m(i+1:9)))/Q(i,i)
            M1(i) = (Q1(i,10)-dot_product(Q1(i, i+1:9),M1(i+1:9)))/Q1(i,i)
    enddo
end program main