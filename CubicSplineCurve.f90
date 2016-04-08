program main
    implicit none
    integer i,j,k
    real*8 Y(0:15) ,X(0:15) ,m(0:15) ,M1(0:15) ,h(0:14)
    real*8 a(0:14) ,b(14) ,a1(0:14) ,b1(14)     !Define a(0) to default ,which is of no use in calculation ,while calculating Q(i,j) for QM=b
    real*8 Q(14,0:15) ,Q1(14,0:15)                 !Define Q(i,0) & Q(i,15) to default ,in case of over boundary
    real*8 xi ,result ,result1 ,alpha(0:15) ,beta(0:15)
    m(0)    = -(1.0+25.0)**(-2)*2.0*(-5.0)
    m(15)  = -(1.0+25.0)**(-2)*2.0*5.0
    M1(0)  = 2.0*(1.0+25.0)**(-3)*4.0*25.0 - (1.0+25.0)**(-2) *8.0*(-5.0)
    M1(15)= 2.0*(1.0+25.0)**(-3)*4.0*25.0 - (1.0+25.0)**(-2) *8.0*(5.0)
    write(*,*) "Cubic Spline Curve (Both San-Wan-Ju & San-Zhuan-Jiao) :"
    write(*,*) "Reading the x and y from file..."
    !Read data knew
    open(55,file='dotKnew.txt')
    do i=0,15
        read(55,*) X(i) ,Y(i)       !Read X,Y
    enddo
    close(55)
    !Calculate the h(i)
    write(*,*) "Calculating the h(i)..."
    do i=0,14
        h(i) = X(i+1)-X(i)
    enddo
    !Calculate the a(i) & b(i)
    write(*,*) "Calculating the a(i) & b(i)..."
    do k=1,14
        ! San Zhuan Jiao Equation
        a(k) = h(k-1)/(h(k-1)+h(k))
        b(k) = 3.0*( (1.0-a(k))*(Y(k)-Y(k-1))/h(k-1) + a(k)*(Y(k+1)-Y(k))/h(k) )
        ! San Wan Ju Equation
        a1(k) = h(k)/(h(k-1)+h(k))
        b1(k) = 6.0/(h(k)+h(k-1))*( (Y(k+1)-Y(k))/h(k) - (Y(k)-Y(k-1))/h(k-1) )
    enddo
    !Build new b ,for Qm=b & QM=b
    b(1)  = b(1) - (1-a(1))*m(0)
    b(14)= b(14) - a(14)*m(15)
    b1(1)  = b1(1) - (1-a1(1))*M1(0)
    b1(14)= b1(14) - a1(14)*M1(15)
    !Build the matrix Q ,for Qm=b & QM=b
    write(*,*) "Building the Q(i,i)..."
    Q=0
    Q1=0
    do i=1,14
        Q(i,i-1) = 1-a(i-1)         !Q(1,0) is already defined to default in case of over boundary
        Q(i,i) = 2
        Q(i,i+1) = a(i)
        Q1(i,i-1) = 1-a1(i-1)
        Q1(i,i) = 2
        Q1(i,i+1) = a1(i)           !Q(14,15) is already defined to default
    enddo

    !Slove equations by Chasing Elimination
    !Build Triangular matrix Q
    write(*,*) "Building the upper Triangular matrix..."
    do i=1,13
            j=i+1
            Q(j,i:i+1) = Q(j,i:i+1) - Q(i,i:i+1)*Q(j,i)/Q(i,i)
            b(j) = b(j) - b(i)*Q(j,i)/Q(i,i)
            Q1(j,i:i+1) = Q1(j,i:i+1) - Q1(i,i:i+1)*Q1(j,i)/Q1(i,i)
            b1(j) = b1(j) - b1(i)*Q1(j,i)/Q1(i,i)
    enddo
    !Slove Equations QM=b
    Q(14,15)=0
    Q1(14,15)=0
    write(*,*) "Calculating the m(i) & M(i)..."
    do i = 14 ,1 ,-1
            m(i)   = (b(i)  -(Q(i, i+1)   *m(i+1)))  /Q(i,i)
            M1(i) = (b1(i)-(Q1(i, i+1)*M1(i+1)))/Q1(i,i) 
    enddo

    write(*,*) "Calculating the S(x)..."
    open(1,file='San-Zhuan-Jiaoout.txt')
    open(2,file='San-Wan-Juout.txt')
    xi=-5.0
    do while (xi<=5.0d0)
        result=0     ! Initialize result
        result1=0   ! Initialize result1
        !Calculate alpha & beta ,in addition to
        !S(x) = Y(i)*alpha(i) +m(i)*beta(i)
        !WRONG PRINT in TEXTBOOK RETORT :
        !In textbook ,alpha(k)=( (xi-X(k-1))/(X(k)-X(k-1)) )**2 * (1+2*(xi-X(k))/(X(k-1)-X(k)))**2
        do k=0,15
            if (k>0 .and. X(k-1)<xi .and. xi <=X(k)) then
                alpha(k)=( (xi-X(k-1))/(X(k)-X(k-1)) )**2 * (1+2*(xi-X(k))/(X(k-1)-X(k)))
                beta(k)  =( (xi-X(k-1))/(X(k)-X(k-1)) )**2 * (xi-X(k))
            else if (k<15 .and. X(k)<xi .and. xi <=X(k+1)) then
                alpha(k)=( (xi-X(k+1))/(X(k)-X(k+1)) )**2 * (1+2*(xi-X(k))/(X(k+1)-X(k)))
                beta(k)  =( (xi-X(k+1))/(X(k)-X(k+1)) )**2 * (xi-X(k))
            else
                alpha(k)=0
                beta(k)=0
            endif
        enddo
        result=dot_product(Y(0:15),alpha(0:15)) + dot_product(m(0:15),beta(0:15))
        write(1,*) xi ,result ,result-1/(1+xi*xi)

        !Calculate the S(i)
        !S(x)=M(k)*(X(k+1)-xi)**3/(6*h(k)) + M(k+1)*(xi-X(k))**3/h(k) + (Y(k)-M(k)*h(k)**2/6)*(X(k+1)-xi)/h(k) + (Y(k+1)-M(k+1)*h(k)**2/6)*(xi-X(k))/h(k)
        do k=0,14
            if (X(k)<=xi .and. xi <X(k+1)) then
                !The equation is too big to complie ,so express it step by step
                result1=M(k) * ((X(k+1)-xi)**3.0) / (6.0*h(k))
                !WRONG PRINT in TEXTBOOK RETORT :
                !result1=result1 + M(k+1) * ((xi-X(k))**3.0) / (h(k)) in textbook!
                result1=result1 + M(k+1) * ((xi-X(k))**3.0) / (6.0*h(k))
                result1=result1 + (Y(k)-M(k) * (h(k)**2.0) / 6.0) * (X(k+1)-xi)/h(k)
                result1=result1 + (Y(k+1)-M(k+1) * (h(k)**2.0) / 6.0) * (xi-X(k))/h(k)
                exit
            endif
        enddo
        write(2,*) xi ,result1 ,result1-1/(1+xi*xi)
        xi=xi+0.01
    enddo
    close(1)
    close(2)
    write(*,*) "Succeed .The data is written in 'San-Zhuan-Jiaoout.txt' and 'San-Wan-Juout.txt' "
end program main