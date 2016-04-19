program main
    implicit none
    integer x ,y ,i ,j
    real*8 Map(0:1800/20 , 0:1500/20)
    real*8,external::Cubic
    Map=0
    !Calculate the Map
    do x=0,1800/20
        ! Get the X-axis value 
        write(*,*) "Calculating the X-axis value",x
        ! Get f(x,0)
        Map(x,0)=Cubic( Map(0,0),Map(600/20,0),Map(1200/20,0),Map(1800/20,0))
        ! Get f(x,500)
        Map(x,500/20)=Cubic( Map(0,500/20),Map(600/20,500/20),Map(1200/20,500/20),Map(1800/20,500/20))
        ! Get f(x,1000)
        Map(x,1000/20)=Cubic( Map(0,1000/20),Map(600/20,1000/20),Map(1200/20,1000/20),Map(1800/20,1000/20))
        ! Get f(x,1500)
        Map(x,1500/20)=Cubic( Map(0,1500/20),Map(600/20,1500/20),Map(1200/20,1500/20),Map(1800/20,1500/20))
        
        ! Get the Y-axis value
        write(*,*) "Calculating the Y-axis value ,y="
        do y=0,1500/20
            ! Get f(x,y)
            Map(x,y)=Cubic( Map(x,0),Map(x,500/20),Map(x,1000/20),Map(x,1500/20))
            write(*,*) "x = " ,x*20 ,"y = " ,y*20
        enddo
    enddo

    ! Write the result
    write(*,*) "write the result..."
    open(1,file='Map.txt')
    do i=0,1800/20
        do j=0,1500/20
            write(1,*) i*20,Map(i,j*20)
        enddo
    enddo
    close(1)
end program main



real*8 function Cubic(num0,num1,num2,num3)
    implicit none
    integer i,j,k
    real*8 :: num0 ,num1 ,num2 ,num3
    real*8 Y(0:3) ,X(0:3) ,m(0:3) ,M1(0:3) ,h(0:2)
    real*8 a(0:2) ,b(2) ,a1(0:2) ,b1(2)     !Define a(0) to default ,which is of no use in calculation ,while calculating Q(i,j) for QM=b
    real*8 Q(2,0:3) ,Q1(2,0:3)                 !Define Q(i,0) & Q(i,3) to default ,in case of over boundary
    real*8 xi ,result ,result1 ,alpha(0:3) ,beta(0:3)
    m(0)    = -(1.0+25.0)**(-2)*2.0*(-5.0)
    m(3)    = -(1.0+25.0)**(-2)*2.0*5.0
    M1(0)  = 2.0*(1.0+25.0)**(-3)*4.0*25.0 - (1.0+25.0)**(-2) *8.0*(-5.0)
    M1(3)  = 2.0*(1.0+25.0)**(-3)*4.0*25.0 - (1.0+25.0)**(-2) *8.0*(5.0)
    write(*,*) "Cubic Spline Curve (Both San-Wan-Ju & San-Zhuan-Jiao) :"
    write(*,*) "Reading the x and y..."
    !Read data knew
    do i=0,3
        X=(/0,1,2,3/)       !Read X,Y
        Y=(/num0,num1,num2,num3/)
    enddo
    !Calculate the h(i)
    write(*,*) "Calculating the h(i)..."
    do i=0,2
        h(i) = X(i+1)-X(i)
    enddo
    !Calculate the a(i) & b(i)
    write(*,*) "Calculating the a(i) & b(i)..."
    do k=1,2
        ! San Zhuan Jiao Equation
        a(k) = h(k-1)/(h(k-1)+h(k))
        b(k) = 3.0d0*( (1.0d0-a(k))*(Y(k)-Y(k-1))/h(k-1) + a(k)*(Y(k+1)-Y(k))/h(k) )
        ! San Wan Ju Equation
        a1(k) = h(k)/(h(k-1)+h(k))
        b1(k) = 6.0d0/(h(k)+h(k-1))*( (Y(k+1)-Y(k))/h(k) - (Y(k)-Y(k-1))/h(k-1) )
    enddo
    !Build new b ,for Qm=b & QM=b
    b(1)  = b(1) - (1-a(1))*m(0)
    b(2)= b(2) - a(2)*m(3)
    b1(1)  = b1(1) - (1-a1(1))*M1(0)
    b1(2)= b1(2) - a1(2)*M1(3)
    !Build the matrix Q ,for Qm=b & QM=b
    write(*,*) "Building the Q(i,i)..."
    Q=0
    Q1=0
    do i=1,2
        Q(i,i-1) = 1-a(i-1)         !Q(1,0) is already defined to default in case of over boundary
        Q(i,i) = 2.0d0
        Q(i,i+1) = a(i)
        Q1(i,i-1) = 1-a1(i-1)
        Q1(i,i) = 2.0d0
        Q1(i,i+1) = a1(i)           !Q(2,3) is already defined to default
    enddo

    !Slove equations by Chasing Elimination
    !Build Triangular matrix Q
    write(*,*) "Building the upper Triangular matrix..."
    do i=1,1
            j=i+1
            Q(j,i:i+1) = Q(j,i:i+1) - Q(i,i:i+1)*Q(j,i)/Q(i,i)
            b(j) = b(j) - b(i)*Q(j,i)/Q(i,i)
            Q1(j,i:i+1) = Q1(j,i:i+1) - Q1(i,i:i+1)*Q1(j,i)/Q1(i,i)
            b1(j) = b1(j) - b1(i)*Q1(j,i)/Q1(i,i)
    enddo
    !Slove Equations QM=b
    Q(2,3)=0
    Q1(2,3)=0
    write(*,*) "Calculating the m(i) & M(i)..."
    do i = 2 ,1 ,-1
            m(i)   = (b(i)  -(Q(i, i+1)   *m(i+1)))  /Q(i,i)
            M1(i) = (b1(i)-(Q1(i, i+1)*M1(i+1)))/Q1(i,i)
    enddo

    write(*,*) "Calculating the S(x)..."
    result=0     ! Initialize result
    result1=0   ! Initialize result1
    !Calculate alpha & beta ,in addition to
    !S(x) = Y(i)*alpha(i) +m(i)*beta(i)
    do k=0,3
        if (k>0 .and. X(k-1)<xi .and. xi <=X(k)) then
            alpha(k)=( (xi-X(k-1))/(X(k)-X(k-1)) )**2 * (1+2*(xi-X(k))/(X(k-1)-X(k)))
            beta(k)  =( (xi-X(k-1))/(X(k)-X(k-1)) )**2 * (xi-X(k))
        else if (k<3 .and. X(k)<xi .and. xi <=X(k+1)) then
            alpha(k)=( (xi-X(k+1))/(X(k)-X(k+1)) )**2 * (1+2*(xi-X(k))/(X(k+1)-X(k)))
            beta(k)  =( (xi-X(k+1))/(X(k)-X(k+1)) )**2 * (xi-X(k))
        else
            alpha(k)=0
            beta(k)=0
        endif
    enddo
    result=dot_product(Y(0:3),alpha(0:3)) + dot_product(m(0:3),beta(0:3))
    write(*,*) "San-Zhuan-Jiao : ",result

    !Calculate the S(i)
    !S(x)=M(k)*(X(k+1)-xi)**3/(6*h(k)) + M(k+1)*(xi-X(k))**3/h(k) + (Y(k)-M(k)*h(k)**2/6)*(X(k+1)-xi)/h(k) + (Y(k+1)-M(k+1)*h(k)**2/6)*(xi-X(k))/h(k)
    do k=0,2
        if (X(k)<=xi .and. xi <X(k+1)) then
            !The equation is too big to complie ,so express it step by step
            result1=M(k) * ((X(k+1)-xi)**3.0) / (6.0*h(k))
            result1=result1 + M(k+1) * ((xi-X(k))**3.0) / (6.0*h(k))
            result1=result1 + (Y(k)-M(k) * (h(k)**2.0) / 6.0) * (X(k+1)-xi)/h(k)
            result1=result1 + (Y(k+1)-M(k+1) * (h(k)**2.0) / 6.0) * (xi-X(k))/h(k)
            exit
        endif
    enddo
    write(*,*) "San-Wan-Ju",result1
    write(*,*) "Succeed ."
    Cubic=result
    !Cubic=result1
    return
end function Cubic