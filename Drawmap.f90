program main
    implicit none
    integer x ,y ,i ,j, tempX ,tempY
    real*8 Map(0:1800/20 , 0:1500/20) ,tempH
    real*8,external::Cubic
    ! Initialize the Map data
    Map=0
    open(55,file='grids.txt')
    do
        read(55,*) tempX ,tempY ,tempH
        if (tempX == -9999 ) then        ! Define -9999 as EOF
            exit
        endif
        Map(tempX/20,tempY/20)=tempH
    enddo
    close(55)
    !Calculate the Map
    do x=0,1800/20
        ! Get the X-axis value 
        ! Get f(x,0)
        Map(x,0)=Cubic( x*20d0,0d0,600d0,1200d0,1800d0,&
                                &Map(0,0),Map(600/20,0),Map(1200/20,0),Map(1800/20,0))
        ! Get f(x,500)
        Map(x,500/20)=Cubic( x*20d0,0d0,600d0,1200d0,1800d0,&
                                        &Map(0,500/20),Map(600/20,500/20),Map(1200/20,500/20),Map(1800/20,500/20))
        ! Get f(x,1000)
        Map(x,1000/20)=Cubic( x*20d0,0d0,600d0,1200d0,1800d0,&
                                          &Map(0,1000/20),Map(600/20,1000/20),Map(1200/20,1000/20),Map(1800/20,1000/20))
        ! Get f(x,1500)
        Map(x,1500/20)=Cubic( x*20d0,0d0,600d0,1200d0,1800d0,&
                                          &Map(0,1500/20),Map(600/20,1500/20),Map(1200/20,1500/20),Map(1800/20,1500/20))
        
        ! Get the Y-axis value
        do y=0,1500/20
            ! Get f(x,y)
            Map(x,y)=Cubic( y*20d0,0d0,500d0,1000d0,1500d0,Map(x,0),Map(x,500/20),Map(x,1000/20),Map(x,1500/20))
            write(*,*) x,y,Map(x,y)
        enddo
    enddo

    ! Write the result
    write(*,*) "------------------------- All Calculation Finished -------------------------"
    write(*,*) "Write the result into 'Map.txt'..."
    open(1,file='Map.txt')
    do i=0,1800/20
        do j=0,1500/20
            !print *,i*20 ,j*20 ,Map(i,j)
            write(1,*) i*20 ,j*20 ,Map(i,j)
        enddo
        write(1,*)      ! Set a blank lines between each X-axis data ,in order to fomat to plot gnuplot pm3d
    enddo
    close(1)
    write(*,*) "Succeed."
end program main



real*8 function Cubic(xi,x0,x1,x2,x3,num0,num1,num2,num3)
    implicit none
    integer i,j,k
    real*8 :: num0 ,num1 ,num2 ,num3,x0,x1,x2,x3
    real*8 Y(0:3) ,X(0:3) ,m(0:3) ,M1(0:3) ,h(0:2)
    real*8 a(0:2) ,b(2) ,a1(0:2) ,b1(2)     !Define a(0) to default ,which is of no use in calculation ,while calculating Q(i,j) for QM=b
    real*8 Q(2,0:3) ,Q1(2,0:3)                 !Define Q(i,0) & Q(i,3) to default ,in case of over boundary
    real*8 xi ,result ,result1 ,alpha(0:3) ,beta(0:3)
    m(0)    = 0
    m(3)    = 0
    M1(0)  = 0
    M1(3)  = 0
    !Read data knew
    X=(/x0,x1,x2,x3/)                       !Read X ,Y
    Y=(/num0,num1,num2,num3/)
    !Calculate the h(i)
    do i=0,2
        h(i) = X(i+1)-X(i)
    enddo
    !Calculate the a(i) & b(i)
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
    do i = 2 ,1 ,-1
            m(i)   = (b(i)  -(Q(i, i+1)   *m(i+1)))  /Q(i,i)
            M1(i) = (b1(i)-(Q1(i, i+1)*M1(i+1)))/Q1(i,i)
    enddo

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
    Cubic=result        ! Use the first derivatives
    !Cubic=result1     ! Use the second derivatives
    return
end function Cubic