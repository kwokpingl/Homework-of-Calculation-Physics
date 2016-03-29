program main

        real*8 A(9,10) ,X(9)
        write(*,*) "Gauss Elimination Method : "
        !Read matrix A ,treat A(i,10) as B(i)
        open(55,file='Matrix.txt',status='old')
        do i = 1 ,9
                read (55,*) A(i,1:10)
        enddo
        close(55)
        !Build Triangular matrix A
        do i=1,8
                do j=i+1,9
                        A(j,i:10) = A(j,i:10) - A(i,i:10)/A(i,i)*A(j,i)
                enddo
        enddo
        !Slove Equations
        do i = 9 ,1 ,-1
                X(i) = (A(i,10)-dot_product(A(i, i+1:9),X(i+1:9)))/A(i,i)
                write(*,*) "X" ,i ," = " ,X(i)
        enddo
        !Mode for check
        write(*,*) "    | Dot_product |","              |Matrix B|"
        do i=1,9
                write (*,*) dot_product(A(i,1:9),X(1:9)),A(i,10)
        enddo
end program main