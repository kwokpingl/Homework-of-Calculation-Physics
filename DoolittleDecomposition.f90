program main

        real*8 X(9) ,Y(9) ,A(9,10) ,L(9,9) ,U(9,9)
        write(*,*) "Doolittle Decomposition Method : "
        !Read matrix A ,treat A(i,10) as B(i)
        open(unit=55,file='Matrix.txt',status='old')
        do i=1,9
                read(55,*) A(i,1:10)
        enddo
        close(55)
        !Build matrix L & matrix U from matrix A
        do k=1,9
                L(k,k)=1                !Initalize matrix L
                U(k,k)=A(k,k)        !Initalize matrix U
                do j=k,9
                        i=j+1
                        U (k,j)=A(k,j)-dot_product(L(k,1:k-1),U(1:k-1,j))                               !Build upper matrix U
                        L (i,k)=( A(i,k)-dot_product(L(i,1:k-1),U(1:k-1,k)) )/U(k,k)               !Build lowwer matrix L
                enddo
        enddo
        !Slove Euqations
        do i=1,9
                Y(i)=A(i,10)-dot_product(L(i,1:i-1),Y(1:i-1))               !Slove LY=b
        enddo
        do i=9,1,-1
                X(i)=( Y(i)-dot_product(U(i,i+1:9),X(i+1:9)) )/U(i,i)   !Slove UX=Y
                write(*,*) "X",i," = ",X(i)
        enddo
        !Mode for check
        write(*,*) "    | Dot_product |","              |Matrix B|"
        do i=1,9
                write (*,*) dot_product(A(i,1:9),X(1:9)),A(i,10)
        enddo
end program main