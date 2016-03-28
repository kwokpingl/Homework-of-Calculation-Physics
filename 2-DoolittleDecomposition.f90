program main

        real*8 X(9) ,Y(9) ,A(9,10) ,L(9,9) ,U(9,9) ,SIZE ,Xtemp(9) ,flag ,loopcount
        data X/9*5/                !Initalize a Xtemp for loop
        data Xtemp/9*5/       !Initalize a Xtemp for loop
        SIZE=0.0001               !Set the accuracy
        loopcount=0               !Use it to count the loop-times
        write(*,*) "Doolittle Decomposition Method : "
        !Read matrix A ,treat A(i,10) as B(i)
        open(unit=55,file='Matrix.txt',status='old')
        do i=1,9
                read(55,*) A(i,1:10)
        enddo
        close(55)
        !Build matrix L & matrix U from matrix A
        do i=1,9
                L(i,i)=1
                do j=1,9
                        if (j>=i) then    !Upper triangular matrix U
                                U(i,j)=A(i,j)-dot_product(L(i,1:i-1),U(1:i-1,j))
                        else                  !Lowwer triangular matrix L
                                L(i,j)=( A(i,j)-dot_product(L(i,1:j-1),U(1:j-1,j)) )/U(j,j)
                        endif
                enddo
        enddo
        !Slove Equations
        do while( flag /= 9 )
                loopcount = loopcount+1
                write(*,*)
                write(*,*) "--------","Loop-times :",loopcount ,"--------"
                !Start a new loop for roots
                do i=1,9
                        Y(i)=A(i,10)-dot_product(L(i,1:i-1),Y(1:i-1))               !Slove LY=b
                        X(i)=( Y(i)-dot_product(U(i,i+1:9),X(i+1:9)) )/U(i,i)   !Slove UX=Y
                        write(*,*) "Y"," = ",Y(i),"X"," = ",X(i)
                enddo
                !Count the number of roots whose accuracy is in SIZE
                flag = 0              !Use as a flag  for showing if the accuracy is in SIZE
                do i=1,9
                        if ( abs(X(i)-Xtemp(i)) < SIZE ) then
                                flag = flag +1
                        endif
                enddo
                Xtemp = X         !Update the Xtemp
                write(*,*) "--------","Flag number :" ,flag ,"--------"
        enddo
end program main
