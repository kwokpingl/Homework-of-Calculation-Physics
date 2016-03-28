program main
        real*8 A(9,10) ,X(9) ,Xtemp(9) ,w ,SIZE ,flag ,loopcount
        data X/9*5/                !Initalize a Xtemp for loop
        data Xtemp/9*5/       !Initalize a Xtemp for loop
        SIZE=0.0001               !Set the accuracy
        loopcount=0               !Use it to count the loop-times
        w=1.55                        !Set the relaxational index
        write(*,*) "Overrelaxation Method : "
        !Read matrix A ,treat A(i,10) as B(i)
        open(55,file='Matrix.txt')
        do i=1,9
                read(55,*) A(i,1:10)
        enddo
        close(55)
        !Slove the Equations
        do while( flag /= 9 )
                loopcount = loopcount+1
                write(*,*)
                write(*,*) "--------","Loop-times :",loopcount ,"--------"
                !Start a new loop for roots
                do i=1,9
                        X(i)=Xtemp(i)+w*(A(i,10)-dot_product(A(i,1:i-1),X(1:i-1))-dot_product(A(i,i:9),Xtemp(i:9)))/A(i,i)
                        write(*,*) "X"," = ",X(i)
                enddo
                !Count the number of roots whose accuracy is in SIZE
                flag = 0          !Use as a flag  for showing if the accuracy is in SIZE
                do i=1,9
                        if ( abs(X(i)-Xtemp(i)) < SIZE ) then
                                flag = flag +1
                        endif
                enddo
                Xtemp = X         !Update the Xtemp
                write(*,*) "--------","Flag number :" ,flag ,"--------"
        enddo
end program main
