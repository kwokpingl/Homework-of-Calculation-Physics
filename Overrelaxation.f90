program main
        real*8 A(9,10) ,X(9) ,Xtemp(9) ,w ,SIZE ,flag ,loopcount
        data X/9*5/                !Initalize a X for loop
        data Xtemp/9*5/       !Initalize a Xtemp for loop
        SIZE=0.000000001               !Set the accuracy
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
                flag = 0               !Use as a flag  to Count the number of roots whose accuracy is in SIZE
                write(*,*)
                write(*,*) "--------","Loop-times :",loopcount ,"--------"
                !Start a new loop for roots
                do i=1,9
                        X(i)=Xtemp(i)+w*(A(i,10)-dot_product(A(i,1:i-1),X(1:i-1))-dot_product(A(i,i:9),Xtemp(i:9)))/A(i,i)
                        if ( abs(X(i)-Xtemp(i)) < SIZE ) then
                                flag = flag +1          !If the number of roots in accuracy is satisfied ,flag is N ,then success
                        endif
                        write(*,*) "X",i," = ",X(i)
                enddo
                Xtemp = X         !Update the Xtemp
                write(*,*) "--------","Flag number :" ,flag ,"--------"
                write(*,*) "----The Accuracy is : ",SIZE,"-------"
        enddo
        !Mode for check
        write(*,*) "    | Dot_product |","            | Matrix B| "
        do i=1,9
                write (*,*) dot_product(A(i,1:9),X(1:9)),A(i,10)
        enddo
end program main