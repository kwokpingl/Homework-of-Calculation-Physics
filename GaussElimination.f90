program main

        real*8 A(9,10) ,X(9) ,Xtemp(9) ,SIZE ,flag ,loopcount
        real*8 sumTemp
        data X/9*5/                !Initalize a X for loop
        data Xtemp/9*5/       !Initalize a Xtemp for loop
        SIZE=0.0001               !Set the accuracy
        loopcount=0               !Use it to count the loop-times
        write(*,*) "Gauss Elimination Method : "
        !Read matrix A ,treat A(i,10) as B(i)
        open(55,file='Matrix.txt',status='old')
        do i = 1 ,9
                read (55,*) A(i,1:10)
        enddo
        close(55)
        !Slove Equations
        do while( flag /= 9 )
                loopcount = loopcount+1
                flag = 0               !Use as a flag  to Count the number of roots whose accuracy is in SIZE
                write(*,*)
                write(*,*) "--------","Loop-times :",loopcount ,"--------"
                !Start a new loop for roots
                do i = 1 ,9
                        sumTemp = dot_product(A(i, i+1:9),X(i+1:9))
                        X(i) = (A(i,10)-sumTemp)/A(i,i)
                        if ( abs(X(i)-Xtemp(i)) < SIZE ) then
                                flag = flag +1          !If the number of roots in accuracy is satisfied ,flag is N ,then success
                        endif
                        write(*,*) "X" ,i ," = " ,X(i)
                enddo
                Xtemp = X         !Update the Xtemp
                write(*,*) "--------","Flag number :" ,flag ,"--------"
        enddo
end program main