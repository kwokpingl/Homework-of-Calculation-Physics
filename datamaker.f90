program main
    real*8 x,y
    x=0
    open(1,file='dotKnew.txt')
    do i=0,15
        x=-5.0+10.0/15.0*i
        y=1/(1+x*x)
        write(1,*)  x ,y
    enddo
    close(1)
end program main