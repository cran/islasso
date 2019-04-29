subroutine check_out(theta,cov,xm,xse,p,intercept)
implicit none
integer :: p,intercept
double precision :: theta(p),cov(p,p),xm(p),xse(p)
integer :: i,j
theta = theta / xse
if(intercept.eq.1) theta(1) = theta(1) - dot_product(theta(2:p),xm(2:p))
do i = 1, p
    cov(i,i) = cov(i,i) / (xse(i) * xse(i))
    do j = (i+1), p
        cov(i,j) = cov(i,j) / (xse(i) * xse(j))
        cov(j,i) = cov(i,j)
    end do
end do
if(intercept.eq.1) then
    cov(1,:) = cov(1,:) - MATMUL(xm(2:p), cov(2:p,:))
    cov(:,1) = cov(1,:)
    cov(1,1) = cov(1,1) - dot_product(cov(1,2:p), xm(2:p))
end if
end subroutine check_out
