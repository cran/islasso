subroutine standardize(X,xm,xse,n,p,intercept)
implicit none
integer :: n,p,intercept
double precision :: X(n,p),xm(p),xse(p)
integer :: j,start
xm = 0.d0
xse = 1.d0
start = 1
if((intercept.eq.1).and.(p.gt.1)) start = 2
do j = start, p
    !if(intercept.eq.1) xm(j) = sum(X(:, j)) / n
    xm(j) = sum(X(:, j)) / n
    xse(j) = sqrt(dot_product(X(:, j), X(:, j)) / n - xm(j)**2)
    X(:, j) = (X(:, j) - xm(j)) / xse(j)
end do
end subroutine standardize

subroutine check_out(theta,cov,xm,xse,p,intercept)
implicit none
integer :: p,intercept
double precision :: theta(p),cov(p,p),xm(p),xse(p)
integer :: i,j
theta = theta / xse
if((intercept.eq.1).and.(p.gt.1)) theta(1) = theta(1) - dot_product(theta(2:p), xm(2:p))
do i = 1, p
    cov(i,i) = cov(i,i) / (xse(i) * xse(i))
    do j = (i+1), p
        cov(i,j) = cov(i,j) / (xse(i) * xse(j))
        cov(j,i) = cov(i,j)
    end do
end do
if((intercept.eq.1).and.(p.gt.1)) then
    cov(1,:) = cov(1,:) - MATMUL(xm(2:p), cov(2:p,:))
    cov(:,1) = cov(1,:)
    cov(1,1) = cov(1,1) - dot_product(cov(1,2:p), xm(2:p))
end if
end subroutine check_out
