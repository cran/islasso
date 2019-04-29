subroutine standardize(X,xm,xse,n,p,intercept)
implicit none
integer :: n,p,intercept
double precision :: X(n,p),xm(p),xse(p)
integer :: j,start
xm = 0.d0
xse = 1.d0
start = 1
if(intercept.eq.1) start = 2
do j = start, p
    xm(j) = sum(X(:, j)) / n
    X(:, j) = X(:, j) - xm(j)
    xse(j) = sqrt(dot_product(X(:, j), X(:, j)) / n)
    X(:, j) = X(:, j) / xse(j)
end do
end subroutine standardize
