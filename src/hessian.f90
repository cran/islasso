subroutine hessian(theta,se,lambda,xtx,pi,p,hess,alpha)
implicit none
integer :: p,i
double precision :: theta(p),se(p),lambda(p),xtx(p,p),pi(p),hess(p,p)
double precision :: dnm,temp1(p),alpha
hess = xtx
do i = 1, p
    temp1(i) = theta(i)/se(i)
    hess(i,i) = hess(i,i) + 2.d0 * lambda(i) * alpha * ( pi(i) * dnm(temp1(i), 0.d0, 1.d0) + &
& (1.d0 - pi(i)) * dnm(temp1(i), 0.d0, 0.01d0) ) / se(i) + (1.d0 - alpha) * lambda(i)
end do
end subroutine hessian
