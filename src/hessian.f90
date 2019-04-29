subroutine hessian(theta,se,lambda,xtx,pi,p,hess)
implicit none
integer :: p,i
double precision :: theta(p),se(p),lambda(p),xtx(p,p),pi(p),hess(p,p)
double precision :: dnm,dnm2,temp1(p)
hess = xtx
do i = 1, p
    temp1(i) = theta(i)/se(i)
    hess(i,i) = hess(i,i) + 2.d0 * lambda(i) * ( pi(i) * dnm(temp1(i)) + (1.d0 - pi(i)) * dnm2(temp1(i)) ) / se(i)
end do
end subroutine hessian
