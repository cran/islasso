subroutine gradient(theta,se,lambda,xtw,res,pi,n,p,grad)
implicit none
integer :: n,p,i
double precision :: theta(p),se(p),lambda(p),xtw(p,n),res(n),pi(p),grad(p)
double precision :: pnm,pnm2,temp1(p),xtw2(n,p)
xtw2 = transpose(xtw)
do i = 1, p
    !grad(i) = - dot_product(xtw2(i,:),res)
    grad(i) = - dot_product(xtw2(:,i),res)
    temp1(i) = theta(i)/se(i)
    grad(i) = grad(i) + lambda(i) * ( pi(i) * (2.d0 * pnm(temp1(i)) - 1.d0) + (1.d0 - pi(i)) * (2.d0 * pnm2(temp1(i)) - 1.d0) )
end do
end subroutine gradient
