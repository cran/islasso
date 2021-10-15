subroutine deviance(weights,res,lambda,alpha,theta,n,p,dev)
implicit none
integer :: n,p
double precision :: weights(n),res(n),lambda(p),alpha,theta(p),dev
dev = sum(weights * (res**2)) + sum(lambda * (alpha * abs(theta) + 0.5d0 * (1 - alpha) * (theta**2)))
end subroutine deviance

!subroutine gradient(theta,se,lambda,xtw,res,pi,n,p,grad,alpha)
!implicit none
!integer :: n,p,i
!double precision :: theta(p),se(p),lambda(p),xtw(p,n),res(n),pi(p),grad(p)
!double precision :: pnm,temp1,alpha,k
!!xtw2 = transpose(xtw)
!grad = 0.d0
!call DGEMV('N', p, n, 1.d0, xtw, p, res, 1, 0.d0, grad, 1)
!grad = -grad
!k = 0.d0
!if (alpha.eq.1.0d0) k = 0.05d0 / n
!do i = 1, p
!    !grad(i) = - dot_product(xtw2(i,:),res)
!!    grad(i) = -dot_product(xtw2(:,i), res)
!    temp1 = theta(i) / se(i)
!    grad(i) = grad(i) + lambda(i) * alpha * ( pi(i) * (2.d0 * pnm(temp1, 0.d0, 1.d0) - 1.d0) + (1.d0 - pi(i)) * &
!        & (2.d0 * pnm(temp1, 0.d0, 0.001d0) - 1.d0) ) + lambda(i) * (1.d0 - alpha) * theta(i) + k * lambda(i) * theta(i)
!end do
!end subroutine gradient

subroutine gradient(theta,se,lambda,xtw,res,pi,n,p,grad,alpha)
implicit none
integer :: n,p,i
double precision :: theta(p),se(p),lambda(p),xtw(p,n),res(n),pi(p),grad(p)
double precision :: pnm,temp1,alpha
!xtw2 = transpose(xtw)
grad = 0.d0
call DGEMV('N', p, n, 1.d0, xtw, p, res, 1, 0.d0, grad, 1)
grad = -grad
do i = 1, p
    temp1 = theta(i) / se(i)
    grad(i) = grad(i) + lambda(i) * alpha * ( pi(i) * (2.d0 * pnm(temp1, 0.d0, 1.d0) - 1.d0) + (1.d0 - pi(i)) * &
        & (2.d0 * pnm(temp1, 0.d0, 0.00001d0) - 1.d0) ) + lambda(i) * (1.d0 - alpha) * theta(i)
end do
end subroutine gradient

!subroutine hessian(theta,se,lambda,xtx,pi,p,n,hess,alpha)
!implicit none
!integer :: p,i,n
!double precision :: theta(p),se(p),lambda(p),xtx(p,p),pi(p),hess(p,p)
!double precision :: dnm,temp1,alpha,k
!hess = xtx
!k = 0.d0
!if (alpha.eq.1.0d0) k = 0.05d0 / n
!do i = 1, p
!    temp1 = theta(i) / se(i)
!    hess(i,i) = hess(i,i) + 2.d0 * lambda(i) * alpha * ( pi(i) * dnm(temp1, 0.d0, 1.d0) + &
!        & (1.d0 - pi(i)) * dnm(temp1, 0.d0, 0.001d0) ) / se(i) + (1.d0 - alpha) * lambda(i) + k * lambda(i)
!end do
!end subroutine hessian

subroutine hessian(theta,se,lambda,xtx,pi,p,hess,alpha)
implicit none
integer :: p,i
double precision :: theta(p),se(p),lambda(p),xtx(p,p),pi(p),hess(p,p)
double precision :: dnm,temp1,alpha
hess = xtx
do i = 1, p
    temp1 = theta(i) / se(i)
    hess(i,i) = hess(i,i) + 2.d0 * lambda(i) * alpha * ( pi(i) * dnm(temp1, 0.d0, 1.d0) + &
        & (1.d0 - pi(i)) * dnm(temp1, 0.d0, 0.00001d0) ) / se(i) + (1.d0 - alpha) * lambda(i)
end do
end subroutine hessian
