subroutine inv(n, A, invA, info)
integer :: info, n
double precision :: A(n, n), invA(n, n)
integer :: i, j
invA = A
call DPOTRF('U', n, invA, n, info)
!if (info.ne.0) call islasso_trace1_3_2(info)
call DPOTRI('U', n, invA, n, info)
!if (info.ne.0) call islasso_trace1_3_2(info)
do j = 1, n
    do i = j + 1, n
        invA(i, j) = invA(j, i)
    end do
end do
end subroutine inv

subroutine solve(A, b, n, info)
integer :: info, n
double precision :: A(n, n), b(n, 1)
call DPOSV('U', n, 1, A, n, b, n, info)
!if (info.ne.0) call islasso_trace1_4_2(info)
end subroutine solve

subroutine crossp(A, xtx, n, p)
integer :: n, p
double precision :: A(n, p), xtx(p, p)
call DGEMM('T', 'N', p, p, n, 1.d0, A, n, A, n, 0.d0, xtx, p)
end subroutine crossp

subroutine tcrossp(x, txx, n, p)
integer :: n,p
double precision :: x(n,p), txx(n,n)
call DGEMM('N', 'T', n, n, p, 1.d0, x, n, x, n, 0.d0, txx, n)
!do i = 1, n
!    txx(i,i) = dot_product(x(i,:),x(i,:))
!    do j = (i+1), n
!        txx(i,j) = dot_product(x(i,:),x(j,:))
!        txx(j,i) = txx(i,j)
!    end do
!end do
end subroutine tcrossp

!subroutine prod1(x,w,xtw,xtx,n,p)
!integer :: n,p
!double precision :: x(n,p),w(n),xtw(p,n),xtw2(n,p),xtx(p,p)
!integer :: i,j
!do i = 1, p
!    !xtw(i,:) = x(:,i) * w
!    xtw2(:,i) = x(:,i) * w
!    !xtx(i,i) = dot_product(xtw(i,:),x(:,i))
!    xtx(i,i) = dot_product(xtw2(:,i),x(:,i))
!    do j = (i+1), p
!        !xtx(i,j) = dot_product(xtw(i,:),x(:,j))
!        xtx(i,j) = dot_product(xtw2(:,i),x(:,j))
!        xtx(j,i) = xtx(i,j)
!    end do
!end do
!xtw = transpose(xtw2)
!end subroutine prod1

subroutine prod1(x,w,xtw,xtx,n,p)
integer :: n,p
double precision :: x(n,p),w(n),xtw(p,n),xtx(p,p)
integer :: i
do i = 1, p
    xtw(i,:) = x(:,i) * w
end do
!xtw = transpose(xtw2)
call DGEMM('N', 'N', p, p, n, 1.d0, xtw, p, x, n, 0.d0, xtx, p)
end subroutine prod1

!subroutine prod2(xtx,tempMat,invH,cov1,hi,p)
!integer :: p
!double precision :: xtx(p,p),tempMat(p,p),invH(p,p),cov1(p,p),hi(p)
!double precision :: tempMat2(p,p),invH2(p,p)
!integer :: i,j,k
!invH2 = transpose(invH)
!do i = 1, p
!    do j = 1, p
!        !tempMat(i,j) = dot_product(invH(i,:),xtx(:,j))
!        tempMat2(j,i) = dot_product(invH2(:,i),xtx(:,j))
!        !tempMat(j,i) = tempMat(i,j)
!    end do
!    hi(i) = tempMat2(i,i)
!    cov1(i,i) = dot_product(tempMat2(:,i),invH(:,i))
!    do k = (i+1), p
!    cov1(i,k) = dot_product(tempMat2(:,i),invH(:,k))
!        cov1(k,i) = cov1(i,k)
!    end do
!end do
!tempMat = transpose(tempMat2)
!end subroutine prod2

subroutine armijo(beta, se, dir, f0, alpha, h, x, y, offset, n, p, lambda, &
& eta, res, pi, xtw)
integer :: n, p
double precision :: beta(p), se(p), dir(p), f0, alpha, h, x(n, p), y(n), offset(n)
double precision :: lambda(p), eta(n), res(n), pi(p), xtw(p, n)
double precision :: fn, betan(p), fac, grad(p)
h = 1.d0
fac = 0.75d0
eta = 0.d0
betan = beta - h * dir
call DGEMV('N', n, p, 1.d0, x, n, betan, 1, 0.d0, eta, 1)
eta = eta + offset
res = y - eta
call gradient(betan, se, lambda, xtw, res, pi, n, p, grad, alpha)
fn = sqrt(sum(grad**2))
do while (fn.gt.f0)
    h = fac * h
    betan = beta - h * dir
    call DGEMV('N', n, p, 1.d0, x, n, betan, 1, 0.d0, eta, 1)
    eta = eta + offset
    res = y - eta
    call gradient(betan, se, lambda, xtw, res, pi, n, p, grad, alpha)
    fn = sqrt(sum(grad**2))
    if(h.le.0.00000001d0) exit
end do
beta = betan
f0 = fn
end subroutine armijo

subroutine armijo2(beta, se, cov, cov1, s2, f0, alpha, h, n, p, lambda, res, pi, xtw)
integer :: n, p
double precision :: beta(p), se(p), cov(p,p), cov1(p,p), s2, f0, alpha, h
double precision :: lambda(p), res(n), pi(p), xtw(p,n)
double precision :: fac, covn(p,p), sen(p), grad(p), fn
integer :: j, err
h = 1.d0
fac = 0.75d0
covn = cov + h * (cov1 - cov)
err = 0
do j = 1, p
    sen(j) = sqrt(s2 * covn(j,j))
    if((s2 * covn(j,j)).lt.0) err = 1
end do
if(err.eq.1) then
    sen = se
    fn = f0
else
    call gradient(beta, sen, lambda, xtw, res, pi, n, p, grad, alpha)
    fn = sqrt(sum(grad**2))
end if
do while ((fn.gt.f0).or.(err.eq.1))
    h = fac * h
    covn = cov + h * (cov1 - cov)
    err = 0
    do j = 1, p
        sen(j) = sqrt(s2 * covn(j,j))
        if((s2 * covn(j,j)).lt.0) err = 1
    end do
    if(err.eq.1) then
        sen = se
        fn = f0
    else
        call gradient(beta, sen, lambda, xtw, res, pi, n, p, grad, alpha)
        fn = sqrt(sum(grad**2))
    end if
    if(h.le.0.00000001d0) exit
end do
cov = covn
se = sen
f0 = fn
end subroutine armijo2



subroutine prod2(xtx,tempMat,invH,cov1,hi,p)
integer :: p
double precision :: xtx(p,p),tempMat(p,p),invH(p,p),cov1(p,p),hi(p)
integer :: i
call DGEMM('N', 'N', p, p, p, 1.d0, invH, p, xtx, p, 0.d0, tempMat, p)
call DGEMM('N', 'N', p, p, p, 1.d0, tempMat, p, invH, p, 0.d0, cov1, p)
do i = 1, p
    hi(i) = tempMat(i,i)
end do
end subroutine prod2


!subroutine linear_predictor(x,beta,eta,offset,n,p)
!integer :: n,p
!double precision :: x(n,p),beta(p),eta(n),offset(n)
!integer :: i
!eta = offset
!eta = 0.d0
!do i = 1, n
!    do j = 1, p
!        eta(i) = eta(i) + sum(x(i,:) * beta)
!    end do
!end do
!end subroutine linear_predictor

subroutine linear_predictor(x,beta,eta,offset,n,p)
integer :: n,p
double precision :: x(n,p),beta(p),eta(n),offset(n)
eta = 0.d0
call DGEMV('N', n, p, 1.d0, x, n, beta, 1, 0.d0, eta, 1)
eta = eta + offset
end subroutine linear_predictor

subroutine setdiff(p,iprofile,ind_noprofile)
integer :: p,iprofile,ind_noprofile(p-1)
integer :: i,k
k = 0
do i = 1, p
    if(i.ne.iprofile) then
        k = k + 1
        ind_noprofile(k) = i
    end if
end do
end subroutine setdiff
