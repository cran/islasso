subroutine family(fam,link,func,x,n,y)
implicit none
integer :: fam,link,func,n
double precision :: x(n),y(n)
select case(fam)
    case(1)                                                 ! 1) Binomial ( link )
        select case(link)
            case(1)                                         ! 1.1) link = "logit"
                select case(func)
                    case(1)                                 ! 1 = linkfun
                        call logitlink(x,n,y)
                    case(2)                                 ! 2 = linkinv
                        call logitlinkinv(x,n,y)
                    case(3)                                 ! 3 = mu.eta
                        call logitmueta(x,n,y)
                    case(4)                                 ! 4 = variance
                        call binomial_variance(x,n,y)
                end select
            case(2)                                         ! 1.2) link = "probit"
                select case(func)
                    case(1)                                 ! 1 = linkfun
                        call probitlink(x,n,y)
                    case(2)                                 ! 2 = linkinv
                        call probitlinkinv(x,n,y)
                    case(3)                                 ! 3 = mu.eta
                        call probitmueta(x,n,y)
                    case(4)
                        call binomial_variance(x,n,y)
                end select
        end select
    case(2)                                                 ! 2) Poisson ( link )
        select case(link)
            case(1)                                         ! 2.1) link = "log"
                select case(func)
                    case(1)                                 ! 1 = linkfun
                        call loglink(x,n,y)
                    case(2)                                 ! 2 = linkinv
                        call loglinkinv(x,n,y)
                    case(3)                                 ! 3 = mu.eta
                        call logmueta(x,n,y)
                    case(4)                                 ! 4 = variance
                        call poisson_variance(x,n,y)
                end select
        end select
    case(3)                                                 ! 3) Gamma ( link )
        select case(link)
            case(1)                                         ! 3.1) link = "inverse"
                select case(func)
                    case(1)                                 ! 1 = linkfun
                        call inverselink(x,n,y)
                    case(2)                                 ! 2 = linkinv
                        call inverselinkinv(x,n,y)
                    case(3)                                 ! 3 = mu.eta
                        call inversemueta(x,n,y)
                    case(4)                                 ! 4 = variance
                        call gamma_variance(x,n,y)
                end select
            case(2)                                         ! 3.2) link = "log"
                select case(func)
                    case(1)                                 ! 1 = linkfun
                        call loglink(x,n,y)
                    case(2)                                 ! 2 = linkinv
                        call loglinkinv(x,n,y)
                    case(3)                                 ! 3 = mu.eta
                        call logmueta(x,n,y)
                    case(4)                                 ! 4 = variance
                        call gamma_variance(x,n,y)
                end select
            case(3)                                         ! 3.3) link = "identity"
                select case(func)
                    case(1)                                 ! 1 = linkfun
                        call identitylink(x,n,y)
                    case(2)                                 ! 2 = linkinv
                        call identitylinkinv(x,n,y)
                    case(3)                                 ! 3 = mu.eta
                        call identitymueta(x,n,y)
                    case(4)                                 ! 4 = variance
                        call gamma_variance(x,n,y)
                end select
        end select
end select
end subroutine family
