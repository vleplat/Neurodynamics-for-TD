function [Xnext, Vnext] = PSO_Code(V, Xbar, Xp, Xg)
varrho = 0.1 ;
eta1 = 0.2 ;
eta2 = 0.2 ;
e1 = rand(1) ;
e2 = rand(1) ;

Vnext =  varrho*V + eta1*e1*(Xp-Xbar) + eta2*e2*(Xg-Xbar) ;
Xnext = Xbar + Vnext ;