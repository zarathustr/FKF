function [xk , Pk ] = kalman_update( xk_1,yk,Pk_1,Phi_k,Xi_k,Eps_k )

    x_=Phi_k*xk_1;
    Pk_=Phi_k*Pk_1*Phi_k'+Xi_k;
    Gk=Pk_*(inv(Pk_+Eps_k));
    Pk=(eye(4)-Gk)*Pk_;
    xk=x_+Gk*(yk-x_);
    
end

