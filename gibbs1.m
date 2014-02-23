%=========================================================================
%                  GIBBS SAMPLER UNDER PRIOR 1
%=========================================================================

ttheta1dp1 = zeros(nsim,1);
ttheta2dp1 = zeros(nsim,1);




for i = 2:nsim
    
    % Draw from conditional density of theta1|theta2
    
    yvector_t1 = yvector - ttheta2dp1(i-1)*xvector(:,2);
    
    ttheta1var     = (xvector(:,1)'*xvector(:,1)+(1/ttau))^(-1);
    ttheta1post    = (ttheta1var)*(xvector(:,1)'*yvector_t1);
    
    ttheta1dp1(i)    = normrnd(ttheta1post,ttheta1var);
    
    
    % Draw from conditional density of theta2|theta1
    
    yvector_t2 = yvector - ttheta1dp1(i-1)*xvector(:,1);
    
    ttheta2var  = (xvector(:,2)'*xvector(:,2)+...
        (ttheta1dp1(i-1)/ttau))^(-1);
    ttheta2post = (ttheta2var)*(xvector(:,2)'*yvector_t2);
    
    ttheta2dp1(i) = normrnd(ttheta2post,ttheta2var);
    
end


