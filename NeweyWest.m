function NWSE = NeweyWest(Z)
% Newey-West standard errors using plug-in procedure

[T,k] = size(Z);

L = floor(4*((T/100)^(2/9)));

S = 0;
for l = 0:L
    w_l = 1-l/(L+1);
    for t = l+1:T
        if (l==0)   
            S = S  + Z(t, :)' * Z(t,:);
        else        
            S = S + w_l *...
                (Z(t, :)' * Z(t-l,:) + Z(t-l, :)' * Z(t,:));
        end
    end
end
S = (1/(T-k)) .*S;

NWSE = sqrt(diag(T.*((Z'*Z)\S/(Z'*Z))));

end
