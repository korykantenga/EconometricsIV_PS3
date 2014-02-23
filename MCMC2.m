
%==========================================================================
%                    MONTE CARLO EXPERIMENT
%==========================================================================

msim    = 100;
nsim    = 1000;
nburnp  = 0.4;

mcmcp1_a1 = zeros(msim,2);
mcmcp2_a1 = zeros(msim,2);
mcmcp1_a0 = zeros(msim,2);
mcmcp2_a0 = zeros(msim,2);

yvector = yvector_a1; %#ok<NASGU>

for e=1:msim
    
    gibbs1
    gibbs2
    
    
    % mcmcp Draws under Prior 1
    fg = ttheta2dp1./ttheta1dp1;
    fg(isnan(fg)) = [];
    fg2= fg.^2;
    mcmcp1_a1(e,:) = [mean(fg(round(nburnp*length(fg)):length(fg),1)),...
        mean(fg2(round(nburnp*length(fg)):length(fg)))];
    
    if e==1
        acfb1  = fg(round(nburnp*length(fg)):length(fg));
        acf2b1 = fg2(round(nburnp*length(fg)):length(fg));
        ttheta111 = ttheta1dp1;
        ttheta211 = ttheta2dp1;
    end
    
    % mcmcp Draws under Prior 2
    fg = ttheta2dp2./ttheta1dp2;
    fg(isnan(fg)) = [];
    fg2= fg.^2;
    mcmcp2_a1(e,:) = [mean(fg(round(nburnp*length(fg)):length(fg))),...
        mean(fg2(round(nburnp*length(fg)):length(fg)))];
    
    if e==1
        acfb2  = fg(round(nburnp*length(fg)):length(fg));
        acf2b2 = fg2(round(nburnp*length(fg)):length(fg));
        ttheta121 = ttheta1dp2;
        ttheta221 = ttheta2dp2;
    end
    
end

yvector = yvector_a0;

for e=1:msim
    
    gibbs1
    gibbs2
    
    
    % mcmcp Draws under Prior 1
    fg = ttheta2dp1./ttheta1dp1;
    fg(isnan(fg)) = [];
    fg2= fg.^2;
    mcmcp1_a0(e,:) = [mean(fg(round(nburnp*length(fg)):length(fg))),...
        mean(fg2(round(nburnp*length(fg)):length(fg)))];
    
    if e==1
        acfb3  = fg(round(nburnp*length(fg)):length(fg));
        acf2b3 = fg2(round(nburnp*length(fg)):length(fg));
        ttheta110 = ttheta1dp1;
        ttheta210 = ttheta2dp1;
    end
    
    % mcmcp Draws under Prior 2
    fg = ttheta2dp2./ttheta1dp2;
    fg(isnan(fg)) = [];
    fg2= fg.^2;
    mcmcp2_a0(e,:) = [mean(fg(round(nburnp*length(fg)):length(fg))),...
        mean(fg2(round(nburnp*length(fg)):length(fg)))];
    
    if e==1
        acfb4  = fg(round(nburnp*length(fg)):length(fg));
        acf2b4 = fg2(round(nburnp*length(fg)):length(fg));
        ttheta120 = ttheta1dp2;
        ttheta220 = ttheta2dp2;
    end
    
end



%==========================================================================
%                   MONTE CARLO VARIANCE
%==========================================================================

var(mcmcp1_a0)
var(mcmcp1_a1)
var(mcmcp2_a0)
var(mcmcp2_a1)

%==========================================================================
%                    PLOT ACFS OF FIRST CHAIN
%==========================================================================



figure;
subplot(2,2,1)
autocorr(acfb1)
title('Prior 1, \alpha = 0.9')
subplot(2,2,2)
autocorr(acfb2)
title('Prior 2, \alpha = 0.9')
subplot(2,2,3)
autocorr(acfb3)
title('Prior 1, \alpha = 0.001')
subplot(2,2,4)
autocorr(acfb4)
title('Prior 2, \alpha = 0.001')
suptitle('E[\beta|Y]')
print -depsc2 acfTAU2.eps

figure;
subplot(2,2,1)
autocorr(acf2b1)
title('Prior 1, \alpha = 0.9')
subplot(2,2,2)
autocorr(acf2b2)
title('Prior 2, \alpha = 0.9')
subplot(2,2,3)
autocorr(acf2b3)
title('Prior 1, \alpha = 0.001')
subplot(2,2,4)
autocorr(acf2b4)
title('Prior 2, \alpha = 0.001')
suptitle('E[\beta^2|Y]')
print -depsc2 acf2TAU2.eps



%==========================================================================
%                    NEWEY WEST STANDARD ERRORS
%==========================================================================

NeweyWest(acfb1)
NeweyWest(acfb2)
NeweyWest(acfb3)
NeweyWest(acfb4)

NeweyWest(acf2b1)
NeweyWest(acf2b2)
NeweyWest(acf2b3)
NeweyWest(acf2b4)

%==========================================================================
%                    MONTE CARLO STANDARD ERRORS
%==========================================================================

v1 = sqrt(mcmcp1_a0(1,2)-(mcmcp1_a0(1,1)^2));
v2 = sqrt(mcmcp2_a0(1,2)-(mcmcp2_a0(1,1)^2));
v3 =sqrt(mcmcp1_a1(1,2)-(mcmcp1_a1(1,1)^2));
v4 = sqrt(mcmcp2_a1(1,2)-(mcmcp2_a1(1,1)^2));

sqrt(2*v1*(2*(mcmcp1_a0(1,1))^2+v1))
sqrt(2*v2*(2*(mcmcp2_a0(1,1))^2+v2))
sqrt(2*v3*(2*mcmcp1_a1(1,1)^2+v3))
sqrt(2*v4*(2*mcmcp2_a1(1,1)^2+v4))



