%==========================================================================
%                 Local Identification Failure
%
% Author: Kory Kantenga
% Date: 2/19/2014
%
% Based on VEC Matlab Code obtained from Frank Schorfheide written by
% Luigi Bocola.
%==========================================================================

%==========================================================================
%                        HOUSEKEEPING
%==========================================================================

tic
close all
clear all
clc


%==========================================================================
%                        CONTROL PARAMETERS
%==========================================================================

gendata  = 1;
plotdata = 1;
plotllhd = 1;
plotdens = 1;
gibbsrun = 1;

%==========================================================================
%                        DATA PARAMETERS
%==========================================================================

rrho1 = 0.2;
rrho2 = 0.2;

x1_0  = 0  ;
x2_0  = 0  ;

aalpha0 = 0.001;
aalpha1 = 0.9  ;

bbeta   = 0.5;


%=========================================================================
%                           PRIOR PARAMETERS
%=========================================================================

ttau = 0.1;

%==========================================================================
%                        GENERATE THE DATA
%=========================================================================

if gendata==1
    
    T = 100;
    
    xerrors = randn(100,2);
    yerrors = randn(100,1);
    
    xvector    = zeros(T,2);
    yvector_a0 = zeros(T,1);
    yvector_a1 = zeros(T,1);
    
    %Initialise at x_0
    xvector(1,1) = rrho1*x1_0 + xerrors(1,1);
    xvector(1,2) = rrho2*x2_0 + xerrors(1,2);
    yvector_a0(1)   = aalpha0*xvector(1,1) + aalpha0*bbeta*xvector(1,2)+...
        + yerrors(1);
    yvector_a1(1)   = aalpha1*xvector(1,1) + aalpha1*bbeta*xvector(1,2)+...
        + yerrors(1);
    
    %Recusively fill in values
    for t = 2:T
        xvector(t,1) = rrho1*xvector(t-1,1) + xerrors(t,1);
        xvector(t,2) = rrho1*xvector(t-1,2) + xerrors(t,2);
        
        yvector_a0(t) = aalpha0*xvector(t,1) + aalpha0*bbeta*xvector(t,2)+...
            + yerrors(t);
        
        yvector_a1(t) = aalpha1*xvector(t,1) + aalpha1*bbeta*xvector(t,2)+...
            + yerrors(t);
    end
    
    save dataTAU2.mat
    
else
    load('dataTAU2.mat')
end

%=========================================================================
%                 PLOT TIME SERIES
%=========================================================================

if plotdata==1
    
    figure;
    subplot(2,1,1)
    plot(1:T,xvector)
    xlabel('Period')
    title('X data')
    legend('x1','x2')
    
    subplot(2,1,2)
    hold on
    plot(1:T,yvector_a0,'-r')
    plot(1:T,yvector_a1,'--')
    xlabel('Period')
    title('Y data')
    legend('\alpha = 0.001','\alpha = 0.9')
    hold off
    
    print -depsc2 dataplotTAU2.eps
    
end


%=========================================================================
%                 PLOT LIKELIHOOD FUNCTION
%=========================================================================

if plotllhd==1
    
    % Discretize Parameter Space
    aalphagrid = linspace(-3,3);
    bbetagrid  = linspace(-3,3);
    
    llvalue1 = zeros(size(aalphagrid,2),size(aalphagrid,2));
    llvalue2 = zeros(size(bbetagrid,2),size(bbetagrid,2));
    
    for i = 1:100
        for j = 1:100
            llvalue1(j,i) = log(gaussllhd(yvector_a0,xvector,...
                [aalphagrid(j) bbetagrid(i)*aalphagrid(j)]',1));
            llvalue2(j,i) = log(gaussllhd(yvector_a0,xvector,...
                [aalphagrid(j) bbetagrid(i)*aalphagrid(j)]',1));
        end
    end
    
    figure;
    subplot(2,1,1)
    surf(aalphagrid,bbetagrid,llvalue1)
    xlabel('\alpha')
    ylabel('\beta')
    zlabel('logLikelihood')
    title('logLikelihood \alpha = 0.001')
    
    subplot(2,1,2)
    surf(aalphagrid,bbetagrid,llvalue2)
    xlabel('\alpha')
    ylabel('\beta')
    zlabel('logLikelihood')
    title('logLikelihood \alpha = 0.9')
    
    print -depsc2 llhdplotTAU2.eps
    
end

%=========================================================================
%                           GIBBS PARAMETERS
%=========================================================================

nsim    = 10000;
nburn   = 0.4*nsim;

%=========================================================================
%                           PARAMETERIZATION 1
%=========================================================================

% Set value for yvector to parameterization 1
yvector = yvector_a1; %#ok<NASGU>
display(' ');
display('Parameterization 1');
display(' ');

if gibbsrun==1
    gibbs1
    gibbs2
else
end


% Posterior Draws under Prior 1
posterior1_a1 = [ttheta1dp1(nburn+1:nsim), ...
    ttheta2dp1(nburn+1:nsim)./ttheta1dp1(nburn+1:nsim)];

% Posterior Draws under Prior 2
posterior2_a1 = [ttheta1dp2(nburn+1:nsim), ...
    ttheta2dp2(nburn+1:nsim)./ttheta1dp2(nburn+1:nsim)];


%=========================================================================
%                           PARAMETERIZATION 0
%=========================================================================

% Set value for yvector to parameterization 0
yvector = yvector_a0;
display(' ');
display('Parameterization 0');
display(' ');

if gibbsrun==1
    gibbs1
    gibbs2
else
    display('')
end

% Posterior Draws under Prior 1
posterior1_a0 = [ttheta1dp1(nburn+1:nsim), ...
    ttheta2dp1(nburn+1:nsim)./ttheta1dp1(nburn+1:nsim)];

% Posterior Draws under Prior 2
posterior2_a0 = [ttheta1dp2(nburn+1:nsim), ...
    ttheta2dp2(nburn+1:nsim)./ttheta1dp2(nburn+1:nsim)];



%=========================================================================
%                       POSTERIOR PLOTS
%=========================================================================

if plotdens==1
    
    % Marginal Density Plots (Parameterization 0)
    
    figure;
    subplot(2,2,1)
    [dens x] = ksdensity(posterior1_a0(:,1));
    plot(x,dens)
    title('Prior 1,\alpha = 0.001')
    xlabel('\alpha')
    
    subplot(2,2,2)
    [dens x] = ksdensity(posterior2_a0(:,1));
    plot(x,dens)
    title('Prior 2,\alpha = 0.001')
    xlabel('\alpha')
    
    subplot(2,2,3)
    [dens x] = ksdensity(posterior1_a0(:,2));
    plot(x,dens)
    title('Prior 1,\alpha = 0.001')
    xlabel('\beta')
    
    subplot(2,2,4)
    [dens x] = ksdensity(posterior2_a0(:,2));
    plot(x,dens)
    title('Prior 2,\alpha = 0.001')
    xlabel('\beta')
    
    print -depsc2 marginala0plotTAU2.eps
    
    
    % Marginal Density Plots (Parameterization 1)
    
    figure;
    subplot(2,2,1)
    [dens x] = ksdensity(posterior1_a1(:,1));
    plot(x,dens)
    title('Prior 1,\alpha = 0.9')
    xlabel('\alpha')
    
    subplot(2,2,2)
    [dens x] = ksdensity(posterior2_a1(:,1));
    plot(x,dens)
    title('Prior 2,\alpha = 0.9')
    xlabel('\alpha')
    
    subplot(2,2,3)
    [dens x] = ksdensity(posterior1_a1(:,2));
    plot(x,dens)
    title('Prior 1,\alpha = 0.9')
    xlabel('\beta')
    
    subplot(2,2,4)
    [dens x] = ksdensity(posterior2_a1(:,2));
    plot(x,dens)
    title('Prior 2,\alpha = 0.9')
    xlabel('\beta')
    
    print -depsc2 marginala1plotTAU2.eps
    
    
    % Joint Density Plot
    
    figure;
    plotmatrix(posterior2_a1)
    title('Joint Posterior, Prior 2, \alpha = 0.9')
    
    print -depsc2 posterior2_a1TAU2.eps
    
    figure;
    plotmatrix(posterior1_a1)
    title('Joint Posterior, Prior 1, \alpha = 0.9')
    print -depsc2 posterior1_a1TAU2.eps
    
    figure;
    plotmatrix(posterior2_a0)
    title('Joint Posterior, Prior 2, \alpha = 0.001')
    print -depsc2 posterior2_a0TAU2.eps
    
    figure;
    plotmatrix(posterior1_a0)
    title('Joint Posterior, Prior 1, \alpha = 0.001')
    print -depsc2 posterior1_a0TAU2.eps
    
end

MCMC