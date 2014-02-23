%==========================================================================
%                   VECTOR ERROR CORRECTION MODEL
%
%
%                      See README_VEC for details
%
%
% Author: Luigi Bocola     lbocola@sas.upenn.edu
% Date  : 06/20/2010
%
% Modified by Kory Kantenga
% Date  : 02/19/2014
%
%==========================================================================


%=========================================================================
%                              HOUSEKEEPING
%=========================================================================

tic
close all
clear all
clc

%Turn on/off plots & sampler
plotdata  = 0;
plottrend = 0;
plotavg   = 0;
plotkernel= 0;
gibbsrun  = 0;

%=========================================================================
%                        LOAD AND PREPARE THE DATA
%=========================================================================

load data

[t,n] = size(gdp);

p = 3;

T = t - p;

Dy = [investment(4:t)-investment(3:t-1), gdp(4:t) - gdp(3:t-1)];

X  = [investment(3:t-1),gdp(3:t-1)];

W = [ones(T,1), investment(2:t-2)-investment(1:t-3),gdp(2:t-2)-gdp(1:t-3),...
    investment(3:t-1)-investment(2:t-2), gdp(3:t-1)-gdp(2:t-2)];


nsim  = 10000;       % Number of draws from Posterior Density
nburn = 0.4*nsim;    % Number of draws to discart


%=========================================================================
%                 PLOT TIME SERIES
%=========================================================================

if plotdata == 1
    figure;
    subplot(2,1,1)
    time=1964:0.25:2006.75;
    [AX,H1,H2] = plotyy(time',gdp,time',investment,'plot');
    xlabel('Year')
    set(get(AX(1),'Ylabel'),'String','log Nominal GDP')
    set(get(AX(2),'Ylabel'),'String','log Investment')
    legend('GDP','Investment','Location','Best')
    set(H1,'LineStyle','--')
    set(H2,'LineStyle','-')
    
    subplot(2,1,2)
    plot(time',investment-gdp)
    xlabel('Year')
    ylabel('log(INVESTMENT/GDP)')
    
    print -dpsc2 dataplot.eps
end

%=========================================================================
%                 STARTING VALUES FOR GIBBS SAMPLER
%=========================================================================

lambda = 0.1;            % Parameter For the Prior
Sigmap = zeros(nsim,2,2);
Gammap = zeros(nsim,5,2);
alphap = ones(2,nsim);
alphap(:,1) = [-.14,0.02];
bp     = -ones(nsim,1);
p      = 2;
counter= 0;
%=========================================================================
%             DRAWS FROM POSTERIOR DENSITY: GIBBS SAMPLING
%=========================================================================
if gibbsrun==1
    disp('                                                                  ');
    disp('        BAYESIAN ESTIMATION OF VEC: GIBBS SAMPLING...             ');
    disp('                                                                  ');
    
    
    for j=2:nsim
        
        % Draws from the density Sigma | Y, Gamma(j-1), alpha (j-1), b (j-1)
        resid = (Dy - X*[1;bp(j-1)]*alphap(:,j-1)' - W*squeeze(Gammap(j-1,:,:)));
        Sigma = resid'*resid;
        
        Sigma = iwishrnd(Sigma,T);
        
        Sigmap(j,:,:) = Sigma;
        
        % Draws from the density vec(Gamma) |Y, Sigma(j), alpha(j-1), b(j-1)
        
        Dytilde  = Dy - X*[1;bp(j-1)]*alphap(:,j-1)';
        
        Gammahat = (W'*W)\(W'*Dytilde);
        
        Gammanew = mvnrnd(reshape(Gammahat,2*(2*p+1),1)',kron(squeeze(Sigmap(j,...
            :,:)),inv(W'*W)))';
        
        Gammap(j,:,:) = reshape(Gammanew,5,2);
        
        % Draws from the density alpha |Y, Sigma(j), Gamma(j), b(j-1)
        
        Zeta   = Dy - W*squeeze(Gammap(j,:,:));
        
        Xtilde = X*[1;bp(j-1)];
        
        alphahat = inv(Xtilde'*Xtilde)'*Xtilde'*Zeta;
        
        alpha = mvnrnd(alphahat,kron(squeeze(Sigmap(j,:,:)),...
            inv(Xtilde'*Xtilde)))';
        
        alphap(:,j) = alpha;
        
        % Draws from the density beta |Y, Sigma(j), Gamma(j), alpha(j)
        
        Zeta = Dy - W*squeeze(Gammap(j,:,:)) - X(:,1)*alphap(:,j)';
        
        C = [(alphap(:,j)/((alphap(:,j)'*alphap(:,j)))), [1,-alphap(1,j)/alphap(2,j)]'];
        
        Zetatilde = Zeta*C;
        
        Sigmatilde  = C'*squeeze(Sigmap(j,:,:))*C;
        
        Zeta12      = Zetatilde(:,1) - Sigmatilde(1,2)*(Sigmatilde(2,2)^(-1))*...
            Zetatilde(:,2);
        
        betaols     = (X(:,2)'*X(:,2))\(X(:,2)'*Zeta12);
        
        invvarols   = inv(((Sigmatilde(1,1) - (Sigmatilde(1,2)*Sigmatilde(2,1)...
            /Sigmatilde(2,2))))/((X(:,2)'*X(:,2))))  ;
        
        varbetapost = inv(invvarols + lambda^(-1));
        
        betahatpost = varbetapost*(invvarols*betaols + (lambda^(-1))*(-1.00)); %#ok<MINV>
        
        bp(j,1)     = betahatpost + ((varbetapost)^(1/2))*randn(1,1);
        
        counter     = counter + 1;
        
        if counter==3000
            disp(['         DRAW NUMBER:   ', num2str(j)]);
            disp('                                                                  ');
            disp(['     REMAINING DRAWS:   ', num2str(nsim-j)]);
            disp('                                                                  ');
            
            counter = 0;
            
        end
        
        
    end
    
    
    %=========================================================================
    %           FIGURE 1: RECURSIVE AVERAGES, SELECTED PARAMETERS
    %=========================================================================
    
    if plotavg==1
        
        Theta       = [bp,alphap',Sigmap(:,1,1),Sigmap(:,2,2)];
        
        pnames      = char('\beta','\alpha_{1}', '\alpha_{2}','\sigma^{2}_{1,1}',...
            '\sigma^{2}_{2,2}');
        
        
        
        figure('Position',[20,20,900,600],'Name',...
            'Recursive Averages','Color','w')
        
        rmean       = zeros(0.2*nsim,5);
        
        for i=1:0.2*nsim
            rmean(i,:) = mean(Theta(1:i,:),1);
        end
        
        for i=1:5
            
            subplot(3,2,i), plot(rmean(:,i),'LineStyle','-','Color','b',...
                'LineWidth',2.5),
            
            title(pnames(i,:),'FontSize',12,'FontWeight','bold');
        end
        
    end
    
    %=========================================================================
    %              STOCHASTIC TREND AND STATIONARY COMPONENT
    %=========================================================================
    
    % Defining some objects
    
    Strend = zeros(nsim-nburn +1,2,T);
    Sstat  = zeros(nsim-nburn +1,2,T);
    
    
    % Estimating the Stochastic Trend and the Stationary Component
    disp('                                                                  ');
    disp('   ESTIMATING STOCHASTIC TREND AND STATIONARY COMPONENT...        ');
    disp('                                                                  ');
    
    for j=nburn:nsim
        % Stochastic Trend
        resid        = (Dy - X*[1;bp(j)]*alphap(:,j)' - W*squeeze(Gammap(j,:,:)));
        Bnull        = [1;(-1/bp(j,1))];
        Anull        = [1,-alphap(1,j)/alphap(2,j)]';
        PHI0         = squeeze(Gammap(j,1,:));
        PHI1         = squeeze(Gammap(j,2:3,:))';
        PHI2         = squeeze(Gammap(j,4:5,:))';
        C            = Bnull*((Anull'*(eye(2)-PHI1-PHI2)*Bnull)^(-1))*Anull';
        Strend(j-nburn+1,:,:) = C*cumsum(resid)' + C*PHI0*cumsum(ones(1,T));
        
        % (De-meaned) Stationary Component
        
        %y(:,1) = resid(1,:)';
        %y(:,2) = resid(2,:)'+PHI1*resid(1,:)';
        
        %for i=3:T
        %   y(:,i) = resid(i,:)' + PHI1*y(:,i-1) + PHI2*y(:,i-2);
        %end
        
        % Sstat(j-nburn+1,:,:) = y;
        
        Sstat(j-nburn+1,:,:) = X'- squeeze(Strend(j-nburn+1,:,:));
        Sstat(j-nburn+1,:,:) = squeeze(Sstat(j-nburn+1,:,:)) - ...
            mean(Sstat(j-nburn+1,:,:),3)'*ones(1,T);
        
    end
    
    % Constructing 90% Credible Sets (HPD) for stochastic trend and stationary
    % component
    
    [nsim,n,T]= size(Strend);
    
    trendmed  = squeeze(median(Strend,1));
    
    statmed   = squeeze(median(Sstat,1));
    
    trend_sort= sort(Strend,1);
    
    stat_sort = sort(Sstat,1);
    
    confidence = 0.9;
    
    M          = floor((1-confidence)*nsim);
    
    h1         = zeros(M,T,n);
    h2         = zeros(M,T,n);
    trendlower = zeros(n,T);
    trendupper = zeros(n,T);
    statlower  = zeros(n,T);
    statupper  = zeros(n,T);
    
    for l=1:2
        for i=1:T
            for j=1:M
                
                h1(j,i,l) = abs(trend_sort(j,l,i) - trend_sort(j + ...
                    round(confidence*nsim),l,i));
                h2(j,i,l) = abs(stat_sort(j,l,i) - stat_sort(j + ...
                    round(confidence*nsim),l,i));
                
            end
            
            [~,k1] = min(h1(:,i,l));
            [~,k2] = min(h2(:,i,l));
            
            trendlower(l,i) = trend_sort(k1,l,i);
            trendupper(l,i) = trend_sort(k1 + round(confidence*nsim),l,i);
            statlower(l,i)  = stat_sort(k2,l,i);
            statupper(l,i)  = stat_sort(k2 + round(confidence*nsim),l,i);
            
        end
    end
    
    
    %=========================================================================
    %           FIGURE 2: STOCHASTIC TREND AND STATIONARY COMPONENT
    %=========================================================================
    
    if plottrend==1
        
        time = 1964.75:0.25:2006.75;
        
        % NBER recession date
        rec1 = 1969.75:0.25:1970.75;
        rec2 = 1973.75:0.25:1975.00;
        rec3 = 1980.00:0.25:1980.50;
        rec4 = 1981.50:0.25:1982.75;
        rec5 = 1990.50:0.25:1991.00;
        rec6 = 2001.00:0.25:2001.75;
        
        
        
        pnames      = char('Stochastic Trend (Inv)','Stochastic Trend (Output)',...
            'Stationary Component (Inv.)','Stationary Component (Output)');
        
        
        
        figure('Position',[20,20,900,600],'Name',...
            'Stochastic Trend and Stationary Components','Color','w')
        
        
        for i=1:2
            ymax = (1.1)*max(trendmed(i,:));
            ymin = 0;
            subplot(2,2,i),
            area(rec1,ymax*ones(length(rec1),1),ymin,'FaceColor', [.5 .5 .5]), hold on
            area(rec2,ymax*ones(length(rec2),1),ymin,'FaceColor', [.5 .5 .5])
            area(rec3,ymax*ones(length(rec3),1),ymin,'FaceColor', [.5 .5 .5])
            area(rec4,ymax*ones(length(rec4),1),ymin,'FaceColor', [.5 .5 .5])
            area(rec5,ymax*ones(length(rec5),1),ymin,'FaceColor', [.5 .5 .5])
            area(rec6,ymax*ones(length(rec6),1),ymin,'FaceColor', [.5 .5 .5])
            plot(time,trendmed(i,:),'LineStyle','-','Color','b',...
                'LineWidth',2), hold on
            plot(time,trendlower(i,:),'LineStyle',':','Color','r',...
                'LineWidth',1.5), hold on
            plot(time,trendupper(i,:),'LineStyle',':','Color','r',...
                'LineWidth',1.5), hold on
            xlabel('Year')
            axis([time(1) 2009 ymin ymax])
            title(pnames(i,:),'FontSize',12,'FontWeight','bold');
        end
        
        
        for i=1:2
            ymax =4*max(statmed(i,:));
            ymin =4*min(statmed(i,:));
            subplot(2,2,2+i),
            area(rec1,ymax*ones(length(rec1),1),ymin,'FaceColor', [.5 .5 .5]), hold on
            area(rec2,ymax*ones(length(rec2),1),ymin,'FaceColor', [.5 .5 .5])
            area(rec3,ymax*ones(length(rec3),1),ymin,'FaceColor', [.5 .5 .5])
            area(rec4,ymax*ones(length(rec4),1),ymin,'FaceColor', [.5 .5 .5])
            area(rec5,ymax*ones(length(rec5),1),ymin,'FaceColor', [.5 .5 .5])
            area(rec6,ymax*ones(length(rec6),1),ymin,'FaceColor', [.5 .5 .5])
            plot(time,statmed(i,:),'LineStyle','-','Color','b',...
                'LineWidth',2), hold on
            plot(time,statlower(i,:),'LineStyle',':','Color','r',...
                'LineWidth',1.5), hold on
            plot(time,statupper(i,:),'LineStyle',':','Color','r',...
                'LineWidth',1.5), hold on
            xlabel('Year')
            axis([time(1) 2009 ymin ymax])
            title(pnames(2+i,:),'FontSize',12,'FontWeight','bold');
        end
        
        disp(['         ELAPSED TIME:   ', num2str(toc)]);
        
    end
    
    %=========================================================================
    %           FIGURE 3: KERNEL DENSITY OF B
    %=========================================================================
    
    if plotkernel==1
        
        [density, x] = ksdensity(bp);
        
        figure;
        hold on
        plot(x,density)
        xlabel('b')
        h = line([-1 -1], [0 max(density)+1]);
        set(h,'Color','r','LineStyle','--')
        ylim([0 max(density)+1])
        hold off
        print -dpsc2 Bkernel.eps
        
    end
    
end

elapsedtime=toc;
