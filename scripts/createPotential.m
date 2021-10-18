function model=createPotential(type,varargin)
% creates a energy potential for Langevin dynamic simulations
%
% custom energy profile
% model=createPotential('custom',xF,xU,xTS,GTS_F,dG,depth)
%
% harmonic energy profile
% model=createPotential('harmonic',x0,kappa)
% 
% quartic energy profile
% model=createPotential('custom',xF,xU,GTS_F)
%
% sinusoidal double well energy profile
% model=createPotential('custom',xF,xU,GTS_F,'d')
%
% sinusoidal repetitive energy profile
% model=createPotential('custom',xF,xU,GTS_F,'r')
%
% Andreas Hartmann | 10.15.2021
    
    if strcmp(type,'custom')
        
        xF=varargin{1};
        xU=varargin{2};
        xTS=varargin{3};
        GTS_F=varargin{4};
        dG=varargin{5};
        depth=varargin{6};        
        
        xx=(xF-5:0.001:xU+5);
    
        % skewed normal distribution
        f=@(x,mu,sig,alpha) 2*normpdf(x,mu,sig).*normcdf(alpha*(x-mu));  
        
        % first guess for fitting
        opt_x=[0.5 0.5 xF xU 0.25*(xU-xF) 0.25*(xU-xF)];
        
        options = optimset('Display','iter','TolX',1e-12,'TolFun',1e-12);

        chi2=1000; % default
        numIter=0;

        % fitting loop with updated boundaries and starting point
        while (chi2>0.001)&&(numIter<10)

            % updated starting point
            x0=[opt_x(1) opt_x(2) opt_x(3) opt_x(4) 0.25*(opt_x(4)-opt_x(3)) 0.25*(opt_x(4)-opt_x(3))];
            % lower boundary
            lb=[0 0 xF-5 xTS 0 0];
            % upper boundary
            ub=[Inf Inf xTS xU+5 0.5*(opt_x(4)-opt_x(3)) 0.5*(opt_x(4)-opt_x(3))];

            % fitting of the potential features
            opt_x=fminsearchbnd(@myPotential,x0,lb,ub,options);
            
            % chi-square update
            [chi2,~]=myPotential(opt_x);

            numIter=numIter+1;
        end
        
        % optimized custom potential
        [~,optG]=myPotential(opt_x);
        
        offset=min(optG);

        G=@(x,xdata) -log(x(1).*(xdata<x(4)).*2.*normpdf(xdata,x(3),x(5)).*normcdf(x(7).*(xdata-x(3)))+x(2).*(xdata>x(3)).*2.*normpdf(xdata,x(4),x(6)).*normcdf(-x(7).*(xdata-x(4))))-x(8);
        grad_G=@(x,xdata) (G(x,xdata+1e-8)-G(x,xdata-1e-8))./(2e-8);
        
        model{1}='custom';
        model{2}=G;
        model{3}=grad_G;
        model{4}=[opt_x(:);depth;offset];
        model{5}=varargin;
    
    elseif strcmp(type,'harmonic')
        
        x0=varargin{1}; % center position
        kappa=varargin{2}; % curvature   
        
        G=@(x,xdata) (x(2)/2).*(xdata-x(1)).^2;
        grad_G=@(x,xdata) x(2).*(xdata-x(1));
        
        model{1}='harmonic';
        model{2}=G;
        model{3}=grad_G;
        model{4}=[x0;kappa];
        model{5}=varargin;
        
    elseif strcmp(type,'quartic')
        
        xF=varargin{1};
        xU=varargin{2};
        GTS_F=varargin{3};
        
        G=@(x,xdata) x(3).*(4.*((xdata-0.5.*(x(2)+x(1)))./(x(2)-x(1))).^2-1).^2;
        grad_G=@(x,xdata) (16*x(3)/(x(2)-x(1))).*(4.*((xdata-0.5*(x(2)+x(1)))./(x(2)-x(1))).^3-(xdata-0.5*(x(2)+x(1)))./(x(2)-x(1)));
        
        model{1}='quartic';
        model{2}=G;
        model{3}=grad_G;
        model{4}=[xF;xU;GTS_F];
        model{5}=varargin;
        
    elseif strcmp(type,'sinusoidal')
        
        xF=varargin{1};
        xU=varargin{2};
        GTS_F=varargin{3};
        border=varargin{4};      
        
        if strcmp(border,'d')

            G=@(x,xdata) (xdata<=x(1)).*((2*(pi^2)*x(3))/((x(2)-x(1))^2)/2).*(xdata-x(1)).^2+(xdata>x(1)&xdata<x(2)).*(0.5*x(3)).*(1-cos(2*pi*(xdata-x(1))./(x(2)-x(1))))+(xdata>=x(2)).*((2*(pi^2)*x(3))/((x(2)-x(1))^2)/2).*(xdata-x(2)).^2;
            grad_G=@(x,xdata) (xdata<=x(1)).*((2*(pi^2)*x(3))/((x(2)-x(1))^2)).*(xdata-x(1))+(xdata>x(1)&xdata<x(2)).*(pi*x(3)/(x(2)-x(1))).*sin(2*pi.*(xdata-x(1))./(x(2)-x(1)))+(xdata>=x(2)).*((2*(pi^2)*x(3))/((x(2)-x(1))^2)).*(xdata-x(2));
            
        elseif strcmp(border,'r')
            
            G=@(x,xdata) (0.5*x(3)).*(1-cos(2*pi*(xdata-x(1))./(x(2)-x(1))));
            grad_G=@(x,xdata) (pi*x(3)/(x(2)-x(1))).*sin(2*pi.*(xdata-x(1))./(x(2)-x(1)));
        end
        
        model{1}='sinusoidal';
        model{2}=G;
        model{3}=grad_G;        
        model{4}=[xF;xU;GTS_F];        
        model{5}=varargin;
    end
        
    function [chi2,G]=myPotential(x)

        A1_=x(1);
        A2_=x(2);
        mu1_=x(3);
        mu2_=x(4);
        sig1_=x(5);
        sig2_=x(6);
        
        f1=zeros(length(xx),1);
        f2=zeros(length(xx),1);
        f1(xx<mu2_)=f(xx(xx<mu2_),mu1_,sig1_,depth);
        f2(xx>mu1_)=f(xx(xx>mu1_),mu2_,sig2_,-depth);

        G=-log(A1_.*f1+A2_.*f2);

        xxL=xx(xx<0.5*(mu1_+mu2_));
        GL=G(xx<0.5*(mu1_+mu2_));

        [GF_,idxF]=min(GL);
        xF_=xxL(idxF);

        xxR=xx(xx>0.5*(mu1_+mu2_));
        GR=G(xx>0.5*(mu1_+mu2_));

        [GU_,idxU]=min(GR);
        xU_=xxR(idxU);

        xxB=xx(xF_<xx&xx<xU_);
        GB=G(xF_<xx&xx<xU_);

        [GTS_,idxM]=max(GB);
        xTS_=xxB(idxM);
        
        GTS_F_=GTS_-GF_;
        
        dG_=GU_-GF_;
        
        chi2=(xF-xF_)^2+(xU-xU_)^2+(xTS-xTS_)^2+(dG-dG_)^2+(GTS_F-GTS_F_)^2;
    end
end