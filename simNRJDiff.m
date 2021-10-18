%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LANGEVIN DYNAMIC 
% BROWNIAN MOTION IN AN ENERGY LANDSCAPE
%
% xF (nm) -> folded state reaction coordinate
% xU (nm) -> unfolded state reaction coordinate
% xTS (nm) -> transition state reaction coordinate
% x0 (nm) -> center coordinate in the case of a harmonic potentail
% dGTS_F (kBT) -> energy barrier of the folded state
% dG=2 (kBT) -> change in free Gibbs energy between folded and unfolded state 
% depth -> steepness of the energy barrier borders
% qL, qR (nm) -> left and right reaction coordinate for transition path time measurements
%                qL < qR!!!
%
% Andreas Hartmann | 17.10.2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
addpath('scripts');

% selection of the energy landscape model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model=createPotential('custom',2,8,3.5,3,2,3); % (xF,xU,xTS,dGTS_F,dG,depth)
% model=createPotential('harmonic',5,10); % (x0,kappa)
model=createPotential('quartic',2,8,3); % (xF,xU,dGTS_F)
% model=createPotential('sinusoidal',2,8,3,'d'); % (xF,xU,dGTS_F,type) type=='d' -> double well or type=='r' -> repetitive

% left and right reaction coordinate for transition path time measurements
qL=3.5; % (nm) qL < qR!!!
qR=6.5; % (nm)

% simulation paramters
D=1; % (nm^2/탎) diffusion constant in the energy landscape
dt=0.01; % (탎) time step size
numDP=100000; % data points per trace
numT=1000; % number of traces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% loading and plotting of the energy landscape model
type=model{1};
G=model{2};
grad_G=model{3};
params=model{4};
input=model{5};

if strcmp(type,'custom')

    % reaction coordinate
    xx=(input{1}-5:0.001:input{2}+5);
    edges=(input{1}-5:0.1:input{2}+5);

    % energy landscape
    fig1=figure();
    plot(xx,G(params,xx),'LineWidth',2);
    hold on;
    plot(qL,G(params,qL),'or');
    plot(qR,G(params,qR),'or');
    xlim([input{1}-2 input{2}+2]);
    ylim([0 1.5*input{4}]);
    xlabel('\itx\rm (nm)');
    ylabel('\itG\rm(\itx\rm) (\itk\rm_B\itT)');
    legend(['\it' type '\rm - energy profile'],'\itq\rm_L & \itq\rm_R');

elseif strcmp(type,'harmonic')
    
    % reaction coordinate
    xx=(input{1}-5:0.001:input{1}+5);
    edges=(input{1}-5:0.1:input{1}+5);

    % energy landscape
    fig1=figure();
    plot(xx,G(params,xx),'LineWidth',2);
    hold on;
    plot(qL,G(params,qL),'or');
    plot(qR,G(params,qR),'or');
    xlim([input{1}-2 input{1}+2]);
    ylim([0 10]);
    xlabel('\itx\rm (nm)');
    ylabel('\itG\rm(\itx\rm) (\itk\rm_B\itT)');
    legend(['\it' type '\rm - energy profile'],'\itq\rm_L & \itq\rm_R');    
    
elseif strcmp(type,'quartic')||strcmp(type,'sinusoidal')
    
    % reaction coordinate
    xx=(input{1}-5:0.001:input{2}+5);
    edges=(input{1}-5:0.1:input{2}+5);

    % energy landscape
    fig1=figure();
    plot(xx,G(params,xx),'LineWidth',2);
    hold on;
    plot(qL,G(params,qL),'or');
    plot(qR,G(params,qR),'or');
    xlim([input{1}-2 input{2}+2]);
    ylim([0 1.5*input{3}]);
    xlabel('\itx\rm (nm)');
    ylabel('\itG\rm(\itx\rm) (\itk\rm_B\itT)');
    legend(['\it' type '\rm - energy profile'],'\itq\rm_L & \itq\rm_R');    
end

xl=xlim;

%% Simulation
histX=zeros(length(edges),1);

arrTPT=zeros(10000,1);
ptrTPT=0; % pointer for TPT field

% starting position
if strcmp(type,'harmonic')
    
    x_arr=input{1}; % center of the potential
else
    x_arr=input{2}; % unfolded state
end

for iterOS=1:numT
    
    x=x_arr(end);
    
    boolU=x>qR;
    boolTS=(x>=qL)&&(x<=qR);
        
    x_arr=zeros(numDP,1);
    
    arrRND=randn(numDP,1);
    
    lenTS=0;
    
    for iter=1:numDP
        
        % Langevin step update
        dx=dt*(-D*grad_G(params,x)+sqrt(2*D/dt)*arrRND(iter));
        x=x+dx;
        
        x_arr(iter)=x;
        
        % folding transition path time analysis
        if (x>qL)&&(x<qR)&&((boolU==1)||(boolTS==1))
            
            boolTS=1;
            boolU=0;
            
            lenTS=lenTS+1;
        end
        
        if (x<qL)&&(boolTS==1)
            
            ptrTPT=ptrTPT+1;
            arrTPT(ptrTPT)=lenTS;
            
            lenTS=0;
            boolTS=0;
            boolU=0;
        end
        
        if x>qR
            
            boolU=1;
            boolTS=0;
            lenTS=0;
        end
        
        if mod(iter,10000)==0
        
            disp([num2str(iterOS) ' - ' num2str(iter) '/' num2str(numDP)]);
        end
    end

    hX=histc(x_arr,edges);
    histX=histX+hX(:);
end

%% Plotting of the simulation results
movegui(fig1,'north');

fig2=figure();
subplot(2,2,1);
b1=bar(edges,histX./trapz(edges,histX),'histc');
set(b1,'EdgeColor','non','FaceColor',[0.5 0.5 0.5]);
hold on;
plot(edges,exp(-G(params,edges))./trapz(edges,exp(-G(params,edges))),'-k','LineWidth',2);
xlim(xl);
ylim([0 1.1.*max(histX./trapz(edges,histX))]);
xlabel('\itx\rm (nm)');
ylabel('\itp\rm (\itx\rm)');

subplot(2,2,2);
histogram(arrTPT(1:ptrTPT).*dt,50);
yL=get(gca,'YLim');
hold on;
plot([mean(arrTPT(1:ptrTPT).*dt) mean(arrTPT(1:ptrTPT).*dt)],[0 yL(2)],'--r');
xlabel('\itt\rm_{TP} (탎)');
ylabel('Counts');
legend('data',['<\itt\rm_{TP}>=' num2str(mean(arrTPT(1:ptrTPT).*dt)) '탎']);

subplot(2,2,[3 4]);
plot((1:1:numDP).*dt,x_arr);
ylabel('\itx\rm (nm)');
xlabel('\itt\rm (ns)');

movegui(fig2,'south');

