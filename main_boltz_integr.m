%% MECH 578 Project Part 2
% Numerically integrate Boltzmann equation for 1D gas as per 
% Boltzmann equations for a binary one-dimensional ideal gas
% Vassili Korotkine Nov. 19 2018
close all;clear all;clc

%% Initialize
%Variables

N=10; %num particles 
L=3; %length of wire 
ma=1; %mass a
mb=4; % mass b

delta_t=0.01; %time step size 
delta_v=0.2; % velocity step size
v_bound=10; %define velocity bound (symmetric)
Nv=2*v_bound/delta_v+1; %size of velocity vec
v=linspace(-v_bound, v_bound, Nv);

ideal_pA=arrayfun(@(x) 1/sqrt(pi)*exp(-x^2), v);
ideal_pB=arrayfun(@(x) 2/sqrt(pi)*exp(-4*x^2), v);

pA=zeros(1,Nv);
pB=zeros(1, Nv);
%Initial condition
pA(floor(Nv/4))=0.5/delta_v;
pB(floor(Nv/4))=0.5/delta_v;
pA(floor(Nv*3/4))=0.5/delta_v;
pB(floor(Nv*3/4))=0.5/delta_v;

pA=rand(1,Nv);
pB=rand(1,Nv);

for i = 1:Nv
    if abs(v(i))>1
        pA(i)=0;
        pB(i)=0;
    end
    
    
    %initial conditions corresponding to the analytical paper
    
    %B (smaller particle) initially at rest
    if v(i)==0
        pB(i)=1;
    else
        pB(i)=0;
    end
    
    %A (bigger particle) two spikes, one on right, one on left, at v0. see
    %above eq15 in analytical paper (Mohazabbi, Maxwellian relaxation)
    if abs(abs(v(i))-1)<0.1
        pA(i)=1;
    else
        pA(i)=0;
    end
    
end

pA=pA/(delta_v*sum(pA));
pB=pB/(delta_v*sum(pB));


n_steps=140;

A.avg=zeros(1,n_steps);
B.avg=zeros(1,n_steps);
A.std_dev=zeros(1,n_steps);
B.std_dev=zeros(1,n_steps);


for t_idx=1:n_steps
    %% Plot
    A.avg(t_idx)=sum(v.*pA*delta_v);
    B.avg(t_idx)=sum(v.*pB*delta_v);
    
    A.std_dev(t_idx)=sum((v-A.avg(t_idx)).^2 ...
                        .*pA*delta_v);
    B.std_dev(t_idx)=sum((v-B.avg(t_idx)).^2 ...
                        .*pB*delta_v);                    
    if rem(t_idx, 10)==0
    figure(1)
    bar(v,pA);hold on
    plot(v, ideal_pA); hold off
    xlabel('Velocity'); 
    ylabel('Probability dist value');
    legend('Current', 'Analytical')
    saveas(gcf,'pAiter','epsc')
    %title(sprintf('Iter %d, particle A Velocity Probability Distribution', t_idx))
    
    figure(2)
    bar(v, pB); hold on
    plot(v, ideal_pB); hold off;
    xlabel('v'); 
    ylabel('pB'); 
    legend('Current', 'Analytical')
    saveas(gcf,'pBiter','epsc')
    %title(sprintf('Iter %d, particle B Velocity Probability Distribution', t_idx))
    
    figure(3)
    plot(delta_t*(1:t_idx), A.std_dev(1:t_idx))
    xlabel('Time'); 
    ylabel('pA std dev');
    ylim([0, max([A.std_dev, 1])]);
    xlim([0,delta_t*n_steps]);
    saveas(gcf,'pASTD', 'epsc')
    %title(sprintf('Standard deviation of probability distribution, particle A', t_idx))
    
    figure(4)
    plot(delta_t*(1:t_idx), B.std_dev(1:t_idx))
    xlabel('time'); 
    ylabel('pB std dev');   
    ylim([0, max([B.std_dev,1])]);
    xlim([0,delta_t*n_steps]);
    saveas(gcf,'pBSTD','epsc')
    disp(t_idx)
    end 
%     subplot(3,2,5)
%     plot(delta_t*(1:t_idx), A.avg(1:t_idx))
%     xlabel('time'); 
%     ylabel('pA avg');
%     xlim([0,delta_t*n_steps]);
%     
%     subplot(3,2,6)
%     plot(delta_t*(1:t_idx), B.avg(1:t_idx))
%     xlabel('time'); 
%     ylabel('pB avg');
%     xlim([0,delta_t*n_steps]);
    
    delta_pA=zeros(1,Nv);
    for i=1:Nv
        %disp(i)
        %summing loop, compute RHS integral for single value of v_i. 
        temp=0;
        for j=1:(Nv)
            w_a=((ma-mb)*v(i)+2*mb*v(j))/(ma+mb);  %majorly illegal, to be checked
            w_b=((mb-ma)*v(j)+2*ma*v(i))/(ma+mb);
            
            if abs(w_a)>v_bound
                w_a=w_a/abs(w_a)*v_bound;
            end
            
            if abs(w_b)>v_bound
                w_b=w_b/abs(w_b)*v_bound;
            end
            
            
            
            temp=temp+ ...
                abs(v(i)-v(j))* ...
                (interp1(v',pA',w_a,'linear')*interp1(v',pB',w_b,'linear') - ...
                interp1(v',pA',v(i),'linear')*interp1(v',pB',v(j),'linear'));
            if isnan(delta_pA(i))
                print('abtin has shitty jokes')
            end
            
        end 
        delta_pA(i)=temp;
    end
    delta_pA=delta_pA*delta_v*delta_t*N/L;
    
    delta_pB=zeros(1,Nv);
    for i=1:Nv
        temp=0;
        for j=1:(Nv)
            w_a=((ma-mb)*v(j)+2*mb*v(i))/(ma+mb);
            w_b=((mb-ma)*v(i)+2*ma*v(j))/(ma+mb);
            
            
            % due to round off error 
            if abs(w_a)>v_bound
                w_a=w_a/abs(w_a)*v_bound;
            end
            
            if abs(w_b)>v_bound
                w_b=w_b/abs(w_b)*v_bound;
            end
            
            temp=temp+ ...
                abs(v(i)-v(j))* ...
                (interp1(v',pA',w_a, 'linear')*interp1(v',pB',w_b, 'linear') - ...
                interp1(v',pA',v(j), 'linear')*interp1(v',pB',v(i), 'linear'));
        end 
        delta_pB(i)=temp;
    end
    delta_pB=delta_pB*delta_v*delta_t*N/L;
    
    pA=pA+delta_pA;
    pB=pB+delta_pB;
    
    
    

    
    
%     subplot(2,3,2)
%     bar(v, delta_pA)
%     xlabel('v'); 
%     ylabel('delta_pA'); 
%     
%     subplot(2,2,4)
%     bar(v, delta_pB)
%     xlabel('v'); 
%     ylabel('delta_pB'); 
    
       
    
    pause(0.001);
        
end
    
    
    


















