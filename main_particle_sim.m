%% MECH 578 Project Part 2
% Direct particle simulation, 1D particles on a wire
close all;clear all;clc

global radius;
global m_a;
global m_b;
global L;

N=5;
L=10;

m_a=4;
m_b=1;

radius=L/20; %how close particles have to be before we consider collision to have occured. 

%A IS THE LARGE PARTICLE
A_part=zeros(1,N);
[pos,vel, m]=vel_init(A_part);
 %values of particle mes w.r.t. particle idx, alternating arr 
%pos=sort(rand(1,N))*L;
right_col.rel_pos=zeros(1,N);
right_col.rel_vel=zeros(1,N);
right_col.tc=zeros(1,N);



n_t=40;
%COLLISIONS

t_arr=zeros(1,n_t); %times between collision
energy=zeros(1,n_t); %plot energy for confirmation 
for t_idx=1:n_t
    %update relative positions, velocities
    [min_tc,min_idx,left_flag, right_flag]=parsim_get_coll(pos, vel);
    [pos, vel] = parsim_solve_collision(pos,vel,min_tc,min_idx,left_flag,right_flag);
    

    
   energy(t_idx)=sum((m.*vel.^2))/2;
   
   
   
   % PLOTS
   figure(1)
   N_plots=6;
   subplot(N_plots,1,1)
   histogram(vel, 10)
   title('vel hist')
   subplot(N_plots,1,2)
   plot(1:t_idx, energy(1:t_idx))
   title('energy vs time')
   
   %figure(2)
   %N_plots=2;
   subplot(N_plots,1,3)
   plot(pos)
   title('pos')
   subplot(N_plots,1,4)
   plot(vel)
   title('vel')
   
   
   subplot(N_plots,1,5)
   %MAKES SIZE WITH PARTICLE RADIUS, BOILERPLATE
   title('particles')
    A_part=logical(A_part);
    s=radius; %particle size

    h=scatter(pos(A_part), zeros(1,length(pos(A_part)))); hold on
    ha=scatter(pos(~A_part), zeros(1,length(pos(~A_part)))); hold on
    hb=scatter(pos(min_idx), 0, 'filled'); hold off
    currentunits = get(gca,'Units');
    set(gca, 'Units', 'Points');
    axpos = get(gca,'Position');
    set(gca, 'Units', currentunits);
    markerWidth = s/diff(xlim)*axpos(3); % Calculate Marker width in points
    set(h, 'SizeData', markerWidth^2)
    set(ha, 'SizeData', markerWidth^2)
    set(hb, 'SizeData', markerWidth^2)
    xlim([0,L])
   %legend('A type', 'B type')
   %min_idx
   subplot(N_plots,1,6)
   plot(min_tc*vel)
   pause(0.5)
%    subplot(N_plots,1,2)
%    plot(1:t_idx, energy(1:t_idx))
%     % debug plots
%     figure(2)
%     N_plots=6;
%     subplot(N_plots,1,1)
%     plot(A_part); title('type')
%     subplot(N_plots,1,2)
%     plot(vel); title('vel')
%     subplot(N_plots,1,3)
%     plot(right_col.rel_vel); title('rel_vel,right')
%     subplot(N_plots,1,4)
%     plot(pos); title('pos')
%     subplot(N_plots,1,5)
%     plot(right_col.rel_pos); title('rel_pos, right')
%     subplot(N_plots,1,6)
%     plot(right_col.tc); title('tc, right')
% 
%     figure(3)
%     N_plots=4;
%     subplot(N_plots,1,1)
%     plot(right_col.rel_pos); title('rel_pos, right')    
%     subplot(N_plots,1,2)
%     plot(right_col.rel_vel); title('rel_vel,right')
%     subplot(N_plots,1,3)
%     plot(right_col.tc); title('tc, right')
%     subplot(N_plots,1,4)
%     plot(left_col.tc); title('tc, left')
x=0;
end






function [pos,vel, m]=vel_init(A_part)

global radius;
global m_a;
global m_b;
%INITIALIZATION
N=length(A_part);
vel_alter=0; %allows to flip velocities for the smaller B particles. 
m=zeros(1,N);


for i=1:N
    if i==1
        pos(i)=2*radius;
    else
        pos(i)=pos(i-1)+2*radius+rand;
    end
    
    if rem(i,2) ==0
        A_part(i)=1;
        m(i)=m_a;
    else
        A_part(i)=0;
        m(i)=m_b;
    end
    
    if A_part(i)
        vel(i)=0;
    else
        if vel_alter
            vel(i)=1;
            vel_alter=0;
        else
            vel(i)=-1;
            vel_alter=1;
        end
    end    
end

end
