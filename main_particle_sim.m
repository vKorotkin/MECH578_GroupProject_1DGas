%% MECH 578 Project Part 2
% Direct particle simulation, 1D particles on a wire
close all;clear all;clc

N=10;
L=10;

m_a=4;
m_b=1;
coll_tol=L/30; %how close particles have to be before we consider collision to have occured. 

%A IS THE LARGE PARTICLE
A_part=zeros(1,N);
pos=sort(rand(1,N))*L;
vel=zeros(1,N);
right_col.rel_pos=zeros(1,N);
right_col.rel_vel=zeros(1,N);
right_col.tc=zeros(1,N);
left_col.rel_pos=zeros(1,N);


%INITIALIZATION
vel_alter=0; %allows to flip velocities for the smaller B particles. 
for i=1:N
    if rem(i,2) ==0
        A_part(i)=1;
    else
        A_part(i)=0;
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


n_t=3;
%COLLISIONS

t_arr=zeros(1,n_t); %times between collision
for t_idx=1:n_t
    %update relative positions, velocities
    for i=1:N
        %temp vars, computation
        if i ~=N
            rel_p=pos(i)-pos(i+1);
            rel_v=vel(i)-vel(i+1);
            tc=-rel_p/rel_v; %minus is intentional, write this out with the particles. will see is correct. 
            %storage

        else %rightmost particle vs wall 
            rel_p=L-pos(N);
            rel_v=vel(i);
            tc=rel_p/rel_v; %no minus sign
        end
        
        if tc<0
            tc=5;
        end
        % storage
        right_col.rel_pos(i)=rel_p;
        right_col.rel_vel(i)=rel_v;
        right_col.tc(i)=tc;
        
        if i ~=1
            rel_p=pos(i)-pos(i-1);
            rel_v=vel(i)-vel(i-1);
            tc=-rel_p/rel_v; %same rationale as above.  write this out. 
            
        else
            rel_p=pos(i);
            rel_v=vel(i);
            tc=-rel_p/rel_v; %leftmost particle, again think about the minus sign.
        end
        
        if tc<0
            tc=5;
        end
        
        
        %storage
        left_col.rel_pos(i)=rel_p;
        left_col.rel_vel(i)=rel_v;
        left_col.tc(i)=tc;
    end
    
    
    [min_tc,min_idx]=min([right_col.tc]); %this wasnt well thought through, probably didnt need left_coll at all. 
    
    if left_col.tc(1)<min_tc
        min_tc=left_col.tc(1); 
        min_idx=1;
    end
    
    
    x=0;
    % debug plots
   figure(1)
    N_plots=6;
    subplot(N_plots,1,1)
    plot(A_part); title('type')
    subplot(N_plots,1,2)
    plot(vel); title('vel')
    subplot(N_plots,1,3)
    plot(right_col.rel_vel); title('rel_vel,right')
    subplot(N_plots,1,4)
    plot(pos); title('pos')
    subplot(N_plots,1,5)
    plot(right_col.rel_pos); title('rel_pos, right')
    subplot(N_plots,1,6)
    plot(right_col.tc); title('tc, right')

    figure(2)
    N_plots=4;
    subplot(N_plots,1,1)
    plot(right_col.rel_pos); title('rel_pos, right')    
    subplot(N_plots,1,2)
    plot(right_col.rel_vel); title('rel_vel,right')
    subplot(N_plots,1,3)
    plot(right_col.tc); title('tc, right')
    subplot(N_plots,1,4)
    plot(left_col.tc); title('tc, left')
x=0;
end



% for i=1:N
%     coll_table.right.rel_pos(i)=1
%     coll_table.right.rel_vel(i)=1
%     
%     
%     coll
%     
% end
