function [min_tc,min_idx,left_flag, right_flag]=parsim_get_coll(pos, vel)
global radius;
global L;

N=length(pos);
    for i=1:N
        %temp vars, computation
        if i ~=N
            rel_p=pos(i)-pos(i+1);
            rel_v=vel(i)-vel(i+1);
            tc=(rel_p+2*radius)/rel_v; %minus is intentional, write this out with the particles. will see is correct. 
            %storage

        else %rightmost particle vs wall 
            rel_p=L-pos(N);
            rel_v=vel(i);
            tc=(rel_p+2*radius)/rel_v; %no minus sign
        end
        
        if tc<0
            tc=5;
        end
        % storage
        right_col.rel_pos(i)=rel_p;
        right_col.rel_vel(i)=rel_v;
        right_col.tc(i)=tc;
    end
    
    
    
    left_flag=0; %spaghetti. but we know all particles will collide with the one on the right, except the leftmost. 
    right_flag=0;
    %this flag tells us if its the leftmost particle colliding with wall.
    [min_tc,min_idx]=min([right_col.tc]); %this wasnt well thought through, probably didnt need left_coll at all. 
    
    left_part_wall_tc=-pos(1)/vel(1);
    if left_part_wall_tc<0
        left_part_wall_tc=100;
    end
    
    if left_part_wall_tc<min_tc
        min_tc=left_part_wall_tc; 
        min_idx=1;
        left_flag=1;
    end
    
    if min_idx==1 && vel(1) < 0 
        right_flag=1;
    end
    
    
end