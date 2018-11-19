function [pos, vel] = parsim_solve_collision(pos,vel,tc,min_idx,left_flag, right_flag)
    %SOLVE FOR COLLISION. 
    global radius;
    global m_a;
    global m_b;
    pos=pos+tc*vel;
    
    
    if left_flag>0.5 || right_flag>0.5
        pos(min_idx)=radius;
        vel(min_idx)=-vel(min_idx);
    else
        %boozer eq 1.
        i=min_idx; %reuse variable, out of for loop. 
        %does not actually change size.. 
        try
        vel(i)=((m(i)-m(i+1))*vel(i)+2*m(i+1)*vel(i+1))/(m_a+m_b);
        vel(i)=((m(i+1)-m(i))*vel(i+1)+2*m(i+1)*vel(i+1))/(m_a+m_b);
        catch
            x=0
        end
        
    end
end

