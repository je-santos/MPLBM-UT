function cube = p_oct_label(octant, label, cube)

% check if there are points in the octant with value 1
if( octant==1 )
    
    % set points in this octant to current label
    % and recurseive labeling of adjacent octants
    idx = cube(:,1) == 1;
    if any(idx)
        cube(idx,1) = label(idx);
    end
    
    idx = cube(:,2) == 1;
    if any(idx)
        cube(idx,2) = label(idx);
        cube(idx,:) = p_oct_label(2,label(idx),cube(idx,:));
    end
    
    idx = cube(:,4) == 1;
    if any(idx)
        cube(idx,4) = label(idx);
        cube(idx,:) = p_oct_label(3,label(idx),cube(idx,:));
    end
    
    idx = cube(:,5) == 1;
    if any(idx)
        cube(idx,5) = label(idx);
        cube(idx,:) = p_oct_label(2,label(idx),cube(idx,:));
        cube(idx,:) = p_oct_label(3,label(idx),cube(idx,:));
        cube(idx,:) = p_oct_label(4,label(idx),cube(idx,:));
    end
    
    idx = cube(:,10) == 1;
    if any(idx)
        cube(idx,10) = label(idx);
        cube(idx,:) = p_oct_label(5,label(idx),cube(idx,:));
    end
    
    idx = cube(:,11) == 1;
    if any(idx)
        cube(idx,11) = label(idx);
        cube(idx,:) = p_oct_label(2,label(idx),cube(idx,:));
        cube(idx,:) = p_oct_label(5,label(idx),cube(idx,:));
        cube(idx,:) = p_oct_label(6,label(idx),cube(idx,:));
    end
    
    idx = cube(:,13) == 1;
    if any(idx)
        cube(idx,13) = label(idx);
        cube(idx,:) = p_oct_label(3,label(idx),cube(idx,:));
        cube(idx,:) = p_oct_label(5,label(idx),cube(idx,:));
        cube(idx,:) = p_oct_label(7,label(idx),cube(idx,:));
    end
    
end

if( octant==2 )
    
    idx = cube(:,2) == 1;
    if any(idx)
        cube(idx,2) = label(idx);
        cube(idx,:) = p_oct_label(1,label(idx),cube(idx,:));
    end

    idx = cube(:,5) == 1;
    if any(idx)
        cube(idx,5) = label(idx);
        cube(idx,:) = p_oct_label(1,label(idx),cube(idx,:));
        cube(idx,:) = p_oct_label(3,label(idx),cube(idx,:));
        cube(idx,:) = p_oct_label(4,label(idx),cube(idx,:));
    end

    idx = cube(:,11) == 1;
    if any(idx)
        cube(idx,11) = label(idx);
        cube(idx,:) = p_oct_label(1,label(idx),cube(idx,:));
        cube(idx,:) = p_oct_label(5,label(idx),cube(idx,:));
        cube(idx,:) = p_oct_label(6,label(idx),cube(idx,:));
    end

    idx = cube(:,3) == 1;
    if any(idx)
        cube(idx,3) = label(idx);
    end

    idx = cube(:,6) == 1;
    if any(idx)
        cube(idx,6) = label(idx);
        cube(idx,:) = p_oct_label(4,label(idx),cube(idx,:));
    end
    
    idx = cube(:,12) == 1;
    if any(idx)
        cube(idx,12) = label(idx);
        cube(idx,:) = p_oct_label(6,label(idx),cube(idx,:));
    end

    idx = cube(:,14) == 1;
    if any(idx)
        cube(idx,14) = label(idx);
        cube(idx,:) = p_oct_label(4,label(idx),cube(idx,:));
        cube(idx,:) = p_oct_label(6,label(idx),cube(idx,:));
        cube(idx,:) = p_oct_label(8,label(idx),cube(idx,:));
    end

end

if( octant==3 )
    
    idx = cube(:,4) == 1;
    if any(idx)
        cube(idx,4) = label(idx);
        cube(idx,:) = p_oct_label(1,label(idx),cube(idx,:));
    end

    idx = cube(:,5) == 1;
    if any(idx)
        cube(idx,5) = label(idx);
        cube(idx,:) = p_oct_label(1,label(idx),cube(idx,:));
        cube(idx,:) = p_oct_label(2,label(idx),cube(idx,:));
        cube(idx,:) = p_oct_label(4,label(idx),cube(idx,:));
    end

    idx = cube(:,13) == 1;
    if any(idx)
        cube(idx,13) = label(idx);
        cube(idx,:) = p_oct_label(1,label(idx),cube(idx,:));
        cube(idx,:) = p_oct_label(5,label(idx),cube(idx,:));
        cube(idx,:) = p_oct_label(7,label(idx),cube(idx,:));
    end

    idx = cube(:,7) == 1;
    if any(idx)
        cube(idx,7) = label(idx);
    end

    idx = cube(:,8) == 1;
    if any(idx)
        cube(idx,8) = label(idx);
        cube(idx,:) = p_oct_label(4,label(idx),cube(idx,:));
    end
    
    idx = cube(:,15) == 1;
    if any(idx)
        cube(idx,15) = label(idx);
        cube(idx,:) = p_oct_label(7,label(idx),cube(idx,:));
    end

    idx = cube(:,16) == 1;
    if any(idx)
        cube(idx,16) = label(idx);
        cube(idx,:) = p_oct_label(4,label(idx),cube(idx,:));
        cube(idx,:) = p_oct_label(7,label(idx),cube(idx,:));
        cube(idx,:) = p_oct_label(8,label(idx),cube(idx,:));
    end
    
end

if( octant==4 )
    
    idx = cube(:,5) == 1;
    if any(idx)
        cube(idx,5) = label(idx);
        cube(idx,:) = p_oct_label(1,label(idx),cube(idx,:));
        cube(idx,:) = p_oct_label(2,label(idx),cube(idx,:));
        cube(idx,:) = p_oct_label(3,label(idx),cube(idx,:));
    end

    idx = cube(:,6) == 1;
    if any(idx)
        cube(idx,6) = label(idx);
        cube(idx,:) = p_oct_label(2,label(idx),cube(idx,:));
    end

    idx = cube(:,14) == 1;
    if any(idx)
        cube(idx,14) = label(idx);
        cube(idx,:) = p_oct_label(2,label(idx),cube(idx,:));
        cube(idx,:) = p_oct_label(6,label(idx),cube(idx,:));
        cube(idx,:) = p_oct_label(8,label(idx),cube(idx,:));
    end
    
    idx = cube(:,8) == 1;
    if any(idx)
        cube(idx,8) = label(idx);
        cube(idx,:) = p_oct_label(3,label(idx),cube(idx,:));
    end

    idx = cube(:,16) == 1;
    if any(idx)
        cube(idx,16) = label(idx);
        cube(idx,:) = p_oct_label(3,label(idx),cube(idx,:));
        cube(idx,:) = p_oct_label(7,label(idx),cube(idx,:));
        cube(idx,:) = p_oct_label(8,label(idx),cube(idx,:));
    end

    idx = cube(:,9) == 1;
    if any(idx)
        cube(idx,9) = label(idx);
    end

    idx = cube(:,17) == 1;
    if any(idx)
        cube(idx,17) = label(idx);
        cube(idx,:) = p_oct_label(8,label(idx),cube(idx,:));
    end

end

if( octant==5 )
    
    idx = cube(:,10) == 1;
    if any(idx)
        cube(idx,10) = label(idx);
        cube(idx,:) = p_oct_label(1,label(idx),cube(idx,:));
    end

    idx = cube(:,11) == 1;
    if any(idx)
        cube(idx,11) = label(idx);
        cube(idx,:) = p_oct_label(1,label(idx),cube(idx,:));
        cube(idx,:) = p_oct_label(2,label(idx),cube(idx,:));
        cube(idx,:) = p_oct_label(6,label(idx),cube(idx,:));
    end
    
    idx = cube(:,13) == 1;
    if any(idx)
        cube(idx,13) = label(idx);
        cube(idx,:) = p_oct_label(1,label(idx),cube(idx,:));
        cube(idx,:) = p_oct_label(3,label(idx),cube(idx,:));
        cube(idx,:) = p_oct_label(7,label(idx),cube(idx,:));
    end

    idx = cube(:,18) == 1;
    if any(idx)
        cube(idx,18) = label(idx);
    end

    idx = cube(:,19) == 1;
    if any(idx)
        cube(idx,19) = label(idx);
        cube(idx,:) = p_oct_label(6,label(idx),cube(idx,:));
    end

    idx = cube(:,21) == 1;
    if any(idx)
        cube(idx,21) = label(idx);
        cube(idx,:) = p_oct_label(7,label(idx),cube(idx,:));
    end

    idx = cube(:,22) == 1;
    if any(idx)
        cube(idx,22) = label(idx);
        cube(idx,:) = p_oct_label(6,label(idx),cube(idx,:));
        cube(idx,:) = p_oct_label(7,label(idx),cube(idx,:));
        cube(idx,:) = p_oct_label(8,label(idx),cube(idx,:));
    end

end

if( octant==6 )
    
    idx = cube(:,11) == 1;
    if any(idx)
        cube(idx,11) = label(idx);
        cube(idx,:) = p_oct_label(1,label(idx),cube(idx,:));
        cube(idx,:) = p_oct_label(2,label(idx),cube(idx,:));
        cube(idx,:) = p_oct_label(5,label(idx),cube(idx,:));
    end

    idx = cube(:,12) == 1;
    if any(idx)
        cube(idx,12) = label(idx);
        cube(idx,:) = p_oct_label(2,label(idx),cube(idx,:));
    end

    idx = cube(:,14) == 1;
    if any(idx)
        cube(idx,14) = label(idx);
        cube(idx,:) = p_oct_label(2,label(idx),cube(idx,:));
        cube(idx,:) = p_oct_label(4,label(idx),cube(idx,:));
        cube(idx,:) = p_oct_label(8,label(idx),cube(idx,:));
    end
    
    idx = cube(:,19) == 1;
    if any(idx)
        cube(idx,19) = label(idx);
        cube(idx,:) = p_oct_label(5,label(idx),cube(idx,:));
    end


    idx = cube(:,22) == 1;
    if any(idx)
        cube(idx,22) = label(idx);
        cube(idx,:) = p_oct_label(5,label(idx),cube(idx,:));
        cube(idx,:) = p_oct_label(7,label(idx),cube(idx,:));
        cube(idx,:) = p_oct_label(8,label(idx),cube(idx,:));
    end
    
    idx = cube(:,20) == 1;
    if any(idx)
        cube(idx,20) = label(idx);
    end

    idx = cube(:,23) == 1;
    if any(idx)
        cube(idx,23) = label(idx);
        cube(idx,:) = p_oct_label(8,label(idx),cube(idx,:));
    end
 
end

if( octant==7 )
    
    idx = cube(:,13) == 1;
    if any(idx)
        cube(idx,13) = label(idx);
        cube(idx,:) = p_oct_label(1,label(idx),cube(idx,:));
        cube(idx,:) = p_oct_label(3,label(idx),cube(idx,:));
        cube(idx,:) = p_oct_label(5,label(idx),cube(idx,:));
    end

    idx = cube(:,15) == 1;
    if any(idx)
        cube(idx,15) = label(idx);
        cube(idx,:) = p_oct_label(3,label(idx),cube(idx,:));
    end

    idx = cube(:,16) == 1;
    if any(idx)
        cube(idx,16) = label(idx);
        cube(idx,:) = p_oct_label(3,label(idx),cube(idx,:));
        cube(idx,:) = p_oct_label(4,label(idx),cube(idx,:));
        cube(idx,:) = p_oct_label(8,label(idx),cube(idx,:));
    end

    idx = cube(:,21) == 1;
    if any(idx)
        cube(idx,21) = label(idx);
        cube(idx,:) = p_oct_label(5,label(idx),cube(idx,:));
    end

    idx = cube(:,22) == 1;
    if any(idx)
        cube(idx,22) = label(idx);
        cube(idx,:) = p_oct_label(5,label(idx),cube(idx,:));
        cube(idx,:) = p_oct_label(6,label(idx),cube(idx,:));
        cube(idx,:) = p_oct_label(8,label(idx),cube(idx,:));
    end

    idx = cube(:,24) == 1;
    if any(idx)
        cube(idx,24) = label(idx);
    end
    
    idx = cube(:,25) == 1;
    if any(idx)
        cube(idx,25) = label(idx);
        cube(idx,:) = p_oct_label(8,label(idx),cube(idx,:));
    end
end

if( octant==8 )
    
    idx = cube(:,14) == 1;
    if any(idx)
        cube(idx,14) = label(idx);
        cube(idx,:) = p_oct_label(2,label(idx),cube(idx,:));
        cube(idx,:) = p_oct_label(4,label(idx),cube(idx,:));
        cube(idx,:) = p_oct_label(6,label(idx),cube(idx,:));
    end

    idx = cube(:,16) == 1;
    if any(idx)
        cube(idx,16) = label(idx);
        cube(idx,:) = p_oct_label(3,label(idx),cube(idx,:));
        cube(idx,:) = p_oct_label(4,label(idx),cube(idx,:));
        cube(idx,:) = p_oct_label(7,label(idx),cube(idx,:));
    end
    
    idx = cube(:,17) == 1;
    if any(idx)
        cube(idx,17) = label(idx);
        cube(idx,:) = p_oct_label(4,label(idx),cube(idx,:));
    end
    
    idx = cube(:,22) == 1;
    if any(idx)
        cube(idx,22) = label(idx);
        cube(idx,:) = p_oct_label(5,label(idx),cube(idx,:));
        cube(idx,:) = p_oct_label(6,label(idx),cube(idx,:));
        cube(idx,:) = p_oct_label(7,label(idx),cube(idx,:));
    end
    
    idx = cube(:,17) == 1;
    if any(idx)
        cube(idx,17) = label(idx);
        cube(idx,:) = p_oct_label(4,label(idx),cube(idx,:));
    end
    
    idx = cube(:,23) == 1;
    if any(idx)
        cube(idx,23) = label(idx);
        cube(idx,:) = p_oct_label(6,label(idx),cube(idx,:));
    end
    
    idx = cube(:,25) == 1;
    if any(idx)
        cube(idx,25) = label(idx);
        cube(idx,:) = p_oct_label(7,label(idx),cube(idx,:));
    end
    
    idx = cube(:,26) == 1;
    if any(idx)
        cube(idx,26) = label(idx);
    end
end
end
