function out = farthest(x,c)
% Return the index of c corresponding to the farther Euclidean distance
% from x

    if( abs(x-c(1)) > abs(x-c(2)) )
        out = 1;
    else
        out = 2;
    end

end