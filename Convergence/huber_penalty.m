function [H, d_H]  =  huber_penalty(d, d_max)

for i = 1:3
    
    if abs(d(i)) >= d_max
        H(i,1)   = d(i) ^ 2 ;
        d_H(i,1) = 2 * d(i) ;
    else
        H(i,1)   = d_max * ( 2 * abs(d(i)) - d_max ) ;
        d_H(i,1) = 2 * d_max * sign(d(i)) ;  
    end
    
end

end