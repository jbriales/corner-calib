function rho = chooseTrihedronNormals( rho, OmU, xi, K )

Dproj = DprojectionP2( xi.c.X );
v_im  = { xi.v(1).X , xi.v(2).X , xi.v(3).X };

% Two possible signs solutions:
rho_{1} = rho;
rho_{2} = diag([-1 1]) * rho;

% Check projection sign and rotation determinant
for k=1:2
    rho = rho_{k};
    n = cell(1,3);
    for i=1:3
        s = sign( v_im{i}' * Dproj * K * OmU{i} * rho(:,i) );
        rho(:,i) = s * rho(:,i);
        n{i} = OmU{i} * rho(:,i);
    end
    R_tri = cell2mat( n );
    
    if abs( det(R_tri) - 1 ) < 0.5
        break % Exit loop with current rho output value
    end
end