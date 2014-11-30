q = [co.LRF_q{:}]
V = co.cam_R_c_w
cam_c = co.cam_c_ray
[R, t_] = P3P_LRF( q, V );

tri_n = co.cam_R_c_w(:,co.thereis_LRF_v)
LRF_v = [co.LRF_v{:}]
for k=1:8
    err2(k) = sum( dot( tri_n, R{k}(:,1:2)*LRF_v, 1 ).^2 );
    disp(k)
    disp( err2(end) )
end
[~,I] = min( err2 );
R = { R{I}, R{I}*diag([-1 -1 +1]) };
t_ = { t_{I}, -t_{I} };

% Check with translation
bp_n = co.cam_reprN(:,co.thereis_LRF_q)
LRF_q = [co.LRF_q{:}]
lam = zeros(1,2);
for k=1:2
    lam(k) = -(bp_n(:,1)'*(R{k}(:,1:2)*LRF_q(:,1) + t_{k})) / ...
        (bp_n(:,1)'*cam_c);
    fprintf('lambda value for %i: %f\n',k,lam(k))
end
id = find( lam > 0 );
lam = lam(id);
t = t_{id} + lam * cam_c;
R = R{id};
