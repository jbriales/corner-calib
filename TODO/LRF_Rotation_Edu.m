% Closed form for LRF calibration rotation

% Generate synthetic data:
for k=1:3
    v{k} = snormalize( rand(2,1) );
    w{k} = snormalize( rand(2,1) );
    a{k} = snormalize( rand(2,1) );
    b{k} = snormalize( rand(2,1) );
end