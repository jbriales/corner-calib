function h = plotTranslation_3D_CostFunction( obj, R, t )
% TODO: Weights

dist  = 0.01;    % Simulation distance (m)
inc   = 0.0005;   % Increments

gv  = -dist:inc:+dist;
Ngv = length(gv);
[w_i,w_j] = meshgrid( gv );

% Linearize code
inc_eps = [ w_i(:), w_j(:), zeros(Ngv^2,1) ]';

weights = obj.FWeights_3D_PlaneDistance( R, t );

labels = { 'x', 'y', 'z'
    'y', 'z', 'x' };
titles = { 'X-Y err^2', 'Y-Z err^2', 'Z-X err^2'
           'X-Y errW' , 'Y-Z errW' , 'Z-X errW' };
figure, hold on
for k=1:3
    Err2 = zeros(1,Ngv^2);
    ErrW = zeros(1,Ngv^2);
    inc_eps_ = circshift( inc_eps, k-1 ); % Shift inc_eps rows
    for i=1:numel(Err2)
        t_ = inc_eps_(:,i) + t;
        residual = obj.FErr_3D_PlaneDistance( R, t_ );
        Err2(i) = residual' * residual;
        ErrW(i) = residual' * weights * residual;
    end
    Err2 = reshape(Err2,Ngv,Ngv);
    ErrW = reshape(ErrW,Ngv,Ngv);
    Err = {Err2, ErrW};
    
    res_GT = obj.FErr_3D_PlaneDistance( R, t );
    err2_GT = { res_GT'*res_GT, res_GT'*weights*res_GT };
    for k_sub=1:2
        subplot(2,3,k +(k_sub-1)*3 ); hold on;
        xlabel(labels{1,k});
        ylabel(labels{2,k});
        title( titles{k_sub,k} );
        surf(w_i,w_j, Err{k_sub});
        % Plot GT point
        plot3(0,0,err2_GT{k_sub}, '.y', 'LineWidth', 3);
        axis([-dist +dist -dist +dist]);
        shading interp;
        view([90 90]);
    end
end
h = []; % TODO
end