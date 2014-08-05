function [residual, J] = Fun( rot_input, R )
% [Res, J] = Fun( R, N, L )

% N_obs = length(rot_input);

% residual = zeros( size([rot_input.N],2), 1);
% J = zeros( size([rot_input.N],2), 3); % If only rotation SO(3)

% cont = 1;
% for i=1:N_obs
%     N_cam = rot_input(i).N;
%     LINE_LRF = rot_input(i).l;
%     
%     for k=1:size(N_cam,2)
%         n_cam = N_cam(:,k);
%         line_LRF = LINE_LRF(:,k);
%         
%         residual(cont) = n_cam' * R(1:3,1:2) * line_LRF;
%         J(i,:) = cross( R(1:3,1:2) * line_LRF, n_cam )';
%         
%         cont = cont + 1;
%     end
% end
% Parallel
n_cam = [rot_input.N];
line_LRF = [rot_input.l];
residual = dot( n_cam, R(1:3,1:2) * line_LRF, 1 )';
J = cross( R(1:3,1:2) * line_LRF, n_cam, 1 )';

end