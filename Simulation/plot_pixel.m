function plot_pixel( results, sd_pixel, N, sd_laser )

% Define the variables
K = size(sd_pixel,2);
mean_err = zeros(K,6);   %[R_w R_nw t_w t_nw t_w3 t_nw3]
var_err  = zeros(K,6);

% Assign the GT values 
R_c_s = results(1,1,1,1).gt.R_c_s;
t_c_s = results(1,1,1,1).gt.t_c_s;

n_n   = results(N,1,sd_laser,1).param.N;
s_l   = results(N,1,sd_laser,1).param.sd_laser;

% Estimates the errors
for i = 1:size(sd_pixel,2)
    for k = 1:size(results,4)
        err_(k,1)   = abs(angularDistance(results(N,i,sd_laser,k).rot.R_c_s_w, R_c_s));
        err_(k,2)   = abs(angularDistance(results(N,i,sd_laser,k).rot.R_c_s_nw, R_c_s));
        err_(k,3)   = norm(results(N,i,sd_laser,k).trans2D.t_c_s_w - t_c_s);
        err_(k,4)   = norm(results(N,i,sd_laser,k).trans2D.t_c_s_nw - t_c_s);
        err_(k,5)   = norm(results(N,i,sd_laser,k).trans3D.t_c_s_w - t_c_s);
        err_(k,6)   = norm(results(N,i,sd_laser,k).trans3D.t_c_s_nw - t_c_s);
    end
    mean_err(i,:)   = mean(err_);
    var_err(i,:)    = var(err_);
end

% Plot the rotation
h_r = figure;
h_r = barwitherr(var_err(:,1:2),mean_err(:,1:2),'BarWidth',0.6);
% %set(h_r,'FaceColor',[0.1 0.6 0.2; 0.1 0.6 0.5] );
legend('Weighted', 'Non-weighted');
title_str = strcat('N = ',int2str(n_n),32, 32, 32, '\sigma_{laser} = ',num2str(s_l,'%f') );
title(title_str);
xlabel('\sigma_{pixel}');
ylabel('Rotation error (º)');
x_label = {};
for i = 1:size(sd_pixel,2)
    x_label{i} = num2str(sd_pixel(i),'%f');
end
set(gca,'XTickLabel',x_label)


% Plot the 2D translation
h_t2 = figure;
h_t2 = barwitherr(var_err(:,3:4),mean_err(:,3:4),'BarWidth',0.6);
% %set(h_r,'FaceColor',[0.1 0.6 0.2; 0.1 0.6 0.5] );
legend('Weighted', 'Non-weighted');
title(title_str);
xlabel('\sigma_{pixel}');
ylabel('Translation error 2D solution (m)');
% x_label = {};
% for i = 1:size(sd_pixel,2)
%     x_label{i} = int2str(sd_pixel(i));
% end
set(gca,'XTickLabel',x_label)


% Plot the 3D translation
h_t3 = figure;
h_t3 = barwitherr(var_err(:,5:6),mean_err(:,5:6),'BarWidth',0.6);
% %set(h_r,'FaceColor',[0.1 0.6 0.2; 0.1 0.6 0.5] );
legend('Weighted', 'Non-weighted');
title(title_str);
xlabel('\sigma_{pixel}');
ylabel('Translation error 3D solution (m)');
% x_label = {};
% for i = 1:size(sd_pixel,2)
%     x_label{i} = int2str(sd_pixel(i));
% end
set(gca,'XTickLabel',x_label)







