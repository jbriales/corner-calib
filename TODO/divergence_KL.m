function [D_KL, rel_norm] = divergence_KL(mu_ref, cov_ref, mu_comp, cov_comp)
% Returns the Kullback-Liebler divergence between the reference
% distribution and the compared distribution, and the relative error with
% the norm of the matrices
%
%[D_KL, rel_norm] = divergence_KL(mu_ref, cov_ref, mu_comp, cov_comp)

% k = size(mu_ref,1);
% 
% D_KL = (-log(det(cov_comp)/det(cov_ref))+trace(cov_ref \ cov_comp)-k+(mu_comp-mu_ref)'*(cov_comp \ (mu_comp-mu_ref)))/2;
% 
% rel_norm = norm(cov_ref-cov_comp)/norm(cov_ref);

k = size(mu_ref,1);

D_KL = (-log(det(cov_comp)/det(cov_ref))+trace(inv(cov_ref)*cov_comp)-k+(mu_comp-mu_ref)'*inv(cov_comp)*(mu_comp-mu_ref))/2;

rel_norm = norm(cov_ref-cov_comp)/norm(cov_ref);