% the function was made by KS Lee

function [res,lambda_ind,recon_re,recon_im,res_freq_depen] ...
    = compute_residue(Z_re,Z_im,A_re,A_im,x,lambda,dataType)
% function which calculate the residue and reconstruct the Z (or C, or Y)
% using the DRT
% res : residue , resRatio : ratio of residue to measurement
% lambda_ind : index of the lambda of which has minimum res
% lambda_ind_ratio : index of the lambda of which has minimum resRatio
% recon_re, recon_im : reconstructed Z (or C, or Y) using DRT data, which has minimum res
% recon_re_ratio, recon_im_ratio : reconstructed Z (or C, or Y) using DRT data,
% which has minimum resRatio
% x is x_comb, lambda is vector of lambdas

res = [];
switch dataType{2}
    case 'absolute'
        for i = 1:numel(lambda)
            res_re = (Z_re - A_re{i}*x(:,i));
            res_im = (Z_im - A_im{i}*x(:,i));
            res = [res; sqrt(res_re'*res_re + res_im'*res_im)];
%             res = [res; sqrt(res_re'*res_re)+sqrt(res_im'*res_im)];
%             res = [res; sqrt((res_re+res_im)'*(res_re+res_im))];
        end
        lambda_ind = max(find(res==min(res)));
        res_re = (Z_re - A_re{lambda_ind}*x(:,lambda_ind));
        res_im = (Z_im - A_im{lambda_ind}*x(:,lambda_ind));
        
    case 'ratio'
        for i = 1:numel(lambda)
            res_re = (Z_re - A_re{i}*x(:,i))./Z_re;
            res_im = (Z_im - A_im{i}*x(:,i))./Z_im;
            res = [res; sqrt(res_re'*res_re + res_im'*res_im)];
%             res = [res; sqrt(res_re'*res_re)+sqrt(res_im'*res_im)];
%             res = [res; sqrt((res_re+res_im)'*(res_re+res_im))];
        end
        lambda_ind = max(find(res==min(res)));
        res_re = (Z_re - A_re{lambda_ind}*x(:,i))./Z_re;
        res_im = (Z_im - A_im{lambda_ind}*x(:,lambda_ind))./Z_im;
end
%% reconstruct using DRT
recon_re = A_re{lambda_ind}*x(:,lambda_ind);
recon_im = A_im{lambda_ind}*x(:,lambda_ind);
res_freq_depen = sqrt(res_re.^2 + res_im.^2);  %%%%%%%% 고쳐볼 수도....?
% res_freq_depen = sqrt(res_re.^2) + sqrt(res_im.^2);
% res_freq_depen = sqrt((res_re+res_im).^2);
end
