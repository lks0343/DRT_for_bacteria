% the function was made by KS Lee
function [re_im,discrepancy,GCV,mGCV,rGCV,x_re,x_im,x_comb,A_re,A_im,normal] = ...
    hyperDRT(freq_original,data_re,data_im,lambda_0,beta,coeff,shape_control,rbf_type,data_type,hyper_type)

%% reconstruct data as a form that frequencies aligned descending
k_max = 1e3;
% (# of analyzed frequencies) = 5(# of measure frequencies)
multi = 1;

% eps = 1e-6;
if freq_original(1,1) < freq_original(end,1)
    freq_original = fliplr(freq_original')';
    data_re = fliplr(data_re')';
    data_im = fliplr(data_im')';
end

%% parameters
Nf = length(freq_original);
mean_diff_freq = mean(diff(log10(freq_original)));

freq_anal = log10(freq_original(1))+2:mean_diff_freq/multi:log10(freq_original(end))-2;
freq_anal = 10.^freq_anal';
Nfexpand = length(freq_anal);

tau_max = ceil(max(log10(1./freq_anal))) + .5;
tau_min = floor(min(log10(1./freq_anal))) - .5;
freq_fine = logspace(-tau_min,-tau_max,10*Nfexpand);
lb = zeros(Nfexpand+2,1);
ub = Inf * ones(Nfexpand+2,1);
x0 = ones(size(lb));

W_re = ones(size(data_re));
W_re = diag(W_re);
W_im = ones(size(data_im));
W_im = diag(W_im);
b_re = data_re;
b_im = data_im;

epsilon = compute_epsilon(freq_anal,coeff,rbf_type,shape_control);     

%% simple regression to obtain the initial x
A_re = assemble_A_re(freq_original,freq_anal,epsilon,rbf_type,data_type);
A_im = assemble_A_im(freq_original,freq_anal,epsilon,rbf_type,data_type);
switch data_type
    case 'Z'
        A_re(:,3) = 1;
        A_im(:,2) = -1./(2*pi*freq_original);
        eps = 1e-3;
        options = optimset('algorithm','interior-point-convex','Display','off', ...
            'TolFun',1e-15,'TolX',1e-10,'MaxFunEvals',1e5);
    case 'Y'
        A_re(:,3) = 1;
        A_im(:,2) = 1./(2*pi*freq_original);
        A_im(:,1) = -2*pi*freq_original;
        eps = 1e-9;
        options = optimset('algorithm','interior-point-convex','Display','off', ...
            'TolFun',1e-17,'TolX',1e-12,'MaxFunEvals',1e5);
    case 'C'
        A_re(:,3) = 1;
        A_im(:,1) = 1./(2*pi*freq_original);
        eps = 1e-17;
        options = optimset('algorithm','interior-point-convex','Display','off', ...
            'TolFun',1e-17,'TolX',1e-17,'MaxFunEvals',1e5);
end

WA_re = W_re*A_re;
WA_im = W_im*A_im;

M = assemble_M_2(freq_anal,epsilon,rbf_type);      % using 2nd order derivative
L = zeros(Nfexpand+3,Nfexpand+3);
L(4:end,4:end) = chol(M(4:end,4:end));

[H_im,f_im] = quad_format(WA_im,b_im,M,1e-7);
[H_re,f_re] = quad_format(WA_re,b_re,M,1e-7);
[H_comb,f_comb] = quad_format_combined(WA_re,WA_im,b_re,b_im,M,1e-7);

x_re = quadprog(H_re,f_re,[],[],[],[],lb,ub,x0,options);
x_im = quadprog(H_im,f_im,[],[],[],[],lb,ub,x0,options);
x_comb = quadprog(H_comb,f_comb,[],[],[],[],lb,ub,x0,options);

%% hyper ridge regression
k = 0;
continue_loop = 1;
while continue_loop
    
    Lx = L*x_re;
    LAMBDA = assemble_LAMBDA(Lx,beta,lambda_0,hyper_type);
    M_hyper = L'*LAMBDA*L;
    x_re_temp = x_re;
    [H_re,f_re] = quad_format(WA_re,b_re,M_hyper,1);
    x_re = quadprog(H_re,f_re,[],[],[],[],lb,ub,x_re_temp,options);
    if (k>=k_max)||~sum(abs(x_re-x_re_temp)>eps)
        continue_loop = 0;
    end
    k = k+1;
end

k = 0;
continue_loop = 1;
while continue_loop
    
    Lx = L*x_im;
    LAMBDA = assemble_LAMBDA(Lx,beta,lambda_0,hyper_type);
    M_hyper = L'*LAMBDA*L;
    x_im_temp = x_im;
    [H_im,f_im] = quad_format(WA_im,b_im,M_hyper,1);
    x_im = quadprog(H_im,f_im,[],[],[],[],lb,ub,x_im_temp,options);
    if (k>=k_max)||~sum(abs(x_im-x_im_temp)>eps)
        continue_loop = 0;
    end
    k = k+1;
end

k = 0;
continue_loop = 1;
while continue_loop
    
    Lx = L*x_comb;
    LAMBDA = assemble_LAMBDA(Lx,beta,lambda_0,hyper_type);
    M_hyper = L'*LAMBDA*L;
    x_comb_temp = x_comb;
    [H_comb,f_comb] = quad_format_combined(WA_re,WA_im,b_re,b_im,M_hyper,1);
    x_comb = quadprog(H_comb,f_comb,[],[],[],[],lb,ub,x_comb_temp,options);
    if (k>=k_max)||~sum(abs(x_comb-x_comb_temp)>eps)
        continue_loop = 0;
    end
    k = k+1;
end

res_ri = (b_re - A_re*x_im)./sqrt(b_re.^2+b_im.^2);
res_ir = (b_im - A_im*x_re)./sqrt(b_re.^2+b_im.^2);
re_im = norm(res_ri)^2 + norm(res_ir)^2;
normal = normalized_factor(epsilon,rbf_type,x_comb(4:end));

normal_re = normalized_factor(epsilon,rbf_type,x_re(4:end));
normal_im = normalized_factor(epsilon,rbf_type,x_im(4:end));
discrepancy = norm(x_re/normal_re - x_im/normal_im)^2;

A = [A_re; A_im];
Z = [data_re; data_im];
K_matrix = A*inv(A'*A + M_hyper)*A';
GCV = (1/(2*Nfexpand))*norm((eye(size(K_matrix))-K_matrix)*Z)^2/ ...
    (1/(2*Nfexpand)*trace(eye(size(K_matrix))-K_matrix))^2;

if Nfexpand<50
    rho = 1.3;
else
    rho = 2;
end
mGCV = (1/(2*Nfexpand))*norm((eye(size(K_matrix))-K_matrix)*Z)^2/ ...
    (1/(2*Nfexpand)*trace(eye(size(K_matrix))-rho*K_matrix))^2;

mu2 = 1/(2*Nfexpand)*trace(K_matrix^2);
if Nfexpand<50
    zeta = 0.2;
else
    zeta = 0.3;
end
rGCV = (zeta + (1-zeta)*mu2)*GCV;

end
