% this function is modified DRTtools codes by KS Lee 
function [discrepancy,crossValidation,res_rr,res_ii,res_ri,res_ir,x_re,x_im,x_comb,A_re,A_im] = ...
    simpleDRT(freq,data_re,data_im,lambda,coeff,shape_control,rbf_type,data_type)

%% reconstruct data as a form that frequencies aligned descending
if freq(1,1) < freq(end,1)
    freq = fliplr(freq')';
    data_re = fliplr(data_re')';
    data_im = fliplr(data_im')';
end

%% parameters
Nf = length(freq);
inter = log10(freq(end)/freq(1))/(Nf-1);
freqExpand = [log10(freq(1))+2:inter:log10(freq(end))-2]';
freqExpand = 10.^freqExpand;
Nfexpand = length(freqExpand);

tau_max = ceil(max(log10(1./freqExpand))) + 0.5;
tau_min = floor(min(log10(1./freqExpand))) - 0.5;
freq_fine = logspace(-tau_min,-tau_max,10*Nfexpand);
lb = zeros(Nfexpand+2,1);
ub = Inf * ones(Nfexpand+2,1);
x0 = ones(size(lb));
options = optimset('algorithm', 'interior-point-convex', 'Display', ...
    'off', 'TolFun', 1e-15, 'TolX', 1e-10, 'MaxFunEvals', 1E5);

W_re = ones(size(data_re))/max(abs(data_re));
W_re = diag(W_re);

W_im = ones(size(data_im))/max(abs(data_im));
W_im = diag(W_im);

b_re = W_re * data_re;
b_im = W_im * data_im;

epsilon = compute_epsilon(freq,coeff,rbf_type,shape_control);

%% regression
A_re = assemble_A_re(freqExpand,epsilon,rbf_type,data_type);
A_im = assemble_A_im(freqExpand,epsilon,rbf_type,data_type);

switch data_type
    case 'Z'
        A_re(:,2) = 1;
        A_im(:,1) = -1./(2*pi*freqExpand);    % assume that circuit is connected with a capacitor (double layer) in series
    case 'Z_noCap'
        A_re(:,2) = 1;          % assume that there's no capacitance connected in seires
        A_im(:,1) = 2*pi*freqExpand;  % assume that there's coil connected in series
    case 'Y'
        A_re(:,2) = 1;
        A_im(:,1) = -2*pi*freqExpand;
    case 'C'
        A_re(:,2) = 1;
        A_im(:,1) = 1./(2*pi*freqExpand);
end

% M = assemble_M_1(freqExpand,epsilon,rbf_type);
M = assemble_M_2(freqExpand,epsilon,rbf_type);      % using 2nd order derivative

Nerase = Nfexpand - Nf;        
A_re(1:Nerase/2,:) = []; 
A_re(end+1-Nerase/2:end,:) = []; 
A_im(1:Nerase/2,:) = []; 
A_im(end+1-Nerase/2:end,:) = []; 

WA_re = W_re*A_re;
WA_im = W_im*A_im;

[H_im,f_im] = quad_format(WA_im,b_im,M,lambda);
[H_re,f_re] = quad_format(WA_re,b_re,M,lambda);
[H_comb,f_comb] = quad_format_combined(WA_re,WA_im,b_re,b_im,M,lambda);

x_re = quadprog(H_re,f_re,[],[],[],[],lb,ub,x0,options);
x_im = quadprog(H_im,f_im,[],[],[],[],lb,ub,x0,options);
x_comb = quadprog(H_comb,f_comb,[],[],[],[],lb,ub,x0,options);

res_rr = A_re*x_re - b_re;
res_ii = A_im*x_im - b_im;
res_ri = A_re(:,3:end)*x_im(3:end) - b_re + A_re(:,1:2)*x_re(1:2);
res_ir = A_im(:,3:end)*x_re(3:end) - b_im + A_im(:,1:2)*x_im(1:2);

discrepancy = (x_re(3:end)-x_im(3:end))' * (x_re(3:end)-x_im(3:end));
crossValidation = res_ri'*res_ri + res_ir'*res_ir;

end
