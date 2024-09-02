% the function was made by KS Lee
function LAMBDA = assemble_LAMBDA(Lx,beta,lambda_0,hyper_type)
j_leng = length(Lx);
lambda_j = zeros(j_leng,1);

switch hyper_type
    case 'Exponential'
        eta = 1/lambda_0;
        lambda_j = 1./((Lx).^2 + eta);
    case 'Gaussian'
        eta = 1/(2*lambda_0^2);
        lambda_j = (sqrt(Lx.^4+8*eta) - Lx.^2)/(4*eta);
    case 'Gamma'
        eta = (beta-1)/lambda_0;
        lambda_j = (beta-1)./(Lx.^2 + eta);
    case 'Inverse Gamma'
        eta = lambda_0*(beta+1);
        lambda_j = (sqrt((beta+1)^2 + 4*eta.*Lx.^2) - (beta+1))./(2.*Lx.^2);
        lambda_j(1:2) = 0;
end

LAMBDA = diag(lambda_j);
        
