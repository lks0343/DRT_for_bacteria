% the function is directly taken from the DRTtools code.
function out_A_im = assemble_A_im(freq,freq_anal,epsilon,rbf_type,data_type)

%   this function assembles the A_im matrix
%   we will assume that the tau vector is identical to the freq vector

%   first get number of frequencies
    N_freq = numel(freq);
    N_anal = length(freq_anal);  % # of analyzed frequencies

%   we compute if the frequencies are sufficiently log spaced
    std_diff_freq = std(diff(log10(freq)));
    mean_diff_freq = mean(diff(log10(freq)));

    std_diff_freq_anal = std(diff(log10(freq_anal)));
    mean_diff_freq_anal = mean(diff(log10(freq_anal)));

%   the define the A_re output matrix
    out_A_im_temp = zeros(N_freq,N_anal);
    out_A_im = zeros(N_freq, N_anal+3);

%   if they are, we apply the toeplitz trick
    toeplitz_trick = std_diff_freq/mean_diff_freq<0.01 ...
        && mean_diff_freq_anal==mean_diff_freq;
    
%   if terms are evenly distributed & do not use PWL & diff_N_freq=diff_N_anal
%   then we do compute only N terms
%   else we compute all terms with brute force
    if toeplitz_trick && ~strcmp(rbf_type,'Piecewise linear')

        % define vectors R and C
        R = zeros(1,N_anal);
        C = zeros(N_freq,1);

        % for clarity the C and R computations are separated
        for iter_freq_n = 1: N_freq

            freq_n = freq(iter_freq_n);
            freq_m = freq_anal(1);
            C(iter_freq_n, 1) = g_ii(freq_n, freq_m, epsilon, rbf_type,data_type);

        end  

        for iter_freq_m = 1: N_anal

            freq_n = freq(1);
            freq_m = freq_anal(iter_freq_m);
            R(1, iter_freq_m) = g_ii(freq_n, freq_m, epsilon, rbf_type,data_type);

        end

        out_A_im_temp = toeplitz(C,R);

    else 
    % compute using brute force

        for iter_freq_n = 1: N_freq

            for iter_freq_m = 1: N_anal

                freq_n = freq(iter_freq_n);
                freq_m = freq_anal(iter_freq_m);

                % this is the usual PWL approximation
                if strcmp(rbf_type,'Piecewise linear')

                    if iter_freq_m == 1

                        freq_m_plus_1 = freq_anal(iter_freq_m+1);
                        out_A_im_temp(iter_freq_n, iter_freq_m) = -0.5*(2*pi*freq_n/freq_m)/(1+((2*pi*freq_n/freq_m))^2)*log((1/freq_m_plus_1)/(1/freq_m));

                    elseif iter_freq_m == N_anal

                        freq_m_minus_1 = freq_anal(iter_freq_m-1);
                        out_A_im_temp(iter_freq_n, iter_freq_m) = -0.5*(2*pi*freq_n/freq_m)/(1+((2*pi*freq_n/freq_m))^2)*log((1/freq_m)/((1/freq_m_minus_1)));

                    else

                        freq_m_plus_1 = freq_anal(iter_freq_m+1);
                        freq_m_minus_1 = freq_anal(iter_freq_m-1);
                        out_A_im_temp(iter_freq_n, iter_freq_m) = -0.5*(2*pi*freq_n/freq_m)/(1+((2*pi*freq_n/freq_m))^2)*log((1/freq_m_plus_1)/(1/freq_m_minus_1));

                    end

                else

                   % compute all RBF terms
                   out_A_im_temp(iter_freq_n, iter_freq_m) = g_ii(freq_n, freq_m, epsilon, rbf_type, data_type);

                end

            end
            
        end
    
    end

%   the first and second row are reserved for L and R respectively
    out_A_im(:, 4:end) = out_A_im_temp;

end
