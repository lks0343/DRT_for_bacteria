% clear all; close all
global time col freq indNC indCon1 indCon2 indCon3 indCon4
% p = parpool(10);
col = [153 153 153; 200 0 0; 0 0 102; 000 051 000; 204 204 051; 153 204 255; ...
    102 000 102; 000 102 000; 204 153 102; 153 153 204; 153 153 102; 000 000 000; ...
    000 000 255; 102 051 153; 000 102 102; 255 102 000; 255 000 255; 051 204 051; ...
    102 051 000; 051 102 204; 102 153 051; 255 153 153; 153 255 204; 255 204 102; ...
    255 204 255; 204 255 204; 255 204 204; 255 255 204; 000 204 204; 204 255 000; 
    50 50 50; 100 100 100; 150 150 150; 200 200 200]/255;
numChan = 16;

indNC = [15];     % Negative Control channel index
indCon1 = [3,4];      % 10^4 CFU/ml channel index
indCon2 = [6,7,8];      % 10^3 CFU/ml channel index
indCon3 = [9,10,11];   % 10^2 CFU/ml channel index
indCon4 = [12,13,14];   % 10^1 CFU/ml channel index
% 
% foldPath = './20231017_polystyrene/';
shape_control = 'FWHM Coefficient';
rbf_type = 'Gaussian';
hyper_type = 'Gamma';
% anal_index = 64;
multi = 1;
% ì²˜ë¦¬?ï¿½ï¿½ ?ï¿½ï¿½?ï¿½ï¿½?ï¿½ï¿½?ï¿½ï¿½ ?ï¿½ï¿½?ï¿½ï¿½?ï¿½ï¿½?ï¿½ï¿½?ï¿½ï¿½ï¿??? ?ï¿½ï¿½?ï¿½ï¿½?ï¿½ï¿½ï¿??? ?ï¿½ï¿½?ï¿½ï¿½ vector
indCheck = [indNC, indCon1, indCon2, indCon3, indCon4];
%% import the data file
dListCGb = dir([foldPath, '*CG*']);
dListZDb = dir([foldPath '*ZD*']);

lengB = zeros(numChan,1);

for i = 1:numChan
    CGb{i} = dlmread([foldPath dListCGb(i).name]);
    ZDb{i} = dlmread([foldPath dListZDb(i).name]);
    lengB(i) = length(CGb{i});
end
freqN = find(CGb{1}(:,5));
freqN = freqN(2)-freqN(1);

lengB = min(lengB);     % the length of the bacterial data
totTime = dlmread([foldPath 'timeLog.txt']);
totTime = totTime(:,1)*3600 + totTime(:,2)*60 + totTime(:,3);
totTime = totTime(2) - totTime(1);
time = totTime/(lengB/freqN)/60/60; % time interval between datas
time = 0:time:(totTime/60/60);  % unit : hour
time = time';
time(end) = [];

startPoint = find(time<1/6);     % ì¸¡ì •?ï¿½ï¿½?ï¿½ï¿½ ?ï¿½ï¿½ 30ï¿??? ?ï¿½ï¿½?ï¿½ï¿½ï¿??? analysis
startPoint = 4;
time(end-startPoint+1:end) = [];
time(end) = [];

freq = CGb{1}(1:freqN,1);   % measuring frequency value 

% import the bacterial data along frequencies
C{freqN} = [];
G{freqN} = [];
Z{freqN} = [];
D{freqN} = [];

for i = 1:freqN
    for j = 1:numChan
        index = i:freqN:lengB;
        index = index';
        C{i} = [C{i}, CGb{j}(index,2)];
        G{i} = [G{i}, CGb{j}(index,4)];
        Z{i} = [Z{i}, ZDb{j}(index,2)];
        D{i} = [D{i}, ZDb{j}(index,4)];
    end
    C{i}(1:startPoint,:) = [];
    G{i}(1:startPoint,:) = [];
    Z{i}(1:startPoint,:) = [];
    D{i}(1:startPoint,:) = [];
    
    C{i}(end,:) = [];
    G{i}(end,:) = [];
    Z{i}(end,:) = [];
    D{i}(end,:) = [];
end
C = C';
D = D';
G = G';
Z = Z';

C_spectrum{numChan} = [];
G_spectrum{numChan} = [];
Z_spectrum{numChan} = [];
D_spectrum{numChan} = [];

for j = 1:numChan
    tmpC = [];
    tmpD = [];
    tmpG = [];
    tmpZ = [];
    for i = 1:freqN
        tmpC = [tmpC; C{i}(:,j)'];
        tmpG = [tmpG; G{i}(:,j)'];
        tmpZ = [tmpZ; Z{i}(:,j)'];
        tmpD = [tmpD; D{i}(:,j)'];
    end
    C_spectrum{j} = tmpC;
    G_spectrum{j} = tmpG;
    Z_spectrum{j} = tmpZ;
    D_spectrum{j} = tmpD;
end

% for anal_index = 32:52
%for anal_index = 35

C_ch = zeros(length(freq),length(time_index),numChan);
wG_ch = zeros(length(freq),length(time_index),numChan);
C_ch_cor = zeros(length(freq),length(time_index),numChan);
wG_ch_cor = zeros(length(freq),length(time_index),numChan);

parfor kk = 1:length(time_index)
    anal_index = time_index(kk);

folderName = [foldPath 'figure_DRT_again' num2str(anal_index) '/'];
if ~isfolder(folderName)
    mkdir(folderName);
end


% load([foldPath 'CNLS_dat/CNLS_delta_index' num2str(anal_index) '.mat'])
w = 2*pi*freq;
% parameters which optimize the lambda
lambda = -10:1/10:3;
lambda = 10.^lambda';

for j = 1:numChan
%     if (j==1)|(j==2)|(j==5)|(j==16)
%         continue
%     end

%%% do not skip the for loop if channel is involved in chan vector
loop_check = 1;
start_loop = 1;
while loop_check<(length(channel_index)+1)
    if channel_index(loop_check)==j
        start_loop = 0;
        break
    end
    loop_check = loop_check+1;
end

if start_loop
    continue
end

%     load([folderName 'Parameters_15_1_channel' num2str(j) '.mat'])
    S = load([folderName 'Parameters_gamma_open_15_1_channel' num2str(j) '.mat'],'re_im');
%     max_re_im = max(S.re_im(41:61));
%     max_re_im = find(S.re_im==max_re_im);
%     min_re_im = min(S.re_im(max_re_im:61));    % only find the optimzed lambda in the range [1e-7 1e-5]
%     min_re_im = find(S.re_im==min_re_im);
    min_re_im = min(S.re_im(45:48));    % only find the optimzed lambda in the range [1e-7 1e-5]
    min_re_im = find(S.re_im==min_re_im);
    % min_re_im = min_vector(j);
%     min_re_im = 47;
    min_re_im = 31;
%     if (j==1)||(j==9)
%         min_re_im = 47;
%     else
%         min_re_im = 46;
%     end3
%     min_re_im = 46;
%     min_re_im = 51;        %%%%%%%%%%t%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if anal_index<17
%     min_re_im = 46;
% else
%     min_re_im = 96;
% end
% % % if (j==9)||(j==12)
% % %     min_re_im = min(re_im(46:61));
% % % elseif j==16
% % %     min_re_im = min(re_im(51:61));
% % % else
% % %     min_re_im = min(re_im(43:61));
% % % end
% % % min_re_im = find(re_im==min_re_im);
%     if (j==1)||(j==2)
%         min_re_im  = 50;
%     elseif (j==4)||(j==12)||(j==15)
%         min_re_im = 48;
%     elseif j==16
%         min_re_im = 56;
%     else
%         min_re_im = 47;
%     end
% % % % % % if anal_index<17
% % % % % % %     min_re_im = 55;
% % % % % % % % min_re_im = min(re_im);
% % % % % % % % min_re_im = find(re_im==min_re_im);
% % % % % % % else
% % % % % % %     min_re_im = 105;
% % % % % % % end
% % % % % %     del = deltaCellOpen{j};
% % % % % %     del(26:end,:) = [];     % analyze the case that high boudary is less than 26
% % % % % % %     if j>4
% % % % % % %         del(13:end,:) = [];
% % % % % % %     end
% % % % % %     minDelta = min(min(del(5:25,5:20)));       % find the index that has minimum delta
% % % % % % %     minDelta = min(min(del));
    %%%%%%%%%%%%%%% for media concentration data %%%%%%%%%%%%%%%%%%%%%%%%%%
%     if j>2
%         minDelta = min(min(del(1:20,5:20)));
%     else
%         minDelta = min(min(del(5:18,1:8)));
%     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % %     [indr,indc] = find(del==minDelta);    
    
    if ~isempty(indLH)              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        indr = indLH(j,2);
        indc = indLH(j,1);
    end
    
% % % 
% % %     Z_re = Z_spectrum{j}(:,anal_index).*cos(D_spectrum{j}(:,anal_index)*pi/180);
% % %     Z_im = Z_spectrum{j}(:,anal_index).*sin(D_spectrum{j}(:,anal_index)*pi/180);
% % %     
% % %     C_ch(:,kk,j) = C_spectrum{j}(:,anal_index);
% % %     wG_ch(:,kk,j) = -G_spectrum{j}(:,anal_index)./w;
% % %     
% % %     fittedZ = warburgOpen(vestimatedCellOpen{j}(:,indr,indc),w);
% % %     Z_re_correct = Z_re - fittedZ(:,1) + vestimatedCellOpen{j}(1,indr,indc);
% % %     Z_im_correct = Z_im - fittedZ(:,2);
% % %     
% % %     C_ch_cor(:,kk,j) = imag(1./(w.*complex(Z_re_correct,Z_im_correct)));
% % %     wG_ch_cor(:,kk,j) = -real(1./(w.*complex(Z_re_correct,Z_im_correct)));
    %%%%%%%%%
    Z_re = Z_spectrum{j}(:,anal_index).*cos(D_spectrum{j}(:,anal_index)*pi/180);
    Z_im = Z_spectrum{j}(:,anal_index).*sin(D_spectrum{j}(:,anal_index)*pi/180);

    v_ini = [500; 1; .0001; 0.48];     % initialize the fitting variables

    result_cnls_open = CNLS(freq,Z_re,Z_im,indLH(j,1),indLH(j,2),v_ini,'open');
    
    fittedZ = warburgOpen(result_cnls_open{1},w);
    Z_re_correct = Z_re - fittedZ(:,1) + result_cnls_open{1}(1);
    Z_im_correct = Z_im - fittedZ(:,2);

    %%%%%%%%%%%
    posiStart = find(Z_im_correct>0);
    if isempty(posiStart)
        posiStart = 0;
    end
    posiStart = posiStart(end)+1;
    if posiStart<5
        posiStart = 5;
    end
     Z_re_diff = Z_re_correct(2:end)-Z_re_correct(1:end-1);
    posiStart2 = find(Z_re_diff>0);
    posiStart2(find(posiStart2>25)) = [];
    posiStart2 = posiStart2(end); 
    posiStart = max(posiStart,posiStart2);
%     posiStart = 8;      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % use the point that has negative reactance (capacitive behavior)
    Z_re_anal = Z_re_correct(posiStart:end);
    Z_im_anal = Z_im_correct(posiStart:end);
    freq2 = freq(posiStart:end);
    if freq2(1) < freq2(end)
        freq2 = fliplr(freq2')';
    end
    % parameters for DRT
    Nf = length(freq2);
    mean_diff_freq = mean(diff(log10(freq2)));
    %%%%%%%%%%%%%fitting frequency range extension%%%%%%%%%%%%%%%%%%%%
    freq_anal = log10(freq2(1))+2:mean_diff_freq/multi:log10(freq2(end))-2;
    freq_anal = 10.^freq_anal';
    Nfexpand = length(freq_anal);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tau_max = ceil(max(log10(1./freq_anal))) + .5;
    tau_min = floor(min(log10(1./freq_anal))) - .5;
    freq_fine = logspace(-tau_min,-tau_max,10*Nfexpand);
    epsilon = compute_epsilon(freq_anal,.5,rbf_type,shape_control);
    
    x_comb = zeros(Nfexpand+3);
    A_re = [];
    A_im = [];
    
        [~,~,~,~,~,~,~,x_comb,A_re, ...
            A_im,~] = hyperDRT(freq(posiStart:end),Z_re_anal, ...
            Z_im_anal,lambda(min_re_im),2,.5,shape_control,rbf_type,'Z',hyper_type);
        
        [gamma_comb,freq_comb] = map_array_to_gamma(freq_fine,freq_anal, ...
            x_comb(4:end),epsilon,rbf_type);
    dlmwrite([folderName '/DRT_gamma_open_15_1_channel' num2str(j) '.dat'],[1./freq_comb',gamma_comb])
end

end