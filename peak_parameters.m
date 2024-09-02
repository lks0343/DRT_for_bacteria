%% import data
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

startPoint = find(time<1/6);     % 痢≪??????????????????????? ?????????? 30?????? ?????????????????????????? analysis
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

%% find peaks & peak parameters
startPoint = [1; -5.8; .1; 28; -5.2; .2; 10; -4.1; .2];   % media : Muller Hinton
peak_position_range = [1e-8, 5e-7; 5e-7, 2.7e-6; 2.7e-6, 3e-5; 3e-5, 2e-4];
fitting_parameters = zeros(length(time_index),12,length(channel_index))*NaN;
fitting_parameters_gauss = zeros(length(time_index),9,length(channel_index))*NaN;
peak_parameter = zeros(length(time_index),8,length(channel_index));
large_peak = zeros(length(time_index),2,length(channel_index));
drtArea = zeros(length(time_index),length(channel_index));
asymmetry_index = zeros(length(time_index),4,length(channel_index));
norm_factor = zeros(length(time_index),1,length(channel_index));

peak_num = zeros(length(time_index),length(channel_index));
widLR = zeros(length(time_index),8,length(channel_index));

maxVal = zeros(length(time_index),1);

for k = 1:length(channel_index)
    j = channel_index(k);
    for kk = 1:length(time_index)
        anal_index = time_index(kk);
        folderName = [foldPath 'figure_DRT_again' num2str(anal_index) '/'];
        dat = dlmread([folderName 'DRT_gamma_open_15_1_channel' num2str(j) '.dat']);
        
        % to determine ylim of figure
        xlow = find(dat(:,1)<5e-7);     
        xhigh = find(dat(:,1)>1e-3);
        xlow = xlow(end);
        xhigh = xhigh(1);
        maxVal(kk) = max(maxVal(kk),max(dat(xlow:xhigh,2)));    
        
        % find the peaks of which position <1e-3
        start_ind = 1;
        finish_ind = find(log10(dat(:,1))>-3.3);
        finish_ind = finish_ind(1);
        
        inter = log10(dat(2,1))-log10(dat(1,1));
        norm_factorTMP = find(dat(:,1)>1e-8 & dat(:,1)<3e-3);
        norm_factor(kk,1,k) = sum(dat(norm_factorTMP,2));
        [pks,locs,w,p] = findpeaks(dat(start_ind:finish_ind,2),log10(dat(start_ind:finish_ind,1)), ...
            'WidthReference','halfheight','MinPeakWidth',15*inter);
        paek_num(kk,k) = length(pks);
        if length(pks)==4
            fitting_parameters(kk,:,k) = [pks(1), locs(1), w(1), ...
                pks(2), locs(2), w(2), ...
                pks(3), locs(3), w(3), ...
                pks(4), locs(4), w(4)];
            
            [asymmetry_index(kk,1,k),widLR(kk,1,k),widLR(kk,2,k)] = ...
                asymmetry(dat(:,2),log10(dat(:,1)),pks(1),locs(1),5);
            [asymmetry_index(kk,2,k),widLR(kk,3,k),widLR(kk,4,k)] = ...
                asymmetry(dat(:,2),log10(dat(:,1)),pks(2),locs(2),5);
            [asymmetry_index(kk,3,k),widLR(kk,5,k),widLR(kk,6,k)] = ...
                asymmetry(dat(:,2),log10(dat(:,1)),pks(3),locs(3),5);
            [asymmetry_index(kk,4,k),widLR(kk,7,k),widLR(kk,8,k)] = ...
                asymmetry(dat(:,2),log10(dat(:,1)),pks(4),locs(4),5);
            
        elseif length(pks)==3   % some pictures donot show smallest peak
            fitting_parameters(kk,:,k) = [pks(1), locs(1), w(1), ...
                NaN, NaN, NaN, ...
                pks(2), locs(2), w(2), ...
                pks(3), locs(3), w(3)];
            
            [asymmetry_index(kk,1,k),widLR(kk,1,k),widLR(kk,2,k)] = ...
                asymmetry(dat(:,2),log10(dat(:,1)),pks(1),locs(1),5);
            asymmetry_index(kk,2,k) = NaN;
            widLR(kk,3,k) = NaN;
            widLR(kk,4,k) = NaN;
            [asymmetry_index(kk,3,k),widLR(kk,5,k),widLR(kk,6,k)] = ...
                asymmetry(dat(:,2),log10(dat(:,1)),pks(2),locs(2),5);
            [asymmetry_index(kk,4,k),widLR(kk,7,k),widLR(kk,8,k)] = ...
                asymmetry(dat(:,2),log10(dat(:,1)),pks(3),locs(3),5);
            
        end
        
        % 1st peak
        start_ind = find(dat(:,1)<peak_position_range(1,1));
        start_ind = start_ind(end);
        tot_start = start_ind;
        finish_ind = find(dat(:,1)>peak_position_range(1,2));
        finish_ind = finish_ind(1);
        
        p0A = max(dat(start_ind:finish_ind,2));     % peak amplitude
        p0P = find(dat(:,2)==p0A);                  % peak position 
        p0P = dat(p0P,1);
        if length(p0P)>1
            p0P = p0P(1);
        end
        
        

        % 2nd peak
        start_ind = find(dat(:,1)<peak_position_range(2,1));
        start_ind = start_ind(end);
        finish_ind = find(dat(:,1)>peak_position_range(2,2));
        finish_ind = finish_ind(1);
        
        p1A = max(dat(start_ind:finish_ind,2));     % peak amplitude
        p1P = find(dat(:,2)==p1A);                  % peak position 
        p1P = dat(p1P,1);
        if length(p1P)>1
            p1P = p1P(1);
        end
        
        % 3rd peak
        start_ind = find(dat(:,1)<peak_position_range(3,1));
        start_ind = start_ind(end);
        finish_ind = find(dat(:,1)>peak_position_range(3,2));
        finish_ind = finish_ind(1);
        
        p2A = max(dat(start_ind:finish_ind,2));
        p2P = find(dat(:,2)==p2A);
        p2P = dat(p2P,1);
        if length(p2P)>1
            p2P = p2P(end);
        end
        
        % 4th peak
        start_ind = find(dat(:,1)<peak_position_range(4,1));
        start_ind = start_ind(end);
        finish_ind = find(dat(:,1)>peak_position_range(4,2));
        finish_ind = finish_ind(1);
        tot_finish = finish_ind;
        
        p3A = max(dat(start_ind:finish_ind,2));
        p3P = find(dat(start_ind:finish_ind,2)==p3A);
        p3P = dat(p3P,1);
        if length(p3P)>1
            p3P = p3P(end);
        end
        drtArea(kk,k) = sum(dat(1:tot_finish,2));
        peak_parameter(kk,:,k) = [p0A, p0P, p1A, p1P, p2A, p2P, p3A, p3P];
        
        p1P = find(dat(:,2)==max(dat(:,2)));
        p1P = dat(p1P,1);
        if length(p1P)>1
            p1P = p1P(1);
        end
        p1A = max(dat(:,2));
        large_peak(kk,:,k) = [p1A,p1P];

        

        start_ind = find(dat(:,1)<4e-7);
        start_ind = start_ind(end);
        finish_ind = find(dat(:,1)>5e-4);
        finish_ind = finish_ind(1);

        result_fitting = fit(log10(dat(start_ind:finish_ind,1)),dat(start_ind:finish_ind,2), ...
             'gauss3','Lower',[0 log10(4e-7) 0 0 log10(2.7e-6) 0 0 log10(3e-5) .0], ...
            'Upper', [100 log10(9e-5) 1/inter 100 log10(3e-5) 1/inter 100 log10(2e-4) 1/inter], ...
            'StartPoint',startPoint);
        fitting_result = [result_fitting.a1, 10^result_fitting.b1, result_fitting.c1, ...
            result_fitting.a2, 10^result_fitting.b2, result_fitting.c2, ...
            result_fitting.a3, 10^result_fitting.b3, result_fitting.c3];
        sorting_el = fitting_result([2,5,8]);
        sorting_max = max(sorting_el);
        bi = find(sorting_el==sorting_max);
        sorting_min = min(sorting_el);
        sm = find(sorting_el==sorting_min);
        mid = [1 2 3];
        mid(find(mid==bi)) = [];
        mid(find(mid==sm)) = [];


        fitting_parameters_gauss(kk,:,k) = [fitting_result((sm-1)*3+1), ...
            fitting_result((sm-1)*3+2), fitting_result((sm-1)*3+3), ...
            fitting_result((mid-1)*3+1), fitting_result((mid-1)*3+2), fitting_result((mid-1)*3+3), ...
            fitting_result((bi-1)*3+1), fitting_result((bi-1)*3+2), fitting_result((bi-1)*3+3)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  new fitting  %%%%%%%%%%%%%%%%%%%%%%%%%%%
standTmp = 10.^locs(find(locs<log10(5e-7)));
standTmp = find(dat(:,1)<standTmp);
standLow = standTmp(end);

if length(pks)==4
    standTmp =  10.^locs(2);
    standTmp = find(dat(:,1)<standTmp);
    standHigh = standTmp(end);
    standTmp = find(dat(standLow:standHigh,2)==min(dat(standLow:standHigh,2)));
    stand(1) = standTmp(1) + standLow;

    standTmp = 10.^locs(3);
    standTmp = find(dat(:,1)<standTmp);
    standLow = standTmp(end);
    standMid = min(dat(standHigh:standLow,2));  % local minima point between peak1, peak2
    standMid = find(dat(standHigh:standLow,2)==standMid)+standHigh;        %
    standMid = log10(dat(standMid(1),1));
    standTmp = 10.^locs(4);
    standTmp = find(dat(:,1)<standTmp);
    standHigh = standTmp(end);
    standTmp = find(dat(standLow:standHigh,2)==min(dat(standLow:standHigh,2)));
    stand(2) = standTmp(1) + standLow;
elseif length(pks)==3
    standTmp = 10.^locs(2);
    standTmp = find(dat(:,1)<standTmp);
    standHigh = standTmp(end);
    standTmp = find(dat(standLow:standHigh,2)==min(dat(standLow:standHigh,2)));
    stand(1) = standTmp(1)+standLow;

    standTmp = 10.^locs(2);
    standTmp = find(dat(:,1)<standTmp);
    standLow = standTmp(end);
    standTmp = 10.^locs(3);
    standTmp = find(dat(:,1)<standTmp);
    standHigh = standTmp(end);
    standTmp = find(dat(standLow:standHigh,2)==min(dat(standLow:standHigh,2)));
    stand(2) = standTmp(1)+standLow-1;
    standMid = locs(3)-.3;;
end

 lowLim = log10(dat(stand(1),1));
highLim = log10(dat(stand(2),1));
maxAmp = max(dat(stand(1):stand(2),2));
maxWidth = max(w);
maxRelax = log10(dat(stand(2),1));
result_fitting = fit(log10(dat(stand(1):stand(2),1)),dat(stand(1):stand(2),2),'gauss3', ...
    'Lower',[0 lowLim 0 0 lowLim 0 0 lowLim 0], ...
    'Upper',[maxAmp standMid maxWidth maxAmp maxRelax maxWidth maxAmp maxRelax maxWidth]);

fitting_result = [result_fitting.a1, result_fitting.b1, result_fitting.c1, ...
    result_fitting.a2, result_fitting.b2, result_fitting.c2, ...
    result_fitting.a3, result_fitting.b3, result_fitting.c3];
sorting_el = fitting_result([2,5,8]);
sorting_max = max(sorting_el);
bi = find(sorting_el==sorting_max);
sorting_min = min(sorting_el);
sm = find(sorting_el==sorting_min);
mid = [1 2 3];
mid(find(mid==bi))=[];
mid(find(mid==sm)) = [];

fitting_parameters_3gauss(kk,:,k) = [fitting_result((sm-1)*3+1), ...
            fitting_result((sm-1)*3+2), fitting_result((sm-1)*3+3), ...
            fitting_result((mid-1)*3+1), fitting_result((mid-1)*3+2), fitting_result((mid-1)*3+3), ...
            fitting_result((bi-1)*3+1), fitting_result((bi-1)*3+2), fitting_result((bi-1)*3+3)];


    end
end
        
parameters = fitting_parameters;
rr = isnan(parameters(:,4,:));
for i = 1:size(parameters,3)
    parameters(:,6,i) = fitting_parameters_3gauss(:,3,i)/0.575;
    parameters(:,5,i) = fitting_parameters_3gauss(:,2,i);
    parameters(:,4,i) = fitting_parameters_3gauss(:,1,i);
end

function [asy,leftP,rightP] = asymmetry(amp_data,pos_data,amp_peak,pos_peak,max_ratio)
fwhm_amp = amp_peak/max_ratio;
index_peak = find(pos_data<pos_peak);
index_peak = index_peak(end);
% find left 1/(max_ratio) of maximum position index
leftP = find(amp_data(1:index_peak)<fwhm_amp);
leftP = leftP(end);
x1 = pos_data(leftP);
x2 = pos_data(leftP+1);
y1 = amp_data(leftP);
y2 = amp_data(leftP+1);
leftP = (x2-x1)/(y2-y1)*(fwhm_amp-y1) + x1;

% find right 1/(max_ratio) of maximum position index
rightP = find(amp_data(index_peak+1:end)<fwhm_amp);
rightP = rightP(1)+index_peak-1;
x1 = pos_data(rightP);
x2 = pos_data(rightP-1);
y1 = amp_data(rightP);
y2 = amp_data(rightP-1);
rightP = (x2-x1)/(y2-y1)*(fwhm_amp-y1) + x1;

% calculate the asymmetry index
leftW = pos_peak - leftP;    % left width
rightW = rightP - pos_peak;  % right width

asy = log10(rightW/leftW);
end
