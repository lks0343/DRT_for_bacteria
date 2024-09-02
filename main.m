global time col freq
%% import the bacterial data along frequencies%% import the data file
dListCGb = dir([foldPath, '*gg*CG*']);
dListZDb = dir([foldPath '*gg*ZD*']);

folderName = [foldPath 'CNLS_dat/'];
if ~isfolder(folderName)
    mkdir(folderName);
end

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

startPoint = find(time<1/6);     % start analyzing after 10 min due to the noise of equipment
time(end-startPoint+1:end) = [];
time(end) = [];

freq = CGb{1}(1:freqN,1);   % measuring frequency value 

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



ind_vec = zeros(25,20)*NaN;     % matrix for the indexing

for kk = 1:length(time_index)
    anal_index = time_index(kk);
	WaitMessage = parfor_wait(25*20*numChan,'Waitbar',true);
	vestimatedMatrixShort = zeros(4,25,20,numChan)*NaN;
	vestimatedMatrixOpen = zeros(4,25,20,numChan)*NaN;
	deltaMatrixShort = zeros(25,20,numChan)*NaN;
	deltaMatrixOpen = zeros(25,20,numChan)*NaN;

parfor ind = 1:25*20*numChan
    [k,i,j] = ind2sub(size(ind_vec),ind);

	loop_check = 1;
	start_parfor = 1;
	while loop_check<(length(channel_index)+1)
	    if channel_index(loop_check)==j
	        start_parfor = 0;
	        break
	    end
	    loop_check = loop_check+1;
	end
	
	if start_parfor
	    WaitMessage.Send;
	    continue
	end



% for i = 1:20

%     for k = 1:30
        lowBound = i;
        highBound = k;

        if highBound-lowBound<4
            WaitMessage.Send;
            continue 
        end
        if lowBound<4
            WaitMessage.Send;
            continue
        end
        if highBound>25
            WaitMessage.Send;
            continue
        end

        v_ini = [500; 1; .0001; 0.48];     % initialize the fitting variables
        
        Z_re = Z_spectrum{j}.*cos(D_spectrum{j}*pi/180);
        Z_im = Z_spectrum{j}.*sin(D_spectrum{j}*pi/180);
        Z_re_anal = Z_re(:,anal_index);
        Z_im_anal = Z_im(:,anal_index);


        result_cnls_open = CNLS(freq,Z_re_anal,Z_im_anal,lowBound,highBound,v_ini,'open');
        vestimatedMatrixOpen(:,ind) = result_cnls_open{1};
        check0 = sum(result_cnls_open{1});
        if check0
            Z_re_subtract = Z_re_anal - result_cnls_open{2}(:,1) + result_cnls_open{1}(1);
            Z_im_subtract = Z_im_anal - result_cnls_open{2}(:,2);
            posiStart = find(Z_im_subtract>0);
            if isempty(posiStart)
                posiStart = 0;
            end
            posiStart = posiStart(end)+1;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Z_re_diff = Z_re_subtract(2:end)-Z_re_subtract(1:end-1);
            posiStart2 = find(Z_re_diff>0);
            posiStart2(find(posiStart2>25)) = [];
            posiStart2 = posiStart2(end); 
            posiStart = max(posiStart,posiStart2);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if posiStart>25
                WaitMessage.Send;
                continue
            end
            [~,~,~,~,~,~,~,x_comb,A_re,A_im,~] = ...
                hyperDRT(freq(posiStart:end),Z_re_subtract(posiStart:end),Z_im_subtract(posiStart:end), ...
                1e-15,2,0.5,'FWHM Coefficient','Gaussian','Z','Gaussian');
            Z_re_rec = A_re*x_comb;
            Z_im_rec = A_im*x_comb;

            W = 1./sqrt(Z_re_subtract.^2 + Z_im_subtract.^2);
            delta_re = (Z_re_subtract(posiStart:end)-Z_re_rec).*W(posiStart:end)/length(Z_re_rec);
            delta_re = delta_re'*delta_re;
            delta_im = (Z_im_subtract(posiStart:end)-Z_re_rec).*W(posiStart:end)/length(Z_im_rec);
            delta_im = delta_im'*delta_im;
            delta_reim = log10(delta_re + delta_im);
            deltaMatrixOpen(ind) = delta_reim;
        end

WaitMessage.Send;
end

close all
for j = 1:numChan
    vestimatedCellShort{j} = vestimatedMatrixShort(:,:,:,j);
    vestimatedCellOpen{j} = vestimatedMatrixOpen(:,:,:,j);
    deltaCellShort{j} = deltaMatrixShort(:,:,j);
    deltaCellOpen{j} = deltaMatrixOpen(:,:,j);
end

saveName = [folderName 'CNLS_delta_index' num2str(anal_index) '.mat'];
save(saveName,'vestimatedCellOpen','vestimatedCellShort','deltaCellOpen','deltaCellShort')
clear deltaCellOpen deltaCellShort vestimatedCellOpen vestimatedCellShort
 
end
