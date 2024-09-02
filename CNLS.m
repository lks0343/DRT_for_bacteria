% the function was made by KS Lee
function yout = CNLS(freq,Z_re,Z_im,lowBound,highBound,v_ini,warburg)
fitData = [Z_re(lowBound:highBound), Z_im(lowBound:highBound)];

w = freq(lowBound:highBound)*2*pi;


options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt', ...
    'FunctionTolerance',1e-15,'StepTolerance',1e-15,'ScaleProblem','jacobian', ...
    'MaxFunctionEvaluations',100000,'MaxIterations',5000);

if strcmp(warburg,'short')
    fittingFunc = @(v,w) warburgShort(v,w);
elseif strcmp(warburg,'open')
    fittingFunc = @(v,w) warburgOpen(v,w);
end

% complex non-linear square fitting 
[vestimated,resnorm,residuals,exitflag,output] = ...        
    lsqcurvefit(fittingFunc,v_ini,w,fitData,[0;0;0;0],[Inf,Inf,Inf,1],options);


fittedZ = fittingFunc(vestimated,freq*2*pi);    % fitting result
correctedZ_re = Z_re(:) - fittedZ(:,1) + vestimated(1);     % (raw data) - (fitting data)
correctedZ_im = Z_im(:) - fittedZ(:,2);

% fittedZ2 = fittingFunc(vestimated,w);
% 
% % find the apex index 
% find_apex = Z_im(1:end-1) - Z_im(2:end);
% find_apex = find_apex>0;
% find_apex = find_apex(1:end-1) - find_apex(2:end);
% find_apex = find(find_apex==-1);
% find_apex = find_apex + 1;
% % check that the fitted data around apex has larger Z_re
% check_large = fittedZ(:,2) - Z_im(find_apex);
% k = 1;
% while k<length(fittedZ)
%     if (check_large(k)<0)&(check_large(k+1)>0)
%         break
%     end
%     k = k+1;
% end
% 
% %  if fitted data around apex is existed on left, then we don't use the results
% if Z_re(find_apex)>fittedZ(k,1)
%     vestimated = [0; 0; 0; 0];
%     fittedZ = [];
% end


yout{1} = vestimated;
yout{2} = fittedZ;

% output
% 
% figure
% scatter(Z_re(:),Z_im(:))    
% hold on
% plot(fittedZ(:,1),fittedZ(:,2))
% scatter(correctedZ_re,correctedZ_im)
% 
% grid on
% text(350,-5000,num2str(vestimated))
% 
% title(['lowBound : ' num2str(lowBound) ',     highBound : ' num2str(highBound)])
end




