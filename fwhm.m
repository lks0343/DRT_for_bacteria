function width = fwhm(data,posinega,ratio)

switch posinega
    case '+'
        posinega = 1;
    case '-'
        posinega = -1;
    otherwise
        error('posinega should have only + or -')
end

dat = data*posinega;
threshold = max(dat);
peakIndex = find(dat==threshold);
peakIndex = peakIndex(1);

threshold = threshold*ratio;
cross = dat - threshold;
cross = cross(1:end-1) .* cross(2:end);         % points that cross the threshold will have negative sign

crossPoint = find(cross<=0);                    % = exact same data point with threshold value
upper = min(crossPoint(crossPoint>peakIndex));     % upper index of crossing point
lower = max(crossPoint(crossPoint<peakIndex));     % lower index of crossing point


if isempty(upper)
    upperV = length(dat);
else
    upperV = 1/(dat(upper+1)-dat(upper))*(threshold - dat(upper)) + upper;     % linear interpolation 
end

if isempty(lower)
    lowerV = 1;
else
    lowerV = 1/(dat(lower+1)-dat(lower))*(threshold - dat(lower)) + lower;     % linear interpolation 
end

width = upperV - lowerV;

end



        


