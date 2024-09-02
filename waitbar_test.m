
startTime = datetime;
startTime.Format = 'HH:mm:ss';
startTime = string(startTime);
gg = waitbar(0,['Start at : ' startTime '  progress : 0 %'])
afterEach(q,gg)
for i = 1:100
    for j = 1:100
        for k = 1:100
            send(p,
        end
        waitbar((i+j/100)/100,gg,['Start at : ' startTime  num2str(i+j/100) ' %'])
    end
end
