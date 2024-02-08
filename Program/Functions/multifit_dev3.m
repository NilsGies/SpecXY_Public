function [gaussianParameters,GOF]= multifit_dev2(x,y_in,NumPeaks,peakshape,extra,Delta,startingGuesses)


NumTrials=1;
TrialError=100000000;

% options = optimset('Display','iter','PlotFcns',@optimplotfval);

options.TolFun = 1e-5;
 %options.Display = 'iter';
options.Display = 'off';
options.TolX = 1e-4;
options.MaxIter = 10000;
NORM = startingGuesses; 
startingGuesses = startingGuesses./NORM;
[parameter, fval, flag, output] = fminsearch(@multi_fit_NG,startingGuesses, options, NORM, x', y_in,NumPeaks,peakshape,extra,Delta,startingGuesses, NumTrials,TrialError);
parameter = parameter.*NORM;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 pos=parameter(:,1);
 wid=parameter(:,2);

for k=1:size(y_in,2)
A=DeconPlotFun(peakshape,x,pos',wid',1,extra,Delta);
%A=gaussplot(pos',1,wid',x);
c = A \ y_in(:,k);
z = A * c;
y_calc=z;
SStot=sum((y_in(:,k)-mean(y_in(:,k))).^2);
SSres=sum((y_in(:,k)-y_calc).^2);
Rsquared(k)=min(SSres./SStot);
MeanFitError(k)=norm(y_in(:,k)-y_calc)./(sqrt(length(x))*max(y_in(:,k)));
end
%plot(theError)

MeanFitError=sum(MeanFitError);
Rsquared=sum(Rsquared)/numel(Rsquared);

 gaussianParameters = [c parameter ];
 % Now sort parameters in order of increasing mean
 gaussianParameters = sortrows(gaussianParameters, 2);

GOF=[MeanFitError Rsquared];

end

function [theError, NumTrials,TrialError]=multi_fit_NG(lambda, NORM, x, y,NumPeaks,peakshape,extra,Delta,startingGuesses, NumTrials,TrialError)
lambda = lambda.*NORM;

pos=lambda(:,1);
wid=lambda(:,2);


theError=ones(1,size(y,2))*10000;
for k=1:size(y,2)

%A=DeconPlotFun(peakshape,x,pos,wid,amp,extra,Delta)
A=DeconPlotFun(peakshape,x,pos,wid,1,extra,Delta);
%A=gaussplot(pos,1,wid,x);
c = A' \ y(:,k);
z = A' * c;
theError(k) = norm(z - y(:,k));

% Penalty so that heights don'x become negative.
if any(max(A<0)) || any(c<0) || isempty(c)
    theError = 10000000*NumTrials;
    return
end

if x(1)>x(end)
    x_test=flip(x);
    y_test=flip(y);
else
x_test=(x);
    y_test=(y);
end
if  any(trapz(x_test,A'))>trapz(x_test,y_test)
    % Penalty so that integration is more than spectrum
    theError = 10000000*NumTrials;
     return
 end
 if  any(pos>max(x)| pos<min(x))
    % Penalty so that integration is more than spectrum
    theError = 10000000*NumTrials;
     return
 end
 if  any(wid>abs(x(1)-x(end)))
    % Penalty so that integration is more than spectrum
    theError = 10000000*NumTrials;
     return
 end
end


theError=sum(theError)/numel(theError);
NumTrials = NumTrials + 1;
TrialError=theError;
end


