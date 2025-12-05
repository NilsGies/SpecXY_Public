function [gaussianParameters,GOF]= NG_fit(x,y,NumPeaks,startingGuesses)
format long g;
format compact;
warning off all
NumTrials=0;
TrialError=[];
iteration_max=2;
for iteration=1:iteration_max

if not(exist('startingGuesses','var'))|| isempty(startingGuesses)
    if x(1)>x(end)
        [amp,locs,wid,~] = findpeaks(flip(y),flip(x));
    else
        [amp,locs,wid,~] = findpeaks(y,x);
    end

     for m=1:length(locs)
         g(:,m)=gaussplot(locs(m),amp(m),wid(m),x);
     end

    if not(isequal(numel(locs),NumPeaks))
        n=(max(x)-min(x));
        peakpos=linspace(n/(NumPeaks+1),n-(n/(NumPeaks+1)),NumPeaks)+min(x);
        peakwid=repmat((n/(3.*NumPeaks)),size(peakpos));
        %peakwid=repmat((n/NumPeaks*0.5),size(peakpos));
        startingGuesses=[peakpos' peakwid'];
    %    'ew'
    else
        peakpos =locs;
        peakwid=wid;%repmat(23,size(wid))%repmat((n/NumPeaks*0.5),size(startingGuessesX));
        startingGuesses=[peakpos' peakwid'];
        'ps'
    end
end


%-------------------------------------------------------------------------------------------------------------------------------------------
% Perform an iterative fit using the FMINSEARCH function to optimize the height, width and center of the multiple Gaussians.
%options = optimset('TolX', 1e-5, 'MaxFunEvals', 10^5);  % Determines how close the model must fit the data
% First, set some options for fminsearch().
options.TolFun = 1e-4;
% options.Display = 'iter';
options.Display = 'off';
options.TolX = 1e-5;
options.MaxIter = 100000;
% options.LargeScale='on';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HEAVY LIFTING DONE RIGHT HERE:
% Run optimization
%[parameter, fval, flag, output] = fminsearch(@(lambda)(fitgauss(lambda, x, y, NumTrials)), startingGuesses, options);
NORM = startingGuesses; 
startingGuesses = startingGuesses./NORM;
%[parameter, fval, flag, output] = fminsearch(@fitgauss,startingGuesses, options, NORM, x, y, NumTrials);
[parameter, fval, flag, output] = fminsearch(@fitgauss,startingGuesses, options, NORM, x, y, NumTrials,TrialError);

parameter = parameter.*NORM;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------------------------------------------------------------------------------------------------------

pos=parameter(:,1);
wid=parameter(:,2);
A=gaussplot(pos,1,wid,x);
c = A' \ y';

gaussianParameters = [c parameter ];
% Now sort parameters in order of increasing mean
gaussianParameters = sortrows(gaussianParameters, 2);
pos=gaussianParameters(:,2);
wid=gaussianParameters(:,3);
amp=gaussianParameters(:,1);
y_calc=sum(gaussplot(pos,amp,wid,x));

SStot=sum((y-mean(y)).^2);
SSres=sum((y-y_calc).^2);
Rsquared=min(SSres./SStot);
MeanFitError=norm(y-y_calc)./(sqrt(length(x))*max(y));
GOF=[MeanFitError Rsquared];

if iteration<iteration_max
startingGuesses=[pos wid];
end


end
end
%=======================================================================================================================================================
function [theError, NumTrials,TrialError]= fitgauss(lambda, NORM, x, y, NumTrials,TrialError)
% Fitting function for multiple overlapping Gaussians, with statements
% added (lines 18 and 19) to slow the progress and plot each step along the
% way, for educational purposes.
% Author: T. C. O'Haver, 2006


lambda = lambda.*NORM;

pos=lambda(:,1);
wid=lambda(:,2);

%g=gaussplot(pos,1,wid,x)

A=gaussplot(pos,1,wid,x);


c = A' \ y';
z = A' * c;
theError = norm(z - y');


if any(c<0)
    theError = theError + 1000000*NumTrials;
    return
end

% Penalty so that heights don'x become negative.
if sum(c < 0) > 0 || any(max(A<0)) || any(c<0) || isempty(c)
    theError = theError + 1000000;
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
     theError = theError + 1000000;
     return
 end
 if  any(pos>max(x)| pos<min(x))
    % Penalty so that integration is more than spectrum
     theError = theError + 1000000;
     return
 end
 if  any(wid>abs(x(1)-x(end)))
    % Penalty so that integration is more than spectrum
     theError = theError + 1000000;
     return
 end

% figure(1)
% cla
% hold on
% plot(x,y,'LineWidth',2)
% A2=gaussplot(pos,c,wid,x);
% plot(x,A2)
% plot(x,sum(A2),':','LineWidth',2)

%plot(x,sum(gaussplot(pos,c,wid,x)))

NumTrials = NumTrials + 1;
	TrialError(NumTrials) = theError;
end % of fitgauss()

