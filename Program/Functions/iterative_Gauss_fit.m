function [gaussianParameters,GOF]= iterative_Gauss_fit(x,y,NumPeaks,plot_figures,startingGuesses)
format long g;
format compact;
global c TrialError NumTrials
warning off all
NumTrials=0;
fontSize = 20;

tFit=x;
if x(1)>x(end)
    [amp,locs,wid,~] = findpeaks(flip(y),flip(x));
else
    [amp,locs,wid,~] = findpeaks(y,x);
end

for m=1:length(locs)
    g(:,m)=gaussplot(locs(m),amp(m),wid(m),x);
end

if not(exist('startingGuesses','var'))
    if not(isequal(numel(locs),NumPeaks))
        startingGuesses=[];
        peakpos=[(max(x)-min(x))/(NumPeaks+1):(max(x)-min(x))/(NumPeaks+1):(max(x)-min(x))-((max(x)-min(x))/(NumPeaks+1))]+min(x);
        for peak=1:NumPeaks
            markx=peakpos(peak);
            startingGuesses=[startingGuesses markx (max(x)-min(x))/ (3.*NumPeaks)];
        end
    else
        startingGuessesX =locs;% (n/(NumPeaks+1):n/(NumPeaks+1):n-(n/(NumPeaks+1)))+min(x);
        startingGuessesY=wid;%repmat(23,size(wid))%repmat((n/NumPeaks*0.5),size(startingGuessesX));
        startingGuesses=[startingGuessesX;startingGuessesY];
        startingGuesses=startingGuesses(:)';

    end
    startingGuesses=startingGuesses;
end
%-------------------------------------------------------------------------------------------------------------------------------------------
% Perform an iterative fit using the FMINSEARCH function to optimize the height, width and center of the multiple Gaussians.
options = optimset('TolX', 1e-4, 'MaxFunEvals', 10^12);  % Determines how close the model must fit the data
% First, set some options for fminsearch().
% options.TolFun = 1e-4;
% options.TolX = 1e-4;
% options.MaxIter = 10000;
options.TolFun = 1e-4;
% options.Display = 'iter';
options.TolX = 1e-4;
options.MaxIter = 100000;
% options.LargeScale='on'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HEAVY LIFTING DONE RIGHT HERE:
% Run optimization
%[parameter, fval, flag, output] = fminsearch(@(lambda)(fitgauss(lambda, tFit, y, NumTrials)), startingGuesses, options);
NORM = startingGuesses;
startingGuesses = startingGuesses./NORM;
%[parameter, fval, flag, output] = fminsearch(@fitgauss,startingGuesses, options, NORM, tFit, y, NumTrials);
[parameter, fval, flag, output] = fminsearch(@fitgauss,startingGuesses, options, NORM, tFit, y, NumTrials);

parameter = parameter.*NORM;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------------------------------------------------------------------------------------------------------
% Now plot results.
yhat = PlotComponentCurves(x, y, tFit, c, parameter,plot_figures);
% Compute the residuals between the actual y and the estimated y and put that into the graph's title.
estimatedMuSigma = reshape(parameter, 2, [])';
gaussianParameters = [c, estimatedMuSigma];
% Now sort parameters in order of increasing mean
gaussianParameters = sortrows(gaussianParameters, 2);

if plot_figures==true
    meanResidual = mean(abs(y - yhat));
    numGaussians = length(c);

    fprintf('The mean of the absolute value of the residuals is %f.\n', meanResidual);
    caption = sprintf('Estimation of %d Gaussian Curves that will fit data.  Mean Residual = %f.', numGaussians, meanResidual);
    title(caption, 'FontSize', fontSize, 'Interpreter', 'none');
    drawnow;
    % Make table for the fitted, estimated results.
    % First make numGaussians row by 3 column matrix: Column 1 = amplitude, column 2 = mean, column 3 = width.
    % 	parameter % Print to command window.

    % Display actual table in the command window.
    % Create table of the output parameters and display it below the actual, true parameters.
    tEstimate = table((1:numGaussians)', c(:), estimatedMuSigma(:, 1), estimatedMuSigma(:, 2), 'VariableNames', {'Number', 'Amplitude', 'Mean', 'Width'});
    % Plot the error as a function of trial number.
    hFigError = figure();
    hFigError.Name = 'Errors';
    plot(TrialError, 'b-');
    % hFigError.WindowState = 'maximized';
    grid on;
    xlabel('Trial Number', 'FontSize', fontSize)
    ylabel('Error', 'FontSize', fontSize)
    set(gca, 'YScale', 'log')
    caption = sprintf('Errors for all %d trials.', length(TrialError));
    title(caption, 'FontSize', fontSize, 'Interpreter', 'none');
    % message = sprintf('Done!\nHere is the result!\nNote: there could be multiple ways\n(multiple sets of Gaussians)\nthat you could achieve the same sum (same test curve).');
    % fprintf('Done running %s.m.\n', mfilename);
    % msgbox(message);
end

SStot=sum((y-mean(y)).^2);
SSres=sum((y-yhat).^2);
Rsquared=min(SSres./SStot);

MeanFitError=norm(y-yhat)./(sqrt(length(x))*max(y));

GOF=[MeanFitError Rsquared];
end
%=======================================================================================================================================================
function yhat = PlotComponentCurves(x, y, t, c, parameter,plot_figures)
fontSize = 20;
% Get the means and widths.
means = parameter(1 : 2 : end);
widths = parameter(2 : 2 : end);
% Now plot results.
if plot_figures==true

    hFig2 = figure;
    hFig2.Name = 'Fitted Component Curves';
    % 	plot(x, y, '--', 'LineWidth', 2)
    hold on;
    set(gca,'XDir','reverse')

end

yhat = zeros(1, length(t));
numGaussians = length(c);
legendStrings = cell(numGaussians + 2, 1);
for k = 1 : numGaussians
    % Get each component curve.
    thisEstimatedCurve = c(k) .* gaussian(t, means(k), widths(k));
    % Plot component curves.
    if plot_figures==true

        plot(x, thisEstimatedCurve, '-', 'LineWidth', 2);
        hold on;
    end
    % Overall curve estimate is the sum of the component curves.
    yhat = yhat + thisEstimatedCurve;
    legendStrings{k} = sprintf('Estimated Gaussian %d', k);
end
if plot_figures==true

    % Plot original summation curve, that is the actual curve.
    plot(x, y, 'r-', 'LineWidth', 1)
    % Plot estimated summation curve, that is the estimate of the curve.
    plot(x, yhat, 'k--', 'LineWidth', 2)
    grid on;
    xlabel('X', 'FontSize', fontSize)
    ylabel('Y', 'FontSize', fontSize)
    set(gca,'XDir','reverse')

    caption = sprintf('Estimation of %d Gaussian Curves that will fit data.', numGaussians);
    title(caption, 'FontSize', fontSize, 'Interpreter', 'none');
    grid on
    legendStrings{numGaussians+1} = sprintf('Actual original signal');
    legendStrings{numGaussians+2} = sprintf('Sum of all %d Gaussians', numGaussians);
    legend(legendStrings);
    xlim(sort([x(1) x(end)]));
    hFig2.WindowState = 'maximized';
    drawnow;
end

end % of PlotComponentCurves
%=======================================================================================================================================================
function [theError, NumTrials]= fitgauss(lambda, NORM, t, y, NumTrials)
% Fitting function for multiple overlapping Gaussians, with statements
% added (lines 18 and 19) to slow the progress and plot each step along the
% way, for educational purposes.
% Author: T. C. O'Haver, 2006
global c NumTrials TrialError


lambda = lambda.*NORM;

A = zeros(length(t), round(length(lambda) / 2));
for j = 1 : length(lambda) / 2
    A(:,j) = gaussian(t, lambda(2 * j - 1), lambda(2 * j))';



    if max(A(:,j))<0 || lambda(2 * j - 1)>max(t) || lambda(2 * j - 1) < min(t)
        theError=1e15;
        return
    end

end
c = A \ y';
z = A * c;
theError = norm(z - y');

% Penalty so that heights don't become negative.
if sum(c < 0) > 0
    theError = theError + 1000000;
end

NumTrials = NumTrials + 1;
TrialError(NumTrials) = theError;
end % of fitgauss()
%=======================================================================================================================================================
function g = gaussian(x, peakPosition, width)
%  gaussian(x,pos,wid) = gaussian peak centered on pos, half-width=wid
%  x may be scalar, vector, or matrix, pos and wid both scalar
%  T. C. O'Haver, 1988
% Examples: gaussian([0 1 2],1,2) gives result [0.5000    1.0000    0.5000]
% plot(gaussian([1:100],50,20)) displays gaussian band centered at 50 with width 20.
g = exp(-((x - peakPosition) ./ (0.60056120439323 .* width)) .^ 2);
end % of gaussian()

% function g=gaussplot(pos,amp,wid,x)
% g = amp.*(exp(-((x-pos)./(0.60056120439323.*wid)).^2));
% end

%=======================================================================================================================================================