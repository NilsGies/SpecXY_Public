function [FitResults,GOF]=get_best_start_guess(signal,center,window,NumPeaks,peakshape,extra,NumTrials,start,autozero,fixedparameters,plots,bipolar,minwidth)
rmsq=zeros(size(start,2),2);
parfor i=1:size(start,2)
    [ FR{i},GF{i},~,~,~,~,~]=peakfit(signal,center,window,NumPeaks,peakshape,extra,NumTrials,start(:,i),autozero,fixedparameters,plots,bipolar,minwidth);
    rmsq(i,:)=GF{i}
end
[~,I]=max(rmsq(:,2));
GOF=GF{I};
FitResults=FR{I};
end