function [Base, Corrected_Spectrum]=baseline_raman(Spectrum,tlimit)
% Raman Spectrum Baseline Removal (https://www.mathworks.com/matlabcentral/fileexchange/69649-raman-spectrum-baseline-removal), MATLAB Central File Exchange. Retrieved December 20, 2023. 

%Input
%-------
%Spectrum: vector of size (N*1)
%Output
%-------
%Base: Identified Baseline vector of size (N*1)
%Corrected_Spectrum: Corrected Spectrum vector of size (N*1)

l=length(Spectrum);  
lp=ceil(0.5*l);
initial_Spectrum=[ones(lp,1)*Spectrum(1) ; Spectrum ; ones(lp,1)*Spectrum(l)];
l2=length(initial_Spectrum);
S=initial_Spectrum;
n=1;
flag1=0;
if exist('tlimit',"var")
    time0 = tic;
end
while flag1==0
    if exist('tlimit',"var") && toc(time0)>tlimit
        error('too long')
    end
    n=n+2;

    i=(n-1)/2;

    [Baseline, stripping]=peak_stripping(S,n);
    A(i)=trapz(S-Baseline);
    Stripped_Spectrum{i}=Baseline;
S=Baseline;
if i>3
    if A(i-1)<A(i-2) && A(i-1)<A(i)
        i_min=i-1;
        flag1=1;
    end 
end
end
Base=Stripped_Spectrum{i_min}; 
Corrected_Spectrum=initial_Spectrum-Base; Corrected_Spectrum=Corrected_Spectrum(lp+1:lp+l);
Base=Base(lp+1:lp+l);
end

function [Baseline, stripping]=peak_stripping(Spectrum,Window)
stripping=0;
y=sgolayfilt(Spectrum,0,Window);
n=length(Spectrum);
Baseline=zeros(n,1);
for i=1:1:n
   if Spectrum(i)>y(i)
       stripping=1;
       Baseline(i)=y(i);
   else
       Baseline(i)=Spectrum(i);
   end
       
    
    
end
end