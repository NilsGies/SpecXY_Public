function output = wn_spec_OH_calc_BC(x,y,calibration,density,E)
thickness=10000;

if not(size(y,1)==size(x,1))
    y=y';
end

output=zeros(1,size(y,2));
tic

if x(1)>x(end)
    x=flipud(x);
    y=flipud(y);
end
y=linCorrect(y);

toc
    M_calc=18.01528;

switch calibration

    case 'Libowitzky and Rossman (1997)'
        if not(exist('molar_absorption_coefficient','var'))
            molar_absorption_coefficient=[246.6.*E.*(3753 - x)']; %'Libowitzky and Rossman (1997) 246.6·[3753 – ν]',
        end
output= (M_calc.*trapz(x,y)'...
    ./((thickness/10000).*(density/0.001)'.*molar_absorption_coefficient)).*10^6;

  
    case 'Balan et al. (2008)'
        if not(exist('molar_absorption_coefficient','var'))
            molar_absorption_coefficient=[382*(3782 - x)']; % 'Balan et al. 2008 382·(3832−ν)'}];
        end

output= (M_calc.*trapz(x,y)'...
    ./((thickness/10000).*(density/0.001)'.*molar_absorption_coefficient)).*10^6;

    case 'Paterson (1982)'
%        %         output2=trapz(x,(((y./thickness)*10000))./(150*E.*(3780-x)));
   %     molar_absorption_coefficient=(150*E.*(3780-mean(x)));
        molar_absorption_coefficient=(150.*E.*(3780-x))';
output= (M_calc.*trapz(x,y)'...
    ./((thickness/10000).*(density/0.001)'.*molar_absorption_coefficient)).*10^6;

    case 'Paterson (1982) simplified'
molar_absorption_coefficient= (150*(1/3)*(3780-3400))';
output= (M_calc.*trapz(x,y)'...
    ./((thickness/10000).*(density/0.001)'.*molar_absorption_coefficient)).*10^6;
% 
    case 'Int. QuantMap [wt.% H2O]'
         output= (M_calc.*trapz(x,y)'...
            ./((thickness/10000).*(density/0.001)'.*E)).*10^6;
 %CH2O= (M_H2O.*Map./((thickness/10000).*(density/0.001)'.*E)).*10^6;
 %CH2O= (M_H2O.*Map./((density/0.001)'.*E)).*10^6;

     %CH2O water content [μg/g H2O]
     % M_H2O molecular weight of H2O [18.01528 g mol-1], 
     %density of the mineral [g/cm^3], 
     % E integral/absolut molar absorption coefficient [L mol(H2O)-1 cm-2 or L mol(H2O)-1 cm-1].
     % thickness= 10000 micron


    case 'Int. QuantMap [μg/g H2O]'
         output= (M_calc.*trapz(x,y)'...
            ./((thickness/10000).*(density/0.001)'.*E)).*10^6;
end
%output=output*10000;
end

%%

%Paterson
%E=1;%/3;%orientation factor for unpolarized 1/3 and 1 for polarized OKAY BECAUSE HERE A+B+C== 3xA






