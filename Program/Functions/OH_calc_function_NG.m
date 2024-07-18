function [OH_ppm,results]=OH_calc_function_NG(intensity,varargin)
%% [OH_ppm,results]=OH_calc_function_NG
% calculates the OH content of FTIR Spectra for
% different minerals and Callibrations
% [OH_ppm,results]=OH_calc_function_NG(intensity,wavenumber,thickness,mineral,density,int_range,lincor,custom_I,sample,full_report)
%
%% %Input:
%   -intensity can be the integrated Absorbance [a], [a,c] or [a,b,c] or [n,n,n,n,...]  where a/b/c or it can be a matrix containing the spectra alpha,
%   beta and gamma (if the input is: one column - A=A=A / two columns -
%   A=A=C / tree columns -  A=B=C) if more than 3 input columns the average
%   of all is taken
%% %Optional inputs:
%   -'wavenumber' [a], [a,c] or [a,b,c],[n,n,n,n,...]  same as Intensity
%   however for n only the first one is considered and should be the same
%   for each n
%   -'thickness' [a], [a,c], [a,b,c]  [n,n,n,n,...] (thickness correction can be either a unique value or a
%   value for each n)
%   -'mineral' default: 'CPX' ['Cpx','Opx','Grt' 'Rutil']
%   -'density' [g/cm^3] default 3.3 or depending on mineral
%   -'int_range' [min max] default int_range=[3050 3800]
%   -'lincor' 0,1 linear background correction default lincor=1
%   -'custom_I' 
%   -'full_report' [0]-simple OH table;[1]- ABS OH table;[2]- extended table;[3]- full struc
%   -'sample' sets the sample name in the detailed report table
%% %Data generation for Examples
%     wavenumber=flipud(3000:3750);
%     pos=[3368 3462 3546 3634];
%     wid=[58 75 80 70];
%     amp=[  0.0119*0.20	0.0675*0.42	0.0377*0.25	0.0759*0.31;....
%     0.0119*0.27	0.0675*0.29	0.0377*0.24	0.0759*0.55; ...
%     0.0119*0.53	0.0675*0.29	0.0377*0.51	0.0759*0.14];
%     c=lines(3);
%     intensity=zeros(length(wavenumber),3);
%     for n=1:3
%     g=[];
%
%     for m=1:length(pos)
%     g(:,m)=gaussplot(pos(m),amp(n,m),wid(m),wavenumber);
%     end
%     intensity(:,n)=sum(g,2);
%     end
%     alpha = intensity(:,1);
%     beta = intensity(:,2);
%     gamma = intensity(:,3);
%% %Exmple 1
%     [OH_ppm,results]=OH_calc_function_NG(10,'mineral','cpx','density',3.2885,'thickness',150)
%
%% %Exmple 2
%     [OH_ppm,results]=OH_calc_function_NG([alpha beta gamma],'wavenumber',wavenumber,'thickness',[250 405 405],'int_range',[3300 3750])
%
%% %Exmple 3
%     [OH_ppm(:,1),OH_table]=OH_calc_function_NG([alpha beta gamma],'wavenumber',wavenumber,'density',3.3,'thickness',[250 405 405],'mineral','cpx','int_range',[3050 3750],'lincor',true)
%     [OH_ppm(:,2),OH_table]=OH_calc_function_NG(intensity,'wavenumber',wavenumber,'density',3.3,'thickness',[250 405 405],'mineral','cpx','int_range',[3300 3750])
%     [OH_ppm(1:4,3),OH_table]=OH_calc_function_NG(trapz(wavenumberwavenumber>3300 &wavenumber<3750,:),[alpha(wavenumber>3300 &wavenumber<3750) beta(wavenumber>3300 &wavenumber<3750) gamma(wavenumber>3300 &wavenumber<3750)]),'thickness',[250 405 405])
%% %Exmple 4
%     [OH_ppm,results]=OH_calc_function_NG([alpha beta gamma],'wavenumber',wavenumber,'density',3.3,'thickness',[250 405 405],'mineral','cpx','custom_i',[34000 50000 106000 145000 160000],'int_range',[3050 3750],'lincor',true)
%% %Exmple 5
%     [OH_ppm,results]=OH_calc_function_NG([alpha beta gamma],'wavenumber',wavenumber,'density',3.3,'thickness',[250 405 405],'mineral','cpx','int_range',[3050 3750],'lincor',true,'full_report',0,'sample','test_sample')
%     [OH_ppm,results]=OH_calc_function_NG([alpha beta gamma],'wavenumber',wavenumber,'density',3.3,'thickness',[250 405 405],'mineral','cpx','int_range',[3050 3750],'lincor',true,'full_report',1,'sample','test_sample')
%     [OH_ppm,results]=OH_calc_function_NG([alpha beta gamma],'wavenumber',wavenumber,'density',3.3,'thickness',[250 405 405],'mineral','cpx','int_range',[3050 3750],'lincor',true,'full_report',2,'sample','test_sample')
%     [OH_ppm,results]=OH_calc_function_NG([alpha beta gamma],'wavenumber',wavenumber,'density',3.3,'thickness',[250 405 405],'mineral','cpx','int_range',[3050 3750],'lincor',true,'full_report',3,'sample','test_sample')

% if length(varargin) > 18
%     error('ff_mat2tab:TooManyOptionalParameters', ...
%         'Too many input arguments. Allowed options are: mineral, density, thickness, lincor');
% end


mineral='Unknown';
sample='Unknown';
thickness=10000;
lincor=1;
M_calc=18.01528;
int_range=[3050 3800];
Cust_molar_absorption_coefficient=[];
full_report=0;
CPX_test={'cpx','clinopyroxenes','clinopyroxene','klinopyroxen','diopside','diopsid','augite','augit','omphacite','omph'};
OPX_test={'opx','orthopyroxenes','orthopyroxene','orthopyroxen','enstatite','enstatit','ferrosilite','ferrosilite'};
Garnet_test={'grt','garnet','garnets','almandine','almandin','grossular','hydrogrossular','pyrope','pyrop','spessartine','spessartin'};
Rutile_test= {'rutile','rutil'};
Quartz_test= {'qz','qtz','quartz','quarz'};
Olivine_test= {'ol','olivines','olivine','olivin'};
Ky_test={'kyanite','disthen','ky'};
        min_test=[Garnet_test OPX_test CPX_test Rutile_test Quartz_test Olivine_test Ky_test];

if (~isempty(varargin))
    for c=1:length(varargin)
        if ischar(varargin{c})

            switch lower(varargin{c})
                case {'density'}
                    density=varargin{c+1};
                case {'thickness'}
                    thickness=varargin{c+1};
                case {'mineral'}
                    mineral=lower(varargin{c+1});
                    mineral_in=varargin{c+1};
                case {'lincor'}
                    lincor=varargin{c+1};
                case {'int_range'}
                    int_range=varargin{c+1};
                case {'wavenumber'}
                    wavenumber=varargin{c+1};
                case {'custom_i'}
                    Cust_molar_absorption_coefficient=varargin{c+1}(:)';
                case {lower(sample)}
                    
                case {'sample'}
                    sample=char(varargin{c+1});
                case {'full_report'}
                    full_report(:,1)=varargin{c+1}(:);
                case {'m_calc'}
                    M_calc=varargin{c+1};
                case min_test
                    %if not in varargin density int_range updaten
%                 otherwise
%                     disp(['Invalid optional argument, ', ...
                %        varargin{c}]);
            end % switch
        end % if ischar
    end % for
end % if isempty
if ~exist("wavenumber","var")
    idx_range=ones(size(intensity));
    wavenumber=1:size(idx_range,1);
end

%ATH find fix
if length(intensity(:,1))<length(intensity(1,:)') && not(numel(thickness)==size(intensity,2))
    intensity=intensity';
end

if size(wavenumber,1)<size(wavenumber,2) && not(numel(thickness)==size(wavenumber,2))
    wavenumber=wavenumber';
end

if size(wavenumber,2)==1
    wn_alpha=wavenumber(:);
    wn_beta=wavenumber(:);
    wn_gamma=wavenumber(:);
elseif size(wavenumber,2)==2
    wn_alpha=wavenumber(:,1);
    wn_beta=wn_alpha;
    wn_gamma=wavenumber(:,2);
else
    wn_alpha=wavenumber(:,1);
    wn_beta=wavenumber(:,2);
    wn_gamma=wavenumber(:,3);
end

idx_range_alpha=wn_alpha>min(int_range) & wn_alpha<(max(int_range));
idx_range_beta=wn_beta>min(int_range) & wn_beta<(max(int_range));
idx_range_gamma=wn_gamma>min(int_range) & wn_gamma<(max(int_range));


switch mineral
    case {'omphacite','omph'}
        if ~exist("density","var") || density==0
            density=3.4;
        end
        mineral='cpx';
    case {'augite','augit'}
        if ~exist("density","var") || density==0
            density=3.35;
        end
        mineral='cpx';
    case {'cpx','clinopyroxenes','clinopyroxene','klinopyroxen','diopside','diopsid'}
        if ~exist("density","var") || density==0
            density=3.3;
        end
        mineral='cpx';
    case OPX_test
        if ~exist("density","var") || density==0
            density=3.32;
        end
        mineral='opx';
        table_row_names={'Bell et al. (1995) OPX 80600'};
        molar_absorption_coefficient= [Cust_molar_absorption_coefficient 80600]';

    case Olivine_test
        if ~exist("density","var") || density==0
            density=3.32;
        end
        mineral='ol';
        table_row_names={'Withers et al. (2012) Olivine 45200 ± 2340 ','Balan (2011) Olivine 43000',' Thomas et al. (2009) Olivine 47,000 ± 1000','Koch- Müller et al. (2006) Olivine 37500 + 5000','Bell et al. (2003) Olivine 28450 + 1830 [3650-3450 cm-1]'};
        molar_absorption_coefficient= [Cust_molar_absorption_coefficient 45200 43000 47000,37500,28450]';

    case Garnet_test
        if ~exist("density","var") || density==0
            switch mineral
                case {'spessartine','spessartin'}
                    density=4.1;
                case {'grossular','hydrogrossular'} 
                    density=3.6;
                case {'pyrope','pyrop'}
                    density=3.74;
                case {'almandine','almandin'}
                    density=4.1;
                otherwise
                    density=3.7;
            end
        end
        mineral='grt';
        % table_row_names={'Maldener et al. (2003) Grt-Sps 4880','Maldener et al. (2003) Grt 14400','Bell et al. (1995) Grt 6700'};
        % molar_absorption_coefficient= [Cust_molar_absorption_coefficient 4880 14400 6700]';
        table_row_names={'Maldener et al. (2003) Grt-Sps 4880','Maldener et al. (2003) Grt 14400','Bell et al. (1995) Grt 6700', 'Rossman and Aines (1991) Grossular 35745','Rossman and Aines (1991) Hydrogrossular 18955','Rossman (1988) Spessartine 35152'};
        molar_absorption_coefficient= [Cust_molar_absorption_coefficient 4880 14400 6700 35745 18955 35152]';

    case Rutile_test
        mineral='rutile';
        molar_absorption_coefficient= [Cust_molar_absorption_coefficient 38000,26160]';
        table_row_names={'Maldener et al. (2001) Rutile 38000+-4000 ',' Hammer and Beran (1991) Rutile y = 1/4-> 26160'};
        if ~exist("density","var") || density==0
            density=4.23;
        end
            if numel(intensity)==1
                intensity(1)=intensity(1);
                intensity(2)=intensity(1);
                intensity(3)=0;
            elseif numel(intensity)==2
                intensity(1)=intensity(1);
                intensity(2)=intensity(2);
                intensity(3)=0;
            elseif    numel(intensity)==3
                if intensity(3)~=0
                    intensity(1)=mean(intensity(1));
                    intensity(2)=mean(intensity(1));
                    intensity(3)=0;
                end
            elseif not(exist("wavenumber","var"))
                intensity=mean(intensity,"omitnan");
                intensity(2)=intensity(1);
                intensity(3)=0;
            end

        
    case Quartz_test
        mineral='quartz';
        molar_absorption_coefficient= [Cust_molar_absorption_coefficient,89000]';
        table_row_names={'Thomas et al. 2009 89000'};
         if ~exist("density","var") || density==0 
            density=2.65;
         end
    case Ky_test
        molar_absorption_coefficient= [Cust_molar_absorption_coefficient,32900]';
        table_row_names={'Bell et al. 2004 32900'};

        if ~exist("density","var") || density==0
            density=3.61;
        end
        mineral='kyanite';
         
    otherwise
        if ~exist("density","var") || density==0
            density=3.3;
        end
                molar_absorption_coefficient= [Cust_molar_absorption_coefficient]';
        table_row_names={};

        
end

switch mineral
    case 'cpx'
        molar_absorption_coefficient=[Cust_molar_absorption_coefficient 38300 65000 110000 160000 83400 46103 ]';
        table_row_names={'Bell et al. (1995) CPX 38300','Koch-Müller et al. (2007) CPX 65000 ','Koch-Müller (2010) Cpx from density 110000 ','Koch-Müller (2010) from molar volume 160000 ','Katayama et al. (2006) CPX 83400','Aubaud et al. (2009) CPX 46103'};
end


if numel(thickness)==1
    thickness_alpha=thickness(1);
    thickness_beta=thickness(1);
    thickness_gamma=thickness(1);
elseif numel(thickness)==2
    thickness_alpha=thickness(1);
    thickness_beta=thickness(1);
    thickness_gamma=thickness(2);
elseif numel(thickness)==3
    thickness_alpha=thickness(1);
    thickness_beta=thickness(2);
    thickness_gamma=thickness(3);
elseif numel(thickness)==size(intensity,2)
    intensity_cor= intensity./repmat((thickness).*10000,size(intensity,1),1);
    alpha=mean(intensity_cor(idx_range_alpha,:),2);
    beta=alpha;
    gamma=alpha;
    thickness_alpha=10000;
    thickness_beta=10000;
    thickness_gamma=10000;
end




cust_table_row_names=string(Cust_molar_absorption_coefficient );
if numel(cust_table_row_names)>=1
    table_row_names=[cust_table_row_names table_row_names];
end

if numel(intensity)==1
    alpha_wt_pct = round((M_calc.*intensity)./(thickness_alpha/10000.*(density./0.001).*molar_absorption_coefficient).*10^6./10000,6);
    intensity(2:3)=intensity(1);
     beta_wt_pct = alpha_wt_pct;
     gamma_wt_pct = alpha_wt_pct;
    
    [OH_ppm results]=get_results(full_report,table_row_names,lincor,wavenumber,intensity,molar_absorption_coefficient,int_range,Cust_molar_absorption_coefficient,thickness,sample,mineral_in,density,intensity(1),intensity(2),intensity(3),[],[],[],thickness_alpha,thickness_beta,thickness_gamma,alpha_wt_pct,beta_wt_pct,gamma_wt_pct);
    return
elseif numel(intensity)==2
    alpha_wt_pct = round((M_calc.*intensity(1))./(thickness_alpha/10000.*(density./0.001).*molar_absorption_coefficient).*10^6./10000,6);
    gamma_wt_pct = round((M_calc.*intensity(2))./(thickness_gamma/10000.*(density./0.001).*molar_absorption_coefficient).*10^6./10000,6);
    beta_wt_pct = alpha_wt_pct;

    intensity(3)=intensity(2);
    intensity(2)=intensity(1);

    [OH_ppm results]=get_results(full_report,table_row_names,lincor,wavenumber,intensity,molar_absorption_coefficient,int_range,Cust_molar_absorption_coefficient,thickness,sample,mineral_in,density,intensity(1),intensity(2),intensity(3),[],[],[],thickness_alpha,thickness_beta,thickness_gamma,alpha_wt_pct,beta_wt_pct,gamma_wt_pct);
    return

elseif numel(intensity)==3
    alpha_wt_pct = round((M_calc.*intensity(1))./(thickness_alpha/10000.*(density./0.001).*molar_absorption_coefficient).*10^6./10000,6);
    beta_wt_pct = round((M_calc.*intensity(2))./(thickness_beta/10000.*(density./0.001).*molar_absorption_coefficient).*10^6./10000,6);
    gamma_wt_pct = round((M_calc.*intensity(3))./(thickness_gamma/10000.*(density./0.001).*molar_absorption_coefficient).*10^6./10000,6);
  
    [OH_ppm results]=get_results(full_report,table_row_names,lincor,wavenumber,intensity,molar_absorption_coefficient,int_range,Cust_molar_absorption_coefficient,thickness,sample,mineral_in,density,intensity(1),intensity(2),intensity(3),[],[],[],thickness_alpha,thickness_beta,thickness_gamma,alpha_wt_pct,beta_wt_pct,gamma_wt_pct);
    return
end

if strcmp(mineral,'rutile') && exist('wavenumber','var')
    if numel(thickness)==1 || numel(thickness)==size(intensity,2) 
        thickness=thickness;
    else
        error('Error - Wrong rutile thickness format ')
    end
    thickness_alpha=10000;
    thickness_beta=10000;
    thickness_gamma=10000;

    if size(intensity,2)==1
        intensity(:,1)=intensity(:,1)./(thickness./10000);
        intensity(:,2)=intensity(:,1);
        intensity(:,3)=zeros(size(intensity(:,1)));
    elseif size(intensity,2)==2
        intensity(:,1)=intensity(:,1)./(thickness./10000);
        intensity(:,2)=intensity(:,2)./(thickness./10000);
        intensity(:,3)=zeros(size(intensity(:,1)));
    elseif size(intensity,2)==3
        if any(intensity(:,3)~=0,'all')
            intensity(:,1)=mean(intensity,2,"omitnan")./(thickness./10000);
            intensity(:,2)=intensity(:,1);
            intensity(:,3)=zeros(size(intensity(:,1)));
        end
    else
        intensity=mean(intensity,2,"omitnan")./(thickness./10000);
        intensity(:,2)=intensity(:,1);
        intensity(:,3)=zeros(size(intensity(:,1)));
    end

alpha=intensity(idx_range_alpha,1);
    beta=intensity(idx_range_beta,2);
    gamma=intensity(idx_range_gamma,3);


else

if size(intensity,2)==1
    alpha=intensity(idx_range_alpha);
    beta=alpha;
    gamma=alpha;
elseif size(intensity,2)==2
    alpha=intensity(idx_range_alpha,1);
    beta=alpha;
    gamma=intensity(idx_range_gamma,2);
elseif size(intensity,2)==3
    alpha=intensity(idx_range_alpha,1);
    beta=intensity(idx_range_beta,2);
    gamma=intensity(idx_range_gamma,3);
elseif  size(intensity,2)>3
    if not(isequal(numel(wn_alpha),numel(wn_beta),numel(wn_gamma)))
        error('Error - Wrong intensity format')
    end
    if numel(thickness)==1
        alpha=mean(intensity(idx_range_alpha,:),2);
        beta=alpha;
        gamma=alpha;
   thickness_alpha=thickness;
        thickness_beta=thickness;
        thickness_gamma=thickness;
    elseif  numel(thickness)==size(intensity,2)
        intensity_cor= intensity./repmat((thickness./10000),size(intensity,1),1);
        alpha=mean(intensity_cor(idx_range_alpha,:),2);
        beta=alpha;
        gamma=alpha;
        thickness_alpha=10000;
        thickness_beta=10000;
        thickness_gamma=10000;
    else
        error('Error - Wrong thickness format')
    end
else
    error('Error - Wrong intensity format')
end
end

%%
wn_alpha=wn_alpha(idx_range_alpha);
wn_beta=wn_beta(idx_range_beta);
wn_gamma=wn_gamma(idx_range_gamma);
% wn_beta
% wn_gamma
if wn_alpha(1)>wn_alpha(end)
    wn_alpha=flipud(wn_alpha);
    alpha=flipud(alpha);
end

if wn_beta(1)>wn_beta(end)
    wn_beta=flipud(wn_beta);
    beta=flipud(beta);
end

if wn_gamma(1)>wn_gamma(end)
    wn_gamma=flipud(wn_gamma);
    gamma=flipud(gamma);
end

%%

%Paterson
Xi=(18/(2*density*1000))*10^6; % density factor
E=1;%/3;%orientation factor for unpolarized 1/3 and 1 for polarized OKAY BECAUSE HERE A+B+C== 3xA


molar_absorption_coefficient_alpha=[repmat(molar_absorption_coefficient,(size(246.6*(3753 - wn_alpha)'))); 246.6*(3753 - wn_alpha)'; 382*(3782 - wn_alpha)'];
molar_absorption_coefficient_beta=[repmat(molar_absorption_coefficient,(size(246.6*(3753 - wn_beta)'))); 246.6*(3753 - wn_beta)'; 382*(3832 - wn_beta)'];
molar_absorption_coefficient_gamma=[repmat(molar_absorption_coefficient,(size(246.6*(3753 - wn_gamma)'))); 246.6*(3753 - wn_gamma)'; 382*(3832 - wn_gamma)'];

if lincor==false
alpha_wt_pct = trapz(wn_alpha,(M_calc.*alpha...
    ./((thickness_alpha/10000).*(density/0.001)'.*molar_absorption_coefficient_alpha')).*10^6)'/10000;
beta_wt_pct = trapz(wn_beta,(M_calc.*beta...
    ./((thickness_beta/10000).*(density/0.001)'.*molar_absorption_coefficient_beta')).*10^6)'/10000;
gamma_wt_pct = trapz(wn_gamma,(M_calc.*gamma...
    ./((thickness_gamma/10000).*(density/0.001)'.*molar_absorption_coefficient_gamma')).*10^6)'/10000;
else
alpha_wt_pct = trapz_BC(wn_alpha,(M_calc.*alpha...
    ./((thickness_alpha/10000).*(density/0.001)'.*molar_absorption_coefficient_alpha')).*10^6)'/10000;
beta_wt_pct = trapz_BC(wn_beta,(M_calc.*beta...
    ./((thickness_beta/10000).*(density/0.001)'.*molar_absorption_coefficient_beta')).*10^6)'/10000;
gamma_wt_pct = trapz_BC(wn_gamma,(M_calc.*gamma...
    ./((thickness_gamma/10000).*(density/0.001)'.*molar_absorption_coefficient_gamma')).*10^6)'/10000;
end

if lincor==false
alpha_wt_Paterson=(Xi/(150*E))*trapz(wn_alpha,(alpha./thickness_alpha/10000)./(3780-wavenumber(idx_range_alpha)))*10000;
beta_wt_Paterson=(Xi/(150*E))*trapz(wn_beta,(beta./thickness_beta/10000)./(3780-wavenumber(idx_range_beta)))*10000;
gamma_wt_Paterson=(Xi/(150*E))*trapz(wn_gamma,(gamma./thickness_gamma/10000)./(3780-wavenumber(idx_range_gamma)))*10000;
else
alpha_wt_Paterson=(Xi/(150*E))*trapz_BC(wn_alpha,(alpha./thickness_alpha/10000)./(3780-wavenumber(idx_range_alpha)))*10000;
beta_wt_Paterson=(Xi/(150*E))*trapz_BC(wn_beta,(beta./thickness_beta/10000)./(3780-wavenumber(idx_range_beta)))*10000;
gamma_wt_Paterson=(Xi/(150*E))*trapz_BC(wn_gamma,(gamma./thickness_gamma/10000)./(3780-wavenumber(idx_range_gamma)))*10000;
end

 
[OH_ppm results]=get_results(full_report,table_row_names,lincor,wavenumber,intensity,molar_absorption_coefficient,int_range,Cust_molar_absorption_coefficient,thickness,sample,mineral_in,density,alpha,beta,gamma,wn_alpha,wn_beta,wn_gamma,thickness_alpha,thickness_beta,thickness_gamma,alpha_wt_pct,beta_wt_pct,gamma_wt_pct,alpha_wt_Paterson,beta_wt_Paterson,gamma_wt_Paterson);
 
function [OH_ppm results]=get_results(full_report,table_row_names,lincor,wavenumber,intensity,molar_absorption_coefficient,int_range,Cust_molar_absorption_coefficient,thickness,sample,mineral_in,density,alpha,beta,gamma,wn_alpha,wn_beta,wn_gamma,thickness_alpha,thickness_beta,thickness_gamma,alpha_wt_pct,beta_wt_pct,gamma_wt_pct,alpha_wt_Paterson,beta_wt_Paterson,gamma_wt_Paterson)
if numel(alpha)==1
    OH_wt_pct=round((alpha_wt_pct+beta_wt_pct+gamma_wt_pct),6);
    OH_abc_ppm=round([alpha_wt_pct beta_wt_pct gamma_wt_pct]*10000,2);
    table_row_names=[table_row_names];
elseif exist("gamma_wt_Paterson","var")
    OH_wt_pct=round([(alpha_wt_pct+beta_wt_pct+gamma_wt_pct);(alpha_wt_Paterson+beta_wt_Paterson+gamma_wt_Paterson)],6);
    OH_abc_ppm=round([[alpha_wt_pct;alpha_wt_Paterson] [beta_wt_pct;beta_wt_Paterson] [gamma_wt_pct;gamma_wt_Paterson]]*10000,2);
    table_row_names=[table_row_names {'Libowitzky and Rossman (1997) 246.6·[3753 – ν]','Paterson (1982) (3780−ν)','Balan et al. 2008 382·(3832−ν)'}];

else
    OH_wt_pct=round((alpha_wt_pct+beta_wt_pct+gamma_wt_pct),6);
    OH_abc_ppm=round([alpha_wt_pct beta_wt_pct gamma_wt_pct]*10000,2);
    table_row_names=[table_row_names {'Libowitzky and Rossman (1997) 246.6·[3753 – ν]','Balan et al. (2008)^ 382·(3832−ν)'}];
end

OH_ppm=OH_wt_pct.*10000;

OH_results=[ OH_wt_pct OH_ppm];
if numel(alpha)>1
    if lincor==false
    abs_alpha=trapz(wn_alpha,alpha/thickness_alpha*10000);
    abs_beta=trapz(wn_beta,beta/thickness_beta*10000);
    abs_gamma=trapz(wn_gamma,gamma/thickness_gamma*10000);
    else
        abs_alpha=trapz_BC(wn_alpha,alpha/thickness_alpha*10000);
    abs_beta=trapz_BC(wn_beta,beta/thickness_beta*10000);
    abs_gamma=trapz_BC(wn_gamma,gamma/thickness_gamma*10000);
    end
else
    abs_alpha=alpha/thickness_alpha*10000;
    abs_beta=beta/thickness_beta*10000;
    abs_gamma=gamma/thickness_gamma*10000;
end
%%
results_0=array2table(OH_results,"RowNames",table_row_names,"VariableNames",{'OH [wt%]','OH [ppm]'});
%%
    OH_abc_ppm_cell=num2cell(OH_ppm)';
    abs_abc_sum=sum([abs_alpha abs_beta abs_gamma])';
    results_1=cell2table([{string(sample)},{string(mineral_in)},density,min(int_range),max(int_range),abs_alpha, abs_beta, abs_gamma,abs_abc_sum,OH_abc_ppm_cell],"VariableNames",[{'Sample','Mineral','Density','Min Wn','Max Wn','a','b','c','Total Int. Abs./cm'}, cellstr(table_row_names)]);
%%
abs_abc=num2cell(sum([abs_alpha abs_beta abs_gamma])');

    for n=1:size(abs_abc,2)
        abs_abc{n}= array2table([[abs_alpha(n,:) abs_beta(n,:) abs_gamma(n,:)], abs_abc{n}],'VariableNames',{'a','b','c','Tot. [Abs/cm]'});
    end

    OH_abc_ppm_cell=num2cell(OH_ppm)';
    for n=1:size(OH_abc_ppm_cell,2)
        OH_abc_ppm_cell{n}= array2table([OH_abc_ppm(n,:), OH_abc_ppm_cell{n}],'VariableNames',{'a','b','c','Tot. [ppm]'});
    end
    results_2=cell2table([{string(sample)},{string(mineral_in)},density,min(int_range),max(int_range),abs_alpha, abs_beta, abs_gamma,abs_abc,OH_abc_ppm_cell],"VariableNames",[{'Sample','Mineral','Density','Min Wn','Max Wn','a','b','c','Total Int. Abs./cm'}, cellstr(table_row_names)]);
%%
    abs_abc_sum=sum([abs_alpha abs_beta abs_gamma])';
    results.abs_abc_table=array2table([abs_alpha abs_beta abs_gamma abs_abc_sum],'VariableNames',{'a','b','c','Tot. [Abs/cm]'});
    results.OH_abc_ppm_table=array2table([OH_abc_ppm, OH_ppm],'VariableNames',{'a','b','c','Tot. [ppm]'},'RowNames',cellstr(table_row_names));

    results.data.OH_abc_ppm_mat=[OH_abc_ppm, OH_ppm];
    results.data.OH_abc_ppm_vars={'a','b','c','Tot. [ppm]'};
    results.data.abs_abc_ppm_mat=[abs_alpha abs_beta abs_gamma abs_abc_sum];
    results.data.abs_abc_ppm_vars={'a','b','c','Tot. [Abs/cm]'};
    results.data.thickness=thickness;
    results.data.sample=sample;
    results.data.density=density;
    results.data.intensity=intensity;
    results.data.wavenumber=wavenumber;

    %% added

    %     results.data.baseline_y=intensity(wavenumber>min(int_range) & wavenumber<(max(int_range)))-linCorrect(intensity(wavenumber>min(int_range) & wavenumber<(max(int_range)))); %Baseline
    %     results.data.baseline_x=wavenumber(wavenumber>min(int_range) & wavenumber<(max(int_range))); %Baseline

    if numel(wavenumber)>1
        results.data.baseline_y=intensity(any(wavenumber>min(int_range) & wavenumber<(max(int_range)),2),:)-linCorrect(intensity(any(wavenumber>min(int_range) & wavenumber<(max(int_range)),2))); %Baseline
        results.data.baseline_x=wavenumber(any(wavenumber>min(int_range) & wavenumber<(max(int_range)),2)); %Baseline
    else
        results.data.baseline_y=nan; %Baseline
        results.data.baseline_x=nan; %Baseline
    end

    %% added end



    results.data.alpha=alpha;
    results.data.beta=beta;
    results.data.gamma=gamma;
    results.data.mineral=mineral_in;
    results.data.lincor=lincor;
    results.data.int_range=int_range;
    results.data.Cust_molar_absorption_coefficient=Cust_molar_absorption_coefficient;
    results.data.molar_absorption_coefficient=molar_absorption_coefficient;
    results.data.table_row_names=table_row_names;
    
    results.results_table_0=results_0;
    results.results_table_1=results_1;
    results.results_table_2=results_2;

    if full_report==1
    results=results_1;
    elseif full_report==2
    results=results_2;
    elseif full_report==0
    results=results_0;
    end


