function [matrix,MapType,MapName]= CalcMatrixFromSpec(x_signal,y_signal,maptype2generate,region,mapsize,density,E,new_name)
    M_calc=18.01528;  % 
    thickness=10000; % 10000 micron

    if isempty(x_signal) || isempty(y_signal)
        matrix=nan(size(mapsize));
        MapType=0;
        MapName='Error in map generation';
        return
    end
switch maptype2generate
    case 'Integration'

        if new_name==true
            MapName=['IntMap_' num2str(min(round(region))) '_' num2str(max(round(region)))];
        end
        region=sort(region);
        integration=zeros(1,size(y_signal,2));
        for m= 1:size(y_signal,2)
            try
            integration(m)=trapz_BC(x_signal(x_signal>min(region)&x_signal<max(region)),y_signal(x_signal>min(region)&x_signal<max(region),m)');
            catch
            integration(m)=nan;
            end
        end

%         if ~any(mapsize(:) == 1)
            if (numel(integration)/mapsize(2))/round(numel(integration)/mapsize(2))==1
                matrix=reshape(integration,[],mapsize(2))';
            else
                matrix=[];
                MapType=[];
                MapName='Error in map generation - Check mapsize';
                return
            end

%         elseif mapsize(1) == 1
%             matrix=reshape(integration,[],1);
%         elseif mapsize(2) == 1
%             matrix=fliplr(reshape(integration,[],1));
% 
%         end
        MapType=1; %intrangemap

        disp('new integration matrix generated')

    case 'Absolute'
        abs_wn=region;
        [ ~, ix ] = min( abs(x_signal-abs_wn ) );
        % wn(ix)
        abs_val=y_signal(ix,:);
        matrix=reshape(abs_val,[],mapsize(2))';

        if new_name==true
            MapName=['AbsMap_' num2str(min(round(abs_wn)))];
        end
        MapType=10; %absolute val map
    case 'Ratio'
        abs_wn=region(1);
        [ ~, ix ] = min( abs(x_signal-abs_wn ) );
        abs_val1=y_signal(ix,:);
        abs_wn=region(2);
        [ ~, ix ] = min( abs(x_signal-abs_wn ) );
        abs_val2=y_signal(ix,:);

        abs_val=abs_val1./abs_val2;
        matrix=reshape(abs_val,[],mapsize(2))';
        if new_name==true
            MapName=['RatioMap_' num2str(round(region(1))) '_' num2str(max(round(region(2))))];
        end

        MapType=30; %ratio map
    case 'Ratio log/log'
        abs_wn=region(1);
        [ ~, ix ] = min( abs(x_signal-abs_wn ) );
        abs_val1=y_signal(ix,:);
        abs_val1(abs_val1<0)=0;

        abs_wn=region(2);
        [ ~, ix ] = min( abs(x_signal-abs_wn ) );
        abs_val2=y_signal(ix,:);
        abs_val2(abs_val2<0)=0;

        abs_val=log(abs_val1)./log(abs_val2);
        matrix=reshape(abs_val,[],mapsize(2))';
        if new_name==true
            MapName=['RatioMap_log_log_' num2str(round(region(1))) '_' num2str(max(round(region(2))))];
        end
        MapType=31; %ratio loglog map
    case 'Ratio Norm'
        abs_wn=region(1);
        [ ~, ix ] = min( abs(x_signal-abs_wn ) );
        abs_val1=y_signal(ix,:);
        % abs_val1(abs_val1<0)=0;

        abs_wn=region(2);
        [ ~, ix ] = min( abs(x_signal-abs_wn ) );
        abs_val2=y_signal(ix,:);
        %abs_val2(abs_val2<0)=0;
        %abs_val= clr([abs_val1' abs_val2']);
        %  abs_val=  clr([ones(size(abs_val1'))*1 abs_val1' abs_val2'])
        abs_val=log(normalize(abs_val1,'range'))./log(normalize(abs_val2,'range'));


        matrix=reshape(abs_val,[],mapsize(2))';
        if new_name==true
            MapName=['RatioMap_Norm_' num2str(round(region(1))) '_' num2str(max(round(region(2))))];
        end
        MapType=32; %ratio loglog map
    case 'Ratio NormLog'
        abs_wn=region(1);
        [ ~, ix ] = min( abs(x_signal-abs_wn ) );
        abs_val1=y_signal(ix,:);
        % abs_val1(abs_val1<0)=0;

        abs_wn=region(2);
        [ ~, ix ] = min( abs(x_signal-abs_wn ) );
        abs_val2=y_signal(ix,:);
        %abs_val2(abs_val2<0)=0;
        %abs_val= clr([abs_val1' abs_val2']);
        %  abs_val=  clr([ones(size(abs_val1'))*1 abs_val1' abs_val2'])
        abs_val=log(normalize(abs_val1,'range')./normalize(abs_val2,'range'));


        matrix=reshape(abs_val,[],mapsize(2))';
        if new_name==true
            MapName=['RatioMap_NormLog_' num2str(round(region(1))) '_' num2str(max(round(region(2))))];
        end
        MapType=33; %ratio normlog/normlog map
    case 'Entropy'
        B=1;
        MSE=zeros(1,size(y_signal,2));
        for k=1:size(y_signal,2)
            MSE(:,k)=immse(x_signal(x_signal>min(region)&x_signal<max(region)),y_signal(x_signal>min(region)&x_signal<max(region),k));
            %immse(in(:,k),[1:numel(in(:,k))]')
        end
        abs_val=20*log(((2^B)-1)./sqrt(MSE));
        %abs_val=sqrt(MSE);
        abs_val=1./abs_val;%sqrt(MSE);
        matrix=reshape(abs_val,[],mapsize(2))';

        matrix=(matrix-min(matrix(:)))./max((matrix(:)-min(matrix(:))))*100;
        %MapValue=(E(id)-min(E)./max((E-min(E)))*100;
%         a=0;
%         b=1;
%         matrix2=a+matrix*(b-a);
% isequal(matrix,matrix2)
        if new_name==true
            MapName=['EntropyMap_' num2str(round(region(1))) '_' num2str(max(round(region(2))))];
        end
        MapType=1000; %absolute val map
        
end

if isempty(matrix)
            MapType=[];
                    MapType=0;
            MapName='Error in map generation';
            return
end

end

