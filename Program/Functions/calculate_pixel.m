function[peakmat_hight,peakmat_area,peakmat_width,peakerrormat1_area,peakerrormat2_area,peakstart,shift] =  calculate_pixel(I,starttest,mask,mask_sort,id,x_signal,y_signal,center,window,NumPeaks,peakshape,extra,NumTrials,start,BaselineMode,fixedparameters,plots,bipolar,minwidth,DELTA,clipheight,fixedparameters_shift)
            shift=zeros(NumPeaks,1);

if numel(mask)==1
    Calculate_ID=true;
elseif numel(mask)>1 && not(mask_sort(id)==0)% when mask existst and is active
    Calculate_ID=true;
else
    Calculate_ID=false;
end

if Calculate_ID==true
    signal=[x_signal y_signal(:,id)];
    multi_start=true;
    if multi_start==true && starttest  %ATH Add multiple start solutions depending on ROIs / UIStarts / class  or Last solution
        [FitResults,GOF]=get_best_start_guess(signal,center,window,NumPeaks,peakshape,extra,NumTrials,start,BaselineMode,fixedparameters,plots,bipolar,minwidth,DELTA,clipheight,fixedparameters_shift);
    else
        try
            if (isempty(start)||start==0) && peakshape==1
                [FitResults,~,~,~,~,~,~]=peakfit(signal,center,round(window),NumPeaks,16,extra,0,start,BaselineMode,fixedparameters,plots,bipolar,minwidth,DELTA);
                s=[FitResults(:,2)';FitResults(:,4)'];
                start=s(:);
            elseif isempty(start) && peakshape==2
                [FitResults,~,~,~,~,~,~]=peakfit(signal,center,round(window),NumPeaks,17,extra,0,start,BaselineMode,fixedparameters,plots,bipolar,minwidth,DELTA);
                s=[FitResults(:,2)';FitResults(:,4)'];
                start=s(:);
            end

            [FitResults,GOF,~,~,~,~,~]=peakfit(signal,center,window,NumPeaks,peakshape,extra,NumTrials,start,BaselineMode,fixedparameters,plots,bipolar,minwidth,DELTA,clipheight,fixedparameters_shift);
            
            if sum(FitResults(:,3))==0 && not(sum(signal(:,2))==0)
                disp('no_result')
                [FitResults,GOF,~,~,~,~,~]=peakfit(signal,center,window,NumPeaks,peakshape,extra,NumTrials,start,BaselineMode,fixedparameters,plots,bipolar,minwidth,0,clipheight);

          
            end
        catch ME
            disp(ME.message)
            FitResults=repmat([nan nan nan nan nan],NumPeaks,1); %NANANANANAANNAANNAAN BATMAAAAAN
            GOF=[nan nan];
        end
    end



        [shift,FitResults]=sort_decon_result( FitResults,fixedparameters); %start
       % FitResults=sort_decon_result( FitResults,start);
   
else
    FitResults=repmat([nan nan nan nan nan],NumPeaks,1); %NANANANANAANNAANNAAN BATMAAAAAN
    GOF=[nan nan];
end


% try
%     
% 
%     [~,idx] = sort(FitResults(:,2)); % sort just the first column
%     FitResults = FitResults(idx,:); % sort the matrix
%     
% catch ME
%     disp(ME.message)
% end




peakmat_hight=FitResults(:,3)';
peakmat_area=FitResults(:,5)';
peakmat_width=FitResults(:,4)';

peakerrormat1_area=GOF(:,1)';
peakerrormat2_area=GOF(:,2)';
peakstart=I;

    function [shift,FitResults_best1]=sort_decon_result( FitResults_best1,pos_ini)
        pos=FitResults_best1(:,2)';
        pos_ini=pos_ini';
        pos_ini2=(pos_ini);
        pos_t2=(pos);
        idxxx=1:size(pos_ini,2);
        idyyy=1:size(pos_ini,2);
        sort1=zeros(1,size(pos_ini,2));
        sort2=zeros(1,size(pos_ini,2));
        for n=1:size(pos_ini,2)

            ix_pos = abs(pos_ini2(:) - pos_t2(:)');
            [xy, col2delete] = min(ix_pos,[],1);

            [~, row2delete ]=min(xy,[],2);

            pos_ini2(col2delete(row2delete))=[];

            sort1(n)=idxxx(col2delete(row2delete));
            idxxx(col2delete(row2delete))=[];
            try
                sort2(n)=idyyy((row2delete));
            catch ME
                disp(ME.message)
            end
            idyyy((row2delete))=[];

            pos_t2(row2delete)=[];
        end


        sorted_ini=pos_ini(sort1);
        sorted_pos_r=pos(sort2);
        FRsort=FitResults_best1(sort2,:);

        [~,idx] = sort(pos_ini);
        [~,idx(idx)] = sort(sorted_ini);
        FitResults_best1=FRsort(idx,:);
        pos=FitResults_best1(:,2);
        shift= pos'-pos_ini;

    end
end
