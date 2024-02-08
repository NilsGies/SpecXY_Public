function [y, baseline]=GetBaselineFun(x,y,signal_split_n,value1,value2,method,interp_method,ignore_co2)
baseline=zeros(size(y));
% if uique y ==0
for m=1:size(y,2)
    if ignore_co2==true 
        try
            co2_low=2245;
        co2_high=2450;

        [ ~, closest_2450 ] = min( abs(x-co2_high) );
        closest_2450=x(closest_2450);

        [ ~, closest_2245 ] = min( abs(x-co2_low) );
        closest_2245=x(closest_2245);

       
        [y(x<closest_2245,m),baseline(x<closest_2245,m)]=BaselineFun(0,x(x<closest_2245),y(x<closest_2245,m),signal_split_n,value1,value2,method,interp_method);
        [y(x>closest_2450,m),baseline(x>closest_2450,m)]=BaselineFun(0,x(x>closest_2450),y(x>closest_2450,m),signal_split_n,value1,value2,method,interp_method);
        baseline(x<=closest_2450 & x>=closest_2245,m) =y (x<=closest_2450 & x>=closest_2245,m);
        
        y(x<=closest_2450 & x>=closest_2245,m)=0;
        catch
        [y(:,m),baseline(:,m)]=BaselineFun(0,x,y(:,m),signal_split_n,value1,value2,method,interp_method);
        end

    else
        [y(:,m),baseline(:,m)]=BaselineFun(0,x,y(:,m),signal_split_n,value1,value2,method,interp_method);
    end
end