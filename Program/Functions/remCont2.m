% Continuum removal for spectra in a hyperspectral image
% The input parameter HShc is a hypercube object representing the hyperspectral image.
% If the Hyperspectral Image Processing Toolbox is not available, input can
% be changed to consist of the image HS_im (numeric 3D array of reflectances)
% and wl (numeric array of wavelengths) objects.
% Scaling of the output is arbitrary, but uniform throughout the image.
% The author takes no responsibility of possibly existing bugs in the code.
% Author: Johanna Torppa
% Date: September 9, 2022

function HS_scale=remCont2(wl,HS_im)
%     HS_im=double(HShc.DataCube);
%     wl=HShc.Wavelength;
    sz=size(HS_im);
    nwl=size(wl,1);
    HS_im(HS_im==0)=NaN;
    HS_scale=HS_im;
    chu=zeros(nwl);
        for j=1:sz(2)
            if (~all(isnan(HS_im(:,j))))
                chu(1)=1; % The first convex hull point
                nch=1; % sThe total number of convex hull points found so far
                HS_scale(1,j)=0.0; % Define the cont removed as hull-original
                while chu(nch)~=nwl
                    nch=nch+1; % index of the new hull point to be found
                    % Set the point next to the last found hull point as
                    % the first trial point and set the corresponding slope
                    % as the starting slmax and the point index as the next
                    % hull point
                    slmax=(HS_im(chu(nch-1)+1,j)-HS_im(chu(nch-1),j))/(wl(chu(nch-1)+1)-wl(chu(nch-1)));
                    chu(nch)=chu(nch-1)+1;
                    for k=(chu(nch-1)+2):nwl % Go through all the following points in spectrum to find the one producing the max slope
                        sl=(HS_im(k,j)-HS_im(chu(nch-1),j))/(wl(k)-wl(chu(nch-1))); % trial slope
                        if sl>slmax % Check if trial is the maximum slope
                            slmax=sl;
                            chu(nch)=k;
                        end
                    end
                    for kint=(chu(nch-1)+1):chu(nch) % Interpolate ch points to all wavelengths between the ch points
                        HS_scale(kint,j)=HS_im(kint,j)-(HS_im(chu(nch-1),j)+slmax*(wl(kint)-wl(chu(nch-1))));
                    end
                end
            end
        end
    
    HS_scale=HS_scale-min(HS_scale,[],'all')+1;
end