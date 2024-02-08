function [Spectra, Wavenumbers, SpectraTitles, pathname, ...
    SpectraComments] = spa_load(pathname)
% Modified after
% LoadSpectra.m by Kurt Oldenburg - 06/28/16
% Imports the absorbance data in .SPA spectrum files 
% 
fid=fopen(pathname,'r');
    

    
    fseek(fid,30,'bof');
    SpectraTitles={char(nonzeros(fread(fid,255,'uint8'))')};
    fseek(fid,564,'bof');
    Spectrum_Pts=fread(fid,1,'int32');
    fseek(fid,576,'bof');
    Max_Wavenum=fread(fid,1,'single');
    Min_Wavenum=fread(fid,1,'single');
    % The Wavenumber values are assumed to be linearly spaced between
    % between the Min and Max values. The array needs to be flipped 
    % around to get the order lined up with the absorbance data.
    
    Wavenumbers(:,1)=flipud(linspace(Min_Wavenum,...
        Max_Wavenum,Spectrum_Pts).')';
    
    % The starting byte location of the absorbance data is stored in the
    % header. It immediately follows a flag value of 3:
    
    Flag=0; 
    
    fseek(fid,288,'bof');
    
    while Flag ~= 3         
        Flag = fread(fid,1,'uint16');   
    end;
    
    DataPosition=fread(fid,1,'uint16')';
    fseek(fid,DataPosition,'bof');
    
    Spectra(:,1)=fread(fid,Spectrum_Pts,'single');
    % Same story goes for the Comments section with a flag of 4.
    % The size of the section is the difference between the two.
    
    Flag=0;
    
    fseek(fid,288,'bof');
    
    while Flag ~= 4         
        Flag = fread(fid,1,'uint16');   
    end 
    
    CommentPosition=fread(fid,1,'uint16')';
    SpectraComments={char(nonzeros(fread(fid, ...
        (DataPosition-DataPosition), 'uint8'))')};
    
    fclose(fid);
    
