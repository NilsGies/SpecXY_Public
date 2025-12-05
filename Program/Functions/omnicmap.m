function FileInfo = omnicmap(file)
% Read OMNIC map

% Read file as 8-bit unsigned integer
fid = fopen(file);
data = fread(fid,inf,'uint8');
fclose (fid);
clear fid

% Retrieve start and stop position, Beam Size and Map Size
strdata= char(data');
locidx = regexp(strdata,'Position ');
mapInfoOffset = locidx(1)-90;
mapDescription = typecast(uint8(strdata(mapInfoOffset:mapInfoOffset +23)),'single');
FileInfo.StartPosition = [mapDescription(1),mapDescription(4)];
FileInfo.EndPosition = [mapDescription(2),mapDescription(5)];
FileInfo.BeamSize = [mapDescription(3),mapDescription(6)];
FileInfo.MapSize = (FileInfo.EndPosition - ...
    FileInfo.StartPosition)./FileInfo.BeamSize + 1;

% Read acquisition information
fmt = 'uint32';
offset = 373;
fudgefactor = 203;
offset = typecast(uint8(data(offset:offset+3)),fmt) - fudgefactor;
IBX = (typecast(uint8(data(373:373+3)),fmt) - fudgefactor)/4;
nVal = 13;
Info = typecast(uint8(data(offset:offset + 4*nVal - 1)),'uint32');
FileInfo.nPoints = Info(1); % number of points / spectrum
FileInfo.ScanPoints = Info(7);
FileInfo.InteroPeakPosition = Info(8);
FileInfo.nSampleScans = Info(9); % number of spectra / pixel
FileInfo.FFTPoints = Info(11);
FileInfo.nBackgroundScans = Info(13); % number of spectra used to compose BKG

% Read additional info
offset = (IBX + 3)*4+1;
fmt = 'single';
nVal = 47;
Info = typecast(uint8(data(offset:offset + 4*nVal - 1)),fmt);
FileInfo.HighWN = Info(1); % Upper wavenumber
FileInfo.LowWN = Info(2); % Lower wavenumber
FileInfo.WaveNumber = double(linspace(Info(2),Info(1),FileInfo.nPoints)); % wn vector
FileInfo.bkgGain = Info(11); % Background gain
FileInfo.startidx = Info(15); % Identifier for start indices of spectra
FileInfo.LaserFrequency = Info(17); % Laser Frequency
%%
% Get spectra
specidx = regexp(strdata,'Spectrum ');
if not(isempty(specidx))

%%
FileInfo.nSpectra = numel(specidx);
%plot(dxdata)
dg=zeros(numel(FileInfo.WaveNumber),FileInfo.nSpectra);
for n=1:(FileInfo.nSpectra)-1
    dg(:,n)=typecast(uint8([data(specidx(n)+4*21:specidx(n+1)-1-4*4)])','single');
end
%     dg = reshape(dxdata,FileInfo.nSpectra,...
%     [])';
else
    specidx =regexp(strdata,'Pixel ');
    if not(isequal(numel(specidx),double(FileInfo.MapSize(1)*FileInfo.MapSize(2))))
        filter=diff(specidx)==mode(diff(specidx));
        specidx=specidx([filter(1:end) true]);
    end

    FileInfo.nSpectra = numel(specidx);

    dg=zeros(numel(FileInfo.WaveNumber),FileInfo.nSpectra);
    for n=1:(FileInfo.nSpectra)-1
        dg(:,n)=typecast(uint8([data(specidx(n)+4*21:specidx(n+1)-1-4*4)])','single');
    end
end

FileInfo.Datagrid = dg;

if FileInfo.BeamSize(1)*FileInfo.BeamSize(2) == 0
    %     FileInfo.BeamSize = [1 1];
    %     %  suggest default x and y
    %     deffac = factor(FileInfo.nSpectra);
    %     def(1) = cellstr(num2str(FileInfo.nSpectra/deffac(end)));
    %     def(2) = cellstr(num2str(deffac(end)));
    %     d = zeros(2,1);
    %     prompt = {'# x points','# y points'};
    %     dlg_title= 'Map Size';
    %     d = inputdlg(prompt,dlg_title,1,def);
    %     d = str2num(char(d));
%     FileInfo.MapSize = d;
FileInfo.MapSize=[1 FileInfo.nSpectra];
end

% Get Background Spectrum, if in file
bkg_offset = regexp(strdata,'Background');
if ~isempty(bkg_offset)
    start = (bkg_offset- FileInfo.nPoints*4 - 1)/4;
    stop = start + FileInfo.nPoints - 1;
    BKG = typecast(uint8(data),'single');
    BKG = BKG(start:stop);
    FileInfo.BackgroundSpectrum = BKG;
end

end
