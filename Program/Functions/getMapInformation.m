function ddict = getMapInformation(data)
    % getMapInformation Extracts map geometrical acquisition parameters from data
    %
    % Parameters:
    % -----------
    % data : The contents of the .map file (as a byte array)
    %
    % Returns:
    % --------
    % ddict : A struct with map geometrical acquisition parameters

    % Determine the chain 'Position'
    chain = 'Position';
    if verLessThan('matlab', '9.0') % Equivalent of checking for Python 3.0
        chain = 'Position';
    else
        chain = unicode2native(chain, 'UTF-8');
    end
    
    % Locate the 'Position' string in the data
    offset = strfind(data, chain);
    positions = offset;

    % Loop to find all occurrences of the chain
    while true
        try
            a = strfind(data((offset(end) + 1):end), chain);
            if isempty(a)
                break;
            end
            offset = offset(end) + a(1);
            positions(end + 1) = offset; %#ok<AGROW>
        catch
            break;
        end
    end

    ddict = struct();

    % Map description position
    if (positions(2) - positions(1)) == 66  % reverse engineered magic number :-)
        mapDescriptionOffset = positions(1) - 90;
        mapDescription = typecast(data(mapDescriptionOffset+1 : mapDescriptionOffset + 24), 'single');
        y0 = mapDescription(1);
        y1 = mapDescription(2);
        deltaY = mapDescription(3);
        x0 = mapDescription(4);
        x1 = mapDescription(5);
        deltaX = mapDescription(6);

        ddict.FirstMapLocation = [x0, y0];
        ddict.LastMapLocation = [x1, y1];
        ddict.MappingStageXStepSize = deltaX;
        ddict.MappingStageYStepSize = deltaY;
        ddict.NumberOfSpectra = abs((1 + ((y1 - y0) / deltaY)) * (1 + ((x1 - x0) / deltaX)));
    end

    % Optional: Add debug information
    fields = fieldnames(ddict);
    for i = 1:numel(fields)
        fieldName = fields{i};
        disp([fieldName ': ' num2str(ddict.(fieldName))]);
    end
end
