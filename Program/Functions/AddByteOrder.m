function AddByteOrder(hdrfn)
%Description: Function to add the line 'byte order = 0' to a header file of 
%             an ENVI format image.  The line is only added if the keyword
%             'byte order' does not already exist in the header file.
%             The line 'byte order = 0' is inserted in the line preceding
%             'interleave = ...'
%
%             The function creates a backup copy of the original header file
%             in C:\Path\to\Data\filename.hdr.bak
%
%Date : 2019-08-22
%Usage: AddByteOrder(hdrfn)
%       hdrfn is a string describing the filename of a header file
%       i.e. hdrfn = 'C:\Path\to\Data\filename.hdr'
%
%AddByteOrder uses the isempty(strfind(a,b)) check instead of contains()
%contains() was only introduced in Matlab R2016b, thus AddByteOrder is
%compatible with older Matlab versions.

  [pathstr,name,ext] = fileparts(hdrfn);
  if ~strcmpi(ext,'.hdr')
    %Check if input is a header file, otherwise return error
    error('.hdr file expected as input')
  else
    hdrtxt = lower(textread(hdrfn,'%s','delimiter','\n'));
    if any(~cellfun(@isempty,strfind(hdrtxt,'byte order')))
      %Check if header file already contains 'byte order' keyword
    %  fprintf('Header file contains the "byte order" keyword\n');
      return;
    else
      %File is a header file and does not contain the 'byte order' keyword
      bckupfn = hdrfn+'.bak';
      copyfile(hdrfn, bckupfn)
      if ~exist(bckupfn, 'file')
        %Make a backup of the header file and see if the backup exists
        error('Could not create backup of headerfile.')
      end
      fidread = fopen(bckupfn, 'rt');
      fidwrite = fopen(hdrfn, 'wt');
      while ~feof(fidread)
        tline = fgetl(fidread);
        if ~isempty(strfind(lower(tline),'interleave'))
          %Check if the current line contains the 'interleave' keyword
          %and insert 'byte order' keyword one line above
          fprintf(fidwrite, 'byte order = 0\n');
        end
        fprintf(fidwrite, [tline '\n']);
      end
      fclose(fidread);
      fclose(fidwrite);
    end
  end
end