function [svgpts, sppts] = svgImport(fname, pathid)
%SVG_IMPORT Reads the points from an SVG file and returns them
% in in x-y coordinates. if name/id of the svg-path of interest 
% (=drawing object!) is unknown and empty, function will assume that
% there is only one path in svg file
% INPUT ARGUMENTS:
% fpath     - target file path
% fname     - target file name
% pts       - outline, vector of points [x, y]
% pathid    - id/name of the svg path object (string)
% OUTPUT ARGUMENT:
% svgpts    - array of outline points (X,Y)

%% read svg file
file=fileread(fname);

%% find and extract all data from all path objects
pathidx = strfind(file, '<path');
for ii=1:numel(pathidx)
    endidx = strfind(file(pathidx(ii):end), '/>');
    pdat{ii}=file(pathidx(ii):pathidx(ii)+endidx(1));
    
    %% find + extract ID
    ididx = strfind(pdat{ii}, 'id="');
    ididx2 = strfind(pdat{ii}(ididx+4:end), '"');
    path{ii}.id= pdat{ii}(ididx+4:ididx+3+ididx2(1));
    %% find + extract data points from path
    pidx = strfind(pdat{ii}, 'd="');
    pidx2 = strfind(pdat{ii}(pidx(1)+3:end), '"');
    dbuf = pdat{ii}(pidx(1)+3:pidx+3+pidx2(1));
    % get rid of non-value path chars for easier readout (only straight
    % lines allowed, no bezier etc)
    nvalc={',', '"', 'C', 'S', 'Q', 'T', 'A', 'Z'};
    for jj=1:numel(nvalc)
        dbuf =  strrep(dbuf, nvalc{jj}, ' ');
    end
    % search and index valid path commands 
    valc = {'M', 'L', 'H', 'V'};
    fidx = [];
    for jj=1:numel(valc)
        fidx = [fidx strfind(dbuf, valc{jj})];
    end
    fidx = [fidx numel(dbuf)];
    fidx=sort(fidx);
    
    pts = [];
    for jj=1:numel(fidx)
       switch dbuf(fidx(jj))
           case 'M'
               pts = [pts; sscanf(dbuf(fidx(jj)+1:fidx(jj+1)), '%f')];
           case 'L'
               pts = [pts; sscanf(dbuf(fidx(jj)+1:fidx(jj+1)), '%f')];
           case 'H' %horizontal lines, only horizontal coordinates given
               ptsbuf = sscanf(dbuf(fidx(jj)+1:fidx(jj+1)), '%f');
               cmplt(1:2:2*numel(ptsbuf)) = ptsbuf;
               cmplt(2:2:2*numel(ptsbuf)) = pts(end); %use last known y coordinate
               cmplt=cmplt';
               pts = [pts; cmplt];
           case 'V' %vertical lines, only vertical coordinates given
               ptsbuf = sscanf(dbuf(fidx(jj)+1:fidx(jj+1)), '%f');
               cmplt(1:2:2*numel(ptsbuf)) = pts(end-1); %use last known x coordinate
               cmplt(2:2:2*numel(ptsbuf)) = ptsbuf;
               cmplt=cmplt';
               pts = [pts; cmplt];           
           otherwise
               % nothing, just for indexing purposes jj+1
       end
    end
     
    % read in data points and store (x,y)
    path{ii}.data(:,1)=pts(1:2:end);
    path{ii}.data(:,2)=pts(2:2:end);
    
   if ~isempty(strfind(path{ii}.id, pathid)) || isempty(pathid)
       svgpts = path{ii}.data;
       sppts = size(path{ii}.data,1);
   end
 
end



end

