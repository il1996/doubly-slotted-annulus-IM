function A = readmatrix(fname, varargin)
%READMATRIX  Octave shim: numeric matrix from a delimited text file, skipping
%  'NumHeaderLines' header lines (default 1), tab/space delimited.
  nh = 1; delim = '';
  for k = 1:2:numel(varargin)
    switch lower(varargin{k})
      case 'numheaderlines', nh = varargin{k+1};
      case 'delimiter', delim = varargin{k+1};
    end
  end
  if strcmp(delim, '\t'), delim = sprintf('\t'); end
  A = dlmread(fname, delim, nh, 0);
end
