function tf = contains(str, pat, varargin)
%CONTAINS  Octave shim: true if pat occurs in str (char or cellstr).
  ignore = false;
  for k = 1:2:numel(varargin)
    if strcmpi(varargin{k}, 'IgnoreCase'), ignore = varargin{k+1}; end
  end
  if ischar(str), str = {str}; single = true; else, single = false; end
  if ischar(pat), pat = {pat}; end
  tf = false(size(str));
  for i = 1:numel(str)
    s = str{i};
    for j = 1:numel(pat)
      p = pat{j};
      if ignore, ok = ~isempty(strfind(lower(s), lower(p))); else, ok = ~isempty(strfind(s, p)); end
      tf(i) = tf(i) || ok;
    end
  end
  if single, tf = tf(1); end
end
