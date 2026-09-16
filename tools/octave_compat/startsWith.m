function tf = startsWith(str, pat)
  if ischar(str), str = {str}; single=true; else, single=false; end
  if ischar(pat), pat = {pat}; end
  tf = false(size(str));
  for i=1:numel(str), for j=1:numel(pat), tf(i) = tf(i) || strncmp(str{i}, pat{j}, numel(pat{j})); end, end
  if single, tf=tf(1); end
end
