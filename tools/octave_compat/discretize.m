function idx = discretize(x, edges)
%DISCRETIZE  Octave shim for MATLAB's discretize(x, edges): bin index of x
%  in the half-open intervals [edges(i), edges(i+1)); the last bin is closed.
%  Values outside return NaN.
  edges = edges(:).';
  idx = nan(size(x));
  n = numel(edges) - 1;
  for i = 1:n
    if i < n
      m = (x >= edges(i)) & (x < edges(i+1));
    else
      m = (x >= edges(i)) & (x <= edges(i+1));
    end
    idx(m) = i;
  end
end
