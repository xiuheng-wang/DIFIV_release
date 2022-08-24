
function save_results( filename, X )
  fid = fopen(filename, 'w');
%   fprintf(fid, '%s', X{1,1});
%   for i = 2:size(X,2)
%       fprintf(fid, '\t%s', X{1,i});
%   end
%   fprintf(fid, '\n');
  for i=2:size(X,1)
    fprintf(fid, '%s', X{i,1});
    fprintf(fid, '\t%2.3f', X{i,2:end});
    fprintf(fid, '\n');
  end
  fclose(fid);
end