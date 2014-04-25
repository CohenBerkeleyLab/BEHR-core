function list_names( structure )
%list_names: Lists the names of a data structure passed to this function.
%   Pass this function a data structure that has an attribute 'Name' and it
%   will list the name of each element of that structure

for a=1:length(structure)
   fprintf('%u : %s \n', a, structure(a).Name); 
end

end

