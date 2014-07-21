function list_all_fields( structure )
%list_names: Lists the names of a data structure passed to this function.
%   Pass this function a data structure that has an attribute 'Name' and it
%   will list the name of each element of that structure

fields = fieldnames(structure);
for a=1:length(structure)
    fprintf('%u: \n',a);
    for b=1:length(fields)
        if isempty(eval(sprintf('structure(a).%s',fields{b})));
            fprintf('\t%s: {} \n',fields{b});
        elseif isstruct(eval(sprintf('structure(a).%s',fields{b})));
            fprintf('\t%s: %s \n', fields{b}, '[struct]'); 
        else
            fprintf('\t%s: %s \n', fields{b}, eval(sprintf('structure(a).%s',fields{b}))); 
        end
    end
end

end

