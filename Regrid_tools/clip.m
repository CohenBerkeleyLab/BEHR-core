%%clip
%%arr 11/5/2007 

function value=clip(value, minvalue, maxvalue);

if value<minvalue;
    value=minvalue;
elseif value>=maxvalue;
    value=maxvalue;
end
