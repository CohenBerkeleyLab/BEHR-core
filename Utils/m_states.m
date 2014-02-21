function m_states( color, varargin )
%m_states.m     Josh Laughner <joshlaugh5@gmail.com> 21 Feb 2014
%Add states: adds state outlines to a map of the US using the m_map package
%for MATLAB. Call this function after the coastlines have already been
%drawn. Takes 2 arguments:
%   color: a normal MATLAB color string that will set the border color
%   varargin: state abbreviations (upper or lower case) that will define
%   the states to be drawn.  Set to 'all' or leave blank to draw all state
%   borders. Set to 'cont' to only draw the 48 continental states.

directory = '/Users/Josh/Documents/MATLAB/BEHR/Utils/mapping_tools'; %Set this directory to wherever the list of 'xxx2pts.mat' files is stored.

if isempty(varargin) 
    states = 'all';
else
    states = varargin;
end

state_abbrev = {'al','ak','az','ar','ca','co','ct','de','fl','ga','hi','id','il','in','ia','ks','ky','la','me','md','ma','mi','mn','ms','mo','mt','ne','nv','nh','nj','nm','ny','nc','nd','oh','ok','or','pa','ri','sc','sd','tn','tx','ut','vt','va','wa','wv','wi','wy'};
state_mats = {'alabama2pts.mat','alaska2pts.mat','arizona2pts.mat','arkansas2pts.mat','california2pts.mat','colorado2pts.mat','connecticut2pts.mat',...
    'delaware2pts.mat','florida2pts.mat','georgia2pts.mat','hawaii2pts.mat','idaho2pts.mat','illinois2pts.mat','indiana2pts.mat','iowa2pts.mat','kansas2pts.mat',...
    'kentucky2pts.mat','louisiana2pts.mat','maine2pts.mat','maryland2pts.mat','massachusetts2pts.mat','michigan2pts.mat','minnesota2pts.mat','mississippi2pts.mat',...
    'missouri2pts.mat','montana2pts.mat','nebraska2pts.mat','nevada2pts.mat','new_hampshire2pts.mat','new_jersey2pts.mat','new_mexico2pts.mat','new_york2pts.mat','north_carolina2pts.mat','north_dakota2pts.mat',...
    'ohio2pts.mat','oklahoma2pts.mat','oregon2pts.mat','pennsylvania2pts.mat','rhode_island2pts.mat','south_carolina2pts.mat','south_dakota2pts.mat','tennessee2pts.mat','texas2pts.mat','utah2pts.mat','vermont2pts.mat','virginia2pts.mat','washington2pts.mat','west_virginia2pts.mat',...
    'wisconsin2pts.mat','wyoming2pts.mat'};

if strcmpi(states,'all')
   to_plot = state_abbrev;
elseif strcmpi(states,'cont')
   to_plot = state_abbrev(~strcmpi('ak',state_abbrev) & ~strcmpi('hi',state_abbrev));
else
   to_plot = states;
end

for s = to_plot
    index = find(strcmp(state_abbrev,s));
    matfile = state_mats{index};
    load(matfile);
    m_line(bndry_lon,bndry_lat,'color',color)
end

end

