function [sp,ts,ms,t_ms,ts_ms,mc,t] = XOVER(t,m,ivar,tvars,hzlimit)
%XOVER Performs a cross-over analysis between cruise datasets.
% Documentation and other supporting material is freely available online:
%  http://dx.doi.org/10.13140/RG.2.1.1629.0405
% Please cite as: Humphreys, Matthew P. (2015). "Cross-over analysis of
%  hydrographic variables: XOVER v1.0". Ocean and Earth Science, University
%  of Southampton, UK. doi:10.13140/RG.2.1.1629.0405
% Correspondence: m.p.humphreys@soton.ac.uk
% Last updated 2017-03-08 [v1.0.1]
% 
% === INPUTS ===
%       t: The test cruise data. [table]
%          Required columns:
%           t.cruise  = identifier that is unique to each cruise
%                        [double|cell of strings]
%           t.station = identifier that is unique to each station within
%                        each cruise (can be duplicated between cruises).
%                        [double|cell of strings]
%           t.lat     = latitude (-90 to +90), in degrees N. [double]
%           t.lon     = longitude (-180 to +180), in degrees E. [double]
%           t.ndate   = sample date and time, in MATLAB datenum() format
%                       [double]
%       m: The master cruise data to compare test data with. [table]
%          Required columns are the same as for input t.
%          If you wish to only use a subset of these data for the analysis,
%           e.g. only below a certain depth, you should carry out that
%           selection manually, and input this subset as m.
%    ivar: Variable name to use as interpolant. [string]
%   tvars: Variable names to compare between cruises. [cell of strings]
% hzlimit: Maximum distance between cross-over stations, in km. [double]
% 
% === OUTPUTS ===
% Each output contains results of the cross-over analysis, organised as
%  follows:
%      sp: List of matching station pairs [table]
%      ts: List of stations in t [table]
%      ms: List of stations in m local to t [table]
%    t_ms: matrices of t v. ms [structure]
%   ts_ms: matrices of ts v. ms [structure]
%      mc: List of matching cruises and statistics [table]
%       t: Updated version of input t including only data used in the
%           cross-over analysis; rows match rows of t_ms [table]
% 

disp('XOVER: performing cross-over analysis...');

%% Select only <m> and <t> stations with data for the variables of interest
tf = t;
m = m(any(~isnan(table2array(m(:,tvars))),2) & ~isnan(m.(ivar)),:);
t = t(any(~isnan(table2array(t(:,tvars))),2) & ~isnan(t.(ivar)),:);
if height(t) == 0
    disp('XOVER: error - input <t> does not contain sufficient data!');
    sp = [];
    ts = [];
    ms = [];
    t_ms = [];
    ts_ms = [];
    mc = [];
    t = tf;
else
    
%% Convert 'cruise' and 'station' fields to cells of strings, if required
t.cruise  = XOVER_num2txt(t.cruise);
t.station = XOVER_num2txt(t.station);
m.cruise  = XOVER_num2txt(m.cruise);
m.station = XOVER_num2txt(m.station);

%% Generate 1-per-station table for <t>: <ts>
ts = unique(t(:,{'cruise' 'station'}));
tsvars = {'lat' 'lon' 'ndate'};
for V = 1:length(tsvars)
    ts.(tsvars{V}) = NaN(height(ts),1);
    for S = 1:height(ts)
        SL = strcmp(t.cruise,ts.cruise{S}) ...
            & strcmp(t.station,ts.station{S});
        ts.(tsvars{V})(S) = mean(t.(tsvars{V})(SL));
    end %for S
end %for V
ts = sortrows(ts,'ndate');
clear tsvars V S SL

%% Convert co-ordinates
t.sxyz0 = XOVER_ll2s(t.lat,t.lon);
ts.sxyz0 = XOVER_ll2s(ts.lat,ts.lon);
m.sxyz0 = XOVER_ll2s(m.lat,m.lon);

%% Cut <m> to within <hzlimit> of <ts>
tsknn = fitcknn(ts.sxyz0,(1:height(ts))');
m.nearest_test = predict(tsknn,m.sxyz0);
m.hzdistance = sqrt(sum((m.sxyz0 - ts.sxyz0(m.nearest_test,:)).^2,2));
m = m(m.hzdistance < hzlimit,:);
clear tsknn

%% Generate 1-per-station table for <m>: <ms>
ms = unique(m(:,{'cruise' 'station'}));
msfields = {'lat' 'lon' 'ndate' 'sxyz0'};
for F = 1:length(msfields)
    ms.(msfields{F}) = NaN(height(ms),size(m.(msfields{F}),2));
    for T = 1:height(ms)
        TL = strcmp(m.cruise,ms.cruise{T}) ...
            & strcmp(m.station,ms.station{T});
        ms.(msfields{F})(T,:) = mean(m.(msfields{F})(TL,:));
    end %for T
end %for F
clear msfields F T TL

%% Get nearby cruises: <ts_ms>
% Rows of <ts_ms> fields correspond to rows of <ts>; columns correspond to
%  rows of <ms>
ts_ms.dist = NaN(height(ts),height(ms));
for S = 1:height(ts)
    ts_ms.dist(S,:) = sqrt(sum((ms.sxyz0 - repmat(ts.sxyz0(S,:), ...
        height(ms),1)).^2,2));
end %for S
ts_ms.match = ts_ms.dist <= hzlimit;
clear S

%% Generate <ms> interpolants
for V = 1:length(tvars)
    ms.(tvars{V}) = cell(height(ms),1);
    ms.([tvars{V} '_interp']) = cell(height(ms),1);
    VL = ~isnan(m.(tvars{V}));
    for S = 1:height(ms)
        SVL = VL & strcmp(m.station,ms.station{S}) ...
            & strcmp(m.cruise,ms.cruise{S});
        
      % Sanitise interpolation inputs
        SVindexVar = unique(m.(ivar)(SVL));
        SVtestVar = NaN(size(SVindexVar));
        for M = 1:length(SVtestVar)
            SVtestVar(M) = mean(m.(tvars{V})(SVL ...
                & m.(ivar) == SVindexVar(M)));
        end %for M
            
        ms.(tvars{V}){S} = [SVindexVar SVtestVar];
        
      % Generate interpolator, if there are sufficient data
        if length(SVtestVar) > 1
            ms.([tvars{V} '_interp']){S} = pchip(SVindexVar,SVtestVar);
        end %if
        
    end %for S
end %for V

clear SV* S M V VL

%% Set up station pairs table <sp>
[t_match,m_match] = ind2sub(size(ts_ms.match),find(ts_ms.match));
spfields = {'cruise' 'station' 'lat' 'lon' 'ndate'};
for F = 1:length(spfields)
    sp.(['t_' spfields{F}]) = ts.(spfields{F})(t_match);
end %for F
for F = 1:length(spfields)
    sp.(['m_' spfields{F}]) = ms.(spfields{F})(m_match);
end %for F
sp.hzdist = sqrt(sum((ms.sxyz0(m_match,:) - ts.sxyz0(t_match,:)).^2,2));
sp = XOVER_s2t(sp);
clear *_match spfields F

%% Generate 1-per-cruise table for master data: <mc>
mc.cruise = unique(ms.cruise);
mc = XOVER_s2t(mc);
mc.ndate = NaN(height(mc),1);
mc.ndate_range = NaN(height(mc),2);
for C = 1:height(mc)
    CL = strcmp(mc.cruise{C},ms.cruise);
    mc.ndate(C) = mean(ms.ndate(CL));
    mc.ndate_range(C,:) = minmax(ms.ndate(CL)');
end %for C
mc = sortrows(mc,'ndate');
clear C CL

%% Interpolate
for V = 1:length(tvars)
    
  % Carry out all interpolations
    t_ms.(tvars{V}) = NaN(height(t),height(ms));
    t_ms.([tvars{V} '_resid']) = NaN(height(t),height(ms));
    for M = 1:height(ms)
        Minterp = ms.([tvars{V} '_interp']){M};
        if ~isempty(Minterp)
            ML = t.(ivar) >= min(Minterp.breaks) ...
                & t.(ivar) <= max(Minterp.breaks) ...
                & ismember(t.station,ts.station(ts_ms.match(:,M))) ...
                & ismember(t.cruise,ts.cruise(ts_ms.match(:,M)));
            t_ms.(tvars{V})(ML,M) = ppval(Minterp,t.(ivar)(ML));
            t_ms.([tvars{V} '_resid'])(ML,M) = t.(tvars{V})(ML) ...
                - ppval(Minterp,t.(ivar)(ML));
        end %if
    end %for M
       
  % Put station pair-by-pair mean resids in <ts_ms>
    ts_ms.([tvars{V} '_resid']) = NaN(height(ts),height(ms));
    ts_ms.([tvars{V} '_std']) = NaN(height(ts),height(ms));
    ts_ms.([tvars{V} '_nobs']) = NaN(height(ts),height(ms));
    for T = 1:height(ts)
        TL = strcmp(t.station,ts.station{T}) ...
            & strcmp(t.cruise,ts.cruise{T});
        for M = 1:height(ms)
            MTL = TL & ~isnan(t_ms.(tvars{V})(:,M));
            ts_ms.([tvars{V} '_resid'])(T,M) ...
                = mean(t_ms.([tvars{V} '_resid'])(MTL,M));
            ts_ms.([tvars{V} '_std'])(T,M) ...
                = std(t_ms.([tvars{V} '_resid'])(MTL,M));
            ts_ms.([tvars{V} '_nobs'])(T,M) = sum(MTL);
        end %for M
    end %for T
    
    sp.([tvars{V} '_resid']) = ts_ms.([tvars{V} '_resid'])(ts_ms.match);
    sp.([tvars{V} '_std']) = ts_ms.([tvars{V} '_std'])(ts_ms.match);
    sp.([tvars{V} '_nobs']) = ts_ms.([tvars{V} '_nobs'])(ts_ms.match);
    
  % Put resids into <ms>
    ms.([tvars{V} '_resid']) = nanmean(t_ms.([tvars{V} '_resid']),1)';
    ms.([tvars{V} '_std']) = nanstd(t_ms.([tvars{V} '_resid']),1)';
    ms.([tvars{V} '_nobs']) = sum(~isnan(t_ms.([tvars{V} '_resid'])),1)';
    
  % Put resids into <mc>
    mc.([tvars{V} '_resid']) = NaN(height(mc),1);
    mc.([tvars{V} '_std']) = NaN(height(mc),1);
    mc.([tvars{V} '_nobs']) = NaN(height(mc),1);
    
    for C = 1:height(mc)
        CL = strcmp(ms.cruise,mc.cruise{C});
        Cresid = t_ms.([tvars{V} '_resid'])(:,CL);
        mc.([tvars{V} '_resid'])(C) = nanmean(Cresid(:));
        mc.([tvars{V} '_std'])(C) = nanstd(Cresid(:));
        mc.([tvars{V} '_nobs'])(C) = sum(~isnan(Cresid(:)));
    end %for C
    
end %for V

clear V M ML Minterp MTL T TL C CL Cresid

disp('XOVER: complete!');

end %if height(t) == 0

end %function XOVER

function sxyz = XOVER_ll2s(lat,lon)
r = 6371; % assume constant Earth radius of 6371 km
[sx,sy,sz] = sph2cart(degtorad(lon),degtorad(lat),r);
sxyz = [sx sy sz];
end %function XOVER_llz2s

function tab = XOVER_s2t(st)
if isa(st,'struct'), tab = struct2table(st);
else, tab = st;
end %if
end %function XOVER_s2t

function txt = XOVER_num2txt(num)
if isa(num,'double'), txt = cellstr(num2str(num,'%-u'));
else, txt = num; end
end %function XOVER_num2txt