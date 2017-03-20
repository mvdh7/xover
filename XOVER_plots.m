%XOVER_plots v1.0
% Example script to performs a cross-over analysis between cruise datasets,
%  using the function XOVER v1.0.
% Also requires the function loadGLODAPv1, available at the File Exchange:
%  http://uk.mathworks.com/matlabcentral/fileexchange/51961-loading-glodap-into-matlab
% Documentation: Humphreys, Matthew P. (2015). "Cross-over analysis of
%  hydrographic variables: XOVER v1.0". Ocean and Earth Science, University
%  of Southampton, UK. doi:10.13140/RG.2.1.1629.0405
% Correspondence: m.p.humphreys@soton.ac.uk
% Last updated 2015-07-21 [v1.0.0.2]

%% Example - load GLODAP data as master data <m>
% Function loadGLODAPv1 is available from:
%  http://uk.mathworks.com/matlabcentral/fileexchange/51961-loading-glodap-into-matlab
% Documentation: http://dx.doi.org/10.13140/RG.2.1.2681.9687
m = loadGLODAPv1('ATL');
m = m(all(~isnan([m.lat m.lon m.ndate]),2),:);

% Select GLODAPv1 cruise #28 as test data <t>
tcruise = 28;
t = m(m.cruise == tcruise,:);
% Remove this cruise from master data <m> (so that a cross-over is not 
%  performed against itself!)
m = m(m.cruise ~= tcruise,:);

% Select only test and master data deeper than 1500m
m = m(m.z > 1500,:);
t = t(t.z > 1500,:);

% Define XOVER inputs
ivar = 's4'; % index variable ('s4' = potential density at 4000 dbar)
tvars = {'dic' 'ta'}; % test variables ('dic' = dissolved inorganic carbon,
                      %  'ta' = total alkalinity)
hzlimit = 100; % max. distance between matching stations, in km

% Run XOVER
[sp,ts,ms,t_ms,ts_ms,mc,t] = XOVER(t,m,ivar,tvars,hzlimit);

%% 1 Individual interpolations
% Choose which one of <tvars> to show
f1tvar = 'dic';
% Choose which master station to show
f1ms_station = 1 + f1ms_station;
% Make figure
figure(1); clf; hold on;
% Set logicals etc. for figure
f1L = strcmp(ms.station(f1ms_station),sp.m_station);
f1Lix = find(f1L);
f1clr = {'r' 'b' 'g' 'm' [1 0.5 0] 'c' 'y' 'k'};
f1clr = repmat(f1clr,1,50);
f1marker = {'x' '+' 'o' 'sq' 'd' '^' 'v'};
f1marker = repmat(f1marker,1,50);
% Plot test data from stations matching the master station
for T = 1:sum(f1L)
    % Plot test data values, as coloured markers (different colour for each
    %  test station)
    TL = strcmp(t.station,sp.t_station{f1Lix(T)});
    scatter(t.(ivar)(TL),t.(f1tvar)(TL),30,f1clr{T}, ...
        'marker',f1marker{T}, 'linewidth',2);
    % Plot matches, as coloured dashed lines
    if ~isempty(ms.([f1tvar '_interp']){f1ms_station})
    Tx = [t.(ivar)(TL) t.(ivar)(TL) NaN(sum(TL),1)]';
        Ty = [t.(f1tvar)(TL) ...
            ppval(ms.([f1tvar '_interp']){f1ms_station},t.(ivar)(TL)) ...
            NaN(sum(TL),1)]';
        Tx(Tx < min(ms.([f1tvar '_interp']){f1ms_station}.breaks) ...
            | Tx > max(ms.([f1tvar '_interp']){f1ms_station}.breaks)) ...
            = NaN;
        plot(Tx(:),Ty(:), 'color',f1clr{T}, 'linestyle',':');
    end %if
end %for T
% Plot master data and interpolation, as black points and line
if ~isempty(ms.([f1tvar '_interp']){f1ms_station})
scatter(ms.(f1tvar){f1ms_station}(:,1), ...
    ms.(f1tvar){f1ms_station}(:,2),30,'ko','filled');
f1x = linspace(min(ms.([f1tvar '_interp']){f1ms_station}.breaks), ...
    max(ms.([f1tvar '_interp']){f1ms_station}.breaks),1000);
    plot(f1x,ppval(ms.([f1tvar '_interp']){f1ms_station},f1x),'k', ...
        'linewidth',2);
end %if
% Axis labels & title
xlabel(ivar);
ylabel(f1tvar);
title(['Master cruise ' ms.cruise{f1ms_station} ', station ' ...
    num2str(ms.station{f1ms_station})]);
% Grid
grid on

%% 2 Map of cross-overs
figure(2); clf; hold on;
% Generate map
worldmap(minmax(t.lat')+[-3 3],minmax(t.lon')+[-3 3]);
land = shaperead('landareas', 'usegeocoords',true);
geoshow(land, 'facecolor','k', 'edgecolor','none');
% Master stations which have data for <tvars>, as grey crosses
f2logic = any(~isnan(table2array(m(:,tvars))),2);
scatterm(m.lat(f2logic),m.lon(f2logic),30,0.7*[1 1 1], 'marker','x');
% All test stations, as light blue plusses
scatterm(ts.lat,ts.lon,30,[0.4 0.7 1], 'marker','+');
% Master stations within <hzlimit> of test stations, as black crosses
scatterm(sp.m_lat,sp.m_lon,30,'kx');
% Test stations involved in cross-over, as dark blue plusses
scatterm(sp.t_lat,sp.t_lon,30,'b+');
% Plot matching pairs, as red lines
f2lats = [sp.m_lat sp.t_lat NaN(size(sp.t_lon))]';
f2lons = [sp.m_lon sp.t_lon NaN(size(sp.t_lon))]';
plotm(f2lats(:),f2lons(:),'r');

%% 3 Histograms of residuals
% Choose which one of <tvars> to show
f3var = 'ta';
% Make figure
figure(3); clf;
% Subplot 1: mean of all master station residuals at each test data point
subplot(1,3,1); hold on;
    % Get data for subplot
    f3data1 = nanmean(t_ms.([f3var '_resid']),2);
    % Create histogram
    histogram(f3data1);
    % Axis labels & title
    xlabel([f3var ' residual']);
    ylabel('Frequency');
    title('Mean residual at each test data point');
    % Stats
    [~,f3pval1] = ttest(f3data1);
    text(0.02,0.96,{['Mean ± SD = ' num2str(nanmean(f3data1)) ' ± ' ...
        num2str(nanstd(f3data1))] ['\itt\rm-test \itp\rm-value = ' ...
        num2str(f3pval1)] ['No. of test data points = ' ...
        num2str(sum(~isnan(f3data1)))]}, 'units','normalized');
    % Extras
	grid on
    plot([0 0],get(gca,'ylim'),'k');
% Subplot 2: mean of all master station residuals at each test station
subplot(1,3,2); hold on;
    % Get data for subplot
    f3data2 = nanmean(ts_ms.([f3var '_resid']),2);
    % Create histogram
    histogram(f3data2);
    % Axis labels & title
    xlabel([f3var ' residual']);
    ylabel('Frequency');
    title('Mean residual at each test station');
    % Stats
    [~,f3pval2] = ttest(f3data2);
    text(0.02,0.96,{['Mean ± SD = ' num2str(nanmean(f3data2)) ' ± ' ...
        num2str(nanstd(f3data2))] ['\itt\rm-test \itp\rm-value = ' ...
        num2str(f3pval2)] ['No. of test stations = ' ...
        num2str(sum(~isnan(f3data2)))]}, 'units','normalized');
    % Extras
    grid on
    plot([0 0],get(gca,'ylim'),'k');
% Subplot 2: mean of all test data point residuals at each master station    
subplot(1,3,3); hold on;
    % Get data for subplot
    f3data3 = nanmean(t_ms.([f3var '_resid']),1);
    % Create histogram
    histogram(f3data3);
    % Axis labels & title
    xlabel([f3var ' residual']);
    ylabel('Frequency');
    title('Mean residual at each master station');
    % Stats
    [~,f3pval3] = ttest(f3data3);
    text(0.02,0.96,{['Mean ± SD = ' num2str(nanmean(f3data3)) ' ± ' ...
        num2str(nanstd(f3data3))] ['\itt\rm-test \itp\rm-value = ' ...
        num2str(f3pval3)] ['No. of master stations = ' ...
        num2str(sum(~isnan(f3data3)))]}, 'units','normalized');    
    % Extras
    grid on
    plot([0 0],get(gca,'ylim'),'k');
    
%% 4 Cruise residuals
% Choose which one of <tvars> to show
f4var = 'ta';
% Calculate regression
f4reg = regstats(mc.([f4var '_resid']),mc.ndate, ...
    'linear',{'beta' 'rsquare'});
% Make figure
figure(4); clf; hold on;
% Scatter master cruise mean residuals
errorbar(mc.ndate,mc.([f4var '_resid']),mc.([f4var '_std']),'b', ...
    'linestyle','none', 'marker','o');
% Date ticks
datetick x
% Plot test cruise date
plot(mean(ts.ndate)*[1 1],get(gca,'ylim'),'r--');
% Plot secular change regression
plot(get(gca,'xlim'),get(gca,'xlim')*f4reg.beta(2)+f4reg.beta(1),'r');
% Extras
grid on
plot(get(gca,'xlim'),[0 0],'k');
% Axis labels & legend
xlabel('Master cruise date');
ylabel([f4var ' residual (test ' endash ' master)']);
l=legend('Mean ± SD of cruise residuals','Date of test cruise', ...
    ['OLS regression: rate = ' num2str(f4reg.beta(2)*365.25) ...
    ' yr^{-1}; \itr\rm^2 = ' num2str(f4reg.rsquare,'%.3f') ...
    '; test cruise intercept = ' num2str(mean(ts.ndate)*f4reg.beta(2) ...
    + f4reg.beta(1)) '; \itn\rm = ' ...
    num2str(sum(~isnan(mc.([f4var '_resid']))))], ...
    'location','northoutside');

%% 5 Variation of test data point residuals with other test data variables
% Choose which one of <tvars> to show
f5var = 'dic';
% Choose which other variables to compare it with
f5others = {'lat' 'lon' 'z' f5var};
% Make figure
figure(5); clf;
set(gcf, 'color','w');
for S = 1:length(f5others)
% Do OLS regression
Sreg = regstats(nanmean(t_ms.([f5var '_resid']),2), ...
    t.(f5others{S}),'linear',{'beta' 'rsquare'});
% Make subplots
subplot(2,2,S); hold on;
    % Scatter data
    scatter(t.(f5others{S}),nanmean(t_ms.([f5var '_resid']),2),20,'b+');
    scatter(t.(f5others{S}),nanmean(t_ms.([f5var '_resid']),2),20, ...
        t.ndate,'filled','markeredgecolor','k');
    % Plot regression
    plot(get(gca,'xlim'),get(gca,'xlim')*Sreg.beta(2)+Sreg.beta(1),'r');
    % Labels and stats
    xlabel(f5others{S});
    ylabel([f5var ' residual']);
    text(0.02,0.9,{[f5var ' resid. = (' num2str(Sreg.beta(2)) ...
        '×' f5others{S} ')' num2str(Sreg.beta(1),'%+f')] ...
        ['\itr\rm^2 = ' num2str(Sreg.rsquare,'%.3f')]}, ...
        'units','normalized', 'color','r');
    % Grid
    grid on
end %for S