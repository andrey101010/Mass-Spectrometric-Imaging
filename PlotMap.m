addpath('..\subroutines\'); % add in Acq subroutines

%% data definition

% Sample: Bradykinin
% where in water
% Sample (UV 2-A): 16kV
% Extraction (UV 2-B): 11kV

filename = '\\win.desy.de\group\cfel\4all\mpsd_drive\massspec\VacuumDIVE\ExperimentalData\TOFData\\2020-10-06\03.signal.div';
logfilename = '\\win.desy.de\group\cfel\4all\mpsd_drive\massspec\VacuumDIVE\ExperimentalData\TOFData\2020-10-06\03.log';
times = [19.9, 19.9, 19.9, 19.9]; % 22.5 for angiotensin
windows = [0.3, 0.2, 0.2, 0.3];
signalBounds = [];

%% data readout
% read logfile
hLogFile = fopen(logfilename, 'r');
logdate = fscanf(hLogFile, 'Date: %s %s\n');
logtriggersperframe = fscanf(hLogFile, 'Triggers per step: %u\n');
logcoordinates = fscanf(hLogFile, '%f %f\n', [2, Inf])';
fclose(hLogFile);

mapwidth = find(abs(logcoordinates(:,2)-logcoordinates(1,2)) > 5, 1, 'first')-1; % finde length of one scan row
mapheight = size(logcoordinates, 1)/mapwidth;

% read signal file
M = length(times);
hFile = fopen(filename, 'r');
header = AcqReadHeader(hFile);
t = AcqGetTime(header);

% determine search window bounds
dx = t(2)-t(1);
windowsIdx = [-1;+1].*ceil(windows./dx);

% find center mass for each peak search window
timeIdx = zeros(1,M);
for k=1:M
    timesIdx(k) = find(times(k)<=t, 1, 'first');
end

% determine peak heights (i,k) where i is the position index, k the mass index
N = floor(header.nFramesCount/logtriggersperframe);
peakHeights = zeros(N, M);
for i=1:N
    frameBounds = [1+(i-1)*logtriggersperframe, i*logtriggersperframe];
    [~, yy] = AcqReadAverage(hFile, header, frameBounds, signalBounds, [], []);
    %figure, plot(t,yy)
    baseline = mean(yy(1:100));
    for k=1:M
        peakHeights(i, k) = max(0, max(yy(timesIdx(k)+(windowsIdx(1,k):windowsIdx(2,k))))-baseline);
    end
    
end

fclose(hFile);

%% reshape data into square format that can be used with imagesc
X = rot90(reshape(logcoordinates(:,1), [mapwidth, mapheight]), -1);
Y = rot90(reshape(logcoordinates(:,2), [mapwidth, mapheight]), -1);
Z = zeros(M, size(X,1), size(X,2));
for k=1:M
    Z(k,:,:) = rot90(reshape(peakHeights(:,k), [mapwidth, mapheight]), -1);
end

% rescale x,y data
X = X-X(end,end);
Y = Y-Y(end,end);

% flip rows from serpentine pattern
for i = 1:size(X,1)
    if X(i,1) > X(i,end)
        X(i,:) = fliplr(X(i,:));
        Y(i,:) = fliplr(Y(i,:));
        Z(:,i,:) = fliplr(squeeze(Z(:,i,:)));
    end
end

%% style definitions
colorMap = getColormap(12);
colorMapFade = zeros(size(colorMap,1), 128, 3);
for i=1:size(colorMap,1)
    mm = max(colorMap(i,:));
    colorMapFade(i,:,:) = min(1,linspace(0,0.2,128)' + sqrt(linspace(0, 1/mm.^2, 128)').*colorMap(i,:));
end
fontSizeL = 16;
fontSizeS = 14;

%% plotting
figure('Position', [20,20,700,600]);

marg_h = 0.05;
marg_w = 0.06;
gap = [0.05, 0.03];

i = 1;
hax(i) = subtightplot(2,2,2, gap, marg_h, marg_w);
imagesc(squeeze(Z(2, :, :)));
colormap(hax(i), squeeze(colorMapFade(1,:,:)));
title(hax(i), '[Bradykinin+H]^+ 1061 m/z', 'FontSize', fontSizeL); %%23 m/z
set(hax(i), 'yticklabel',[]);
set(hax(i), 'xticklabel',[]);
hc(i) = colorbar('FontSize', fontSizeS);
axis equal tight
ylabel(hc(i), 'signal average (a. u.)', 'FontSize', fontSizeL+2);

i = 2;
hax(i) = subtightplot(2,2,4, gap, marg_h, marg_w);
imagesc(squeeze(Z(4, :, :)));
colormap(hax(i), squeeze(colorMapFade(2,:,:)));
title(hax(i), '[Bradykinin+H]^+ 1061 m/z', 'FontSize', fontSizeL); %%
set(hax(i), 'yticklabel',[]);
set(hax(i), 'xticklabel',[]);
hc(i) = colorbar('FontSize', fontSizeS);
axis equal tight
ylabel(hc(i), 'signal average (a. u.)', 'FontSize', fontSizeL+2);

i = 3;
hax(i) = subtightplot(2,2,1, gap, marg_h, marg_w);
imagesc(squeeze(Z(1, :, :)));
colormap(hax(i), squeeze(colorMapFade(3,:,:)));
title(hax(i), '[Bradykinin+H]^+ 1061 m/z','FontSize', fontSizeL); %%[H]^+
set(hax(i), 'yticklabel',[]);
set(hax(i), 'xticklabel',[]);
hc(i) = colorbar('FontSize', fontSizeS);
axis equal tight

i = 4;
hax(i) = subtightplot(2,2,3, gap, marg_h, marg_w);
imagesc(squeeze(Z(3, :, :)));
colormap(hax(i), squeeze(colorMapFade(4,:,:)));
title(hax(i), '[Bradykinin+H]^+ 1061 m/z', 'FontSize', fontSizeL); %%
set(hax(i), 'yticklabel',[]);
set(hax(i), 'xticklabel',[]);
hc(i) = colorbar('FontSize', fontSizeS);
axis equal tight

hax(5) = axes('Position', get(hax(4), 'Position'));
%ht = text(hax(5), 7, 0, '200 µm/pixel', 'FontSize', fontSizeL);
hold on
hl = plot(hax(5), [1, 6], [0, 0], 'LineWidth', 10, 'Color', [0 0 1]);
hold off
axis(hax(5), [get(hax(4), 'XLim'), get(hax(4), 'YLim')]);
set(hax(5), 'Clipping', 'off');
set(hax(5),'Visible','off');


% addpath('\\win.desy.de\group\cfel\4all\mpsd_drive\massspec\Users\Frederik\MATLAB Libraries\export_fig\');
% export_fig(gcf, [filename(1:end-length('.signal.div')), ' map'], '-jpg', '-transparent', '-m2');