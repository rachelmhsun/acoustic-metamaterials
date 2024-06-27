%% ACOUSTIC Simulations 
clc; clear; 
close all

%% plot flags/options 
pfacFlag = 0; % plot participation factors ON = 1, OFF = 0
multiDispRel = 0; % multiple dispersion relations ON = 1, OFF = 0
IBZ = 'tetrag'; % for reduced wavevector plot labels 
bandgapsFlag = 1; % plot bandgaps ON = 1, OFF = 0

%% Read data
filename = '1_0d_cen_5_gap_tetragonalFourFold'; 
veldata = readmatrix([filename,'.xlsx']);
tstr = 'Adjusted Modulus BC'; % title string across all masses

%% sort data if multiple dispersion relations 
% (for multiple geometrical parameters) 
if multiDispRel == 1
    % sort full data matrix from csv into separate cells for each mass
    id = veldata(:,1) ; % column index = 1 for center/corner only, 2 for resonators
    [c,ia,ib] = unique(id); 
    NRad = length(c);
    matcell = cell(NRad,1);
    for i = 1:NRad
      matcell{i} = veldata(ib==i,:); 
    end  
    tIndex = 2; % t parameter column index 
    fIndex = 3; % frequency column index 
else 
    % only 1 dispersion relation 
    NRad = 1;
    matcell{1} = veldata; 
    tIndex = 1; 
    fIndex = 2; 
end 

%% Assemble and calculate data 
% initialize cells for real and imaginary data
matdata = {}; 
imagmatdata = {};

% assemble one data structure with everything
for NRadIndex = 1:NRad
    alldata = [];
    ccell = matcell{NRadIndex};
     
    % store fundamental values 
    alldata(:,1) = ccell(:,tIndex); % wave vector: t value, column index = 2 for center/corner only, 1 for resonators
    alldata(:,2) = ccell(:,fIndex)*1e-6; % frequency f in MHz
    if pfacFlag == 1
        % calculate participation metric and store 
        % (vol int of velocity in x/y/z dir) / (vol int of velocity magnitude) * 1/2
        % [m^4/s]/[m^4/s]*[1]
        pmetricx = (1/2)*abs(ccell(:,3))./(ccell(:,6));
        pmetricy = (1/2)*abs(ccell(:,4))./(ccell(:,6));
        pmetricz = (1/2)*abs(ccell(:,5))./(ccell(:,6));
        alldata(:,3) = pmetricx; 
        alldata(:,4) = pmetricy; 
        alldata(:,5) = pmetricz;
    end 
    % take real and imaginary parts of data matrix and store in cells
    imagdata = imag(alldata);
    alldata = real(alldata);
    matdata{1,NRadIndex} = alldata; 
    imagmatdata{1,NRadIndex} = imagdata(:,2); 
end 


%% Plot Dispersion Relation ONLY 

RDstring = {'15','20','25','30'}; 



for NRadIndex = 1:NRad % l index controls which mass(es) you plot 
    
    % FIND BANDGAPS 
    master_arr = [matdata{NRadIndex}(:,1), matdata{NRadIndex}(:,2)]; % assemble a matrix with t and freq 
    master_arr = sortrows(master_arr); % sort the ef_vals in order from small to big
    sort_arr = [master_arr(:,1), master_arr(:,2)]; % this is the arr with sorted items
    % making them the sorted ones
    t = sort_arr(:,1); % t (wavevector) 
    ef = sort_arr(:,2); % frequencies 

    ef_uni = unique(ef); % removing duplicate values 
    bandgaps = []; % holds all the ef_value ranges where there is a bandgap
    index = 1;
    BGW = 70000; % width of smallest bandgap (identify by eye first) 
    for j = 1:size(ef_uni)-1
        if abs(ef_uni(j+1)-ef_uni(j)) > BGW*1e-6 % eyeballing, min range where there is a bandgap
            bandgaps(index,1) = ef_uni(j);
            bandgaps(index,2) = ef_uni(j+1);
            index = index+1;
        end
    end

    % x values for graphing fill box 
    min_t = min(t);
    max_t = max(t);
    
    x = [min_t:max_t, fliplr(min_t:max_t)];
%     bandgaps = bandgaps(end-2:end-1,:) % trim  off extra bandgaps
    [xlen, ylen] = size(bandgaps); % find size of bandgap matrix 
    
    
    % calculate phase and group velocity
    gvels = []; % initialize group velocity matrix 
    NEfreqs = 30; % number of eigenfrequencies per t 
    for b = 1:NEfreqs
        % take every 30th value (since there are 30 eigenvalues per t)
        arrpart = [matdata{NRadIndex}(b:30:end,1), matdata{NRadIndex}(b:30:end,2)]; 
        ef = arrpart(:,2);
        kvec = arrpart(:,1);
        dt = kvec(2)-kvec(1);
        phvel = ef./kvec; % calculate phase velocity w/k
        grvel = gradient(ef,dt); % calculate group velocity dw/dk 
        grvelsign = sign(grvel); % take sign of group velocity 
        gvels = [gvels; kvec ef grvel grvelsign phvel]; % compile into matrix 
    end 
    

%% PLOT DISPERSION RELATION 
    figure
    set(gcf,'Position',[750,0,500,831])
%     subplot(1,4,1);
    hold on
    box on
%     scatter(gvels(:,1),gvels(:,2)*10^6/nfac,13,gvels(:,3),'filled') % plot group velocities with normalized frequencies 
%     scatter(gvels(:,1),gvels(:,2),13,gvels(:,3),'filled') % plot group velocities 
    scatter(matdata{1,NRadIndex}(:,1),matdata{1,NRadIndex}(:,2),13,'k','filled') % plot plain disp rel w/o group vel
    if pfacFlag == 1
        pfacstrings = {'X' 'Y' 'Z'}; % ends of labels for participation factor color bars 
        colors = {[1 0 0],[0 1 0],[0 0 1],[0 0 0]}; % colors for participation factor figures for each possible mass (r,g,b,black)

        % plot participation factor (matdata column 5 = z, 4 = y, 3 = x)
        scatter(matdata{1,NRadIndex}(:,1),matdata{1,NRadIndex}(:,2),...
                13,...%point size, biggest will be 5
                matdata{1,NRadIndex}(:,5),'filled')% color values, 1 or -1 basically, change matdata column for diff direction
        red = [1, 0, 0];
        green = [0, 1, 0];
        pink = [255, 192, 203]/255;
        lblue = [0.5, 0.5, 1]; 
        dblue = [0,   0,   1];
    
        % set colors
        dcol = dblue; 
        lcol = colors{NRadIndex}; 
        white = [0.9, 0.9, 0.9];
        % colorbar and colormap 
        colors_p = [linspace(white(1),lcol(1))', linspace(white(2),lcol(2))', linspace(white(3),lcol(3))'];
        colormap(colors_p);
        caxis([0 1e-1])
        c = colorbar;
        ylabel(c,['Z Participation Factor'])
    end 

    % figure settings 
    ylabel('Frequency (MHz)')
    % set x axis labels and ticks 
    xlabel('Reduced Wavevector')
    xticks([0,1,2,3,4,5,6,7,8])
    if strcmp(IBZ,'cubic')
        xticklabels({'\Gamma','X', 'M','\Gamma','R','M','X','R','\Gamma'}) % cubic 
    else 
        xticklabels({'\Gamma','X', 'R','A','Z','\Gamma','M','X','\Gamma'}) % tetragonal 
    end 
    set(gca,'FontName','Helvetica Neue','fontsize',15,'LineWidth',1)


    % plot bandgaps on disp rel
    if bandgapsFlag == 1
        for i = 1:xlen
            yline(bandgaps(i,1)) 
            yline(bandgaps(i,2))
            X = ones(min_t+1, max_t+1);
            % X = ones(1,9);
            y1 = bandgaps(i,1)*X;
            y2 = bandgaps(i,2)*X;
            y = [y1, fliplr(y2)];
            fill(x,y,'k','FaceAlpha',0.2)
        end
        % store bandgap information in cell
        bgvec{NRadIndex} = [bandgaps, bandgaps(:,2)-bandgaps(:,1), (bandgaps(:,2)+bandgaps(:,1))/2];
    end 
    hold off
end 
ylim([0,12])
%% save dispersion relation plot 
print(gcf,'-vector','-dsvg',['images/',filename,'.svg']) % svg
saveas(gca, ['images/',filename],'fig')
saveas(gca, ['images/',filename],'png')
