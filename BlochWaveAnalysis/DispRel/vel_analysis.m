%% ACOUSTIC Simulations GtoX on Euler
clc; clear; close all
% Read data
veldata = dlmread('BC20_cornerResonator_vel.csv',',',1,0); % read in data from csv file
tstr = 'BC20 Corner Resonator M'; % title name for plots across all masses
pfacs = 0; % 1 to plot participation factors, 0 to not plot
% sort full data matrix from csv into separate cells for each mass
id = veldata(:,2) ; % column index = 1 for center/corner only, 2 for resonators
[c,ia,ib] = unique(id); 
NRad = length(c) ;
matcell = cell(NRad,1);
for i = 1:NRad
  matcell{i} = veldata(ib==i,:); 
end 

% initialize cells for real and imaginary data
matdata = {}; 
imagmatdata = {};

% assemble one data structure with everything
for l = 1:NRad
    alldata = [];
    ccell = matcell{l};
     
    % store fundamental values 
    alldata(:,1) = ccell(:,1); % wave vector: t value, column index = 2 for center/corner only, 1 for resonators
    alldata(:,2) = ccell(:,3)*1e-6; % frequency f in MHz
%     solidVolume = ccell(1,8); % volume in m3
    
    % calculate participation metric and store 
    % (vol int of velocity in x/y/z dir) / (vol int of velocity magnitude) * 1/2
    % [m^4/s]/[m^4/s]*[1]
    pmetricx = (1/2)*abs(ccell(:,4))./abs(ccell(:,7));
    pmetricy = (1/2)*abs(ccell(:,5))./abs(ccell(:,7));
    pmetricz = (1/2)*abs(ccell(:,6))./abs(ccell(:,7));
    alldata(:,3) = pmetricx; 
    alldata(:,4) = pmetricy; 
    alldata(:,5) = pmetricz;
    
    % take real and imaginary parts of data matrix and store in cells
    imagdata = imag(alldata);
    alldata = real(alldata);
    matdata{1,l} = alldata; 
    imagmatdata{1,l} = imagdata(:,2); 
end 

pfacstrings = {'X' 'Y' 'Z'}; % ends of labels for participation factor color bars 
colors = {[1 0 0],[0 1 0],[0 0 1],[0 0 0]}; % colors for participation factor figures for each possible mass (r,g,b,black)
bmat = []; % initialize cell for border points matrix 

for l = 1:NRad
    
    % FIND BANDGAPS 
    master_arr = [matdata{l}(:,1), matdata{l}(:,2)]; % assemble a matrix with t and freq 
    master_arr = sortrows(master_arr); % sort the ef_vals in order from small to big
    sort_arr = [master_arr(:,1), master_arr(:,2)]; % this is the arr with sorted items
    % making them the sorted ones
    t = sort_arr(:,1); % t (wavevector) 
    ef = sort_arr(:,2); % frequencies 

    ef_uni = unique(ef); % removing duplicate values 
    bandgaps = []; % holds all the ef_value ranges where there is a bandgap
    index = 1; % counter index for bandgaps matrix
    BGW = 120000; % width of smallest bandgap in Hz (identify by eye first) 
    for j = 1:size(ef_uni)-1
        if abs(ef_uni(j+1)-ef_uni(j)) > BGW*1e-6 % eyeballing, min range where there is a bandgap
            bandgaps(index,1) = ef_uni(j); % lower bound of bandgap
            bandgaps(index,2) = ef_uni(j+1); % upper bound of bandgap
            index = index+1;
        end
    end

    % x values for graphing fill box 
    min_t = min(t); % min t values
    max_t = max(t); % max t values 
    x = [min_t:max_t, fliplr(min_t:max_t)];
%     bandgaps = bandgaps(end-3:end,:) % trim  off extra bandgaps
    [xlen, ylen] = size(bandgaps); % find size of bandgap matrix 
    
    % bandgap metric compiler: 
    % collects and stores participation factor data for points that border bandgaps 
    bgvec = [];
    for n = 1:xlen
    bgvec = [bgvec bandgaps(n,:)];
    end 
    for i = 1:length(bgvec)
        [row, col] = find(bgvec(i)==matdata{1,l}(:,2)); % find row indices of bandgap values
        row = unique(row); % remove duplicate row indices
        border = matdata{1,l}(row,:); % store border points
        if pfacs 
            vals = border(:,3:5); % store participation factor values
    %         border(isinf(border)) = max(vals(~isinf(vals)));
            border(isinf(border)) = 0; % set inf values to 0
            for s = 1:length(row)
                normb = norm([(border(s,3)) (border(s,4)) (border(s,5))]); % normalize participation factor values
                bavg = [l bgvec(i) (border(s,3))./normb (border(s,4))./normb (border(s,5))./normb]; % store normalized values in matrix
                bmat = [bmat; bavg]; 
            end 
        end 
        
    end 
%     figure 
%     scatter(matdata{1,l}(:,1),imagmatdata{l})
    
    % plot original disp rel without color 
    figure
    set(gcf,'Position',[750,0,500*4,831]) % set figure size and location
    subplot(1,4,1); 
    hold on
    scatter(matdata{1,l}(:,1),matdata{1,l}(:,2),13,'k','filled')
    xlabel('Reduced Wavevector')
    ylabel('Frequency (MHz)')
    % xlim([0,1])
    xticks([0,1,2,3,4,5,6,7,8]) 
    xticklabels({'\Gamma','X','R','A','Z','\Gamma','M','X','\Gamma'})
    set(gca,'FontName','Arial','fontsize',15,'LineWidth',2)
    mass = string(l);
    title(strcat(tstr, mass))
    
    % plot bandgaps on disp rel 
    for i = 1:xlen
        yline(bandgaps(i,1)) 
        yline(bandgaps(i,2))
        X = ones(min_t+1, max_t+1);
        y1 = bandgaps(i,1)*X;
        y2 = bandgaps(i,2)*X;
        y = [y1, fliplr(y2)];
        fill(x,y,'k','FaceAlpha',0.2)
    end
    hold off
    
    % then plot same disp rel with participation factors in x, y, z dirs 
    for i = 1:3
        subplot(1,4,i+1);
        hold on 
        % maxPointSize = max(abs(pmetric)); 
        %use this for size-dependent participation
        % (pmetric./maxPointSize)*100

        scatter(matdata{1,l}(:,1),matdata{1,l}(:,2),...
                13,...%point size, biggest will be 5
                matdata{1,l}(:,i+2),'filled')% color values, 1 or -1 basically
        red = [1, 0, 0];
        green = [0, 1, 0];
        pink = [255, 192, 203]/255;
        lblue = [0.5, 0.5, 1]; 
        dblue = [0,   0,   1];

        % set colors
        dcol = dblue; 
        lcol = colors{l}; 
        white = [1 1 1];
        % white = dblue;
    %     colors_p = [linspace(dcol(1),white(1))', linspace(dcol(2),white(2))', linspace(dcol(3),white(3))';
    %         linspace(white(1),lcol(1))', linspace(white(2),lcol(2))', linspace(white(3),lcol(3))'];

        colors_p = [linspace(white(1),lcol(1))', linspace(white(2),lcol(2))', linspace(white(3),lcol(3))'];

        colormap(colors_p);
        caxis([0 1e-1])
        c = colorbar;
        ylabel(c,[pfacstrings{i} ' Participation Factor'])
        xlabel('Reduced Wavevector')
        ylabel('Frequency (MHz)')
        % xlim([0,1])
        xticks([0,1,2,3,4,5,6,7,8])
        xticklabels({'\Gamma','X','R','A','Z','\Gamma','M','X','\Gamma'})
        set(gca,'FontName','Arial','fontsize',15,'LineWidth',2)
        mass = string(l);
        title(strcat(tstr, mass))
        
        % plot bandgaps on disp rel
        for k = 1:xlen
            yline(bandgaps(k,1)) 
            yline(bandgaps(k,2))
            X = ones(min_t+1, max_t+1);
            y1 = bandgaps(k,1)*X;
            y2 = bandgaps(k,2)*X;
            y = [y1, fliplr(y2)];
            fill(x,y,'k','FaceAlpha',0.2)
        end
        hold off
    end 
%     saveas(gcf,strcat(filebasename, mass, '.png'))
end 


%% plots participation factor data for points that border bandgaps 
figure
BCstl = stlread('BC30/BC30.stl');
trimesh(BCstl,'FaceColor','k','FaceAlpha',0.1,'EdgeColor','k','EdgeAlpha',0.1)

hold on 
fac = 100; 
lw = 2.5; 
for i = 1:length(bmat)
    plot3([0 bmat(i,3)]*fac,[0 bmat(i,4)]*fac,[0 bmat(i,5)]*fac,'Color',colors{bmat(i,1)},'LineWidth',lw)
end 
maxlim = 60; 
xlim([0,maxlim])
ylim([0,maxlim])
zlim([0,maxlim])
xlabel('x')
ylabel('y')
zlabel('z')
title(strcat(tstr,'ass'))

ncolors = length(unique(bmat(:,1)));
h = zeros(ncolors, 1);
for j = 1:ncolors
h(j) = plot(NaN,NaN, '-', 'LineWidth', lw, 'Color', colors{j});
end 
% legend(h, 'Mass 3');
legend(h, 'Mass 1','Mass 2','Mass 3','Mass 4');

view(30,20)

%% Plot Dispersion Relation ONLY 

for l = 1:3 % l index controls which mass(es) you plot 
    
    % FIND BANDGAPS 
    master_arr = [matdata{l}(:,1), matdata{l}(:,2)]; % assemble a matrix with t and freq 
    master_arr = sortrows(master_arr); % sort the ef_vals in order from small to big
    sort_arr = [master_arr(:,1), master_arr(:,2)]; % this is the arr with sorted items
    % making them the sorted ones
    t = sort_arr(:,1); % t (wavevector) 
    ef = sort_arr(:,2); % frequencies 

    ef_uni = unique(ef); % removing duplicate values 
    bandgaps = []; % holds all the ef_value ranges where there is a bandgap
    index = 1;
    BGW = 120000; % width of smallest bandgap (identify by eye first) 
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
%     bandgaps = bandgaps(end-3:end,:) % trim  off extra bandgaps
    [xlen, ylen] = size(bandgaps); % find size of bandgap matrix 
    
    
    % calculate phase and group velocity by sorting eigenvalues into longi, shear, etc. bands 
    gvels = []; % initialize group velocity matrix 
    NEigenfreqs = 30; % number of eigenvalues per t
    for b = 1:NEigenfreqs
        % take every 30th value (since there are 30 eigenvalues per t)
        arrpart = [matdata{l}(b:30:end,1), matdata{l}(b:30:end,2)]; % assemble a matrix with t and freq
        ef = arrpart(:,2); % frequency values 
        kvec = arrpart(:,1); % reduced wavevector values 
        phvel = ef./kvec; % calculate phase velocity w/k
        grvel = gradient(ef); % calculate group velocity dw/dk 
        grvelsign = sign(grvel); % take sign of group velocity 
        gvels = [gvels; kvec ef grvel grvelsign phvel]; % compile into matrix 
    end 
    

    % plot original disp rel without color 
    figure
    set(gcf,'Position',[750,0,500,831])
%     subplot(1,4,1);
    hold on
    scatter(gvels(:,1),gvels(:,2),13,gvels(:,3),'filled') % plot group velocities 
%     scatter(matdata{1,l}(:,1),matdata{1,l}(:,2),13,'k','filled') % plot plain disp rel w/o group vel
    xlabel('Reduced Wavevector')
    ylabel('Frequency (MHz)')
    xticks([0,1,2,3,4,5,6,7,8])
    xticklabels({'\Gamma','X','R','A','Z','\Gamma','M','X','\Gamma'})
    set(gca,'FontName','Arial','fontsize',15,'LineWidth',2)
%     ylim([0,16])
    mass = string(l);
%     mass = '1';
    title(strcat(tstr, mass))
    
    % plot bandgaps on disp rel
    for i = 1:xlen
        yline(bandgaps(i,1)) 
        yline(bandgaps(i,2))
        X = ones(min_t+1, max_t+1);
        y1 = bandgaps(i,1)*X;
        y2 = bandgaps(i,2)*X;
        y = [y1, fliplr(y2)];
        fill(x,y,'k','FaceAlpha',0.2)
    end
    hold off
    
    % color bar values for group velocity plot 
    dcol = dblue; 
    lcol = red; 
    white = [1 1 1];
    colors_p = [linspace(dcol(1),white(1))', linspace(dcol(2),white(2))', linspace(dcol(3),white(3))';
        linspace(white(1),lcol(1))', linspace(white(2),lcol(2))', linspace(white(3),lcol(3))'];
    colormap(colors_p);
    c = colorbar;
%     caxis([-0.05 0.05])
    ylabel(c,'Group Velocity')
end 
