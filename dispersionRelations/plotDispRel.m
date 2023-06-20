%% ACOUSTIC Simulations 
clc; clear; close all
% Read data
veldata = readmatrix('TBD/BC6p3CorM13p8.xlsx');
tstr = 'Adjusted Modulus BC'; % title string across all masses

% sort full data matrix from csv into separate cells for each mass
id = veldata(:,1) ; % column index = 1 for center/corner only, 2 for resonators
[c,ia,ib] = unique(id); 
NRad = length(c);
NRad = 1;
matcell = cell(NRad,1);
% for i = 1:NRad
%   matcell{i} = veldata(ib==i,:); 
% end  
matcell{1} = veldata; 

% initialize cells for real and imaginary data
matdata = {}; 
imagmatdata = {};

% assemble one data structure with everything
for l = 1:NRad
    alldata = [];
    ccell = matcell{l};
     
    % store fundamental values 
    alldata(:,1) = ccell(:,1); % wave vector: t value, column index = 2 for center/corner only, 1 for resonators
    alldata(:,2) = ccell(:,2)*1e-6; % frequency f in MHz
    
    % take real and imaginary parts of data matrix and store in cells
    imagdata = imag(alldata);
    alldata = real(alldata);
    matdata{1,l} = alldata; 
    imagmatdata{1,l} = imagdata(:,2); 
end 

pfacstrings = {'X' 'Y' 'Z'}; % ends of labels for participation factor color bars 
colors = {[1 0 0],[0 1 0],[0 0 1],[0 0 0]}; % colors for participation factor figures for each possible mass (r,g,b,black)
bmat = []; % initialize cell for border points matrix 

%% Plot Dispersion Relation ONLY 

RDstring = {'15','20','25','30'}; 



for l = 1 % l index controls which mass(es) you plot 
    
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
    BGW = 100000; % width of smallest bandgap (identify by eye first) 
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
    for b = 1:30
        % take every 30th value (since there are 30 eigenvalues per t)
        arrpart = [matdata{l}(b:30:end,1), matdata{l}(b:30:end,2)]; 
        ef = arrpart(:,2);
        kvec = arrpart(:,1);
        dt = kvec(2)-kvec(1);
        phvel = ef./kvec; % calculate phase velocity w/k
        grvel = gradient(ef,dt); % calculate group velocity dw/dk 
        grvelsign = sign(grvel); % take sign of group velocity 
        gvels = [gvels; kvec ef grvel grvelsign phvel]; % compile into matrix 
    end 
    

    % PLOT DISPERSION RELATION 
    figure
    set(gcf,'Position',[750,0,500,831])
%     subplot(1,4,1);
    hold on
%     scatter(gvels(:,1),gvels(:,2)*10^6/nfac,13,gvels(:,3),'filled') % plot group velocities with normalized frequencies 
%     scatter(gvels(:,1),gvels(:,2),13,gvels(:,3),'filled') % plot group velocities 
    scatter(matdata{1,l}(:,1),matdata{1,l}(:,2),13,'k','filled') % plot plain disp rel w/o group vel
    xlabel('Reduced Wavevector')
    ylabel('Frequency (MHz)')
    xticks([0,1,2,3,4,5,6,7,8])
    xticklabels({'\Gamma','X','R','A','Z','\Gamma','M','X','\Gamma'}) % LABEL IBZ POINTS 
    set(gca,'FontName','Helvetica Neue','fontsize',15,'LineWidth',1)

    mass = RDstring{l};

    
    % plot bandgaps on disp rel
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
    hold off
    bgvec = [bandgaps, bandgaps(:,2)-bandgaps(:,1), (bandgaps(:,2)+bandgaps(:,1))/2];

end 



