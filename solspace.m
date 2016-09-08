% Function that calculates space of solutions via det(M). The inputs xrange
% and yrange are the ranges of the two components of the response
% coefficient to be analyzed. The output is a list with the xy coordinates
% of the minima. The function also saves a figure with the det-landscape
% and stores the xyz as well as seeds on a .dat file.

function space = solspace(xrange,yrange)
% We load the parameters
global motor chib bc Sp kp xihead k xi

%% CALCULATE DETERMINANT SURFACE
% Depending on the motor model we generate xyz-data and calculate minima  
switch motor
    case 'sliding'
        % Determinant of bvp for sliding control 
        determinant = @(res) det(bcmat([res(1),res(2),0,0]));
        
        % We compute the determinant on the given range
        detdat = zeros(length(xrange),length(yrange));
        for i=1:length(xrange)
            for j=1:length(yrange)
                detdat(i,j)=abs(determinant([xrange(i),yrange(j)]));
            end
        end
        
        % Data is stored in log-scale as xyz
        [x,y]   = meshgrid(xrange,yrange);
        stringx = 'Sliding response \chi''';
        stringy = 'Sliding response \chi''''';
        z       = log10(detdat);
        
        % And calculate the minima
        [~,~,~,imin] = extrema2(z');
        
    case 'dyn-curvature'
        % Determinant of bvp for sliding control 
        determinant = @(res) det(bcmat([res(1),0,0,res(2)]));
        
        % We compute the determinant on the given range
        detdat = zeros(length(xrange),length(yrange));
        for i=1:length(xrange)
            for j=1:length(yrange)
                detdat(i,j)=abs(determinant([xrange(i),yrange(j)]));
            end
        end
        
        % Data is stored in log-scale as xyz
        [x,y]   = meshgrid(xrange,yrange);
        stringx = 'Sliding response \chi''';
        stringy = 'Curvature response \beta''''';
        z       = log10(detdat);
        
        % And calculate the minima
        [~,~,~,imin] = extrema2(z'); 
        
    case 'curvature'
        % Determinant of bvp for sliding control 
        determinant = @(res) det(bcmat([k,xi,res(1),res(2)]));
        
        % We compute the determinant on the given range
        detdat = zeros(length(xrange),length(yrange));
        for i=1:length(xrange)
            for j=1:length(yrange)
                detdat(i,j)=abs(determinant([xrange(i),yrange(j)]));
            end
        end
        
        % Data is stored in log-scale as xyz
        [x,y]   = meshgrid(xrange,yrange);
        stringx = 'Curvature response \beta''';
        stringy = 'Curvature response \beta''''';
        z       = log10(detdat);
        
        % And calculate the minima
        [~,~,~,imin] = extrema2(z');              
end

%% OBTAIN SEEDS
% We find the points that are xy null
xmin=x(imin);
ymin=y(imin);
nonull=((xmin~=0)+(ymin~=0))~=0;

% And select the rest
xmin = xmin(nonull);
ymin = ymin(nonull);
space.seeds = [xmin';ymin'];
space.xyz   = {x,y,z};

%% PLOT HEATMAP
% Generate heat map
figure;
colormap(winter);
%imagesc([xrange(1),xrange(end)],[yrange(1),yrange(end)],z);
contourf(x,y,z');
daspect([1 1 1]);
xlabel(stringx,'FontSize',12,'FontName','Helvetica');
ylabel(stringy,'FontSize',12,'FontName','Helvetica');
h = colorbar;
ylabel(h, 'log_{10}(det(M))');
hold on

% Add dots and numbers for modes
plot(xmin,ymin,'r.')
[~,indexseeds]=sort(xmin.^2+ymin.^2);
text(xmin(indexseeds)-0*abs(xrange(end)-xrange(1)),...
     ymin(indexseeds)-0*abs(xrange(end)-xrange(1)),...
     strread(num2str(1:length(indexseeds)),'%s'),'Color',[1 1 1]);

% Anotate the plot
inset={['Motor model: ' motor],['Boundary conditions: ' bc],...
       ['S_p=' num2str(Sp) ';   \chi_b=' num2str(chib)],...
       ['k_p=' num2str(kp) ';   \xi_h=' num2str(xihead)], ...
       ['k_i=' num2str(k)  ';   \xi_i=' num2str(xi)]};

annotation('textbox',...
    [.15 .7 .33 .17],...
    'String',inset,...
    'FontSize',12,...
    'FontName','Helvetica',...
    'EdgeColor',[1 1 1],...
    'FitBoxToText','off',...
    'Color',[1 1 1],...
    'FitHeightToText','on');
hold off