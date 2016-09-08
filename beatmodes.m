% Function that calculates the beating modes given a list of seeds.

function sol = beatmodes(seeds)
%Load parameters
global Sp motor k xi D0;

%% STRUCTURE
nseeds = size(seeds,2);
sol(nseeds).res = [0,0];
sol(nseeds).k   = [0,0,0,0];
sol(nseeds).A   = [0,0,0,0,0];
sol(nseeds).err =  0;
titles = cell(1,nseeds);

%% CALCULATE SOLUTIONS
switch motor
    case 'sliding'
        % Reverse condition number for sliding control
        revcond = @(res) rcond(bcmat([res(1),res(2),0,0]));
        
        % Set options for the solver
        options = optimset('Display','off','TolFun',1e-15,'TolX',1e-15);
        
        % Loop over all seeds
        for i=1:nseeds
            % Solve the singular system using singular value decomposition
            x0         = [seeds(1,i),seeds(2,i)];
            sol(i).res = fmincon(revcond,x0,[[1,0];[0,1]],[0 0],[],[],...
                               [],[],[],options);
            M          = bcmat([sol(i).res(1),sol(i).res(2),0,0]);
            [~,~,V]    = svd(M,'econ');
            
            % Assign amplitudes, modes, and errors
            sol(i).A = V(:,size(M,2));
            sol(i).k = roots([1,0,-(sol(i).res(1)+1i*sol(i).res(2)),0,...
                1i*Sp]);
            sol(i).err = sum(abs(M*sol(i).A));
            
            % Create plot title
            titles{i}=['Mode \chi='...
                        num2str(sol(i).res(1)+1i*sol(i).res(2))...
                       '  with error   '  num2str(sol(i).err)];
        end
        
    case 'dyn-curvature'
        % Reverse condition number for dynamic curvature control
        revcond = @(res) rcond(bcmat([res(1),0,0,res(2)]));
        
        % Set options for the solver
        options = optimset('Display','off','TolFun',1e-15,'TolX',1e-15);

        % Loop over all seeds
        for i=1:nseeds
            % Solve the singular system using singular value decomposition
            x0       = [seeds(1,i),seeds(2,i)];
            sol(i).res=fmincon(revcond,x0,[-1,0],0,[],[],[],[],[],options);
            M        = bcmat([sol(i).res(1),0,0,sol(i).res(2)]);
            [~,~,V]  = svd(M,'econ');
            
            % Assign amplitudes, modes, and errors
            sol(i).A = V(:,size(M,2));
            sol(i).k =roots([1,-1i*sol(i).res(2),-sol(i).res(1),0,1i*Sp]);
            sol(i).err = sum(abs(M*sol(i).A));
            
            % Create plot title
            titles{i}=['Mode \chi''=' num2str(sol(i).res(1))...
                       '     and \beta''''=' num2str(1i*sol(i).res(2))...
                       '  with error   '  num2str(sol(i).err)];
        end
        
    case 'curvature'
        % Reverse condition number for dynamic curvature control
        revcond = @(res) rcond(bcmat([k,xi,res(1),res(2)]));
        
        % Set options for the solver
        options = optimset('Display','off','TolFun',1e-15,'TolX',1e-15);
        
        % Loop over all seeds
        for i=1:nseeds
            % Solve the singular system using singular value decomposition
            x0       = [seeds(1,i),seeds(2,i)];
            sol(i).res=fmincon(revcond,x0,[-1,0],10^7,[],[],[],[],[],options);
            M        = bcmat([k,xi,sol(i).res(1),sol(i).res(2)]);
            [~,~,V]  = svd(M,'econ');
            
            % Assign amplitudes, modes, and errors
            sol(i).A = V(:,size(M,2));
            sol(i).k =roots([1,-(sol(i).res(1)+1i*sol(i).res(2)),...
                      -(k+xi),0,1i*Sp]);
            sol(i).err = sum(abs(M*sol(i).A));
            
            % Create plot title
            titles{i}=['Mode \beta='...
                        num2str(sol(i).res(1)+1i*sol(i).res(2))...
                       '  with error   '  num2str(sol(i).err)];
        end
end

%% PLOTS
for i=1:nseeds
    % Sample a solution in arc-length "s" and time "t"
    ds  =0.01; dt = 0.1;
    s   = 0:ds:1;
    t   = 0:dt:1;
    psi = (sol(i).A(1)*exp(sol(i).k(1)*s)+sol(i).A(2)*exp(sol(i).k(2)*s)...
        +  sol(i).A(3)*exp(sol(i).k(3)*s)+sol(i).A(4)*exp(sol(i).k(4)*s));
    psi = psi/(max([0.01 abs(psi(1)) abs(psi(end))]))*0.5;
    
    % Plot amplitude, phase, real, and imaginary parts
    figure;
    subplot(3,2,1),plot(s,abs(psi),'k');
    xlabel('Arc-length, s/L'); ylabel('|\psi|');
     
    subplot(3,2,3),plot(s,unwrap(angle(psi)),'k');
    xlabel('Arc-length, s/L'); ylabel('arg(\psi)');
    
    subplot(3,2,5),plot(s,real(psi),'b',s,imag(psi),'r');
    xlabel('Arc-length, s/L');
    ylabel('{\color{blue}Re(\psi)} and {\color{red}Im(\psi)}');

    % Plot angle time-series
    subplot(3,2,2),plot(s,2*real(psi'*exp(2*pi*1i*t)));
    xlabel('Arc-length, s/L');
    ylabel('\psi(s,t)');
    
    % Plot positional representation
    subplot(3,2,[4 6])
    plot(ds*cumtrapz(cos(D0*(s'*ones(1,length(t)))+2*real(psi'*exp(2*pi*1i*t)))),...
         ds*cumtrapz(sin(D0*(s'*ones(1,length(t)))+2*real(psi'*exp(2*pi*1i*t)))));
    %axis equal 
    xlabel('x-position');
    ylabel('y-position');
    
    % Generate title
    set(gcf,'NextPlot','add');
    axes; 
    set(gca,'Visible','off'); 
    h = title(titles{i},...
        'FontWeight','b',...
        'FontSize',14,...
        'FontName','Helvetica');
    set(h,'Visible','on');
end