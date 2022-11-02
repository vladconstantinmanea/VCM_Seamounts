close all;
clear all;

fid = fopen('pde2d.rdm');
% Read solution from pde2d.rdm
NSAVE = fscanf(fid,'%g',1);
NEQN = fscanf(fid,'%g',1);
NX = fscanf(fid,'%g',1);
NY = fscanf(fid,'%g',1);
k = 0;
if (NX<1 || NY<1)
   error('NX or NY < 1')
end

L0 = 0;
T = zeros(NSAVE+1,1);
X = zeros(NY+1,NX+1);
Y = zeros(NY+1,NX+1);
U = zeros(NY+1,NX+1,NSAVE+1,3,NEQN);
HF = zeros(NY+1,NX+1,NSAVE+1,3,NEQN);
X_PDE2_HF = zeros(NX+1,NSAVE+1);
PDE2D_HF = zeros(NX+1,NSAVE+1);

for i=0:NX
for j=0:NY
   X(j+1,i+1) = fscanf(fid,'%g',1);
   Y(j+1,i+1) = fscanf(fid,'%g',1);     
end
end

% Read measured heat flow data
filename_1 = 'GrizzlyBear_Heat_Flow_data_Left_side.dat';
fid_1 = fopen(filename_1);
hf_data_left = textscan(fid_1,'%f %f','HeaderLines',1) ;
hf_data_left = cell2mat(hf_data_left);
X_HF_left=hf_data_left(:,1);
OBSERVED_HF_left=hf_data_left(:,2);
% get the X data bounds
[ min_X_HF_left , max_X_HF_right ] = bounds( X_HF_left , 'all' );
fclose(fid_1);

% Read measured heat flow data
filename_2 = 'GrizzlyBear_Heat_Flow_data_Right_side.dat';
fid_2 = fopen(filename_2);
hf_data_right = textscan(fid_2,'%f %f','HeaderLines',1) ;
hf_data_right = cell2mat(hf_data_right);
X_HF_right=hf_data_right(:,1);
OBSERVED_HF_right=hf_data_right(:,2);
fclose(fid_2);

% Read PDE2D multiple output heat flow files
for L=1:NSAVE
    filename_3 = ['heat_flow_PDE2D_' num2str(L)  '.out'];
    fid_3 = fopen(filename_3);
    hf_PDE2D = textscan(fid_3,'%f %f %f %f','HeaderLines',7) ;
    hf_PDE2D = cell2mat(hf_PDE2D);
    X_PDE2_HF(:,L) = hf_PDE2D(:,1);
    % heat flow in mW/m2
    PDE2D_HF(:,L) = -1000 * hf_PDE2D(:,3);
    fclose(fid_3);
end

% Read model limits 
% Surface
filename_4 = 'Grizzly_Bear_surface.dat';
fid_4 = fopen(filename_4);
XY_top = textscan(fid_4,'%f %f','HeaderLines',1) ;
XY_top = cell2mat(XY_top);
X_top=XY_top(:,1);
Y_top=XY_top(:,2);
fclose(fid_4);
% Sediments-permeable basalt limit
filename_5 = 'Grizzly_Bear_sediments_basalt_limit.dat';
fid_5 = fopen(filename_5);
XY_sed_bas = textscan(fid_5,'%f %f','HeaderLines',1) ;
XY_sed_bas = cell2mat(XY_sed_bas);
X_sed_bas=XY_sed_bas(:,1);
Y_sed_bas=XY_sed_bas(:,2);
fclose(fid_5);
% Top of the basement
filename_6 = 'Grizzly_Bear_top_basement.dat';
fid_6 = fopen(filename_6);
XY_bottom = textscan(fid_6,'%f %f','HeaderLines',1) ;
XY_bottom = cell2mat(XY_bottom);
X_bottom=XY_bottom(:,1);
Y_bottom=XY_bottom(:,2);
fclose(fid_6);


for L=L0:NSAVE
T(L+1) = fscanf(fid,'%g',1);
for ieq=1:NEQN
for ider=1:3
for i=0:NX
for j=0:NY
    
   U(j+1,i+1,L+1,ider,ieq) = fscanf(fid,'%g',1);
      
end
end
end
end
end
xmin = min(min(X(:,:)));
xmax = max(max(X(:,:)));
ymin = min(min(Y(:,:)));
ymax = max(max(Y(:,:)));
hx = 0.1*(xmax-xmin);
hy = 0.1*(ymax-ymin);


% ******* Choose variables for scalar plots
            for IEQ=1:NEQN
% Plot U_IEQ, if IDER=1
%      A_IEQ, if IDER=2
%      B_IEQ, if IDER=3
IDER = 1;
umin = min(min(min(U(:,:,L0+1:NSAVE+1,IDER,IEQ))));
umax = max(max(max(U(:,:,L0+1:NSAVE+1,IDER,IEQ))));
if (umax==umin); umax = umin+1; end

%Polt pressure
if(IEQ == 1)
for L=L0:NSAVE
   h1=figure;
   hold on;
% Pressure in MPa
   U(:,:,L+1,IDER,IEQ) = U(:,:,L+1,IDER,IEQ)/1E6;
   pcolor(X,Y,U(:,:,L+1,IDER,IEQ))
   [c, h] = contour(X,Y,U(:,:,L+1,IDER,IEQ),'LineColor','black');
   clabel(c)
   shading interp;
   axis equal tight;
   c = colorbar;
   xlabel('Distance X (m)')
   ylabel('Depth Y (m)')
    title(['Pressure field, time = ', num2str(T(L+1)/(365*24*3600),'%.2f'),' years'], ...
        'Units', 'normalized','Position', [0.5, 1.1, 1.0]) 
    c.Label.String = 'P (MPa)';
    colormap(flipud(copper));
    
   % Plot model surfaces
   plot(X_top,Y_top, 'Color', 'k','LineWidth',1.5);
   hold on;
   plot(X_sed_bas,Y_sed_bas, '--', 'Color', 'w','LineWidth',1.5); 
   hold on;
   plot(X_bottom,Y_bottom, '--', 'Color', 'w','LineWidth',1.5); 
   hold off;
   % Save plot for each NSAVE 
   fig_name_1 = strcat('Pressure_PDE2D_',num2str(L),'.png');
   print(h1,fig_name_1,'-dpng','-r300'); 
end
end

% Plot temperature & heat flow (observed and from PDE2D)
if(IEQ == 2)
for L=L0:NSAVE
   h2=figure;
   subplot(2,1,1)
% align the two plots
   tmp=get(gca,'position');
   set(gca,'position',[tmp(1) tmp(2) 0.87*tmp(3) tmp(4)])
% plot the left side observed heat flow data
   plot(X_HF_left,OBSERVED_HF_left,'.','markersize', 8, 'Color', 'blue');
   hold on;
   % plot the right side observed heat flow data
   plot(X_HF_right,OBSERVED_HF_right,'.','markersize', 8, 'Color', 'red');
   axis([0 xmax 0 500]);
   hold on;
   if(L>0)   
       
% plot PDE2D output heat flow whithin the observed X bounds 
   tolerance = 5 * eps(1E12);
   index_1 = find(ismembertol(X_PDE2_HF(:,1), min_X_HF_left, tolerance));
   index_2 = find(ismembertol(X_PDE2_HF(:,1), max_X_HF_right, tolerance));   
   plot(X_PDE2_HF(index_1+1:index_2,L), PDE2D_HF(index_1+1:index_2,L), ...
        'Color', 'black','LineWidth',1.1);
% plot PDE2D output heat flow along the entire profile
%        plot(X_PDE2_HF(:,L),PDE2D_HF(:,L), 'Color', 'black',...
%            'LineWidth',1.1);

        xlabel('Distance X (m)')
        ylabel('Heat flow (mW/m^2)')
        title(['Heat flow, time = ',num2str(T(L+1)/(365*24*3600),'%.2f'),...
            ' years'], 'Units', 'normalized','Position', [0.5, 1.05, 1.0])
   end
   legend('Observed heat flow (left side)' , ...
       'Observed heat flow (right side)', ...
       'Calculated heat flow (PDE2D)')
   hold on;
   subplot(2,1,2)   
% plot the temperature (in degC)
   pcolor(X,Y,U(:,:,L+1,IDER,IEQ)-273)
   hold on;
% plot temperature contours (in degC)
   [c, h] = contour(X,Y,U(:,:,L+1,IDER,IEQ)-273,'LineColor','black');
   clabel(c)
   shading interp;
   axis([0 xmax 0 (ymax+100)]);
%   axis equal tight;
   c = colorbar;
   xlabel('Distance X (m)')
   ylabel('Depth Y (m)')
   title(['Temperature field, time = ',num2str(T(L+1)/(365*24*3600),'%.2f'),' years'], ...
        'Units', 'normalized','Position', [0.55, 1.05, 1.0])
    c.Label.String = 'T (\circC)'; 
    colormap(turbo);
    brighten(turbo,.5) 
   % Plot model surfaces
   plot(X_top,Y_top, 'Color', 'k','LineWidth',1.5);
   hold on;
   plot(X_sed_bas,Y_sed_bas,'--','Color', 'k','LineWidth',1.5); 
   hold on;
   plot(X_bottom,Y_bottom,'--','Color', 'k','LineWidth',1.5); 
% Save plot for each NSAVE 
   fig_name_2 = strcat('Temperature_PDE2D_',num2str(L),'.png');
   print(h2,fig_name_2,'-dpng','-r300'); 
   hold off;
end
end
            end
            
% ******* Choose variables for vector plots
% IEQi,IDERi indicates:
%      U_IEQi, if IDERi=1
%      A_IEQi, if IDERi=2
%      B_IEQi, if IDERi=3
IEQ1 = 1;
IDER1 = 2;
IEQ2 = 1;
IDER2 = 3;
for L=L0:NSAVE
   Umax = max(max(abs(U(:,:,L+1,IDER1,IEQ1))));
   Vmax = max(max(abs(U(:,:,L+1,IDER2,IEQ2))));
   h3=figure;
   % Plot model surfaces
   plot(X_top,Y_top, 'Color', 'k','LineWidth',1.5);
   hold on;
   plot(X_sed_bas,Y_sed_bas,'--','Color', 'k','LineWidth',1.5); 
   hold on;
   plot(X_bottom,Y_bottom,'--', 'Color', 'k','LineWidth',1.5); 
   hold on;
% For plotting only velocity vectors replace streamslice with quiver
   density = 5;
   l = streamslice(X,Y,U(:,:,L+1,IDER1,IEQ1), ...
              U(:,:,L+1,IDER2,IEQ2),density);
   set(l, 'LineWidth' , 0.5)
   set(l, 'Color' , 'k' ); 
   axis equal tight;   
   axis([0 xmax 0 (ymax+100)]);
   xlabel('Distance X (m)')
   ylabel('Depth Y (m)')
    title(['Streamlines, time = ',num2str(T(L+1)/(365*24*3600),'%.2f'),' years'], ...
        'Units', 'normalized','Position', [0.5, 1.05, 1.0])
% Save plot for each NSAVE
   fig_name_3 = strcat('Streamlines_PDE2D_',num2str(L),'.png');
   print(h3,fig_name_3,'-dpng','-r300'); 
   hold off;
end
