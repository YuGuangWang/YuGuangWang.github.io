function fig2 = pltfunc_rbf_loc_2(F, X,tstr)
% [K, fh1, fh2, fig1,fig2] = pltfunc(F, X, tstr, ifig, scale, ipr)
% Plot the function F at the m points X on the unit sphere S^2 in R^3
% Title string tstr using figure windows ifig and ifig+1 and scaling scale
% are all optional

% Default arguments
% if nargin < 6
%     ipr = 0;
% end;
% if nargin < 5
%     scale = 0.5;
% end;
% if nargin < 4
%     ifig = 1;
% end;
if nargin < 3
    tstr = [];
end;

scale = 0.5;

% Default view
%vw = [45, 10];
% vw = [-60 25];

% Use Matlab 6 routine based on QHULL
K = convhulln(X');

X1=[X(1,K(:,1)); X(1,K(:,2)); X(1,K(:,3))];
Y1=[X(2,K(:,1)); X(2,K(:,2)); X(2,K(:,3))];
Z1=[X(3,K(:,1)); X(3,K(:,2)); X(3,K(:,3))];
C = F(K');

% if ifig > 0
%     fig1 = figure(ifig); clf;
%     fh1 = patch(X1, Y1, Z1, C);
%     colormap(jet(255));
%     axis vis3d
%     axis equal tight
%     view([160 10]);
%     grid off
%     set(gca, 'Visible', 'off')
% %     colorbar('SouthOutside');
%     cbh1 = colorbar('Location','EastOutside');
%     if ipr > 10
%         hold on
%         plot3(X(1,:), X(2,:), X(3,:), 'k.', 'MarkerSize', 16);
%         hold off
%     else
%         set(fh1, 'EdgeColor', 'none')
%     end;
%     title(tstr);
% end;

[Fmax, imax] = max(F);
%xmax = X(:,imax);
[Fmin, imin] = min(F);
%xmin = X(:,imin);
fprintf('Minimum function value = %.6f, Maximum function value = %.6f\n', Fmin, Fmax);

FS = 1 + (scale/(Fmax-Fmin))*(C-Fmin);
%set(gca,'CLim',[0 0.23]); %%%%%%
%FS = C;

fig2 = figure;
clf
fh2 = patch(X1.*FS, Y1.*FS, Z1.*FS, C);
set(fh2, 'EdgeColor', 'none');
%colormap(jet(255));
colormap(jet(1023));
%colormap(gray(1023));
% cbh = colorbar('SouthOutside');
cbh2 = colorbar('Location','EastOutside');
%set(gca,'CLim',[0 0.23]); %%%%%%
%cbh = colorbar;
cbp = get(cbh2, 'Position');
% position = [left bottom width height]
%cbp(4) = 0.5*cbp(4);
% cbp(3) = 1*cbp(3);
% set(cbh, 'Position', cbp);
%set(cbh, 'FontSize', 6);
%view(90,0);
axis vis3d
axis equal tight
%view(xmax);
    AZ=160;
    EL=10;
    view([AZ EL]) % modified view
%grid on
axis off
