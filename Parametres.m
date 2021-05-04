%% Parametres %% (A MODIFIER SELON VOS BESOINS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%donnees
lw = 2; fs = 16; nx = 64;
nom = '3.3_a)_fixed.gif';

repertoire = ''; % Chemin d'acces au code compile
executable = 'Exercice7'; % Nom de l'executable
input = 'configuration.in'; % Nom du fichier d'entree

%% Simulations et analyse %%
%%%%%%%%%%%%%%%%%

cmd = sprintf('%s %s %s %s%s %s','set', 'path=%path:C:\Program Files\MATLAB\R2020b\bin\win64;=%', '&', repertoire, executable, input);
system(cmd);
disp('Done.')

data = load("output_mesh.out");
x = data(1,:);
y = data(2,:);
%% XY
data = load("output_E.out");
t = data(:,1); 
%%E(t)

data = load("output_f.out");
f = data(:,2:end);
h = figure;
for i = 1:size(t)
    %Dessin de la fonction
    f = data(1+(i-1)*nx:i*nx,2:end); 
    %%f (1+ n_x .(k-1), k . nx)
    surf(x,y,f);
    %%dimensionsd e XYZ
    xlim([0 10])
    ylim([0 6])
    zlim([-3 3])
    xlabel('X', 'fontsize', fs)
    ylabel('Y', 'fontsize', fs)
    zlabel('Z', 'fontsize', fs)

     drawnow
%      pause(10);

  frame = getframe(h); 
  im = frame2im(frame); 
  [imind,cm] = rgb2ind(im,256); 

  if i == 1 
      imwrite(imind,cm,nom,'gif', 'Loopcount',inf); 
  else 
      imwrite(imind,cm,nom,'gif','WriteMode','append'); 
  end 

end