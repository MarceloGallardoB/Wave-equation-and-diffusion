%% Parametres %% (A MODIFIER SELON VOS BESOINS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lw = 2; fs = 16; Nx = 64;
nom = '1_a)_2D.gif';

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

data = load("output_E.out");
t = data(:,1);

data = load("output_f.out");
f = data(:,2:end);

% size(x)
% size(y)
% size(f)
h = figure;
for i = 1:size(t)
    % Calculs de valeur
    f = data(1+(i-1)*Nx:i*Nx,2:end);
    maxf = 0;
    for j = 1:Nx
        for l = 1:Nx
           if abs(f(j,l)) > maxf
              maxf =  abs(f(j,l));
           end
        end  
    end
    c = sprintf('t = %.6g s;  max(|f|) = %.6g', t(i), maxf);
    
    % Dessin d'un pas de temps
    plot(x,f);
    ylim([-3 3])
    xlim([0 10])
    xlabel('X', 'fontsize', fs)
    ylabel('Z', 'fontsize', fs)
    title(c)

    drawnow
    
    
    % Sauvegarde la simulation en un .gif
    frame = getframe(h); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 

    if i == 1 
        imwrite(imind,cm,nom,'gif', 'Loopcount',inf); 
    else 
        imwrite(imind,cm,nom,'gif','WriteMode','append'); 
    end 
end