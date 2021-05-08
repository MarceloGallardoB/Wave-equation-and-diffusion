%% Parametres %% (A MODIFIER SELON VOS BESOINS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lw = 2; fs = 16;
u = 4; nx = 64; ny = 64; Lx = 10; Ly = 6; m = 2; n = 2; hx = Lx/(nx-1); hy = Ly/(ny-1); kx = pi*m/Lx; ky = pi*n/Ly; omega = 4*pi*sqrt(m^2/Lx^2 + n^2/Ly^2);
tfin = 10;

repertoire = ''; % Chemin d'accès au code compilé
executable = 'Exercice7'; % Nom de l'exécutable
input = 'configuration.in'; % Nom du fichier d'entrée

nsimul = 10;

nper = round(logspace(2,3,nsimul)); % Nombre d'itérations entier de 10^2 à 10^4  
CFL = (63./nper);
dt = tfin./nper;

paramstr = 'CFL'; % Nom du paramètre à scanner  
param = CFL; % Valeurs du paramètre à scanner 

%% Simulations %%
%%%%%%%%%%%%%%%%%
output = cell(1, nsimul); % Tableau de cellules contenant les noms des fichiers de sortie
for i = 1:nsimul
    output{i} = [paramstr, '=', num2str(param(i)), '.out'];
    % Execution du programme en lui envoyant la valeur a scanner en argument
    para = sprintf('%s%s%s', paramstr, '=', num2str(param(i)));
    cmd = sprintf('%s %s %s %s%s %s CFL=%.15g output_file=%s', 'set', 'path=%path:C:\Program Files\MATLAB\R2020b\bin\win64;=%', '&', repertoire, executable, input, CFL(i), output{i});
    system(cmd);
    disp('Done.')
end

%% Analyse %%
%%%%%%%%%%%%%%%%%
data = load("output_mesh.out");
x = data(1,:);
y = data(2,:);

fana = zeros(nx,ny);
for j = 1:nx
    for l = 1:ny
        fana(j,l) = cos(x(j)*kx - y(l)*ky) - cos(x(j)*kx + y(l)*ky);
    end
end

error = zeros(1,nsimul);
for i = 1:nsimul % Parcours des resultats de toutes les simulations
    data = load(output{i});
    fnum = data(end-nx+1:end,2:end);

    result = 0;
    for j = 1:nx-1
        for l = 1:ny-1
            result = result + (fana(j,l) - fnum(j,l))^2 + (fana(j+1,l) - fnum(j+1,l))^2 + (fana(j,l+1) - fnum(j,l+1))^2 + (fana(j+1,l+1) - fnum(j+1,l+1))^2;
        end
    end
    
    error(i) = 0.5*sqrt(hx*hy*result);
end

%%Polyfit : approximation linéaire
%% Figures
p = polyfit(CFL,error, 1)
fit_f = polyfit(log(CFL), log(error), 1);

figure
loglog(CFL,error,'k+','MarkerSize',12)
hold on
loglog(CFL, exp(polyval(fit_f, log(CFL))), 'k--','linewidth',lw);
set(gca,'fontsize',fs)
xlabel('$\beta_{CFL}$', 'interpreter', 'latex','fontsize', fs)
ylabel('Error')
legend('Data', strcat("fit : slope = ", string(fit_f(1))),'location', 'nw','interpreter', 'latex');
grid on