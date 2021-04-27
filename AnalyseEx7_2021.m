<pre><div class="text_to_html">%% Chargement des resultats %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ny=5; % Doit etre consistant avec l&#039;input du code C++
fichier = &#039;output&#039;;
fid=fopen([fichier,&#039;_mesh.out&#039;]); % ouverture du fichier
x=str2num(fgetl(fid)); % lit la 1e ligne du fichier, qui contient les x_i
y=str2num(fgetl(fid)); % lit la 2e ligne du fichier, qui contient les y_j
u = load([fichier,&#039;_u.out&#039;]); % lit tout le fichier des u(x_i,t_j)
data = load([fichier,&#039;_E.out&#039;]); % fichier de (t,E)
t = data(:,1);
E = data(:,2);
data = load([fichier,&#039;_f.out&#039;]);
% Le fichier contient, ligne par ligne, {t_k, {f(i,j),j=1,Nx}, i=i,Ny},
% k=1,nsteps
[ii,jj]=size(data);
Nx=jj-1;
nsteps=ii/Ny;

% Exemple 1: extraire le tableau des f(i,j) au temps final
istart=(nsteps-1)*Ny +1; % indice de la premiere ligne a lire du tableau global (data)
iend=istart+Ny-1; % indice de la derniere ligne a lire du tableau global (data)\
f=(data(istart:iend,2:jj))&#039;; % extrait le f(i,j); le &#039; transpose pour avoir f(x_i,y_j)
fs=16; % font size
% Contour plot de f(x,y)
figure
contourf(x,y,f&#039;) % lignes de niveaux avec coloriage (le &#039; transpose)
% voir aussi contour, surf, surfc, ...
set(gca,&#039;fontsize&#039;,fs) % set font size for labels, title, etc
xlabel(&#039;x [m]&#039;)
ylabel(&#039;y [m]&#039;)
colorbar

% Exemple 2: une coupe a y=y0=const de f(x,y,t) --&gt; fcut(x_i,t_k)
% extraire le tableau des f(i,j,k) : f(x,y=y0,t) a y0 fixe.
i0=2; % indice pour y_0 fixe (ici on prend le 2e point de maillage en y)
fcut=(data(i0:Ny:end,2:jj))&#039;; % toutes les Ny lignes, on est au meme y
t=data(i0:Ny:end,1); % la 1e colonne contient le temps

% contour plot de f(x,y0,t)

figure
contourf(x,t,fcut&#039;)
set(gca,&#039;fontsize&#039;,fs) % set font size for labels, title, etc
xlabel(&#039;x [m]&#039;)
ylabel(&#039;t [s]&#039;)
colorbar


</div></pre>