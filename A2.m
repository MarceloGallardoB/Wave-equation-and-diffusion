lw = 2; fs = 16; Nx = 64;
nom = 'try.gif';

x = [0,10];
y = [0,6];
[X,Y] = meshgrid(0:0.1:10,0:0.1:6);
t = [];
for i=1:129
    t = [t; i/100];
end

size(t)
h = figure;
for i = 1:size(t)
    %Dessin de la fonction
    surf(X,Y,f(X,Y,t(i)));
    t(i)
    xlim([0 10])
    ylim([0 6])
    zlim([-3 3])
    xlabel('X', 'fontsize', fs)
    ylabel('Y', 'fontsize', fs)
    zlabel('Z', 'fontsize', fs)

     drawnow
     pause(0.1);

  frame = getframe(h); 
  im = frame2im(frame); 
  [imind,cm] = rgb2ind(im,256); 

  if i == 1 
      imwrite(imind,cm,nom,'gif', 'Loopcount',inf); 
  else 
      imwrite(imind,cm,nom,'gif','WriteMode','append'); 
  end 

end


function ave = f(x,y,t)
    ave = real(exp(i*pi*(x/5 + y/3 + (sqrt(34)/5)*4*t)));
end