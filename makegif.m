function makegif(frames,numframes,str)

delay=0.03;
for i = 1:numframes
im = frame2im(frames(i)); 
[imind,cm] = rgb2ind(im,256);
      
      % Write to the GIF File 
      if i == 1 
          imwrite(imind,cm,str,'gif','DelayTime', delay, 'Loopcount',inf); 
      else 
          imwrite(imind,cm,str,'gif','WriteMode','append','DelayTime', delay); 
      end 
      
end