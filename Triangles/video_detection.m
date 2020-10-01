
clc
clear 


Z = peaks;
surf(Z); 
axis tight manual 
set(gca,'nextplot','replacechildren'); 

% Create a video writer object for the output video file and open the object for writing.

v = VideoWriter('peaks.avi', 'Uncompressed AVI');
open(v);

% Generate a set of frames, get the frame from the figure, and then write each frame to the file.

for k = 1:20 
    
   surf(sin(2*pi*k/20)*Z,Z)
   frame = getframe(gcf);
   writeVideo(v,frame);
   
%    ss(:,:,k) = sin(2*pi*k/20)*Z;
%    Ss(:,:,k) = uint8(ss(:,:,k));
end

%close(v);

%% EXAMPLE:

% load mristack
% 
% implay(mristack)
% 
% implay(ss)
% implay(frame)
implay('peaks.avi');



