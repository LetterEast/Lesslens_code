function [IllumPaSet] = getMultiAngleIllum(InputImg, MNZ_result, PixelSize, WaveLength)
% this function calculates the illumination patterns of angle lights.
k0           = 2*pi/WaveLength;              % wave number
[mRow, nCol] = size(InputImg);
x            = -nCol/2*(PixelSize):PixelSize:(nCol/2-1)*PixelSize;
y            = -mRow/2*PixelSize:PixelSize:(mRow/2-1)*PixelSize;


%% 
M = MNZ_result.M;
N = MNZ_result.N;

LEDX = M * PixelSize;
LEDY = N * PixelSize;
D_LED2Sample = MNZ_result.Z;
[Object_X, Object_Y] = meshgrid(x,y);


%% Sphere Wave

r = sqrt((Object_X - LEDX).^2 + (Object_Y - LEDY).^2 + (D_LED2Sample)^2);
IllumPaSet = exp(1i*k0.*r) ./ r;

%% imshow
% phi = angle(IllumPaSet);
% amp = abs(IllumPaSet);
% figure;
% subplot(1,2,1);imshow(mat2gray(amp));xlabel("Illum amp")
% subplot(1,2,2);imshow(mat2gray(phi));xlabel("Illum phs");
% drawnow;

end