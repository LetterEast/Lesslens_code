function outputImg = subpixelshift(inputImg, shift_X_in_Pixel, shift_Y_in_Pixel)
% Applying lateral shifting to the input image by muliplying the FT of the
% input image with an additional phase term
[mRow, nCol, nLayers]=size(inputImg);
outputImg = zeros(mRow, nCol, nLayers);

%[xGrid, yGrid] = meshgrid(1:nCol, 1:mRow);
[xGrid, yGrid] = meshgrid([0:1:nCol-1]/nCol, [0:1:mRow-1]/mRow);
%      xGrid  = ifftshift(xGrid);
%      yGrid  = ifftshift(yGrid);
for iLayer = 1 : nLayers
    AdditionalPhase = xGrid.*(shift_X_in_Pixel(1,iLayer)) + ...
        yGrid.*(shift_Y_in_Pixel(1,iLayer));

    H = exp(-1*1i*2*pi*AdditionalPhase);

    thisInputImgFT = fftshift(fft2(inputImg));
    thisOutputImgFT = thisInputImgFT .* H;
    thisOutputImg = abs(ifft2(ifftshift(thisOutputImgFT)));

    outputImg(:,:,iLayer) = thisOutputImg;
end
end