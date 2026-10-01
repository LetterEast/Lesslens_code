function writeMaskedPNG(pixels, filename, mask)
%WRITEMASKEDPNG Save unsupported pixels as black, without transparency.
assert(isequal(size(mask),[size(pixels,1),size(pixels,2)]), ...
    'Output mask must match the image dimensions.');
invalid = repmat(~logical(mask),1,1,size(pixels,3));
pixels(invalid) = 0;
imwrite(pixels,filename);
end
