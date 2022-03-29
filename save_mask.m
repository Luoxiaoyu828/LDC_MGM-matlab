function save_mask(mask,mask_name)
% % ±£´æÑÚÄ¤
if exist(mask_name,'file')
    delete(mask_name);
end


import matlab.io.*
fptr = fits.createFile(mask_name);
[size_x,size_y,size_z] = size(mask);
fits.createImg(fptr,'int32',[size_x,size_y,size_z]);
fits.writeImg(fptr,mask);
fits.closeFile(fptr);
end