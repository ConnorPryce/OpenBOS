function image_ndgrid = image_to_ndgrid(image)
    image_ndgrid = pagetranspose(flip(image,1));
end

