%The aim of this code is to write a [2 4 8 16] bits tiff files based on a
%tiff object, a dataset(movie) and a number of bit .
%
    %this function is similar to nBTiff but allow to control the writing
    %outside of the function which is convenient to write data by small
    %chunks instead of giving the full dataset as input.

function [tifObject] = writeTiff(tifObject,mov,bit)
    assert(isa(tifObject,'Tiff'),'Input 1 should be a tiff');
    assert(ismember(length(size(mov)),[2 3]), 'The data you are trying to save has an unexpected dimension');
    assert(ismember(bit,[2 4 8 16 32]), 'Unexpected number of bit for saving tiff');

    setTag(tifObject,'ImageLength',size(mov,1))
    setTag(tifObject,'ImageWidth',size(mov,2))
    setTag(tifObject,'Photometric',Tiff.Photometric.MinIsBlack)
    setTag(tifObject,'BitsPerSample',bit)
    setTag(tifObject,'SamplesPerPixel',1)
    setTag(tifObject,'PlanarConfiguration',Tiff.PlanarConfiguration.Chunky)
    tifObject.write(mov(:,:,1))

    for i = 2:size(mov,3)
        tifObject.writeDirectory
        setTag(tifObject,'ImageLength',size(mov,1))
        setTag(tifObject,'ImageWidth',size(mov,2))
        setTag(tifObject,'Photometric',Tiff.Photometric.MinIsBlack)
        setTag(tifObject,'BitsPerSample',bit)
        setTag(tifObject,'SamplesPerPixel',1)
        setTag(tifObject,'PlanarConfiguration',Tiff.PlanarConfiguration.Chunky)
        tifObject.write(mov(:,:,i))
    end    
end