%The aim of this code is to save [2 4 8 16] bits tiff files based on a
%filename, a dataset(movie) and a number of bit

%There is no output to this function.

function nBTiff(tifName,mov,bit)
    assert(ischar(tifName),'Filename needs to be a char');
    assert(ismember(length(size(mov)),[2 3]), 'The data you are trying to save has an unexpected dimension');
    assert(ismember(bit,[2 4 8 16 32]), 'Unexpected number of bit for saving tiff');

    t = Tiff(tifName, 'w');
    setTag(t,'ImageLength',size(mov,1))
    setTag(t,'ImageWidth',size(mov,2))
    setTag(t,'Photometric',Tiff.Photometric.MinIsBlack)
    setTag(t,'BitsPerSample',bit)
    setTag(t,'SamplesPerPixel',1)
    setTag(t,'PlanarConfiguration',Tiff.PlanarConfiguration.Chunky)
    t.write(mov(:,:,1))

    for i = 2:size(mov,3)
        t.writeDirectory
        setTag(t,'ImageLength',size(mov,1))
        setTag(t,'ImageWidth',size(mov,2))
        setTag(t,'Photometric',Tiff.Photometric.MinIsBlack)
        setTag(t,'BitsPerSample',bit)
        setTag(t,'SamplesPerPixel',1)
        setTag(t,'PlanarConfiguration',Tiff.PlanarConfiguration.Chunky)
        t.write(mov(:,:,i))
    end    
    end