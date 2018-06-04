%The aim of this code is to save binary tiff files based on a binary
%image and a filename

%There is no output to this function.

function saveBinaryTiff(tifName,BW)

    t = Tiff(tifName, 'w');
    setTag(t,'ImageLength',size(BW,1))
    setTag(t,'ImageWidth',size(BW,2))
    setTag(t,'Photometric',Tiff.Photometric.MinIsBlack)
    setTag(t,'BitsPerSample',1)
    setTag(t,'SamplesPerPixel',1)
    setTag(t,'PlanarConfiguration',Tiff.PlanarConfiguration.Chunky)
    t.write(BW(:,:,1))

    for i = 2:size(BW,3)
        t.writeDirectory
        setTag(t,'ImageLength',size(BW,1))
        setTag(t,'ImageWidth',size(BW,2))
        setTag(t,'Photometric',Tiff.Photometric.MinIsBlack)
        setTag(t,'BitsPerSample',1)
        setTag(t,'SamplesPerPixel',1)
        setTag(t,'PlanarConfiguration',Tiff.PlanarConfiguration.Chunky)
        t.write(BW(:,:,i))
    end
    t.close 
end