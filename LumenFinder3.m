function [ GraftVesselDensity, MyocardialVesselDensity ] = LumenFinder3( FileName )
% LumenFinder funct>>ion detects luminal structures and counts vessel density on
% fluorescently stained tissue sections (image scale 

n = 2; %this is the number of images in the directory

for i = 1:n
    FileName = sprintf('.jpg',i); % Previous version: FileName = input('Open file: ', 's');
    % Load file and define graft and myocardial regions
    fileIn = imread (FileName); %read file
    scale = 0.628; %For 1024x1024 20X images: 0.628 microns/pixel

    [graftMask,GX,GY] = (roipoly(fileIn));%draw region of interest
    
    Graft_area = polyarea(GX, GY); %calculates area of ROI
    graftMask = uint8(cat(3,graftMask,graftMask,graftMask));%apply mask to all levels of colored image
    Graft = rgb2gray(fileIn.*graftMask); %convert to grayscale

    [myocardiumMask,MX,MY] = (roipoly(fileIn)); %repeat for myocardial region
    Myocardium_area = polyarea(MX, MY);
    myocardiumMask = uint8(cat(3,myocardiumMask,myocardiumMask,myocardiumMask));
    Myocardium = rgb2gray(fileIn.*myocardiumMask);

    Mask = Graft + Myocardium; %create display mask with both ROIs
    a = figure; hold on 
    imshow(Mask); title('Graft and Myocardium ROI'); saveas(a, ['ROI_' FileName], 'tif')

    %Graft thresholding and particle count
    Graft_threshold = im2bw(Graft,0.08);
    b = figure; hold on 
    imshow(Graft_threshold); title('Graft threshold'); saveas(b, ['Graft_thresh_' FileName], 'tif')

    Graft_threshold_filter = bwareaopen(Graft_threshold, 80);
    c = figure; hold on
    imshow(Graft_threshold_filter); title('Filtered Graft threshold'); saveas(c, ['Graft_thresh_filter_' FileName], 'tif')

    [B,L,N,A] = bwboundaries(Graft_threshold_filter,8,'noholes');

    d = figure; imshow(label2rgb(L, @jet, [.5 .5 .5]))
    hold on
    title('Graft Lumens')
        for k = 1:length(B)
            boundary = B{k};
            plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2)
            graft_LA (k) = polyarea(boundary(:,1),boundary(:,2));
        end
    saveas(d, ['Graft_lumens_' FileName], 'tif')

    %Myocardium thresholding and particle count
    Myocardium_threshold = im2bw(Myocardium,0.08);
    e = figure; hold on 
    imshow(Myocardium_threshold); title('Myocardium threshold'); saveas(e, ['Myocard_thresh_' FileName], 'tif')

    Myocardium_threshold_filter = bwareaopen(Myocardium_threshold, 80);
    f = figure; hold on
    imshow(Myocardium_threshold_filter); title('Filtered Myocardium threshold'); saveas(f, ['Myocard_thresh_filter_' FileName], 'tif')

    [B2,L2,N2,A2] = bwboundaries(Myocardium_threshold_filter,8,'noholes');

    g = figure; imshow(label2rgb(L2, @jet, [.5 .5 .5]))
    hold on
    title('Myocardium Lumens')
        for j = 1:length(B2)
            boundary2 = B2{j};
            plot(boundary2(:,2), boundary2(:,1), 'w', 'LineWidth', 2)
            myocardium_LA (j) = polyarea(boundary2(:,1),boundary2(:,2));
        end
    saveas(g, ['Myocardium_lumens_' FileName], 'tif')

    %converst graft area from pixels^2 to microns^2 to mm^2
    GraftArea = Graft_area*(scale^2)/(1000^2); 
    MyocardiumArea = Myocardium_area*(scale^2)/(1000^2);

    %find average lumen size and convert from pixels^2 to microns^2 
    %vessel area (area of graft vessels as percent of total graft area
    if isempty(B2)== 1
        myocardium_lumen_size  = 0;
        Myocardium_percent_area = 0;
    else
        myocardium_lumen_size = mean(myocardium_LA)*(scale^2);
        Myocardium_percent_area = 100*sum(myocardium_LA)/Myocardium_area; 
    end

    if isempty(B)== 1
        graft_lumen_size = 0;
        Graft_percent_area = 0;
    else
        graft_lumen_size = mean(graft_LA)*(scale^2);
        Graft_percent_area = 100*sum(graft_LA)/Graft_area;
    end

    Graft_vessels = N;
    GraftVesselDensity = Graft_vessels/GraftArea; %Number of vessels/mm^2

    VascularDensity = [Graft_vessels, graft_lumen_size, Graft_percent_area, GraftVesselDensity];
   
    disp(VascularDensity)

end








