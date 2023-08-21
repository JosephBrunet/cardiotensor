% MAIN SCRIPT FOR CALCULATING THE FIBRE ORIENTATION IN THE WHOLE HEART

close all
clear, clc
% ===============================================================
% Images path
% ===============================================================
% path_ima = '/Volumes/PGC2/FETAL HEARTS/New samples - 2017/2537/images/2537_scaled_crop_flip/';
path_ima = uigetdir;
path_ima = [path_ima,'/'];
ndir_ima = dir([path_ima,'*.tif']);
ndir_ima = ndir_ima(arrayfun(@(x) ~strcmp(x.name(1),'.'),ndir_ima));

% ===============================================================
% Fibers path
% ===============================================================
% path_fibers = '/Volumes/PGC2/FETAL HEARTS/New samples - 2017/2537/fibers/';
path_fibers = uigetdir;
path_fibers = [path_fibers,'/'];

if exist([path_fibers,'FI/'], 'dir')
   warning('Fiber FI folder already exist!');
else
   mkdir([path_fibers,'FI/']);
end

if exist([path_fibers,'FA/'], 'dir')
   warning('Fiber VF folder already exist!');
else
   mkdir([path_fibers,'FA/']);
end

if exist([path_fibers,'VF/'], 'dir')
   warning('Fiber VF folder already exist!');
else
   mkdir([path_fibers,'VF/']);
end

% ===============================================================
% Load fibers files: LV center and LV vector (vec1n)
% ===============================================================
files_fibers = dir([path_fibers,'*.mat']) ;
files_fibers = files_fibers(arrayfun(@(x) ~strcmp(x.name(1),'.'),files_fibers));
if (numel(files_fibers) < 2)
    error('One fiber file missing!');
else
    load([path_fibers,files_fibers(1).name]);
    load([path_fibers,files_fibers(2).name]);
end

% ===============================================================
% Number of consecutive images to calculate structure tensor. 
% Nz mut be an odd number (9, 11, 13 etc.); Minimum Nz = 9;
% Definition of first and last image to compute fibre orientation
% ===============================================================
Nz = 3;
ini_ima = center_LV_i(1,3);
fin_ima = center_LV_i(end,3);
N = fin_ima - ini_ima - (Nz-1);

% ===============================================================
% Definition of input variables
% ===============================================================
% Usually, if first image is the apex and last is the base, 
% flip = 0, but check it first in 1 image before runngin the whole volume.
bool_flip = 0;
bool_correct = 1;
FilterType = 1;

% ===============================================================
% START COMPUTATION
% ===============================================================
fprintf ('Starting fibre analysis...\n');
for j = 1 : 1 : N
   if j ~= 155
      continue
   end
    cima = j - 1 + ini_ima + (Nz-1)*0.5;
    [pathstr,name,ext] = fileparts([path_ima,ndir_ima(cima).name]);
    k=1;
    for i = cima-(Nz-1)*0.5 : cima+(Nz-1)*0.5
        I3d(:,:,k)=imread([path_ima,ndir_ima(i).name]);
        k = k+1;
    end
    center = center_LV_i(j+(Nz-1)*0.5,1:2);
    [VF,FI,EV,~] = Volume2ST2VF(I3d, center, vec1n, bool_correct, bool_flip, FilterType);
    FI = single(FI);  
    EV = single(EV);
    VF = single(VF);
    l1 = EV(:,:,1);
    l2 = EV(:,:,2);
    l3 = EV(:,:,3);
    lm = (l1 + l2 + l3)./3;
    FA = sqrt(3/2)*(sqrt((l1-lm).^2 +(l2-lm).^2 + (l3-lm).^2)./sqrt(l1.^2 + l2.^2 +l3.^2));
    
% ===============================================================
%  Saving the data    
% ===============================================================
    save([path_fibers,'FI/',name,'_FI.mat'],'FI');
    save([path_fibers,'FA/',name,'_FA.mat'],'FA');
    save([path_fibers,'VF/',name,'_VF.mat'],'VF');
    
    fprintf ('Image %d of %d done!\n', j, N);
    
% ===============================================================
%  Plot FI   
% ===============================================================
   if j == 155
      close(figure(1)), figure(1), imagesc(rad2deg(FI)); hold on, 
      cmap = [0 0 0; hsv(180)]; colormap(cmap);
      plot(center(1,1),center(1,2),'xr','MarkerSize',20);
   end
   clear I3d FI EV FA VF center


end
fprintf ('Fibre analysis done!\n');  


