%------------------------------------------------------------------
% Introduce path of the images:
path_dir_ima = uigetdir;
ndir = dir([path_dir_ima,'/*.tif']);
ndir = ndir(arrayfun(@(x) ~strcmp(x.name(1),'.'),ndir));
num_images = numel(ndir)

sep = strfind(path_dir_ima,'/');
path_dir0 = path_dir_ima(1:sep(end-1));
sampleID = path_dir_ima(sep(end-2)+1:sep(end-1)-1);

%------------------------------------------------------------------
% Introduce the x,y,z coordinates of apical and mitral valve point:
prompt = {'Enter x-coordinate:','Enter y-coordinate:','Enter Z-coordinate:'};
title = 'Apical point';
answer = inputdlg(prompt,title);
apex_point(1,1) = str2num(answer{1});
apex_point(1,2) = str2num(answer{2});
apex_point(1,3) = str2num(answer{3});

prompt = {'Enter x-coordinate:','Enter y-coordinate:','Enter Z-coordinate:'};
title = 'Mitral valve point';
answer = inputdlg(prompt,title);
mitral_point(1,1) = str2num(answer{1});
mitral_point(1,2) = str2num(answer{2});
mitral_point(1,3) = str2num(answer{3});

center_LV(1,1:3) = apex_point;
center_LV(2,1:3) = mitral_point;

%------------------------------------------------------------------
% Interpolation of points.
cz_i = 1:1:num_images;
fx = fit(center_LV(:,3),center_LV(:,1),'poly1');
fy = fit(center_LV(:,3),center_LV(:,2),'poly1');
cx_i = feval(fx,cz_i)';
cy_i = feval(fy,cz_i)';
center_LV_i = [cx_i;cy_i;cz_i].';

p1=center_LV(1,:); %centre of the appex
p2=center_LV(2,:); %center of the mitral vave
vec1 = [p2(1)-p1(1) p2(2)-p1(2) p2(3)-p1(3)];
vec1n = vec1/norm(vec1);

%------------------------------------------------------------------
% Saving the data
path_dir = [path_dir0,'fibers/'];
if exist(path_dir, 'dir')
   warning('Fiber folder already exist!');
else
   mkdir(path_dir);
end
save([path_dir,sampleID,'_center.mat'],'center_LV_i');
save([path_dir,sampleID,'_vec1n.mat'],'vec1n');
fprintf('Centre of LV calculated ans saved!\n');