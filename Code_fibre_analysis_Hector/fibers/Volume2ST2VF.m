function [VF, FI, EV, FT] = Volume2ST2VF(I3d, center, vec1n, correct, flip, FilterType)
% AUTHOR: Vedrana B & Patricia G.
%
% FUNCTION DESCRIPTION:
%   Function Volume2ST2VF calculates a vector field for the central slice
%   in the volumetric image I3d, using a structure tensor method and PCA.
%   It also returns FI as an angle between vector components Z and X, with 
%   respect to the given center of cylindrical coord. system of the image,
%   also known as helical angle, EV as the eigenvalues for each pixel, and 
%   FT as the transverse angle. 
%
% INPUTS:
%   - I3d = volumetric image, size (MxNxL), which L should be bigger than
%   filter size. Minimum L = 9;
%   - center = coordinate of the center of the LV.
%   - vec1n = coordinats of the normalized LV axis vector.
%   - correct = 1 if you want to correct the missaligment between
%     LV long axis (vec1n) and the image system vertical axis [0 0 1]. 
%     correct = 0 if you don't.
%   - flip = 1: flip the Z coordinates when calculate the strcuture tensor. 
%     FI in endocardium shoul be positive, and negative in epicardiu, so if
%     is in the other way just select flip = 1; Otherwise keep flip = 0;
%   - FilterType = 1 or 2. FilterType = 1 uses Finite difference filter 
%     to calculate the granit FlterType = 2 used Solber filter. Defasult is 
%      FilterType = 1;
%
% EXAMPLE CALL:
%   [VF, FI, EV, FT] = Volume2ST2VF(I3d, center, vec1n, 1, 0, 1)
%   [VF, FI, EV, FT] = Volume2ST2VF(I3d, [700 700 5], [0.0 0.0 1.0], 0, 0, 1)

%%

[h, w, d] = size(I3d);
wd = w*d;
hw = h*w;
k  = (d+1)/2;

% ===================================================
% CALUCLATE GRADIENTS
% ===================================================

nT = d - 2;

if(FilterType==1) % Finite difference
    [gx,gy,gz] = imgradientxyz(I3d,'central');
    if(flip==0)
        gz = -gz;
    end
else % Sober filter
    [gx,gy,gz] = imgradientxyz(I3d,'sobel');
    if(flip==0)
        gz = -gz;
    end
end

% ===================================================
% CALUCLATE SUMS OF GRADIENT PRODUCTS
% ===================================================

sumMask3  = ones(1,1,nT);
sumMask22 = ones(nT,nT);

tmp2d = reshape(gx.*gx,[h wd]);
tmp2d = conv2(tmp2d,sumMask22,'same');
tmp3d = reshape(tmp2d,[h w d]);
tmp3d = convn(tmp3d,sumMask3,'same');
Sgxgx = reshape(tmp3d(:,:,k),[hw 1]);

tmp2d = reshape(gy.*gy,[h wd]);
tmp2d = conv2(tmp2d,sumMask22,'same');
tmp3d = reshape(tmp2d,[h w d]);
tmp3d = convn(tmp3d,sumMask3,'same');
Sgygy = reshape(tmp3d(:,:,k),[hw 1]);

tmp2d = reshape(gz.*gz,[h wd]);
tmp2d = conv2(tmp2d,sumMask22,'same');
tmp3d = reshape(tmp2d,[h w d]);
tmp3d = convn(tmp3d,sumMask3,'same');
Sgzgz = reshape(tmp3d(:,:,k),[hw 1]);

tmp2d = reshape(gx.*gy,[h wd]);
tmp2d = conv2(tmp2d,sumMask22,'same');
tmp3d = reshape(tmp2d,[h w d]);
tmp3d = convn(tmp3d,sumMask3,'same');
Sgxgy = reshape(tmp3d(:,:,k),[hw 1]);

tmp2d = reshape(gx.*gz,[h wd]);
tmp2d = conv2(tmp2d,sumMask22,'same');
tmp3d = reshape(tmp2d,[h w d]);
tmp3d = convn(tmp3d,sumMask3,'same');
Sgxgz = reshape(tmp3d(:,:,k),[hw 1]);

tmp2d = reshape(gy.*gz,[h wd]);
tmp2d = conv2(tmp2d,sumMask22,'same');
tmp3d = reshape(tmp2d,[h w d]);
tmp3d = convn(tmp3d,sumMask3,'same');
Sgygz = reshape(tmp3d(:,:,k),[hw 1]);

clear gx gy gz tmp2d tmp3d

% ===================================================
% FIND SMALLEST EIGENVALUE - SOLVE 3rd ORDER EQUATION
% ===================================================

d = Sgxgx.*(- Sgygz.^2 + Sgygy.*Sgzgz) + Sgxgz.*(Sgxgy.*Sgygz - Sgxgz.*Sgygy) + Sgxgy.*(Sgxgz.*Sgygz - Sgxgy.*Sgzgz);
c = Sgxgy.^2 + Sgxgz.^2 + Sgygz.^2 - Sgxgx.*(Sgygy + Sgzgz) - Sgygy.*Sgzgz;
b = Sgxgx + Sgygy + Sgzgz;
a = -ones(size(d));

ind_a0 = a == 0;
ind_b0 = b == 0;
ind_c0 = c == 0;

Eigens = zeros(hw,3);
nn = zeros(hw,1) * NaN;

ind_ord0 = ind_a0 & ind_b0 &  ind_c0; 
ind_ord1 = ind_a0 & ind_b0 & ~ind_c0;
ind_ord2 = ind_a0 - ind_ord0 - ind_ord1;

delta = b.*b - 4*a.*c;
sqdelta = sqrt(delta);
ind_ord2_1 = ind_ord2 & (delta==0);
ind_ord2_2 = ind_ord2 & (delta~=0);
inva = 0.5./a(ind_ord2_2);

ind_ord3 = ~ind_a0;

Eigens(ind_ord0,:) = [nn(ind_ord0) nn(ind_ord0) nn(ind_ord0)];
Eigens(ind_ord1,:) = [-d(ind_ord1)./c(ind_ord1) nn(ind_ord1) nn(ind_ord1)];
Eigens(ind_ord2_1,:) = -0.5* [b(ind_ord2_1)./a(ind_ord2_1) b(ind_ord2_1)./a(ind_ord2_1) nn(ind_ord2_1)];
Eigens(ind_ord2_2,:) = [(-b(ind_ord2_2)+sqdelta(ind_ord2_2)).*inva (-b(ind_ord2_2)-sqdelta(ind_ord2_2)).*inva nn(ind_ord2_2)];
Eigens(ind_ord3,:) = solveOrd3Eigens( a(ind_ord3), b(ind_ord3), c(ind_ord3), d(ind_ord3));

Em = min(Eigens,[],2); %(min(Eigens'))';

% ===================================================
% RESHAPE EIGENVALUES FROM VECTORIZED TO 2D FORMAT
% ===================================================

EV(:,:,3) = reshape(Eigens(:,3), [h,w]);
EV(:,:,2) = reshape(Eigens(:,2), [h,w]);
EV(:,:,1) = reshape(Eigens(:,1), [h,w]);

% ===================================================
% FIND CORRESPONDING EIGENVECTOR
% ===================================================

num = (Em.*Sgygz + Sgxgy.*Sgxgz - Sgxgx.*Sgygz);
den =  Em.^2 - Em.*Sgxgx - Em.*Sgygy -Sgxgy.^2 + Sgxgx.*Sgygy;
v3 = ones(size(Em));
v2 = num./den;
v1 = (Sgxgy.*v2 + Sgxgz)./(Em-Sgxgx);
Ampl = sqrt(v1.^2 + v2.^2 + v3.^2);
posdef = sum(Eigens<0,2) == 0;

FV(:,3) = v3./Ampl.*posdef;
FV(:,2) = v2./Ampl.*posdef;
FV(:,1) = v1./Ampl.*posdef;

% =====================================================
% ROTATE THE LONG LV AXIS TO ALINEATE WITH ACQ. AXIS
% =====================================================
if(correct)
    ssc = @(v) [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
    RU = @(A,B) eye(3) + ssc(cross(A,B)) + ssc(cross(A,B))^2*(1-dot(A,B))/(norm(cross(A,B))^2);

    RUU = RU([0 0 1]',vec1n');
    FV=RUU*FV'; %el nuevo eje que coincide con el que tenemos
    FV = FV';
end

% ===================================================
% RESHAPE FROM VECTORIZED TO 2D FORMAT
% ===================================================

VF(:,:,3) = reshape(FV(:,3), [h,w]);
VF(:,:,2) = reshape(FV(:,2), [h,w]);
VF(:,:,1) = reshape(FV(:,1), [h,w]);

% =====================================================
% CHANGE THE COORDINATE SYSTEM FROM CARTESIAN TO RADIAL
% =====================================================

[Cx, Cy] = deal(center(1),center(2)); 
[X, Y]   = meshgrid((1:w)-Cx,(1:h)-Cy);
R = sqrt((X.^2 + Y.^2));
gama = -atan2(Y(:),X(:));
cg   = cos(gama);
sg   = sin(gama);

FVrot(:,3) = FV(:,3);
FVrot(:,1) = cg.*FV(:,1) - sg.*FV(:,2);
FVrot(:,2) = sg.*FV(:,1) + cg.*FV(:,2);

VFrot(:,:,3) = reshape(FVrot(:,3), [h,w]);
VFrot(:,:,1) = reshape(FVrot(:,1), [h,w]);
VFrot(:,:,2) = reshape(FVrot(:,2), [h,w]);

% =====================================================
% CALCULATE VECTOR PROJECTIONS IN RADIAL COOR. SYSTEM
% =====================================================

FT = atan(VFrot(:,:,1)./VFrot(:,:,2));
FI = atan(VFrot(:,:,3) ./ (VFrot(:,:,2)));
FI_unwrapped = atan2(VFrot(:,:,3), (VFrot(:,:,2)));
ind_wrapped = (FI_unwrapped-FI)~=0;
VF_ind = repmat(ind_wrapped,[1 1 3]);
VF(VF_ind) = VF(VF_ind)*(-1);

end

% =====================================================
% FUNCTION FOR DEFINING GRADIENT MASK
% =====================================================
function maskDiff = centralDiff(nG)
    if nG==1 
        maskDiff = [-1 1];
    elseif nG==2 
        maskDiff = [ -1/2 0 1/2];
    elseif nG==4 
        maskDiff = [1/12 -2/3 0 2/3 -1/12];
    elseif nG==6 
        maskDiff = [-1/60  3/20 -3/4 0 3/4 -3/20 1/60];
    elseif nG==8 
        maskDiff = [1/280 -4/105 1/5  -4/5 0 4/5 -1/5  4/105 1/280]; 
    end
end
% =====================================================
% FUNCTION FOR DEFINING CODE PARAMETERS
% =====================================================

function codePar = aboutCode()

codePar.nT = 7; %7
codePar.nG = 2; %2
codePar.d = codePar.nT + codePar.nG - 1;

end 

% =====================================================
% FUNCTION FOR SOLVING 3rd ORDER POLYNOMIAL
% =====================================================

function x = solveOrd3Eigens(a,b,c,d)

numElem = length(a);
%  /* calculate the normalized form x^3 + a2 * x^2 + a1 * x + a0 = 0 */
inv_a = 1 ./ a;
b_a = inv_a .* b; b_a2 = b_a .* b_a;
c_a = inv_a .* c;
d_a = inv_a .* d;

%  /* solve the cubic equation */
Q = (3 * c_a - b_a2) / 9;
R = (9 * b_a .* c_a - 27 * d_a - 2 * b_a .* b_a2) / 54;
Q3 = Q .* Q .* Q;
D = Q3 + R .* R;
b_a_3 = (1 ./ 3) * b_a;

th  = acos(R ./ sqrt(-Q3));
sqQ = sqrt(-Q);
x   = zeros(numElem,3);
nn  = zeros(numElem,1) * NaN;

ind_Q0R0  = (Q==0) & (R==0);
ind_Q0RX  = (Q==0) & ~(R==0);
x(ind_Q0R0,:)  = [  -b_a_3(ind_Q0R0) ...
                    -b_a_3(ind_Q0R0) ...
                    -b_a_3(ind_Q0R0)];
x(ind_Q0RX,:)  = [  (2*R(ind_Q0RX)).^(1./3) - b_a_3(ind_Q0RX) ... 
                    nn(ind_Q0RX) ...
                    nn(ind_Q0RX)];
                
                
ind_Dleq0 = (D<=0) & ~(ind_Q0R0 | ind_Q0RX) ;
x(ind_Dleq0,:) = [  2 * sqQ(ind_Dleq0) .* cos(th(ind_Dleq0)/ 3.0) - b_a_3(ind_Dleq0), ...
                    2 * sqQ(ind_Dleq0) .* cos( (th(ind_Dleq0)+2*pi)/ 3.0 ) - b_a_3(ind_Dleq0),...
                    2 * sqQ(ind_Dleq0) .* cos( (th(ind_Dleq0)+4*pi)/ 3.0 ) - b_a_3(ind_Dleq0)];

indAllPrev = ~(ind_Dleq0 | ind_Q0R0 | ind_Q0RX);
ind_Rg0   =  R>0  & indAllPrev;
ind_Rl0   =  R<0  & indAllPrev;
ind_R0    =  R==0 & indAllPrev;
ind_R = ind_Rg0 | ind_Rl0  | ind_R0;

multi = ind_Rg0*1 + ind_Rl0*(-1) + ind_R0*0; 
AD = ((abs(R)+sqrt(D)) .^ (1/3)) .* multi;
BD = AD==0 .* 0 + AD~=0 .* (-Q ./ AD);

x(ind_R,:) = [AD(ind_R) + BD(ind_R) - b_a_3(ind_R) ...
     nn(ind_R),...
     nn(ind_R)];

end
