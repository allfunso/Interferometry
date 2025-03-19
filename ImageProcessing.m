clear
clc

% Esfera
% x1 = 1200;
% x2 = 3100;
% y1 = 900;
% y2 = 2700;
% imagen = "Esfera";
% h = 1.6; % ALtura real
% x = 4.47;  %1700p = 4 cm 
% y = 4.23;

% Engrane
x1 = 100;
x2 = 1000;
y1 = 100;
y2 = 1000;
imagen = "Engrane";
h = 1;
x = 6.46; %2135p = 6 cmo
y = 6.46;

% Audífono
% x1 = 1200;
% x2 = 3500;
% y1 = 400;
% y2 = 2400;
% imagen = "Emblema";
% h = 0.1;  %2680 = 8 cm
% x = 6.86;
% y = 5.97;

xs = linspace(0,x,(x2-x1)+1);
ys = linspace(0,y,(y2-y1)+1);

im11 = im2double(rgb2gray(imread(imagen+"1.JPG")));
im12 = im2double(rgb2gray(imread(imagen+"2.JPG")));
im13 = im2double(rgb2gray(imread(imagen+"3.JPG")));
im14 = im2double(rgb2gray(imread(imagen+"4.JPG")));

im11 = im11(y1:y2, x1:x2); %y, z
im12 = im12(y1:y2, x1:x2);
im13 = im13(y1:y2, x1:x2);
im14 = im14(y1:y2, x1:x2);

im21 = im2double(rgb2gray(imread(imagen+"5.JPG")));
im22 = im2double(rgb2gray(imread(imagen+"6.JPG")));
im23 = im2double(rgb2gray(imread(imagen+"7.JPG")));
im24 = im2double(rgb2gray(imread(imagen+"8.JPG")));

im21 = im21(y1:y2, x1:x2); %y, z
im22 = im22(y1:y2, x1:x2);
im23 = im23(y1:y2, x1:x2);
im24 = im24(y1:y2, x1:x2);


im15 = im2double(rgb2gray(imread("Plano"+imagen+"1.JPG")));
im16 = im2double(rgb2gray(imread("Plano"+imagen+"2.JPG")));
im17 = im2double(rgb2gray(imread("Plano"+imagen+"3.JPG")));
im18 = im2double(rgb2gray(imread("Plano"+imagen+"4.JPG")));

im15 = im15(y1:y2, x1:x2); %y, z
im16 = im16(y1:y2, x1:x2);
im17 = im17(y1:y2, x1:x2);
im18 = im18(y1:y2, x1:x2);

im25 = im2double(rgb2gray(imread("Plano"+imagen+"5.JPG")));
im26 = im2double(rgb2gray(imread("Plano"+imagen+"6.JPG")));
im27 = im2double(rgb2gray(imread("Plano"+imagen+"7.JPG")));
im28 = im2double(rgb2gray(imread("Plano"+imagen+"8.JPG")));

im25 = im25(y1:y2, x1:x2); %y, z
im26 = im26(y1:y2, x1:x2);
im27 = im27(y1:y2, x1:x2);
im28 = im28(y1:y2, x1:x2);

w1 = 0.5;
w2 = 1-w1;

im1 = w1*im11 + w2*im21;
im2 = w1*im12 + w2*im22;
im3 = w1*im13 + w2*im23;
im4 = w1*im14 + w2*im24;
im5 = w1*im15 + w2*im25;
im6 = w1*im16 + w2*im26;
im7 = w1*im17 + w2*im27;
im8 = w1*im18 + w2*im28;

fase_objeto = (atan2(-(im2-im4),(im1-im3)))+pi;
fase_plano = (atan2(-(im6-im8),(im5-im7)))+pi;


figure(1)
imagesc(xs,ys,fase_objeto); colormap gray, axis equal;
figure(2)
imagesc(xs,ys,fase_plano); colormap gray, axis equal;

filtro = 4;

[unwri1, aunwri1, k1] = unwrapDet(fase_objeto, filtro);
[unwri2, aunwri2, k2] = unwrapDet(fase_plano, filtro);

aUnwri0 = aunwri2 - aunwri1;
Unwri0 = unwri2 - unwri1;
aUnwri0 = 4/3*(aUnwri0 + max(aUnwri0(:))/2);
aUnwri0 (aUnwri0<0) =0;
% aUnwri0 = aUnwri0*h/max(aUnwri0(:));


% Correción de perspectiva

a = 0.14;
alfa = -24*pi/180;
beta = 0;
aUnwri0 = aUnwri0./(2*pi)*(a/(tan(alfa)+tan(beta)));

figure(3)
imagesc(xs,ys,aUnwri0), colormap gray; axis equal;
%imagesc(-1*Unwri0), colormap gray;

figure(4)
surf(xs,ys,aUnwri0), colormap gray; 
%surf(-1*Unwri0), colormap gray;
shading interp
axis equal



function [psi, apsi, K] = unwrapDet(phi, sigma)

% Input
% phi :   wrapped phase map with values in the interval [-pi,pi)
% sigma : (Optional) if present, a Gaussian filter with parameter
%         'sigma' is applied to reduce the effect of noise. 
%
% Output
% psi :   "exact" solution for unwrapped phase
% apsi :  approximate solution for unwrapped phase
% K :     integer field that unwraps the phase
%
% Reference
% + Vyacheslav V. Volkov and Yimei Zhu, "Deterministic phase
%   unwrapping in the presence of noise", Opt. Lett. 28, 2156-2158,
%   2003

% This software is in the Public Domain   
% Initial version, Ulf Griesmann, NIST, April 2016

    % continuous function Z(x,y)) from wrapped phase
    Z = exp(1i*phi); % Eq. 4
    
    % apply filter
    if ~isempty(sigma)
        Z = gaussfilterZ(Z, sigma);
    end
    
    % gradient of unwrapped psi(x,y)
    [DZx, DZy] = gradient(Z);
    Dpsix = real(-1i*DZx ./ Z);
    Dpsiy = real(-1i*DZy ./ Z);
    
    % integrate with Frankot-Chellappa --> approximate solution
    apsi = intgrad(Dpsix, Dpsiy); % Eq. 3
    
    % calculate K(x,y) and "exact" solution
    K = round(0.5 * (apsi - phi) / pi);
    psi = 2*pi * K + phi; % Eq. 1
    
end

function z = intgrad(dzdx,dzdy)

   
    if ~all(size(dzdx) == size(dzdy))
       error('Gradient matrices must match');
    end

    [rows,cols] = size(dzdx);

    
    [wx, wy] = meshgrid(([1:cols]-(fix(cols/2)+1))/(cols-mod(cols,2)), ...
			([1:rows]-(fix(rows/2)+1))/(rows-mod(rows,2)));

    wx = ifftshift(wx); 
    wy = ifftshift(wy);

    DZDX = fft2(dzdx);   % Fourier transforms of derivatives
    DZDY = fft2(dzdy);

    Z = (-j*wx.*DZDX -j*wy.*DZDY)./(wx.^2 + wy.^2 + eps);  % Equation 21
    
    z = 0.5 * real(ifft2(Z)) / pi;  % Reconstruction (divided by 2pi; U.G. 16apr2016)

end

function fzmap = gaussfilterZ(zmap, sigma)


    % calculate filter size
    sze = ceil(6*sigma);  
    if ~mod(sze,2)    % Ensure filter size is odd
        sze = sze+1;
    end
    sze = max(sze,1); % and make sure it is at least 1
    
    % calculate Gaussian kernel
    h = fspecial('gaussian', [sze sze], sigma);
    
    % apply filter
    fzmap = complex(filter2(h, real(zmap)), filter2(h, imag(zmap)));
    
end


% Explicación del código

% Se le introduce la fase y el parámetro del filtro Gaussiano,
% de ahí se genera una función Z, posteriormente se le aplica el filtro
% con lla función de fspecial para generar el kernel de un filtro gaussiano
% y que se le aplica de manera separada el Kernel a la parte real y la imaginaria
% de ahí se  obtine el gradiente de Z y se integrá con el étodo de Frankot-Chellappa
% que consiste en que con el gradiente de la fase envuelta en x y el gradiente en y
% se obtenga la solución aproximada; para dicho método se calculan las
% trnaformadas de fourier de dzdx y dzdy se construye la función de
% frecuencia en base a la ecuación 21 (quien sabe cual sea, el paper solo llega a 4)
% y se saca la transformada inversa. Ya con la solución aproximada se
% calcula el campo entero que desenvuelve la fase (K) y con la ecuación 1
% se desenvuelve la fase al multiplicar 2\pi * K + phi(la fase envuelta)
% obteniendo la fase desenvuelta exacta
