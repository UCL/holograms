clear 
close all

% AUTHOR:  Filipe M. Ferreira (filipef@ieee.org)

% Based on:
% Carmelo Rosales-Guzmán, Andrew Forbes,
% "How to Shape Light with Spatial Light Modulators," https://doi.org/10.1117/3.2281295

disp("by far the best algorithm for intensity search")

en = [0     6    12    18    24    30    36    42    48    54    60    66    72    78    84    90    96]*1e4;
ep = [0.1977    0.0932    0.0726    0.0619    0.0559    0.0517    0.0485    0.0462    0.0440    0.0422    0.0409    0.0396    0.0384    0.0375    0.0367    0.0360    0.0352]; % 64
es = [0.4242    0.7526    0.8303    0.8713    0.8869    0.9034    0.9142    0.9199    0.9270    0.9336    0.9376    0.9409    0.9445    0.9465    0.9483    0.9507    0.9522]; % 64

downSamp = 8; % power of 2 only

%% Load target
rawData1 = importdata('mandrill.tiff');
mandrill = double(rgb2gray(rawData1));
mandrill = mandrill./256;
mandrill = mandrill(1:downSamp:512,1:downSamp:512);
T=single(mandrill);
T=T./max(max(T));
avgT=mean(mean(T));
sz = size(T);

%% Find Hologram
H    = (ifft2((T)));
H    = (floor((angle(H) + pi)/2/pi*2^8)/2^8*2*pi-pi); % 8-bits SLM
R    = 1/sz(1)*(fft2((exp(1j.*(H)))));
avgR = mean(mean(abs(R)));
R    = (R./avgR).*avgT;
rmse = mean((abs(R(:))-T(:)).^2).^0.5;

T1   = T.*exp(1j*angle(R)); % new target

perf(1,:) = [1 rmse ssim(abs(R)./max(abs(R(:))),T)]; % performance

ii = 2;
for n = 0:99; % iterations to optimize the phase hologram
    H    = (ifft2((T1)));
    H    = (floor((angle(H) + pi)/2/pi*2^8)/2^8*2*pi-pi); % 8-bits SLM
    R    = 1/sz(1)*(fft2((exp(1j.*(H)))));
    avgR = mean(mean(abs(R)));
    R    = (R./avgR).*avgT;
    rmse = mean((abs(R(:))-T(:)).^2).^0.5;
    
    % new target
    T1   = T.*exp(1j*angle(R));

    if 1%rem(n,1)
        perf(ii,:) = [ii rmse ssim(abs(R)./max(abs(R(:))),T)]; % performance
    
        figure(123)
        subplot(1,3,1)
        semilogy(perf(1:ii,1),perf(1:ii,2:3), '.-')%,en,ep,'b--',en,es,'r--')
        axis([0 perf(ii,1) 10^-3 1])
        subplot(1,3,2)
        imshow(abs(R./max(max(abs(R)))));%imshow([abs(R./max(max(abs(R))));T]);
        title(['Reconst image SSIM=',num2str(ssim(abs(R)./max(abs(R(:))),T),4)])
        subplot(1,3,3)
        imshow([mod(angle(R),2*pi)/2/pi]);
        title(['replay phase'])
        drawnow
        tic    
        
        ii = ii + 1;
    end
end
hold off


    
