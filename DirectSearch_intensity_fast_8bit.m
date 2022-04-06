clear 
close all

% AUTHOR:  Filipe M. Ferreira (filipef@ieee.org)

% Based on:
% Peter J. Christopher, Jamie D. Lake, Daoming Dong, Hannah J. Joyce, and Timothy D. Wilkinson,
% "Improving holographic search algorithms using sorted pixel selection," J. Opt. Soc. Am. A 36, 1456-1462 (2019)

%% maximum expectable convergence performance - for direct search
en = [0     6    12    18    24    30    36    42    48    54    60    66    72    78    84    90    96]*1e4;
ep = [0.1977    0.0932    0.0726    0.0619    0.0559    0.0517    0.0485    0.0462    0.0440    0.0422    0.0409    0.0396    0.0384    0.0375    0.0367    0.0360    0.0352]; % 64
es = [0.4242    0.7526    0.8303    0.8713    0.8869    0.9034    0.9142    0.9199    0.9270    0.9336    0.9376    0.9409    0.9445    0.9465    0.9483    0.9507    0.9522]; % 64

%% downsampling for faster testing
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
Tr = abs(T).*exp(1i*(rand(size(T))-0.5)*2*pi); % randomise the phase 
H  = angle(ifft2(Tr));
H  = 1/sz(1)*exp(1j.*(floor((H + pi)/2/pi*2^8)/2^8*2*pi-pi)); % 8-bits SLM
R  = (fft2((H)));
avgR = mean(mean(abs(R)));
erro = mean((abs(R(:))-T(:)).^2);%/avgR*avgT
    
perf(1,:) = [0 erro.^0.5 ssim(abs(R)./max(abs(R(:))),T)]; % performance

%% start search
ii = 2;
[yy,xx] = meshgrid(0:sz(1)-1,0:sz(2)-1);
tic
for n=1:1e6 % iterations to optimize the phase hologram
    %% Select random pixel
    rp = ceil(rand(1,1)*sz(1)^2);
    rp_x = rem(rp-1,sz(1))+1;
    rp_y = (rp-rp_x)/sz(1) + 1;
    
    nn = (floor(rand(1,1)*2^8)/2^8);
    np = abs(H(rp))*exp(1i*(nn)*2*pi); % 8-bits SLM
%     if abs((np)-(H(rp))) < 1e-10 % avoid selecting the same pixel
%         np = abs(H(rp))*exp(1i*(mod(nn+1,2^8))*2*pi);
%     end 
    
    %% Calculate Error  
    Rt = R+((+(np)-(H(rp))).*exp( -1j*2*pi/(sz(1)) * (yy*(rp_y-1)+xx*(rp_x-1)))); % just update the replay for the pixel that changed
    % Ht = H; Ht(rp) = np;
    % Rt   = fft2(Ht); % re-calculate the whole replay

    errot = mean((abs(Rt(:))-T(:)).^2);

    %% check error
    if errot < erro
        H(rp) = np;
        erro = errot;
        R = Rt;
    end
    
    if rem(n,1e3) == 0 % hologram not 100% efficient, obviously - boost the amplitude boost of the incident beam  
        R   = fft2(H);
        avgR = mean(abs(R(:)));
        
        H   = ifft2(R * 1/avgR*avgT);
        R   = fft2(H);
        erro = mean((abs(R(:))-T(:)).^2);
    end

    if rem(n,10000) == 0
        R   = fft2(H);        
        perf(ii,:) = [n erro.^0.5 ssim(abs(R)./max(abs(R(:))),T)];  % performance
        
        toc
        figure(123)
        subplot(1,3,1)
        semilogy(perf(1:ii,1),perf(1:ii,2:3), '.-',en,ep,'b--',en,es,'r--')
        axis([0 perf(ii,1) min(ep) max(es)])
        ylabel('error metric')
        xlabel('#iterations')
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

