clear 
% close all

% AUTHOR:  Filipe M. Ferreira (filipef@ieee.org)

% Based on: 
% Peter J. Christopher, Ralf Mouthaan, Miguel El Guendy, Timothy D. Wilkinson,
% "Linear-time algorithm for phase-sensitive holography," Opt. Eng. 59(8) 085104 (28 August 2020) https://doi.org/10.1117/1.OE.59.8.085104

en = [0     6    12    18    24    30    36    42    48    54    60    66    72    78    84    90    96]*1e4;
ep = [0.1977    0.0932    0.0726    0.0619    0.0559    0.0517    0.0485    0.0462    0.0440    0.0422    0.0409    0.0396    0.0384    0.0375    0.0367    0.0360    0.0352]; % 64
es = [0.4242    0.7526    0.8303    0.8713    0.8869    0.9034    0.9142    0.9199    0.9270    0.9336    0.9376    0.9409    0.9445    0.9465    0.9483    0.9507    0.9522]; % 64

thrP = 45; % number of points of the Fourier transform of the mask to be used
method = 'fastest';
downSamp = 16; % power of 2 only

%% Load target
rawData1 = importdata('mandrill.tiff');
mandrill = double(rgb2gray(rawData1));
mandrill = mandrill./256;
mandrill = mandrill(1:downSamp:512,1:downSamp:512);
mandrill = double(mandrill./max(max(mandrill)));
TarAmp=(mandrill);
our_data = TarAmp;
our_data_pad = double(rand(2*size(our_data)));
our_data_pad(size(our_data,1)/2+1:1.5*size(our_data),size(our_data,1)/2+1:1.5*size(our_data)) = our_data;
TarAmp = our_data_pad;
TarAmp=TarAmp./max(max(TarAmp));
avgTarAmp=mean(mean(TarAmp));
avgT = avgTarAmp;
sz = size(TarAmp);

rawData1 = importdata('peppers.tiff');
peppers = double(rgb2gray(rawData1));
peppers = peppers./256;
peppers = peppers(1:downSamp:512,1:downSamp:512);
peppers = double(peppers./max(max(peppers)));
TarPha=(peppers);
our_data = TarPha;
our_data_pad = double(rand(2*size(our_data)));
our_data_pad(size(our_data,1)/2+1:1.5*size(our_data),size(our_data,1)/2+1:1.5*size(our_data)) = our_data;
TarPha = our_data_pad;
TarPha=TarPha./max(max(TarPha));
avgTarPha=mean(mean(TarPha));

mask = zeros(size(TarPha));
mask(size(mandrill,1)/2+1:1.5*size(mandrill),size(mandrill,1)/2+1:1.5*size(mandrill)) = 1;
mask_inds = find(mask == 1);

%% Find Hologram
T = abs(TarAmp).*exp(1i*(TarPha)*2*pi);
H = ifft2(T);%phase_quantisation(ifft2(T));
H = angle(H);
H = exp(1j.*(floor((H + pi)/2/pi*2^8)/2^8*2*pi-pi));

R = 1/sz(1)*fft2(H);

avgR = mean(abs(R(mask_inds)));
avgT = mean(abs(T(mask_inds)));
rmse = mean( mask(:).*abs( R(:) - T(:) ).^2 );%/avgR*avgT

output    = R(size(our_data,1)/2+1:1.5*size(our_data),size(our_data,1)/2+1:1.5*size(our_data));
perf(1,:) = [0 rmse ssim(abs(output)./max(abs(output(:))),mandrill) ssim(mod(angle(output),2*pi)/2/pi,peppers)];%];ssim(abs(output)./max(abs(output(:))),mandrill./max(abs(mandrill(:)))) ssim(mod(angle(output),2*pi)/2/pi,peppers./max(abs(peppers(:))))];%];

%% Mask and Fourier pre-calculation and Mask Fourier Simplification to thrP
M  = double(mask);
L  = 1/sz(1)*ifft2(M);
F  = sz(1)*ifft2(M.*T);
K  = sz(1)*ifft2(M.*R);

[~ ,Linds] = maxk(abs(L(:)),thrP);
Lz = zeros(size(L));
Lz(Linds) = L(Linds);
[LindsX,LindsY] = ind2sub(size(L),Linds);
LindsX = LindsX-1;
LindsY = LindsY-1;
Lzi = Lz(Linds);

%%
tic
ii = 2;
for n=1:1e9 % iterations to optimize the phase hologram
    %% Select random pixel
    rp = ceil(rand(1,1)*sz(1)^2);
    % [rp_x,rp_y] = ind2sub(sz,rp); rp_x = rp_x - 1; rp_y = rp_y - 1;
    nn = (floor(rand(1,1)*2^8)/2^8);
    np = abs(H(rp))*exp(1i*nn*2*pi); % 8-bits SLM
    if abs(phase(np)-phase(H(rp))) < 1e-10 % avoid selecting the same pixel
        np = abs(H(rp))*exp(1i*mod(nn+1,2^8)*2*pi);
    end   
    rp_x = rem(rp-1,sz(1));
    rp_y = (rp-(rp_x+1))/sz(1);
  
    switch method%'naive'%method
        case 'naive'
            Ht = H;
            Ht(rp) = np;
            Rt   = 1/sz(1)*fft2(Ht);

            avgRt = mean(abs(Rt(mask_inds)));
            rmset_naive = mean( mask(:).*abs( Rt(:) - T(:) ).^2 ); % /avgRt*avgT
            delta_actual = rmset_naive - rmse;
                   
            % check error
            if rmset_naive < rmse
                H(rp) = np;
                rmse = rmset_naive;
                R = Rt;
            end
        case 'fastest' % error calculated in the hologram plane - no FFTs required to evaluate the errror
            nx = (mod(LindsX+rp_x,sz(1)));
            ny = (mod(LindsY+rp_y,sz(1)));
            dinds = (nx+1)+(ny)*sz(1);

            Fir = real(F(dinds));
            Fii = imag(F(dinds));
            Kir = real(K(dinds));
            Kii = imag(K(dinds));
            dKi = Lzi*(np-H(rp))*sz(1);
            dKir = real(dKi);
            dKii = imag(dKi);
            deltaE = sum(2*(-Fir).*(dKir)+2*(-Fii).*(dKii) + 2*(dKir).*(Kir) + 2*(dKii).*(Kii) + (dKir).^2 + (dKii).^2)/sz(1)^2;
            %rmset_fastest = rmse + deltaE;
            %((rmset_fastest - rmse)/delta_actual-1)*100   

            if deltaE < 0
                H(rp) = np;
                rmse = rmse + deltaE;
                K(dinds) = K(dinds) + dKi;
            end
        case 'fast' % actually slower than naive - just an intermediate step between 'naive' and 'fastest'
            Ht = H;
            Ht(rp) = np;
            Rt = 1/sz(1)*fft2(Ht);
            Kt = sz(1)*ifft2(M.*Rt);
            rmset_fastnaive = mean(mean( F.*conj(F) - F.*conj(Kt) - conj(F).*Kt + Kt.*conj(Kt) ));
            % ((rmset_fastnaive - rmse)/delta_actual-1)*100
            
            % check error
            if rmset_fastnaive < rmse
                H(rp) = np;
                rmse = rmset_fastnaive;
            end
        case 'vfast' % actually no much faster than naive - just a step closer to 'fastest'
            deltaK = circshift(Lz,[rp_x rp_y])*(np-H(rp))*sz(1);
            deltaE_fast = mean(mean( -F.*conj(deltaK) - conj(F).*deltaK + deltaK.*conj(K) + K.*conj(deltaK) + deltaK.*conj(deltaK) ));%*sz(1)^4;
            %rmset_fast = rmse + deltaE_fast;
            %((rmset_fast - rmse)/delta_actual-1)*100   
            
            if deltaE_fast < 0
                H(rp) = np;
                rmse = rmse + deltaE_fast;
                K = K + deltaK;
            end
    end
    
    if rem(n,1e3) == 0 % hologram not 100% efficient, obviously - boost the amplitude boost of the incident beam  
        R   = 1/sz(1)*fft2(H);
        avgT = mean(abs(T(mask_inds)));
        avgR = mean(abs(R(mask_inds)));
        
        Hg   = sz(1)*ifft2(R * 1/avgR*avgT);
        Rg   = 1/sz(1)*fft2(Hg);
        
        K = sz(1)*ifft2(M.*Rg);
        H = Hg;
        
        if strcmp(method,'naive') || strcmp(method,'fast') 
            R = Rg;
            rmse = mean( mask(:).*abs( R(:) - T(:) ).^2 ); 
        end
    end
    
    % Plot progress
    if rem(n,1e5) == 0
        R   = 1/sz(1)*fft2(H);
        output = R(size(our_data,1)/2+1:1.5*size(our_data),size(our_data,1)/2+1:1.5*size(our_data));
        perf(ii,:) = [n rmse ssim(abs(output)./max(abs(output(:))),mandrill) ssim(mod(angle(output),2*pi)/2/pi,peppers)];%];       

        toc
        figure(123)
        subplot(1,3,1)
        semilogy(perf(1:ii,1),perf(1:ii,2:4), '.-',en,es,'r--')%,en,ep,'b--'
        axis([0 perf(ii,1) 10^-2 1])
        hold on
        subplot(1,3,2)
        imshow([mod(abs(R)/max(abs(R(mask_inds))),1);TarAmp]);
        title(['monke SSIM=',num2str(perf(ii,3),4)])
        subplot(1,3,3)
        imshow([mod(angle(R),2*pi)/2/pi;TarPha]);
        title(['pepe SSIM=',num2str(perf(ii,4),4)])
        drawnow
        tic  
        
        if perf(ii,3) > 0.99 || perf(ii,4) > 0.99
            break;
        end
        
        ii = ii + 1;
    end
end
ii = ii - 1;
hold off
