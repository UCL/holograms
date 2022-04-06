clear
close all

sz = [16 16];
[mm,nn] = meshgrid(0:sz(1)-1,0:sz(2)-1);

H = rand(sz)*2*pi;
H = exp(1i*H);

Hfft = zeros(size(H));
for u = 0:sz(1)-1
    for v = 0:sz(1)-1
        Hfft(u+1,v+1) = sum(sum( H.*exp( -1j*2*pi/(sz(1)) * (mm*(v)+nn*(u))) ));
    end
end

A = Hfft;
B = (fft2(H));
figure()
subplot(3,2,1)
surf(abs(A),'EdgeColor','none'),view([0 0 1])
subplot(3,2,2)
surf(abs(B),'EdgeColor','none'),view([0 0 1])
subplot(3,2,3)
surf(abs(A./B)-1,'EdgeColor','none')%,view([0 0 1])
subplot(3,2,4)
surf(angle(A./B),'EdgeColor','none')%,view([1 1 0])
subplot(3,2,5)
surf((angle(A))-(angle(B)),'EdgeColor','none')%,view([0 0 1])
subplot(3,2,6)
surf(unwrap(angle(B)),'EdgeColor','none')%,view([0 0 1])

%% IFFT

sz = [16 16];
[mm,nn] = meshgrid(0:sz(1)-1,0:sz(2)-1);

H = B;

Hifft = zeros(size(H));
for u = 0:sz(1)-1
    for v = 0:sz(1)-1
        Hifft(u+1,v+1) = sum(sum( H.*exp( 1j*2*pi/(sz(1)) * (mm*(v)+nn*(u))) ));
    end
end

A = 1/sz(1)^2*Hifft;
B = (ifft2(H));
figure()
subplot(3,2,1)
surf(abs(A),'EdgeColor','none'),view([0 0 1])
subplot(3,2,2)
surf(abs(B),'EdgeColor','none'),view([0 0 1])
subplot(3,2,3)
surf(abs(A./B)-1,'EdgeColor','none')%,view([0 0 1])
subplot(3,2,4)
surf(angle(A./B),'EdgeColor','none')%,view([1 1 0])
subplot(3,2,5)
surf((angle(A))-(angle(B)),'EdgeColor','none')%,view([0 0 1])
subplot(3,2,6)
surf(unwrap(angle(B)),'EdgeColor','none')%,view([0 0 1])
