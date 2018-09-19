function fZ = fzeta_loop_weights(Ntheta,Nrho,betas,rhos,a,osthlarge)
krho=(-Nrho/2:Nrho/2-1);
Nthetalarge=osthlarge*Ntheta;
thsplarge=(-Nthetalarge/2:Nthetalarge/2-1)/Nthetalarge*betas;
fZ=zeros(Nthetalarge,Nrho);
h=ones(size(thsplarge));
% correcting=1+[-3 4 -1]/24;correcting(1)=2*(correcting(1)-0.5);
% correcting=1+[-23681 55688 -66109 57024 -31523 9976 -1375]/120960;correcting(1)=2*(correcting(1)-0.5);
correcting=1+[-216254335 679543284 -1412947389 2415881496 -3103579086 2939942400 -2023224114 984515304 -321455811 63253516 -5675265]/958003200;correcting(1)=2*(correcting(1)-0.5);
h(1)=h(1)*(correcting(1));
for j=2:length(correcting)
    h(j)=h(j)*correcting(j);
    h(end-j+2)=h(end-j+2)*(correcting(j));
end
for j=1:numel(krho)
    fcosa=cos(thsplarge).^(-2*pi*1i*krho(j)/rhos-1-a);
    fZ(:,j)=fftshift(fft(fftshift(h.*fcosa)));
end
fZ=fZ(end/2-Ntheta/2+1:end/2+Ntheta/2,:);
fZ=fZ*(thsplarge(2)-thsplarge(1));
%put imag to 0 for the border
fZ(1,:)=0;%real(fZ(1,:));
fZ(:,1)=0;%real(fZ(:,1));
end