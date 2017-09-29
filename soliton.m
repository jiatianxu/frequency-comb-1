% kdv.m - Solve KdV equation by Fourier spectral/ETDRK4 scheme
%         A.-K. Kassam and L. N. Trefethen 4/03
%
% This code solves the Korteweg-de Vries eq. u_t+uu_x+u_xxx=0
% with periodic BCs on [-pi,pi] and initial condition given by
% a pair of solitons.  The curve evolves up to t=0.005 and at
% the end u(x=0) is printed to 6-digit accuracy.  Changing N
% to 384 and h to 2.5e-7 improves this to 10 digits but takes
% four times longer.

% Set up grid and two-soliton initial data:
  N = 256; x = (2*pi/N)*(-N/2:N/2-1)';
  A = 25; B = 16;
  u = 3*A^2*sech(.5*(A*(x+2))).^2+3*B^2*sech(.5*(B*(x+1))).^2;
  p = plot(x,u,'linewidth',3);
  axis([-pi pi -200 2200]), grid on

% Precompute ETDRK4 scalar quantities (Kassam-Trefethen):
  h = 1e-6;                               % time step
  k = [0:N/2-1 0 -N/2+1:-1]';             % wave numbers
  L = 1i*k.^3;                            % Fourier multipliers
  E = exp(h*L); E2 = exp(h*L/2);
  M = 64;                                 % no. pts for complex means
  r = exp(2i*pi*((1:M)-0.5)/M);           % roots of unity
  LR = h*L(:,ones(M,1))+r(ones(N,1),:);
  Q  = h*mean(                  (exp(LR/2)-1)./LR   ,2);
  f1 = h*mean((-4-LR+exp(LR).*(4-3*LR+LR.^2))./LR.^3,2);
  f2 = h*mean(    (4+2*LR+exp(LR).*(-4+2*LR))./LR.^3,2);
  f3 = h*mean((-4-3*LR-LR.^2+exp(LR).*(4-LR))./LR.^3,2);
  g = -.5i*k;

% Time-stepping by ETDRK4 formula (Cox-Matthews):
  set(gcf,'doublebuffer','on')
  disp('press <return> to begin'), pause  % wait for user input
  t = 0; step = 0; v = fft(u);
  while t+h/2 < .005
    step = step+1;
    t = t+h;
    Nv = g.*fft(real(ifft(v)).^2);
    a = E2.*v+Q.*Nv;        Na = g.*fft(real(ifft(a)).^2);
    b = E2.*v+Q.*Na;        Nb = g.*fft(real(ifft(b)).^2);
    c = E2.*a+Q.*(2*Nb-Nv); Nc = g.*fft(real(ifft(c)).^2);
    v = E.*v+(Nv.*f1+(Na+Nb).*f2+Nc.*f3);
    if mod(step,25)==0
      u = real(ifft(v)); set(p,'ydata',u)
      title(sprintf('t = %7.5f',t),'fontsize',18), drawnow
    end
  end
  text(-2.4,900,sprintf('u(0) = %11.7f',u(N/2+1)),...
               'fontsize',18,'color','r')

