function A=DeconPlotFun(peakshape,x,pos,wid,amp,extra,Delta)
switch peakshape
    case 1
        A=gaussian(x,pos,wid);
    case 2
        A=lorentzian(x,pos,wid);
    case 3
        A=logistic(x,pos,wid);
    case 4
        A=pearson(x,pos,wid,extra);
    case 5
        A=expgaussian(x,pos,wid,-extra)';
    case 6
        A=gaussian(x,pos,wid);
    case 7
        A=lorentzian(x,pos,wid);
    case 8
        A=expgaussian(x,pos,wid,-extra)';
    case 9
        A=exppulse(x,pos,wid);
    case 10
        A=upsigmoid(x,pos,wid);
    case 11
        A=gaussian(x,pos,wid);
    case 12
        A=lorentzian(x,pos,wid);
    case 13
        A=GL(x,pos,wid,extra);
    case 14
        A=BiGaussian(x,pos,wid,extra);
    case 15
        A=BWF(x,pos,wid,extra);
    case 16
        A=gaussian(x,pos,wid);
    case 17
        A=lorentzian(x,pos,wid);
    case 18
        A=explorentzian(x,pos,wid,-extra)';
    case 19
        A=alphafunction(x,pos,wid);
    case 20
        A=voigt(x,pos,wid,extra);
    case 21
        A=triangular(x,pos,wid);
    case 22
        A=[];%peakfunction(shapesvector,x,pos,wid,extra);
    case 23
        A=downsigmoid(x,pos,wid);
    case 24
        A=nbinpdf(x,pos,wid);
    case 25
        A=lognormal(x,pos,wid);
    case 26
        A=linslope(x,pos,wid);
    case 27
        A=d1gauss(x,pos,wid);
    case 28
        A=polynomial(x,pos,wid);
    case 29
        A=segmented(x,yy,PEAKHEIGHTS);
    case 30
        A=voigt(x,pos,wid,extra);
    case 31
        A=expgaussian(x,pos,wid,extra);
    case 32
        A=pearson(x,pos,wid,extra);
    case 33
        A=GL(x,pos,wid,extra);
    case 34

        A=voigt(x,pos, abs(wid),extra);
    case 35
        A=GL(x,pos,wid,extra);
    case 36
        A=expgaussian(x,pos,wid,extra);
    case 37
        A=pearson(x,pos,wid,extra);
    case 38
        A=explorentzian(x,pos,wid,extra);
    case 39
        A=ExpGaussian2(x,pos,wid,extra);
    case 40
        A=GL(x,pos,wid);
    case 41
        A=rectangle(x,pos,wid);
    case 42
        A=ngaussian(x,pos,wid,extra);
    case 43
        A=Gompertz(x,pos,wid,extra);
    case 44
        A=OneMinusExp(x,pos,wid);
    case 45
        A=FourPL(x,pos,wid,extra);
    case 46
        A=quadslope(x,pos,wid);
    case 47
        A=blackbody(x,pos);
    case 48
        A=exppulse(x,pos,wid);
    case 49
        A=[];
        A=PearsonIV(x,pos,wid,extra,Delta);
end % switch
% for

A=amp.*A;
% ----------------------------------------------------------------------
function g = gaussian(x,pos,wid)
%  gaussian(X,pos,wid) = gaussian peak centered on pos, half-width=wid
%  X may be scalar, vector, or matrix, pos and wid both scalar
% Examples: gaussian([0 1 2],1,2) gives result [0.5000    1.0000    0.5000]
% plot(gaussian([1:100],50,20)) displays gaussian band centered at 50 with width 20.
g = exp(-((x-pos)./(0.60056120439323.*wid)).^2);


% ----------------------------------------------------------------------
function g = ngaussian(x,pos,wid,n)
%  ngaussian(x,pos,wid) = flattened Gaussian centered on x=pos, half-width=wid
%  x may be scalar, vector, or matrix, pos and wid both scalar
% Shape is Gaussian when n=1. Becomes more rectangular as n increases.
%  T. C. O'Haver, 1988, revised 2014
% Example: ngaussian([1 2 3],1,2,1) gives result [1.0000    0.5000    0.0625]
if n>0
    g = 1-(10.^-(n.*gaussian(x,pos,wid)));
    g=g./max(g);
else
    g = gaussian(x,pos,wid);
end
% ----------------------------------------------------------------------

% ----------------------------------------------------------------------
function g = lorentzian(x,position,width)
% lorentzian(x,position,width) Lorentzian function.
% where x may be scalar, vector, or matrix
% position and width scalar
% T. C. O'Haver, 1988
% Example: lorentzian([1 2 3],2,2) gives result [0.5 1 0.5]
g=ones(size(x))./(1+((x-position)./(0.5.*width)).^2);

% ----------------------------------------------------------------------
function g = logistic(x,pos,wid)
% logistic function.  pos=position; wid=half-width (both scalar)
% logistic(x,pos,wid), where x may be scalar, vector, or matrix
% pos=position; wid=half-width (both scalar)
% T. C. O'Haver, 1991
n = exp(-((x-pos)/(.477.*wid)) .^2);
g = (2.*n)./(1+n);

% ----------------------------------------------------------------------
function g = triangular(x,pos,wid)
%Triangle function.  pos=position; wid=half-width (both scalar)
%trianglar(x,pos,wid), where x may be scalar or vector,
%pos=position; wid=half-width (both scalar)
% T. C. O'Haver, 1991
% Example
% x=[0:.1:10];plot(x,trianglar(x,5.5,2.3),'.')
g=1-(1./wid) .*abs(x-pos);
for i=1:length(x)
    if g(i)<0,g(i)=0;end
end
% ----------------------------------------------------------------------

% ----------------------------------------------------------------------
function g = rectangle(x,pos,wid)
%rectangle function.  pos=position; wid=half-width (both scalar)
%rectangle(x,pos,wid), where x may be scalar or vector,
%pos=position; wid=half-width (both scalar)
% T. C. O'Haver, 2016
% Example
% x=[0:.1:10];plot(x,rectangle(x,5.5,2.3),'.')
g=zeros(size(x));
hw=wid./2;
for i=1:length(x)
    if x(i)<pos-hw,g(i)=0;end
    if x(i)>pos-hw,g(i)=1;end
    if x(i)>pos+hw,g(i)=0;end
end
% ----------------------------------------------------------------------
function g = pearson(x,pos,wid,m)
% Pearson VII function.
% g = pearson(x,pos,wid,m) where x may be scalar, vector, or matrix
% pos=position; wid=half-width (both scalar)
% m=some number
%  T. C. O'Haver, 1990
g=ones(size(x))./(1+((x-pos)./((0.5.^(2/m)).*wid)).^2).^m;


% ----------------------------------------------------------------------
function g = expgaussian(x,pos,wid,timeconstant)
%  Exponentially-convoluted gaussian(x,pos,wid) = gaussian peak centered on pos, half-width=wid
%  x may be scalar, vector, or matrix, pos and wid both scalar
%  T. C. O'Haver, 2006
g = exp(-((x-pos)./(0.60056120439323.*wid)) .^2);
g = ExpBroaden(g',timeconstant);
% ----------------------------------------------------------------------
function emg=ExpGaussian2(t,mu,s,lambda)
%   Fitting functions for multiple exponentially-modified Gaussian bands signal.
%  T. C. O'Haver, May 21, 2017.
% Same shape as expgaussian except parameterized differently
%if tau~=0,
%    lambda=1/tau;
% disp([mu,s,lambda]) % <<<<<<<<<<<<<<<<<<<<
EMG=s.*lambda.*sqrt(pi/2).*exp(0.5.*(s.*lambda).^2-lambda.*(t-mu)).*erfc((1/sqrt(2)).*(s.*lambda-((t-mu)./s)));
emg=EMG./max(EMG);
%end
% ----------------------------------------------------------------------
function g = explorentzian(x,pos,wid,timeconstant)
%  Exponentially-broadened lorentzian(x,pos,wid) = lorentzian peak centered on pos, half-width=wid
%  x may be scalar, vector, or matrix, pos and wid both scalar
%  T. C. O'Haver, 2013
g = ones(size(x))./(1+((x-pos)./(0.5.*wid)).^2);
g = ExpBroaden(g',timeconstant);
% ----------------------------------------------------------------------
function yb = ExpBroaden(y,t)
% ExpBroaden(y,t) zero pads y and convolutes result by an exponential decay
% of time constant t by multiplying Fourier transforms and inverse
% transforming the result.
hly=round(length(y)./2);
ey=[y(1).*ones(1,hly)';y;y(length(y)).*ones(1,hly)'];
fy=fft(ey);
a=exp(-(1:length(fy))./t);
fa=fft(a);
fy1=fy.*fa';
ybz=real(ifft(fy1))./sum(a);
yb=ybz(hly+2:length(ybz)-hly+1);
% ----------------------------------------------------------------------

% ----------------------------------------------------------------------
function g = exppulse(x,t1,t2)
% Exponential pulse of the form
% g = (x-spoint)./pos.*exp(1-(x-spoint)./pos);
e=(x-t1)./t2;
p = 4*exp(-e).*(1-exp(-e));
p=p .* (p>0);
g = p';
% ----------------------------------------------------------------------


% ----------------------------------------------------------------------
function g = alphafunction(x,pos,spoint)
% alpha function.  pos=position; wid=half-width (both scalar)
% alphafunction(x,pos,wid), where x may be scalar, vector, or matrix
% pos=position; wid=half-width (both scalar)
% Taekyung Kwon, July 2013
g = (x-spoint)./pos.*exp(1-(x-spoint)./pos);
for m=1:length(x);if g(m)<0;g(m)=0;end;end
% ----------------------------------------------------------------------


% ----------------------------------------------------------------------
function g=downsigmoid(x,t1,t2)
% down step sigmoid
g=.5-.5*erf(real((x-t1)/sqrt(2*t2)));
% ----------------------------------------------------------------------
function g=upsigmoid(x,t1,t2)
% up step sigmoid
g=1/2 + 1/2* erf(real((x-t1)/sqrt(2*t2)));
% ----------------------------------------------------------------------


% ----------------------------------------------------------------------
function g = GL(x,pos,wid,m)
% Gaussian/Lorentzian blend. m = percent Gaussian character
% pos=position; wid=half-width
% m = percent Gaussian character.
%  T. C. O'Haver, 2012
% sizex=size(x)
% sizepos=size(pos)
% sizewid=size(wid)
% sizem=size(m)
g=2.*((m/100).*gaussian(x,pos,wid)+(1-(m(1)/100)).*lorentzian(x,pos,wid))/2;
% ----------------------------------------------------------------------


% ----------------------------------------------------------------------
function v=voigt(x,pos,Gausswidth,voigtalpha)
% Unit height Voigt profile function. x is the independent variable
% (energy, wavelength, etc), Gausswidth is the Gaussian(Doppler) width,
% and voigtalpha is ratio of the Gausswidth to the LorentzWidth (pressure
% width). Version 3, August, 2019
LorentzWidth=Gausswidth.*voigtalpha;
if LorentzWidth<0, LorentzWidth=-LorentzWidth;end
if Gausswidth<0, Gausswidth=-Gausswidth;end
dx=x(2)-x(1); % x increment
ex=[x-max(x)-dx x x+max(x)]; % Extended x
gau=gaussian(ex,0,Gausswidth);
lor=lorentzian(ex,pos,LorentzWidth);
VoigtConv=ifft(fft(gau).*fft(lor))./sum(lor);
g=VoigtConv./max(VoigtConv);
oex=ex-max(x);
outrange=val2ind(oex,0):val2ind(oex,max(x))-dx;
% minoutrange=min(outrange)
% maxoutrange=max(outrange)
v=g(outrange+1);
% ----------------------------------------------------------------------


% ----------------------------------------------------------------------
function g = BiGaussian(x,pos,wid,m)
% BiGaussian (different width on leading edge and trailing edge).
% pos=position; wid=width
% m = ratio of widths of right-hand to left-hand halves.
% If m=1, it becomes identical to a Gaussian.
% Verison 2, T. C. O'Haver, 2012
% Example: Plots Gaussian and BiGaussian with m=3, both with halfwidth=20.
% x=[1:100];
% g=gaussian(x,50,20);
% bg=bigaussian(x,50,20,3);
% plot(x,g,x,bg)
%
lx=length(x);
hx=val2ind(x,pos);
hwid=2.*wid;
g(1:hx)=gaussian(x(1:hx),pos,hwid/(m+1));
g(hx+1:lx)=gaussian(x(hx+1:lx),pos,m.*hwid/(m+1));
% ----------------------------------------------------------------------


% ----------------------------------------------------------------------
function g = BWF(x,pos,wid,m)
% BWF (Breit-Wigner-Fano) http://en.wikipedia.org/wiki/Fano_resonance
% pos=position; wid=width; m=Fano factor
%  T. C. O'Haver, 2014
y=((m*wid/2+x-pos).^2)./(((wid/2).^2)+(x-pos).^2);
% y=((1+(x-pos./(m.*wid))).^2)./(1+((x-pos)./wid).^2);
g=y./max(y);
% ----------------------------------------------------------------------

function g = lognormal(x,pos,wid)
% lognormal function.  pos=position; wid=half-width (both scalar)
% lognormal(x,pos,wid), where x may be scalar, vector, or matrix
% pos=position; wid=half-width (both scalar)
% T. C. O'Haver, 1991
g = exp(-(log(x/pos)/(0.01.*wid)) .^2);
% ----------------------------------------------------------------------
function err = fitewGLblend(lambda,t,y,extra)
% Fitting function for multiple Gaussian peaks with equal peak widths.
global PEAKHEIGHTS BaselineMode BIPOLAR LOGPLOT
numpeaks=round(length(lambda)-1);
A = zeros(length(t),numpeaks);
for j = 1:numpeaks
    A(:,j) = GL(t,lambda(j),lambda(numpeaks+1),extra)';
end
if BaselineMode==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
% function g=sine(t,f,phase)
% % Sine wave (alpha test)
% g=sin(2*pi*f*(t+phase));
% % ----------------------------------------------------------------------


% ----------------------------------------------------------------------
function y=d1gauss(x,p,w)
% First derivative of Gaussian (alpha test)
y=-(5.54518.*(x-p).*exp(-(2.77259.*(p-x).^2)./w^2))./w.^2;
y=y./max(y);
% ----------------------------------------------------------------------

function y=polynomial(t,coeff)
y=polyval(coeff,t);
% ----------------------------------------------------------------------

function yi=segmented(x,y,segs)
global PEAKHEIGHTS
clear yy
for n=1:length(segs)
    yind=val2ind(x,segs(n));
    yy(n)=y(yind(1));
end
yi=interp1(segs,yy,x);
PEAKHEIGHTS=segs;
% ----------------------------------------------------------------------


% ----------------------------------------------------------------------
function y=linslope(x,slope,intercept)
y=x.*slope+intercept;
% y=y./max(y);
% ----------------------------------------------------------------------


% ----------------------------------------------------------------------
function y=quadslope(x,boa,coa) % normalized quadratic
y=(x.^2+(boa).*x+coa);
y=y./max(y);
% ----------------------------------------------------------------------

% ----------------------------------------------------------------------
function y = PearsonIV(x,a1,a2,a3,a4)
% Unit-height Pearson IV probability distrubtion function as described by
% the list of functions from PeakFit(TM) software documentation on page 261.
% Suggested by Chris George. Code by T. C. O'Haver, 2021
% Example:
% x=1:.1:10;y=PearsonIV(x,6,3,4,2);plot(x,y)
%
num=(x-((a2.*a4)./(2.*a3))-a1);
y=(((1+num.^2./a2.^2).^-a3) .* exp(-a4.*(atan(num./a2)+atan(a4./(2.*a3))))) ./ (1+((a4.^2)./(4.*a3.^2))).^-a3;
y=y./max(y);
% ----------------------------------------------------------------------
function y=Gompertz(t,Bo,Kh,L)
% A Gompertz curve or Gompertz function, named after Benjamin Gompertz, is
% a sigmoid function. It is a type of mathematical model for a time series,
% where growth is slowest at the start and end of a time period. The
% right-hand or future value asymptote of the function is approached much
% more gradually by the curve than the left-hand or lower valued asymptote,
% in contrast to the simple logistic function in which both asymptotes are
% approached by the curve symmetrically. It is a special case of the
% generalized logistic function.
%
% Example:
% x=1:.1:10;y=gompertz(x,6,3,4);plot(x,y)
%
y=Bo*exp(-exp((Kh*exp(1)/Bo)*(L-t) +1));

