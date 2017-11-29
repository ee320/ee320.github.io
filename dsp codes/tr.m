function [classif]=tr(start,finish,filename)
start=0;
classif = [];

[Y1,FS,NBITS,OPTS]=wavread(filename); % y1 is sampled data, sample rate(FS) in Hz, bitsper sample(NBITS) 
duration = finish-start;
kkk = size(Y1,2);
arxh = start*FS;
m = 100;
timeread = 1;

Vt = 0*[1:50];Vta = [];Vtz = 0*[1:50];Vtz1 = 0*[1:50];Wtaz = [];Wtaz1 = [];shma = [];temp = [];
y = 0;
py = [1:1:50];

scp = 0;
sc = 0;
Vtp = [1:1:50];
i=1;
[PYs,FS,NBITS,OPTS]=wavread(filename,[arxh+1+(i-1)*timeread*FS  arxh+i*timeread*FS]);
sf11 = 0;
aVt = [];

%READING - Computation of RMS-ZC-RMS*ZC per 20 msec
for i=1:floor(duration/timeread),
   [Ys,FS,NBITS,OPTS]=wavread(filename,[arxh+1+(i-1)*timeread*FS  arxh+i*timeread*FS]);
   if kkk==2,
      Ys = Ys(:,1)+Ys(:,2);
   end
   for j=1:timeread*50,
      y = Ys(round((j-1)*0.02*FS+1):floor(j*0.02*FS));
      rms = sqrt(sum(y.*y));
      Vtp = Vt;
      Vt = [Vt(2:50) rms];
      l = length(y);
      y1 = y(1:l-1);
      y2 = y(2:l);
      mm1 = y1.*y2;
      zc = 0.5*sum((abs(sign(mm1))-sign(mm1)));
      Vtz1 = [Vtz1(2:50) zc];
      zc = zc*rms;
      Vtz = [Vtz(2:50) zc];
   end
   
   m = mean(Vt);
   v = var(Vt);
   Vta = [Vta Vt];
   Wtaz = [Wtaz Vtz];
   Wtaz1 = [Wtaz1 Vtz1];
   
   if (v ~= 0 & m ~= 0)
      b(i) = v/m;
      a(i) = (m/b(i))-1;
   else
      b(i) = 1000000;
      a(i) = -1;
   end
end

%Computation of  KVRMS, calculating nRMS varience
l = length(Vta);
ftm(1) = mean(Vta(1:50));
ftv(1) = var(Vta(1:50));
for i=2:l-50,
  ftm(i) = ftm(i-1)+(Vta(i+49)-Vta(i-1))/50;
  ftv(i) = ftv(i-1) + ftm(i-1)*ftm(i-1) + ((Vta(i+49)*Vta(i+49)-Vta(i-1)*Vta(i-1))/50) - ftm(i)*ftm(i);  
  
  if  ftm(i)~= 0,
     kvrms(i) =   ftv(i)/(ftm(i)*ftm(i));
  else
     kvrms(i) = 0;
  end
end

l = floor(duration/timeread);
m = mean(kvrms);
v = var(kvrms);
b1111 = v/m; % exact NRMS
a111 = (m/b1111)-1;

%computation of similarity between windows  based on RMS distribution 
for i=2:l,
   k = (a(i)+a(i-1)+2)/2;
   somiot1(i-1) = ((2/(b(i-1)+b(i)))^k)*(b(i-1)^((a(i)+1)/2))*(b(i)^((a(i-1)+1)/2))*gamma(k)/(sqrt(gamma(a(i-1)+1)*gamma(a(i)+1)));     
   if somiot1(i-1) ~= somiot1(i-1),
      somiot1(i-1) = 0;
   end 
end

for i=3:l,
   k = (a(i)+a(i-2)+2)/2;
   somiot2(i-2) = ((2/(b(i-2)+b(i)))^k)*(b(i-2)^((a(i)+1)/2))*(b(i)^((a(i-2)+1)/2))*gamma(k)/(sqrt(gamma(a(i-2)+1)*gamma(a(i)+1)));     
   if somiot2(i-2) ~= somiot2(i-2),
      somiot2(i-2) = 0;
   end
end

P = [];

%Computation of probability  P of change  (at window i)
for i=2:l-1,
   P(i-1) =(1-somiot2(i-1));
end

P = [0 P 0];
 
M = mean(P);
V1 = mean(abs(P(1:l-3)-P(2:l-2)));
V2 = mean(abs(P(1:l-4)-P(3:l-2)));

V = (V1+V2)*0.25;

V = 0.25*median(abs((P(1:l-4)-P(3:l-2))));

Pd(1) = P(1);
Pd(2) = P(2);

for i=3:l-2,
   m = mean(P([i-2:i-1 i+1:i+2]));
   m1 = max(P([i-2:i-1 i+1:i+2]));
   
   if (m < 0.1*M)
      m1 = 0.1*M;
   end
   d = (P(i)-m)/m1;
     
   if d < 0,
      d = 0;
   end
   Pd(i) = P(i)*d/V;
end
Pd(l) = P(l);
Pd(l-1) = P(l-1);

%window size = 1 sec 
%Selection of candicate windows where there is a change 
j = 0;
for i=3:l-3,
   if (Pd(i) > 5 & P(i) > P(i+1) & P(i) > P(i-1)),
      j = j+1;
      pos(j) = i;
   end
end
if j == 0,
   return;
end

%Computation of change with 20 msec accurancy 
k2 = 1;
c = 25;
for k1=1:length(pos),
   i = 50*(pos(k1)-2);  % msec
   for j=-c:50+c,
      if (ftv(i+j) ~= 0 & ftm(i+j) ~= 0),
         b1 = ftv(i+j)/ftm(i+j);
         a1 = ftm(i+j)/b1-1;
      else
         b1 = 10000;
         a1 = -1;
      end
      if (ftv(i+j+50) ~= 0 & ftm(i+j+50) ~= 0)
         b2 = ftv(i+j+50)/ftm(i+j+50);
         a2 = ftm(i+j+50)/b2-1;
      else
         b2 = 10000;
         a2 = -1;
      end
      k = (a2+a1+2)/2;
      Om(j+c+1) = ((2/(b1+b2))^k)*(b1^((a2+1)/2))*(b2^((a1+1)/2))*gamma(k)/(sqrt(gamma(a1+1)*gamma(a2+1)));     
      if Om(j+c+1)~= Om(j+c+1),
         Om(j+c+1) = 0;
      end
   end
   [x posms(k1)] = min(Om);
   h = 1-Om;
   if mean(h(max(posms(k1)-25,1):min(2*c+51,posms(k1)+25))) >  0.1,
      msec(k2) = posms(k1)-1-c;
      mpos(k2) = pos(k1);
      k2 = k2+1;
   end
end

posms = [];posms = msec;pos = [];pos = mpos;
fpos = start+pos-1+0.02*posms;
gfpos = floor(50*pos+posms-50);

bvoise = 0.14339738781836;bmusic = 0.04399754659151;amusic = 1.66349817725076;avoise = 2.32677887950291;

gfpos = [1 gfpos length(kvrms)];
l = length(gfpos)-1;
V = [];
Pzero = [];
Fsp1 = [];

%Classification for each segment 
for i=1:l,
   d = 0;
   x1 = mean(kvrms([gfpos(i):gfpos(i+1)]));
   x = x1;
   y = Wtaz([gfpos(i)+d:gfpos(i+1)-d]);
   y = y/(2*max(Vta(gfpos(i)+d:gfpos(i+1)-d)) - min(Vta(gfpos(i)+d:gfpos(i+1)-d))-median(Vta(gfpos(i)+d:gfpos(i+1)-d)));
   thor(i) = 50*mean(exp(-y));
   PPmusic(i) = (x^amusic)*exp(-x/bmusic)/(gamma(amusic+1)*bmusic^(amusic+1));
   PPvoise(i) = 0.5*(x^avoise)*exp(-x/bvoise)/(gamma(avoise+1)*bvoise^(avoise+1));
   sm = 0.7*median(ftm(gfpos(i):gfpos(i+1)))+0.3*mean(ftm(gfpos(i):gfpos(i+1)));
   Pspace(i) =  6*exp((-sm*sm)/(2*0.6*0.6));
   zpos(i) = sum(exp(-10*Wtaz1([gfpos(i)+d:gfpos(i+1)-d])))/length(Wtaz1([gfpos(i)+d:gfpos(i+1)-d]));
   
   if zpos(i) > 0.08 | x > 3
     PPvoise(i)  = 10;
  else
      PPmusic(i) = (x^amusic)*exp(-x/bmusic)/(gamma(amusic+1)*bmusic^(amusic+1));
      PPvoise(i)  = 0.5*(x^avoise)*exp(-x/bvoise)/(gamma(avoise+1)*bvoise^(avoise+1));
   end
   
   V(i,:) = [PPmusic(i)  PPvoise(i)];
   [Pmax(i) type(i)] = max(V(i,:));
      
  thor2(i) = 50*mean(exp(-4*y));
  PPmusicg(i) = 0;
   if thor2(i) > 3.5,
      PPvoiseg(i) = 10;
   else
      PPmusicg(i) = (x^amusic)*exp(-x/bmusic)/(gamma(amusic+1)*bmusic^(amusic+1));
      PPvoiseg(i)  = 0.5*(x^avoise)*exp(-x/bvoise)/(gamma(avoise+1)*bvoise^(avoise+1));
   end
   
   Vg(i,:) = [PPmusicg(i)  PPvoiseg(i)];
   [Pmaxg(i) typeg(i)] = max(Vg(i,:));
   Vk = (sign(Wtaz1([gfpos(i):gfpos(i+1)])));
   tempV = [];
   tempV = Vta([gfpos(i):gfpos(i+1)]);
   rr = max(tempV)+0.001;
   tempV1 = tempV/rr;
   Vk0 = Vk;
   
   for j = 1:length(Vk)
      if (tempV1(j) < 0.1 & tempV(j) < 1.5) | tempV(j) < 1
         Vk(j) = 0;
      end
   end
   Vk1 = Vk;
   
   for j = 1:length(Vk)
      if (tempV1(j) < 0.4 | tempV(j) < 2 | tempV(j) < mean(tempV(j))) 
         Vk1(j) = 0;
      end
   end
   
   Zca = [];
   Zca = Vk1.*(Wtaz1([gfpos(i):gfpos(i+1)]));
   Pol = (ones(1,length(Vk)) -  sign(Vk)).*Vk0;        
   VZC = Pol.*Wtaz1([gfpos(i):gfpos(i+1)]);   
   dVk = abs(Vk(2:length(Vk))-Vk(1:length(Vk)-1));
   Freq(i,1) = 50*sum(dVk)/(2*length(Vk)); %FSP%
   Freq(i,2) = length(Vk)/50;%duration in sec of the segment
   Freq(i,3) = sum(Freq(:,2)); 
   Freq(i,4) = zpos(i);
   Freq(i,5) = thor2(i);
   Freq(i,6) =   max(Zca);
   PPmusicg(i) = 0;
   PPvoiseg(i) = 0;
   if Freq(i,1) < 0.59 & Freq(i,2) > 2.5,
      PPmusicg(i) = 10;
   else
      PPmusicg(i) = (x^amusic)*exp(-x/bmusic)/(gamma(amusic+1)*bmusic^(amusic+1));
      PPvoiseg(i)  = 0.5*(x^avoise)*exp(-x/bvoise)/(gamma(avoise+1)*bvoise^(avoise+1));
      if x > 3,
         PPvoiseg(i) = PPvoiseg(i)+0.1;         
      end
      if x > 300 ,
         PPmusicg(i) = 9;
      end
      if thor2 > 3.5
         PPvoiseg(i) = 11;
      end
      if max(Zca) > 280,
         PPmusicg(i) = 11;
         max(Zca)  
      end
      if Freq(i,1) > 4.62 & Freq(i,2) > 2.5 & (thor2(i) > 0.1 | zpos(i) > 0.05),
         PPvoiseg(i) = 12;
      end
      if zpos(i) > 0.15 
         PPvoiseg(i)  = 13;
      end
  end
   Vg(i,:) = [PPmusicg(i)  PPvoiseg(i)];
   [Pmaxg(i) typeg(i)] = max(Vg(i,:));

	Fsp1 = [Fsp1 Vk];  
 
   if (Pspace(i) > 4 & thor2(i) < 1)  | (Pspace(i) > 4.5),
      type(i) = 3; %silence
      typeg(i) = 3;
   end 
   if (gfpos(i+1) - gfpos(i))/50 < 0.5 & i >= 1, % short segments
      type(i) = type(i-1);
   end
end
   g = [];g1 = [];tt = [];pt = [];pm = [];ps = [];pv = [];pmf = []; 
  
for i=1:l,
   tt = [tt type(i)*ones(1,gfpos(i+1)-gfpos(i))];
   pt = [pt Pmax(i)*ones(1,gfpos(i+1)-gfpos(i))];
   pm = [pm PPmusic(i)*ones(1,gfpos(i+1)-gfpos(i))];
   pv = [pv PPvoise(i)*ones(1,gfpos(i+1)-gfpos(i))];
   ps = [ps Pspace(i)*ones(1,gfpos(i+1)-gfpos(i))];
end

PosV = 100*sum(exp(-10*abs(tt-2*ones(1,length(tt)))))/length(tt);
PosM = 100*sum(exp(-10*abs(tt-ones(1,length(tt)))))/length(tt);

l1 = l;
l = length(tt);

l = l1;
 g = [];g1 = [];tt = [];pt = [];pm = [];ps = [];pv = [];pmf = [];

for i=1:l,
   tt = [tt typeg(i)*ones(1,gfpos(i+1)-gfpos(i))];
   pt = [pt Pmaxg(i)*ones(1,gfpos(i+1)-gfpos(i))];
   pm = [pm PPmusicg(i)*ones(1,gfpos(i+1)-gfpos(i))];
   pv = [pv PPvoiseg(i)*ones(1,gfpos(i+1)-gfpos(i))];
   ps = [ps Pspace(i)*ones(1,gfpos(i+1)-gfpos(i))];
end
l = length(tt);

PosVg = 100*sum(exp(-10*abs(tt-2*ones(1,length(tt)))))/length(tt);
PosMg = 100*sum(exp(-10*abs(tt-ones(1,length(tt)))))/length(tt);

for k = 1:length(Vta),
   if Vta(k) < 2 | Vta(k) < 0.1*max(Vta),
      Wtaz1(k) = 0;
   end
end
m = 100;
temp = [];
h1 = [];
temp = Wtaz1;
h1 = hist(temp,m);
h1(1) = 0;
h1 = h1/sum(h1);
a1 = [min(temp):(max(temp)-min(temp))/m:max(temp)-(max(temp)-min(temp))/m];

figure(2);

subplot(2,1,1);
plot([1:l]*(finish-start-1)/l+start,tt);
axis([0.5 duration 0 3.5]);
title('classification type (1 : Music, 2: Voice, 3: Silence)');
xlabel('sec');
subplot(2,1,2);
plot(start+([1:length(Vta)]*duration/length(Vta)),Vta);
axis([0 duration 0 (floor(max(Vta))+1)]);
title('RMS');
xlabel('sec');


classif = tt;


a=0; b=0; c=0; z=0;
sz=size(classif);
i=0;
for i=1:sz(2)
    if classif(i)==1
        a=a+1;
    end
    if classif(i)==2
        b=b+1;
    end
    if classif(i)==3
        c=c+1;
    end
    z=z+1;
end
histo=[a b c];
histo=100*histo./(sz(2));
bar(histo);