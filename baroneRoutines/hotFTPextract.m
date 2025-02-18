function  [ctd,iso] = hotFTPextract(hot_n)
%
% function  [ctd,iso] = hotFTPextract(hot_n)
%
%   Extracts CTD data from the FTP site of the Hawaii Ocean Time-series
%   (ftp://mananui.soest.hawaii.edu/pub/hot/) and merges them into MATLAB
%   structure variables on depth intervals and potential density intervals.
%   The routine computes the value of several oceanographic variables
%   in part using the Thermodynamic Equation of Seawater - 2010 
%   (http://www.teos-10.org/index.htm).
%
%   This function requires the Gibbs-SeaWater (GSW) Oceanographic Toolbox
%   available from the website http://www.teos-10.org/software.htm
%
%   This function also requires the function RunMedian.m to compute the
%   running median fluorescence profile and find the depth of the DCM 
%
% INPUT:
%   hot_n: number of the HOT cruise desired
%
% OUTPUT:
%   ctd: structure containing the different oceanographic variables binned
%       on pressure intervals of 2 db and in the 0-1000 db range.
%   iso: structure containing the different oceanographic variables binned
%       on potential density intervals equal to the average potential
%       density at the 2 db pressure intervals of the values of ctd.p
% 
% EXAMPLE:
%   [ctd,iso] = hotFTPextract(210);
%
% VERSION HISTORY:
%   0.1     First version
%   0.11    Fix bugs in cast number and date; transform -9 to NaN,
%           fix bug in MLD computation
%           find DCM characteristics based on running median filter
%   0.12    Fix bug in Brunt Vaisala freuency calculation (complex numbers = NaN)
%   0.13    transform -99 to NaN
%   0.14    Fix bug that prevented to read CTD files for HOT cruises before 10
%   0.15    Change repository to ftp.soest.hawaii.edu
%
% Benedetto Barone - June 2020 - Version 0.15

% Open connection with ftp server
%nftp = ftp('mananui.soest.hawaii.edu');
nftp = ftp('ftp.soest.hawaii.edu');

hn = num2str(hot_n); % string of hot number
hn_ctd = hn; % string for hot ctd files
if hot_n<10
    hn_ctd = ['0' hn_ctd];
end
% Extract time of the casts (start time)
cd(nftp,'/hot/cruise.summaries/');
mget(nftp,['hot' hn '.sum']);
fid = fopen(['hot' hn '.sum']);
A = textscan(fid,'%*s %*s %f %f %*s %2f%2f%2f %2f%2f %*s %*f %*f %*s %*f %*f %*s %*s %*f %*f %*f %*s %*s %*s %*s %*s %*s','Headerlines',4,'collectoutput',1);
fclose(fid);
A = A{1};
A = A(A(:,1)==2,:);
i2000 = A(:,5)<85;
A(i2000,5) = 2000+A(i2000,5);% change year field
A(~i2000,5) = 1900+A(~i2000,5);% change year field
castnumber = A(:,2);
castdates = datenum(A(:,5),A(:,3),A(:,4),A(:,6),A(:,7),0*A(:,7));
castdates = castdates -10/24; % In Hawaii time
clear A i2000
delete(['hot' hn '.sum'])
% Initialize final variables
cd(nftp,['/hot/ctd/hot-' hn]);
files = dir(nftp,['h' hn_ctd 'a02*.ctd']);
lf = length(files);
ctd.cruise = hot_n;
ctd.p = (0:2:1000)';
ctd.depth = -gsw_z_from_p(ctd.p,22.75);
ctd.cast = NaN(1,lf);
ctd.date = NaN(1,lf);
ctd.t = NaN(501,lf);
ctd.sp = NaN(501,lf);
ctd.sa = NaN(501,lf);
ctd.ct = NaN(501,lf);
ctd.sig = NaN(501,lf);
ctd.f = NaN(501,lf);
ctd.o = NaN(501,lf);
ctd.n = NaN(501,lf);
% Extract data from single CTD data files
for i = 1:lf
    filename = files(i).name;
    mget(nftp,filename); % get file
    ctd.cast(i) = str2double(filename((end-5):(end-4)));
    ctd.date(i) = castdates(castnumber == ctd.cast(i));
    castdata = dlmread(filename,'',6,0); % load data
    castdata(castdata(:,1)>1000 | castdata(:,1)<0,:) = [];
    castdata(castdata==-9) = NaN; % assign NaN
    castdata(castdata==-99) = NaN; % assign NaN
    ind = castdata(:,1)/2 +1;
    ctd.t(ind,i) = castdata(:,2);
    ctd.sp(ind,i) = castdata(:,3);
    ctd.o(ind,i) = castdata(:,4);
    ctd.f(ind,i) = castdata(:,6);
    ctd.n(ind,i) = castdata(:,5);
    delete(filename) % delete file
    clear  castdata hdr filename
end
clear castdates files fid i ind
% Sort data 
[ctd.cast,ind] = sort(ctd.cast);
ctd.date = ctd.date(ind); ctd.t = ctd.t(:,ind); ctd.sp = ctd.sp(:,ind); ctd.f = ctd.f(:,ind); ctd.o = ctd.o(:,ind); ctd.n = ctd.n(:,ind);
clear ind
% Compute GSW-10 oceanographic variables
ctd.sa = gsw_SA_from_SP(ctd.sp,ctd.p,-158,22.75); % Absolute salinity
ctd.ct = gsw_CT_from_t(ctd.sa,ctd.t,ctd.p); % Conservative temperature
ctd.sig = gsw_sigma0(ctd.sa,ctd.ct); % Potential density anomaly (0 dbar)
% From mol/kg to mol/L
ctd.o = ctd.o.*(1000+ctd.sig)/1000; % from umol kg-1 to umol L-1
ctd.n = ctd.n.*(1000+ctd.sig)/1000; % from umol kg-1 to umol L-1

% Brunt Vaisala frequency
d_dz = NaN*ctd.sig;
ctd.bvf = NaN*ctd.sig;
for i = 1:lf 
    for j=2:500        
        d_dz(j,i) = (ctd.sig(j+1,i) - ctd.sig(j-1,i))/(ctd.depth(j+1) - ctd.depth(j-1));
    end
    ddz_for_bvf = d_dz(:,i);
    ddz_for_bvf(ddz_for_bvf<0) = NaN;
    ctd.bvf(:,i) = sqrt((9.8./(1000+ctd.sig(:,i))).*ddz_for_bvf);
    clear ddz_for_bvf
end
clear d_dz i j

% Find mixed layer depth (db)
ind_10 = find(ctd.p ==10);
mld = [ctd.sig - repmat(ctd.sig(ind_10,:),[501,1])];
ctd.mld003 = NaN(1,lf);
for i = 1:lf
    if sum(isnan(mld(:,i))==0)>0
        ind_003 = min(find(mld(:,i)>0.03 & ctd.p > 10));
        if isempty(ind_003)
            ctd.mld003(i) = NaN;
        else
            ctd.mld003(i) = ctd.p(ind_003);
        end
    end
end
clear mld ind_10 ind_003
d_1 = ctd.sig(1,:);
d_1(isnan(d_1)) = ctd.sig(2,isnan(d_1));
d_1(isnan(d_1)) = ctd.sig(3,isnan(d_1));
mld = [ctd.sig - repmat(d_1,[501,1])];
ctd.mld0125 = NaN(1,lf);
ctd.mld001 = NaN(1,lf);
ctd.mld0005 = NaN(1,lf);
for i = 1:lf
    if sum(isnan(mld(:,i))==0)>0
        ind_0125 = min(find(mld(:,i)>0.125));
        if isempty(ind_0125)
            ctd.mld0125(i) = NaN;
        else
            ctd.mld0125(i) = ctd.p(ind_0125);
        end
        ind_001 = min(find(mld(:,i)>0.01));
        if isempty(ind_001)
            ctd.mld001(i) = NaN;
        else
            ctd.mld001(i) = ctd.p(ind_001);
        end
        ind_0005 = min(find(mld(:,i)>0.005));        
        if isempty(ind_0005)
            ctd.mld0005(i) = NaN;
        else
            ctd.mld0005(i) = ctd.p(ind_0005);
        end
    end
end
clear ind_0005 ind_001 ind_0125 d_1 i mld

% Find charachteristics of the DCM based on 5 points Running Median filter
% below 30 db
ctd.pcm = NaN(1,lf);
ctd.fcm = NaN(1,lf);
ctd.sigcm = NaN(1,lf); 
for i = 1:lf
    [ctd.fcm(i),ind_max] = nanmax(RunMedian(ctd.f(16:126,i),5));
    ind_max = ind_max+15;
    if isnan(ctd.fcm(i))
        ctd.pcm(i) = NaN;
        ctd.sigcm(i) = NaN;
    else
        ctd.pcm(i) = ctd.p(ind_max);
        ctd.sigcm(i) = ctd.sig(ind_max,i);
    end
    clear ind_max
end

% Cast decimal hour
ctd.decimal_hour = datevec(ctd.date);
ctd.decimal_hour = ctd.decimal_hour(:,4) + ctd.decimal_hour(:,3)/60;

% Isopycnal data structure
iso.cruise = hot_n;
% Data on isopycnals
iso.sig = nanmean(ctd.sig,2);
iso.sig = sort(iso.sig);
% initialize table of isopycnal data
iso.t = NaN*ctd.t;
iso.sp = NaN*ctd.sa;
iso.sa = NaN*ctd.sa;
iso.ct = NaN*ctd.ct;
iso.f = NaN*ctd.f;
iso.o = NaN*ctd.o;
iso.n = NaN*ctd.n;
% interpolate at isopycnal values
for i = 1:lf
    ind_nan = isnan(ctd.sig(:,i));
    if sum(~ind_nan) ~= 0
        iso.t(:,i) = interp1(ctd.sig(~ind_nan,i),ctd.t(~ind_nan,i),iso.sig);
        iso.sp(:,i) = interp1(ctd.sig(~ind_nan,i),ctd.sp(~ind_nan,i),iso.sig);
        iso.sa(:,i) = interp1(ctd.sig(~ind_nan,i),ctd.sa(~ind_nan,i),iso.sig);
        iso.o(:,i) = interp1(ctd.sig(~ind_nan,i),ctd.o(~ind_nan,i),iso.sig);
        iso.f(:,i) = interp1(ctd.sig(~ind_nan,i),ctd.f(~ind_nan,i),iso.sig);
        iso.n(:,i) = interp1(ctd.sig(~ind_nan,i),ctd.n(~ind_nan,i),iso.sig);
        iso.ct(:,i) = interp1(ctd.sig(~ind_nan,i),ctd.ct(~ind_nan,i),iso.sig);
    end
    clear ind_nan
end
clear i
% Units
ctd.units = {'cruise: #';'p: db';'depth: m';'cast: #';'date: HST matlab format';'t: deg C';'sp: PSU';'sa: g kg-1';'ct: deg C'; ...
    'sig: kg m-3';'f: mg m-3';'o: mmol m-3';'n: mmol m-3';'bvf: s-1';'mld: db';'pcm: db';'fcm: mg m-3';'sigcm: kg m-3';'decimal_hour: h';};
iso.units = {'cruise: #';'sig: kg m-3';'t: deg C';'sp: PSU';'sa: g kg-1';'ct: deg C'; ...
    'f: mg m-3';'o: mmol m-3';'n: mmol m-3'};

%{
eval(['CTD_h' hn ' = ctd']);
eval(['ISO_h' hn ' = iso']);
eval(['save CTD_h' hn ' CTD_h' hn ' ISO_h' hn])
%}
    
%{
% Contours of oxygen

subplot(2,1,1)
contourf(ctd.date,ctd.p(1:101),ctd.o(1:101,:),100:1:240,'Edgecolor','none')
hold on, plot(ctd.date,ctd.mld0005,'k-'),hold off
set(gca,'ydir','rev')
datetick('x','HH:MM'); xlim([ctd.date(1) ctd.date(end)])
xlabel('HH:MM (HST)'); ylabel('Pressure (db)')
cb = colorbar; title(cb,'Oxygen mmol m^{-3}');
subplot(2,1,2)
contourf(ctd.date,ctd.p(1:101),iso.o(1:101,:),100:1:240,'Edgecolor','none')
set(gca,'ydir','rev')
datetick('x','HH:MM'); xlim([ctd.date(1) ctd.date(end)])
xlabel('HH:MM (HST)'); ylabel('Average isopycnal pressure (db)')
cb = colorbar; title(cb,'Oxygen mmol m^{-3}');
%}