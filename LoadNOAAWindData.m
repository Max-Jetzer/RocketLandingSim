function data = LoadNOAAWindData(cityCode)

%NM is 'elp'
URLL = 'https://www.aviationweather.gov/windtemp/data?level=low&fcst=06&region=all&layout=on&date=';
URLH = 'https://www.aviationweather.gov/windtemp/data?level=high&fcst=06&region=all&layout=on&date=';
s = urlread(URLL); %lower level winds for the Dallas-Fort Worth area

k1 = findstr(s,cityCode);
if( isempty(k1) )
  error('City code not found.');
end

kn = regexp(s, '[\n]');

kft = findstr(s,'FT');
k3 = kn(kn>kft);
k3 = k3(1);
header = s(kft:k3-1);

k2 = kn(kn>k1);
k2 = k2(1);

line = s(k1:k2-1);

hdrcell = strsplit(header);
linecell = strsplit(line);
dif = length(hdrcell)-length(linecell);
if dif == 1                                             %Making it work with initial altitudes greater than 3 kft
    new={'0000'};
    linecell = [linecell{1}, new, linecell{2:end}];
end
if dif == 2
    new={'0000', '0000'};
    linecell = [linecell{1}, new, linecell{2:end}];
end 
if dif == 3
    new={'0000', '0000', '0000'};
    linecell = [linecell{1}, new, linecell{2:end}];
end 

    
s = urlread(URLH); %Upper level winds too

k1 = findstr(s,cityCode);
if( isempty(k1) )
  fprintf('\nupper level winds unavailable');
else

kn = regexp(s, '[\n]');

kft = findstr(s,'FT');
k3 = kn(kn>kft);
k3 = k3(1);
header = s(kft:k3-1);

k2 = kn(kn>k1);
k2 = k2(1);

line = s(k1:k2-1);

hdrcell1 = strsplit(header);
linecell1 = strsplit(line);

hdrcell = [hdrcell, hdrcell1{2:end}];
linecell = [linecell, linecell1{2:end}];
end

altFt = [];
for i=2:length(hdrcell)
  if( ~isempty(hdrcell{i}) )
    altFt(end+1) = str2num(hdrcell{i});
  end
end

altM = floor(altFt*0.3048); %Truncating the meters to make it work better

windNorth = [];
windEast = [];
dir = [];
spd = [];
for i=2:length(linecell)
  if( ~isempty(linecell{i}) )
    windDir = str2num(linecell{i}(1:2));    %Wind direction is the direction that wind comes FROM.
    windSpd = -str2num(linecell{i}(3:4));   %So to fix, we just make windspeed negative.
    if windDir>36                          %Wind direction is encoded in tens of degrees.
        windDir=(windDir-50);                %When wind speed is greater than 100 knots,
        windSpd=-100-windSpd;               %they add 50 to the direction string and subtract 100 from speed.
    end                                     %That way it stays only 4 digits.
    if windDir==49
        windSpd=0.1;
    end
    windDir=10*windDir;                      
    dir(end+1)=windDir;
    spd(end+1)=windSpd;
    windNorth(end+1) = windSpd*cosd(windDir);
    windEast(end+1)  = windSpd*sind(windDir);
  end
end

knotsToMPS = 0.51444445;

windNorthMPS = windNorth*knotsToMPS;
windEastMPS = windEast*knotsToMPS;

data.altM  = [0 altM];
data.altFt = [0 altFt];
data.windNorthMPS = [0 windNorthMPS];
data.windEastMPS  = [0 windEastMPS];
data.windNorthFPS = [0 windNorthMPS/.3048];
data.windEastFPS  = [0 windEastMPS/.3048];
data.angle = dir;
data.windspd = spd*knotsToMPS;