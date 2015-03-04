%Establish size of domain walls
if (max(fixdistbot) ~= 0 || max(fixdisttop) ~= 0)
    nbinswallbot = round(fixdistbot./binsize);
    nbinswalltop = round(fixdisttop./binsize);
    wallbot = fixdistbot;
    walltop = fixdisttop;
elseif (max(tethdistbot) ~= 0 || max(tethdisttop) ~= 0)
    nbinswallbot = round(tethdistbot./binsize);
    nbinswalltop = round(tethdisttop./binsize);
    wallbot = tethdistbot;
    walltop = tethdisttop;
end

%Establish size of liquid region
nbinsliquid = nbins-nbinswalltop - nbinswallbot;
liquiddomain = globaldomain-walltop - wallbot;
    

%Store string for sliding wall
if ((max(slidedisttop) && max(slidedistbot)) ~= 0)
	sliding_wall = 'both';
	error('Output from MD codes does not allow this yet') 
elseif (max(slidedisttop ) ~= 0)
    sliding_wall = 'top';
    topwallslidev = wallslidev;
    botwallslidev = 0.0;
elseif (max(slidedistbot) ~= 0)
    sliding_wall = 'bottom';
    topwallslidev = 0.0;
    botwallslidev = wallslidev;
end




