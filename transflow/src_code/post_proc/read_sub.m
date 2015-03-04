%==========================================================================
% Read a subsection from a skiped 3D subdomain file
%==========================================================================
% Usage:
% subregion = read_sub(filename,n1,n2,n3,p1,p2,p3,skip1,skip2,skip3,nvar);
% 
% n1, n1 and n3 are the sizes of the original domain for dims 1, 2 and 3
% p1, p1 and p3 are the region subscripts for dims 1, 2 and 3
% use [] to specify the entire domain range
%
% skip1, skip2 and skip3 are the skips for the subdomain along each dim
% 
% Example:
% 
% files = dir('SubDom_dble.*');
% 
% ngx = 3072 +1;
% ngy =  192 +1;
% ngz =  192 +1;
% 
% px = 1:1500;
% py = 20;
% pz = 1:97;
% 
% nvar = 3;
% V = read_sub(files(10).name,ngz,ngx,ngy,pz,px,py,kskip,iskip,jskip,nvar);
% 
% figure
% subplot(3,1,1)
% imagesc(V{1}); axis xy
% subplot(3,1,2)
% imagesc(V{2}); axis xy
% subplot(3,1,3)
% imagesc(V{3}); axis xy


function subregion = read_sub(filename,n1,n2,n3,p1,p2,p3,skip1,skip2,skip3,nvar)

% Size of data file
ns1 = ((n1-1)/skip1) + 1;
ns2 = ((n2-1)/skip2) + 1;
ns3 = ((n3-1)/skip3) + 1;

if isempty(p1)
    p1 = 1:ns1;
end
   
if isempty(p2)
    p2 = 1:ns2;
end

if isempty(p3)
    p3 = 1:ns3;
end

% Input subscripts
[D1 D2 D3] = ndgrid(p1,p2,p3);

D1 = D1(:);
D2 = D2(:);
D3 = D3(:);

% Find indices of input subscripts
s = sub2ind([ns1,ns2,ns3],D1,D2,D3);

% Find sequential points
s  = sort(s);
ss = find(diff(s)~=1);

assignin('base','s',s)
% if all points are sequential, read once else minimise read operations
if isempty(ss)
    slen = length(s);
    ss = s(1);
else
    slen = mode(diff(ss));
    ss = [0; ss];
    ss = ss + 1;
    ss = s(ss);
end

% ---- Read velocity field ----

% Open file
fid = fopen(filename,'r','ieee-le.l64');

% add fseek here if you have header data to skip

    for nv = 1:nvar
        
        for n = 1:length(ss)
            
            fseek(fid,(nv-1)*8*ns1*ns2*ns3,'bof');
            fseek(fid,8*(ss(n)-1),'cof');
            subregion{nv}(:,n) = fread(fid,slen,'double');
            
        end
    end

fclose(fid);

    % Clean up output
    for nv = 1:nvar
        subregion{nv} = reshape(subregion{nv},[length(p1),length(p2),length(p3)]);
        subregion{nv} = squeeze(subregion{nv});
    end
    
% End of function
    