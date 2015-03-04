function [status] = change_line(codefile,variable_name,variables,variable_type)

%Store Present Working directory
pwdir = pwd;

%Go to code file
cd (codefile)

%Open MD.in file
fid = fopen('./MD.in');
%Find location in file
tline = fgets(fid);
var_in = []; lineno = 1;
while isempty(var_in)
    status = feof(fid); % Test for end of file
    if (status > 0)
        break
    end
    tline = fgets(fid);
    disp(tline)
    disp(variable_name)
    var_in = regexp(tline,variable_name);
    lineno = lineno + 1;
end

%If character not found then return error message and status 1
if (status ~= 0)
    disp(strcat('Character',variable_name,' not found'));
    return
end

%Backup previous MD.in
sedstr = 'cp MD.in MD.in.backup &&';
[status] = system(sedstr);
%Change line
for i = 1:length(variables)
    lineno = lineno + 1;
    if ( strcmp(variable_type,'integer'))
        sedstr = strcat('sed',32,int2str(lineno),'s/.*/',int2str(variables(i)),'/ <MD.in >new', ...
                    32,'&& mv new MD.in');
    else
        sedstr = strcat('sed',32,int2str(lineno),'s/.*/',num2str(variables(i)),'/ <MD.in >new', ...
                    32,'&& mv new MD.in');
    end
    [status] = system(sedstr);
end

cd(pwdir)

