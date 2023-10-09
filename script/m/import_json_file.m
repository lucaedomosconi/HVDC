fname = 'I_displ_file.json';
fid = fopen(fname);
raw = fread(fid,inf);
str = char(raw');
fclose(fid);
DATA = jsondecode(str);
