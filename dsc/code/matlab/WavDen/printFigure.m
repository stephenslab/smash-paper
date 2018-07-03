function printFigure(k,file)

% printFigure
%
% Saves current figure window contents as a colour EPS ( -depsc ) file or etc.  
% This is of general use.
%
% For postscript, the -loose option was supposed to make it look like it does in the
% Matlab figure...  That didn't work. Setting paperpositionmode to auto does.
%
% cPbL, Oct 99: added optional arguments k and file
% cPbL Oct 99: turned off inverthardcopy in the figure
% cPbL, March1998: Created fcn

set(gcf,'paperpositionmode','auto'); % Should preseve wysiwig geometry.
set(gcf,'inverthardcopy','off'); % Should leave black as black and white as white in EPS

qList = {'EPS (colour)', 'EPS (b&w)',  'PS (colour)',  'PICT',   'JPEG',   'Printer (colour)', 'Printer (b&w)'};
devList={'-depsc -loose -adobecset -cmyk ',      '-deps -loose -adobecset ',  '-dpsc2 -loose -adobecset -cmyk ',     '-dpict',   '-djpeg',  '-dpsc2 -cmyk ', '-dps2'};
fileList={'file',        'file',        'file',        'file',   'file',    'printer',          'printer'};
sufList ={'*.eps',       '*.eps',       '*.ps',        '*.pict',  '*.jpeg'};

if nargin < 1,
   k=menu('Send to what device?',qList);
end%if

if strcmp(fileList{k},'file'),
   if nargin < 2
      [fullFileName,path] = uiputfile(sufList{k}, ['Choose a name for your ' qList{k} ' file...']);
      if (fullFileName == 0)
         return; %Cancel button selected, load canceled.
      end%if file cancelled
      file =['''' path fullFileName ''' '];    
   end%if interactive mode
end%if file

eval(['print ' devList{k} ' ' file ' ']);
disp([' Figure printed to ' devList{k} ' on ' fileList{k} ' ' file]); 