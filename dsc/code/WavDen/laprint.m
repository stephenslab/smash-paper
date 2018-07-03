function laprint(figno,filename,varargin)

%  This function prints a figure for inclusion in LaTeX documents.
%  It creates an eps-file and a tex-file. The tex-file contains the
%  annotation of the figure such as titles, labels and texts. The
%  eps-file contains the non-text part of the figure as well as the 
%  position of the text-objects. The packages 'epsfig' and 'psfrag' are 
%  required for the LaTeX run. A postscript driver like 'dvips' is 
%  required for printing.
%
%
% Usage:
%
%      laprint or laprint(figno) or laprint filename or laprint(figno,filename)
%      or laprint(figno,filename,opt1,opt2,..).
%
%
% Input:
%
%    figno            : integer, figure to be printed (default: current figure)
%    filename         : string, basename of files to be created (default:
%                       unnamed). Thus filename.tex and filename.eps are created.
%    opt1,..          : strings, describing optional inputs as follows
%    'width=xx'       : xx is the width of the figure in the tex file
%                       (in cm). Its default is 12. 
%    'factor=xx'      : xx is the factor by which the figure in the LaTeX
%                       document is smaller than the figure in the
%                       postscript file. Thus, if 10pt fonts are active
%                       while LaTeX scans the tex-file, then 10pt fonts 
%                       are inserted in the postscript file (by psfrag).
%                       The figure in the postscript file is scaled by xx,
%                       such that (10*xx)pt fonts are finally used in the 
%                       LaTeX figure. Linethicknesses etc. scale 
%                       accordingly. The default is factor=0.8. A
%                       non-positive number for xx lets the factor be
%                       computed such that the figure  in the postscript
%                       file has the same size as the figure on screen.
%    'asonscreen'     : prints a graphics 'as on screen', retaining
%                       ticks, ticklabels and lims.
%    'verbose'        : verbose mode; asks before overwriting files and
%                       issues some more messages. 
%    'keepticklabels' : keeps the tick labels within the eps-file
%    'mathticklabels' : tick labels are set in LaTeX math mode
%    'keepfontprops'  : tries to translate the MATLAB font properties 
%                       (size, width, angle) into similar LaTeX font
%                       properties. (By default, laprint does not 
%                       introduce any LaTeX font selection commands.)
%    'noscalefonts'   : does not scale the fonts with the figure. 
%    'noextrapicture' : does not add extra picture environments. The 
%                       picture would be empty, but alows to place LaTeX 
%                       objects in arbitrary positions at later times.
%    'nofigcopy'      : directly modifies the figure figno. It will
%                       get messed up by tags.
%    'nohead'         : does not place a commenting head in the tex-file. 
%    'caption=xx'     : adds \caption{xx} and \label{fig:filename} 
%                       entries to the tex-file.
%    'comment=xx'     : places the comment xx into the header of the
%                       tex-file
%    'viewfile=xx'    : creates an additional file xx.tex containing a
%                       LaTeX document which calls the tex-file.
%
%  Examples:
%     Suppose you have created Figure No. 1. Calling 
%     >> laprint(1,'f1')  (or   >> laprint f1, if Figure No. 1 is current)    
%     creates f1.eps and f1.tex. You can input f1.tex into a LaTeX-document:
%     \documentclass{article}\usepackage{epsfig,psfrag}..\input{f1}..
%     This will create a figure of width 12cm. If 10pt fonts are active
%     while scanning f1.tex, then the text in the figure will be of size 8pt.  %
%     Similarly
%     >> laprint(1,'f2','width=8','factor=0.6')
%     will create a 8cm x 6cm graphics with 6pt text.

%  This version of laprint is written and tested under MATLAB 5.3.
%  More information on the usage and the most recent version may be 
%  found under http://www.uni-kassel.de/~linne/.
%  Report bugs, suggestions and comments to linnemann@uni-kassel.de.

%  known problems and limitations: 
%  --  The matlab functions copyobj and plotedit have bugs. 
%      If this is a problem, use option 'nofigcopy'.
%  --  multi-line text is not supported

%  (c) 1999   Arno Linnemann.   All rights reserved. 
%  The author of this program assumes no responsibility for any  errors 
%  or omissions. In no event shall he be liable for damages  arising out of 
%  any use of the software. Redistribution of the unchanged file is allowed.
%  Distribution of changed versions is allowed provided the file is renamed
%  and the source and authorship of the original version is acknowledged in 
%  the modified file.

%  Written by:
%  Arno Linnemann,
%  Control and Systems Theory,
%  Department of Electrical Engineering, 
%  University of Kassel,
%  34109 Kassel,
%  Germany.
%  mailto: linnemann@uni-kassel.de
%  http://www.uni-kassel.de/~linne/

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%% PART 1: Initialize, check inputs, etc.
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize
laprintident = '1.01 (1.10.99)'; 
furtheroptions='';

% default values
verbose=0;
asonscreen=0;
keepticklabels=0;
mathticklabels=0;
keepfontprops=0;
extrapicture=1;
nofigcopy=0;
nohead=0;
noscalefonts=0;
caption=0;
commenttext='';
captiontext='';
width=12;
factor=0.8;
viewfile=0;

% no output
if nargout
  error('No output argument, please.')
end  

% check inputs 1 and 2
if nargin==0
  figno=gcf;
  filename='unnamed';
elseif nargin==1
  if ~isa(figno,'char')
    filename='unnamed';
  else  
    filename=figno;
    figno=gcf;
  end  
end
if ~isa(figno,'double') 
  figno
  error('This is not a figure handle.') 
end
if ~any(get(0,'children')==figno)
  figno
  error('This is not a figure handle.') 
end
if ~isa(filename,'char')
  filename
  error('This is not a file name.') 
end

% read and check options  
if nargin>2  
  for i=1:nargin-2
    if ~isa(varargin{i},'char')
      error('Options must be character arrays.')
    end  
    oriopt=varargin{i}(:)';
    opt=[ lower(strrep(oriopt,' ','')) '                   ' ];
    if strcmp(opt(1:7),'verbose')
      verbose=1;
      furtheroptions=[ furtheroptions ' / ' deblank(opt) ]; 
    elseif strcmp(opt(1:10),'asonscreen')
      asonscreen=1;
      furtheroptions=[furtheroptions ' / ' deblank(opt) ]; 
    elseif strcmp(opt(1:14),'keepticklabels')
      keepticklabels=1;
      furtheroptions=[furtheroptions ' / ' deblank(opt) ]; 
    elseif strcmp(opt(1:14),'mathticklabels')
      mathticklabels=1;
      furtheroptions=[furtheroptions ' / ' deblank(opt) ]; 
    elseif strcmp(opt(1:13),'keepfontprops')
      keepfontprops=1;
      furtheroptions=[furtheroptions ' / ' deblank(opt) ]; 
    elseif strcmp(opt(1:14),'noextrapicture')
      extrapicture=0;
      furtheroptions=[furtheroptions ' / ' deblank(opt) ]; 
    elseif strcmp(opt(1:9),'nofigcopy')
      nofigcopy=1;
      furtheroptions=[furtheroptions ' / ' deblank(opt) ]; 
    elseif strcmp(opt(1:12),'noscalefonts')
      noscalefonts=1;
      furtheroptions=[furtheroptions ' / ' deblank(opt) ]; 
    elseif strcmp(opt(1:6),'nohead')
      nohead=1;
      furtheroptions=[furtheroptions ' / ' deblank(opt) ]; 
    elseif strcmp(opt(1:7),'caption')
      caption=1;
      eqpos=findstr(oriopt,'=');
      if isempty(eqpos)
	furtheroptions=[furtheroptions ' / ' deblank(opt) ];
	captiontext=[];
      else	
	furtheroptions=[furtheroptions ' / ' oriopt ];
	captiontext=oriopt(eqpos+1:length(oriopt));
      end	
    elseif strcmp(opt(1:8),'comment=')
      eqpos=findstr(oriopt,'=');
      furtheroptions=[furtheroptions ' / ' oriopt ];
      commenttext=oriopt(eqpos(1)+1:length(oriopt));
    elseif strcmp(opt(1:9),'viewfile=')
      viewfile=1;
      eqpos=findstr(oriopt,'=');
      furtheroptions=[furtheroptions ' / ' oriopt ];
      viewfilename=oriopt(eqpos(1)+1:length(oriopt));
    elseif strcmp(opt(1:6),'width=')
      eval([ opt ';' ]);
    elseif strcmp(opt(1:7),'factor=')
      eval([ opt ';' ]);
    else
      error([ 'option ' varargin{i} ' not recognized.'])
    end   
  end
end
furtheroptions=strrep(strrep(furtheroptions,'\','\\'),'%','%%');
captiontext=strrep(strrep(captiontext,'\','\\'),'%','%%');
commenttext=strrep(strrep(commenttext,'\','\\'),'%','%%');

if verbose, 
  disp([ 'This is laprint, version ' laprintident '.' ]); 
end  

if mathticklabels
  Do='$';
else  
  Do='';
end  

% eps- and tex- filenames
[epsfullnameext,epsbasenameext,epsbasename,epsdirname]= ...
   getfilenames(filename,'eps',verbose);

[texfullnameext,texbasenameext,texbasename,texdirname]= ...
                 getfilenames(filename,'tex',verbose);
if ~strcmp(texdirname,epsdirname)
  disp('warning: eps-file and tex-file are placed in different directories.')
end  
if viewfile
  [viewfullnameext,viewbasenameext,viewbasename,viewdirname]= ...
                 getfilenames(viewfilename,'tex',verbose);
  if strcmp(texfullnameext,viewfullnameext)
    error('The tex- and view-file coincide. Use different names.')
  end  
  if ~strcmp(texdirname,viewdirname)
    disp([ 'warning: eps-file and view-file are placed '...
	   'in different directories.' ])
  end  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%% PART 2: Create new figure, insert tags, and bookkeep original text
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% open new figure (if required) and set properties

if ~nofigcopy
  figno=copyobj(figno,0);
  set(figno,'Numbertitle','off')
  set(figno,'MenuBar','none')
  %plotedit(figno,'hidetoolsmenue')
  pause(0.5)  
end  

if asonscreen
  xlimmodeauto=findobj(figno,'xlimmode','auto');
  xtickmodeauto=findobj(figno,'xtickmode','auto');
  xticklabelmodeauto=findobj(figno,'xticklabelmode','auto');
  ylimmodeauto=findobj(figno,'ylimmode','auto');
  ytickmodeauto=findobj(figno,'ytickmode','auto');
  yticklabelmodeauto=findobj(figno,'yticklabelmode','auto');
  zlimmodeauto=findobj(figno,'zlimmode','auto');
  ztickmodeauto=findobj(figno,'ztickmode','auto');
  zticklabelmodeauto=findobj(figno,'zticklabelmode','auto');
  set(xlimmodeauto,'xlimmode','manual')
  set(xtickmodeauto,'xtickmode','manual')
  set(xticklabelmodeauto,'xticklabelmode','manual')
  set(ylimmodeauto,'ylimmode','manual')
  set(ytickmodeauto,'ytickmode','manual')
  set(yticklabelmodeauto,'yticklabelmode','manual')
  set(zlimmodeauto,'ylimmode','manual')
  set(ztickmodeauto,'ytickmode','manual')
  set(zticklabelmodeauto,'yticklabelmode','manual')
end  
set(figno,'paperunits','centimeters');
set(figno,'units','centimeters');
oripp=get(figno,'PaperPosition');
orip=get(figno,'Position');

if factor <= 0
  factor=width/orip(3);
end 
latexwidth=width;
epswidth=latexwidth/factor;
epsheight = epswidth*orip(4)/orip(3);

set(figno,'PaperPosition',[1 1 epswidth epsheight ])
set(figno,'Position',[orip(1)+0.5 orip(2)-0.5 epswidth epsheight ])
set(figno,'Name',[ 'To be printed; size: ' num2str(factor,3) ...
      ' x (' num2str(epswidth,3) 'cm x ' num2str(epsheight,3) 'cm)' ])

asonscreen_dummy=0;
if asonscreen_dummy
  set(xlimmodeauto,'xlimmode','auto')
  set(xtickmodeauto,'xtickmode','auto')
  set(xticklabelmodeauto,'xticklabelmode','auto')
  set(ylimmodeauto,'ylimmode','auto')
  set(ytickmodeauto,'ytickmode','auto')
  set(yticklabelmodeauto,'yticklabelmode','auto')
  set(zlimmodeauto,'ylimmode','auto')
  set(ztickmodeauto,'ytickmode','auto')
  set(zticklabelmodeauto,'yticklabelmode','auto')
end
  
% some warnings
if (epswidth<13) | (epsheight<13*0.75)
  disp('warning: The size of the eps-figure is quite small.')
  disp('         The text objects might not be properly set.')
  disp('         Reducing ''factor'' might help.')
end
if latexwidth/epswidth<0.5
  disp([ 'warning: The size of the eps-figure is large compared ' ...
	  'to the latex figure.' ])
  disp('         The text size might be too small.')
  disp('         Increasing ''factor'' might help.')
end  
if (orip(3)-epswidth)/orip(3) > 0.1
  disp('warning: The size of the eps-figure is much smaller than the original')
  disp('         figure on screen. Matlab might save different ticks and')
  disp([ '         ticklabels than in the original figure. See option ' ...
            '''asonsceen''. ' ])
end  

if verbose
  disp('Strike any key to continue.');
  pause
end  

%
% TEXT OBJECTS: modify new figure 
%

% find all text objects
hxl=get(findobj(figno,'type','axes'),'xlabel');
hyl=get(findobj(figno,'type','axes'),'ylabel');
hzl=get(findobj(figno,'type','axes'),'zlabel');
hti=get(findobj(figno,'type','axes'),'title');
hte=findobj(figno,'type','text');
% array of all text handles
htext=[ celltoarray(hxl) celltoarray(hyl) celltoarray(hzl) ...
      celltoarray(hti) celltoarray(hte)];
nt=length(htext);

% generate new strings and store old ones
oldstr=get(htext,'string');
newstr=cell(nt,1);
basestr='str00';
for i=1:nt
  if isa(oldstr{i},'cell')
    if length(oldstr{i})>1
      disp('warning: Annotation in form of a cell is currently not supported.')
      disp('         Ignoring all but first component.')
    end
    % To do: place a parbox here. 
    oldstr{i}=oldstr{i}{1};
  end  
  if size(oldstr{i},1)>1
    disp([ 'warning: Annotation in form of string matrices ' ...
	  'is currently not supported.' ])
    disp('         Ignoring all but first row.')
    % To do: place a parbox here. 
    oldstr{i}=oldstr{i}(1,:);
  end  
  if length(oldstr{i})
    oldstr{i}=strrep(strrep(oldstr{i},'\','\\'),'%','%%');
    newstr{i} = overwritetail(basestr,i);
  else  
    newstr{i}='';    
  end
end

% replace strings in figure
for i=1:nt
    set(htext(i),'string',newstr{i});
    %set(htext(i),'visible','on');
end    

% get alignments
hora=get(htext,'HorizontalAlignment');
vera=get(htext,'VerticalAlignment');
align=cell(nt,1);
for i=1:nt
  align{i}=hora{i}(1);
  if strcmp(vera{i},'top')
    align{i}=[align{i} 't'];
  elseif strcmp(vera{i},'cap')
    align{i}=[align{i} 't'];
  elseif strcmp(vera{i},'middle')
    align{i}=[align{i} 'c'];
  elseif strcmp(vera{i},'baseline')
    align{i}=[align{i} 'B'];
  elseif strcmp(vera{i},'bottom')
    align{i}=[align{i} 'b'];
  end
end  

% get font properties and create commands
[fontsizecmd{1:nt}] = deal('');
[fontanglecmd{1:nt}] = deal('');
[fontweightcmd{1:nt}] = deal('');
selectfontcmd='';
  
if keepfontprops

  % fontsize
  set(htext,'fontunits','points');
  fontsize=get(htext,'fontsize');
  for i=1:nt
    fontsizecmd{i}=[ '\\fontsize{' num2str(fontsize{i}) '}{' ...
	  num2str(fontsize{i}*1.5) '}'  ];
  end
    
  % fontweight
  fontweight=get(htext,'fontweight');
  for i=1:nt
    if strcmp(fontweight{i},'light')
      fontweightcmd{i}=[ '\\fontseries{l}\\mathversion{normal}' ];
    elseif strcmp(fontweight{i},'normal')
      fontweightcmd{i}=[ '\\fontseries{m}\\mathversion{normal}' ];
    elseif strcmp(fontweight{i},'demi')
      fontweightcmd{i}=[ '\\fontseries{sb}\\mathversion{bold}' ];
    elseif strcmp(fontweight{i},'bold')
      fontweightcmd{i}=[ '\\fontseries{bx}\\mathversion{bold}' ];
    else
      disp([ ' warning: unknown fontweight:' fontweight{i} ])
      fontweightcmd{i}=[ '\\fontseries{m}\\mathversion{normal}' ];
    end
  end  

  % fontangle
  fontangle=get(htext,'fontangle');
  for i=1:nt
    if strcmp(fontangle{i},'normal')
      fontanglecmd{i}=[ '\\fontshape{n}' ];
    elseif strcmp(fontangle{i},'italic')
      fontanglecmd{i}=[ '\\fontshape{it}' ];
    elseif strcmp(fontangle{i},'oblique')
      fontangle{i}=[ '\\fontshape{it}' ];
    else
      disp([ ' warning: unknown fontangle:' fontangle{i} ])
      fontanglecmd{i}=[ '\\fontshape{n}' ];
    end
  end  
  selectfontcmd= '\\selectfont ';
   
end

%
% LABELS: modify new figure
%

if ~keepticklabels

  % all axes
  hax=celltoarray(findobj(figno,'type','axes'));
  na=length(hax);

  % try to figure out if we have 3D axes an warn
  issuewarning=0;
  for i=1:na
    issuewarning=max(issuewarning,is3d(hax(i)));
  end
  if issuewarning
    disp('warning: There seems to be a 3D plot. The LaTeX labels are ')
    disp('         possibly incorrect. The option  ''keepticklabels'' might')
    disp('         help. The option ''nofigcopy'' might be wise, too.')
  end

  % try to figure out if we linear scale with extra factor 
  % and determine powers of 10
  powers=NaN*zeros(na,3);  % matrix with powers of 10 
  for i=1:na                    % all axes
    allxyz={ 'x', 'y', 'z' };
    for ixyz=1:3                % x,y,z
      xyz=allxyz{ixyz};
      ticklabelmode=get(hax(i),[ xyz 'ticklabelmode']);
      if strcmp(ticklabelmode,'auto')
        tick=get(hax(i),[ xyz 'tick']);
        ticklabel=get(hax(i),[ xyz 'ticklabel']);	      
	nticks=size(ticklabel,1);
	if nticks==0,
          powers(i,ixyz)=0;
	end  
        for k=1:nticks        % all ticks
	  label=str2num(ticklabel(k,:));
	  if length(label)==0, 
	    powers(i,ixyz)=0;
	    break; 
	  end  
	  if ( label==0 ) & ( abs(tick(k))>1e-10 )
	    powers(i,ixyz)=0;
	    break; 
          end	      
	  if label~=0    
            expon=log10(tick(k)/label);
	    rexpon=round(expon);
	    if abs(rexpon-expon)>1e-10
              powers(i,ixyz)=0;
	      break; 
	    end	
            if isnan(powers(i,ixyz))
	      powers(i,ixyz)=rexpon;
	    else 	
	      if powers(i,ixyz)~=rexpon
        	powers(i,ixyz)=0;
	        break; 
              end		
	    end 
          end  	    
	end % k	    
      else % if 'auto'
        powers(i,ixyz)=0;
      end % if 'auto'
    end % ixyz
  end % i
  
  % replace all ticklabels and bookkeep
  nxlabel=zeros(1,na);
  nylabel=zeros(1,na);
  nzlabel=zeros(1,na);
  allxyz={ 'x', 'y', 'z' };
  for ixyz=1:3
    xyz=allxyz{ixyz};
    k=1;
    basestr=[ xyz '00' ];
    if strcmp(xyz,'y') % 'y' is not horizontally centered! 
      basestr='v00';
    end  
    oldtl=cell(na,1);
    newtl=cell(na,1);
    for i=1:na
      % set(hax(i),[ xyz 'tickmode' ],'manual')
      % set(hax(i),[ xyz 'ticklabelmode' ],'manual')
      oldtl{i}=chartocell(get(hax(i),[ xyz 'ticklabel' ]));
      nlabel(i)=length(oldtl{i});
      newtl{i}=cell(1,nlabel(i));
      for j=1:nlabel(i)
        newtl{i}{j} = overwritetail(basestr,k);
        k=k+1;
        oldtl{i}{j}=deblank(strrep(strrep(oldtl{i}{j},'\','\\'),'%','%%'));
      end
      set(hax(i),[ xyz 'ticklabel' ],newtl{i});
    end  
    eval([ 'old' xyz 'tl=oldtl;' ]);
    eval([ 'new' xyz 'tl=newtl;' ]);
    eval([ 'n' xyz 'label=nlabel;' ]);
  end

  % determine latex commands for font properties
  
  if keepfontprops

    % font size
    afsize=zeros(na,1);
    for i=1:na
      afsize(i)=get(hax(i),'fontsize');
    end          
    if (any(afsize ~= afsize(1) ))
      disp('warning: Different font sizes for axes not supported.')
      disp([ '         All axses will have font size ' ...
	     num2str(afsize(1)) '.' ] )
    end      
    afsizecmd = [ '\\fontsize{' num2str(afsize(1)) '}{' ...
	  num2str(afsize(1)*1.5) '}'  ];

    % font weight
    afweight=cell(na,1);
    for i=1:na
      afweight{i}=get(hax(i),'fontweight');
    end
    if strcmp(afweight{1},'light')
      afweightcmd=[ '\\fontseries{l}\\mathversion{normal}' ];
    elseif strcmp(afweight{1},'normal')
      afweightcmd=[ '\\fontseries{m}\\mathversion{normal}' ];
    elseif strcmp(afweight{1},'demi')
      afweightcmd=[ '\\fontseries{sb}\\mathversion{bold}' ];
    elseif strcmp(afweight{1},'bold')
      afweightcmd=[ '\\fontseries{bx}\\mathversion{bold}' ];
    else
      disp([ ' warning: unknown fontweight:' afweight{1} ])
      afweightcmd=[ '\\fontseries{m}\\mathversion{normal}' ];
    end
    for i=1:na
      if ~strcmp(afweight{i},afweight{1})
        disp('warning: Different font weights for axes not supported.')
        disp([ '         All axes will have font weight ' afweightcmd ] )
      end      
    end      

    % font angle
    afangle=cell(na,1);
    for i=1:na
      afangle{i}=get(hax(i),'fontangle');
    end
    if strcmp(afangle{1},'normal')
      afanglecmd=[ '\\fontshape{n}' ];
    elseif strcmp(afangle{1},'italic')
      afanglecmd=[ '\\fontshape{it}' ];
    elseif strcmp(afangle{1},'oblique')
      afanglecmd=[ '\\fontshape{it}' ];
    else
      disp([ ' warning: unknown fontangle:' afangle{1} ])
      afanglecmd=[ '\\fontshape{n}' ];
    end
    for i=1:na
      if ~strcmp(afangle{i},afangle{1})
        disp('warning: Different font angles for axes not supported.')
        disp([ '         All axes will have font angle ' afanglecmd ] )
      end      
    end      
  
  end

end

%
% extra picture environment
%

if extrapicture
  %%%%%%%%%%%%%%%%%%% AJOUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % all axes
  hax=celltoarray(findobj(figno,'type','axes'));
  na=length(hax);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  unitlength=zeros(na,1);
  ybound=zeros(na,1);
  for i=1:na
    if ~is3d(hax(i))
      xlim=get(hax(i),'xlim');
      ylim=get(hax(i),'ylim');
      axes(hax(i));
      hori=text(ylim(1),ylim(1),[ 'origin' int2str(i) ]);
      set(hori,'VerticalAlignment','bottom');
      set(hori,'Fontsize',2);
      pos=get(hax(i),'Position');
      unitlength(i)=pos(3)*epswidth;
      ybound(i)=(pos(4)*epsheight)/(pos(3)*epswidth);
    else
      disp('warning: Option ''extrapicture'' for 3D axes not supported.')
    end
  end 
end  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%% PART 3: save eps and tex files
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% save eps file
%
cmd=[ 'print(''-deps'',''-f' int2str(figno) ''',''' epsfullnameext ''')' ];
if verbose
  disp([ 'executing: '' ' cmd ' ''' ]);
end
eval(cmd);

%
% create latex file
%
if verbose
  disp([ 'writing to: '' ' texfullnameext ' ''' ])
end
fid=fopen(texfullnameext,'w');

% head
if ~nohead
  fprintf(fid,[ '%% This file is generated by the MATLAB m-file laprint.m.' ...
       ' It can be included\n']);
  fprintf(fid,[ '%% into LaTeX documents using the packages epsfig and ' ...
       'psfrag. It is accompanied\n' ]);
  fprintf(fid,  '%% by a postscript file. A sample LaTeX file is:\n');
  fprintf(fid, '%%    \\documentclass{article} \\usepackage{epsfig,psfrag}\n');
  fprintf(fid,[ '%%    \\begin{document}\\begin{figure}\\input{' ...
	texbasename '}\\end{figure}\\end{document}\n' ]);
  fprintf(fid, [ '%% See http://www.uni-kassel.de/~linne/ for recent ' ...
	'versions of laprint.m.\n' ]);
  fprintf(fid,  '%%\n');
  fprintf(fid,[ '%% created by:           ' 'laprint version ' ...
	laprintident '\n' ]);
  fprintf(fid,[ '%% created on:           ' datestr(now) '\n' ]);
  fprintf(fid,[ '%% options used:        ' furtheroptions '\n' ]);
  fprintf(fid,[ '%% latex width:          ' num2str(latexwidth) ' cm\n' ]);
  fprintf(fid,[ '%% factor:               ' num2str(factor) '\n' ]);
  fprintf(fid,[ '%% eps file name:        ' epsbasenameext '\n' ]);
  fprintf(fid,[ '%% eps bounding box:     ' num2str(epswidth) ...
      ' cm x ' num2str(epsheight) ' cm\n' ]);
  fprintf(fid,[ '%% comment:              ' commenttext '\n' ]);
  fprintf(fid,'%%\n');
else 
  fprintf(fid,[ '%% generated by laprint.m\n' ]);
  fprintf(fid,'%%\n');
end

% go on
fprintf(fid,'\\begin{psfrags}%%\n');
%fprintf(fid,'\\fontsize{10}{12}\\selectfont%%\n');
fprintf(fid,'\\psfragscanon%%\n');

% text strings

numbertext=0;
for i=1:nt
  numbertext=numbertext+length(newstr{i});
end
if numbertext>0,
  fprintf(fid,'%%\n');
  fprintf(fid,'%% text strings:\n');
  for i=1:nt
    if length(newstr{i})
      alig=strrep(align{i},'c','');
      fprintf(fid,[ '\\psfrag{' newstr{i} '}[' alig '][' alig ']{' ...
        fontsizecmd{i} fontweightcmd{i} fontanglecmd{i} selectfontcmd ...
        oldstr{i} '}%%\n' ]);
    end
  end
end

% labels

if ~keepticklabels
  if keepfontprops
    fprintf(fid,'%%\n');
    fprintf(fid,'%% axes font properties:\n');
    fprintf(fid,[ afsizecmd afweightcmd '%%\n' ]);
    fprintf(fid,[ afanglecmd '\\selectfont%%\n' ]);
  end  
  nxlabel=zeros(1,na);
  nylabel=zeros(1,na);
  nzlabel=zeros(1,na);
  for i=1:na
    nxlabel(i)=length(newxtl{i});
    nylabel(i)=length(newytl{i});
    nzlabel(i)=length(newztl{i});
  end    
      
  allxyz={ 'x', 'y', 'z' };
  for ixyz=1:3
    xyz=allxyz{ixyz};
    eval([ 'oldtl=old' xyz 'tl;' ]);
    eval([ 'newtl=new' xyz 'tl;' ]);
    eval([ 'nlabel=n' xyz 'label;' ]);
    if sum(nlabel) > 0
      fprintf(fid,'%%\n');
      fprintf(fid,[ '%% ' xyz 'ticklabels:\n']);
      if xyz=='x'
        poss='[t][t]';
      else
        poss='[r][r]';
      end  
      for i=1:na
        if nlabel(i)
          if strcmp(get(hax(i),[ xyz 'scale']),'linear')
	    % lin scale
	    % all but last
            for j=1:nlabel(i)-1
              fprintf(fid,[ '\\psfrag{' newtl{i}{j} '}' poss '{' ...
		  Do oldtl{i}{j} Do '}%%\n' ]);
            end 
            % last
            rexpon=powers(i,ixyz);
	    if rexpon
	      if xyz=='x'
	        fprintf(fid,[ '\\psfrag{' newtl{i}{nlabel(i)} ...
                   '}' poss '{\\shortstack{' ... 
                   Do oldtl{i}{nlabel(i)} Do '\\\\$\\times 10^{'...
		   int2str(rexpon) '}\\ $}}%%\n' ]);
	       else
	         fprintf(fid,[ '\\psfrag{' newtl{i}{nlabel(i)} ...
		     '}' poss '{' Do oldtl{i}{nlabel(i)} Do ...
                   '\\setlength{\\unitlength}{1ex}%%\n' ...
		   '\\begin{picture}(0,0)\\put(0.5,1.5){$\\times 10^{' ...
		   int2str(rexpon) '}$}\\end{picture}}%%\n' ]);
	       end
            else
	      fprintf(fid,[ '\\psfrag{' newtl{i}{nlabel(i)} '}' poss '{' ...
                Do oldtl{i}{nlabel(i)} Do '}%%\n' ]);
            end
          else
          % log scale
            for j=1:nlabel
              fprintf(fid,[ '\\psfrag{' newtl{i}{j} '}' poss '{$10^{' ...
                oldtl{i}{j} '}$}%%\n' ]);
            end
          end
        end   
      end
    end
  end
end  

% extra picture
if extrapicture
  %%%%%%%%%%%%%%%%%%% AJOUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % all axes
  hax=celltoarray(findobj(figno,'type','axes'));
  na=length(hax);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  fprintf(fid,'%%\n');
  fprintf(fid,'%% extra picture(s):\n');
  for i=1:na
    fprintf(fid,[ '\\psfrag{origin' int2str(i) '}[lb][lb]{' ...
                  '\\setlength{\\unitlength}{' ...
		  num2str(unitlength(i),'%5.5f') 'cm}%%\n' ]);
    fprintf(fid,[ '\\begin{picture}(1,' ...
		  num2str(ybound(i),'%5.5f') ')%%\n' ]);
    %fprintf(fid,'\\put(0,0){}%% lower left corner\n');
    %fprintf(fid,[ '\\put(1,' num2str(ybound(i),'%5.5f') ...
    %	          '){}%% upper right corner\n' ]);
    fprintf(fid,'\\end{picture}%%\n');
    fprintf(fid,'}%%\n');
  end
end  

% figure
fprintf(fid,'%%\n');
fprintf(fid,'%% Figure:\n');
if caption
  fprintf(fid,[ '\\parbox{' num2str(latexwidth) 'cm}{\\centering%%\n' ]);
end  
if noscalefonts
  fprintf(fid,[ '\\epsfig{file=' epsbasenameext ',width=' ...
	num2str(latexwidth) 'cm}}%%\n' ]);
else  
  fprintf(fid,[ '\\resizebox{' num2str(latexwidth) 'cm}{!}' ...
      '{\\epsfig{file=' epsbasenameext '}}%%\n' ]);
end
if caption
  if isempty(captiontext)
    captiontext=[ texbasenameext ', ' epsbasenameext ];
  end  
  fprintf(fid,[ '\\caption{' captiontext '}%%\n' ]);
  fprintf(fid,[ '\\label{fig:' texbasename '}%%\n' ]);
  fprintf(fid,[ '}%%\n' ]);
end  
fprintf(fid,'\\end{psfrags}%%\n');
fprintf(fid,'%%\n');
fprintf(fid,[ '%% End ' texbasenameext '\n' ]);
fclose(fid);

set(figno,'Name','Printed by laprint')
if ~nofigcopy
  if verbose
    disp('Strike any key to continue.');
    pause
  end  
  close(figno)
end

%
% create view file
%
if viewfile
  if verbose
    disp([ 'writing to: '' ' viewfullnameext ' ''' ])
  end
  fid=fopen(viewfullnameext,'w');

  if ~nohead
    fprintf(fid,[ '%% This file is generated by laprint.m.\n' ]);
    fprintf(fid,[ '%% It calls ' texbasenameext ...
		  ', which in turn  calls ' epsbasenameext '.\n' ]);
    fprintf(fid,[ '%% Process this file using\n' ]);
    fprintf(fid,[ '%%   latex ' viewbasenameext '\n' ]);
    fprintf(fid,[ '%%   dvips -o' viewbasename '.ps ' viewbasename '.dvi' ...
		'\n']);
    fprintf(fid,[ '%%   ghostview ' viewbasename '.ps&\n' ]);
  else 
    fprintf(fid,[ '%% generated by laprint.m\n' ]);
  end

  fprintf(fid,[ '\\documentclass{article}\n' ]);
  fprintf(fid,[ '\\usepackage{epsfig,psfrag,a4}\n' ]);
  fprintf(fid,[ '\\usepackage[latin1]{inputenc}\n' ]);
  if ~strcmp(epsdirname,viewdirname)
    %disp([ 'warning: The view-file has to be supplemented by '...
    %	   'path information.' ])
    fprintf(fid,[ '\\graphicspath{{' epsdirname '}}\n' ]);
  end  
  fprintf(fid,[ '\\begin{document}\n' ]);
  fprintf(fid,[ '\\pagestyle{empty}\n' ]);
  fprintf(fid,[ '\\begin{figure}[ht]\n' ]);
  fprintf(fid,[ '  \\begin{center}\n' ]);
  if strcmp(texdirname,viewdirname) 
    %fprintf(fid,[ '    \\fbox{\\input{' texbasenameext '}}\n' ]);
    fprintf(fid,[ '    \\input{' texbasenameext '}\n' ]);
  else
    %fprintf(fid,[ '    \\fbox{\\input{' texdirname texbasenameext '}}\n' ]);
    fprintf(fid,[ '    \\input{' texdirname texbasenameext '}\n' ]);
  end
  fprintf(fid,[ '    %% \\caption{A laprint figure}\n' ]);
  fprintf(fid,[ '    %% \\label{fig:' texbasename '}\n' ]);
  fprintf(fid,[ '  \\end{center}\n' ]);
  fprintf(fid,[ '\\end{figure}\n' ]);
  fprintf(fid,[ '\\vfill\n' ]);
  fprintf(fid,[ '\\begin{flushright}\n' ]);
  fprintf(fid,[ '\\tiny printed with laprint on ' ...
		datestr(now) '\\\\\n' ]);
  fprintf(fid,[  viewdirname viewbasenameext '\\\\\n' ]);
  fprintf(fid,[  '( ' texdirname texbasenameext ' )\\\\\n' ]);
  fprintf(fid,[  '( ' epsdirname epsbasenameext ' )\n' ]);
  fprintf(fid,[ '\\end{flushright}\n' ]);
  fprintf(fid,[ '\\end{document}\n' ]);
  fclose(fid);
  if verbose
    yn=input([ 'Perform LaTeX run on ' viewbasenameext '? (y/n) '],'s');
    if strcmp(yn,'y') 
      cmd=[ '!latex ' viewbasenameext ];
      disp([ 'executing: '' ' cmd ' ''' ]);
      eval(cmd);
      yn=input([ 'Perform dvips run on ' viewbasename '.dvi? (y/n) '],'s');
      if strcmp(yn,'y') 
        cmd=[ '!dvips -o' viewbasename '.ps ' viewbasename '.dvi' ];
        disp([ 'executing: '' ' cmd ' ''' ]);
        eval(cmd);
        yn=input([ 'Call ghostview on ' viewbasename '.ps? (y/n) '],'s');
        if strcmp(yn,'y') 
          cmd=[ '!ghostview ' viewbasename '.ps&' ];
          disp([ 'executing: '' ' cmd ' ''' ]);
          eval(cmd);
        end
      end
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%% PART 4: functions used
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fullnameext,basenameext,basename,dirname]= getfilenames(...
    filename,extension,verbose);
% appends an extension to a filename (as '/home/tom/tt') and determines  
%   fullnameext: filename with extension with dirname, as '/home/tom/tt.tex'
%   basenameext: filename with extension without dirname, as 'tt.tex'
%   basename   : filename without extension without dirname, as 'tt'
%   dirname    : dirname without filename, as '/home/tom/'
% In verbose mode, it asks if to overwrite or to modify.
%
[dirname, basename] = splitfilename(filename);
fullnameext = [ dirname basename '.' extension ];
basenameext = [ basename '.' extension ];
if verbose
  quest = (exist(fullnameext)==2);
  while quest
    yn=input([ fullnameext ' exists. Overwrite? (y/n) '],'s');
    if strcmp(yn,'y') 
      quest=0;
    else
      filename=input( ...
	  [ 'Please enter new filename (without extension .' ...
	  extension '): ' ],'s');
      [dirname, basename] = splitfilename(filename);
      fullnameext = [ dirname basename '.' extension ];
      basenameext = [ basename '.' extension ];
      quest = (exist(fullnameext)==2);
    end
  end
end
if ( exist(dirname)~=7 & ~strcmp(dirname,[ '.' filesep ]) ...
      & ~strcmp(dirname,filesep) )
  error([ 'Directory ' dirname ' does not exist.' ] )
end  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dirname,basename]=splitfilename(filename);
% splits filename into dir and base
slashpos=findstr(filename,filesep);
nslash=length(slashpos);
nfilename=length(filename);
if nslash
  dirname = filename(1:slashpos(nslash));
  basename = filename(slashpos(nslash)+1:nfilename);
else  
  dirname = [ pwd filesep ] ;
  basename = filename;
end  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function yesno=is3d(haxes);
% tries to figure out if axes is 3D
yesno=0;
CameraPosition=get(haxes,'CameraPosition');
CameraTarget=get(haxes,'CameraTarget');
CameraUpVector=get(haxes,'CameraUpVector');
if CameraPosition(1)~=CameraTarget(1)
  yesno=1;
end  
if CameraPosition(2)~=CameraTarget(2)
  yesno=1;
end  
if any(CameraUpVector~=[0 1 0])
  yesno=1;
end  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function b=celltoarray(a);
% converts a cell of doubles to an array
if iscell(a),
  b=[];
  for i=1:length(a),
    b=[b a{i}]; 
  end  
else, 
  b=a(:)';
end  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function b=chartocell(a)
% converts a character array into a cell array of characters

% convert to cell 
if isa(a,'char')
  n=size(a,1);
  b=cell(1,n);
  for j=1:n
    b{j}=a(j,:); 
  end  
else
  b=a;
end  
% convert to char
n=length(b);
for j=1:n
  if isa(b{j},'double')
    b{j}=num2str(b{j});
  end  
end	

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function b=overwritetail(a,k)
% overwrites tail of a by k
% a,b: strings
% k: integer
ks=int2str(k);
b = [ a(1:(length(a)-length(ks))) ks ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 