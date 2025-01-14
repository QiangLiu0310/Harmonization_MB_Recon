function filestr = mrir_sysutil__wildfile(wildfile, varargin)
%MRIR_SYSUTIL__WILDFILE
%
% filestr = mrir_sysutil__wildfile(wildstr, varargin)


% Copyright © 2006-2012 Jonathan R. Polimeni and
%   The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense


% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2009/dec/28
% $Id: mrir_sysutil__wildfile.m,v 1.4 2012/05/26 09:58:27 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.4 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  FLAG__return_all = 0;
  if ( nargin >= 3 ),
    FLAG__return_all = varargin{2};
  end;

  if ( iscellstr(wildfile) ),
    wildfile = cell2mat(wildfile);
  end;



  % expand wild cards
  [pathstr, basestr, ext] = fileparts(wildfile);

  % overwrite pathstr if input argument provided
  if ( nargin >= 2 ),
    if ( ~isempty(varargin{1}) ),
      pathstr = varargin{1};
    end;
  end;

  if ( isempty(ls(fullfile(pathstr, [basestr, ext]))) ),
    disp(sprintf('pattern [%s] does not match existing file\n', wildfile));
    filestr = '';
    return;
  end;

  wildstr = fullfile(pathstr, [basestr, ext]);

  [status, listing] = unix(sprintf('ls --color=never %s', wildstr));

  filestr = deblank(listing);

  if ( ~exist(filestr, 'file') ),
    error('pattern [%s] does not match existing file\n', wildfile);
  end;

  % TODO: add flag for 1st hit versus all hits a la "wildfile.m"


  return;



  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_sysutil__wildfile.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:


