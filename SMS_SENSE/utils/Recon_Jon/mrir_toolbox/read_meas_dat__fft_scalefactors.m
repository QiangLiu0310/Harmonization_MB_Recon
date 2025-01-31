function [fft_scale_sparse, channels] = read_meas_dat__fft_scalefactors(header);
%READ_MEAS_DAT__FFT_SCALEFACTORS
%
% fft_scale = read_meas_dat__fft_scalefactors(header);
%
% NOTE: to match ICE, FFT scale factors should ONLY be applied to data
% that will be Fourier transformed, i.e., **NOT** noise data!

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/may/17
% $Id$
%**************************************************************************%

  VERSION = '$Revision: 1.5 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %------------------------------------------------------------------------%
  % extract individual FFT scale factors for each coil in an array

  fft_scalefactor_str = 'aFFT_SCALE';
  fft_scalefactor_value_match   = regexp(header, ['(?<param>' fft_scalefactor_str ')\[(?<channel>\d*)\]\.flFactor\s*=\s*(?<value>\d*\.?\d*)'], ...
                                         'names');

  indices = str2double({fft_scalefactor_value_match.channel});
  values  = str2double({fft_scalefactor_value_match.value});

  % if scale factors repeated in header, take value of first instance
  [indices, instance] = unique(indices, 'first');

  % channel index is 0-based
  fft_scale(indices+1) = values(instance);


  %------------------------------------------------------------------------%
  % extract corresponding channel labels; in some cases (e.g., 7T host),
  % channel labels are not consecutive.

  fft_scalefactor_channel_match = regexp(header, ['(?<param>' fft_scalefactor_str ')\[(?<channel>\d*)\]\.lRxChannel\s*=\s*(?<value>\d*\.?\d*)'], ...
                                         'names');

  indices = str2double({fft_scalefactor_channel_match.channel});
  values  = str2double({fft_scalefactor_channel_match.value});

  % if scale factors repeated in header, take value of first instance
  [indices, instance] = unique(indices, 'first');

  % channel index is 0-based
  channels(indices+1) = values(instance);


  %------------------------------------------------------------------------%
  % map scale factors to sparse array as a lazy way of dealing with
  % non-consecutive channel labels

  fft_scale_sparse(channels) = fft_scale;


  return;



  %************************************************************************%
  %%% $Source: /home/jonp/cvsroot/dotfiles/emacs,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
