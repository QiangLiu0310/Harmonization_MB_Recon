function noisepower = mrir_array_SNR_noisepower(sensitivity, covmtx, varargin)
%MRIR_ARRAY_SNR_NOISEPOWER
%
% noisepower = mrir_array_SNR_noisepower(sensitivity, covmtx)

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/mar/08
% $Id: mrir_array_SNR_noisepower.m,v 1.1 2007/06/01 21:48:34 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;

  FLAG__PREWHITEN = 0;

  
  %==--------------------------------------------------------------------==%
  
  if ( nargin > 3 ),
    intensities = varargin{1};
  else,
    intensities = sensitivity;
  end;
  
  
  %==--------------------------------------------------------------------==%

  Ncol = size(sensitivity, mrir_DIM_COL);
  Nlin = size(sensitivity, mrir_DIM_LIN);
  Ncha = size(sensitivity, mrir_DIM_CHA);
  Npar = size(sensitivity, mrir_DIM_PAR);
  Nslc = size(sensitivity, mrir_DIM_SLC);

  if ( FLAG__PREWHITEN ),
    % whitening operator for pre-whitening data
    W = mrir_array_whitening_operator(covmtx, 'svd');

    sensitivity_whitened = reshape(sensitivity, [], size(W,1)) * W.';
    sensitivity = reshape(sensitivity_whitened, size(sensitivity));

    intensities_whitened = reshape(intensities, [], size(W,1)) * W.';
    intensities = reshape(intensities_whitened, size(intensities));
    
  else,
    covmtxinv = inv(covmtx);
  end;

  % preallocate         1     2  3  4  5  6  7  8     9     0
  noisepower = zeros(Ncol, Nlin, 1, 1, 1, 1, 1, 1, Npar, Nslc);
  
  for kk = 1:Npar,
  for ii = 1:Ncol,
    for jj = 1:Nlin,
      S = squeeze(sensitivity(ii, jj, :,1,1,1,1,1,kk));
      I = squeeze(intensities(ii, jj, :,1,1,1,1,1,kk));

      if ( FLAG__PREWHITEN ),
	noisepower_inv = S' * I;
      else,
	noisepower_inv = S' * covmtxinv * I;
      end;
      
      noisepower(ii, jj, 1,1,1,1,1,1,kk) = abs( inv(noisepower_inv) );
    end;
  end;
  end;
  
  % recall that oversampling does not affect SNR map, but normalization
  % factor already applied to sensitivity map
  scale_factor = 1 / sqrt(2 * Ncol * Nlin * Npar);
  
  noisepower = scale_factor * noisepower;

  
  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_array_SNR_noisepower.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End: