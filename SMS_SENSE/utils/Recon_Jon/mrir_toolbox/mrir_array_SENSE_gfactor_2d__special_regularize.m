function gfactor = mrir_array_SENSE_gfactor_2d__special(sensitivity, covmtx, R1, R2, varargin)
%MRIR_ARRAY_SENSE_GFACTOR_2D__SPECIAL
%
% gfactor = mrir_array_SENSE_gfactor_2d__special(sensitivity, covmtx, R1, R2)

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/apr/04
% $Id: mrir_array_Gfactor_2d.m,v 1.3 2007/05/14 02:37:50 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.3 $';
  if ( nargin == 0 ), help(mfilename); return; end;

  f = 0.0000005;

  if ( nargin >= 4 ),
    f = varargin{1};
  end;


  %==--------------------------------------------------------------------==%

  Ncol = size(sensitivity, 1); % frequency encoded
  Nlin = size(sensitivity, 2); % phase encoded
  Ncha = size(sensitivity, 3);

  covmtxinv = inv(covmtx);

  aliasmap1 = mrir_array_accelerated_aliasmap(Nlin, R1);
  aliasmap2 = mrir_array_accelerated_aliasmap(Ncol, R2);

  if ( prod(size(aliasmap1)) ~= Nlin )
    warning('R1 not compatible with Nlin');
  end;

  if ( prod(size(aliasmap2)) ~= Ncol )
    warning('R2 not compatible with Ncol');
  end;

  % preallocate
  gfactor  = zeros(Ncol, Nlin);

  for jj = 1:floor(Nlin/R1),

    % indices of the aliased pixels for this position from lookup table
    aliased_ind1 = aliasmap1(jj, 1:R1);

    for ii = 1:floor(Ncol/R2),

      % indices of the aliased pixels for this position from lookup table
      aliased_ind2 = aliasmap2(ii, 1:R2);

      % R2 x R1 x Nchan
      b = sensitivity(aliased_ind2, aliased_ind1, :);

      % encoding matrix: Nchan x R1*R2, or Nchan x Nalias
      E = permute(reshape(b, [R1*R2, Ncha]), [2, 1]);

      errcovmtx_inv = E' * covmtxinv * E;

      e1 = (max(eig(errcovmtx_inv)));
      e2 = (min(eig(errcovmtx_inv)));

      errcovmtx = inv( errcovmtx_inv + f*e1*eye(size(errcovmtx_inv)) );

      U = errcovmtx * E' * covmtxinv;
      rho = U * E;
      A = diag(1./diag(rho));

      nu = A * rho * errcovmtx * A;

      gfactor(aliased_ind2,aliased_ind1) = reshape(sqrt(abs(diag(nu).*diag(errcovmtx_inv))), [R2,R1]);

      eigmax(aliased_ind2,aliased_ind1) = e1;
      eigmin(aliased_ind2,aliased_ind1) = e2;

    end;
  end;


  return;



  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_array_Gfactor_2d.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End: