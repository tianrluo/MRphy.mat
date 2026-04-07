function [d, weight] = target_ovs(cube, beta_iv, beta_ov, iv, ov, varargin)
  % Compute the target steady-state magnetization profile for an OVS SPGR sequence.
  % The target is the desired SS magnetization right before the alpha pulse,
  % specified via target beta flip angles for inner- and outer-volume regions.
  %
  % Usage:
  %   [d, weight] = mrphy.steady_state.spgr.target_ovs(cube, 0, 90, iv, ov, 'alpha', 20, 'TR', 80e-3)
  %
  %INPUTS:
  % - cube    (1,) mrphy.SpinCube
  % - beta_iv (1,) deg, target inner-volume beta flip angle
  % - beta_ov (1,) deg, target outer-volume beta flip angle
  % - iv      (*Nd) logical, inner-volume mask
  % - ov      (*Nd) logical, outer-volume mask
  %OPTIONALS:
  % - weight_iv (1,) [1.0], loss weight for iv region
  % - weight_ov (1,) [1.0], loss weight for ov region
  % - doEmbed   [T/f], return embedded or compact form
  % - alpha     (1,) deg [0],     uniform excitation flip angle
  % - TR        (1,) Sec [55e-3], repetition time
  %OUTPUTS:
  % - d      (nM, xyz) | (*Nd, xyz), target SS magnetization
  % - weight (nM, 1)   | (*Nd, 1),  voxel weights
  import attr.*

  [arg.weight_iv, arg.weight_ov] = deal(1.0, 1.0);
  [arg.doEmbed, arg.alpha, arg.TR] = deal(true, 0, 55e-3);
  arg = attrParser(arg, varargin);

  iv_ = cube.extract(iv);  % (nM, 1)
  ov_ = cube.extract(ov);  % (nM, 1)

  M0 = cube.M_(:, 3);             % (nM, 1), equilibrium z-magnetization
  E1 = exp(-arg.TR ./ cube.T1_);  % (nM, 1)

  beta_rad  = deg2rad(beta_iv .* iv_ + beta_ov .* ov_);  % (nM, 1)
  alpha_rad = deg2rad(arg.alpha);

  denom = 1 - cos(beta_rad) .* cos(alpha_rad) .* E1;
  scale = M0 .* (1 - E1) ./ denom;

  d_ = zeros(cube.nM, 3);
  d_(:, 3) = scale .* cos(beta_rad);  % Mz before alpha
  d_(:, 1) = scale .* sin(beta_rad);  % Mxy before alpha (along x)
  d_(isnan(d_)) = 0;

  weight_ = arg.weight_iv .* iv_ + arg.weight_ov .* ov_;  % (nM, 1)

  if arg.doEmbed
    d      = cube.embed(d_);
    weight = cube.embed(weight_);
  else
    d      = d_;
    weight = weight_;
  end
end
