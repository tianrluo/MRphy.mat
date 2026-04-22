function [M_ss, Mhst_] = spgr_ovs(spin, pulse, varargin)
  % Compute the steady-state magnetization for an outer-volume-suppression (OVS) SPGR sequence.
  % Assumes the following sequence structure:
  %   [beta - alpha - readout] x N_rep
  % where beta (the input pulse) is a saturation preparation pulse and alpha is
  % a non-selective excitation. If alpha=0, the pulse is treated as a selective
  % excitation rather than a saturation preparation.
  % The output M_ss is the steady-state magnetization right before the alpha pulse.
  %
  % Usage:
  %   [M_ss, Mhst_] = mrphy.steady_state.spgr.spgr_ovs(spin, pulse, 'alpha', 20, 'TR', 80e-3)
  %
  %INPUTS:
  % - spin  (1,) mrphy.SpinArray or mrphy.SpinCube
  % - pulse (1,) mrphy.Pulse
  %OPTIONALS:
  % - TR      (1,) Sec [55e-3], repetition time
  % - alpha   (1,) deg [0],     uniform excitation flip angle
  % - loc_   ^ loc   (nM, xyz)       ^ (*Nd, xyz), XOR; ignored for SpinCube
  % - b0Map_ | b0Map (nM, 1, nCoils) ^ (*Nd, 1, nCoils); ignored for SpinCube
  % - b1Map_ | b1Map (nM, 1, nCoils) ^ (*Nd, 1, nCoils)
  % - doCim    [T/f]
  % - doEmbed  [t/F]
  % - doUpdate [t/F], update spin.M_ to the SS magnetization
  %OUTPUTS:
  % - M_ss  (nM, xyz) | (*Nd, xyz)
  % - Mhst_ (nM, xyz, nT) | (*Nd, xyz, nT), evolving history
  import attr.*

  [arg.loc,   arg.loc_]   = deal([], []);
  [arg.b0Map, arg.b0Map_] = deal([], []);
  [arg.b1Map, arg.b1Map_] = deal([], []);
  [arg.doCim, arg.doEmbed, arg.doUpdate] = deal(true, false, false);
  [arg.alpha, arg.TR] = deal(0, 55e-3);
  arg = attrParser(arg, varargin);

  if isa(spin, 'mrphy.SpinCube')
    kw = {  'loc_',spin.loc_, 'b0Map_',spin.b0Map_ ...
          , 'b1Map',arg.b1Map, 'b1Map_',arg.b1Map_ ...
          , 'doEmbed',false};
    sa = spin.spinarray;
  else
    kw = {  'loc',arg.loc,     'loc_',arg.loc_ ...
          , 'b0Map',arg.b0Map, 'b0Map_',arg.b0Map_ ...
          , 'b1Map',arg.b1Map, 'b1Map_',arg.b1Map_ ...
          , 'doEmbed',false};
    sa = spin;
  end

  beff_ = sa.pulse2beff(pulse, kw{:});
  [Mo_, Mhst_] = mrphy.sims.blochsim(sa.M_, beff_, sa.T1_, sa.T2_ ...
                                     , pulse.dt, sa.gam_ ...
                                     , arg.doCim);

  Mo_(isnan(Mo_)) = 0;
  beta_ = atan2((Mo_(:,1).^2 + Mo_(:,2).^2).^0.5, Mo_(:,3));  % saturation angle (nM,1)

  E1 = exp(-arg.TR ./ sa.T1_);
  denom = 1 - cos(beta_) .* cosd(arg.alpha) .* E1;
  scale = sa.M_(:,3) .* (1 - E1) ./ denom;

  M_ss = zeros(size(Mo_));
  M_ss(:,3) = scale .* cos(beta_);  % steady-state Mz before alpha
  M_ss(:,1) = scale .* sin(beta_);  % steady-state Mxy before alpha (along x)

  if arg.doUpdate, sa.M_ = M_ss; end
  if arg.doEmbed
    [M_ss, Mhst_] = deal(spin.embed(M_ss), spin.embed(Mhst_));
  end
end
