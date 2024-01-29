ENTANGLEMENT_VARS = {
  'Rchsh'        : { 'threshold' : 1, 'symbol' : r'\mathfrak{m}_{12}[C]' },
  'concurrence'  : { 'threshold' : 0, 'symbol' : r'\mathcal{C}[\rho]' },
  'steerability' : { 'threshold' : 1, 'symbol' : r'S' },
}

CENTRAL_CHOICES = [ 'nominal', 'median' ]

AXIS_CHOICES = [ 'beam', 'higgs' ]

SPIN_ANALYZERS = {
  'asymmetry'          : 'FB asymm.',
  'differentialXsec1d' : '1d distr.',
  'differentialXsec2d' : '2d distr.',
  'mlfit'              : 'ML-fit',
  'summation'          : 'Exp. value',
}

DECAY_MODES = {
  'pi_pi'   : { 'label' : r'$\pi^+\pi^-$',       'marker' : 'o' },
  'pi_rho'  : { 'label' : r'$\pi^\pm\rho^\mp$',  'marker' : 'v' },
  'pi_a1'   : { 'label' : r'$\pi^\pm a_1^\mp$',  'marker' : '^' },
  'rho_rho' : { 'label' : r'$\rho^+\rho^-$',     'marker' : 's' },
  'rho_a1'  : { 'label' : r'$\rho^\pm a_1^\mp$', 'marker' : 'X' },
  'a1_a1'   : { 'label' : r'$a_1^+a_1^-$',       'marker' : 'd' },
  'comb'    : { 'label' : 'Combination',         'marker' : '*' },
#  'had_had' : { 'label' : 'All hadronic',        'marker' : 'r$\club$' },
}

INIT_METHODS = {
  'gen'      : 'MC-truth level',
#  'startPos' : 'Analytic',
  'kinFit'   : 'KF after smearing',
  'svFit'    : 'SVfit after smearing',
}

EXTENSIONS = [ 'pdf', 'png' ]
