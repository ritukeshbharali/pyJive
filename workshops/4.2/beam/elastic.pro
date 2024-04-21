init =
{
  nodeGroups = [ left, right, mid ];

  mesh = 
  {
    type = geo;
    file = beam.geom;
  };

  left = [0];
  mid = [1];
  right = [2];
};

model =
{
  type = Multi;

  models = [ frame, diri, neum ];

  frame =
  {
    type = Frame;

    elements = all;

    subtype = nonlin;

    EA = 2000;
    GAs = 1000;
    EI = 1;

    shape =
    {
      type = Line2;
      intScheme = Gauss1;
    };
  };

  diri =
  {
    type = Dirichlet; 

    groups = [ left, right, right ];
    dofs   = [ dy, dx, dy ];
    values = [ 0.0, 0.0, 0.0 ];
  };

  neum = 
  {
    type = Neumann;

    groups = [ left, mid ];
    dofs = [ dx, dy ];
    values = [ 0., 0. ];
    loadIncr = [ 0.999, 0.001 ];
  };
};

nonlin =
{
  type = Arclen;
  nsteps = 50;
  itermax = 10;
  tolerance = 1e-6;
  beta = 0.1;
  dl = 0.05;
};

loaddisp = 
{
  type = LoadDisp;
  groups = [ left, mid, right ];
};

frameview =
{
  type = FrameView;
  deform = 1.;
  interactive = True;
  plotStress = M;
};

