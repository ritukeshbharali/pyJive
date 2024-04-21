init =
{
  nodeGroups = [ left, right, mid ];

  mesh = 
  {
    type = geo;
    file = euler.geom;
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

    subtype = linear;

    EA = 20000;
    GAs = 10000;
    EI = 10;

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

    groups = [ left ];
    dofs = [ dx ];
    values = [ 1. ];
  };
};

linbuck =
{
  type = LinBuck;
};

frameview =
{
  type = FrameView;
  deform = 10.;
  label = Mode;
  maxStep = 20;
  step0 = 0;
};
