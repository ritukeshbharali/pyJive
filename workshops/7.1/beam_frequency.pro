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

  models = [ frame, diri ];

  frame =
  {
    type = Frame;

    elements = all;

    subtype = linear;

    EA = 8e9;
    GAs = 4e9;
    EI = 2.667e9;
    rhoA = 1000;
    rhoI = 333.3;

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
};

modeshape =
{
  type = ModeShape;
};

frameview =
{
  type = FrameView;
  deform = 5.;
  label = Mode;
  maxStep = 20;
  step0 = 0;
};
