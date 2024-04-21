init =
{
  nodeGroups = [ botleft, topleft, topright, botright ];

  mesh = 
  {
    type = geo;
    file = frame.geom;
  };

  botleft = [0];
  topleft = [1];
  topright = [2];
  botright = [3];
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

    groups = [ botleft, botleft, botright, botright ];
    dofs   = [ dx, dy, dx, dy ];
    values = [ 0.0, 0.0, 0.0, 0.0 ];
    dispIncr = [ 0.0, 0.0, 0.0, 0.0 ];
  };

  neum = 
  {
    type = Neumann;

    groups = [topright];
    dofs = [dy];
    values = [-1.];
    loadIncr = [0.];
  };
};

linbuck =
{
  type = LinBuck;
};

frameview =
{
  type = FrameView;
  deform = 1.;
  step0 = 0;
};

