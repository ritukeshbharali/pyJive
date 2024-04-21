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

  models = [ frame, diri ];

  frame =
  {
    type = Frame;

    elements = all;

    subtype = nonlin;

    EA = 10000;
    GAs = 5000;
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

    groups = [ botleft, botleft, botright, botright, topright ];
    dofs   = [ dx, dy, dx, dy, dy ];
    values = [ 0.0, 0.0, 0.0, 0.0, 0.0 ];
    dispIncr = [ 0.0, 0.0, 0.0, 0.0, -0.00005 ];
  };
};

nonlin =
{
  type = Nonlin;
  nsteps = 150;
  itermax = 10;
  tolerance = 1e-6;
};

frameview =
{
  type = FrameView;
  deform = 1.;
  interactive = True;
  step0 = 150;
};

loaddisp = 
{
  type = LoadDisp;
  groups = [ topright ];
};

graph = 
{
  type = Graph;
  xData = [loaddisp.topright.disp.dy];
  yData = [loaddisp.topright.load.dy];
  xlabel = u;
  ylabel = F;
};
