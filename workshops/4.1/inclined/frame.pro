init =
{
  nodeGroups = [ bl, br, top ];

  mesh = 
  {
    type = geo;
    file = frame.geom;
  };

  bl = [2];
  top = [1];
  br = [0];
};

model =
{
  type = Multi;

  models = [ column, spring, diri, neum ];

  column =
  {
    type = Frame;

    elements = member0;

    subtype = nonlin;

    EA = 100000;
    GAs = 20000;
    EI = 1.e6;

    shape =
    {
      type = Line2;
      intScheme = Gauss1;
    };
  };

  spring =
  {
    type = Frame;

    elements = member1;

    subtype = nonlin;

    EA = 14.14213562;
    GAs = 0;
    EI = 0;

    shape =
    {
      type = Line2;
      intScheme = Gauss1;
    };
  };

  diri =
  {
    type = Dirichlet; 

    groups = [ bl, bl, bl, br, br ];
    dofs   = [ dx, dy, phi, dx, dy ];
    values = [ 0.0, 0.0, 0.0, 0.0, 0.0 ];
    dispIncr = [ 0.0, 0.0, 0.0, 0.0, 0.0 ];
  };

  neum = 
  {
    type = Neumann;

    groups = [ top ];
    dofs = [ dy ];
    values = [ -0.0 ];
    loadIncr = [ -0.2 ];
  };
};

arclen =
{
  type = Arclen;
  nsteps = 100;
  itermax = 10;
  tolerance = 1e-6;
  beta = 0.2;
  dl = 0.02;
};

loaddisp = 
{
  type = LoadDisp;
  groups = [ top ];
};

graph =
{
  type = Graph;
  xData = [loaddisp.top.disp.dx,loaddisp.top.disp.dy];
  yData = [loaddisp.top.load.dy,loaddisp.top.load.dy];
};

frameview = 
{
  type = FrameView;
  interactive = True;
  deform = 1.;
};
