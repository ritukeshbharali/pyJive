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
    plastic = True;

    EA = 20000;
    GAs = 10000;
    EI = 10;
    Mp = 0.4;

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
    values = [ 0, 0, 0, 0 ];
  };
  
  neum = 
  {
    type = Neumann;
    groups = [ topleft, topleft, topright ];
    dofs = [ dx, dy, dy ];
    values = [ 0, 0, 0 ];
    loadIncr = [ 1, -5, -5];
  };
};

nonlin =
{
  type = Arclen;
  nsteps = 100;
  itermax = 10;
  tolerance = 1e-6;
  beta = 0.1;
  dl = 0.02;
};

frameview =
{
  type = FrameView;
  deform = 1.;
  interactive = True;
};

loaddisp = 
{
  type = LoadDisp;
  groups = [ topleft ];
};

graph = 
{
  type = Graph;
  xData = [loaddisp.topleft.disp.dx];
  yData = [loaddisp.topleft.load.dx];
  xlabel = u;
  ylabel = F;
};
