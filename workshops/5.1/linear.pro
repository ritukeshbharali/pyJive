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

    subtype = linear;
    plastic = False;

    EA = 100000;
    GAs = 50000;
    EI = 10;
    Mp = 1;

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
    
    groups = [ topleft, topleft, topright ];
    dofs = [ dx, dy, dy ];
    values = [ 1, -5, -5 ];
    loadIncr = [ 1, -5, -5 ];
  };
};

nonlin =
{
  type = Nonlin;
  nsteps = 50;
  itermax = 10;
  tolerance = 1e-6;
};

frameview =
{
  type = FrameView;
  deform = 1.;
  interactive = True;
  plotStress = M;
  step0 = 49;
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
