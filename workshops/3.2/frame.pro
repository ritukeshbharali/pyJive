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

    groups = [ botleft, botleft, botright, botright];
    dofs   = [ dx, dy, dx, dy ];
    values = [ 0.0, 0.0, 0.0, 0.0 ];
  };
  
  neum = 
  {
    type = Neumann;
    
    groups = [ topleft ];
    dofs = [ dx ];
    values = [ 1. ];
  };
};

solver =
{
  type = Solver;
  nsteps = 1;
};

frameview =
{
  type = FrameView;
  deform = 1.;
  interactive = False;
  plotStress = M;
};
