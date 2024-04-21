init =
{
  nodeGroups = [ botleft, top, botright ];

  mesh = 
  {
    type = geo;
    file = vierendeel.geom;
  };

  botleft = [0];
  botright = [6];
  top = [7,8,9,10,11,12,13];
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

    groups = [ botleft, botleft, botright];
    dofs   = [ dx, dy, dy ];
    values = [ 0.0, 0.0, 0.0 ];
  };
  
  neum = 
  {
    type = Neumann;
    
    groups = [ top ];
    dofs = [ dy ];
    values = [ -1. ];
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
};
