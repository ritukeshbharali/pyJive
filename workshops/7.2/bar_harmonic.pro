init =
{
  nodeGroups = [ left, right, all ];

  mesh = 
  {
    type = geo;
    file = beam.geom;
  };

  left = [0];
  right = [1];
  all = {};
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

    EA = 10;
    GAs = 1;
    EI = 1;
    rhoA = 1;
    rhoI = 2;

    shape =
    {
      type = Line2;
      intScheme = Gauss1;
    };
  };

  diri =
  {
    type = Dirichlet; 

    groups = [ left, all, all ];
    dofs   = [ dx, dy, phi ];
    values = [ 0.0, 0.0, 0.0 ];
  };
  
  neum = 
  {
    type = Neumann;
    
    groups = [ right ];
    dofs = [ dx ];
    values = [ 1 ];
    deltaTime = 0.01;
    timeSignal = np.sin(t)**2;
  };
};

stepper =
{
  type = Newmark;
  deltaTime = 0.01;
  gamma = 0.5;
  beta = 0.25;
  nsteps = 2000;
};

lodi = 
{
  type = LoadDisp;
  groups = [left, right];
};

frameview =
{
  type = FrameView;
  deform = 5.;
  plotDirichlet = False;
  plotStress = N;
  vectorSize = 0.5;
};
