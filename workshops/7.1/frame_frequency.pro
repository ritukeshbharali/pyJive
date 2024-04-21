init =
{
  nodeGroups = [ bot, top ];

  mesh = 
  {
    type = geo;
    file = house.geom;
  };
  
  bot = [0, 4];   
  top = [2];
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

    EA = 20000;
    GAs = 10000;
    EI = 10;
    rhoA = 2000;
    rhoI = 1;

    shape =
    {
      type = Line2;
      intScheme = Gauss1;
    };
  };

  diri =
  {
    type = Dirichlet; 

    groups = [ bot, bot, bot ];
    dofs   = [ dy, dx, phi ];
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
  deform = 1.;
  label = Mode;
  maxStep = 20;
  step0 = 0;
};
