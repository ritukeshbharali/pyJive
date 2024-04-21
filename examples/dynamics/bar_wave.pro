init =
{
  nodeGroups = [ left, right ];

  mesh = 
  {
    type = geo;
    file = bar_5elem.geom;
  };

  left = [0];
  right = [1];
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

    groups = [ left, left, left ];
    dofs   = [ dx, dy, phi ];
    values = [ 0.0, 0.0, 0.0 ];
  };
  
  neum =
  {
    type = Neumann;
    
    groups = [right];
    dofs = [ dx ];
    values = [1];
    loadIncr = [ 0.01 ];
  };

};


explicittime =
{
    type = ExplicitTime;
    deltaTime = 0.05;
    nsteps = 200;
    storeMatrix = True;
    storeMassMatrix = True;
    tolerance = 1;
};

frameview =
{
  type = FrameView;
  deform = 1.;
  step0 = 0;
  maxStep = 200;
  vectorSize = 0.5;
};
