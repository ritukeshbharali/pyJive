init =
{
  nodeGroups = [ botleft, botright, topleft ];

  mesh = 
  {
    type = geo;
    file = portal.geom;
  };

  botleft = [0];
  botright = [3];
  topleft = [1];
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

    EA = 1;
    GAs = 1;
    EI = 1;
    rhoA = 1;
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

    groups = [ botleft, botleft, botleft, botright, botright, botright];
    dofs   = [ dx, dy, phi, dx, dy, phi ];
    values = [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ];
  };
  
  neum =
  {
    type = Neumann;
    
    groups = [topleft];
    dofs = [ dx ];
    values = [0.1];
    loadIncr = [0];
  };

};


explicittime =
{
    type = ExplicitTime;
    deltaTime = 0.05;
    nsteps = 400;
    storeMatrix = True;
    storeMassMatrix = True;
};

frameview =
{
  type = FrameView;
  deform = 1.;
  step0 = 0;
  maxStep = 400;
  vectorSize = 0.5;
};
