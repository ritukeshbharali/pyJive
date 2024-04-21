init =
{
  nodeGroups = [ left, right ];

  mesh =
  {
    type = manual;
    file = timoshenko.mesh;
  };

  left =
  {
    xtype = min;
  };

  right =
  {
    xtype = max;
  };
};

model =
{
  type = Multi;

  models = [ timoshenko, diri, neum ];

  timoshenko =
  {
    type = Timoshenko;

    elements = all;

    EI = 20000;

    GAs = 100000;

    shape =
    {
      type = Line2;
      intScheme = Gauss2;
    };
  };

  diri =
  {
    type = Dirichlet; 

    groups = [ left, left ];
    dofs   = [ phi, dy ];
    values = [ 0.0, 0.0 ];
  };

  neum = 
  {
    type = Neumann;

    groups = [ right ];
    dofs = [ dy ];
    values = [ 1.0 ];
  };
};

solver =
{
  nsteps = 1;
};
