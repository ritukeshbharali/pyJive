init =
{
  nodeGroups = [ left, right ];

  mesh = 
  {
    type = manual;
    file = bar.mesh;
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

  models = [ bar, diri, neum ];

  bar =
  {
    type = Bar;

    elements = all;

    EA = 1000.0;

    shape =
    {
      type = Line2;
      intScheme = Gauss2;
    };
  };

  diri =
  {
    type = Dirichlet; 

    groups = [ right ];
    dofs   = [ dx ];
    values = [ 0.0 ];
  };

  neum = 
  {
    type = Neumann;

    groups = [ left ];
    dofs = [ dx ];
    values = [ -1.0 ];
  };
};

solver =
{
  nsteps = 1;
  storeMatrix = True;
};
