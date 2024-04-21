init =
{
  type = Init;

  mesh =
  {
    type = gmsh;
    file = bar.msh;
  };

  nodeGroups = [ left, right ];

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

  models = [ solid, diri, neum ];

  solid =
  {
    type = Solid;

    elements = all;

    area = 1.0;

    material =
    {
      type = J2;
      rank = 1;
      anmodel = bar;

      E = 1000.;

      yield = 50.0 + 50.0*kappa;
      maxIter = 1000;
      tolerance = 1.e-7;
    };


    shape =
    {
      type = Line2;
      intScheme = Gauss1;
    };
  };

  diri =
  {
    type = Dirichlet;

    groups = [ left ];
    dofs   = [ dx ];
    values = [ 0. ];
  };

  neum =
  {
    type = Neumann;
    groups = [ right ];
    dofs = [ dx ];
    values = [ 11.0 ];
    timeSignal = t;
    deltaTime = 1.0;
  };
};

nonlin =
{
  type = Nonlin;
  nsteps = 10;
  itermax = 1000;
  tolerance = 1e-6;
  lenient = False;
};

lodi =
{
  type = LoadDisp;
  groups = [right];
};

graph =
{
  xData = [lodi.right.disp.dx];
  yData = [lodi.right.load.dx];
};
