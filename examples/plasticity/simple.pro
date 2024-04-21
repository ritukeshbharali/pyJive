init =
{
  type = Init;

  mesh =
  {
    type = gmsh;
    file = simple.msh;
  };

  nodeGroups = [ left, right, leftb, rightb ];

  left =
  {
    xtype = min;
  };

  right =
  {
    xtype = max;
  };

  leftb =
  {
    xtype = min;
    ytype = min;
  };

  rightb =
  {
    xtype = max;
    ytype = min;
  };
};

model =
{
  type = Multi;

  models = [ solid, diri ];

  solid =
  {
    type = Solid;

    elements = all;

    material =
    {
      type = J2;
      rank = 2;
      anmodel = plane_strain;

      E = 1000.;
      nu = 0.2;

      yield = 64.8 - 33.6*exp(kappa/-0.003407);
      maxIter = 1000;
      tolerance = 1.e-7;
    };

    thickness = 1.0;

    shape =
    {
      type = Triangle3;
      intScheme = Gauss1;
    };
  };

  diri =
  {
    type = Dirichlet;

    groups = [ left, right, leftb, rightb ];
    dofs   = [ dx  , dx   , dy   , dy     ];
    values = [ 0.  , .01  , 0.   , 0.     ];
    dispIncr=[ 0.  , .01  , 0.   , 0.     ];
  };
};

nonlin =
{
  nsteps = 20;
  itermax = 10;
  tolerance = 1e-6;
};

lodi =
{
  type = LoadDisp;
  groups = [ right ];
};

