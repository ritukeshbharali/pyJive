init =
{
  type = Init;

  mesh =
  {
    type = gmsh;
    file = voids.msh;
  };

  nodeGroups = [ right, ???, ??? ];

  right =
  {
    xtype = max;
  };

  ??? =
  {
    ??? = ???;
  };
  
  ??? =
  {
    ??? = ???;
  };
};

model =
{
  type = Multi;

  models = [ solid, dispcontrol ];

  solid =
  {
    type = Solid;

    elements = all;

    material =
    {
      type = J2;
      rank = 2;
      anmodel = plane_stress;

      E = 1000.;
      nu = 0.2;

      yield = ???;
      maxIter = 1000;
      tolerance = 1.e-10;
    };

    thickness = 1.0;

    shape =
    {
      type = Triangle3;
      intScheme = Gauss1;
    };
  };

  dispcontrol =
  {
    type = ???;

    groups = [ ???, ???, right ];
    dofs   = [ ???, ???, dx    ];
    values = [ ???, ???, 1.0   ];
    timeSignal = ???;
    deltaTime = 0.01;
  };
};

nonlin =
{
  type = Nonlin;
  nsteps = 100;
  itermax = 100;
  tolerance = 1e-6;
  lenient = False;
};

lodi =
{
  type = LoadDisp;
  groups = [???];
};

graph =
{
  xData = [lodi.???.disp.???];
  yData = [lodi.???.load.???];
};
