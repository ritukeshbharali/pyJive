init =
{
  nodeGroups = [ lb, rb ];

  mesh = 
  {
    type = gmsh;
    file = solid.msh;
  };

  lb =
  {
    xtype = min;
    ytype = min;
  };

  rb =
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
      type = Elastic;
      E = 20e9;
      nu = 0.2;
      rank = 2;
      rho = 2500;
      anmodel = plane_stress;
    };

    thickness = 0.2;

    shape =
    {
      type = Triangle3;
      intScheme = Gauss1;
    };
  };

  diri =
  {
    type = Dirichlet; 

    groups = [ lb, rb, rb ];
    dofs   = [ dy, dx, dy ];
    values = [ 0.0, 0.0, 0.0 ];
  };
};

modeshape =
{
  type = ModeShape;
};

view = 
{
  type = View;
  deform = 2;
  plot = solution[dy];
  interactive = True;
  maxStep = 50;
  colorMap = plasma_r;
};
  
