init =
{
  nodeGroups = [ lb, rb, tm ];

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
  
  tm =
  {
    xtype = mid;
    ytype = max;
  };
};

model =
{
  type = Multi;

  models = [ solid, dirichlet ];

  solid =
  {
    type = Solid;

    elements = all;
    
    material = 
    {
      type = Elastic;
      E = 10000.;
      nu = 0.2;
      rank = 2;
      rho = 1;
      anmodel = plane_stress;
    };
    
    thickness = 0.2;

    shape =
    {
      type = Triangle3;
      intScheme = Gauss1;
    };
  };

  dirichlet =
  {
    type = Dirichlet; 

    groups = [ lb, lb, rb ];
    dofs   = [ dx, dy, dy ];
    values = [ 0., 0., 0. ];
  };
  

};

solver =
{
  nsteps = 1;
  storeMatrix = True;
};
