init =
{
  nodeGroups = [ left, right, mid ];

  mesh = 
  {
    type = geo;
    file = euler.geom;
  };

  left = [0];
  mid = [1];
  right = [2];
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

    EA = 20000;
    GAs = 10000;
    EI = 10;

    shape =
    {
      type = Line2;
      intScheme = Gauss1;
    };
  };

  diri =
  {
    type = Dirichlet; 

    groups = [ left, left, right, right ];
    dofs   = [ dy, dx, dy, dx ];
    values = [ 0.0, 0.0, 0.0, 0.0 ];
    dispIncr = [ 0.0, 0.0001, 0.0, 0.0 ];
  };

  neum = 
  {
    type = Neumann;

    groups = [ mid ];
    dofs = [ dy ];
    values = [ 0.001 ];
    loadIncr = [ 0.0 ];
  };
};

nonlin =
{
  type = Nonlin;
  nsteps = 100;
  itermax = 20;
  tolerance = 1e-6;
};

loaddisp = 
{
  type = LoadDisp;
  groups = [ left, mid, right ];
};

frameview =
{
  type = FrameView;
  deform = 1.;
  interactive = True;
  vectorSize = 1.;
};

graph = 
{
  type = Graph;
  xData = [loaddisp.left.disp.dx,loaddisp.mid.disp.dy];
  yData = [loaddisp.left.load.dx,loaddisp.left.load.dx];
  legend = [u_end, v_mid];
  xlabel = displacement;
  ylabel = force;
};
