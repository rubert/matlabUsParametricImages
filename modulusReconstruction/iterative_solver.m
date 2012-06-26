function mu = iterative_solver(tempDp, inAlpha)


global mu_0
global alpha
global dp_m
global residArray

residArray = zeros(100,1);

mu_0 = 1;
dp_m = tempDp;
alpha = inAlpha;


%work out number of four node rectangular elements
[rows, cols] = size(dp_m);  
rows_el = rows - 1; cols_el = cols - 1;
num_elems =  rows_el*cols_el;
init_modulus = ones(rows_el,cols_el);  %%if I write the modulus distribution as mu_0*f, this is f


lb = (1E-6)*ones(num_elems,1);   %lower bound on modulus
ub = (1E6)*ones(num_elems,1);   %upper bound on modulus

options = optimset('algorithm', 'trust-region-reflective', 'hessian', 'user-supplied', 'gradobj', 'on', 'display', 'iter', 'TolFun', 1E-8);

mu = fmincon('computeAll', init_modulus, [], [], [], [], lb,ub,[], options);