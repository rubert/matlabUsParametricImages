function [resid, grad, hess] = computeAll(func)
%%%%going to work out the functional, gradient, and Hessian in this code.

global alpha;  
global dp_m;


%%%flattened and arranged in alternating array, going x,y,x,y
[mesh_rows, mesh_cols] = size(dp_m);
num_nodes = mesh_rows*mesh_cols;


el_rows = mesh_rows - 1;
el_cols = mesh_cols - 1;
num_elems = el_rows*el_cols;


ind = 1;
nodeNo = zeros(mesh_rows,mesh_cols);


for xx = 1:mesh_cols
    for yy = 1:mesh_rows
        
           nodeNo(yy, xx) = ind;
           ind = ind + 1;
        
    end
end



nodes = zeros(num_elems, 4);

%start in upper left corner go around counter-clockwise
for ex = 1:el_cols
    for ey = 1:el_rows

        e = ey + (ex-1)*el_rows;
       
        %start in upper left corner go around counter-clockwise
        nodes(e,1) = nodeNo( ey + 1, ex  );
        nodes(e,2) = nodeNo( ey+1, ex+1 );
        nodes(e,3) = nodeNo(  ey, ex + 1 );
        nodes(e,4) = nodeNo(  ey, ex );
        
    
    end
end


el_xDOF = (nodes-1)*2 + 1;
el_yDOF = nodes*2;
%pre-allocating enormous global stiffness matrix
K = zeros(2*num_nodes, 2*num_nodes);


load lagrange_ke



%nodes are counted counter-clockwise:
%
%       4         3
%
%
%       1         2
%now to create the global stiffness matrix
for e = 1:num_elems
   
    
    tempKe = func(e)*ke;
    
    for i = 1:4
        for j = 1:4
          %function values
          K(el_yDOF(e,i), el_yDOF(e,j)  ) = K( el_yDOF(e,i), el_yDOF(e,j)  ) + tempKe(2*i,2*j );
          K(el_xDOF(e,i), el_xDOF(e,j)  ) = K( el_xDOF(e,i), el_xDOF(e,j)  ) + tempKe(1+2*(i-1),1 + 2*(j-1) );
          K(el_yDOF(e,i), el_xDOF(e,j)  ) = K( el_yDOF(e,i), el_xDOF(e,j)  ) + tempKe(2*i,1 +2*(j-1) );
          K(el_xDOF(e,i), el_yDOF(e,j)  ) = K( el_xDOF(e,i), el_yDOF(e,j)  ) + tempKe(1+2*(i-1),2*j );
        end
    end
   
   
end





%now adjust the global stiffness matrix and force vector to enforce 
%bndry conditions, change rows to 0 save diagonal entry, which is one.
%let force vector be equal to dp.  Creates equation Q = dp_m

yDOF = 2*nodeNo;   %for implementing boundary conditions
xDOF = 1 + 2*(nodeNo-1);

bc = [yDOF(1,:), yDOF(end,:), yDOF(2:end-1,1)', yDOF(2:end-1,end)', xDOF(end,end) ];
bcDp = [dp_m(1,:), dp_m(end,:), dp_m(2:end-1,1)', dp_m(2:end-1,end)', 0 ];
   


K(bc, :) = 0;
f = zeros(num_nodes, 1);

for t = 1:length(bc )
    DOF = bc(t);
    K(DOF, DOF) = 1;
    f(DOF) = bcDp(t);            
end


%%%remove every degree of freedom not associated with axial displacements

K = sparse(K);

%now solve
%for the predicted displacement field
dp_p = K\f;

%%%%%%%%


%%%%%now calculate gradient and Hessian simultaneously
%%%%%Hessian will be approximated by [J]^T[J], where J is Jacobian
delta_u = zeros(2*length(dp_m(:)),1);
temp =   reshape(dp_p(2:2:end), mesh_rows, mesh_cols) - dp_m;    %dp_m is possibly a matrix.  Flatten it.
delta_u(2:2:end) = temp;
w = -K'\delta_u(:);


grad = zeros(num_elems,1);


dK_dE = sparse(2*num_nodes, 2*num_nodes, 60);
jacobian = zeros(2*num_nodes, num_elems);

for el = 1:num_elems


     for i = 1:4
        for j = 1:4
            
           dK_dE(el_yDOF(el,i), el_yDOF(el,j)  ) = dK_dE( el_yDOF(el,i), el_yDOF(el,j)  ) + ke(2*i,2*j );
           dK_dE(el_xDOF(el,i), el_xDOF(el,j)  ) = dK_dE( el_xDOF(el,i), el_xDOF(el,j)  ) + ke(1+2*(i-1),1 + 2*(j-1) );
           dK_dE(el_yDOF(el,i), el_xDOF(el,j)  ) = dK_dE( el_yDOF(el,i), el_xDOF(el,j)  ) + ke(2*i,1 +2*(j-1) );
           dK_dE(el_xDOF(el,i), el_yDOF(el,j)  ) = dK_dE( el_xDOF(el,i), el_yDOF(el,j)  ) + ke(1+2*(i-1),2*j );
        end
     end
     
     
    dK_dE(bc, :) = 0;

    %%%%first the gradient
    grad(el) = w'*(dK_dE*dp_p );
   
    %%%%now one column of Jacobian matrix
     jacobian(:, el) = K\(dK_dE*dp_p);
   
    dK_dE(:) = 0;

    
end



[norm, gradNorm, hessNorm] = calculate_L1_regularization(func);
hess = jacobian'*jacobian + alpha*hessNorm;


hess = sparse(hess);
%%%%%%%%%
grad = grad(:) + alpha*gradNorm(:);

%%%%%%%%%%%%%%5
resid = delta_u.^2;
resid = .5*sum(resid(:));
resid = resid + alpha*norm;

