function [norm, grad, hess] = calculate_L1_regularization(f)

beta =.1; %Huber switching parameter

[rows, cols] = size(f);


%first work out matrixes that give the gradient when they multiply
%the image, f when f is treated as a vector

Dx = zeros(length(f(:)));
Dy = zeros(length(f(:)));



for I = 1:length(f(:));
    
    [i,j] = ind2sub(size(f), I);    
        
        
        if i ~= rows
            Dy(I, I) = 1;
            Dy(I, I + 1) = -1;
        end
        
        if j ~= cols
            Dx(I, I) = 1;
            Dx(I, I + rows) = -1;
        end
        
        
end

Dx = sparse(Dx);
Dy = sparse(Dy);

xDeriv = Dx*f(:);
yDeriv = Dy*f(:);
innerTerm = xDeriv.^2 + yDeriv.^2;

deriv_vec = zeros(length(f(:)), 1);
deriv_II_vec = zeros(length(f(:)), 1);
norm = 0;

%%%These formulas can be found in the monograph Computational Methods
%%%for inverse problems by Vogel
for I = 1:length(f(:))
     if innerTerm(I)  > beta^2 
        norm = norm +  2*sqrt( innerTerm(I) ) - beta;
        deriv_vec(I) = 1/ sqrt(  innerTerm(I) );
        deriv_II_vec(I) = -(1/2)*1/innerTerm(I)^(3/2);
     else
         norm = norm + ( innerTerm(I) )/beta;
         deriv_vec(I) = 1/beta;
         deriv_II_vec(I) = 0;
     end
end

%%%%%%%%gradient calculation%%%%%%%%%%%%%%
L = (Dx'*diag(deriv_vec)*Dx + Dy'*diag(deriv_vec)*Dy );
grad = L*f(:);

%%%%%%%%%%%Hessian calculation%%%%%%%%%%%%%%%%
L_prime = Dx'*diag( xDeriv.^2 .* deriv_vec ) *Dx +...
          Dx'*diag( yDeriv.*xDeriv .* deriv_vec ) *Dy +...
          Dy'*diag( xDeriv.*yDeriv.* deriv_vec ) *Dx +...
          Dy'*diag( yDeriv.^2 .* deriv_vec ) *Dy;
    
hess = L + L_prime;