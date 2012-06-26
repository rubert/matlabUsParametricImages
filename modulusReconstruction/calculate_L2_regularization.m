function [norm, grad, hess] = calculate_L2_regularization(f)


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


%%%These formulas can be found in the monograph Computational Methods
%%%for inverse problems by Vogel
norm = sum(yDeriv.^2) + sum(xDeriv.^2);
grad = (Dx'*Dx + Dy'*Dy)*f(:);
hess = (Dx'*Dx + Dy'*Dy);