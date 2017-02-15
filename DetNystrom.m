%% Function for Fredholm determinant

function deter = DetNystrom(K,z,a,b,m)
[w,y] = quad(m,a,b); % w = quad(m,a,b)
w = sqrt(w); 
[yj,yi] = meshgrid(y,y);
deter = det(eye(m)+z*(w'*w).*K(yi,yj));
end
