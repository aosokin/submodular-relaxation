
function [f, g] = testfunc1(x)
  f = 0;
  for j = 1:1000
    f = f + x(j)*x(j);
    g(j) = 2*x(j);
  end
end
