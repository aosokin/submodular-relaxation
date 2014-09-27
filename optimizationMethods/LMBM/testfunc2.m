
function [f, g] = testfunc2(x)

  f = 0.0;
  g(1) = 0.0;
  n = length(x);

  for i=1:(n-1)
    g(i+1) = 0.0;
    a = -x(i)-x(i+1);
    b = -x(i)-x(i+1)+(x(i)*x(i)+x(i+1)*x(i+1)-1.0);

    if (a >= b)
      f = f+a;
      g(i) = g(i)-1.0;
      g(i+1) = -1.0;
    else
      f = f+b;  
      g(i) = g(i)-1.0+2.0*x(i);
      g(i+1) = -1.0+2*x(i+1);
    end
  end

end
