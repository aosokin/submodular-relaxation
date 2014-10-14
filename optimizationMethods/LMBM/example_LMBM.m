function example_LMBM

[x,fval,niter,nfeval,term,time]=lmbm_driver('testfunc1',10*ones(1000,1), 1000, 1);

rpar = [ 0, 0, 0, 1e-5, 0, 0.5, 0, 0 ];
ipar = [ 0, 5000000, 5000000, 0, -1, 0, 0 ];
lmbm_driver('testfunc2', -0.5*ones(1000,1), 1000, 1, 300, 10, 7, 7, rpar, ipar);

end
