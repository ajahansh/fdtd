function return_energy=pml_optimize(k,m,f)
return_args=thin_wall(0,k,m,f,1);
return_energy=return_args.E;