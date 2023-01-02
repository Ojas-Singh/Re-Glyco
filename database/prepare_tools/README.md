#Tools to add new structure to database from MD simulation.

# step 1
```
parm structure1.ion.prm7
trajin prod1.nc  
trajin prod2.nc
autoimage
strip :WAT,Cl-,Na+
strip @H* 
trajout prodfull.dry.pdb pdb
```

# step 2

