:Begin:
:Function:	aa_to_cart
:Pattern: 	AAToCart[J_List, th_List]
:Arguments: 	{J, th}
:ArgumentTypes: {RealList, RealList}
:ReturnType: 	Manual
:End:
:Evaluate:	AAToCart::usage = "AAToCart[{Lz, L, Jr}, {th1, th2, th3}] takes the three actions and their corresponding angles in the isochrone potential and calculates the corresponding Cartesian position and velocity, returning a six-vector w={x,y,z,vx,vy,vz}. For purely radial orbits the inclination angle is undefined; these orbits are placed in the z=0 plane and a warning is printed. The program returns a zero vector and prints an error message if the orbit is unbound or has Lz>L."