:Begin:
:Function:	cart_to_aa
:Pattern: 	CartToAA[Miso_Real, biso_Real, x_List, v_List]
:Arguments: 	{Miso, biso, x, v}
:ArgumentTypes: {Real, Real, RealList, RealList}
:ReturnType: 	Manual
:End:
:Evaluate:	CartToAA::usage = "CartToAA[Miso, biso, {x,y,z}, {vx,vz,vz}] returns the three actions {Lz,L,Jr} and the conjugate angles, in the isochrone potential with total mass Miso (in kpc^3/Myr^2) and scale radius biso (in kpc), given the Cartesian phase-space coordinates x and v. If the orbit is unbound the program returns the nonsense values {-100,-100,-100} for the actions."