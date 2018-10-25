function [u,uz] = LogarithmicVelocityProfile(us,kappa,z0,z)

u = us/kappa*log((z+z0)/z0);
uz = (us/kappa)*1./(z+z0);
