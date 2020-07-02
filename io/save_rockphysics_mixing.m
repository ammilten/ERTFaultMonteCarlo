function save_rockphysics_mixing(fname,rho)

rho.Rb_g = mixing(rho.Ra_g,rho.phi_g,rho.Rw);
rho.Rb_s = mixing(rho.Ra_s,rho.phi_s,rho.Rw);
rho.Rb_ms = mixing(rho.Ra_ms,rho.phi_ms,rho.Rw);
rho.Rb_m = mixing(rho.Ra_m,rho.phi_m,rho.Rw);

fid = fopen(fname,'w');
fprintf(fid, '%f,,,\n', rho.h);
fprintf(fid,'%f, %f, %f, %f\n',rho.Ra_g,rho.Ra_s,rho.Ra_ms,rho.Ra_m);
fprintf(fid,'%f, %f, %f, %f',rho.Rb_g,rho.Rb_s,rho.Rb_ms,rho.Rb_m);
fclose(fid);
