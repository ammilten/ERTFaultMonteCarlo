function save_rockphysics_archie(fname,rho)

rho.Rt_g = rho.Rt_g + archie(rho.R0_g,rho.S_g);
rho.Rt_s = rho.Rt_s + archie(rho.R0_s,rho.S_s);
rho.Rt_ms = rho.Rt_ms + archie(rho.R0_ms,rho.S_ms);
rho.Rt_m = rho.Rt_m + archie(rho.R0_m,rho.S_m);

fid = fopen(fname,'w');
fprintf(fid, '%f,,,\n', rho.h);
fprintf(fid,'%f, %f, %f, %f\n',rho.Rt_g,rho.Rt_s,rho.Rt_ms,rho.Rt_m);
fprintf(fid,'%f, %f, %f, %f',rho.R0_g,rho.R0_s,rho.R0_ms,rho.R0_m);
fclose(fid);
