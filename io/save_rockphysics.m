function save_rockphysics(fname,rho)

fid = fopen(fname,'w');
fprintf(fid, '%f,,,\n', rho.h);
fprintf(fid,'%f, %f, %f, %f\n',rho.g_a,rho.s_a,rho.ms_a,rho.m_a);
fprintf(fid,'%f, %f, %f, %f',rho.g_b,rho.s_b,rho.ms_b,rho.m_b);
fclose(fid);

