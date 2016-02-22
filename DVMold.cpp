

/*
 
 ////Delete vortices that enter the body
 void DVMBase::absorb()
 {
 int size=m_vortex.size();
 std::vector<unsigned> closest_panel;
 std::vector<double> min_dist;
 closest_panel.resize(m_vortex.size());
 min_dist.resize(m_vortex.size());
 std::vector<int> inside_index;
 inside_index.resize(m_vortex.size());
 
 double x_init,z_init,x_0,z_0,x_1,z_1;
 unsigned count;
 double dx, dz, dr, min_prev;
 count=0;
 
 for(unsigned i=0;i<m_vortex.size();i++)
 {
 min_dist[i]=0;
 closest_panel[i]=0;
 inside_index[i]=inside_body(m_vortex.x[i],m_vortex.z[i]);
 if(inside_index[i]==1)
 {
 min_prev=10E6;
 // Find which panel is closer to the vortex
 for(unsigned j=0;j<m_vortsheet.size();j++)
 {
 dx=m_vortsheet.xc[j]-m_vortex.x[i];
 dz=m_vortsheet.zc[j]-m_vortex.z[i];
 dr=std::sqrt(std::pow(dx,2)+std::pow(dz,2));
 if(dr<min_prev)
 {
 closest_panel[i]=j;
 min_prev=dr;
 }
 }
 // Assign the vorticity to the new closes panel;
 m_Gamma_abs[closest_panel[i]]=m_Gamma_abs[closest_panel[i]]+m_vortex.circ[i];
 
 //Now we can delete the vortex blob
 
 m_vortex.ID[i]=m_vortex.ID.back();m_vortex.ID.pop_back();
 m_vortex.x[i]=m_vortex.x.back();m_vortex.x.pop_back();
 m_vortex.z[i]=m_vortex.z.back();m_vortex.z.pop_back();
 m_GammaDel= m_GammaDel + m_vortex.circ[i];
 m_vortex.circ[i]=m_vortex.circ.back();m_vortex.circ.pop_back();
 m_vortex.sigma[i]=m_vortex.sigma.back();m_vortex.sigma.pop_back();
 m_vortex.u[i]=m_vortex.u.back();m_vortex.u.pop_back();
 m_vortex.w[i]=m_vortex.w.back();m_vortex.w.pop_back();
 m_vortex.omega[i]=m_vortex.omega.back();m_vortex.omega.pop_back();
 m_vortex.uvs[i]=m_vortex.uvs.back();m_vortex.uvs.pop_back();
 m_vortex.wvs[i]=m_vortex.wvs.back();m_vortex.wvs.pop_back();
 count=count+1;
 
 if(size!= (int)m_vortex.size())
 {
 --i;
 size=m_vortex.size();
 }
 
 }
 }
 
 std::cout<<"Number of vortex blob that have crossed the body :  "<<count<<std::endl;
 }
 
 */

/*
 
 ////Delete vortices that enter the body
 void DVMBase::delete_vort()
 {
 std::cout<<"Checking whether a vortex blob has crossed the body or not ..."<<std::endl;
 std::vector<int> inside_index;
 inside_index.resize(m_vortex.size());
 int size=m_vortex.size();
 int count=0;
 m_GammaDel=0;
 for(unsigned i=0;i<size;i++)
 {
 inside_index[i]=inside_body(m_vortex.x[i],m_vortex.z[i]);
 //if 0 is out of the domain and if 1 it is inside
 if(inside_index[i]==1)
 {
 m_vortex.ID[i]=m_vortex.ID.back();m_vortex.ID.pop_back();
 m_vortex.x[i]=m_vortex.x.back();m_vortex.x.pop_back();
 m_vortex.z[i]=m_vortex.z.back();m_vortex.z.pop_back();
 m_GammaDel= m_GammaDel + m_vortex.circ[i];
 m_vortex.circ[i]=m_vortex.circ.back();m_vortex.circ.pop_back();
 m_vortex.sigma[i]=m_vortex.sigma.back();m_vortex.sigma.pop_back();
 m_vortex.u[i]=m_vortex.u.back();m_vortex.u.pop_back();
 m_vortex.w[i]=m_vortex.w.back();m_vortex.w.pop_back();
 m_vortex.omega[i]=m_vortex.omega.back();m_vortex.omega.pop_back();
 m_vortex.uvs[i]=m_vortex.uvs.back();m_vortex.uvs.pop_back();
 m_vortex.wvs[i]=m_vortex.wvs.back();m_vortex.wvs.pop_back();
 count=count+1;
 
 if(size!= (int)m_vortex.size())
 {
 --i;
 size=m_vortex.size();
 }
 }
 }
 std::cout<<"Number of vortex blob that have crossed the body :  "<<count<<std::endl;
 std::cout<<"Total circulation absorbed by the body is : "<<m_GammaDel<<std::endl;
 }
 
 */

/*
 
 void DVMBase::imagebc()
 {
 double pi=4.0*atan(1.0);
 double rsigmasqr;
 double rpi2=1.0/(2.0*pi);
 double dx_ij,dz_ij,dr_ij2,threshold,xkernel;
 double dK_ij, zkernel;
 double x_i,z_i,u_i,w_i,c_i;
 
 std::vector<double> ximage,zimage,circI;
 uI.resize(m_vortex.size());
 wI.resize(m_vortex.size());
 ximage.resize(m_vortex.size());
 zimage.resize(m_vortex.size());
 circI.resize(m_vortex.size());
 
 
 
 for(unsigned i=0;i<m_vortex.size();i++)
 {
 ximage[i]=m_vortex.x[i];
 zimage[i]=-m_vortex.z[i];
 circI[i]=-m_vortex.circ[i];
 uI[i]=0.0;
 wI[i]=0.0;
 }
 
 for(unsigned i=0;i<m_vortex.size();i++)
 {
 x_i=m_vortex.x[i];
 z_i=m_vortex.z[i];
 u_i=0.0;
 w_i=0.0;
 c_i=m_vortex.circ[i];
 
 for(unsigned j=1;j<m_vortex.size();j++)
 {
 dx_ij = x_i-ximage[j];
 dz_ij = z_i-zimage[j];
 dr_ij2 = std::pow(dx_ij,2)+std::pow(dz_ij,2);
 
 // Not sure what sigma we need to use here.
 threshold = m_kernel_threshold * std::pow(m_vortex.sigma[j],2);
 rsigmasqr = 1.0/std::pow(m_vortex.sigma[j],2);
 
 if(dr_ij2<threshold)
 {
 dK_ij=(1.0-std::exp(-dr_ij2*rsigmasqr))/dr_ij2;
 }else{
 dK_ij= 1.0/dr_ij2;
 }
 
 xkernel = dK_ij*dz_ij;
 zkernel = dK_ij*dx_ij;
 u_i = - xkernel*circI[j];
 w_i = + zkernel*circI[j];
 uI[i] = uI[i] + u_i;
 wI[i] = wI[i] + w_i;
 }
 }
 
 for(unsigned i=0;i<m_vortex.size();i++)
 {
 uI[i]=rpi2*uI[i];
 wI[i]=rpi2*wI[i];
 }
 }
 
 */