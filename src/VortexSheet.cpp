#include "VortexSheet.hpp"

VortexSheet::VortexSheet()
{
	// Don't put anything in here, it is never called
}

VortexSheet::VortexSheet(const XmlHandler &xml)
{
	m_rho = xml.getValueAttribute("constants", "density");
	m_dt = xml.getValueAttribute("time", "dt");

	m_pi = 4.0 * atan(1.0);
	m_rpi2 = 1.0 / (2.0 * m_pi);
}

void VortexSheet::resize(unsigned size)
{
	gamma.set_size(size);
	x.set_size(size);
	z.set_size(size);
	xc.set_size(size);
	zc.set_size(size);
	theta.set_size(size);
	ds.set_size(size);
	enx.set_size(size);
	enz.set_size(size);
	etx.set_size(size);
	etz.set_size(size);
}

void VortexSheet::vortexsheetbc(VortexBlobs &blobs)
{
	double c1, c2, c5, c6, c7, c8, c9;
	Matrix px, py, qx, qy;

	px.set_size(blobs.size(), size());
	qx.copy_size(px);
	py.copy_size(px);
	qy.copy_size(px);

// compute the coefficients
	for (unsigned i = 0; i < blobs.size(); i++) {

		double xi = blobs.m_x[i];
		double zi = blobs.m_z[i];

		for (unsigned j = 0; j < size(); j++) {

			double xj = xc[j];
			double zj = zc[j];
			double thetaj = theta[j];
			double dsj = ds[j];

			c1 = -(xi - xj) * cos(thetaj) - (zi - zj) * sin(thetaj);
			c2 = std::pow(xi - xj, 2.0) + std::pow(zi - zj, 2.0);
			c5 = (xi - xj) * sin(thetaj) - (zi - zj) * cos(thetaj);
			c6 = log(1.0 + dsj * ((dsj + 2 * c1) / c2));
			c7 = atan2(c5 * dsj, c2 + c1 * dsj);
			c8 =
			    (xi - xj) * sin(-2.0 * thetaj) + (zi - zj) * cos(-2.0 * thetaj);
			c9 =
			    (xi - xj) * cos(-2.0 * thetaj) - (zi - zj) * sin(-2.0 * thetaj);

			qx(i, j) = -sin(thetaj) + 0.5 * c8 * c6 / dsj
			           + (c1 * cos(thetaj) + c5 * sin(thetaj)) * c7 / dsj;
			px(i, j) = -0.5 * c6 * sin(thetaj) - c7 * cos(thetaj) - qx(i, j);

			qy(i, j) = cos(thetaj) + 0.5 * c9 * c6 / dsj
			           + (c1 * sin(thetaj) - c5 * cos(thetaj)) * c7 / dsj;
			py(i, j) = 0.5 * c6 * cos(thetaj) - c7 * sin(thetaj) - qy(i, j);
		}
	}

	// calculate the vortexsheet induced velocities
	unsigned last = size() - 1;

#pragma omp parallel for
	for (unsigned i = 0; i < blobs.size(); i++) {

		blobs.m_uvs[i] = 0.0;
		blobs.m_wvs[i] = 0.0;

		for (unsigned j = 0; j < size(); j++) {
			if (j == last) {
				blobs.m_uvs[i] += px(i, last) * gamma[last]
				                   + qx(i, last) * gamma[0];
				blobs.m_wvs[i] += py(i, last) * gamma[last]
				                   + qy(i, last) * gamma[0];
			} else {
				blobs.m_uvs[i] += px(i, j) * gamma[j]
				                   + qx(i, j) * gamma[j + 1];
				blobs.m_wvs[i] += py(i, j) * gamma[j]
				                   + qy(i, j) * gamma[j + 1];
			}
		}
		blobs.m_uvs[i] *= m_rpi2;
		blobs.m_wvs[i] *= m_rpi2;
	}
}

void VortexSheet::compute_loads()
{
	Vector p = -m_rho / m_dt * (gamma % ds);

	p(0) = 0;

	m_fx = arma::sum(p % enx % ds);
	m_fz = arma::sum(p % enz % ds);

	std::cout << "Fx = " << m_fx << " Fz = " << m_fz << std::endl;
}

unsigned VortexSheet::size()
{
	if ((gamma.size() != xc.size())
	    && (gamma.size() != zc.size())
	    && (gamma.size() != ds.size())
	    && (gamma.size() != theta.size())
	    && (gamma.size() != enx.size())
	    && (gamma.size() != enz.size())
	    && (gamma.size() != etx.size())
	    && (gamma.size() != etz.size())) {
		throw std::string("Size mismatch in VortexSheet");
	}
	return gamma.size();
}

void VortexSheet::print_collocation()
{
	for (unsigned i = 0; i < size(); i++) {
		std::cout << "(xc,zc) = (" << xc(i) << "," << zc(i) << ")" << std::endl;
	}
}

void VortexSheet::print_unit_vectors()
{
	for (unsigned i = 0; i < size(); i++) {
		std::cout << "For panel No " << i << "the normal vectors are en = ("
		          << enx(i) << "," << enz(i) << ") and et = (" << etx(i) << ","
		          << etz(i) << ")" << std::endl;
	}
}

void VortexSheet::print_gamma()
{
	gamma.print();
}

std::tuple<double, double> VortexSheet::get_forces()
{
	return std::make_tuple(m_fx, m_fz);
}
