#include "VortexSheet.hpp"

#include <filesystem>

namespace fs = std::filesystem;

VortexSheet::VortexSheet()
{
    // Don't put anything in here, it is never called
}

VortexSheet::VortexSheet(const XmlHandler &xml, const std::string &stamp)
{
    m_rho = xml.getValueAttribute("constants", "density");
    m_nu = xml.getValueAttribute("constants", "nu");
    m_dt = xml.getValueAttribute("time", "dt");

    m_Ux = xml.getValueAttribute("flow", "ux");
    m_Uz = xml.getValueAttribute("flow", "uz");

    m_maxNumPanelVort = static_cast<unsigned>(xml.getValueAttribute("constants", "max_NumPanelVort"));
    m_cutoff_exp = xml.getValueAttribute("constants", "cutoff_exp");

    auto indir = xml.getStringAttribute("io", "input_dir");
    fs::path inpath(indir);

    std::string domain_file = xml.getStringAttribute("io", "domain_file");
    fs::path domain = inpath / domain_file;

    read_input_coord(domain.string());
    form_vortex_sheet();
    compute_influence_matrix();

    // Initialise the output file
    auto outdir = xml.getStringAttribute("io", "output_dir", true);
    fs::path outpath(outdir);

    auto gamma_path = outpath / (stamp + "_gamma.dat");
    m_gammafile.open(gamma_path);
    if (!m_gammafile) {
        throw dvm::FileIOException::cannotOpen(gamma_path.string());
    }
    m_gammafile << m_xc.n_elem << " # Number of collocation points\n"
                << m_dt << " # Time Step\n"
                << xml.getValueAttribute("time", "steps") << " # Steps\n"
                << "Time [s] Gamma - Vortex sheet strength\n";

    auto force_path = outpath / (stamp + "_loads.dat");
    m_forcefile.open(force_path);
    if (!m_forcefile) {
        throw dvm::FileIOException::cannotOpen(force_path.string());
    }
    m_forcefile << "Time [s]\tF_x [-]\tF_z [s]\n";
}

void VortexSheet::read_input_coord(const std::string& file)
{
    // Make use of the armadillo functionality to read the file
    Matrix coords;
    bool status = coords.load(file, arma::raw_ascii);

    if (!status) {
        throw dvm::FileIOException::cannotRead(file);
    }

    m_x = coords.col(0);
    m_z = coords.col(1);

    std::cout << "Successfully loaded coordinate file with " << m_x.n_elem
              << " points." << std::endl;

    // Assume for the moment that we are looking at a cylinder - get the
    // diameter for the drag coefficient reference area.
    m_A = arma::max(m_x) - arma::min(m_x);
    std::cout << "A = " << m_A << " (assuming cylinder)\n";
}

void VortexSheet::form_vortex_sheet()
{
    // Create the collocation points
    const auto end = m_x.n_elem - 1;
    m_xc = 0.5 * (m_x.rows(0, end - 1) + m_x.rows(1, end));
    m_zc = 0.5 * (m_z.rows(0, end - 1) + m_z.rows(1, end));

    // ds and theta along the vortex sheet
    Vector dx = m_x.rows(1, end) - m_x.rows(0, end - 1);
    Vector dz = m_z.rows(1, end) - m_z.rows(0, end - 1);

    m_ds = arma::sqrt(arma::pow(dx, 2) + arma::pow(dz, 2));
    m_theta = arma::atan2(dz, dx);

    // Outwards facing normals and tangentials
    m_enx = dz / m_ds;
    m_enz = -dx / m_ds;
    m_etx = -m_enz;
    m_etz = m_enx;

    // Make sure the surface vorticity is the correct size
    m_gamma.copy_size(m_xc);

    std::cout << "Created vortex sheet of size " << m_gamma.n_elem << std::endl;
}

void VortexSheet::compute_influence_matrix()
{
    // Compute influence matrix according to coefficients after Morgenthal
    // Follow Morgenthal (2002)
    const unsigned Nl = static_cast<unsigned>(m_gamma.n_elem);

    std::cout << "Nl = " << Nl << "\n";

    m_infM.set_size(Nl + 1, Nl);

    Matrix p(Nl, Nl);
    Matrix q(Nl, Nl);

    for (unsigned i = 0; i < Nl; i++) {
        const double xci = m_xc(i);
        const double zci = m_zc(i);
        const double thetai = m_theta(i);

        for (unsigned j = 0; j < Nl; j++) {
            if (i == j) {
                p(i, j) = -1.0;
                q(i, j) = 1.0;
            } else {
                const double xj = m_x(j);
                const double zj = m_z(j);
                const double thetaj = m_theta(j);
                const double dsj = m_ds(j);

                const double c1 = -(xci - xj) * std::cos(thetaj)
                                  - (zci - zj) * std::sin(thetaj);
                const double c2 = std::pow(xci - xj, 2) + std::pow(zci - zj, 2);
                const double c3 = std::sin(thetai - thetaj);
                const double c4 = std::cos(thetai - thetaj);
                const double c5 = (xci - xj) * std::sin(thetaj)
                                  - (zci - zj) * std::cos(thetaj);
                const double c6 = std::log(1.0 + dsj * ((dsj + 2 * c1) / c2));
                const double c7 = std::atan2((c5 * dsj), (c2 + c1 * dsj));
                const double c9 = (xci - xj) * std::cos(thetai - 2.0 * thetaj)
                                  - (zci - zj) * std::sin(thetai - 2.0 * thetaj);

                q(i, j) = c4 + 0.5 * c9 * c6 / dsj - ((c1 * c3 + c4 * c5) * c7) / dsj;
                p(i, j) = 0.5 * c4 * c6 + c3 * c7 - q(i, j);
            }
        }
    }

    for (unsigned i = 0; i < Nl; i++) {
        m_infM(i, 0) = -math::rpi2 * (p(i, 0) + q(i, Nl - 1));

        for (unsigned j = 1; j < Nl; j++) {
            m_infM(i, j) = -math::rpi2 * (p(i, j) + q(i, j - 1));
        }
    }

    // Enforcing the total circulation
    for (unsigned j = 0; j < Nl; j++) {
        m_infM(Nl, j) = m_ds(j);
    }

    std::cout << "Computed influence matrix of size " << arma::size(m_infM)
              << " after Kuette and Chow" << std::endl;
}

void VortexSheet::solvevortexsheet(VortexBlobs &blobs)
{
    const unsigned Nl = static_cast<unsigned>(m_gamma.n_elem);

    Vector brhs(Nl + 1);

    if (blobs.size() == 0) {
        brhs.rows(0, Nl - 1) = m_Ux * m_enx + m_Uz * m_enz;
        brhs(Nl) = -blobs.totalcirc();

    } else {
        Vector u(Nl);
        Vector w(Nl);

        #pragma omp parallel for
        for (unsigned i = 0; i < Nl; i++) {
            u(i) = 0.0;
            w(i) = 0.0;

            for (unsigned j = 0; j < blobs.size(); j++) {
                const double dx_ij = m_xc(i) - blobs.m_x(j);
                const double dz_ij = m_zc(i) - blobs.m_z(j);
                const double dr_ij2 = std::pow(dx_ij, 2) + std::pow(dz_ij, 2);

                const double threshold = 10 * std::pow(blobs.m_sigma(j), 2);
                const double rsigmasqr = 1.0 / std::pow(blobs.m_sigma(j), 2);

                double dK_ij;
                if (dr_ij2 < threshold) {
                    dK_ij = (1.0 - std::exp(-dr_ij2 * rsigmasqr)) / dr_ij2;
                } else {
                    dK_ij = 1.0 / dr_ij2;
                }

                u(i) -= dK_ij * dz_ij * blobs.m_circ(j);
                w(i) += dK_ij * dx_ij * blobs.m_circ(j);
            }
        }

        u *= math::rpi2;
        w *= math::rpi2;

        brhs.rows(0, Nl - 1) = (m_Ux + u) % m_enx + (m_Uz + w) % m_enz;
        brhs(Nl) = -blobs.totalcirc();
    }

    // Solve system
    m_gamma = arma::solve(m_infM.t() * m_infM, m_infM.t() * brhs);
}

void VortexSheet::vortexsheetbc(VortexBlobs &blobs)
{
    Matrix px(blobs.size(), size());
    Matrix qx(blobs.size(), size());
    Matrix py(blobs.size(), size());
    Matrix qy(blobs.size(), size());

    // compute the coefficients
    for (unsigned i = 0; i < blobs.size(); i++) {
        const double xi = blobs.m_x(i);
        const double zi = blobs.m_z(i);

        for (unsigned j = 0; j < size(); j++) {
            const double xj = m_xc(j);
            const double zj = m_zc(j);
            const double thetaj = m_theta(j);
            const double dsj = m_ds(j);

            const double c1 = -(xi - xj) * std::cos(thetaj) - (zi - zj) * std::sin(thetaj);
            const double c2 = std::pow(xi - xj, 2.0) + std::pow(zi - zj, 2.0);
            const double c5 = (xi - xj) * std::sin(thetaj) - (zi - zj) * std::cos(thetaj);
            const double c6 = std::log(1.0 + dsj * ((dsj + 2 * c1) / c2));
            const double c7 = std::atan2(c5 * dsj, c2 + c1 * dsj);
            const double c8 = (xi - xj) * std::sin(-2.0 * thetaj) + (zi - zj) * std::cos(-2.0 * thetaj);
            const double c9 = (xi - xj) * std::cos(-2.0 * thetaj) - (zi - zj) * std::sin(-2.0 * thetaj);

            qx(i, j) = -std::sin(thetaj) + 0.5 * c8 * c6 / dsj
                       + (c1 * std::cos(thetaj) + c5 * std::sin(thetaj)) * c7 / dsj;
            px(i, j) = -0.5 * c6 * std::sin(thetaj) - c7 * std::cos(thetaj) - qx(i, j);

            qy(i, j) = std::cos(thetaj) + 0.5 * c9 * c6 / dsj
                       + (c1 * std::sin(thetaj) - c5 * std::cos(thetaj)) * c7 / dsj;
            py(i, j) = 0.5 * c6 * std::cos(thetaj) - c7 * std::sin(thetaj) - qy(i, j);
        }
    }

    // calculate the vortexsheet induced velocities
    const unsigned last = size() - 1;

    #pragma omp parallel for
    for (unsigned i = 0; i < blobs.size(); i++) {
        blobs.m_uvs(i) = 0.0;
        blobs.m_wvs(i) = 0.0;

        for (unsigned j = 0; j < size(); j++) {
            if (j == last) {
                blobs.m_uvs(i) += px(i, last) * m_gamma(last) + qx(i, last) * m_gamma(0);
                blobs.m_wvs(i) += py(i, last) * m_gamma(last) + qy(i, last) * m_gamma(0);
            } else {
                blobs.m_uvs(i) += px(i, j) * m_gamma(j) + qx(i, j) * m_gamma(j + 1);
                blobs.m_wvs(i) += py(i, j) * m_gamma(j) + qy(i, j) * m_gamma(j + 1);
            }
        }
        blobs.m_uvs(i) *= math::rpi2;
        blobs.m_wvs(i) *= math::rpi2;
    }
}

VortexBlobs VortexSheet::release_nascent_vortices_rw(Random &_rand)
{
    const auto Nvs = size();

    // First we need to determine how many of these vortices will be created at
    // each panel. This is in accordance with Morgenthal PhD 4.3.6 Vortex release
    // algorithm (eq. 4.108)
    Vector PanelCirc = m_gamma % m_ds; // Panel's total circulation

    // Absolute value of the panel vorticity
    Vector AbsPanelCirc = arma::abs(PanelCirc / arma::max(PanelCirc) * m_maxNumPanelVort);

    Vector PanelNewVort = arma::floor(AbsPanelCirc) + arma::ones(Nvs);
    // Make sure the new panel vortices are non-negative
    if (!arma::all(PanelNewVort > 1.0e-15)) {
        throw dvm::VortexException("VortexSheet::release_nascent_vortices_rw: negative PanelNewVort");
    }

    Vector Circ_new = PanelCirc / PanelNewVort;

    // Number of the new vortices.
    // if we set the maximum number of vortices from each panel=1, then Nrv=Nsv
    const unsigned Nrv = static_cast<unsigned>(arma::sum(PanelNewVort));

    // Initialize the vortexblobs ready for release
    VortexBlobs nascentVort(Nrv);

    // Calculating the position of the newly released vortices
    unsigned counter = 0;

    for (unsigned i = 0; i < Nvs; i++) {
        unsigned j = 0;
        while (j < static_cast<unsigned>(PanelNewVort(i))) {
            // Find locations at the ith panel from which the vortices should be released
            const double xm = m_x(i) + m_ds(i) * m_etx(i) / (PanelNewVort(i) + 1);
            const double zm = m_z(i) + m_ds(i) * m_etz(i) / (PanelNewVort(i) + 1);

            // Create Random walk values
            const double R = _rand.rand();

            // This is half normal distribution only in the outwards direction
            const double rrw = std::abs(std::sqrt(4.0 * m_nu * m_dt * std::log(1.0 / R)));

            // Add the released vortex
            nascentVort.m_x(counter) = xm + rrw * m_enx(i);
            nascentVort.m_z(counter) = zm + rrw * m_enz(i);
            nascentVort.m_circ(counter) = Circ_new(i);

            // Now for the cut-off kernel we implement things as suggested by
            // Mirta Perlman 1985 (JCP) using sigma=ds^q where 0.5<q<1
            nascentVort.m_sigma(counter) = std::pow(m_ds(i), m_cutoff_exp);

            counter++;
            assert(counter <= Nrv);
            j++;
        }
    }
    return nascentVort;
}

void VortexSheet::write_step(double time)
{
    for (unsigned i = 0; i < m_xc.n_elem; ++i) {
        m_gammafile << time << "\t " << m_xc(i) << "\t" << m_zc(i) << "\t"
                    << m_gamma(i) << "\t" << m_ds(i) << "\n";
    }

    m_forcefile << time << "\t" << m_fx << "\t" << m_fz << "\n";
}

void VortexSheet::reflect(VortexBlobs &vortex)
{
    Vector_un closest_panel(vortex.size());
    Vector min_dist(vortex.size());

    for (unsigned i = 0; i < vortex.size(); i++) {
        min_dist(i) = 0;
        closest_panel(i) = 0;

        if (inside_body(vortex.m_x(i), vortex.m_z(i))) {
            // Find which panel is closest to the vortex
            double min_prev = 10E6;
            for (unsigned j = 0; j < size(); j++) {
                const double dx = m_xc(j) - vortex.m_x(i);
                const double dz = m_zc(j) - vortex.m_z(i);
                const double dr = std::sqrt(std::pow(dx, 2.0) + std::pow(dz, 2.0));

                if (dr < min_prev) {
                    closest_panel(i) = j;
                    min_prev = dr;
                }
            }

            // Find the mirror image vortex blob
            const double x_init = vortex.m_x(i);
            const double z_init = vortex.m_z(i);

            const double x_0 = m_x(closest_panel(i));
            const double z_0 = m_z(closest_panel(i));
            const double x_1 = m_x(closest_panel(i) + 1);
            const double z_1 = m_z(closest_panel(i) + 1);

            Vector _mirror = mirror(x_init, z_init, x_0, z_0, x_1, z_1);

            vortex.m_x(i) = _mirror(0);
            vortex.m_z(i) = _mirror(1);
            // All other properties of the vortices remain the same
        }
    }
}

bool VortexSheet::inside_body(double xcoor, double zcoor) const
{
    // Using this algorithm, cn == 1 when the point is inside the boundary
    int cn = 0;

    for (unsigned i = 0; i < size(); i++) {
        if (((m_z(i) <= zcoor) && (m_z(i + 1) > zcoor))
            || ((m_z(i) > zcoor) && (m_z(i + 1) <= zcoor))) {

            const auto vt = (zcoor - m_z(i)) / (m_z(i + 1) - m_z(i));
            if (xcoor < m_x(i) + vt * (m_x(i + 1) - m_x(i))) {
                ++cn;
            }
        }
    }

    return (cn == 1);
}

Vector VortexSheet::mirror(double x_init,
                           double z_init,
                           double x_0,
                           double z_0,
                           double x_1,
                           double z_1) const
{
    Vector p2(2);

    const double dx = x_1 - x_0;
    const double dz = z_1 - z_0;

    const double a = (dx * dx - dz * dz) / (dx * dx + dz * dz);
    const double b = 2.0 * dx * dz / (dx * dx + dz * dz);

    p2(0) = a * (x_init - x_0) + b * (z_init - z_0) + x_0;
    p2(1) = b * (x_init - x_0) - a * (z_init - z_0) + z_0;

    return p2;
}

void VortexSheet::compute_loads(double Ur)
{
    Vector P(size(), arma::fill::zeros);
    Vector Dp = -m_rho / m_dt * (m_gamma % m_ds);

    for (unsigned i = 0; i < size() - 1; i++) {
        P(i + 1) = P(i) + Dp(i);
    }

    // Finding the average of the two pressures at the boundary
    for (unsigned i = 0; i < size() - 1; i++) {
        P(i) = 0.5 * (P(i) + P(i + 1));
    }

    P(size() - 1) = 0.5 * (P(size() - 1) + P(0));
    const double pmax = arma::max(P);
    const double pref = 0.5 * m_rho * Ur * Ur;

    P += (pref - pmax) * arma::ones(size(), 1);

    // To find the forces we need to measure on the particle
    // we need to change sign
    m_fx = -arma::sum(P % m_enx % m_ds);
    m_fz = -arma::sum(P % m_enz % m_ds);

    std::cout << "C_D = " << m_fx / (0.5 * m_rho * m_A * Ur * Ur) << "\t"
              << "C_L = " << m_fz / (0.5 * m_rho * m_A * Ur * Ur) << std::endl;
}

unsigned VortexSheet::size() const
{
    if ((m_gamma.size() != m_xc.size()) || (m_gamma.size() != m_zc.size())
        || (m_gamma.size() != m_x.size() - 1)
        || (m_gamma.size() != m_z.size() - 1)
        || (m_gamma.size() != m_ds.size())
        || (m_gamma.size() != m_theta.size())
        || (m_gamma.size() != m_enx.size())
        || (m_gamma.size() != m_enz.size())
        || (m_gamma.size() != m_etx.size())
        || (m_gamma.size() != m_etz.size())) {
        throw dvm::SizeMismatchException("VortexSheet");
    }
    return static_cast<unsigned>(m_gamma.size());
}

void VortexSheet::print_collocation() const
{
    for (unsigned i = 0; i < size(); i++) {
        std::cout << "(xc,zc) = (" << m_xc(i) << "," << m_zc(i) << ")"
                  << std::endl;
    }
}

void VortexSheet::print_unit_vectors() const
{
    for (unsigned i = 0; i < size(); i++) {
        std::cout << "For panel No " << i << " the normal vectors are en = ("
                  << m_enx(i) << "," << m_enz(i) << ") and et = (" << m_etx(i)
                  << "," << m_etz(i) << ")" << std::endl;
    }
}

void VortexSheet::print_gamma() const
{
    m_gamma.print();
}

std::tuple<double, double> VortexSheet::get_forces() const
{
    return std::make_tuple(m_fx, m_fz);
}

void VortexSheet::probe_velocities(const VortexBlobs &blobs, Probe &probe)
{
    // Compute the velocity vector at the probe points
    for (unsigned i = 0; i < probe.size(); i++) {
        probe.m_u(i) = 0.0;
        probe.m_w(i) = 0.0;
    }

    for (unsigned i = 0; i < probe.size(); i++) {
        const double x_i = probe.m_x(i);
        const double z_i = probe.m_z(i);
        double u_i = 0.0;
        double w_i = 0.0;

        for (unsigned j = 1; j < blobs.size(); j++) {
            const double dx_ij = x_i - blobs.m_x(j);
            const double dz_ij = z_i - blobs.m_z(j);
            const double dr_ij2 = std::pow(dx_ij, 2) + std::pow(dz_ij, 2);

            const double threshold = 10 * std::pow(blobs.m_sigma(j), 2);
            const double rsigmasqr = 1.0 / std::pow(blobs.m_sigma(j), 2);

            double dK_ij;
            if (dr_ij2 < threshold) {
                dK_ij = (1.0 - std::exp(-dr_ij2 * rsigmasqr)) / dr_ij2;
            } else {
                dK_ij = 1.0 / dr_ij2;
            }

            const double xkernel = dK_ij * dz_ij;
            const double zkernel = dK_ij * dx_ij;
            u_i = -xkernel * blobs.m_circ(j);
            w_i = +zkernel * blobs.m_circ(j);
            probe.m_u(i) = probe.m_u(i) + u_i;
            probe.m_w(i) = probe.m_w(i) + w_i;
        }
    }

    Matrix px(probe.size(), size());
    Matrix qx(probe.size(), size());
    Matrix py(probe.size(), size());
    Matrix qy(probe.size(), size());

    for (unsigned i = 0; i < probe.size(); i++) {
        probe.m_uvs(i) = 0.0;
        probe.m_wvs(i) = 0.0;
    }

    for (unsigned i = 0; i < probe.size(); i++) {
        for (unsigned j = 1; j < size(); j++) {
            const double c1 = -(probe.m_x(i) - m_xc(j)) * std::cos(m_theta(j))
                              - (probe.m_z(i) - m_zc(j)) * std::sin(m_theta(j));
            const double c2 = std::pow((probe.m_x(i) - m_xc(j)), 2)
                              + std::pow((probe.m_z(i) - m_zc(j)), 2);
            const double c5 = (probe.m_x(i) - m_xc(j)) * std::sin(m_theta(j))
                              - (probe.m_z(i) - m_zc(j)) * std::cos(m_theta(j));
            const double c6 = std::log(1.0 + m_ds(j) * ((m_ds(j) + 2 * c1) / c2));
            const double c7 = std::atan2((c5 * m_ds(j)), (c2 + c1 * m_ds(j)));
            const double c8 = (probe.m_x(i) - m_xc(j)) * std::sin(-2.0 * m_theta(j))
                              + (probe.m_z(i) - m_zc(j)) * std::cos(-2.0 * m_theta(j));
            const double c9 = (probe.m_x(i) - m_xc(j)) * std::cos(-2.0 * m_theta(j))
                              + (probe.m_z(i) - m_zc(j)) * std::sin(-2.0 * m_theta(j));

            qx(i, j) = -std::sin(m_theta(j)) + 0.5 * c8 * c6 / m_ds(j)
                       + (c1 * std::cos(m_theta(j)) + c5 * std::sin(m_theta(j))) * c7 / m_ds(j);
            px(i, j) = -0.5 * c6 * std::sin(m_theta(j)) - c7 * std::cos(m_theta(j)) - qx(i, j);

            qy(i, j) = std::cos(m_theta(j)) + 0.5 * c9 * c6 / m_ds(j)
                       + (c1 * std::cos(m_theta(j)) - c5 * std::sin(m_theta(j))) * c7 / m_ds(j);
            py(i, j) = -0.5 * c6 * std::cos(m_theta(j)) - c7 * std::sin(m_theta(j)) - qy(i, j);

            if (j == size() - 1) {
                probe.m_uvs(i) = (px(i, size() - 1) * m_gamma(size() - 1)
                                  + qx(i, size() - 1) * m_gamma(0));
                probe.m_wvs(i) = (py(i, size() - 1) * m_gamma(size() - 1)
                                  - qy(i, size() - 1) * m_gamma(0));
            } else {
                probe.m_uvs(i) = (px(i, j) * m_gamma(j) + qx(i, j) * m_gamma(j + 1));
                probe.m_wvs(i) = (py(i, j) * m_gamma(j) + qy(i, j) * m_gamma(j + 1));
            }
        }
    }

    for (unsigned i = 0; i < probe.size(); i++) {
        probe.m_u(i) = math::rpi2 * probe.m_u(i) + probe.m_uvs(i) + m_Ux;
        probe.m_w(i) = math::rpi2 * probe.m_w(i) + probe.m_wvs(i) + m_Uz;
    }
}
