# \file blockflow2d.py
# \ingroup mbcns2
# \brief A module to pick up the flow data for mbcns2 block-structured grids.
# \author PAJ
# \version 11-Apr-2007 : first cut after talking about it for so long
"""
Pick up the flow data for mbcns2 block-structured grids.
"""

from copy import copy

try:
    from numpy import array, zeros
except:
    try:
        from Numeric import array, zeros
    except:
        print "Could import neither numpy nor Numeric."


from math import fabs

UFLOWZ = 1.0e-20

def apply_uflowz(val):
    if fabs(val) <= UFLOWZ:
        return 0.0
    else:
        return val
        
class BlockFlow2D(object):
    """
    Storage and service functions for a block of cells within a block-structured grid.
    """
    def __init__(self, ni=None, nj=None, label=None,
                 nsp=1, nvib=0, separate_electron_energy=0,
                 radiation=0, k_omega=0, scalar_pdf=0, debug=0):
        """
        ni: number of cells in the i-index direction
        nj: number of cells in the j-index direction
        label: (optional) string label for the block
        
        If the number of vertices are specified in each direction,
        we actually create the storage arrays now.
        The number of grid vertices in each index-direction
        will be one more than the number of finite-volume cells
        in that direction.
        """
        self.ni = ni
        self.nj = nj
        if debug:
            print "New BlockFlow2D: ni=", ni, "nj=", nj
        if ni != None and nj != None:
            self.init_arrays()
        self.label = label
        self.t = 0.0
        self.nsp = nsp
        self.nvib = nvib
        self.separate_electron_energy = separate_electron_energy
        self.radiation = radiation
        self.k_omega = k_omega
        self.scalar_pdf = scalar_pdf
        return

    def init_arrays(self):
        """
        Create the storage arrays at the previously-specified sizes.

        Note that the mass fractions are stored in an array of lists.
        (Each list contains the species mass fractions for a cell.)
        The same arrangement is used for the vibrational energies
        and temperatures.
        """
        self.x = zeros((self.ni, self.nj), 'd')
        self.y = zeros((self.ni, self.nj), 'd')
        self.rho = zeros((self.ni, self.nj), 'd')
        self.vx = zeros((self.ni, self.nj), 'd')
        self.vy = zeros((self.ni, self.nj), 'd')
        self.vz = zeros((self.ni, self.nj), 'd')
        self.e = zeros((self.ni, self.nj), 'd')
        self.p = zeros((self.ni, self.nj), 'd')
        self.a = zeros((self.ni, self.nj), 'd')
        self.T = zeros((self.ni, self.nj), 'd')
        self.mu = zeros((self.ni, self.nj), 'd')
        self.mu_t = zeros((self.ni, self.nj), 'd')
        self.S = zeros((self.ni, self.nj), 'd')
        self.f = zeros((self.ni, self.nj, self.nsp), 'd')
        if self.nvib > 0:
            self.e_vib = zeros((self.ni, self.nj, self.nvib), 'd')
            self.T_vib = zeros((self.ni, self.nj, self.nvib), 'd')
        if self.radiation:
            self.Q_rE_rad = zeros((self.ni, self.nj), 'd')
        if self.k_omega:
            self.tke = zeros((self.ni, self.nj), 'd')
            self.omega = zeros((self.ni, self.nj), 'd')
        if self.scalar_pdf:
            self.sigma_T = zeros((self.ni, self.nj), 'd')
            self.sigma_c = zeros((self.ni, self.nj), 'd')
        if self.separate_electron_energy:
            self.e_e = zeros((self.ni, self.nj), 'd')
            self.T_e = zeros((self.ni, self.nj), 'd')
            self.p_e = zeros((self.ni, self.nj), 'd')
        return

    def read_solution_for_cell(self, fp, i, j):
        """
        Read the next available flow data and store it in the i,j cell.

        This function must be kept in sync with read_solution_for_cell()
        which is found in lib/fv_core/source/cns_cell.cxx.
        """
        line = fp.readline().strip()
        items = line.split()
        self.x[i][j] = float(items[0])
        self.y[i][j] = float(items[1])
        if len(items) == 3:
            self.volume = float(items[2])
        else:
            self.volume = 1.0
        line = fp.readline().strip()
        items = line.split()
        self.rho[i][j] = float(items[0])
        self.vx[i][j] = float(items[1])
        self.vy[i][j] = float(items[2])
        self.vz[i][j] = float(items[3])
        self.e[i][j] = float(items[4])
        line = fp.readline().strip()
        items = line.split()
        self.p[i][j] = float(items[0])
        self.a[i][j] = float(items[1])
        self.T[i][j] = float(items[2])
        self.mu[i][j] = float(items[3])
        self.mu_t[i][j] = float(items[4])
        self.S[i][j] = float(items[5])
        if self.radiation:
            self.Q_rE_rad[i][j] = float(items[6])
        if self.k_omega:
            line = fp.readline().strip()
            items = line.split()
            self.tke[i][j] = float(items[0])
            self.omega[i][j] = float(items[1])
        if self.scalar_pdf:
            line = fp.readline().strip()
            items = line.split()
            self.sigma_T[i][j] = float(items[0])
            self.sigma_c[i][j] = float(items[1])
        line = fp.readline().strip()
        for isp, item in enumerate(line.split()):
            self.f[i][j][isp] = float(item)
        if self.nvib > 0:
            line = fp.readline().strip()
            items = line.split()
            for ivib in range(len(items)/2):
                self.e_vib[i][j][ivib] = float(items[ivib*2])
                self.T_vib[i][j][ivib] = float(items[ivib*2+1])
        if self.separate_electron_energy:
            line = fp.readline().strip()
            items = line.split()
            self.e_e[i][j] = float(items[0])
            self.T_e[i][j] = float(items[1])
            self.p_e[i][j] = float(items[2])
        return
    
    def read(self, fp, debug=0):
        """
        Read the next available block of flow data.

        The solution file may contain several solution times,
        each with several blocks of data
        """
        line = fp.readline().strip()
        items = line.split()
        self.t = float(items[0])
        line = fp.readline().strip()
        items = line.split()
        self.ni = int(items[0])
        self.nj = int(items[1])
        if debug:
            print "Read block: t=", self.t, "ni=", self.ni, "nj=", self.nj
        self.init_arrays()
        for i in range(self.ni):
            for j in range(self.nj):
                self.read_solution_for_cell(fp, i, j)
        return

    def data_in_dictionary(self):
        data = {}
        data["x"] = copy(self.x)
        data["y"] = copy(self.y)
        data["rho"] = copy(self.rho)
        data["vel_x"] = copy(self.vx)
        data["vel_y"] = copy(self.vy)
        data["vel_z"] = copy(self.vz)
        data["vel"] = (copy(self.vx), copy(self.vy), copy(self.vz))
        data["e"] = copy(self.e)
        data["p"] = copy(self.p)
        data["a"] = copy(self.a)
        data["T"] = copy(self.T)
        data["mu"] = copy(self.mu)
        data["mu_t"] = copy(self.mu_t)
        data["S"] = copy(self.S)
        if self.radiation:
            data["Q_rE_rad"] = copy(self.Q_rE_rad)
        if self.k_omega:
            data["tke"] = copy(self.tke)
            data["omega"] = copy(self.omega)
        if self.scalar_pdf:
            data["sigma_T"] = copy(self.sigma_T)
            data["sigma_c"] = copy(self.sigma_c)
        for isp in range(len(self.f[0][0])):
            data["f-%d" % isp] = copy(self.f[:,:,isp])
        if self.nvib > 0:
            for ivib in range(len(self.T_vib[0][0])):
                data["T_vib-%d" % ivib] = copy(self.T_vib[:,:,ivib])
                data["e_vib-%d" % ivib] = copy(self.e_vib[:,:,ivib])
        if self.separate_electron_energy:
            data["e_e"] = copy(self.e_e)
            data["T_e"] = copy(self.T_e)
            data["p_e"] = copy(self.p_e)

        # Flush small values to zero
        for data_list in data.itervalues():
            if type(data_list) == tuple:
                # We have a vector most likely
                for v_array in data_list:
                    for i in range(len(v_array.flat)):
                        v_array.flat[i] = apply_uflowz(v_array.flat[i])
            else:
                for i in range(len(data_list.flat)):
                    data_list.flat[i] = apply_uflowz(data_list.flat[i])
        
        return data

#------------------------------------------------------------------------

if __name__ == '__main__':
    print "Begin read block data for cone20 test case..."
    file_pointer = open("cone20.s", "r")
    block0 = BlockFlow2D()
    block1 = BlockFlow2D()
    for it in range(4):
        block0.read(file_pointer)
        block1.read(file_pointer)
        print "Have read data for time t=", block0.t
    print "Pressure and mass flux along cone surface:"
    for i in range(block1.ni):
        print block1.x[i][0], block1.y[i][0], block1.p[i][0], block1.f[i][0][0]
    file_pointer.close()
    print "Done"
