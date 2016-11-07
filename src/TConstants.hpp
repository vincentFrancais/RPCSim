// TConstants.hpp
// 
// Copyright 2016 Vincent Français <francais@clermont.in2p3.fr>
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
// MA 02110-1301, USA.
// 
// 

namespace Constants {
	
	static const double Pi = 3.1415926535897932384626433832795;
	static const double TwoPi = 2. * Pi;
	static const double HalfPi = 0.5 * Pi;
	static const double Pi2 = Pi * Pi;
	
	// Elementary particle masses [MeV / c2]
	static const double ElectronMass = 0.510998910e6;
	static const double MuonMass = 105.658367e6;
	static const double ProtonMass = 938.272013e6;
	static const double NeutronMass = 939.56536e6;
	
	static const double VacuumPermittivity = 8.854187817e-12; /* A^2 s^4 / kg m^3 */
	
	static const double ElectronCharge = 1.602176487e-19; /* Coulomb */
	
	/* Conversions to m */
	static const double cm = 1.e-2;
	static const double mm = 1.e-3;
	
}
