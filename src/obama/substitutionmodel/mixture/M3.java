package obama.substitutionmodel.mixture;

public class M3 extends M1 {
	@Override
	public double[][] getEmpiricalRates() {
        double[][] rates = new double[][] {
{0, 0.733336, 0.558955, 0.503360, 4.149599, 1.415369, 1.367574, 1.263002, 0.994098, 0.517204, 0.775054, 0.763094, 1.890137, 0.540460, 0.200122, 4.972745, 1.825593, 0.450842, 0.526135, 3.839269},
{0.733336, 0, 0.597671, 0.058964, 2.863355, 2.872594, 0.258365, 0.366868, 2.578946, 0.358350, 0.672023, 5.349861, 0.691594, 0.063347, 0.032875, 0.821562, 0.580847, 0.661866, 0.265730, 0.395134},
{0.558955, 0.597671, 0, 5.581680, 1.279881, 1.335650, 0.397108, 1.840061, 5.739035, 0.284730, 0.109781, 1.612642, 0.466979, 0.141582, 0.019509, 4.670980, 1.967383, 0.088064, 0.581928, 0.145401},
{0.503360, 0.058964, 5.581680, 0, 0.225860, 0.434096, 2.292917, 1.024707, 0.821921, 0.027824, 0.021443, 0.088850, 0.060820, 0.018288, 0.042687, 1.199607, 0.420710, 0.037642, 0.141233, 0.090101},
{4.149599, 2.863355, 1.279881, 0.225860, 0, 1.043232, 0.209978, 0.823594, 3.039380, 1.463390, 1.983693, 0.397640, 2.831098, 4.102068, 0.059723, 5.901348, 2.034980, 2.600668, 5.413080, 4.193725},
{1.415369, 2.872594, 1.335650, 0.434096, 1.043232, 0, 4.534772, 0.377181, 4.877840, 0.370939, 1.298542, 3.509873, 2.646440, 0.087872, 0.072299, 1.139018, 0.864479, 0.390688, 0.322761, 0.625409},
{1.367574, 0.258365, 0.397108, 2.292917, 0.209978, 4.534772, 0, 0.496780, 0.532488, 0.232460, 0.169219, 0.755219, 0.379926, 0.020447, 0.023282, 0.503875, 0.577513, 0.109318, 0.153776, 0.696533},
{1.263002, 0.366868, 1.840061, 1.024707, 0.823594, 0.377181, 0.496780, 0, 0.398817, 0.008940, 0.043707, 0.436013, 0.087640, 0.064863, 0.036426, 1.673207, 0.124068, 0.218118, 0.039217, 0.104335},
{0.994098, 2.578946, 5.739035, 0.821921, 3.039380, 4.877840, 0.532488, 0.398817, 0, 0.349195, 0.838324, 0.888693, 0.488389, 1.385133, 0.050226, 0.962470, 0.502294, 1.065585, 8.351808, 0.377304},
{0.517204, 0.358350, 0.284730, 0.027824, 1.463390, 0.370939, 0.232460, 0.008940, 0.349195, 0, 5.102837, 0.561690, 7.010411, 3.054968, 0.039318, 0.204155, 2.653232, 0.564368, 0.854294, 15.559906},
{0.775054, 0.672023, 0.109781, 0.021443, 1.983693, 1.298542, 0.169219, 0.043707, 0.838324, 5.102837, 0, 0.401070, 8.929538, 5.525874, 0.067505, 0.273372, 0.437116, 1.927515, 0.940458, 2.508169},
{0.763094, 5.349861, 1.612642, 0.088850, 0.397640, 3.509873, 0.755219, 0.436013, 0.888693, 0.561690, 0.401070, 0, 1.357738, 0.043394, 0.023126, 0.567639, 1.048288, 0.120994, 0.180650, 0.449074},
{1.890137, 0.691594, 0.466979, 0.060820, 2.831098, 2.646440, 0.379926, 0.087640, 0.488389, 7.010411, 8.929538, 1.357738, 0, 3.135353, 0.012695, 0.570771, 2.319555, 1.856122, 0.975427, 3.404087},
{0.540460, 0.063347, 0.141582, 0.018288, 4.102068, 0.087872, 0.020447, 0.064863, 1.385133, 3.054968, 5.525874, 0.043394, 3.135353, 0, 0.015631, 0.458799, 0.151684, 4.154750, 11.429924, 1.457957},
{0.200122, 0.032875, 0.019509, 0.042687, 0.059723, 0.072299, 0.023282, 0.036426, 0.050226, 0.039318, 0.067505, 0.023126, 0.012695, 0.015631, 0, 0.233109, 0.077004, 0.011074, 0.026268, 0.052132},
{4.972745, 0.821562, 4.670980, 1.199607, 5.901348, 1.139018, 0.503875, 1.673207, 0.962470, 0.204155, 0.273372, 0.567639, 0.570771, 0.458799, 0.233109, 0, 8.113282, 0.377578, 0.429221, 0.260296},
{1.825593, 0.580847, 1.967383, 0.420710, 2.034980, 0.864479, 0.577513, 0.124068, 0.502294, 2.653232, 0.437116, 1.048288, 2.319555, 0.151684, 0.077004, 8.113282, 0, 0.222293, 0.273138, 2.903836},
{0.450842, 0.661866, 0.088064, 0.037642, 2.600668, 0.390688, 0.109318, 0.218118, 1.065585, 0.564368, 1.927515, 0.120994, 1.856122, 4.154750, 0.011074, 0.377578, 0.222293, 0, 4.731579, 0.564762},
{0.526135, 0.265730, 0.581928, 0.141233, 5.413080, 0.322761, 0.153776, 0.039217, 8.351808, 0.854294, 0.940458, 0.180650, 0.975427, 11.429924, 0.026268, 0.429221, 0.273138, 4.731579, 0, 0.681215},
{3.839269, 0.395134, 0.145401, 0.090101, 4.193725, 0.625409, 0.696533, 0.104335, 0.377304, 15.559906, 2.508169, 0.449074, 3.404087, 1.457957, 0.052132, 0.260296, 2.903836, 0.564762, 0.681215, 0},
        };
		return rates;
	}

	@Override
	public double[] getEmpiricalFrequencies() {
        double[] f = new double[] {0.062457,0.066826,0.049332,0.065270,0.006513,0.041231,0.058965,0.080852,0.028024,0.037024,0.075925,0.064131,0.019620,0.028710,0.104579,0.056388,0.062027,0.008241,0.033124,0.050761};
        return f;
	}

}
