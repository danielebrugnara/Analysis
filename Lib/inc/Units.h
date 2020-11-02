#pragma once

namespace UNITS
{
//
// Length [L]
//
static const double meter   = 1.;
static const double meter2  = meter * meter;
static const double meter3  = meter * meter * meter;

static const double millimeter  = 1.E-3*meter;
static const double millimeter2 = millimeter * millimeter;
static const double millimeter3 = millimeter * millimeter * millimeter;

static const double centimeter  = 10. * millimeter;
static const double centimeter2 = centimeter * centimeter;
static const double centimeter3 = centimeter * centimeter * centimeter;


static const double kilometer   = 1000. * meter;
static const double kilometer2  = kilometer * kilometer;
static const double kilometer3  = kilometer * kilometer * kilometer;

static const double parsec = 3.0856775807e+16 * meter;

static const double micrometer  = 1.e-6 * meter;
static const double nanometer   = 1.e-9 * meter;
static const double angstrom    = 1.e-10 * meter;
static const double fermi       = 1.e-15 * meter;

static const double barn        = 1.e-28 * meter2;
static const double millibarn   = 1.e-3 * barn;
static const double microbarn   = 1.e-6 * barn;
static const double nanobarn    = 1.e-9 * barn;
static const double picobarn    = 1.e-12 * barn;

// symbols
static const double nm = nanometer;
static const double um = micrometer;

static const double mm  = millimeter;
static const double mm2 = millimeter2;
static const double mm3 = millimeter3;

static const double cm  = centimeter;
static const double cm2 = centimeter2;
static const double cm3 = centimeter3;

static const double m   = meter;
static const double m2  = meter2;
static const double m3  = meter3;

static const double km  = kilometer;
static const double km2 = kilometer2;
static const double km3 = kilometer3;

static const double pc = parsec;

//
// Angle
//
static const double radian = 1.;
static const double milliradian = 1.e-3 * radian;
static const double degree = (3.14159265358979323846 / 180.0) * radian;

static const double steradian = 1.;

// symbols
static const double rad = radian;
static const double mrad = milliradian;
static const double sr = steradian;
static const double deg = degree;

//
// Time [T]
//
static const double second      = 1.;
static const double millisecond = 1.e-3 * second;
static const double microsecond = 1.e-6 * second;
static const double nanosecond  = 1.e-9 * second;
static const double picosecond  = 1.e-12 * second;
static const double femtosecond = 1.e-15 * second;
static const double aptosecond  = 1.e-18 * second;
static const double zeptosecond = 1.e-21 * second;
static const double yoctosecond = 1.e-24 * second;

static const double minutes = 60 * second;
static const double hours   = 3600 * second;
static const double days    = 24 * minutes;
static const double years   = 365.25 * days;


static const double hertz       = 1. / second;
static const double kilohertz   = 1.e+3 * hertz;
static const double megahertz   = 1.e+6 * hertz;

// symbols
static const double ns  = nanosecond;
static const double us  = microsecond;
static const double s   = second;
static const double ms  = millisecond;
static const double ps  = picosecond;
static const double fs  = femtosecond;
static const double as  = aptosecond;
static const double zs  = zeptosecond;
static const double ys  = yoctosecond;


//
// Electric charge [Q]
//
static const double eplus   = 1.;             // positron charge
static const double e_SI    = 1.602176487e-19; // positron charge in coulomb
static const double coulomb = eplus / e_SI; // coulomb = 6.24150 e+18 * eplus

//
// Energy [E]
//
static const double megaelectronvolt    = 1.;
static const double electronvolt        = 1.e-6 * megaelectronvolt;
static const double kiloelectronvolt    = 1.e-3 * megaelectronvolt;
static const double gigaelectronvolt    = 1.e+3 * megaelectronvolt;
static const double teraelectronvolt    = 1.e+6 * megaelectronvolt;
static const double petaelectronvolt    = 1.e+9 * megaelectronvolt;

static const double joule = electronvolt / e_SI; // joule = 6.24150 e+12 * MeV

// symbols
static const double MeV = megaelectronvolt;
static const double eV  = electronvolt;
static const double keV = kiloelectronvolt;
static const double GeV = gigaelectronvolt;
static const double TeV = teraelectronvolt;
static const double PeV = petaelectronvolt;

//
// Mass [E][T^2][L^-2]
//
static const double kilogram    = joule * second * second / (meter2);
static const double gram        = 1.e-3 * kilogram;
static const double milligram   = 1.e-3 * gram;

// symbols
static const double kg  = kilogram;
static const double g   = gram;
static const double mg  = milligram;

//
// Power [E][T^-1]
//
static const double watt = joule / second; // watt = 6.24150 e+3 * MeV/ns

//
// Force [E][L^-1]
//
static const double newton = joule / meter; // newton = 6.24150 e+9 * MeV/mm

//
// Pressure [E][L^-3]
//
//#define pascal hep_pascal                          // a trick to avoid warnings
static const double hep_pascal  = newton / m2;     // pascal = 6.24150 e+3 * MeV/mm3
static const double pascal      = newton / m2;         // pascal = 6.24150 e+3 * MeV/mm3
static const double bar         = 100000 * pascal;        // bar    = 6.24150 e+8 * MeV/mm3
static const double atmosphere  = 101325 * pascal; // atm    = 6.32420 e+8 * MeV/mm3
static const double torr        = 0.00133322 * bar;

//
// Electric current [Q][T^-1]
//
static const double ampere      = coulomb / second; // ampere = 6.24150 e+9 * eplus/ns
static const double milliampere = 1.e-3 * ampere;
static const double microampere = 1.e-6 * ampere;
static const double nanoampere  = 1.e-9 * ampere;

//
// Electric potential [E][Q^-1]
//
static const double megavolt    = megaelectronvolt / eplus;
static const double kilovolt    = 1.e-3 * megavolt;
static const double volt        = 1.e-6 * megavolt;

//
// Electric resistance [E][T][Q^-2]
//
static const double ohm = volt / ampere; // ohm = 1.60217e-16*(MeV/eplus)/(eplus/ns)

//
// Electric capacitance [Q^2][E^-1]
//
static const double farad       = coulomb / volt; // farad = 6.24150e+24 * eplus/Megavolt
static const double millifarad  = 1.e-3 * farad;
static const double microfarad  = 1.e-6 * farad;
static const double nanofarad   = 1.e-9 * farad;
static const double picofarad   = 1.e-12 * farad;

//
// Magnetic Flux [T][E][Q^-1]
//
static const double weber = volt * second; // weber = 1000*megavolt*ns

//
// Magnetic Field [T][E][Q^-1][L^-2]
//
static const double tesla = volt * second / meter2; // tesla =0.001*megavolt*ns/mm2

static const double gauss       = 1.e-4 * tesla;
static const double kilogauss   = 1.e-1 * tesla;

//
// Inductance [T^2][E][Q^-2]
//
static const double henry = weber / ampere; // henry = 1.60217e-7*MeV*(ns/eplus)**2

//
// Temperature
//
static const double kelvin = 1.;

//
// Amount of substance
//
static const double mole = 1.;

//
// Activity [T^-1]
//
static const double becquerel   = 1. / second;
static const double curie       = 3.7e+10 * becquerel;

//
// Absorbed dose [L^2][T^-2]
//
static const double gray        = joule / kilogram;
static const double kilogray    = 1.e+3 * gray;
static const double milligray   = 1.e-3 * gray;
static const double microgray   = 1.e-6 * gray;

//
// Luminous intensity [I]
//
static const double candela = 1.;

//
// Luminous flux [I]
//
static const double lumen = candela * steradian;

//
// Illuminance [I][L^-2]
//
static const double lux = lumen / meter2;

//
// Miscellaneous
//
static const double perCent     = 0.01;
static const double perThousand = 0.001;
static const double perMillion  = 0.000001;

namespace CONSTANTS
{
static const double pi      = 3.1415926535897932384626433832795028841971693993751;
static const double halfpi  = 3.1415926535897932384626433832795028841971693993751 / 2.;
} // namespace CONSTANTS

namespace  PHYSICS
{
static const double emass_SI    = 9.1093837015E-31 * kg;
static const double emass_c2    = 0.51099895000 * MeV;
static const double echarge_SI  = e_SI;
static const double pmass_SI    = 1.67262192369E-27 * kg;
static const double pmass_c2    = 938.27208816 * MeV;
static const double nmass_SI    = 1.674927351E-27 * kg;
static const double nmass_c2    = 939.565378 * MeV;
static const double amu_c2      = 931.49432 * MeV;
static const double c_light     = 299792458 * m/s;

}

} // namespace UNITS