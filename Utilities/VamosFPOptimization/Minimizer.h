

class Minimizer{
    public:
        Minimizer(double, bool, double,  double, double, double, int, double);
        ~Minimizer();

        inline double Step(double );
        double GetCoefficient(){return _Coefficient;};
    private:
        std::array <double, 2> _Y;
        std::array <double, 2> _X;
        std::array < long double, 2> _Derivative;
        std::array <long double, 2> _Step;
        std::array <long double, 2> _Rate;

        int _NSteps;
        bool _Verbose;
        double _StartingRate;
        double _Coefficient;
        double _Threshold;
        double _Quenching;
        double _MaxSteps;
};