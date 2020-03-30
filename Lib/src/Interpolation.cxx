#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include "Interpolation.h"

//ClassImp(Interpolation);

Interpolation::Interpolation(std::string input_file) : spline(nullptr)
{
    if (!ReadFile(input_file))
        throw std::runtime_error("Unable to generate interpolation from text file\n");
}

Interpolation::Interpolation(TFile *file) : spline(nullptr)
{
    if (!(file->IsOpen()))
        throw std::runtime_error("Interpolation file pointer non initialized\n");

    TIter contents(file->GetListOfKeys());
    TKey *key;
    TObject *obj;
    while ((key = (TKey *)contents()))
    {
        obj = file->Get(key->GetName());
        if (obj->InheritsFrom("TSpline"))
        {
            TSpline *tmp = (TSpline *)obj;
            spline = tmp;
            return;
        }
    }
    if (!spline)
        throw std::runtime_error(std::string("Interpolation not found in TFile: ") +
                                 file->GetName() +
                                 "\n");
}

Interpolation::~Interpolation()
{
    delete spline;
}

bool Interpolation::ReadFile(std::string input_file_name)
{
    std::ifstream in_file(input_file_name);
    if (!in_file.is_open())
        return false;
    std::string line;
    std::string delimiter = ", ";
    int ii = -1;
    std::vector<double> fValues_x;
    std::vector<double> fValues_y;
    while (std::getline(in_file, line))
    {
        ii++;
        if (ii == 0)
            continue;
        std::string x_string = line.substr(2, line.find(", ") - 2);
        std::string y_string = line.substr(line.find(", ") + 2, line.rfind("}") - line.find(", ") - 3);
        double x, y;
        try
        {
            x = std::stof(x_string);
            y = std::stof(y_string);
            if (ii > 1 && fValues_x.back() > x)
                break;
            fValues_x.push_back(x);
            fValues_y.push_back(y);
        }
        catch (std::invalid_argument &err)
        {
            std::cerr << "Invalid arg :" << err.what() << "\n";
            //            return false;
        }
        spline = new TSpline3(input_file_name.c_str(), &fValues_x[0], &fValues_y[0], fValues_x.size());
    }
    return true;
}

TSpline *Interpolation::GetSpline()
{
    return spline;
}

#endif