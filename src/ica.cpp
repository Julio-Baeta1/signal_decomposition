#include "ica.h"

Ica::Ica(std::shared_ptr<Mat> initial_S)
{
    S = initial_S;
}

void Ica::setSource(std::shared_ptr<Mat> new_S)
{
    S = new_S;
}


void Ica::whiten()
{
    int x=1;
}

void Ica::decompose()
{
    int x=1;
}