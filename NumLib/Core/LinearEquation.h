
#pragma once

namespace NumLib
{

template<typename Tmat, typename Tvec>
class LinearEquation
{
private:
    Tmat* _A;
    Tvec* _RHS;
public:

    Tmat* getA() {return _A;};
    Tvec* getRHS() {return _RHS;}

    void initialize() 
    {
        throw std::exception("The method or operation is not implemented.");
    }

    void reset() 
    {
        throw std::exception("The method or operation is not implemented.");
    }

    template<typename T2mat, typename T2vec>
    void add( LinearEquation<T2mat,T2vec> &localEQS, std::vector<size_t> &dofmap ) 
    {
        throw std::exception("The method or operation is not implemented.");
    }

    double* getX() 
    {
        throw std::exception("The method or operation is not implemented.");
    }





};

}
