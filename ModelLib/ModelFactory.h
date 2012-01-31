
#pragma once

#include "IModel.h"
#include "GROUNDWATER_FLOW.h"

namespace ModelLib
{

struct ModelType
{
    enum type {
        GROUNDWATER_FLOW,
        INVALID
    };
};


class ModelFactory
{
public:
    static IModel* createModel(ModelType::type model_type) 
    {
        switch (model_type) {
        case ModelType::GROUNDWATER_FLOW:
            return new GROUNDWATER_FLOW();
        }
        return 0;
    };
};
}
