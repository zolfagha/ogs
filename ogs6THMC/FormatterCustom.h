
#pragma once

#include "logog/include/logog.hpp"
#include "logog/include/formatter.hpp"

/**
 * new formatter for logog
 */
class FormatterCustom : public logog::FormatterGCC
{
    virtual TOPIC_FLAGS GetTopicFlags( const logog::Topic &topic )
    {
        return ( Formatter::GetTopicFlags( topic ) &
                 ~( TOPIC_LEVEL_FLAG | TOPIC_FILE_NAME_FLAG | TOPIC_LINE_NUMBER_FLAG ));
    }
};
