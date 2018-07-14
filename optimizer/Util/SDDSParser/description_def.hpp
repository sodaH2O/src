//
//  Copyright & License: See Copyright.readme in src directory
//

/*!
  Class documentation
*/

#ifndef DESCRIPTION_DEF_HPP_
#define DESCRIPTION_DEF_HPP_

#include "description.hpp"

#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_fusion.hpp>
#include <boost/spirit/include/phoenix_bind.hpp>

namespace SDDS { namespace parser
{
    template <typename Iterator>
    description_parser<Iterator>::description_parser(error_handler<Iterator> & _error_handler):
        description_parser::base_type(start)
    {
        using qi::on_error;
        using qi::fail;
        using phx::function;
        typedef function<error_handler<Iterator> > error_handler_function;

        qi::_1_type _1;
        qi::_3_type _3;
        qi::_4_type _4;
        qi::lexeme_type lexeme;
        qi::lit_type lit;
        qi::char_type char_;
        qi::_val_type _val;

        quoted_string %= lexeme['"' >> +(char_ - '"') >> '"'];
        description_text = lit("text") > '=' > quoted_string;
        description_content = lit("contents") > '=' > quoted_string;

        start =
                lit("&description")
                > ((description_text[phx::at_c<0>(_val) = _1]
                   > -(',' > description_content[phx::at_c<1>(_val) = _1]))
                  |(description_content[phx::at_c<1>(_val) = _1]
                    > -(',' > description_text[phx::at_c<0>(_val) = _1])))
                > lit("&end");

        BOOST_SPIRIT_DEBUG_NODES(
            (start)
        )

        on_error<fail>(start,
            error_handler_function(_error_handler)(
                                                   std::string("Error! Expecting "), _4, _3));
    }
}}
#endif /* DESCRIPTION_DEF_HPP_ */