//
//  Copyright & License: See Copyright.readme in src directory
//

/*!
  Class documentation
*/

#ifndef VERSION_DEF_HPP_
#define VERSION_DEF_HPP_

#include "version.hpp"

#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_fusion.hpp>
#include <boost/spirit/include/phoenix_bind.hpp>

namespace SDDS { namespace parser
{
    template <typename Iterator>
    version_parser<Iterator>::version_parser(error_handler<Iterator> & _error_handler):
        version_parser::base_type(start)
    {
        using qi::on_error;
        using qi::fail;
        using phx::function;
        typedef function<error_handler<Iterator> > error_handler_function;

        // qi::_1_type _1;
        qi::_3_type _3;
        qi::_4_type _4;
        // qi::_val_type _val;
        qi::short_type short_;
        qi::no_skip_type no_skip;
        qi::lit_type lit;

        start %= no_skip[lit("SDDS") > short_];

        BOOST_SPIRIT_DEBUG_NODES(
            (start)
        )

        on_error<fail>(start,
            error_handler_function(_error_handler)(
                  std::string("Error! Expecting "), _4, _3));
    }
}}
#endif /* VERSION_DEF_HPP_ */