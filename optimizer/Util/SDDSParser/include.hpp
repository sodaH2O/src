//
//  Copyright & License: See Copyright.readme in src directory
//

/*!
  Class documentation
*/

#ifndef INCLUDE_HPP_
#define INCLUDE_HPP_

#include "ast.hpp"
#include "skipper.hpp"
#include "error_handler.hpp"

#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_fusion.hpp>
#include <boost/spirit/include/phoenix_bind.hpp>
#include <boost/fusion/include/adapt_struct.hpp>

#include <list>

#define BOOST_SPIRIT_NO_PREDEFINED_TERMINALS
#define BOOST_SPIRIT_QI_DEBUG

namespace SDDS {
    struct include
    {
        enum attributes { FILENAME
                        , INCLUDE
        };

        template <attributes A>
        struct complainUnsupported
        {
            static bool apply()
            {
                std::string attributeString;
                switch(A)
                {
                case FILENAME:
                    attributeString = "filename";
                    break;
                case INCLUDE:
                    attributeString = "include";
                    break;
                default:
                    return true;
                }
                std::cerr << attributeString << " not supported yet" << std::endl;
                return false;
            }
        };
    };

    struct includeList: std::list<include> {};

    inline std::ostream& operator<<(std::ostream& out, const include& ) {
        return out;
    }
}

namespace SDDS { namespace parser
{
    namespace qi = boost::spirit::qi;
    namespace ascii = boost::spirit::ascii;
    namespace phx = boost::phoenix;

    template <typename Iterator>
    struct include_parser: qi::grammar<Iterator, include(), skipper<Iterator> >
    {
        include_parser(error_handler<Iterator> & _error_handler);

        qi::rule<Iterator, include(), skipper<Iterator> > start;
        qi::rule<Iterator, std::string(), skipper<Iterator> > string, quoted_string,
                include_filename;
    };
}}
#endif /* INCLUDE_HPP_ */
