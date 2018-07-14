//
//  Copyright & License: See Copyright.readme in src directory
//

/*!
  Class documentation
*/

#ifndef VERSION_HPP_
#define VERSION_HPP_

#include "ast.hpp"
#include "skipper.hpp"
#include "error_handler.hpp"

#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/fusion/include/adapt_struct.hpp>

#define BOOST_SPIRIT_NO_PREDEFINED_TERMINALS
#define BOOST_SPIRIT_QI_DEBUG

namespace SDDS {
    struct version
    {
        short layoutVersion_m;
    };

    inline std::ostream& operator<<(std::ostream& out, const version& head) {
        out << "layout version is " << head.layoutVersion_m;
        return out;
    }
}

BOOST_FUSION_ADAPT_STRUCT(
    SDDS::version,
    (short, layoutVersion_m)
)

namespace SDDS { namespace parser
{
    namespace qi = boost::spirit::qi;
    namespace ascii = boost::spirit::ascii;
    namespace phx = boost::phoenix;

    template <typename Iterator>
    struct version_parser: qi::grammar<Iterator, version(), skipper<Iterator> >
    {
        version_parser(error_handler<Iterator> & _error_handler);

        qi::rule<Iterator, version(), skipper<Iterator> > start;
    };
}}
#endif /* VERSION_HPP_ */
