#ifndef SWITCHERERROR_H
#define SWITCHERERROR_H

#include "Utilities/ClassicException.h"

class SwitcherError:public ClassicException
{
public:
    SwitcherError(const string &meth, const string &msg);

    SwitcherError(const SwitcherError &);
    virtual ~SwitcherError();

private:

    // Not implemented.
    SwitcherError();    
};

#endif
