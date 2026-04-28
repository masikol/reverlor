#!/usr/bin/env python3

__author__ = 'Maksim Sikolenko'
__email__ = 'sikolenko@mbio.bas-net.by'
__version__ = '0.0.a'
__last_update_date__ = None


# >>> Check python interpreter version >>>

import sys

if sys.version_info.major < 3:
    print(
        "\nYour python interpreter version is " + "%d.%d" % (
            sys.version_info.major,
            sys.version_info.minor
        )
    )
    print("   Please, use Python 3.\a")
    # In python 2 'raw_input' does the same thing as 'input' in python 3.
    # Neither does 'input' in python2.
    if sys.platform.startswith("win"):
        raw_input("Press ENTER to exit:")
    # end if
    sys.exit(1)
# end if

# <<<


# >>> Import public functions >>>
from reverlor_find import reverlor_find
from reverlor_verify import reverlor_verify
# <<<


# >>> Functions >>>

def main():
    reverlor()
# end def

def reverlor():
    reverlor_find()
    reverlor_verify()
# end def

# <<<


# >>> Proceed >>>

if __name__ == '__main__':
    main()
# end if

# <<<

sys.exit(0)
