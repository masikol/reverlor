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


# >>> Import functions >>>
from src.FindArgs import FindArgs
from src.find_repeats import find_repeats as find_repeats_minimap2
from src.find_repeats_RepeatScout import find_repeats as find_repeats_repeatscout
from src.find_repeats_phRAIDER import find_repeats as find_repeats_phraider
# <<<

# >>> Functions >>>

def reverlor_find():
    args = FindArgs.parse_args()
    if args.finder == 'repeatscout':
        repeats_bed = find_repeats_repeatscout(args)
    elif args.finder == 'phraider':
        repeats_bed = find_repeats_phraider(args)
    else:
        repeats_bed = find_repeats_minimap2(args)
    # end if
# end def

# <<<


# >>> Proceed >>>

if __name__ == '__main__':
    reverlor_find()
# end if

# <<<

sys.exit(0)
