#!/usr/bin/env python3

# >>> Import >>>

import sys
import logging

from .src.FindArgs import FindArgs
from .src._version import report_version_and_author
from .src.find_repeats_minimap2 import find_repeats as find_repeats_minimap2

# <<<

# >>> Functions >>>

def main():
    args = FindArgs.parse_args()
    reverlor_find(args)
# end def

def reverlor_find(args: FindArgs):
    log_level = logging.getLogger().level
    if log_level in (logging.DEBUG, logging.INFO):
        report_version_and_author()
        logging.info(args)
    # end if

    logging.info('Repeat search started')
    repeat_bed_fpath = find_repeats_minimap2(args)

    logging.info('Repeat search completed!')
    logging.info(f'Output directory: `{args.output_dir}`')
    logging.info(f'Repeats are listed in this BED file: `{repeat_bed_fpath}`')
# end def

# <<<


# >>> Proceed >>>

if __name__ == '__main__':
    main()
    sys.exit(0)
# end if
