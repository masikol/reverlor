
import sys


__version__ = '0.0.1'
__last_update_date__ = '2026-07-24'

__author__ = 'Maksim Sikolenko'


def report_version_and_author() -> None:
    ascii_art = r"""
        ____  _______    ____________  __    ____  ____ 
       / __ \/ ____/ |  / / ____/ __ \/ /   / __ \/ __ \
      / /_/ / __/  | | / / __/ / /_/ / /   / / / / /_/ /
     / _, _/ /___  | |/ / /___/ _, _/ /___/ /_/ / _, _/ 
    /_/ |_/_____/  |___/_____/_/ |_/_____/\____/_/ |_|  
                                                        
"""
    sys.stdout.write(ascii_art)
    sys.stderr.write('REVERLOR:\n')
    sys.stderr.write('A program for REpeat VERification using LOng Reads\n')
    sys.stderr.write('Version: {} ({} edition)\n'.format(
        __version__,
        __last_update_date__
    ))
    sys.stderr.write('By {}\n'.format(
        __author__
    ))
    sys.stderr.write('https://github.com/masikol/reverlor\n')
    sys.stderr.write('=' * 20 + '\n')
    sys.stderr.flush()
# end def
