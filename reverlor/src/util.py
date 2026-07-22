
import os
import sys
import logging


def rm_files_if_exist(*fpaths: str) -> None:
    for fpath in fpaths:
        if os.path.isfile(fpath):
            try:
                os.unlink(fpath)
            except OSError as err:
                logging.warning('Error: cannot remove temp file `{}`'.format(fpath))
                logging.warning(str(err))
                logging.warning('Ignoring...')
            # end try
        # end if
    # end for
# end def


def rm_empty_dir_if_exists(dir_path: str) -> None:
    try:
        os.rmdir(dir_path)
    except OSError as err:
        logging.warning('Error: cannot remove temp dir `{}`'.format(dir_path))
        logging.warning(str(err))
        logging.warning('Ignoring...')
    # end try
# end def
