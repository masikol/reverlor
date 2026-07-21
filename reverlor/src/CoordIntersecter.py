#!/usr/bin/env python3

import os
import sys
import subprocess as sp
from typing import Optional

import pysam


UID_SEP = '$/$'


class CoordIntersecter:
    def __init__(self,
                 bam_fpath: str,
                 f_flags: Optional[list[int]] = None,
                 F_flags: Optional[list[int]] = None):

        if not os.path.isfile(bam_fpath):
            raise FileNotFoundError(
                'BAM file does not exist: `{}`'.format(bam_fpath)
            )
        # end if

        self._bam_file = pysam.AlignmentFile(bam_fpath, 'rb')

        if f_flags is None or len(f_flags) == 0:
            self._f_flags = None
        else:
            self._f_flags = f_flags
        # end if

        if F_flags is None or len(F_flags) == 0:
            self._F_flags = None
        else:
            self._F_flags = F_flags
        # end if
    # end def

    def __del__(self):
        self._bam_file.close()
    # end def


    def intersect_coords(self,
                         rname: str,
                         pos_1: int,
                         pos_2: int) -> set[str]:

        if pos_1 < 0:
            raise ValueError('pos_1 must be >= 0, got {}'.format(pos_1))
        # end if
        if pos_2 < 0:
            raise ValueError('pos_2 must be >= 0, got {}'.format(pos_2))
        # end if

        uids_1 = self._get_read_uids_at(rname, pos_1)
        uids_2 = self._get_read_uids_at(rname, pos_2)
        common_uids = uids_1 & uids_2

        return frozenset((
            _uid_to_qname(uid) for uid in common_uids
        ))
    # end def

    def _get_read_uids_at(self,
                          rname: str,
                          pos: int) -> frozenset[str]:
        pos_segments = filter(
            self._segment_passes_flag_filter,
            self._bam_file.fetch(rname, pos, pos + 1) # +1 becaulse end is open (exclusive)
        )

        pos_uids = map(
            _sam_line_to_uid,
            pos_segments
        )

        return frozenset(pos_uids)
    # end def

    def _segment_passes_flag_filter(self, seg: pysam.AlignedSegment) -> bool:

        if self._f_flags is not None:
            f_fails = not all(map(
                lambda flag: (flag & seg.flag) == flag,
                self._f_flags
            ))
            if f_fails:
                return False
            # end if
        # end if

        if self._F_flags is not None:
            F_fails = any(map(
                lambda flag: (flag & seg.flag) == flag,
                self._F_flags
            ))
            if F_fails:
                return False
            # end if
        # end if

        return True
    # end def
# end class

def _sam_line_to_uid(read: pysam.AlignedSegment) -> str:
    return UID_SEP.join((
        read.query_name.strip(),
        str(read.flag).strip(),
        str(read.reference_start).strip(),
    ))
# end def

def _uid_to_qname(uid: str) -> str:
    return uid.partition(UID_SEP)[0]
# end def
