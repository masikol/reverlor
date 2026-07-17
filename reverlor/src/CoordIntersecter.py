#!/usr/bin/env python3

import os
import sys
import subprocess as sp
from typing import Optional


UID_SEP = '$/$'


class CoordIntersecter:
    def __init__(self,
                 bam_fpath: str,
                 samtools_fpath: str = 'samtools',
                 f_flags: Optional[list[int]] = None,
                 F_flags: Optional[list[int]] = None):

        if not os.path.isfile(bam_fpath):
            raise FileNotFoundError(
                'BAM file does not exist: `{}`'.format(bam_fpath)
            )
        # end if

        f_str = ' '.join('-f {}'.format(flag) for flag in (f_flags or []))
        F_str = ' '.join('-F {}'.format(flag) for flag in (F_flags or []))
        self._flag_str = ' '.join((f_str, F_str)).strip()
        self._bam_fpath = os.path.abspath(bam_fpath)
        self._samtools_fpath = samtools_fpath
    # end def

    def intersect_coords(self,
                         rname: str,
                         pos_1: int,
                         pos_2: int) -> set[str]:

        if pos_1 < 1:
            raise ValueError('pos_1 must be >= 1, got {}'.format(pos_1))
        # end if
        if pos_2 < 1:
            raise ValueError('pos_2 must be >= 1, got {}'.format(pos_2))
        # end if

        uids_1 = self._get_read_uids_at(rname, pos_1)
        uids_2 = self._get_read_uids_at(rname, pos_2)
        common_uids = uids_1 & uids_2

        return frozenset((
            self._uid_to_qname(uid) for uid in common_uids
        ))
    # end def

    def _get_read_uids_at(self,
                          rname: str,
                          pos: int) -> frozenset[str]:

        region = '{}:{}-{}'.format(rname, pos, pos)
        cmd_items = [
            self._samtools_fpath,
            'view',
            self._flag_str,
            self._bam_fpath,
            region,
        ]
        cmd = ' '.join(cmd_items)

        pipe = sp.Popen(
            cmd,
            shell=True,
            stdout=sp.PIPE,
            stderr=sp.PIPE,
            encoding='utf-8',
        )
        outs, errs = pipe.communicate()

        if pipe.returncode != 0:
            raise RuntimeError(
                'samtools view failed:\nCMD: {}\n{}'.format(cmd, errs)
            )
        # end if

        return frozenset(
            self._sam_line_to_uid(line) for line in outs.splitlines()
        )
    # end def

    @staticmethod
    def _sam_line_to_uid(line: str) -> str:
        vals = line.split('\t')[:4]
        return UID_SEP.join((
            vals[0].strip(),
            vals[1].strip(),
            vals[3].strip(),
        ))
    # end def

    @staticmethod
    def _uid_to_qname(uid: str) -> str:
        return uid.partition(UID_SEP)[0]
    # end def
# end class
