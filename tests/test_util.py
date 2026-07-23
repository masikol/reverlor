import os

import pytest

from reverlor.src.util import rm_files_if_exist, rm_empty_dir_if_exists


# >>> rm_files_if_exist tests >>>

def test_single_existing_file_is_removed(tmp_path):
    fpath = tmp_path / 'a.txt'
    fpath.write_text('hello')
    assert fpath.is_file()
    rm_files_if_exist(str(fpath))
    assert not fpath.is_file()
# end def


def test_multiple_existing_files_are_removed(tmp_path):
    files = [tmp_path / f'file_{i}.txt' for i in range(4)]
    for f in files:
        f.write_text('x')
    rm_files_if_exist(*[str(f) for f in files])
    for f in files:
        assert not f.is_file()
# end def


def test_nonexistent_file_is_silent(tmp_path):
    fpath = tmp_path / 'does_not_exist.txt'
    rm_files_if_exist(str(fpath))
    assert not fpath.is_file()
# end def


def test_mixed_existing_and_nonexistent(tmp_path):
    existing = tmp_path / 'exists.txt'
    missing = tmp_path / 'missing.txt'
    existing.write_text('y')
    rm_files_if_exist(str(existing), str(missing))
    assert not existing.is_file()
    assert not missing.is_file()
# end def


def test_directory_is_not_removed(tmp_path):
    d = tmp_path / 'mydir'
    d.mkdir()
    rm_files_if_exist(str(d))
    assert d.is_dir()
# end def


def test_symlink_to_file_is_removed(tmp_path):
    target = tmp_path / 'target.txt'
    target.write_text('data')
    link = tmp_path / 'link.txt'
    link.symlink_to(target)
    rm_files_if_exist(str(link))
    assert not link.is_symlink()
    assert target.is_file()
# end def


def test_readonly_file_logs_warning(tmp_path, monkeypatch, caplog):
    def boom(path):
        raise OSError('permission denied')
    # end def

    monkeypatch.setattr(os, 'unlink', boom)

    fpath = tmp_path / 'cant_delete.txt'
    fpath.write_text('data')
    rm_files_if_exist(str(fpath))
    assert fpath.is_file()
    assert 'cannot remove temp file' in caplog.text
# end def


def test_no_args_does_nothing(tmp_path):
    rm_files_if_exist()
# end def


def test_unlink_called_for_each_file(tmp_path, monkeypatch):
    calls = []
    original_unlink = os.unlink

    def mock_unlink(path):
        calls.append(path)
        return original_unlink(path)
    # end def

    monkeypatch.setattr(os, 'unlink', mock_unlink)

    a = tmp_path / 'a.txt'
    b = tmp_path / 'b.txt'
    a.write_text('1')
    b.write_text('2')
    rm_files_if_exist(str(a), str(b))
    assert len(calls) == 2
    assert str(a) in calls
    assert str(b) in calls
# end def


def test_unlink_error_is_caught(tmp_path, monkeypatch):
    def boom(path):
        raise OSError('permission denied')
    # end def

    monkeypatch.setattr(os, 'unlink', boom)

    fpath = tmp_path / 'cant_delete.txt'
    fpath.write_text('data')
    rm_files_if_exist(str(fpath))
# end def


# <<< rm_files_if_exist tests <<<


# >>> rm_empty_dir_if_exists tests >>>

def test_empty_dir_is_removed(tmp_path):
    d = tmp_path / 'empty'
    d.mkdir()
    rm_empty_dir_if_exists(str(d))
    assert not d.exists()
# end def


def test_nonexistent_dir_is_silent(tmp_path):
    d = tmp_path / 'nope'
    rm_empty_dir_if_exists(str(d))
    assert not d.exists()
# end def


def test_nonempty_dir_is_not_removed(tmp_path):
    d = tmp_path / 'full'
    d.mkdir()
    (d / 'file.txt').write_text('data')
    rm_empty_dir_if_exists(str(d))
    assert d.is_dir()
# end def


def test_file_path_is_not_removed(tmp_path):
    f = tmp_path / 'not_a_dir.txt'
    f.write_text('hello')
    rm_empty_dir_if_exists(str(f))
    assert f.is_file()
# end def


def test_rmdir_error_is_caught(tmp_path, monkeypatch):
    def boom(path):
        raise OSError('permission denied')
    # end def

    monkeypatch.setattr(os, 'rmdir', boom)

    d = tmp_path / 'cant_delete'
    d.mkdir()
    rm_empty_dir_if_exists(str(d))
    assert d.is_dir()
# end def


def test_warning_message_on_failure(tmp_path, monkeypatch, caplog):
    def boom(path):
        raise OSError('disk full')
    # end def

    monkeypatch.setattr(os, 'rmdir', boom)

    d = tmp_path / 'fail'
    d.mkdir()
    rm_empty_dir_if_exists(str(d))
    assert 'cannot remove temp dir' in caplog.text
    assert 'disk full' in caplog.text
# end def


def test_os_rmdir_called(tmp_path, monkeypatch):
    calls = []
    original_rmdir = os.rmdir

    def mock_rmdir(path):
        calls.append(path)
        return original_rmdir(path)
    # end def

    monkeypatch.setattr(os, 'rmdir', mock_rmdir)

    d = tmp_path / 'target'
    d.mkdir()
    rm_empty_dir_if_exists(str(d))
    assert len(calls) == 1
    assert calls[0] == str(d)
# end def


# <<< rm_empty_dir_if_exists tests <<<
