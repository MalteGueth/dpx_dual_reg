# -*- coding: utf-8 -*-
"""Utility functions.

Authors: José C. García Alanis <alanis.jcg@gmail.com>

License: BSD (3-clause)
"""
import os
import warnings
import re

import numpy as np

import string


class FileNames(object):
    """
    Utility class to manage filenames.

    Notes:
    ------
    Use the `add` method to add new filenames. You specify a short "alias" for
    them, which you can use to retrieve the full filename later:
    
        >>> fname = FileNames()
        >>> fname.add('my_file', '/path/to/file1')
        >>> fname.my_file
        '/path/to/file1'

    Filenames can also be templates that can be used to generate
    filenames for different subjects, conditions, etc.:
    
        >>> fname = FileNames()
        >>> fname.add('epochs', '/data/{subject}/{cond}-epo.fif')
        >>> fname.epochs(subject='sub001', cond='face')
        '/data/sub001/face-epo.fif'
        
    Templates can contain placeholders in the way `string.format` allows,
    including formatting options:
    
        >>> fname = FileNames()
        >>> fname.add('epochs', '/data/sub{subject:03d}/{cond}-epo.fif')
        >>> fname.epochs(subject=1, cond='face')
        '/data/sub001/face-epo.fif'
        
    If a placeholder happens to be the alias of a file that has been added earlier,
    the placeholder is automatically filled:
    
        >>> fname = FileNames()
        >>> fname.add('subjects', '/data/subjects_dir')
        >>> fname.add('epochs', '{subjects}/{subject}/{cond}-epo.fif')
        >>> fname.epochs(subject='sub001', cond='face')
        '/data/subjects_dir/sub001/face-epo.fif'

    If all placeholders could be automatically filled, no brackets () are required
    when accessing it:

        >>> fname = FileNames()
        >>> fname.add('subjects', '/data/subjects_dir')
        >>> fname.add('fsaverage', '{subjects}/fsaverage-src.fif')
        >>> fname.fsaverage
        '/data/subjects_dir/fsaverage-src.fif'

    If computing the file path gets more complicated than the cases above, you can
    supply your own function. When the filename is requested, your function will
    get called with the FileNames object as first parameter, followed by any
    parameters that were supplied along with the request:

        >>> fname = FileNames()
        >>> fname.add('basedir', '/data/subjects_dir')
        >>> def my_function(files, subject):
        ...     if subject == 1:
        ...         return files.basedir + '/103hdsolli.fif'
        ...     else:
        ...         return files.basedir + '/%s.fif' % subject
        >>> fname.add('complicated', my_function)
        >>> fname.complicated(subject=1)
        '/data/subjects_dir/103hdsolli.fif'

    Author: Marijn van Vliet <w.m.vanvliet@gmail.com>
    """   # noqa: E501

    def files(self):
        """Obtain a list of file aliases known to this FileNames object.
        Returns
        -------
        files : list of str
            The list of file aliases.
        """
        files = dict()
        for name, value in self.__dict__.items():
            public_methods = ['list_filenames', 'add']
            if not name.startswith('_') and name not in public_methods:
                files[name] = value
        return files

    def add(self, alias, fname):
        """Add a new filename.
        Parameters
        ----------
        alias : str
            A short alias for the full filename. This alias can later be used
            to retrieve the filename. Aliases can not start with '_' or a
            number.
        fname : str | function
            The full filename. Either a string, with possible placeholder
            values, or a function that will compute the filename. If you
            specify a function, it will get called with the FileNames object as
            first parameter, followed by any parameters that were supplied
            along with the request.
        """
        if callable(fname):
            self._add_function(alias, fname)
        else:
            # Determine whether the string contains placeholders and whether
            # all placeholders can be pre-filled with existing file aliases.
            placeholders = _get_placeholders(fname)
            if len(placeholders) == 0:
                self._add_fname(alias, fname)  # Plain string filename
            else:
                prefilled = _prefill_placeholders(placeholders, self.files(),
                                                  dict())
                if len(prefilled) == len(placeholders):
                    # The template could be completely pre-filled. Add the
                    # result as a plain string filename.
                    self._add_fname(alias, fname.format(**prefilled))
                else:
                    # Add filename as a template
                    self._add_template(alias, fname)

    def _add_fname(self, alias, fname):
        """Add a filename that is a plain string."""
        self.__dict__[alias] = fname

    def _add_template(self, alias, template):
        """Add a filename that is a string containing placeholders."""
        # Construct a function that will do substitution for any placeholders
        # in the template.
        def fname(**kwargs):
            return _substitute(template, self.files(), kwargs)

        # Bind the fname function to this instance of FileNames
        self.__dict__[alias] = fname

    def _add_function(self, alias, func):
        """Add a filename that is computed using a user-specified function."""
        # Construct a function that will call the user supplied function with
        # the proper arguments. We prepend 'self' so the user supplied function
        # has easy access to all the filepaths.
        def fname(**kwargs):
            return func(self, **kwargs)

        # Bind the fname function to this instance of FileNames
        self.__dict__[alias] = fname


def _get_placeholders(template):
    """Get all placeholders from a template string.
    Parameters
    ----------
    template : str
        The template string to get the placeholders for.
    Returns
    -------
    placeholders : list of str
        The list of placeholder names that were found in the template string.

    Author: Marijn van Vliet <w.m.vanvliet@gmail.com>
    """
    return [p[1] for p in string.Formatter().parse(template)
            if p[1] is not None and len(p[1]) > 0]


def _substitute(template, files, user_values):
    """Makes a filename from a template.
    Any placeholders that point to known file aliases will be prefilled. The
    rest is filled given the values provided by the user when requesting the
    filename.
    Parameters
    ----------
    template : str
        The template string for the filename.
    files : list of str
        A list of file aliases that are already known.
    user_values : dict
        The key=value parameters that the user specified when requesting the
        filename.
    Returns
    -------
    filename : str
        The filename, obtained by filling all the placeholders of the template
        string.

    Author: Marijn van Vliet <w.m.vanvliet@gmail.com>
    """
    # Get all placeholder names
    placeholders = _get_placeholders(template)

    # Pre-fill placeholders based on existing file aliases
    placeholder_values = _prefill_placeholders(placeholders, files,
                                               user_values)

    # Add user specified values for the placeholders
    placeholder_values.update(**user_values)

    # Check whether all placeholder values are now properly provided.
    provided = set(placeholder_values.keys())
    needed = set(placeholders)
    missing = needed - provided
    if len(missing) > 0:
        raise ValueError('Cannot construct filename, because the following '
                         'parameters are missing: %s' % missing)

    # Do the substitution
    return template.format(**placeholder_values)


def _prefill_placeholders(placeholders, files, user_values):
    """Search through existing file aliases to pre-fill placeholder values.
    Parameters
    ----------
    placeholders : list of str
        The list of placeholder names that were found in the template string.
    files : list of str
        A list of file aliases that are already known.
    user_values : dict
        The key=value parameters that the user specified when requesting the
        filename. Can be empty if no parameters were specified (yet).
    Returns
    -------
    placeholder_values : dict
        A dictionary containing the values for the placeholders that could be
        pre-filled.

    Author: Marijn van Vliet <w.m.vanvliet@gmail.com>
    """
    placeholder_values = dict()

    for placeholder in placeholders:
        if placeholder in files:
            # Placeholder name is a filename, so get the path
            path = files[placeholder]
            if not isinstance(path, str):
                try:
                    path = path(**user_values)
                except ValueError:
                    # Placeholder could not be pre-filled given the supplied
                    # values by the user.
                    continue

            # Add the path as possible placeholder value
            placeholder_values[placeholder] = path

    return placeholder_values


def validate_sourcedata(path, source_type, pattern='sub-\\d+'):
    """
    This function validates the "sourcedata/" directory provided by user to
    see if it's contents are consistent with the BIDS specification.
    """

    if not path:
        path = './'

    if not source_type:
        source_type = ['eeg']

    # construct relative sourcedata path
    sourcedata_dir = os.path.join(path, 'sourcedata')

    source_name = None

    sub_dirs = []
    n_subs = 0
    subject_names = None

    data_dirs = []
    data_files = []

    # browse through directories in parent_dir
    for root, directories, files in os.walk(path):
        if root == path:
            if 'sourcedata' not in directories:
                raise ValueError('The provided directory does not contain "sourcedata/".')  # noqa: E501
            else:
                source_name = 'OK'

        if root == sourcedata_dir:
            sub_dirs.extend([os.path.join(root, sub) for sub in directories])
            n_subs = len(sub_dirs)

            valid_names = []
            for sub in sub_dirs:
                if re.findall(pattern=pattern,
                              string=os.path.basename(sub)):
                    valid_names.append(True)
                else:
                    valid_names.append(False)

            if len(set(valid_names)) > 1:
                bad_names = np.where(np.logical_not(valid_names))[0]
                bad_names = [sub_dirs[bn] for bn in bad_names]

                warnings.warn('The directory names in "sourcedata/" are not BIDS conform.')  # noqa: E501
                subject_names = 'error in %s' \
                                % (', '.join([str(bn) for bn in bad_names]))
            else:
                subject_names = 'OK'

        elif root in sub_dirs:
            data_dirs.append([data_dir for data_dir in directories
                              if data_dir in source_type])
            n_dirs = [len(d) for d in data_dirs]
            if len(set(n_dirs)) > 1:
                warnings.warn('Subject directories contain different number of subdirectories.')  # noqa: E501
                modal = max(set(n_dirs), key=n_dirs.count)

                inconscistent = [i for i, n in enumerate(n_dirs)
                                 if not n == modal]
                bads = [sub_dirs[i] for i in inconscistent]

        elif root in [os.path.join(s_dir, source)
                      for s_dir in sub_dirs for source in source_type]:
            data_files.extend([os.path.join(root, file) for file in files])

            valid_file_names = []
            file_patterns = [pattern + '_\\w+_' + source + '.*'
                             for source in source_type]

            for file in data_files:
                for file_pattern in file_patterns:
                    if re.findall(pattern=file_pattern,
                                  string=os.path.basename(file)):
                        valid_file_names.append(True)
                    else:
                        valid_file_names.append(False)

            if len(set(valid_file_names)) > 1:
                bad_files = np.where(np.logical_not(valid_file_names))[0]
                bad_files = [data_files[bf] for bf in bad_files]

                warnings.warn('Some file names are not BIDS conform.')
                file_names = 'error in %s' % (', '.join([str(bn)
                                                         for bn in bad_files]))
            else:
                file_names = 'OK'

    data_val = {
        'source_data_path':
            {
                'path': sourcedata_dir,
                'naming': source_name,
                'dirs_in_sourcedata': n_subs
            },
        'subject_directories':
            {
                'name_pattern': pattern,
                'naming': subject_names
            },
        'data_files':
            {
                'file_pattern': (', '.join([f_pat for f_pat in file_patterns])),  # noqa
                'naming': file_names  # noqa
            }
    }

    return data_val
