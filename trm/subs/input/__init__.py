"""
a module for parameter input to scripts

This module handles parameter input from the user and storage and
retrieval of parameters from disk files. This gives scripts a memory
and can save a lot of typing, especially with frequently invoked
scripts.

Classes
=======

Input      -- the main class for parameter input
InputError -- exception class for the package

Functions
=========

clist      -- split up a command string appropriately for Input

Examples of parameter input
===========================

Here are some examples of usage to illustrate this:

A command with inputs 'device' (hidden), 'npoint' and 'output' could be 
invoked variously as

command

(npoint and output will be prompted for)

or 

command device=/ps npoint=20

(output will be prompted for)

or 

command device=/ps \\

(\ indicates take default values for npoint and output, in UNIX shells it must be escaped
hence \\)

or

command 20

(npoint will be set = 20, output will be prompted for). Note that such
unnamed parameters can only set the values of parameters which are by default
prompted for. Hidden parameters must always be explicitly named to specify
them on the command line. This prevents differing behaviour when the prompt flag
'?' is specified.

'?' or 'PROMPT' in the argument list (quotes for emphasis) force prompting for all
variables.  On UNIX shells you may have to escape ? as in \? to get it through or
quote it: '?'. 'NODEFS' bypasses any attempt to read or write the default files.
It is mainly provided as a way to avoid clashes between multiple processes but may also
speed commands a tiny bit.

When you get prompting, <tab> allows you to complete filenames. Entering '?' gives
the parameter range if any has been supplied.
"""
from __future__ import print_function
import os, re, sys, pickle, warnings

# next two lines allow tab completion of file names
import readline
readline.parse_and_bind("tab: complete")

#def complete(text,state):
#    results = ["example",None]
#    return results[state]

#readline.set_completer(complete)

def add_extension(fname, ext):
    """Add extension ext to a file name if it is not already there, and returns
    the revised name

    """
    if len(ext) and not fname.endswith(ext):
        return fname + ext
    else:
        return fname

def clist(command):
    """
    Splits up a command string returning a list suitable for
    constructing Cline objects. The reason for using this rather than
    a simple string split is that it allows you to use double quotes
    to get strings with spaces through. Returns a list of strings.
    """

    cl = re.findall('\"[^"]*\"|\S+', command)
    return [c.lstrip('"').rstrip('"') for c in cl]


class Input:

    """Class to handle command line inputs. In particular this allows storage and
    retrieval of input parameter values from files which allows different
    scripts to communicate parameters to each other through 'global' defaults,
    and commands to have a 'memory' between different invocations. To use the
    class you first create an instance, then register each parameter name, and
    finally get the input, either from the user, default values or disk.

    Here is some example code::

      >> import input as inp
      >>
      >> # Initialize Input. COMM_ENV is an environment
      >> # variable specifying a directory where the files
      >> # are stored. '.comm' is the name of a directory
      >> # under the home directory that will be used by
      >> # default.
      >>
      >> cline = inp.Input('COMM_ENV', '.comm', 'command', args)
      >>
      >> # register parameters
      >> cline.register('device', inp.Input.GLOBAL, inp.Input.HIDE)
      >> cline.register('npoint', inp.Input.LOCAL,  inp.Input.PROMPT)
      >> cline.register('output', inp.Input.LOCAL,  inp.Input.PROMPT)
      >>
      >> try:
      >>    device = cline.get_value('device', 'plot device', '/xs')
      >>    npoint = cline.get_value('npoint', 'number of points', 10, 1, 100)
      >>    output = cline.get_value('output', 'output file', 'save.dat')
      >> except inp.InputError, err:
      >>    print('Error on parameter input:')
      >>    print(err)
      >>    exit(1)
      >>
      >> # rest of program here ...

    :class:`Input` objects define the four static variables GLOBAL, LOCAL,
    PROMPT, HIDE which should be used when registering parameters to define
    their properties

    """

    # All attributes are hidden, so not publicly documented. They are::
    #
    #    _ddir    -- name of defaults directory (string)
    #    _lname   -- name of local defaults file (string)
    #    _gname   -- name of global defaults file (string)
    #    _lpars   -- parameters loaded from local default file (dict)
    #    _gpars   -- parameters loaded from global default file (dict)
    #    _cname   -- command name (string)
    #    _pbynam  -- parameter/value pairs read from arguments (dict)
    #    _pbypos  -- parameter values by position read from arguments (list)
    #    _rpars   -- List of registered parameters. For each one a dictionary
    #                specifies whether they are to be found in the global or
    #                local default and whether they should be prompted for or
    #                not.
    #    _prompt  -- force prompting or not (bool)
    #    _list    -- list the parameter name / value pairs or not (bool)
    #    _nodefs  -- controls whether disk default files will be accessed or not
    #                (True ==> no access) (bool)
    #    _usedef  -- use default rather than prompt from user. (bool)

    GLOBAL   = 1
    LOCAL    = 2
    PROMPT   = 3
    HIDE     = 4

    def __init__(self, direnv, defdir, args):
        """
        Initialize an Input class object

        Arguments::

           direnv : (string)
              environment variable pointing at a directory where default
              files will be stored.

           defdir : (string)
              default directory (sub-directory of HOME) if the enviroment
              variable 'direnv' is not defined

           args : (list of strings)
              command-line arguments. The first one must be the command name.

        """

        # Flag to indicate whether inputs need to be saved
        self._save = False

        # Extract special keywords NODEFS, PROMPT and LIST from argument list

        # Store the command name
        self._cname = args.pop(0)

        if 'NODEFS' in args:
            self._nodefs = True
            args = filter(lambda a: a != 'NODEFS', args)
        else:
            self._nodefs = False

        if 'PROMPT' in args:
            self._prompt = True
            args = filter(lambda a: a != 'PROMPT', args)
        else:
            self._prompt = False

        if 'LIST' in args:
            self._list = True
            args = filter(lambda a: a != 'LIST', args)
        else:
            self._list = False

        if self._list:
            print ('\n' + self._cname)

        if not self._nodefs:

            if direnv is None and defdir is None:
                raise InputError(
                    'input.Input: no default file environment variable or directory name supplied')

            if direnv is not None and direnv in os.environ:
                self._ddir = os.environ[direnv]
            else:
                home = os.environ['HOME']
                self._ddir = os.path.join(home, defdir)

            # read local and global default files
            self._lname = os.path.join(self._ddir, os.path.basename(self._cname) + '.def')
            try:
                with open(self._lname,'rb') as flocal:
                    self._lpars  = pickle.load(flocal)
            except IOError:
                self._lpars = {}
            except (EOFError, pickle.UnpicklingError):
                warnings.warn(
                    'input.Input: failed to read local defaults file {0:s}; possible corrupted file\n'.format(self._lname), UserWarning)
                self._lpars = {}

            self._gname = os.path.join(self._ddir, 'GLOBAL.def')
            try:
                with open(self._gname,'rb') as fglobal:
                    self._gpars  = pickle.load(fglobal)
            except IOError:
                self._gpars = {}
            except (EOFError, pickle.UnpicklingError):
                warnings.warn(
                    'input.Input: failed to read global defaults file {0:s}; possible corrupted file.\n'.format(self._gname),UserWarning)
                self._gpars = {}

        else:
            self._ddir = None
            self._lpars = {}
            self._gpars = {}

        # _pbynam and _pbypos are a dictionary of parameters defined by name
        # and a list of parameters defined by position sent in through the
        # command line arguments
        self._pbynam = {}
        self._pbypos = []
        checker = re.compile('[a-zA-Z0-9]+=')
        for arg in args:
            if checker.match(arg):
                p,v = arg.split('=',1)
                if p in self._pbynam:
                    raise InputError(
                        'input.Input: parameter = ' + p + ' defined more than once in argument list.')
                self._pbynam[p] = v
            else:
                self._pbypos.append(arg)

        self._rpars = {}
        self.narg = 0
        self._usedef = False

    def prompt_state(self):
        """Says whether prompting is being forced or not. Note the propting state does
        not change once an Input is initialized, being fixed by the presence
        of 'PROMPT' on the command line or not.

        Returns True if prompting is being forced.

        """
        return self._prompt

    def register(self, param, g_or_l, p_or_h):
        """Registers a parameter as one to be expected and defines basic
        properties. You must call this once for every parameter that you might
        call 'get_value' for.

        Arguments::

          param : (string)
              parameter name. Must have no spaces, equal signs or quotes.

          g_or_l : (int)
              defines whether the parameter should be global, i.e. stored
              in a file called GLOBAL to allow access from other commands,
              or just local to this command. Use the static variables GLOBAL
              and LOCAL to set this, e.g. input.Input.GLOBAL

          p_or_h : (int)
              defines whether the parameter is prompted for by default or
              hidden. Parameters that are rarely changed are better hidden to
              reduce clutter. The PROMPT command-line keyword forces even
              hidden parameters to be prompted for in the rare cases that they
              need to be changed. Use the static variables PROMPT and HIDE to
              set this.

        Sometimes you may want to set the default values of hidden parameters
        unless you are happy for the set value to be retained.

        """

        if param.find(' ') != -1 or param.find('\t') != -1 or param.find('=') != -1 or param.find('"') != -1 or param.find("'") != -1:
            raise InputError('Parameter = ' + param + ' is illegal.')

        if g_or_l != Input.GLOBAL and g_or_l != Input.LOCAL:
            raise InputError('g_or_l must either be Input.GLOBAL or Input.LOCAL')

        if p_or_h != Input.PROMPT and p_or_h != Input.HIDE:
            raise InputError('p_or_h must either be Input.PROMPT or Input.HIDE')

        if param in self._rpars:
            raise InputError('parameter = ' + param + ' has already been registered.')

        self._rpars[param] = {'g_or_l' : g_or_l, 'p_or_h' : p_or_h}

    def set_default(self, param, defval):
        """Set the default value of a parameter automatically. This is often useful
        for changing hidden parameters on the fly.

        """
        if param not in self._rpars:
            raise InputError('set_default: parameter = "' + param + '" has not been registered.')

        if self._rpars[param]['g_or_l'] == Input.GLOBAL:
            self._gpars[param] = defval
        else:
            self._lpars[param] = defval

    def get_value(self, param, prompt, defval, minval=None, maxval=None,
                  lvals=None, fixlen=True, multipleof=None):
        """Gets the value of a parameter, either from the command arguments, or by
        retrieving default values or by prompting the user as required. This
        is the main function of Input. The value obtained is used to update
        the defaults which, if 'NODEFS' has not been defined, are written to
        disk at the end of the command.

        Arguments::

          param : (string)
             parameter name.

          prompt : (string)
             the prompt string associated with the parameter

          defval : (various)
             default value if no other source can be found (e.g. at start).
             This also defines the data type of the parameter (see below for
             possibilities)

        minval : (same as defval's type)
             the minimum value of the parameter to allow.

        maxval : (same ad defval's type)
             the maximum value of the parameter to allow.

        lvals : (list)
             list of possible values (exact matching used)

        fixlen : (bool)
             for lists or tuples, this insists that the user input has the
             same length

        multipleof : (int)
             specifies a number that the final value must be a multiple of
             (integers only)

        Data types: at the moment, only certain data types are recognised by
        this routine. These are the standard numerical types, 'int', 'float',
        the logical type 'bool' which can be set with any of (case
        insensitively) 'true', 'yes', 'y', '1' (all True), or 'false', 'no',
        'n', '0' (all False), strings, and Fname objects to represent
        filenames with specific extensions, and lists. In the case of tuples,
        it is the default value 'defval' which sets the type.

        """

        # calling this routine indicates that we need to save the inputs at
        # some point
        self._save = True
        if param not in self._rpars:
            raise InputError(
                'input.Input.get_value: parameter = "{:s}" has not been registered.'.format(param.upper())
                )

        if lvals != None and defval not in lvals:
            raise InputError(
                'input.Input.get_value: default = {!s} not in allowed list = {!s}'.format(defval,lvals)
            )

        # Now get the parameter value by one of three methods
        if param in self._pbynam:
            # get value from name/value pairs read from command line arguments
            # of the form param=value
            value = self._pbynam[param]

        elif self.narg < len(self._pbypos) and self._rpars[param]['p_or_h'] == Input.PROMPT:
            # get value from bare values in the command line such as '23' '\\'
            # indicates use the default value and also to use defaults for any
            # other unspecified parameters that come later (_usedef set to
            # True)
            if self._pbypos[self.narg] == '\\':
                if self._rpars[param]['g_or_l'] == Input.GLOBAL and param in self._gpars:
                    value = self._gpars[param]
                elif self._rpars[param]['g_or_l'] == Input.LOCAL and param in self._lpars:
                    value = self._lpars[param]
                else:
                    value = defval
                self._usedef = True
            else:
                value = self._pbypos[self.narg]
            self.narg += 1

        else:
            # load default from values read from file or the initial value
            if self._rpars[param]['g_or_l'] == Input.GLOBAL and param in self._gpars:
                value = self._gpars[param]
            elif self._rpars[param]['g_or_l'] == Input.LOCAL and param in self._lpars:
                value = self._lpars[param]
            else:
                value = defval

            # prompt user for input
            if not self._usedef and \
               (self._prompt or self._rpars[param]['p_or_h'] == Input.PROMPT):
                reply = '?'
                while reply == '?':
                    reply = input(
                        '{:s} - {:s} [{!s}]: '.format(param,prompt,value)
                    )
                    if reply == '\\':
                        self._usedef = True
                    elif reply == '?':
                        print()
                        if minval != None and maxval != None:
                            print('Parameter = "{:s}" must lie from {!s} to {!s}'.format(param,minval,maxval))
                        elif minval != None:
                            print('Parameter = "{:s}" must be greater than {!s}'.format(param,minval))
                        elif maxval != None:
                            print('Parameter = "{:s}" must be less than {!s}'.format(param,maxval))
                        else:
                            print('Parameter = "{:s}" has no restriction on its value'.format(param))

                        print('"{:s}" has data type = {!s}'.format(param,type(defval)))
                        if lvals != None:
                            print('Only the following values are allowed:')
                            print(lvals)
                        if isinstance(defval, (list, tuple)) and fixlen:
                            print('You must enter exactly {:d} values'.format(len(defval)))
                        print()
                    elif reply != '':
                        if  isinstance(defval, (list, tuple)) and fixlen and len(reply.split()) != len(defval):
                            print ('You must enter exactly {:d} values. [You only entered {:d}]'.format(len(defval),len(reply.split())))
                            reply = '?'
                        else:
                            value = reply

        # at this stage we have the value, now try to convert to the right
        # type according to the type of 'defval'
        try:
            if isinstance(defval, Fname):
                value = defval(value)
            elif isinstance(defval, str):
                value = str(value)
            elif isinstance(defval, bool):
                if isinstance(value, str):
                    if value.lower() == 'true' or value.lower() == 'yes' or \
                       value.lower() == '1' or value.lower() == 'y':
                        value = True
                    elif value.lower() == 'false' or value.lower() == 'no' or \
                         value.lower() == '0' or value.lower() == 'n':
                        value = False
                    else:
                        raise InputError(
                            'input.Input.get_value: could not translate "' + value + '" to a boolean True or False.')
            elif isinstance(defval, int):
                value = int(value)
            elif isinstance(defval, float):
                value = float(value)
            elif isinstance(defval, list):
                if isinstance(value, str):
                    value = map(type(defval[0]), value.split())
                else:
                    value = list(value)
            elif isinstance(defval, tuple):
                if isinstance(value, str):
                    value = tuple(map(type(defval[0]), value.split()))
                else:
                    value = tuple(value)
            else:
                raise InputError('input.Input.get_value: did not recognize the data type of the default supplied for parameter {:s} = {!s}'.format(param, type(defval)))

        except ValueError as err:
            raise InputError('input.Input.get_value: {!s}'.format(err))

        # ensure value is within range
        if minval != None and value < minval:
            raise InputError(
                'input.Input.get_value: ' +
                param + ' = ' + str(value) + ' < ' + str(minval))
        elif maxval != None and value > maxval:
            raise InputError(
                'input.Input.get_value: ' +
                param + ' = ' + str(value) + ' > ' + str(maxval))

        # and that it is an OK value
        if lvals != None and value not in lvals:
            raise InputError(
                'input.Input.get_value: ' +
                str(value) + ' is not one of the allowed values = ' + str(lvals))

        if multipleof != None and value % multipleof != 0:
            raise InputError(
                'input.Input.get_value: ' +
                str(value) + ' is not a multiple of ' + str(multipleof))

        # update appropriate set of defaults. In the case of Fnames, strip the
        # extension
        if self._rpars[param]['g_or_l'] == Input.GLOBAL:
            if isinstance(defval, Fname):
                self._gpars[param] = defval.noext(value)
            else:
                self._gpars[param] = value
        else:
            if isinstance(defval, Fname):
                self._lpars[param] = defval.noext(value)
            else:
                self._lpars[param] = value

        if self._list:
            print (param,'=',value)

        return value

    def get_rest(self):
        """
        Returns any unused command-line arguments as a list or None
        if there aren't any
        """
        if self.narg < len(self._pbypos):
            return self._pbypos[self.narg:]
        else:
            return None

    def __del__(self):
        if self._save:
            warnings.warn(
                'input.Input.__del__: an input or inputs has/have not been saved; possible programming error', UserWarning)

    def save(self):
        """Saves parameter values to disk (if NODEFS has not been specified). This
        should always be called after the last input in order to save the
        values for next time.

        """

        if not self._nodefs:

            # make the default directory if need be
            if not os.path.lexists(self._ddir):
                os.mkdir(self._ddir, 0o755)

            # save local defaults
            with open(self._lname, 'wb') as flocal:
                pickle.dump(self._lpars, flocal)

            # save global defaults
            with open(self._gname, 'wb') as fglobal:
                pickle.dump(self._gpars, fglobal)

        self._save = False

class Fname(str):
    """Defines a callable parameter type for the :class:`Cline` to allow for some
    early checks on file names. This is mainly to prevent a whole series of
    parameters being input, only to find that the file name input in the first
    one is incorrect.

    """

    OLD       = 0
    NEW       = 1
    NOCLOBBER = 2

    def __new__(cls, root, ext='', ftype=OLD, exist=True):
        """Constructor distinct from __init__ because str is immutable. In the
        following text items in capitals such as 'OLD' are static variables so
        that one should use input.Fname.OLD or equivalent to refer to
        them.

        Arguments::

          root : (string)
             root name of file (if it ends with 'ext', an extra 'ext' will
             not be added)

          ext : (string)
             extension, e.g. '.dat'

          ftype : (int)
             OLD = existing or possibly existing file; NEW = new file which
             will overwrite anything existing; NOCLOBBER is same as NEW but
             there must not be an existing one of the specified name.

          exist : (bool)
             If exist=True and ftype=OLD, the file :must: exist. If
             exist=False, the file may or may not exist already.

        """

        if ftype != Fname.OLD and ftype != Fname.NEW and ftype != Fname.NOCLOBBER:
            raise Exception(
                'input.Fname.__new__: ftype must be either OLD, NEW or NOCLOBBER')

        # store root with no extension
        if len(ext) and root.endswith(ext):
            fname = super().__new__(cls, root[:-len(ext)])
        else:
            fname = super().__new__(cls, root)

        return fname

    def __init__(self, root, ext='', ftype=OLD, exist=True):
        """Initialiser. In the following text items in capitals such as 'OLD' are
        static variables so that one should use input.Fname.OLD or
        equivalent to refer to them.

        Arguments::

          root : (string)
             root name of file (if it ends with 'ext', an extra 'ext' will
             not be added)

          ext : (string)
             extension, e.g. '.dat'

          ftype : (int)
             If exist=True and ftype=OLD, the file :must: exist. If exist=False, the file may or may
             not exist already.

          exist : (bool)
             If True, the file must exist.

        ext, ftype and exist are stored as identically-named attributes. 'root'
        is stored as the base string.

        """

        self.ext = ext
        self.ftype = ftype
        self.exist = exist

    def __call__(self, fname):
        """Given a potential file name, this first ensures that it has the correct
        extension, and then tests for its existence if need be, depending upon
        the values of `ftype` and `exist` defined at instantiation.

        Arguments::

           fname : (string)

              file name. The extension associated with the :class:`Fname` will
              be added if necessary.

        Returns the file name to use. Raises an exception if there
        are problems.

        """

        # Add extension if not already present.
        fname = add_extension(fname, self.ext)

        if self.exist and self.ftype == Fname.OLD and not os.path.exists(fname):
            raise Exception(
                'input.Fname.__call__: could not find file = ' + fname)

        if self.ftype == Fname.NOCLOBBER and os.path.exists(fname):
            raise Exception(
                'input.Fname.__call__: file = ' + fname + ' already exists')

        return fname

    def noext(self, fname):
        """Returns the suggested file name, `fname`, with the extension removed"""
        if len(self.ext) and fname.endswith(self.ext):
            return fname[:-len(self.ext)]
        else:
            return fname

    def exists(self):
        """
        Checks that the Fname file exists
        """
        return os.path.exists(self)

    def __getnewargs__(self):

        """Enables pickling of :class:`Fname` objects. This returns a tuple of
        arguments that are passed off to __new__

        """
        return (self,self.ext,self.ftype,self.exist)

class InputError(Exception):
    """For throwing exceptions from the input module"""
