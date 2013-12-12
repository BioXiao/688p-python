Lecture 6: Software Best Practices
==================================
Keith Hughitt
2013/12/11

Overview
--------

This lecture describes some of the best practices to follow when writing
software for scientific work. Although the focus is on Python, the general
concepts apply equally to other programming environments.

Topic covered:

* documentation
* coding standards
* version control
* software testing
* debugging
* reproducible research
* data provenance

Documentation
-------------

## Writing documentation

Probably the most important thing you can do to make your code more useful both
to other people and to future self is to document it thoroughly -- Your code 
may still suck, but at least this way people have a good shot at understanding
what it's doing!

To effectively document code, think about creating documentation at several
different levels:

1. project
2. file
3. function
4. line

For the first three, you are basically providing an overview of the project
at an increasingly granular level.

In Python, you generally want to use triple-quotes (docstrings) for most
documentation. Only single-line comments should use "#".

### 1. Project documentation

This will vary depending on the scale of the project you are working on (here
I'm using "project" very loosely -- this could range from a short script you
put together to a full-scale software application.)

For smaller scripts, having a `README.txt` or a `README.md` is probably 
sufficient but for larger projects you will likely want to use a documentation 
generation tool such as `Sphinx`.

Some things to include:

* Project description (what does it do?)
* Contact information
* Date created / last modified
* Installation requirements
* Inputs and outputs
* Usage examples

### 2. File documentation

A couple lines to describe the contents of a file. If you are just working on
a small script then this can also take the form of your project description.

### 3. Function documention

For function documentation you want to include:

* A description of what the function is doing and possibly also a brief
  description of how it works.
* Input parameters and return values
* References, Usage examples, etc.

Example function doc from
[Pandas](https://github.com/pydata/pandas/blob/master/pandas/io/data.py):

```python
def DataReader(name, data_source=None, start=None, end=None,
               retry_count=3, pause=0.001):
    """
    Imports data from a number of online sources.

    Currently supports Yahoo! Finance, Google Finance, St. Louis FED (FRED)
    and Kenneth French's data library.

    Parameters
    ----------
    name : str or list of strs
        the name of the dataset. Some data sources (yahoo, google, fred) will
        accept a list of names.
    data_source: str
        the data source ("yahoo", "google", "fred", or "ff")
    start : {datetime, None}
        left boundary for range (defaults to 1/1/2010)
    end : {datetime, None}
        right boundary for range (defaults to today)

    Examples
    ----------

    # Data from Yahoo! Finance
    gs = DataReader("GS", "yahoo")

    # Data from Google Finance
    aapl = DataReader("AAPL", "google")

    # Data from FRED
    vix = DataReader("VIXCLS", "fred")

    # Data from Fama/French
    ff = DataReader("F-F_Research_Data_Factors", "famafrench")
    ff = DataReader("F-F_Research_Data_Factors_weekly", "famafrench")
    ff = DataReader("6_Portfolios_2x3", "famafrench")
    ff = DataReader("F-F_ST_Reversal_Factor", "famafrench")
    """
```

### 4. Line comments

Line comments are the most granular level of documentation. You don't have to
comment every single line of code, but you should consider including a comment
whenever you are doing something complex or not immediately obvious, or to
describe the basic flow of a function or script.

Examples:

```python
# filter out non-gene entries
gene_rows = [ x for x in lines if 'gene\t' in x]

# gene id
match = re.search('ID=[^;]+;', row['attributes'])
gene_id = match.group()[3:-1]
```

## Documentation generation

Finally, for larger scale projects, you may want to consider using a tool to
help generate your documentation from your code comments.

For Python projects, [Sphinx](http://sphinx-doc.org/) is the preferred tool
and it is what was used to create the documentation for NumPy, Matplotlib,
Pandas, etc.

To use it, you simply format your file and code comments in a specific manner.
Sphinx can then scan through all of your files, extract the comments, and
parse them into elegant-looking HTML or PDF documentation.

Coding Standards
----------------

A coding standard is a set of rules or guidelines to follow when writing
code and documentation for a project which, when followed, help to ensure a
certain level of code quality and consistency.

In addition to setting some basic conventions to ensure consistency within and
across files, many of the guidelines are chosen for specific reasons such
as readability, and ease of debugging.

Whenever you are starting to work with a new programming language, it is a good
idea read a bit about coding standard and best practices for that language.
Each language will have its own conventions and you can learn a lot about a
language from reading about these. Also, the earlier you adopt a standard, the
easier it will be to develop those habits early on.

The Python community largely follows a coding standard referred to as 
[PEP8](http://www.python.org/dev/peps/pep-0008/) and there are [numerous 
tools](&oq=pep8+checker&aqs=chrome..69i57j0l3.2103j0j1&sourceid=chrome&ie=UTF-8)
to help check if your code is following these conventions.

In fact, Spyder has [an option to check for PEP8 compliance built-in](http://www.southampton.ac.uk/~fangohr/blog/spyder-the-python-ide.html#warn-if-pep8-coding-guidelines-are-violated):
enable it!

Version Control
---------------

Version control software ("VCS") are tools which help to track changes made
to a collection of files over time, and to assist in the collaborative editing
of these files.

There are a number of tools that do this, but the main ones in use now are:

* git
* subversion (svn)
* cvs
* mercurial
* bazaar

Some of the reasons for using a VCS are:

### File history

ProjectKeeping track of the history of changes made to each file in a project. 
This way if something breaks, an earlier version of the file can be restored
and used to track down the problematic changes.

### Collaboration

VCS help to provide mechanisms for multiple developers to work on different
parts of a codebase and merge those changes together.

### Backup

VCS tools generally provide a mechanism to store a codebase at some central
location (svn and cvs), or to keep multiple copies of the entire codebase
at different locations (distributed VCS).

Many websites (e.g. Github) provide free hosting for open-source code stored
in a VCS repository.

This way if your laptop dies, or you accidentally delete an important directory
you will still have a backup of the code and all of the changes made to it
available online.

### More info:

* [Git: Getting Started - About Version Control](http://git-scm.com/book/en/Getting-Started-About-Version-Control)
* [What Is Git & Why You Should Use Version Control If You’re a Developer](http://www.makeuseof.com/tag/git-version-control-youre-developer/)
* [Github for beginners](http://readwrite.com/2013/09/30/understanding-github-a-journey-for-beginners-part-1#awesm=~opJ2C5U2izU4OD)

Software Testing
----------------

Thoroughly [testing](http://en.wikipedia.org/wiki/Software_testing) your
software is one way to ensure find and eliminate bugs and prevents issues
from creeping into code as it and the software it depends on changes over
time.

There are many different types of tests for software but the most common kind
and the easiest to begin using is called [unit
testing](http://en.wikipedia.org/wiki/Unit_testing). The basic idea for unit
testing is to write one or more "test" functions for each function you want
to test where you describe some inputs to the function along with know
(expected) outputs. A [test framework](https://wiki.python.org/moin/PythonTestingToolsTaxonomy)
is then used to scan your code for these tests, apply them to your functions, 
and ensure that the output is as expected. Testing frameworks exist for all
different programming languages, and Python is no exception. In addition to
some [builtin](http://docs.python.org/2/library/doctest.html) [methods](http://docs.python.org/2/library/unittest.html) 
for testing, there are several popular test frameworks including:

* [nose](http://nose.readthedocs.org/en/latest/)
* [py.test](http://pytest.org/latest/)

For more information on testing Python code, see:

* http://pythontesting.net/start-here/
* http://docs.python-guide.org/en/latest/writing/tests/

Debugging
---------

Although adopting the above practices will almost certainly save you some
heartache in the future, you will still encounter issues during coding that
you need to resolve. Although learning how to properly debug software may
initially seem daunting and takes some time to get used to, becoming 
comfortable this process will make your life much easier and save you a lot
of time in the long run.

Here are a few references that you may find useful when you encounter problems
during Python coding:

* [Python Development with IPython: Debugging Python using the IPython Shell](http://www.youtube.com/watch?v=zm8tX0e2kyI)
* [IPython reference: using IPython for interactive work](http://ipython.org/ipython-doc/rel-1.1.0/interactive/reference.html)
* [pdb — The Python Debugger](http://docs.python.org/2/library/pdb.html)
* [Stack Overflow - Python debugging tips](http://stackoverflow.com/questions/1623039/python-debugging-tips)


Reproducible Research
---------------------

See: https://github.com/umd-byob/presentations/tree/master/2013/0903-knitr_reproducible_research


Data Provenance
---------------

The basic idea here is to try and include additional metadata with any data
you create in such a way that future users of that data product will be able
to determine where the data came from, how the data was produced, and how they
could recreate the data themselves. This is related to the idea of reproducible
research above. Often, data files (e.g. images, csv tables, etc.) will be
passed along from one location to another, occasionally changing in filename or
contents. If proper provenance information is included in the file, for example
as a comment at the top of the file, then the receiver of the data should have
sufficient information to understand the where the data came from and how is
changed from the time of creation (post-processing steps, etc).

Some information to consider including:

* A brief description of the file contents
* Dates of creation and modification
* Contact information
* How the data was generated, including versions and parameters of any software
  called
* Any changes that have been made to the file since its original creation

