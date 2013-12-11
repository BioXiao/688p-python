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

```
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


